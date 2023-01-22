"""
Contains the main object Fusor for running MaxFuse pipeline
"""

import warnings
from collections.abc import Iterable
from collections import defaultdict
import numpy as np
from sklearn.cross_decomposition import CCA

# import bisect
# import ruptures as rpt

from . import match_utils, graph, utils

import matplotlib.pyplot as plt
from matplotlib.ticker import MaxNLocator


class Fusor:
    """
    Main object for fitting MaxFuse.
    """
    def __init__(
            self, shared_arr1, shared_arr2, active_arr1, active_arr2,
            method='centroid_shrinkage', labels1=None, labels2=None
    ):
        """
        Initialization for Fusor object.

        :param shared_arr1: First dataset with shared features.
        :type shared_arr1: np.ndarray
        :param shared_arr2: Second dataset with shared features.
        :type shared_arr2: np.ndarray
        :param active_arr1: First dataset with all active features.
        :type active_arr1: np.ndarray
        :param active_arr2: Second dataset with all active features.
        :type active_arr2: np.ndarray
        :param method: One of 'centroid_shrinkage' or 'graph_smoothing',
            controlling how the fuzzy smoothing is done.
        :type method: str
        :param labels1: Optional cluster labels for the first dataset.
        :type labels1: Union[np.ndarray, List, None]
        :param labels2: Optional cluster labels for the second dataset
        :type labels2: Union[np.ndarray, List, None]
        """
        # input
        self.shared_arr1 = shared_arr1
        self.shared_arr2 = shared_arr2
        self.active_arr1 = active_arr1
        self.active_arr2 = active_arr2
        assert method in {'centroid_shrinkage', 'graph_smoothing'}
        self.method = method
        self.labels1 = labels1
        self.labels2 = labels2

        # batching
        # linear assignment will be performed on a n1*n2 bipartite graph
        # self.max_outward_size specifies the max n1 in one batch
        self.max_outward_size = None
        # self.max_pivot_size * self.matching_ratio is the max n2 in one batch
        self.matching_ratio = None
        # arr1 will be clustered into len(self.active_arr1) / self.metacell_size many clusters
        # and the cluster centroids are the metacells
        self.metacell_size = None
        # store indices for each batch
        self._batch_to_indices1 = None
        self._batch_to_indices2 = None
        # store the map between batches in each round of matching
        self._batch1_to_batch2 = None

        # graph construction and clustering
        # self._edges1[b] is the graph constructed on the b-th batch of arr1
        # self._edges1[b] = rows, cols, vals, each edge is (rows[i], cols[i]), and its weight is vals[i]
        self._edges1 = None
        self._edges2 = None
        # self._labels1[b] is the clustering labels constructed on the b-th batch of arr1
        self._labels1 = None
        self._labels2 = None
        # self._metacell_labels1[b] is the metacell labels constructed on the b-th batch of arr1
        self._metacell_labels1 = None

        # matching
        # for batch index b, the matching is self._init_matching[b]
        # the indices in the matching is the row indices of the data at this batch
        # if metacells are used, it will be the metacell indices
        self._init_matching = None
        self._refined_matching = None
        # the matching for batch b that remained after filtering is
        # self._refined_matching[b][0][self._remaining_indices_in_refined_matching[b]]
        # self._refined_matching[b][1][self._remaining_indices_in_refined_matching[b]]
        # self._refined_matching[b][2][self._remaining_indices_in_refined_matching[b]]
        self._remaining_indices_in_refined_matching = None
        self._propagated_matching = None
        # store the propagated matching after filtering in the same way as the refined matching
        self._remaining_indices_in_propagated_matching = None

        # hyperparameters regarding refined matching
        self._svd_components1_for_cca_embedding = None
        self._svd_components2_for_cca_embedding = None
        self._randomized_svd_for_cca_embedding = None
        self._svd_runs_for_cca_embedding = None
        self._cca_components = None
        self._cca_max_iter = None

        # memorizing info for dimension reduction
        self._cca_on_pivots = None
        self._rotation_before_cca_on_active_arr1 = None
        self._rotation_before_cca_on_active_arr2 = None

        # output
        # both are defaultdict(list)
        # structure is {(idx2: [(idx1, score1), (idx1', score1'), ...])}
        # i.e., idx2 is matched to multiple indices1, and their similarity scores (not distances) are stored
        self._pivot2_to_pivots1 = None
        self._propidx2_to_propindices1 = None

    def split_into_batches(
            self, max_outward_size=5000, matching_ratio=3, metacell_size=2,
            method='random', batching_scheme='pairwise',
            prebatching_smoothing=False,
            shared_wt1=1, shared_wt2=1,
            active_wt1=1, active_wt2=1,
            shared_svd_components1=None, shared_svd_components2=None,
            active_svd_components1=None, active_svd_components2=None,
            randomized_svd=False, svd_runs=1,
            seed=None, verbose=True
    ):
        """
        Split the data into batches.

        Parameters
        ----------
        max_outward_size: int, default=10000
            Max number of cells to match in arr1.
        matching_ratio: int, default=3
            One cell in arr1 is matched to how many cells in arr2 on average.
        metacell_size: int, default=2
            For arr1, how many cells will be aggregated into one metacell on average.
        method: str, default='random'
            Either 'random', doing random split, or 'binning', doing binning by largest singular vector.
        batching_scheme: str, default='pairwise'
            Either 'cyclic' or 'pairwise'
                if cyclic, pair batches in the two datasets in a cyclic fashion
                if pairwise, all possible combinations of pairs of batches are considered
        prebatching_smoothing: bool, default=False
            Whether to smooth towards the low rank approximation before batching
        shared_wt1: float, default=1
            The shrinkage weight to put on the raw data for shared_arr1.
        shared_wt2: float, default=1
            The shrinkage weight to put on the raw data for shared_arr2.
        active_wt1: float, default=1
            The shrinkage weight to put on the raw data for active_arr1.
        active_wt2: float, default=1
            The shrinkage weight to put on the raw data for active_arr2.
        shared_svd_components1: None or int, default=None
            If not None, perform SVD to reduce the dimension of self.shared_arr1.
        shared_svd_components2: None or int, default=None
            If not None, perform SVD to reduce the dimension of self.shared_arr2.
        active_svd_components1: None or int, default=None
            If not None, perform SVD to reduce the dimension of self.active_arr1.
        active_svd_components2: None or int, default=None
            If not None, perform SVD to reduce the dimension of self.active_arr2.
        randomized_svd: bool, default=False
            Whether to use randomized SVD.
        svd_runs: int, default=1
            Perform multiple runs of SVD and the one with lowest Frobenious reconstruction error is selected.
        seed: None or int, default=None
            Numpy random seed.
        verbose: bool, default=True
            Whether to print the progress.

        Returns
        -------
        None
        """
        def smoothing_towards_low_rank_approximation(arr, wt, n_components, randomized, n_runs):
            """Return wt * arr + (1-wt) * low rank approximation of arr.
            """
            approx = utils.svd_denoise(arr=arr, n_components=n_components, randomized=randomized, n_runs=n_runs)
            return wt * arr + (1-wt) * approx

        if prebatching_smoothing:
            print('Smoothing towards low rank approximations...', flush=True)
            # smoothing shared_arr and active_arr towards their low rank approximation
            self.shared_arr1 = smoothing_towards_low_rank_approximation(
                arr=self.shared_arr1, wt=shared_wt1, n_components=shared_svd_components1,
                randomized=randomized_svd, n_runs=svd_runs
            )
            self.shared_arr2 = smoothing_towards_low_rank_approximation(
                arr=self.shared_arr2, wt=shared_wt2, n_components=shared_svd_components2,
                randomized=randomized_svd, n_runs=svd_runs
            )
            self.active_arr1 = smoothing_towards_low_rank_approximation(
                arr=self.active_arr1, wt=active_wt1, n_components=active_svd_components1,
                randomized=randomized_svd, n_runs=svd_runs
            )
            self.active_arr2 = smoothing_towards_low_rank_approximation(
                arr=self.active_arr2, wt=active_wt2, n_components=active_svd_components2,
                randomized=randomized_svd, n_runs=svd_runs
            )

        self.max_outward_size = max_outward_size
        self.matching_ratio = matching_ratio
        self.metacell_size = metacell_size

        if method == 'random':
            def _split(arr, n_batches, curr_seed):
                indices = list(np.random.RandomState(curr_seed).permutation(arr.shape[0]))
                res = []
                batch_size = int(len(indices) // n_batches)
                for b in range(n_batches):
                    res.append(indices[b * batch_size:(b + 1) * batch_size])
                res[-1].extend(indices[n_batches * batch_size:])
                return res

        elif method == 'binning':
            def _split(arr, n_batches, curr_seed):
                if curr_seed is not None:
                    np.random.seed(curr_seed)
                largest_eigenvec = utils.svd_embedding(arr, n_components=1).flatten()
                # sort the largest eigen vector
                ordered_indices = np.argsort(largest_eigenvec)
                batch_size = int(arr.shape[0] // n_batches)
                # partition ordered_indices into batch_size many parts
                partitions = []
                for i in range(batch_size):
                    partitions.append(ordered_indices[i * n_batches:(i+1) * n_batches])
                    partitions[-1] = list(partitions[-1][np.random.permutation(len(partitions[-1]))])
                extra_part = ordered_indices[batch_size * n_batches:]

                # select one point from each part in partitions
                batch_to_indices = [[] for _ in range(n_batches)]
                for b in range(n_batches):
                    for part in partitions:
                        # print(part)
                        batch_to_indices[b].append(part[b])
                batch_to_indices[-1].extend(extra_part)
                return batch_to_indices

        else:
            raise NotImplementedError

        # cut arr2 into batches
        max_arr2_batch_size = int(max_outward_size * matching_ratio)
        n_batches2 = max(1, int(self.shared_arr2.shape[0] // max_arr2_batch_size))
        arr2_batch_size = int(self.shared_arr2.shape[0] // n_batches2)
        self._batch_to_indices2 = _split(
            arr=self.active_arr2,
            n_batches=n_batches2,
            curr_seed=seed
        )

        # cut arr1 into batches
        max_arr1_batch_size = int(arr2_batch_size // (matching_ratio / metacell_size))
        n_batches1 = max(1, int(self.shared_arr1.shape[0] // max_arr1_batch_size))
        arr1_batch_size = int(self.shared_arr1.shape[0] // n_batches1)
        if seed is not None:
            seed = seed + 1
        self._batch_to_indices1 = _split(
            arr=self.active_arr1,
            n_batches=n_batches1,
            curr_seed=seed
        )

        # construct mapping between batches
        b1 = 0
        b2 = 0
        self._batch1_to_batch2 = []
        if batching_scheme == 'cyclic':
            for i in range(max(n_batches1, n_batches2)):
                self._batch1_to_batch2.append((b1, b2))
                b1 = (b1 + 1) % n_batches1
                b2 = (b2 + 1) % n_batches2
        elif batching_scheme == 'pairwise':
            for b1 in range(n_batches1):
                for b2 in range(n_batches2):
                    self._batch1_to_batch2.append((b1, b2))

        if verbose:
            print(('The first data is split into {} batches, '
                   'average batch size is {}, and max batch size is {}.').format(
                n_batches1, arr1_batch_size, len(self._batch_to_indices1[-1])
            ), flush=True)
            print(('The second data is split into {} batches, '
                   'average batch size is {}, and max batch size is {}.').format(
                n_batches2, arr2_batch_size, len(self._batch_to_indices2[-1])
            ), flush=True)
            print('Batch to batch correspondence is:\n  {}.'.format(
                [str(i) + '<->' + str(j) for i, j in self._batch1_to_batch2]
            ), flush=True)

        if arr1_batch_size <= 1000:
            warnings.warn('Batch size for arr1 is <= 1000: '
                          'consider setting a smaller matching ratio for better performance.')

    def plot_singular_values(
            self, target='shared_arr1', batch=None,
            n_components=None, randomized_svd=False, svd_runs=1
    ):
        """
        Plot the singular values of the target array for selection of SVD components.

        Parameters
        ----------
        target: string or np.array, default='shared_arr1'
            The array to fit SVD on.
            If is a string, then must be one of {'shared_arr1', 'shared_arr2', 'active_arr1', 'active_arr2'};
            otherwise, must be an np.array.
        batch: None or int, default=None
            If None, then randomly select on batch.
            If target array is an np.array, then this argument is disregarded.
        n_components: None or int, default=None
            How many components to plot.
            If None, set it to be either one or min(100, min(target.shape)-1), whichever is larger.
        randomized_svd: bool, default=False
            Whether to use randomized SVD.
        svd_runs: int, default=1,
            How many randomized SVDs will be executed.
            In the end the one with the lowest Frobenious reconstruction error will be taken.
            If randomized=False, this argument is disregarded.

        Returns
        -------
        fig, ax
        """
        if isinstance(target, str):
            assert target in {'shared_arr1', 'shared_arr2', 'active_arr1', 'active_arr2'}
            if batch is None:
                batch = np.random.randint(len(self._batch_to_indices1)) if target[-1] == '1' \
                    else np.random.randint(len(self._batch_to_indices2))
            if target == 'shared_arr1':
                target = self.shared_arr1[self._batch_to_indices1[batch], :]
            elif target == 'shared_arr2':
                target = self.shared_arr2[self._batch_to_indices2[batch], :]
            elif target == 'active_arr1':
                target = self.active_arr1[self._batch_to_indices1[batch], :]
            else:
                target = self.active_arr2[self._batch_to_indices2[batch], :]
        else:
            assert isinstance(target, np.ndarray) and len(target.shape) == 2

        if n_components is None:
            n_components = max(1, min(100, min(target.shape) - 1))

        _, s, _ = utils.robust_svd(arr=target, n_components=n_components, randomized=randomized_svd, n_runs=svd_runs)
        if not randomized_svd:
            s = s[::-1]  # svds orders the singular vectors in decreasing order
        fig, ax = plt.subplots()
        ax.plot(list(range(len(s))), s, '--o')
        ax.set(xlabel='Index', ylabel='Singular value')
        if batch is not None:
            ax.set(title='Singular value v.s. component index for batch {}'.format(batch))
        else:
            ax.set(title='Singular value v.s. component index')
        ax.xaxis.set_major_locator(MaxNLocator(integer=True))
        ax.grid()
        return fig, ax

    def _construct_graphs(
            self, target=1, n_neighbors=15, metacell=False, svd_components=None,
            randomized_svd=False, svd_runs=1, verbose=True, metric='correlation'
    ):
        """Construct neighborhood graphs for all batches of the target array.
        """
        if verbose:
            print('Constructing neighborhood graphs for cells in arr{}...'.format(target), flush=True)

        if target == 1:
            arr = self.active_arr1
            batch_to_indices = self._batch_to_indices1
            self._edges1 = []
            edges = self._edges1
        else:
            arr = self.active_arr2
            batch_to_indices = self._batch_to_indices2
            self._edges2 = []
            edges = self._edges2

        for b, indices in enumerate(batch_to_indices):
            if verbose:
                print('Now at batch {}...'.format(b), flush=True)
            if target == 1 and metacell:
                curr_arr = utils.get_centroids(arr=arr[indices, :], labels=self._metacell_labels1[b])
            else:
                curr_arr = arr[indices, :]

            edges.append(
                graph.construct_graph(
                    arr=curr_arr,
                    randomized_svd=randomized_svd,
                    svd_runs=svd_runs,
                    svd_components=svd_components,
                    n_neighbors=n_neighbors,
                    verbose=False,
                    metric=metric
                )
            )

        if verbose:
            print('Graph construction finished!', flush=True)

    def _cluster_graphs(
            self, target=1, metacell=False, resolution=1, leiden_runs=1, seed=None, verbose=True
    ):
        """Cluster the neighborhood graphs.
        """
        if verbose:
            print('Clustering the graphs for cells in arr{}...'.format(target), flush=True)

        if target == 1:
            batch_to_indices = self._batch_to_indices1
            edges = self._edges1
            self._labels1 = []
            labels = self._labels1
        else:
            batch_to_indices = self._batch_to_indices2
            edges = self._edges2
            self._labels2 = []
            labels = self._labels2

        for b, curr_edges in enumerate(edges):
            if verbose:
                print('Now at batch {}...'.format(b), flush=True)
            if target == 1 and metacell:
                n = len(np.unique(self._metacell_labels1[b]))
            else:
                n = len(batch_to_indices[b])
            labels.append(
                graph.graph_clustering(
                    n=n,
                    edges=curr_edges,
                    resolution=resolution,
                    n_runs=leiden_runs,
                    seed=seed,
                    verbose=False
                )
            )

        if verbose:
            print('Graph clustering finished!', flush=True)

    def construct_graphs(
            self, n_neighbors1=15, n_neighbors2=15,
            svd_components1=None, svd_components2=None,
            resolution1=1, resolution2=1,
            randomized_svd=False, svd_runs=1,
            resolution_tol=0.1,
            leiden_runs=1,
            metric='correlation',
            leiden_seed=None,
            verbose=True
    ):
        """
        Construct neighborhood graphs and cluster them as needed.

        Parameters
        ----------
        n_neighbors1: int, default=15
            Number of neighbors for graph construction for arr1.
        n_neighbors2: int, default=15
            Number of neighbors for graph construction for arr2.
        svd_components1: None or int, default=None
            If not None, perform SVD to reduce the dimension of self.active_arr1 before doing neighborhood search.
        svd_components2
            If not None, perform SVD to reduce the dimension of self.active_arr2 before doing neighborhood search.
        resolution1: int, default=1
            Resolution parameter for Leiden algorithm when clustering the graphs for arr1.
        resolution2: int, default=1
            Resolution parameter for Leiden algorithm when clustering the graphs for arr2.
        randomized_svd: bool, default=False
            Whether to perform randomized SVD.
        svd_runs: int, default=1
            Perform multiple runs of SVD and select the one with lowest Frobenious reconstruction error is selected.
        resolution_tol: float, default=0.1
            Any resolution within the range of plus/minus resolution_tol will not be differentiated.
        leiden_runs: int, default=1
            Perform multiple runs of Leiden algorithm and the one with highest modularity is selected.
        metric: string, default='correlation'
            The metric to use in nearest neighbor search.
        leiden_seed: None or int, default=None
            Random seed for Leiden algorithm. If leiden_runs>1, leiden_seed will be reset to None.
        verbose: bool, default=True
            Whether to print the progress.

        Returns
        -------
        None
        """
        if self.metacell_size > 1:
            if verbose:
                print(
                    'Aggregating cells in arr1 into metacells of average size {}...'.format(self.metacell_size),
                    flush=True
                )
            self._construct_graphs(
                target=1, n_neighbors=n_neighbors1, metacell=False, svd_components=svd_components1,
                randomized_svd=randomized_svd, svd_runs=svd_runs, verbose=verbose, metric=metric
            )
            if verbose:
                print('Clustering into metacells...', flush=True)
            self._metacell_labels1 = []
            for b, curr_edges in enumerate(self._edges1):
                if verbose:
                    print('Now at batch {}...'.format(b), flush=True)
                n = len(self._batch_to_indices1[b])
                self._metacell_labels1.append(
                    graph.graph_clustering(
                        n=n,
                        edges=curr_edges,
                        n_clusters=int(n // self.metacell_size),
                        resolution=None,
                        n_runs=leiden_runs,
                        resolution_tol=resolution_tol,
                        seed=leiden_seed,
                        verbose=False
                    )
                )

                if verbose:
                    print('Metacell clustering finished!', flush=True)

        if self.method == 'graph_smoothing':
            # only need to construct graphs, no need to do clustering
            # note that if we use metacell, then self._edges1 is not None and is filled with the graphs on single cells
            # we will overwrite it with graphs on metacells
            self._construct_graphs(
                target=1, n_neighbors=n_neighbors1, metacell=(self.metacell_size > 1),
                svd_components=svd_components1,
                randomized_svd=randomized_svd, svd_runs=svd_runs, verbose=verbose, metric=metric
            )
            self._construct_graphs(
                target=2, n_neighbors=n_neighbors2, metacell=False, svd_components=svd_components2,
                randomized_svd=randomized_svd, svd_runs=svd_runs, verbose=verbose, metric=metric
            )
        elif self.method == 'centroid_shrinkage':
            # for arr1
            if self.labels1 is not None:
                # separate self.labels1 into batches
                self._labels1 = []
                for b, indices in enumerate(self._batch_to_indices1):
                    if self.metacell_size > 1:
                        # majority voting to determine the metacell cluster labels
                        curr_labels = utils.summarize_clustering(
                            clustering=self._metacell_labels1[b], true_labels=self.labels1[indices]
                        )
                    else:
                        curr_labels = self.labels1[indices]
                    self._labels1.append(utils.recode(curr_labels)[0])
            else:
                # need to construct the graph + cluster
                self._construct_graphs(
                    target=1, n_neighbors=n_neighbors1, metacell=(self.metacell_size > 1),
                    svd_components=svd_components1,
                    randomized_svd=randomized_svd, svd_runs=svd_runs, verbose=verbose, metric=metric
                )
                self._cluster_graphs(
                    target=1, metacell=(self.metacell_size > 1),
                    resolution=resolution1, leiden_runs=leiden_runs, seed=leiden_seed, verbose=verbose
                )
            # for arr2
            if self.labels2 is not None:
                # separate self.labels2 into batches
                self._labels2 = []
                for indices in self._batch_to_indices2:
                    self._labels2.append(utils.recode(self.labels2[indices])[0])
            else:
                # need to construct the graph + cluster
                self._construct_graphs(
                    target=2, n_neighbors=n_neighbors2, metacell=False, svd_components=svd_components2,
                    randomized_svd=randomized_svd, svd_runs=svd_runs, verbose=verbose, metric=metric
                )
                self._cluster_graphs(
                    target=2, metacell=False,
                    resolution=resolution2, leiden_runs=leiden_runs, seed=leiden_seed, verbose=verbose
                )
        else:
            raise ValueError('self.method must be one of \'centroid_shrinkage\' or \'graph_smoothing\'.')

    def find_initial_pivots(
            self,
            wt1=0.3, wt2=0.3,
            svd_components1=None, svd_components2=None,
            randomized_svd=False, svd_runs=1,
            verbose=True
    ):
        """
        Perform initial matching.

        Parameters
        ----------
        wt1: float, default=0.3
            The shrinkage weight to put on the raw data for arr1.
        wt2: float, default=0.3
            The shrinkage weight to put on the raw data for arr2.
        svd_components1: None or int, default=None
            If not None, perform SVD to reduce the dimension of self.shared_arr1.
        svd_components2: None or int, default=None
            If not None, perform SVD to reduce the dimension of self.shared_arr2.
        randomized_svd: bool, default=False
            Whether to use randomized SVD.
        svd_runs: int, default=1
            Perform multiple runs of SVD and the one with lowest Frobenious reconstruction error is selected.
        verbose: bool, default=True
            Whether to print the progress.

        Returns
        -------
        None
        """
        self._init_matching = []
        for b1, b2 in self._batch1_to_batch2:
            if verbose:
                print(
                    'Now at batch {}<->{}...'.format(b1, b2),
                    flush=True
                )
            if self.metacell_size > 1:
                arr1 = utils.get_centroids(
                    arr=self.shared_arr1[self._batch_to_indices1[b1], :],
                    labels=self._metacell_labels1[b1]
                )
            else:
                arr1 = self.shared_arr1[self._batch_to_indices1[b1], :]

            arr2 = self.shared_arr2[self._batch_to_indices2[b2], :]

            edges1, edges2, clust_labels1, clust_labels2 = None, None, None, None
            if self.method == 'centroid_shrinkage':
                clust_labels1 = self._labels1[b1]
                clust_labels2 = self._labels2[b2]
            elif self.method == 'graph_smoothing':
                edges1 = self._edges1[b1]
                edges2 = self._edges2[b2]
            else:
                raise ValueError('self.method must be one of \'centroid_shrinkage\' or \'graph_smoothing\'.')

            self._init_matching.append(
                match_utils.get_initial_matching(
                    arr1=arr1,
                    arr2=arr2,
                    clust_labels1=clust_labels1,
                    clust_labels2=clust_labels2,
                    edges1=edges1,
                    edges2=edges2,
                    wt1=wt1,
                    wt2=wt2,
                    randomized_svd=randomized_svd,
                    svd_runs=svd_runs,
                    svd_components1=svd_components1,
                    svd_components2=svd_components2,
                    verbose=False
                )
            )

        if verbose:
            print('Done!', flush=True)

    def plot_canonical_correlations(
            self, batch=None,
            svd_components1=None, svd_components2=None,
            cca_components=None,
            filter_prop=0.,
            randomized_svd=False,
            svd_runs=1,
            cca_max_iter=2000
    ):
        """
        Perform CCA on active arrays aligned by initial matching, and plot the canonical correlations.

        Parameters
        ----------
        batch: None or int or tuple of two integers
            The arrays to perform CCA on.
            If None, randomly select a batch index,
            if an integer, then the batches are self._batch1_to_batch2[batch],
            if a tuple of two integers, and batch[0] and batch[1] are selected.
        svd_components1: None or int, default=None
            If not None, perform SVD to reduce the dimension of self.active_arr1 before feeding it to CCA.
        svd_components2: None or int, default=None
            If not None, perform SVD to reduce the dimension of self.active_arr2 before feeding it to CCA.
        cca_components: None or int, default=None
            Number of CCA components.
            If None, it is set to 100 or self.active_arr1.shape[1] or self.active_arr2.shape[1], whichever is smaller.
        filter_prop: float, default=0.
            CCA is performed on top 1-filter_prop slice of the data on which the matched distances are smallest.
        randomized_svd: bool, default=False
            Whether to perform randomized SVD.
        svd_runs: int, default=1
            Perform multiple runs of SVD and the one with lowest Frobenious reconstruction error is selected.
        cca_max_iter: int, default=2000
            Maximum iteration number for CCA.

        Returns
        -------
        fig, ax
        """
        if batch is None:
            idx = np.random.randint(len(self._batch1_to_batch2))
            b1, b2 = self._batch1_to_batch2[idx]
        elif isinstance(batch, int):
            idx = batch
            b1, b2 = self._batch1_to_batch2[idx]
        elif isinstance(batch, Iterable):
            idx, b1, b2 = None, None, None
            for curr_idx, (curr_b1, curr_b2) in enumerate(self._batch1_to_batch2):
                input_b1, input_b2 = batch
                if curr_b1 == input_b1 and curr_b2 == input_b2:
                    idx = curr_idx
                    b1 = curr_b1
                    b2 = curr_b2
                    break
            if idx is None:
                input_b1, input_b2 = batch
                raise ValueError('Batch {} in arr1 is not matched to batch {} in arr2'.format(input_b1, input_b2))
        else:
            raise ValueError('batch must be an integer, a tuple of two integers, or None.')

        if self.metacell_size > 1:
            arr1 = utils.get_centroids(
                arr=self.active_arr1[self._batch_to_indices1[b1], :],
                labels=self._metacell_labels1[b1]
            )
        else:
            arr1 = self.active_arr1[self._batch_to_indices1[b1], :]

        arr2 = self.active_arr2[self._batch_to_indices2[b2], :]
        # do SVD if needed
        if svd_components1 is not None:
            arr1 = utils.svd_embedding(arr=arr1, n_components=svd_components1, randomized=randomized_svd,
                                       n_runs=svd_runs)
        if svd_components2 is not None:
            arr2 = utils.svd_embedding(arr=arr2, n_components=svd_components2, randomized=randomized_svd,
                                       n_runs=svd_runs)

        if cca_components is None:
            cca_components = min([
                100,
                arr1.shape[1],
                arr2.shape[1],
                int(arr1.shape[0]//2),
                int(arr2.shape[0]//2)
            ])

        # fit CCA and extract canonical correlations
        cc = utils.cca_embedding(
            arr1=arr1, arr2=arr2,
            init_matching=self._init_matching[idx],
            filter_prop=filter_prop, n_components=cca_components, max_iter=cca_max_iter
        )[2]

        fig, ax = plt.subplots()
        ax.plot(list(range(len(cc))), cc, '--o')
        ax.set(xlabel='Index', ylabel='Canonical correlation')
        ax.set(title='Canonical correlation v.s. component index for batch {}<->{}'.format(b1, b2))
        ax.xaxis.set_major_locator(MaxNLocator(integer=True))
        ax.grid()
        return fig, ax

    def refine_pivots(
            self,
            wt1=0.5, wt2=0.5,
            svd_components1=None, svd_components2=None,
            cca_components=None,
            filter_prop=0,
            n_iters=1,
            randomized_svd=False, svd_runs=1,
            cca_max_iter=2000,
            verbose=True
    ):
        """
        Perform refined matching.

        Parameters
        ----------
        wt1: float, default=0.3
            The shrinkage weight to put on the raw data for arr1.
        wt2: float, default=0.3
            The shrinkage weight to put on the raw data for arr2.
        svd_components1: None or int, default=None
            If not None, perform SVD to reduce the dimension of self.active_arr1 before feeding it to CCA.
        svd_components2: None or int, default=None
            If not None, perform SVD to reduce the dimension of self.active_arr2 before feeding it to CCA.
        cca_components: None or int, default=None
            Number of CCA components.
            If None, it is set to 100 or self.active_arr1.shape[1] or self.active_arr2.shape[1], whichever is smaller.
        filter_prop: float, default=0.
            CCA is performed on top 1-filter_prop slice of the data on which the matched distances are smallest.
        n_iters: int, default=1
            Number of refinement iterations.
        randomized_svd: bool, default=False
            Whether to perform randomized SVD.
        svd_runs: int, default=1
            Perform multiple runs of SVD and the one with lowest Frobenious reconstruction error is selected.
        cca_max_iter: int, default=2000
            Maximum iteration number for CCA.
        verbose: bool, default=True
            Whether to print the progress.

        Returns
        -------
        None
        """
        # save cca parameters for later use
        self._svd_components1_for_cca_embedding = svd_components1
        self._svd_components2_for_cca_embedding = svd_components2
        self._randomized_svd_for_cca_embedding = randomized_svd
        self._svd_runs_for_cca_embedding = svd_runs
        self._cca_components = cca_components
        self._cca_max_iter = cca_max_iter

        self._refined_matching = []
        for batch_idx, (b1, b2) in enumerate(self._batch1_to_batch2):
            if verbose:
                print(
                    'Now at batch {}<->{}...'.format(b1, b2),
                    flush=True
                )
            arr1_init, arr2_init = None, None
            if self.metacell_size > 1:
                arr1 = utils.get_centroids(
                    arr=self.active_arr1[self._batch_to_indices1[b1], :],
                    labels=self._metacell_labels1[b1]
                )
            else:
                arr1 = self.active_arr1[self._batch_to_indices1[b1], :]

            arr2 = self.active_arr2[self._batch_to_indices2[b2], :]
            arr2_init = self.shared_arr2[self._batch_to_indices2[b2], :]

            edges1, edges2, clust_labels1, clust_labels2 = None, None, None, None
            if self.method == 'centroid_shrinkage':
                clust_labels1 = self._labels1[b1]
                clust_labels2 = self._labels2[b2]
            elif self.method == 'graph_smoothing':
                edges1 = self._edges1[b1]
                edges2 = self._edges2[b2]
            else:
                raise ValueError('self.method must be one of \'centroid_shrinkage\' or \'graph_smoothing\'.')

            self._refined_matching.append(
                match_utils.get_refined_matching(
                    init_matching=self._init_matching[batch_idx],
                    arr1=arr1,
                    arr2=arr2,
                    randomized_svd=randomized_svd,
                    svd_runs=svd_runs,
                    svd_components1=svd_components1,
                    svd_components2=svd_components2,
                    clust_labels1=clust_labels1,
                    clust_labels2=clust_labels2,
                    edges1=edges1,
                    edges2=edges2,
                    wt1=wt1,
                    wt2=wt2,
                    n_iters=n_iters,
                    filter_prop=filter_prop,
                    cca_components=cca_components,
                    cca_max_iter=cca_max_iter,
                    verbose=False
                )
            )

        if verbose:
            print('Done!', flush=True)

#     def plot_matching_scores(
#             self, target='pivot',
#             batch=None,
#             detect_changepoint=True,
#             min_score=0.,
#             max_score=0.8,
#             n_bkps=5
#     ):
#         """
#         Plot the similarity scores of matched pairs.
#
#         Parameters
#         ----------
#         target: 'pivot' or 'propagated'
#             Either to plot scores on refined matching or the propagated matching.
#         batch: None or int or tuple of two integers
#             The batch to extract the scores.
#             If None, randomly select a batch index,
#             if an integer, then the batches are self._batch1_to_batch2[batch],
#             if a tuple of two integers, and batch[0] and batch[1] are selected.
#         detect_changepoint: bool, default=True
#             Whether to perform changepoint detect to suggest the cutoff.
#         min_score: float, default=0.
#             Anything <= min_score will definitely be filtered.
#         max_score: float, default=0.8
#             Anything >= max_score will definitely be filtered.
#         n_bkps: int, default=5
#             Number of break points to detect.
#             Changepoint detection will be performed on the consecutive differences of the sorted scores whose values
#             are bounded in (min_score, max_score).
#
#         Returns
#         -------
#         fig, ax if detect_changepoint=False,
#         fig, ax, break_point_locations if detect_changepoint=True.
#         """
#         # for choosing the filtering threshold
#         # anything >= max_score will definitely be retained
#
#         if batch is None:
#             idx = np.random.randint(len(self._batch1_to_batch2))
#             b1, b2 = self._batch1_to_batch2[idx]
#         elif isinstance(batch, int):
#             idx = batch
#             b1, b2 = self._batch1_to_batch2[idx]
#         elif isinstance(batch, Iterable):
#             idx, b1, b2 = None, None, None
#             for curr_idx, (curr_b1, curr_b2) in enumerate(self._batch1_to_batch2):
#                 input_b1, input_b2 = batch
#                 if curr_b1 == input_b1 and curr_b2 == input_b2:
#                     idx = curr_idx
#                     b1 = curr_b1
#                     b2 = curr_b2
#                     break
#             if idx is None:
#                 input_b1, input_b2 = batch
#                 raise ValueError('Batch {} in arr1 is not matched to batch {} in arr2'.format(input_b1, input_b2))
#         else:
#             raise ValueError('batch must be an integer, a tuple of two integers, or None.')
#
#         if target == 'pivot':
#             increasing_scores = 1 - self._refined_matching[idx][2]
#         elif target == 'propagated':
#             increasing_scores = 1 - self._propagated_matching[idx][2]
#         else:
#             raise ValueError('target must be in {\'pivot\', \'propagated\'}')
#
#         increasing_scores = np.array(sorted(increasing_scores))
#         if not detect_changepoint:
#             fig, ax = plt.subplots()
#             ax.plot(np.arange(len(increasing_scores)) / len(increasing_scores), increasing_scores[::-1])
#             ax.set(xlabel='Quantile', ylabel='Matching score')
#             ax.set(title='Matching scores in decreasing order for batch {}<->{}'.format(b1, b2))
#             ax.grid()
#             return fig, ax
#
#         # do changepoint detection
#         assert min_score < max_score
#         signal = increasing_scores[1:] - increasing_scores[:-1]
#         # first time the score is > min_score
#         start = bisect.bisect_right(increasing_scores, min_score)
#         if start == len(increasing_scores):
#             raise ValueError('All the scores are <= min_score, please set min_score to be smaller.')
#         # first time score is >= max_score
#         end = bisect.bisect_left(increasing_scores, max_score)
#         if end == 0:
#             raise ValueError('All the scores are >= max_score, please set max_score to be larger.')
#
#         algo = rpt.Binseg(model='l2').fit(signal[start:end])
#         bkps = np.array(algo.predict(n_bkps=n_bkps)) + start
#
#         # convert locations of min_scores, max_score, and bkps to indices in increasing_scores[::-1]
#         start, end = len(increasing_scores) - end - 1, len(increasing_scores) - start - 1
#         bkps = len(increasing_scores) - bkps - 1
#         start, end, bkps = start/len(increasing_scores), end/len(increasing_scores), bkps/len(increasing_scores)
#         # plot
#         fig, ax = plt.subplots()
#         ax.plot(np.arange(len(increasing_scores)) / len(increasing_scores), increasing_scores[::-1])
#         ax.axvline(start, ls='-.', color='red')
#         ax.text(start, max_score, 'max_score={}'.format(max_score))
#         for bkp in bkps:
#             ax.axvline(bkp, ls='--')
#         ax.axvline(end, ls='-.', color='red')
#         ax.text(end, max_score, 'min_score={}'.format(min_score))
#         ax.set(xlabel='Quantile', ylabel='Matching score')
#         ax.set(title='Matching scores in decreasing order for batch {}<->{}'.format(b1, b2))
#         ax.grid()
#         return fig, ax, bkps[::-1]

    def _fit_svd_on_full_data(self):
        """Perform SVD on full self.active_arr1 and self.active_arr2 and save the functions that reduce the dimension
        of the data.
        """
        if self._svd_components1_for_cca_embedding is not None:
            u1, s1, vh1 = utils.robust_svd(
                arr=self.active_arr1, n_components=self._svd_components1_for_cca_embedding,
                randomized=self._randomized_svd_for_cca_embedding, n_runs=self._svd_runs_for_cca_embedding
            )
            self._rotation_before_cca_on_active_arr1 = lambda arr: arr @ vh1.T
        else:
            self._rotation_before_cca_on_active_arr1 = lambda arr: arr

        if self._svd_components2_for_cca_embedding is not None:
            u2, s2, vh2 = utils.robust_svd(
                arr=self.active_arr2, n_components=self._svd_components2_for_cca_embedding,
                randomized=self._randomized_svd_for_cca_embedding, n_runs=self._svd_runs_for_cca_embedding
            )
            self._rotation_before_cca_on_active_arr2 = lambda arr: arr @ vh2.T
        else:
            self._rotation_before_cca_on_active_arr2 = lambda arr: arr

    def filter_bad_matches(self, target='pivot', filter_prop=0., verbose=True):
        """
        Filter matched pairs of low quality and score each remaining pairs.

        Parameters
        ----------
        target: 'pivot' or 'propagated'
            Which matching to perform filtering on.
        filter_prop: float, default=0.
            The proportion of matching to be filtered out.
        verbose: bool, default=True
            Whether to print the progress.

        Returns
        -------
        None
        """
        # construct mapping from metacell index to single cell indices
        # useful when stitching the batch-wise matchings together
        metacell_idx_to_array_indices = []
        if self.metacell_size > 1:
            for metacell_labels in self._metacell_labels1:
                curr_dict = defaultdict(list)
                for single_idx, meta_idx in enumerate(metacell_labels):
                    curr_dict[meta_idx].append(single_idx)
                metacell_idx_to_array_indices.append(curr_dict)

        # conduct filtering
        n_remaining = 0
        if verbose:
            print('Begin filtering...', flush=True)
        if target == 'pivot':
            matching_to_be_filtered = self._refined_matching
        elif target == 'propagated':
            matching_to_be_filtered = self._propagated_matching
        else:
            raise ValueError('target must be in {\'pivot\', \'propagated\'}.')

        batch_to_remaining_indices_after_filtering = []
        idx2_to_indices1 = defaultdict(set)
        # record the locations that survive the filtering
        for batch_idx, (b1, b2) in enumerate(self._batch1_to_batch2):
            if verbose:
                print('Now at batch {}<->{}...'.format(b1, b2), flush=True)
            rows, cols, vals = matching_to_be_filtered[batch_idx]
            # anything with val <= thresh will be retained
            thresh = np.quantile(vals, 1-filter_prop)
            batch_to_remaining_indices_after_filtering.append([i for i in range(len(vals)) if vals[i] <= thresh])
            n_remaining += len(batch_to_remaining_indices_after_filtering[batch_idx])
            for i in batch_to_remaining_indices_after_filtering[batch_idx]:
                r, c = rows[i], cols[i]
                idx2 = self._batch_to_indices2[b2][c]
                if self.metacell_size > 1:
                    indices1 = np.array(self._batch_to_indices1[b1])[metacell_idx_to_array_indices[b1][r]]
                else:
                    indices1 = [self._batch_to_indices1[b1][r]]
                for idx1 in indices1:
                    idx2_to_indices1[idx2].add(idx1)

        if target == 'pivot':
            self._remaining_indices_in_refined_matching = batch_to_remaining_indices_after_filtering
        else:
            self._remaining_indices_in_propagated_matching = batch_to_remaining_indices_after_filtering

        if verbose:
            print('{}/{} pairs of matched cells remain after the filtering.'.format(
                n_remaining, np.sum([len(per_batch_matching[0]) for per_batch_matching in matching_to_be_filtered])
            ))

        # convert idx2_to_indices1 to a dict of lists for ease of later usage
        idx2_to_indices1 = {idx2: sorted(indices1) for idx2, indices1 in idx2_to_indices1.items()}

        # fit CCA on pivots
        if target == 'pivot':
            if verbose:
                print('Fitting CCA on pivots...', flush=True)
            self._fit_svd_on_full_data()
            arr1, arr2 = [], []
            for idx2, indices1 in idx2_to_indices1.items():
                arr1.append(np.mean(self.active_arr1[indices1, :], axis=0))
                arr2.append(self.active_arr2[idx2, :])
            arr1 = self._rotation_before_cca_on_active_arr1(np.array(arr1))
            arr2 = self._rotation_before_cca_on_active_arr2(np.array(arr2))
            self._cca_on_pivots = CCA(n_components=self._cca_components, max_iter=self._cca_max_iter)
            self._cca_on_pivots.fit(arr1, arr2)

        if verbose:
            print('Scoring matched pairs...', flush=True)
        # transform the whole dataset
        arr1 = self._rotation_before_cca_on_active_arr1(self.active_arr1)
        arr2 = self._rotation_before_cca_on_active_arr2(self.active_arr2)
        arr1, arr2 = self._cca_on_pivots.transform(arr1, arr2)
        arr1 = utils.center_scale(arr1)
        arr2 = utils.center_scale(arr2)

        # add distances to idx2_to_indices1
        all_indices1, all_indices2 = [], []
        for idx2, indices1 in idx2_to_indices1.items():
            for idx1 in indices1:
                all_indices1.append(idx1)
                all_indices2.append(idx2)
        # compute all possible pairs of distances
        pearson_correlations = utils.pearson_correlation(arr1[all_indices1, :], arr2[all_indices2])

        cnt = 0
        for idx2, indices1 in idx2_to_indices1.items():
            indices_and_scores = []
            for idx1 in indices1:
                indices_and_scores.append((idx1, pearson_correlations[cnt]))
                cnt += 1
            idx2_to_indices1[idx2] = indices_and_scores

        if target == 'pivot':
            self._pivot2_to_pivots1 = utils.sort_dict(idx2_to_indices1)
            unique_pivots1 = set()
            for idx2, indices1_and_scores in self._pivot2_to_pivots1.items():
                for idx1, _ in indices1_and_scores:
                    unique_pivots1.add(idx1)
            print('{}/{} cells in arr1 are selected as pivots.'.format(
                len(unique_pivots1), self.active_arr1.shape[0]
            ), flush=True)
            print('{}/{} cells in arr2 are selected as pivots.'.format(
                len(self._pivot2_to_pivots1), self.active_arr2.shape[0]
            ), flush=True)
        else:
            self._propidx2_to_propindices1 = utils.sort_dict(idx2_to_indices1)

        if verbose:
            print('Done!', flush=True)

    def propagate(
            self, wt1=0.7, wt2=0.7, svd_components1=None, svd_components2=None,
            metric='euclidean',
            randomized_svd=False, svd_runs=1, verbose=True
    ):
        """
        For indices not in pivots, find their matches by nearest neighbor search.

        Parameters
        ----------
        svd_components1: None or int, default=None
            If not None, perform SVD to reduce the dimension of self.active_arr1
            before doing internal nearest neighbor search.
        svd_components2: None or int, default=None
            If not None, perform SVD to reduce the dimension of self.active_arr1
            before doing internal nearest neighbor search.
        wt1: float, default=0.7
            Weight to put on raw data of self.active_arr1 when doing smoothing.
        wt2: float, default=0.7
            Weight to put on raw data of self.active_arr2 when doing smoothing.
        metric: string, default='correlation'
            The metric to use in nearest neighbor search.
        randomized_svd: bool, default=False
            Whether to perform randomized SVD.
        svd_runs: int, default=1
            Perform multiple runs of SVD and the one with lowest Frobenious reconstruction error is selected.
        verbose: bool, default=True
            Whether to print the progress.

        Returns
        -------
        None
        """
        self._propagated_matching = []
        for batch_idx, (b1, b2) in enumerate(self._batch1_to_batch2):
            if verbose:
                print('Now at batch {}<->{}...'.format(b1, b2), flush=True)

            curr_propagated_matching = [[], [], []]
            curr_refined_matching = self._refined_matching[batch_idx]
            # get good pivot indices that survived pivot filtering
            existing_indices = curr_refined_matching[0][self._remaining_indices_in_refined_matching[batch_idx]]
            good_indices1 = curr_refined_matching[0][existing_indices]
            good_indices2 = curr_refined_matching[1][existing_indices]

            # get arrays that were used when doing refined matching
            if self.metacell_size > 1:
                curr_arr1 = utils.get_centroids(
                    arr=self.active_arr1[self._batch_to_indices1[b1], :], labels=self._metacell_labels1[b1]
                )
            else:
                curr_arr1 = self.active_arr1[self._batch_to_indices1[b1], :]
            curr_arr2 = self.active_arr2[self._batch_to_indices2[b2], :]

            # do smoothing
            if self.method == 'centroid_shrinkage':
                clust_labels1 = self._labels1[b1]
                clust_labels2 = self._labels2[b2]
                curr_arr1 = utils.shrink_towards_centroids(arr=curr_arr1, labels=clust_labels1, wt=wt1)
                curr_arr2 = utils.shrink_towards_centroids(arr=curr_arr2, labels=clust_labels2, wt=wt2)
            elif self.method == 'graph_smoothing':
                edges1 = self._edges1[b1]
                edges2 = self._edges2[b2]
                curr_arr1 = utils.graph_smoothing(arr=curr_arr1, edges=edges1, wt=wt1)
                curr_arr2 = utils.graph_smoothing(arr=curr_arr2, edges=edges2, wt=wt1)
            else:
                raise ValueError('self.method must be one of \'centroid_shrinkage\' or \'graph_smoothing\'.')

            # get remaining indices
            # propagation will only be done for those indices
            good_indices1_set = set(good_indices1)
            remaining_indices1 = [i for i in range(curr_arr1.shape[0]) if i not in good_indices1_set]
            good_indices2_set = set(good_indices2)
            remaining_indices2 = [i for i in range(curr_arr2.shape[0]) if i not in good_indices2_set]

            # propagate for remaining indices in arr1
            if len(remaining_indices1) > 0:
                # get 1-nearest-neighbors and the corresponding distances
                remaining_indices1_nns, remaining_indices1_nn_dists = graph.get_nearest_neighbors(
                    query_arr=curr_arr1[remaining_indices1, :],
                    target_arr=curr_arr1[good_indices1, :],
                    svd_components=svd_components1,
                    randomized_svd=randomized_svd,
                    svd_runs=svd_runs,
                    metric=metric
                )
                matched_indices2 = good_indices2[remaining_indices1_nns]
                curr_propagated_matching[0].extend(remaining_indices1)
                curr_propagated_matching[1].extend(matched_indices2)
                curr_propagated_matching[2].extend(remaining_indices1_nn_dists)

            # propagate for remaining indices in arr2
            if len(remaining_indices2) > 0:
                # get 1-nearest-neighbors and the corresponding distances
                remaining_indices2_nns, remaining_indices2_nn_dists = graph.get_nearest_neighbors(
                    query_arr=curr_arr2[remaining_indices2, :],
                    target_arr=curr_arr2[good_indices2, :],
                    svd_components=svd_components2,
                    randomized_svd=randomized_svd,
                    svd_runs=svd_runs,
                    metric=metric
                )
                matched_indices1 = good_indices1[remaining_indices2_nns]
                curr_propagated_matching[0].extend(matched_indices1)
                curr_propagated_matching[1].extend(remaining_indices2)
                curr_propagated_matching[2].extend(remaining_indices2_nn_dists)

            self._propagated_matching.append(
                (np.array(curr_propagated_matching[0]),
                 np.array(curr_propagated_matching[1]),
                 np.array(curr_propagated_matching[2]))
            )

        if verbose:
            print('Done!', flush=True)

    def get_matching(self, order=None, target='pivot'):
        """
        Return a copy of the desired matching.

        Parameters
        ----------
        order: None or (1, 2) or (1, 2), default=None
            If (1, 2), then every cell in target arr1 has at least one match
            if (2, 1), then does the other way around,
            if None, then every cell in target arr1 and every cell in target arr2 both have at least one match
        target: 'pivot' or 'full_data'
            If 'pivot', then only return matching on pivots, else return matching on all the data.

        Returns
        -------
        A matching of format dict or list.
        """
        if target not in {'pivot', 'full_data'}:
            raise ValueError('mode must be in {\'pivot_only\', \'full_data\'}.')
        # if return_format not in {'dict', 'list'}:
        #     raise ValueError('return_format must be in {\'dict\', \'list\'}.')
        res = [[], [], []]
        for idx2, indices1_and_scores in self._pivot2_to_pivots1.items():
            for idx1, score in indices1_and_scores:
                res[0].append(idx1)
                res[1].append(idx2)
                res[2].append(score)
        if target == 'pivot':
            return res
        elif target == 'full_data':
            if order == (1, 2):
                # add propagated matching for non-pivot cells in the first dataset
                existing_indices1 = np.unique(res[0])
                remaining_indices1 = [i for i in range(self.active_arr1.shape[0]) if i not in existing_indices1]
                propagated_idx1_to_indices2 = defaultdict(list)
                for idx2, indices1_and_scores in self._propidx2_to_propindices1.items():
                    for idx1, score in indices1_and_scores:
                        propagated_idx1_to_indices2[idx1].append((idx2, score))
                for idx1 in remaining_indices1:
                    if idx1 in propagated_idx1_to_indices2:
                        for idx2, score in propagated_idx1_to_indices2[idx1]:
                            res[0].append(idx1)
                            res[1].append(idx2)
                            res[2].append(score)
            elif order == (2, 1):
                # add propagated matching for non-pivot cells in the second dataset
                existing_indices2 = np.unique(res[1])
                remaining_indices2 = [i for i in range(self.active_arr2.shape[0]) if i not in existing_indices2]
                for idx2 in remaining_indices2:
                    if idx2 in self._propidx2_to_propindices1:
                        for idx1, score in self._propidx2_to_propindices1[idx2]:
                            res[0].append(idx1)
                            res[1].append(idx2)
                            res[2].append(score)
            elif order is None:
                # first do order (1, 2) and then do order (2, 1)
                # add propagated matching for non-pivot cells in the first dataset
                existing_indices1 = np.unique(res[0])
                remaining_indices1 = [i for i in range(self.active_arr1.shape[0]) if i not in existing_indices1]
                propagated_idx1_to_indices2 = defaultdict(list)
                for idx2, indices1_and_scores in self._propidx2_to_propindices1.items():
                    for idx1, score in indices1_and_scores:
                        propagated_idx1_to_indices2[idx1].append((idx2, score))
                for idx1 in remaining_indices1:
                    if idx1 in propagated_idx1_to_indices2:
                        for idx2, score in propagated_idx1_to_indices2[idx1]:
                            res[0].append(idx1)
                            res[1].append(idx2)
                            res[2].append(score)

                # add propagated matching for non-pivot cells in the second dataset
                existing_indices2 = np.unique(res[1])
                remaining_indices2 = [i for i in range(self.active_arr2.shape[0]) if i not in existing_indices2]
                for idx2 in remaining_indices2:
                    if idx2 in self._propidx2_to_propindices1:
                        for idx1, score in self._propidx2_to_propindices1[idx2]:
                            res[0].append(idx1)
                            res[1].append(idx2)
                            res[2].append(score)

            else:
                raise NotImplementedError('order must be None or (1, 2) or (2, 1).')
        else:
            raise NotImplementedError('target must be in {\'pivot\', \'full_data\'}.')

        return match_utils.address_matching_redundancy(matching=res, order=order)

    def get_embedding(
            self,
            active_arr1,
            active_arr2,
            refit=False,
            matching=None,
            order=None,
            cca_components=None,
            cca_max_iter=None
    ):
        """
        Get CCA embedding.

        Parameters
        ----------
        active_arr1: np.ndarray of shape (n_samples_1, n_features_1)
            The first data matrix.
        active_arr2: np.ndarray of shape (n_samples_2, n_features_2)
            The second data matrix.
        refit: bool, default=False
            Whether to refit CCA.
        matching: None or list of length three, default=None
            Must be provided when refit=True,
            rows, cols, vals = matching,
            each matched pair of rows[i], cols[i], their score (the larger, the better) is vals[i].
        order: None or (1, 2) or (2, 1), default=None
            If None, then directly use matching to align the data and fit CCA;
            if (1, 2), then for every cell in arr1, average the corresponding matches in arr2;
            if (2, 1), do the other way around.
        cca_components: None or int, default=None
            Number of CCA components, if None, use self._cca_components.
        cca_max_iter: None or int, default=None
            Maximum number of CCA iteration, if None, use self._cca_max_iter.

        Returns
        -------
        arr1, arr2: np.arrays representing the CCA embedding
        """
        if refit:
            assert matching is not None
            # need to refit
            if cca_components is None:
                cca_components = self._cca_components
            if cca_max_iter is None:
                cca_max_iter = self._cca_max_iter

            arr1, arr2 = [], []
            if order is None:
                arr1, arr2 = self.active_arr1[matching[0], :], self.active_arr2[matching[1], :]
            elif order == (1, 2):
                idx1_to_indices2_and_scores = utils.list_to_dict(matching)
                for idx1, indices2_and_scores in idx1_to_indices2_and_scores.items():
                    indices2 = [idx2 for idx2, _ in indices2_and_scores]
                    arr1.append(self.active_arr1[idx1, :])
                    arr2.append(np.mean(self.active_arr2[indices2, :], axis=0))
            elif order == (2, 1):
                idx2_to_indices1_and_scores = utils.list_to_dict([matching[1], matching[0], matching[2]])
                for idx2, indices1_and_scores in idx2_to_indices1_and_scores.items():
                    indices1 = [idx1 for idx1, _ in indices1_and_scores]
                    arr1.append(np.mean(self.active_arr1[indices1, :], axis=0))
                    arr2.append(self.active_arr2[idx2, :])
            else:
                raise NotImplementedError('order must be None or (1, 2) or (2, 1).')

            arr1 = self._rotation_before_cca_on_active_arr1(np.array(arr1))
            arr2 = self._rotation_before_cca_on_active_arr2(np.array(arr2))
            cca = CCA(n_components=cca_components, max_iter=cca_max_iter)
            cca.fit(arr1, arr2)
        else:
            cca = self._cca_on_pivots

        arr1_cca, arr2_cca = cca.transform(
            self._rotation_before_cca_on_active_arr1(active_arr1),
            self._rotation_before_cca_on_active_arr2(active_arr2)
        )
        arr1_cca = utils.center_scale(arr1_cca)
        arr2_cca = utils.center_scale(arr2_cca)

        return arr1_cca, arr2_cca

