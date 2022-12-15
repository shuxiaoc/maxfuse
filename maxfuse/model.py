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
