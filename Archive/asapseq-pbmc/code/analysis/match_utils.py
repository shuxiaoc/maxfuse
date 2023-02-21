import numpy as np
from scipy.optimize import linear_sum_assignment
import utils


def address_matching_redundancy(matching, order=(1, 2)):
    """
    Make a potentially multiple-to-multiple matching to an one-to-one matching according to order.
    Parameters
    ----------
    matching: list of length three
        rows, cols, vals = matching: list
        Each matched pair of rows[i], cols[i], their score (the larger, the better) is vals[i]
    order: (1, 2) or (2, 1)
        If (1, 2), then the redundancy is addressed by making matching
        an injective map from the first dataset to the second;
        if (2, 1), do the other way around.

    Returns
    -------
    rows, cols, vals: list
        Each matched pair of rows[i], cols[i], their score is vals[i].
    """
    res = [[], [], []]
    if order == (1, 2):
        idx1_to_idx2 = dict()
        idx1_to_score = dict()
        for i, j, score in zip(matching[0], matching[1], matching[2]):
            if i not in idx1_to_idx2:
                idx1_to_idx2[i] = j
                idx1_to_score[i] = score
            elif score > idx1_to_score[i]:
                idx1_to_idx2[i] = j
                idx1_to_score[i] = score
        for idx1, idx2 in idx1_to_idx2.items():
            res[0].append(idx1)
            res[1].append(idx2)
            res[2].append(idx1_to_score[idx1])
    elif order == (2, 1):
        idx2_to_idx1 = dict()
        idx2_to_score = dict()
        for i, j, score in zip(matching[0], matching[1], matching[2]):
            if j not in idx2_to_idx1:
                idx2_to_idx1[j] = i
                idx2_to_score[j] = score
            elif score > idx2_to_score[j]:
                idx2_to_idx1[j] = i
                idx2_to_score[j] = score
        for idx2, idx1 in idx2_to_idx1.items():
            res[0].append(idx1)
            res[1].append(idx2)
            res[2].append(idx2_to_score[idx2])
    else:
        raise NotImplementedError('order must be in {(1, 2), (2, 1)}.')

    return res


def match_cells(arr1, arr2, verbose=True):
    """
    Get matching between arr1 and arr2 using linear assignment, the distance is 1 - Pearson correlation.
    Parameters
    ----------
    arr1: np.array of shape (n_samples1, n_features)
        The first data matrix
    arr2: np.array of shape (n_samples1, n_features)
        The second data matrix
    verbose: bool, default=True
        Whether to print the progress

    Returns
    -------
    rows, cols, vals: list
        Each matched pair of rows[i], cols[i], their distance is vals[i]
    """
    if verbose:
        print('Start the matching process...', flush=True)
        print('Computing the distance matrix...', flush=True)
    dist = utils.cdist_correlation(arr1, arr2)
    if verbose:
        print('Solving linear assignment...', flush=True)
    rows, cols = linear_sum_assignment(dist)
    if verbose:
        print('Linear assignment completed!', flush=True)

    return rows, cols, np.array([dist[i, j] for i, j in zip(rows, cols)])


def get_initial_matching(
        arr1, arr2,
        clust_labels1=None, clust_labels2=None,
        edges1=None, edges2=None,
        wt1=0.3, wt2=0.3,
        randomized_svd=True,
        svd_runs=1,
        svd_components1=None, svd_components2=None,
        verbose=True
):
    """
    Assume the features of arr1 and arr2 are column-wise directly comparable,
        obtain a matching by minimizing the correlation distance between arr1 and arr2.
    Parameters
    ----------
    arr1: np.array of shape (n_samples1, n_features1)
        The first data matrix.
    arr2: np.array of shape (n_samples2, n_features2)
        The second data matrix.
    clust_labels1: None or np.array of shape (n_samples1, )
        If not None, then it is the clustering label of the first data matrix,
        and the smoothing of this matrix will be done via cluster centroid shrinkage.
    clust_labels2: None or np.array of shape (n_samples2, )
        Same as clust_labels1 but for the second data matrix.
    edges1: None or list of length two or three
        If not None, then each edge in the graph is (edges[0][i], edges[1][i]) with weight edges[2][i] (if exists)
        and the smoothing of this matrix will be done via graph smoothing.
    edges2: None or scipy.sparse.csr_matrix of shape (n_samples2, n_samples2)
        Same as edges1 but for the second data matrix.
    wt1: float, default=0.3
        The smoothing of the first data matrix will be wt1 * arr1 + (1-wt1) * shrinkage_targets,
        where the shrinkage_targets will be either the cluster centroids or the average of graph neighbors.
    wt2: float, default=0.3
        Same as wt1 but for the second data matrix.
    randomized_svd: bool, default=False
        Whether to use randomized svd.
    svd_runs: int, default=1
        Randomized SVD will result in different runs,
        so if randomized_svd=True, perform svd_runs many randomized SVDs, and pick the one with the
        smallest Frobenious reconstruction error.
        If randomized_svd=False, svd_runs is forced to be 1.
    svd_components1: None or int
        If None, then do not do SVD,
        else, number of components to keep when doing SVD de-noising for the first data matrix.
    svd_components2: None or int
        Same as svd_components1 but for the second data matrix.
    verbose: bool, default=True
        Whether to print the progress.

    Returns
    -------
    matching: list of length 3
        rows, cols, vals = matching,
        Each matched pair is rows[i], cols[i], their distance is vals[i].
    """
    assert arr1.shape[1] == arr2.shape[1]
    # labels and edges cannot be specified simultaneously
    assert (clust_labels1 is None) or (edges1 is None)
    assert (clust_labels2 is None) or (edges2 is None)

    arr1, arr2 = utils.drop_zero_variability_columns(arr_lst=[arr1, arr2])

    # smoothing and denoising
    if verbose:
        print("Denoising the data...", flush=True)

    if clust_labels1 is not None:
        arr1 = utils.shrink_towards_centroids(arr=arr1, labels=clust_labels1, wt=wt1)
    elif edges1 is not None:
        arr1 = utils.graph_smoothing(arr=arr1, edges=edges1, wt=wt1)
    arr1 = utils.svd_denoise(
        arr=arr1, n_components=svd_components1, randomized=randomized_svd,
        n_runs=svd_runs
    )

    if clust_labels2 is not None:
        arr2 = utils.shrink_towards_centroids(arr=arr2, labels=clust_labels2, wt=wt2)
    elif edges2 is not None:
        arr2 = utils.graph_smoothing(arr=arr2, edges=edges2, wt=wt2)
    arr2 = utils.svd_denoise(
        arr=arr2, n_components=svd_components2, randomized=randomized_svd,
        n_runs=svd_runs
    )

    res = match_cells(arr1=arr1, arr2=arr2, verbose=verbose)
    if verbose:
        print('Initial matching completed!', flush=True)

    return res


def get_refined_matching_one_iter(
        init_matching, arr1, arr2,
        clust_labels1=None, clust_labels2=None,
        edges1=None, edges2=None,
        wt1=0.5, wt2=0.5,
        filter_prop=0,
        cca_components=15,
        cca_max_iter=2000,
        verbose=True
):
    """
    Run one iteration of CCA refinement.
    Parameters
    ----------
    init_matching: list
        init_matching[0][i], init_matching[1][i] is a matched pair,
        and init_matching[2][i] is the distance for this pair
    arr1: np.array of shape (n_samples1, n_features1)
        The first data matrix.
    arr2: np.array of shape (n_samples2, n_features2)
        The second data matrix.
    clust_labels1: None or np.array of shape (n_samples1, )
        If not None, then it is the clustering label of the first data matrix,
        and the smoothing of this matrix will be done via cluster centroid shrinkage.
    clust_labels2: None or np.array of shape (n_samples2, )
        Same as clust_labels1 but for the second data matrix.
    edges1: None or list of length two or three
        If not None, then each edge in the graph is (edges[0][i], edges[1][i]) with weight edges[2][i] (if exists)
        and the smoothing of this matrix will be done via graph smoothing.
    edges2: None or scipy.sparse.csr_matrix of shape (n_samples2, n_samples2)
        Same as edges1 but for the second data matrix.
    wt1: float, default=0.5
        The smoothing of the first data matrix will be wt1 * (cca embedding of arr1) + (1-wt1) * shrinkage_targets,
        where the shrinkage_targets will be either the cluster centroids or the average of graph neighbors.
    wt2: float, default=0.5
        Same as wt1 but for the second data matrix.
    filter_prop: float, default=0
        Proportion of matched pairs to discard before feeding into refinement iterations.
    cca_components: int, default=15
        Number of CCA components.
    cca_max_iter: int, default=2000
        Maximum number of CCA iterations.
    verbose: bool, default=True
        Whether to print the

    Returns
    -------
    rows, cols, vals: list
        Each matched pair of rows[i], cols[i], their distance is vals[i]
    """
    if verbose:
        print('Fitting CCA...', flush=True)
    arr1_cca, arr2_cca, _ = utils.cca_embedding(
        arr1=arr1, arr2=arr2,
        init_matching=init_matching, filter_prop=filter_prop, n_components=cca_components, max_iter=cca_max_iter
    )
    if verbose:
        print('Denoising the data...', flush=True)
    if clust_labels1 is not None:
        arr1_cca = utils.shrink_towards_centroids(arr=arr1_cca, labels=clust_labels1, wt=wt1)
    elif edges1 is not None:
        arr1_cca = utils.graph_smoothing(arr=arr1_cca, edges=edges1, wt=wt1)
    if clust_labels2 is not None:
        arr2_cca = utils.shrink_towards_centroids(arr=arr2_cca, labels=clust_labels2, wt=wt2)
    elif edges2 is not None:
        arr2_cca = utils.graph_smoothing(arr=arr2_cca, edges=edges2, wt=wt2)

    return match_cells(arr1=arr1_cca, arr2=arr2_cca, verbose=verbose)


def get_refined_matching(
        init_matching, arr1, arr2,
        randomized_svd=False, svd_runs=1,
        svd_components1=None, svd_components2=None,
        clust_labels1=None, clust_labels2=None,
        edges1=None, edges2=None,
        wt1=0.5, wt2=0.5,
        n_iters=3, filter_prop=0,
        cca_components=15,
        cca_max_iter=2000,
        verbose=True
):
    """
    Refinement of init_matching.
    Parameters
    ----------
    init_matching: list
        init_matching[0][i], init_matching[1][i] is a matched pair,
        and init_matching[2][i] is the distance for this pair.
    arr1: np.array of shape (n_samples1, n_features1)
        The first data matrix.
    arr2: np.array of shape (n_samples2, n_features2)
        The second data matrix.
    randomized_svd: bool, default=False
        Whether to use randomized SVD
    svd_runs: int, default=1
        Randomized SVD will result in different runs,
        so if randomized_svd=True, perform svd_runs many randomized SVDs, and pick the one with the
        smallest Frobenious reconstruction error.
        If randomized_svd=False, svd_runs is forced to be 1.
    svd_components1: None or int
        If None, then do not do SVD,
        else, number of components to keep when doing SVD de-noising for the first data matrix
        before feeding into CCA.
    svd_components2: None or int
        Same as svd_components1 but for the second data matrix.
    clust_labels1: None or np.array of shape (n_samples1, )
        If not None, then it is the clustering label of the first data matrix,
        and the smoothing of this matrix will be done via cluster centroid shrinkage.
    clust_labels2: None or np.array of shape (n_samples2, )
        Same as clust_labels1 but for the second data matrix.
    edges1: None or list of length two or three
        If not None, then each edge in the graph is (edges[0][i], edges[1][i]) with weight edges[2][i] (if exists)
        and the smoothing of this matrix will be done via graph smoothing.
    edges2: None or scipy.sparse.csr_matrix of shape (n_samples2, n_samples2)
        Same as edges1 but for the second data matrix.
    wt1: float, default=0.5
        The smoothing of the first data matrix will be wt1 * (cca embedding of arr1) + (1-wt1) * shrinkage_targets,
        where the shrinkage_targets will be either the cluster centroids or the average of graph neighbors.
    wt2: float, default=0.5
        Same as wt1 but for the second data matrix.
    n_iters: int, default=3
        Number of refinement iterations.
    filter_prop: float, default=0
        Proportion of matched pairs to discard before feeding into refinement iterations.
    cca_components: int, default=15
        Number of CCA components.
    cca_max_iter: int, default=2000,
        Maximum number of CCA iterations.
    verbose: bool, default=True
        Whether to print the progress.

    Returns
    -------
    matching: list of length 3
        rows, cols, vals = matching,
        Each matched pair is rows[i], cols[i], their distance is vals[i].
    """
    ns = [len(x) for x in init_matching]
    assert ns[0] == ns[1] == ns[2]
    # labels and edges can not be specified simultaneously
    assert (clust_labels1 is None) or (edges1 is None)
    assert (clust_labels2 is None) or (edges2 is None)
    assert isinstance(n_iters, int) and n_iters >= 1
    assert 0 <= int(ns[0] * filter_prop) < ns[0]

    assert 1 <= cca_components <= min(arr1.shape[1], arr2.shape[1])

    arr1 = utils.drop_zero_variability_columns(arr_lst=[arr1])[0]
    arr2 = utils.drop_zero_variability_columns(arr_lst=[arr2])[0]

    if verbose:
        print('Normalizing and reducing the dimension of the data...', flush=True)
    arr1_svd = utils.svd_embedding(
        arr=arr1, n_components=svd_components1,
        randomized=randomized_svd, n_runs=svd_runs
    )
    arr2_svd = utils.svd_embedding(
        arr=arr2, n_components=svd_components2,
        randomized=randomized_svd, n_runs=svd_runs
    )

    cca_matching = init_matching
    # iterative refinement
    for it in range(n_iters):
        if verbose:
            print('Now at iteration {}:'.format(it), flush=True)
        cca_matching = get_refined_matching_one_iter(
            init_matching=cca_matching,
            arr1=arr1_svd, arr2=arr2_svd, clust_labels1=clust_labels1,
            clust_labels2=clust_labels2, edges1=edges1, edges2=edges2,
            wt1=wt1, wt2=wt2, filter_prop=filter_prop,
            cca_components=cca_components, cca_max_iter=cca_max_iter, verbose=verbose
        )

    arr1_cca, arr2_cca, _ = utils.cca_embedding(
        arr1=arr1_svd, arr2=arr2_svd,
        init_matching=cca_matching,
        filter_prop=filter_prop,
        n_components=cca_components,
        max_iter=2000)

    if verbose:
        print('Refined matching completed!', flush=True)
    return cca_matching
