"""
General utility functions
"""

import warnings
from collections import defaultdict, Counter
import numpy as np
import anndata as ad
import scanpy as sc
from scipy.sparse.linalg import svds
from sklearn.cross_decomposition import CCA
from sklearn.utils.extmath import randomized_svd


def center_scale(arr):
    """
    Column-wise center and scale by standard deviation.

    Parameters
    ----------
    arr: np.ndarray of shape (n_samples, n_features)

    Returns
    -------
    Center and scaled version of arr.
    """
    return (arr - arr.mean(axis=0)) / arr.std(axis=0)


def sort_dict(dict_):
    """
    Return a copy of dict_ with both keys and values sorted.

    Parameters
    ----------
    dict_: dict

    Returns
    -------
    A sorted copy of dict_
    """
    res = defaultdict(list)
    for key in sorted(dict_.keys()):
        # sort according to decreasing order in scores
        res[key] = sorted(dict_[key], key=lambda x: -x[1])
    return res


def dict_to_list(dict_):
    """
    Convert dict_ into a list.

    Parameters
    ----------
    dict_: dict
        A dict representing a matching. The structure is {row: [(col1, val1), (col2, val2), ...]}

    Returns
    -------
    res
        A list of length three, say (rows, cols, vals). Each element is further a list and
        each matched pair is (rows[i], cols[i]), and their similarity score is vals[i].
    """
    res = [[], [], []]
    for curr_idx1, indices2_and_scores in dict_.items():
        for curr_idx2, curr_score in indices2_and_scores:
            res[0].append(curr_idx1)
            res[1].append(curr_idx2)
            res[2].append(curr_score)
    return res


def list_to_dict(list_):
    """
    Convert list_ to a dict.

    Parameters
    ----------
    list_: list of length three
        rows, cols, vals = list_, each element is further a list and
        each matched pair is (rows[i], cols[i]), and their similarity score is vals[i].

    Returns
    -------
    A dict representing a matching.
    The structure is {row: [(col1, val1), (col2, val2), ...]}.
    """
    res = defaultdict(list)
    for idx1, idx2, score in zip(list_[0], list_[1], list_[2]):
        res[idx1].append((idx2, score))
    return res


def summarize_clustering(clustering, true_labels):
    """
    Compute the majority cell type for each cluster.

    Parameters
    ----------
    clustering: np.array of shape (n_samples,)
        Clustering labels, coded from 0, 1, ..., n_clusters.
    true_labels: np.array of shape (n_samples,)
        Groundtruth labels.

    Returns
    -------
    np.array of shape (n_clusters,)
        The majority voting results.
    """
    res = []
    cluster_label_to_indices = defaultdict(list)
    for i, l in enumerate(clustering):
        cluster_label_to_indices[l].append(i)
    for i in range(len(cluster_label_to_indices)):
        curr_indices = cluster_label_to_indices[i]
        curr_true_labels = true_labels[curr_indices]
        # majority voting
        # randomly break the ties
        curr_true_labels = np.random.permutation(curr_true_labels)
        counter = Counter(curr_true_labels)
        most_common_element, _ = counter.most_common(1)[0]
        res.append(most_common_element)

    return res


def drop_zero_variability_columns(arr_lst: list, tol=1e-8):
    """
    Drop columns for which its standard deviation is zero in any one of the arrays in arr_list.

    Parameters
    ----------
    arr_lst: list of np.array
        List of arrays
    tol: float, default=1e-8
        Any number less than tol is considered as zero

    Returns
    -------
    List of np.array where no column has zero standard deviation
    """
    assert all([arr_lst[i].shape[1] == arr_lst[0].shape[1] for i in range(len(arr_lst))])
    bad_columns = set()
    for arr in arr_lst:
        curr_std = np.std(arr, axis=0)
        for col in np.nonzero(np.abs(curr_std) < tol)[0]:
            bad_columns.add(col)
    good_columns = [i for i in range(arr_lst[0].shape[1]) if i not in bad_columns]
    return [arr[:, good_columns] for arr in arr_lst]


def recode(labels):
    """
    Recode the cluster labels to 0, 1, ..., num_clusters-1

    Parameters
    ----------
    labels: np.array of shape (n_samples,)
        Cluster labels, each element must be hashable

    Returns
    -------
    new_labels: np.array of shape (n_samples,)
        Recoded cluster labels
    new_to_old: dict
        Dictionary of {new_label: old_label}

    """
    unique_labels = np.unique(labels)
    old_to_new = {old: new for new, old in enumerate(unique_labels)}
    new_to_old = {new: old for new, old in enumerate(unique_labels)}
    new_labels = []
    for l in labels:
        new_labels.append(old_to_new[l])
    return np.array(new_labels), new_to_old


def shrink_towards_centroids(arr, labels, wt):
    """
    For each row of arr, shrink it towards its cluster centroid by taking wt*raw_data + (1-wt)*centroid

    Parameters
    ----------
    arr: np.array of shape (n_samples, n_features)
        Data matrix
    labels: np.array of shape (n_samples,)
        Cluster labels
    wt: float
        Weight for shrinkage

    Returns
    -------
    denoised_arr: np.array of shape (n_samples, n_features)
        Original array after centroid shrinkage
    """
    labels, _ = recode(labels)
    centroids = get_centroids(arr, labels)
    centroids = centroids[labels, :]
    return wt * arr + (1 - wt) * centroids


def graph_smoothing(arr, edges, wt):
    """
    For each row of arr, shrink it towards the average of its neighborhood by taking wt*raw_data + (1-wt)*nhbd_avg

    Parameters
    ----------
    arr: np.array of shape (n_samples, n_features)
        Data matrix
    edges: list of length two or three
        Each edge of the graph is (edges[0][i], edges[1][i]) and the weight is edges[2][i] (if exists)
    wt: float
        Weight for shrinkage

    Returns
    -------
    denoised_arr: np.array of shape (n_samples, n_features)
        Original array after graph_smoothing
    """
    n = arr.shape[0]
    adj_list = [[] for _ in range(n)]
    wt_list = [[] for _ in range(n)]
    for i in range(len(edges[0])):
        adj_list[edges[0][i]].append(edges[1][i])
        if len(edges) > 2:
            wt_list[edges[0][i]].append(edges[2][i])

    # for i in range(len(adj.indptr)-1):
    #     col_indices = adj.indices[adj.indptr[i]:adj.indptr[i+1]]
    #     if len(col_indices) == 0:
    #         adj_list.append(np.array([i]))
    #     else:
    #         adj_list.append(adj.indices[adj.indptr[i]:adj.indptr[i+1]])
    centroids = []
    for i in range(n):
        if len(edges) == 2:
            centroids.append(np.mean(arr[adj_list[i], :], axis=0))
        else:
            centroids.append(np.average(arr[adj_list[i], :], axis=0, weights=wt_list[i]))

    centroids = np.array(centroids)
    return wt * arr + (1-wt) * centroids


def get_centroids(arr, labels):
    """
    Compute the centroids (cluster mean) of arr.

    Parameters
    ----------
    arr: np.array of shape (n_samples, n_features)
        Data matrix
    labels: np.array of shape (n_samples,)
        Cluster labels of each sample, coded from 0, ..., num_clusters-1

    Returns
    -------
    centroids: np.array of shape (n_centroids, n_features)
        Matrix of cluster centroids, the i-th column is the centroid of the i-th cluster
    """
    cluster_label_to_indices = defaultdict(list)
    for i, l in enumerate(labels):
        cluster_label_to_indices[l].append(i)

    unique_labels = sorted(cluster_label_to_indices.keys())
    if not all(i == l for i, l in enumerate(unique_labels)):
        raise ValueError('labels must be coded in integers from 0, ..., n_clusters-1.')

    centroids = np.empty((len(unique_labels), arr.shape[1]))
    for curr_label, indices in cluster_label_to_indices.items():
        centroids[curr_label, :] = arr[indices, :].mean(axis=0)
    return centroids


def robust_svd(arr, n_components, randomized=False, n_runs=1):
    """
    Do deterministic or randomized SVD on arr.

    Parameters
    ----------
    arr: np.array
        The array to do SVD on
    n_components: int
        Number of SVD components
    randomized: bool, default=False
        Whether to run randomized SVD
    n_runs: int, default=1
        Run multiple times and take the realization with the lowest Frobenious reconstruction error

    Returns
    -------
    u, s, vh: np.array
        u @ np.diag(s) @ vh is the reconstruction of the original arr
    """
    if randomized:
        best_err = float('inf')
        u, s, vh = None, None, None
        for _ in range(n_runs):
            curr_u, curr_s, curr_vh = randomized_svd(arr, n_components=n_components, random_state=None)
            curr_err = np.sum((arr - curr_u @ np.diag(curr_s) @ curr_vh) ** 2)
            if curr_err < best_err:
                best_err = curr_err
                u, s, vh = curr_u, curr_s, curr_vh
        assert u is not None and s is not None and vh is not None
    else:
        if n_runs > 1:
            warnings.warn("Doing deterministic SVD, n_runs reset to one.")
        u, s, vh = svds(arr*1.0, k=n_components) # svds can not handle integer values
    return u, s, vh


def svd_denoise(arr, n_components=20, randomized=False, n_runs=1):
    """
    Compute best rank-n_components approximation of arr by SVD.

    Parameters
    ----------
    arr: np.array of shape (n_samples, n_features)
        Data matrix
    n_components: int, default=20
        Number of components to keep
    randomized: bool, default=False
        Whether to use randomized SVD
    n_runs: int, default=1
        Run multiple times and take the realization with the lowest Frobenious reconstruction error

    Returns
    -------
    arr: array_like of shape (n_samples, n_features)
        Rank-n_comopnents approximation of the input arr.
    """
    if n_components is None:
        return arr
    u, s, vh = robust_svd(arr, n_components=n_components, randomized=randomized, n_runs=n_runs)
    return u @ np.diag(s) @ vh


def svd_embedding(arr, n_components=20, randomized=False, n_runs=1):
    """
    Compute rank-n_components SVD embeddings of arr.

    Parameters
    ----------
    arr: np.array of shape (n_samples, n_features)
        Data matrix
    n_components: int, default=20
        Number of components to keep
    randomized: bool, default=False
        Whether to use randomized SVD
    n_runs: int, default=1
        Run multiple times and take the realization with the lowest Frobenious reconstruction error

    Returns
    -------
    embeddings: array_like of shape (n_samples, n_components)
        Rank-n_comopnents SVD embedding of arr.
    """
    if n_components is None:
        return arr
    u, s, vh = robust_svd(arr, n_components=n_components, randomized=randomized, n_runs=n_runs)
    return u @ np.diag(s)


def cdist_correlation(arr1, arr2):
    """Calculate pair-wise 1 - Pearson correlation between X and Y.

    Parameters
    ----------
    arr1: np.array of shape (n_samples1, n_features)
        First dataset.
    arr2: np.array of shape (n_samples2, n_features)
        Second dataset.

    Returns
    -------
    array-like of shape (n_samples1, n_samples2)
        The (i, j)-th entry is 1 - Pearson correlation between i-th row of arr1 and j-th row of arr2.
    """
    n, p = arr1.shape
    m, p2 = arr2.shape
    assert p2 == p

    arr1 = (arr1.T - np.mean(arr1, axis=1)).T
    arr2 = (arr2.T - np.mean(arr2, axis=1)).T

    arr1 = (arr1.T / np.sqrt(1e-6 + np.sum(arr1 ** 2, axis=1))).T
    arr2 = (arr2.T / np.sqrt(1e-6 + np.sum(arr2 ** 2, axis=1))).T

    return 1 - arr1 @ arr2.T


def pearson_correlation(arr1, arr2):
    """Calculate the vector of pearson correlations between each row of arr1 and arr2.

    Parameters
    ----------
    arr1: np.array of shape (n_samples, n_features)
        First dataset.
    arr2: np.array of shape (n_samples, n_features)
        Second dataset.

    Returns
    -------
    np.array of shape (n_samples,), the i-th entry is the pearson correlation between arr1[i, :] and arr2[i, :].
    """
    n, p = arr1.shape
    m, p2 = arr2.shape
    assert n == m and p2 == p

    arr1 = (arr1.T - np.mean(arr1, axis=1)).T
    arr2 = (arr2.T - np.mean(arr2, axis=1)).T

    arr1 = (arr1.T / np.sqrt(1e-6 + np.sum(arr1 ** 2, axis=1))).T
    arr2 = (arr2.T / np.sqrt(1e-6 + np.sum(arr2 ** 2, axis=1))).T

    return np.sum(arr1 * arr2, axis=1)


def filter_bad_matches(matching, filter_prop=0.1):
    """
    Filter bad matches according to the distances of matched pairs.

    Parameters
    ----------
    matching: list
        rows, cols, vals = init_matching, where each matched pair is (rows[i], cols[i]),
        and their distance is vals[i]
    filter_prop: float
        Matched pairs with distance in top filter_prop are discarded
    Returns
    -------
    rows, cols, vals: list
        Each matched pair of rows[i], cols[i], their distance is vals[i]
    """
    init_rows, init_cols, init_vals = matching
    thresh = np.quantile(init_vals, 1 - filter_prop)
    rows = []
    cols = []
    vals = []
    for i, j, val in zip(init_rows, init_cols, init_vals):
        if val < thresh:
            rows.append(i)
            cols.append(j)
            vals.append(val)
    return np.array(rows, dtype=np.int32), np.array(cols, dtype=np.int32), np.array(vals, dtype=np.float32)


def cca_embedding(arr1, arr2, init_matching, filter_prop, n_components, max_iter=2000):
    """
    Filter bad matched pairs, align arr1 and arr2 using init_matching, fit CCA, and get CCA embeddings of arr1 and arr2.

    Parameters
    ----------
    arr1: np.ndarray of shape (n_samples1, n_features1)
        The first data matrix
    arr2: np.ndarray of shape (n_samples2, n_features2)
        The second data matrix
    init_matching: list
        rows, cols, vals = init_matching, where each matched pair is (rows[i], cols[i]),
        and their distance is vals[i]
    filter_prop: float
        Matched pairs with distance in top filter_prop are discarded when fitting CCA
    n_components: int
        Number of components to keep when fitting CCA
    max_iter: int, default=2000
        Maximum number of iterations for CCA

    Returns
    -------
    arr1_cca: np.array of shape (n_samples1, n_components)
    arr2_cca: np.array of shape (n_samples2, n_components)
    canonical_correlations: np.array of shape (n_components,)
    """

    # filter bad matched pairs
    arr1_indices, arr2_indices, _ = filter_bad_matches(init_matching, filter_prop)

    # align
    arr1_aligned = arr1[arr1_indices, :]
    arr2_aligned = arr2[arr2_indices, :]

    # cca
    cca = CCA(n_components=n_components, max_iter=max_iter)
    cca.fit(arr1_aligned, arr2_aligned)
    arr1_aligned_cca, arr2_aligned_cca = cca.transform(arr1_aligned, arr2_aligned)
    arr1_aligned_cca = center_scale(arr1_aligned_cca)
    arr2_aligned_cca = center_scale(arr2_aligned_cca)

    canonical_correlations = np.corrcoef(
        arr1_aligned_cca, arr2_aligned_cca, rowvar=False).diagonal(offset=n_components)
    arr1_cca, arr2_cca = cca.transform(arr1, arr2)
    arr1_cca = center_scale(arr1_cca)
    arr2_cca = center_scale(arr2_cca)

    return arr1_cca, arr2_cca, canonical_correlations


def process_count_data(arr, target_sum=1e4, min_mean=0.0125, max_mean=3, min_disp=0.5, max_value=10):
    """
    Process count data according to scanpy pipeline.

    Parameters
    ----------
    arr: np.array of shape (n_samples, n_features)
        Data matrix
    target_sum: float, default=1e4
        Parameter in sc.pp.normalize_total
    min_mean: float, default=0.0125
        Parameter in sc.pp.higly_variable_genes
    max_mean: float, default=3
        Parameter in sc.pp.higly_variable_genes
    min_disp: float, default=0.5
        Parameter in sc.pp.higly_variable_genes
    max_value: float, default=10
        Parameter in sc.pp.scale

    Returns
    -------
    An np.array representing the processed version of arr
    """
    adata = ad.AnnData(arr, dtype=np.float32)
    sc.pp.normalize_total(adata, target_sum=target_sum)
    sc.pp.log1p(adata)
    sc.pp.highly_variable_genes(adata, min_mean=min_mean, max_mean=max_mean, min_disp=min_disp)
    adata = adata[:, adata.var.highly_variable].copy()
    sc.pp.scale(adata, max_value=max_value)
    return adata.X
