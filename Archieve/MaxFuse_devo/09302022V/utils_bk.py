import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import pynndescent
from scipy.optimize import linear_sum_assignment, nnls
from scipy.sparse import csr_matrix, coo_matrix
from scipy.sparse.csgraph import min_weight_full_bipartite_matching
from sklearn.cross_decomposition import CCA
from sklearn.utils import check_random_state
from sklearn.utils.extmath import randomized_svd
from sknetwork.clustering import Louvain
from umap.umap_ import fuzzy_simplicial_set


def get_active_data(df, df_shared, quantile=0.5):
    """
    Get top (1-quantile) active columns of the data, according to the column-wise standard deviations.
    Parameters
    ----------
    df: pd.DataFrame of shape (n_samples, n_features)
        Dataframe
    quantile: float, default=0.5
        Any columns whose standard deviations are above this quantile will be kept.
    df_shared: pd.DataFrame of shape (n_samples, n_features) of shared features
        Different filtering standard for shared genes
    Returns
    -------
    df: pd.DataFrame of shape (n_samples, n_features_active)
        Dataframe with active features

    """
    std = df.std(axis=0)
    all_std_features = df.columns[std >= max(1e-5, np.quantile(std, quantile))]
    shared_std_features = df_shared[std > 0]
    features_use = all_std_features.union(shared_std_features)

    return df[features_use]


def louvain_clustering(df, n_neighbors, metric='correlation', resolution=1, verbose=False):
    """
    Apply Louvain modularity maximization algorithm to the UMAP graph of df.

    Parameters
    ----------
    df: array_like of shape (n_samples, n_features)
    n_neighbors: int
        Number of neighbors desired
    metric: string, default='correlation'
        Metric used when constructing the initial knn graph
    resolution: float, default=1
        Resolution parameter in Louvain algorithm
    verbose: bool, default=True
        Whether to print progress

    Returns
    -------
    labels: array_like of shape (n_samples,)
        Cluster labels
    rows, cols: list
        Each edge is rows[i], cols[i] and the distance is vals[i]
    """
    n = df.shape[0]
    rows, cols = get_feature_edges(data=df, n_neighbors=n_neighbors, metric=metric, verbose=verbose)
    labels = Louvain(resolution=resolution).fit_transform(
        csr_matrix(
            ([1] * len(rows), (rows, cols)),
            shape=(n, n)
        )
    )
    return labels, rows, cols


def get_feature_edges(data, n_neighbors, metric='correlation', verbose=False):
    """
    Compute mutual k-nearest neighbors of data and return the UMAP graph.

    Parameters
    ----------
    data: array_like of shape (n_samples, n_features)
        Data matrix
    n_neighbors: int
        Number of neighbors desired
    metric: string, default='correlation'
        Metric used when constructing the initial knn graph
    verbose: bool, default=True
        Whether to print progress

    Returns
    -------
    rows, cols: list
        Each edge is rows[i], cols[i]
    """
    data = np.array(data)
    if verbose:
        print("Doing k-NN search...", flush=True)
    index = pynndescent.NNDescent(data, n_neighbors=n_neighbors, metric=metric)
    knn_indices, knn_dists = index.neighbor_graph

    # knn_indices, knn_dists, knn_search_index = nearest_neighbors(
    #     data,
    #     n_neighbors=n_neighbors,
    #     metric=metric,
    #     metric_kwds={},
    #     angular=False,
    #     random_state=random_state,
    #     low_memory=True,
    #     use_pynndescent=True,
    #     n_jobs=1,
    #     verbose=False,
    # )

    if verbose:
        print("Constructing UMAP graph...", flush=True)
    adj, _, _ = fuzzy_simplicial_set(
        X=data,
        n_neighbors=n_neighbors,
        metric=metric,
        random_state=check_random_state(0),
        knn_indices=knn_indices,
        knn_dists=knn_dists,
    )
    adj = coo_matrix(adj)
    rows = adj.row
    cols = adj.col

    if verbose:
        print("Done!", flush=True)

    return rows, cols


def plot_edges(rows, cols, labels, markersize=0.005):
    """
    Plot adjacency matrix.
    Parameters
    ----------
    rows, cols: list
        Each edge is rows[i], cols[i]
    labels: array_like of shape (n_samples, )
        Cluster labels
    markersize: float, default=0.005
        Size of the marker for each edge
    """
    n = len(labels)
    perm = np.argsort(np.array(labels))
    old2new = {old: new for new, old in enumerate(perm)}
    ordered_rows = []
    ordered_cols = []

    for i in range(len(rows)):
        ordered_rows.append(old2new[rows[i]])
        ordered_cols.append(old2new[cols[i]])

    adj = csr_matrix(([1] * len(ordered_rows), (ordered_rows, ordered_cols)), shape=(n, n))
    plt.spy(adj, markersize=markersize)


def get_centroids(df, labels):
    """
    Compute the centroids (cluster mean) of df.

    Parameters
    ----------
    df: pd.DataFrame of shape (n_samples, n_features)
        Dataframe
    labels: array_like of shape (n_samples,)
        Cluster labels of each sample

    Returns
    -------
    centroids: pd.DataFrame of shape (n_centroids, n_features)
        Dataframe of cluster centroids
    """
    labels = np.array(labels)
    unique_labels = np.unique(labels)
    centroids = []
    for each_label in unique_labels:
        mask = labels == each_label
        centroids.append(list(df.iloc[mask].mean(axis=0)))
    centroids = pd.DataFrame(centroids, columns=df.columns, index=unique_labels)
    return centroids


def svd_denoise(arr, n_components=20):
    """
    Compute best rank-n_components approximation of arr by SVD.

    Parameters
    ----------
    arr: array_like of shape (n_samples, n_features)
        Data matrix
    n_components: int, default=20
        Number of components to keep

    Returns
    -------
    arr: array_like of shape (n_samples, n_features)
        Rank-n_comopnents approximation of the input arr.
    """
    arr = np.array(arr)
    u, s, vh = randomized_svd(arr, n_components=n_components, random_state=None)
    return u @ np.diag(s) @ vh


def svd_embedding(arr, n_components=20):
    """
    Compute rank-n_components SVD embeddings of arr.

    Parameters
    ----------
    arr: array_like of shape (n_samples, n_features)
        Data matrix
    n_components: int, default=20
        Number of components to keep

    Returns
    -------
    embeddings: array_like of shape (n_samples, n_components)
        Rank-n_comopnents SVD embedding of arr.
    """
    arr = np.array(arr)
    u, s, vh = randomized_svd(arr, n_components=n_components, random_state=None)
    return u @ np.diag(s)


def shrink_towards_centroids(df, labels, wt):
    """
    For each row of df, shrink it towards its cluster centroid by taking wt*raw_data + (1-wt)*centroid

    Parameters
    ----------
    df: pd.DataFrame of shape (n_samples, n_features)
        Data matrix
    labels: array_like of shape (n_samples,)
        Cluster labels
    wt: float
        Weight for shrinkage

    Returns
    -------
    denoised_df: pd.DataFrame of shape (n_samples, n_features)
        Original df after centroid shrinkage
    """

    centroids = get_centroids(df, labels)
    centroids = centroids.loc[labels]
    return pd.DataFrame(wt * df.to_numpy() + (1 - wt) * centroids.to_numpy(), columns=df.columns)


def cdist_correlation(X, Y):
    """Calculate pair-wise 1 - Pearson correlation between X and Y.

    Parameters
    ----------
    X: array-like of shape (n_samples_of_X, n_features)
        First dataset.
    Y: array-like of shape (n_samples_of_Y, n_features)
        Second dataset.

    Returns
    -------
    array-like of shape (n_samples_of_X, n_samples_of_Y)
        The (i, j)-th entry is 1 - Pearson correlation between i-th row of X and j-th row of Y.
    """
    X = np.array(X)
    Y = np.array(Y)
    n, p = X.shape
    m, p2 = Y.shape
    assert p2 == p

    X = (X.T - np.mean(X, axis=1)).T
    Y = (Y.T - np.mean(Y, axis=1)).T

    X = (X.T / np.sqrt(1e-6 + np.sum(X ** 2, axis=1))).T
    Y = (Y.T / np.sqrt(1e-6 + np.sum(Y ** 2, axis=1))).T

    return 1 - X @ Y.T


def get_matching(df1, df2, sparsity=None, metric='correlation', verbose=True):
    """
    Get matching between df1 and df2 using linear assignment.
    Parameters
    ----------
    df1: array_like of shape (n_samples, n_features)
        The first data matrix
    df2: array_like of shape (n_samples, n_features)
        The second data matrix
    sparsity: int or None, default=None
        if None, compute full correlation distance matrix for linear assignment
        if an integer, keep only kNN graph with k=sparsity for linear assignment
    metric: string, default='correlation'
        if sparsity=None, only metric='correlation' is used
        if sparsity is an integer, then can use available metrics in pynndescent package
    verbose: bool, default=True
        Whether to print the progress

    Returns
    -------
    rows, cols, vals: list
        Each matched pair of rows[i], cols[i], their distance is vals[i]
    """
    if verbose:
        print('Start the matching process...', flush=True)
    df1 = np.array(df1)
    df2 = np.array(df2)
    if sparsity is None:
        if verbose:
            print('Computing the distance matrix...', flush=True)
        dist = cdist_correlation(df1, df2)
        if verbose:
            print('Solving linear assignment...', flush=True)
        rows, cols = linear_sum_assignment(dist)
    else:
        if verbose:
            print('Computing the distance matrix...', flush=True)
        # build an initial knn graph for efficient search
        index = pynndescent.NNDescent(df2, n_neighbors=50, metric=metric)
        index.prepare()
        knn_indices, knn_dists = index.query(df1, k=sparsity)
        dist_rows = np.repeat(list(range(df1.shape[0])), sparsity)
        dist_cols = knn_indices.flatten()
        dist_vals = knn_dists.flatten()
        dist_vals = dist_vals + max(1e-6, -min(dist_vals))
        dist = csr_matrix((dist_vals, (dist_rows, dist_cols)), shape=(df1.shape[0], df2.shape[0]))
        if verbose:
            print('Solving linear assignment...', flush=True)
        rows, cols = min_weight_full_bipartite_matching(dist)

    return rows, cols, np.array([dist[i, j] for i, j in zip(rows, cols)])


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


def cca_embedding(df1, df2, init_matching, filter_prop, n_components):
    """
    Filter bad matched pairs, align df1 and df2 using init_matching, fit CCA, and get CCA embeddings of df1 and df2.
    Parameters
    ----------
    df1: array_like of shape (n_samples_1, n_features_1)
        The first dataframe
    df2: array_like of shape (n_samples_2, n_features_2)
        The second dataframe
    init_matching: list
        rows, cols, vals = init_matching, where each matched pair is (rows[i], cols[i]),
        and their distance is vals[i]
    filter_prop: float
        Matched pairs with distance in top filter_prop are discarded when fitting CCA
    n_components: int
        Number of components to keep when fitting CCA
    wt: float
        Weight for shrinkage

    Returns
    -------
    df1_cca: array_like of shape (n_samples_1, n_components)
    df2_cca: array_like of shape (n_samples_1, n_components)
    """

    # filter bad matched pairs
    df1_indices, df2_indices, _ = filter_bad_matches(init_matching, filter_prop)

    # align
    df1_aligned = df1[df1_indices]
    df2_aligned = df2[df2_indices]

    # cca
    cca = CCA(n_components=n_components, max_iter=2000)
    cca.fit(df1_aligned, df2_aligned)
    df1_cca, df2_cca = cca.transform(df1, df2)
    df1_cca = pd.DataFrame(df1_cca)
    df2_cca = pd.DataFrame(df2_cca)
    return df1_cca, df2_cca


def augment_new_columns(
        df1_shared, df2_shared, df1_additional, matched_pairs, add_thresh=0.6, drop_thresh=0.2, verbose=True
):
    """
    For each column in df1_additional, fit a non-negative least squares regression on df1_shared.
    Then apply the model to df2_shared to get a new column.
    Align this column with the column in df1_additional using an initial matching from matched paris
    and compute Pearson correlation. If the correlation exceeds add_thresh, add this predicted column to df2_shared
    and the original column in df1_additional to df1_shared.
    Now for each pair column of the augmented df1 and df2,
    drop all pairs whose estimated correlation (after aligning by matched_pairs) is below drop_thresh.

    Parameters
    ----------
    df1_shared: pd.DataFrame of shape (n_samples_1, n_features)
        The first data matrix
    df2_shared: pd.DataFrame of shape (n_samples_2, n_features)
        The second data matrix, whose columns are aligned with the columns in df1_shared
    df1_additional: pd.DataFrame of shape (n_samples_1, n_features_additional)
        Additional data matrix
    matched_pairs: list
        df1_indices, df2_indices = matched_pairs, each matched pair is (df1_indices[i], df2_indices[i])
    add_thresh: float, default=0.6
        The threshold for keeping a predicted column
    drop_thresh: float, default=0.2
        The threshold for dropping an existing column
    verbose: bool, default=True
        Whether to print the progress

    Returns
    -------
    df1_augmented: pd.DataFrame of shape (n_samples_1, n_features_augmented)
        Augmented df1_shared
    df2_augmented: pd.DataFrame of shape (n_samples_2, n_features_augmented)
        Augmented df2_shared
    """
    df1_indices, df2_indices = matched_pairs

    # identify bad correpondences
    bad_correspondences = set()
    for j in range(df1_shared.shape[1]):
        predicted_cor = np.corrcoef(
            np.array(df1_shared.iloc[:, j])[df1_indices],
            np.array(df2_shared.iloc[:, j])[df2_indices]
        )[0, 1]
        if predicted_cor < drop_thresh:
            bad_correspondences.add(j)

    df1_colname_to_predicted_df2_columns = {}
    df1_shared_colnames = df1_shared.columns
    df1_shared = np.array(df1_shared)
    df2_shared_colnames = df2_shared.columns
    df2_shared = np.array(df2_shared)

    # delete bad columns from colnames
    df1_shared_colnames = [colname for j, colname in enumerate(df1_shared_colnames) if j not in bad_correspondences]
    df2_shared_colnames = [colname for j, colname in enumerate(df2_shared_colnames) if j not in bad_correspondences]

    for j, colname in enumerate(df1_additional.columns):
        # fit nnl
        beta = nnls(df1_shared, np.array(df1_additional[colname]))[0]
        predicted_df2_column = df2_shared @ beta
        # compute correlation
        predicted_cor = np.corrcoef(
            np.array(df1_additional[colname])[df1_indices], predicted_df2_column[df2_indices]
        )[0, 1]
        # store new columns
        if predicted_cor > add_thresh:
            df1_colname_to_predicted_df2_columns[colname] = predicted_df2_column
        if verbose and j % int(df1_additional.shape[1] / 10) == 0:
            print('Completed {}/{} columns...'.format(j, df1_additional.shape[1]), flush=True)

    # add new columns and delete bad columns
    # the new columns have a suffix of '_aug'
    df1_new_columns = df1_additional[list(df1_colname_to_predicted_df2_columns.keys())]

    df1_augmented = np.concatenate(
        [df1_shared[:, [j for j in range(df1_shared.shape[1]) if j not in bad_correspondences]], df1_new_columns],
        axis=1
    )

    df1_augmented = pd.DataFrame(
        df1_augmented,
        columns=list(df1_shared_colnames) + \
                [str(colname) + '_aug' for colname in df1_colname_to_predicted_df2_columns.keys()]
    )

    df2_new_columns = pd.DataFrame(df1_colname_to_predicted_df2_columns)

    df2_augmented = np.concatenate(
        [df2_shared[:, [j for j in range(df1_shared.shape[1]) if j not in bad_correspondences]], df2_new_columns],
        axis=1
    )

    df2_augmented = pd.DataFrame(
        df2_augmented,
        columns=list(df2_shared_colnames) + \
                [str(colname) + '_aug' for colname in df1_colname_to_predicted_df2_columns.keys()]
    )

    if verbose:
        if len(df1_colname_to_predicted_df2_columns) == 0:
            print('No good augmentation columns found!', flush=True)
            return df1_shared, df2_shared
        else:
            print('Added {} good columns'.format(len(df1_colname_to_predicted_df2_columns)), flush=True)

        print('Deleted {} bad columns.'.format(len(bad_correspondences)), flush=True)
        print('Done!', flush=True)

    return df1_augmented, df2_augmented
