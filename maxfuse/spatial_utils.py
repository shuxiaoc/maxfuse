"""
Utility functions for dealing with spatial data
"""

import numpy as np
from sklearn.neighbors import NearestNeighbors


def bind_spatial(features, nbhd, wt_on_features=0.7):
    """
    Return a new array of form [wt_on_features * features / feature_norm, (1-wt_on_features) * nbhd / nbhd_norm]

    Parameters
    ----------
    features: np.ndarray of shape (n_samples, n_features)
        Feature matrix
    nbhd: np.ndarray of shape (n_samples, n_clusters)
        Cell neighborhood composition matrix
    wt_on_features: float, default=0.7
        Weight to put on the feature matrix.

    Returns
    -------
    res: np.ndarray of shape (n_samples, n_features+n_clusters)

    """
    # normalize two kinds of info for easier tuning of weight
    feature_norm = np.linalg.norm(features)
    nbhd_norm = np.linalg.norm(nbhd)
    res = np.concatenate((
        wt_on_features * features / feature_norm,
        (1-wt_on_features) * nbhd / nbhd_norm
    ), axis=1)
    return res


def get_spatial_knn_indices(locations, n_neighbors=15, method='kd_tree'):
    """
    Compute k-nearest neighbors of locations.

    Parameters
    ----------
    locations: np.ndarray of shape (n_samples, 2)
        Data matrix
    n_neighbors: int
        Number of nearest neighbors
    method: str, default='kd_tree'
        Method to use when computing the nearest neighbors, one of ['ball_tree', 'kd_tree', 'brute']

    Returns
    -------
    knn_indices: np.ndarray of shape (n_samples, n_neighbors)
        Each row represents the knn of that sample
    """
    locations = np.array(locations)
    assert n_neighbors <= locations.shape[0]
    # k-NN indices, may be asymmetric
    _, knn_indices = NearestNeighbors(
        n_neighbors=n_neighbors, algorithm=method
    ).fit(locations).kneighbors(locations)
    return knn_indices


def get_neighborhood_composition(knn_indices, labels, log1p=False):
    """
    Compute the composition of neighbors for each sample.

    Parameters
    ----------
    knn_indices: np.ndarray of shape (n_samples, n_neighbors)
        Each row represents the knn of that sample
    labels: np.ndarray of shape (n_samples, )
        Cluster labels
    log1p: bool, default=False
        Whether to apply log1p transformation

    Returns
    -------
    comp: np.ndarray of shape (n_samples, n_neighbors)
        The composition (in proportion) of neighbors for each sample.
    """
    labels = list(labels)
    n, k = knn_indices.shape
    unique_clusters = np.unique(labels)
    n_clusters = len(unique_clusters)
    label_to_clust_idx = {label: i for i, label in enumerate(unique_clusters)}

    comp = np.zeros((n, n_clusters))
    for i, neighbors in enumerate(knn_indices):
        good_neighbors = [nb for nb in neighbors if nb != -1]
        for nb in good_neighbors:
            comp[i, label_to_clust_idx[labels[nb]]] += 1

    if log1p:
        comp = np.log1p(comp)
    return comp
