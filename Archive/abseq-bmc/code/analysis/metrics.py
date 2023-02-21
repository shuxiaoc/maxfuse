import numpy as np
import itertools
from collections.abc import Iterable

import match_utils


def get_matching_acc(matching, labels1, labels2, order=None):
    """
    Compute the cluster level matching accuracy.
    Parameters
    ----------
    matching: a list of length three.
        The matched pairs are (matching[0][i], matching[1][i]),
        and its score (the higher, the better) is matching[2][i].
    labels1: np.array of shape (n_samples1,)
        The first label vector.
    labels2: np.array of shape (n_samples2,)
        The first label vector.
    order: None or (1, 2) or (2, 1), default=None
        If None, then directly use matching without addressing any redundancy.
        If (1, 2), find one-to-one matching from the first dataset to the second dataset;
        if (2, 1), do the other way around.

    Returns
    -------
    Matching accuracy.
    """
    if order is None:
        return np.mean([labels1[i] == labels2[j] for i, j in zip(matching[0], matching[1])])
    matching = match_utils.address_matching_redundancy(matching=matching, order=order)
    rows, cols, _ = matching
    return np.mean([labels1[i] == labels2[j] for i, j in zip(rows, cols)])


def get_foscttm(dist, true_matching='identity'):
    """
    Compute the fraction of samples closer than true match.
    Parameters
    ----------
    dist: np.ndarray of shape (n1, n2)
        Distance matrix.
    true_matching: 'identity' or Iterable of length n1, default='identity'
        If is a list, then the ground truth matched pairs are (i, true_matching[i])
        If is 'identity', then true_matching = [0, 1..., n1].

    Returns
    -------
    The fraction of samples closer than true match.
    """
    n1, _ = dist.shape
    if true_matching == 'identity':
        true_matching = np.arange(n1)
    elif isinstance(true_matching, Iterable):
        true_matching = [i for i in true_matching]
    else:
        raise NotImplementedError('true_matching must be \'identity\' or Iterable of length dist.shape[0].')
    # mask[i, j] = True iff dist[i, j] < dist[i, true_matching[i]]
    mask = (dist.T < dist[np.arange(n1), true_matching]).T
    return np.mean(np.mean(mask, axis=1))


def get_matching_alignment_score(estimated_matching, n_samples, true_matching='identity'):
    """
    Compute the alignment between the estimated matching and the true_matching
    according to the metric in https://openproblems.bio/neurips_docs/about_tasks/task2_modality_matching/.
    Parameters
    ----------
    estimated_matching: a list of length three.
        The matched pairs are (matching[0][i], matching[1][i]),
        and its score (the higher, the better) is matching[2][i].
    n_samples: int
        The sample size for the first dataset.
    true_matching: 'identity' or Iterable of length n_samples, default='identity'
        If is a list, then the ground truth matched pairs are (i, true_matching[i])
        If is 'identity', then true_matching = [0, 1..., n_samples].

    Returns
    -------
    The alignment score.
    """
    if true_matching == 'identity':
        true_matching = np.arange(n_samples)
    elif isinstance(true_matching, Iterable):
        true_matching = [i for i in true_matching]
    else:
        raise NotImplementedError('true_matching must be \'identity\' or Iterable of length dist.shape[0].')

    idx1_to_indices2_and_scores = dict()
    for i, j, score in zip(estimated_matching[0], estimated_matching[1], estimated_matching[2]):
        if i not in idx1_to_indices2_and_scores:
            idx1_to_indices2_and_scores[i] = [[j], [score]]
        else:
            idx1_to_indices2_and_scores[i][0].append(j)
            idx1_to_indices2_and_scores[i][1].append(score)

    for idx1, indices2_and_scores in idx1_to_indices2_and_scores.items():
        indices2_and_scores[1] = list(np.array(indices2_and_scores[1]) / np.sum(indices2_and_scores[1]))

    res = 0
    for idx1, idx2 in enumerate(true_matching):
        if idx1 in idx1_to_indices2_and_scores:
            for loc in range(len(idx1_to_indices2_and_scores[idx1][0])):
                candidate_idx2 = idx1_to_indices2_and_scores[idx1][0][loc]
                if idx2 == candidate_idx2:
                    res += idx1_to_indices2_and_scores[idx1][1][loc]
    return res / len(idx1_to_indices2_and_scores)


def get_knn_alignment_score(dist, k_max, true_matching='identity'):
    """
    For each 1 <= k <= k_max, obtain knn matching from dist,
    and compute its matching proximity with the true matching.
    The proximity is calculated by:
    for each cell in arr1, claim it is successfully matched when the true match is in the k-nearest-neighborhood;
    then calculate the average success rate.

    Parameters
    ----------
    dist: np.ndarray of shape (n1, n2)
        Distance matrix.
    k_max: int
        Maximum k for knn matching.
    true_matching: 'identity' or Iterable of length n1, default='identity'
        If is a list, then the ground truth matched pairs are (i, true_matching[i])
        If is 'identity', then true_matching = [0, 1..., n1].

    Returns
    -------
    np.ndarray of shape (k_max,) representing the score for each 1<=k<=k_max.
    """
    n1, n2 = dist.shape
    assert k_max <= n2
    knn_indices = np.argsort(dist, axis=1)[:, :k_max]
    # knn_scores = 1 - dist[np.arange(n1)[:, None], knn_indices]

    if true_matching == 'identity':
        true_matching = np.arange(n1)
    elif isinstance(true_matching, Iterable):
        true_matching = [i for i in true_matching]
    else:
        raise NotImplementedError('true_matching must be \'identity\' or Iterable of length dist.shape[0].')

    res = np.zeros(k_max)
    for idx1, idx2 in enumerate(true_matching):
        candidates = knn_indices[idx1, :]
        idx2_location = np.where(candidates == idx2)[0]
        if len(idx2_location) == 0:
            # even k_max-nn matching does not contain the true match
            continue
        # find the first occurrence of idx2
        # every knn matching with k >= idx_location is able to find the true match
        idx2_location = idx2_location[0]
        # curr_scores = knn_scores[idx1, idx2_location] / np.cumsum(knn_scores[idx1, :])
        # res[idx2_location:] = res[idx2_location:] + curr_scores[idx2_location:]
        res[idx2_location:] = res[idx2_location:] + 1
    return res / n1


# this function is specifically designed for methods that does not direcly produce matching information
# but only produces a reduced embedding after integration
def get_knn_matching(dist, k_max = 1):
    """
    For each 1 <= k <= k_max, obtain knn matching from dist.
    Parameters
    ----------
    dist: np.ndarray of shape (n1, n2)
        Distance matrix.
    k_max: int
        Maximum k for knn matching.
    Returns
    -------
    knn matching indices of n1 to n2.
    """
    n1, n2 = dist.shape
    assert k_max <= n2
    knn_indices = np.argsort(dist, axis=1)[:, :k_max]
    knn_list = list(itertools.chain(*knn_indices.tolist()))
    knn_scores = []
    for i in range(n1):
        score = 1 - dist[i,knn_list[i]]
        knn_scores.append(score)
    return knn_list, knn_scores


#merged = list(itertools.chain(*lg_full_match.tolist()))


