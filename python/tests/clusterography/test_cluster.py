import numpy as np

import libcasm.clusterography as clust


def test_cluster():
    cluster = clust.Cluster.from_list(
        [
            [0, 0, 0, 0],
            [0, 1, 0, 0],
        ]
    )
    cluster.append([0, 2, 0, 0])
    assert len(cluster.sites()) == 3
    assert cluster.size() == 3


def test_cluster_list():
    cluster = clust.Cluster.from_list(
        [
            [0, 0, 0, 0],
            [0, 1, 2, 3],
        ]
    )
    assert cluster.to_list() == [[0, 0, 0, 0], [0, 1, 2, 3]]


def test_cluster_translate_add():
    cluster = clust.Cluster.from_list(
        [
            [0, 1, 0, 0],
        ]
    )
    translation = np.array([0, 0, 1])
    transformed = cluster + translation
    assert transformed.to_list() == [[0, 1, 0, 1]]


def test_cluster_translate_iadd():
    cluster = clust.Cluster.from_list(
        [
            [0, 1, 0, 0],
        ]
    )
    translation = np.array([0, 0, 1])
    cluster += translation
    assert cluster.to_list() == [[0, 1, 0, 1]]


def test_cluster_translate_sub():
    cluster = clust.Cluster.from_list(
        [
            [0, 1, 0, 0],
        ]
    )
    translation = np.array([0, 0, 1])
    transformed = cluster - translation
    assert transformed.to_list() == [[0, 1, 0, -1]]


def test_cluster_translate_isub():
    cluster = clust.Cluster.from_list(
        [
            [0, 1, 0, 0],
        ]
    )
    translation = np.array([0, 0, 1])
    cluster -= translation
    assert cluster.to_list() == [[0, 1, 0, -1]]
