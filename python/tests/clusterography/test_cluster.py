import numpy as np

import libcasm.clusterography as clust
import libcasm.sym_info as sym_info
import libcasm.xtal.prims as xtal_prims


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


def test_cluster_compare():
    # construct Cluster
    A = clust.Cluster.from_list(
        [
            [0, 0, 0, 0],  # [b, i, j, k]
            [0, 0, 1, 0],
        ]
    )
    B = clust.Cluster.from_list(
        [
            [0, 0, 0, 0],  # [b, i, j, k]
            [0, 1, 0, 0],
        ]
    )
    assert A < B
    assert A <= B
    assert A <= A
    assert B > A
    assert B >= A
    assert B >= B
    assert A == A
    assert B == B
    assert A != B


def test_cluster_rmul():
    xtal_prim = xtal_prims.FCC(r=1.0, occ_dof=["A", "B", "Va"])
    prim_factor_group = sym_info.make_factor_group(xtal_prim)
    site_factor_group_rep = clust.make_integral_site_coordinate_symgroup_rep(
        prim_factor_group.elements, xtal_prim
    )

    # construct Cluster
    cluster = clust.Cluster.from_list(
        [
            [0, 0, 0, 0],  # [b, i, j, k]
            [0, 0, 1, 0],
        ]
    )
    for site_rep in site_factor_group_rep:
        transformed_cluster = site_rep * cluster
        assert isinstance(transformed_cluster, clust.Cluster)
        assert transformed_cluster.site(0) == site_rep * cluster.site(0)
        assert transformed_cluster.site(1) == site_rep * cluster.site(1)
