import numpy as np

import libcasm.clusterography as clust
import libcasm.sym_info as sym_info
import libcasm.xtal as xtal
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


def test_cluster_iter():
    # construct Cluster
    A = clust.Cluster.from_list(
        [
            [0, 0, 0, 0],  # [b, i, j, k]
            [0, 0, 1, 0],
        ]
    )
    for site in A:
        assert isinstance(site, xtal.IntegralSiteCoordinate)
        assert (site in A) is True


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


def test_cluster_repr():
    import io
    from contextlib import redirect_stdout

    obj = clust.Cluster.from_list(
        [
            [0, 0, 0, 0],
            [0, 1, 0, 0],
        ]
    )

    f = io.StringIO()
    with redirect_stdout(f):
        print(obj)
    out = f.getvalue()
    assert "sites" in out


def test_cluster_distances():
    xtal_prim = xtal_prims.FCC(r=1.0, occ_dof=["A", "B", "Va"])

    # construct Cluster
    cluster_1 = clust.Cluster.from_list(
        [
            [0, 0, 0, 0],  # [b, i, j, k]
            [0, 0, 1, 0],
        ]
    )

    assert np.allclose(cluster_1.distances(xtal_prim=xtal_prim), [2.0])

    # construct Cluster
    cluster_2 = clust.Cluster.from_list(
        [
            [0, 0, 0, 0],  # [b, i, j, k]
            [0, 0, 1, 0],
            [0, 0, 2, 0],
        ]
    )

    assert np.allclose(cluster_2.distances(xtal_prim=xtal_prim), [2.0, 2.0, 4.0])

    # construct Cluster
    phenomenal_cluster = clust.Cluster.from_list(
        [
            [0, 0, -1, 0],  # [b, i, j, k]
        ]
    )

    assert np.allclose(
        cluster_1.phenomenal_distances(
            xtal_prim=xtal_prim,
            phenomenal=phenomenal_cluster,
        ),
        [2.0, 4.0],
    )


def test_IntegralClusterOrbitGenerator():
    xtal_prim = xtal_prims.FCC(r=1.0, occ_dof=["A", "B", "Va"])
    cluster = clust.Cluster.from_list(
        [
            [0, 0, 0, 0],
            [0, 1, 0, 0],
        ]
    )
    orbit_gen = clust.ClusterOrbitGenerator(
        prototype=cluster,
        include_subclusters=True,
    )
    assert isinstance(orbit_gen, clust.ClusterOrbitGenerator)

    data = orbit_gen.to_dict(xtal_prim=xtal_prim)
    assert isinstance(data, dict)

    orbit_gen2 = clust.ClusterOrbitGenerator.from_list([data], xtal_prim=xtal_prim)[0]
    assert orbit_gen2.prototype == orbit_gen.prototype
    assert orbit_gen2.include_subclusters == orbit_gen.include_subclusters

    import io
    from contextlib import redirect_stdout

    f = io.StringIO()
    with redirect_stdout(f):
        print(orbit_gen)
    out = f.getvalue()
    assert "sites" in out
