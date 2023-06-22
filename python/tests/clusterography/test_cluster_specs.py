import numpy as np
import libcasm.xtal as xtal
import libcasm.xtal.prims as xtal_prims
import libcasm.configuration as config
import libcasm.clusterography as clust
import libcasm.sym_info as sym_info


def test_cluster_specs():
    xtal_prim = xtal_prims.FCC(r=1.0, occ_dof=["A", "B", "Va"])
    prim_factor_group = sym_info.make_factor_group(xtal_prim)
    cluster_specs = clust.ClusterSpecs(
        xtal_prim=xtal_prim,
        generating_group=prim_factor_group,
        max_length=[0.0, 0.0, 2.01],
    )
    assert isinstance(cluster_specs, clust.ClusterSpecs)

    orbits = cluster_specs.make_orbits()
    assert len(orbits) == 3
    assert len(orbits[0]) == 1
    assert len(orbits[1]) == 1
    assert len(orbits[2]) == 6


def test_make_periodic_cluster_specs():
    xtal_prim = xtal_prims.FCC(r=1.0, occ_dof=["A", "B", "Va"])
    cluster_specs = clust.make_periodic_cluster_specs(
        xtal_prim=xtal_prim, max_length=[0.0, 0.0, 2.01]
    )
    assert isinstance(cluster_specs, clust.ClusterSpecs)

    orbits = cluster_specs.make_orbits()
    assert len(orbits) == 3
    assert len(orbits[0]) == 1
    assert len(orbits[1]) == 1
    assert len(orbits[2]) == 6


def test_local_cluster_specs():
    xtal_prim = xtal_prims.FCC(r=1.0, occ_dof=["A", "B", "Va"])
    prim_factor_group = sym_info.make_factor_group(xtal_prim)
    symgroup_rep = clust.make_integral_site_coordinate_symgroup_rep(
        prim_factor_group.elements(), xtal_prim
    )
    phenomenal = clust.Cluster(
        [
            xtal.IntegralSiteCoordinate.from_list([0, 0, 0, 0]),
            xtal.IntegralSiteCoordinate.from_list([0, 1, 0, 0]),
        ]
    )
    cluster_group = clust.make_cluster_group(
        cluster=phenomenal,
        group=prim_factor_group,
        lattice=xtal_prim.lattice(),
        integral_site_coordinate_symgroup_rep=symgroup_rep,
    )
    cluster_specs = clust.ClusterSpecs(
        xtal_prim=xtal_prim,
        generating_group=cluster_group,
        max_length=[0.0, 0.0],
        phenomenal=phenomenal,
        cutoff_radius=[0.0, 2.01],
    )
    assert isinstance(cluster_specs, clust.ClusterSpecs)

    orbits = cluster_specs.make_orbits()
    assert len(orbits) == 5
    assert len(orbits[0]) == 1
    assert len(orbits[1]) == 4
    assert len(orbits[2]) == 4
    assert len(orbits[3]) == 8
    assert len(orbits[4]) == 2


def test_make_local_cluster_specs():
    xtal_prim = xtal_prims.FCC(r=1.0, occ_dof=["A", "B", "Va"])
    phenomenal_cluster = clust.Cluster(
        [
            xtal.IntegralSiteCoordinate.from_list([0, 0, 0, 0]),
            xtal.IntegralSiteCoordinate.from_list([0, 1, 0, 0]),
        ]
    )
    cluster_specs = clust.make_local_cluster_specs(
        xtal_prim=xtal_prim,
        phenomenal_cluster=phenomenal_cluster,
        max_length=[0.0, 0.0],
        cutoff_radius=[0.0, 2.01],
    )
    assert isinstance(cluster_specs, clust.ClusterSpecs)

    orbits = cluster_specs.make_orbits()
    assert len(orbits) == 5
    assert len(orbits[0]) == 1
    assert len(orbits[1]) == 4
    assert len(orbits[2]) == 4
    assert len(orbits[3]) == 8
    assert len(orbits[4]) == 2
