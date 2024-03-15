import libcasm.clusterography as clust
import libcasm.sym_info as sym_info
import libcasm.xtal as xtal
import libcasm.xtal.prims as xtal_prims


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

    data = cluster_specs.to_dict()
    print(xtal.pretty_json(data))
    # fmt: off
    expected_data = {
        "generating_group": [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16,
                             17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30,
                             31, 32, 33, 34, 35, 36, 37, 38, 39, 40, 41, 42, 43, 44,
                             45, 46, 47],
        "orbit_branch_specs": {
            "0": { "max_length": 0.0 },
            "1": { "max_length": 0.0 },
            "2": { "max_length": 2.01 }
        },
        "site_filter_method": "dof_sites"
    }
    # fmt: on
    assert data == expected_data

    integral_site_coordinate_symgroup_rep = (
        clust.make_integral_site_coordinate_symgroup_rep(
            group_elements=prim_factor_group.elements,
            xtal_prim=xtal_prim,
        )
    )
    cluster_specs_in = clust.ClusterSpecs.from_dict(
        data=data,
        xtal_prim=xtal_prim,
        prim_factor_group=prim_factor_group,
        integral_site_coordinate_symgroup_rep=integral_site_coordinate_symgroup_rep,
    )
    assert cluster_specs_in.generating_group().head_group_index == list(range(48))
    assert cluster_specs_in.max_length() == cluster_specs.max_length()
    assert cluster_specs_in.cutoff_radius() == cluster_specs.cutoff_radius()


def test_cluster_specs_from_dict_1():
    xtal_prim = xtal_prims.FCC(r=1.0, occ_dof=["A", "B", "Va"])
    prim_factor_group = sym_info.make_factor_group(xtal_prim)
    integral_site_coordinate_symgroup_rep = (
        clust.make_integral_site_coordinate_symgroup_rep(
            group_elements=prim_factor_group.elements,
            xtal_prim=xtal_prim,
        )
    )
    data = {
        "orbit_branch_specs": {
            "2": {"max_length": 3.01},
            "3": {"max_length": 2.01},
        },
        "site_filter_method": "dof_sites",
    }
    cluster_specs_in = clust.ClusterSpecs.from_dict(
        data=data,
        xtal_prim=xtal_prim,
        prim_factor_group=prim_factor_group,
        integral_site_coordinate_symgroup_rep=integral_site_coordinate_symgroup_rep,
    )
    assert cluster_specs_in.max_length() == [0.0, 0.0, 3.01, 2.01]


def test_cluster_specs_from_dict_2():
    xtal_prim = xtal_prims.FCC(r=1.0, occ_dof=["A", "B", "Va"])
    prim_factor_group = sym_info.make_factor_group(xtal_prim)
    integral_site_coordinate_symgroup_rep = (
        clust.make_integral_site_coordinate_symgroup_rep(
            group_elements=prim_factor_group.elements,
            xtal_prim=xtal_prim,
        )
    )
    phenomenal = clust.Cluster.from_list([[0, 0, 0, 0], [0, 1, 0, 0]])
    custom_generators = [
        clust.ClusterOrbitGenerator(
            prototype=clust.Cluster.from_list([[0, 2, 0, 0], [0, 10, 0, 0]]),
            include_subclusters=True,
        )
    ]
    cluster_specs = clust.ClusterSpecs(
        xtal_prim=xtal_prim,
        generating_group=prim_factor_group,
        max_length=[0.0, 0.0, 3.01, 2.01],
        custom_generators=custom_generators,
        phenomenal=phenomenal,
        cutoff_radius=[0.0, 0.0, 5.01, 4.01],
    )
    data = cluster_specs.to_dict()
    print(xtal.pretty_json(data))

    cluster_specs_in = clust.ClusterSpecs.from_dict(
        data=data,
        xtal_prim=xtal_prim,
        prim_factor_group=prim_factor_group,
        integral_site_coordinate_symgroup_rep=integral_site_coordinate_symgroup_rep,
    )
    assert cluster_specs_in.max_length() == [0.0, 0.0, 3.01, 2.01]
    assert len(cluster_specs_in.custom_generators()) == 1
    assert (
        cluster_specs_in.custom_generators()[0].prototype
        == custom_generators[0].prototype
    )
    assert cluster_specs_in.max_length() == [0.0, 0.0, 3.01, 2.01]
    assert cluster_specs_in.cutoff_radius() == [0.0, 0.0, 5.01, 4.01]


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
        prim_factor_group.elements, xtal_prim
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
