import libcasm.clusterography as clust
import libcasm.sym_info as sym_info
import libcasm.xtal as xtal
import libcasm.xtal.prims as xtal_prims


def test_make_periodic_orbit():
    xtal_prim = xtal_prims.FCC(r=1.0, occ_dof=["A", "B", "Va"])
    prim_factor_group = sym_info.make_factor_group(xtal_prim)
    factor_group_site_rep = clust.make_integral_site_coordinate_symgroup_rep(
        prim_factor_group.elements, xtal_prim
    )
    cluster = clust.Cluster.from_list(
        [
            [0, 0, 0, 0],
            [0, 1, 0, 0],
        ]
    )

    ## Check orbit ##
    orbit = clust.make_periodic_orbit(
        orbit_element=cluster,
        integral_site_coordinate_symgroup_rep=factor_group_site_rep,
    )

    assert len(orbit) == 6

    ## Check equivalence map ##
    equivalence_map = clust.make_periodic_equivalence_map(
        orbit=orbit,
        symgroup=prim_factor_group,
        lattice=xtal_prim.lattice(),
        integral_site_coordinate_symgroup_rep=factor_group_site_rep,
    )

    assert len(equivalence_map) == 6
    prototype = orbit[0]
    for i_equiv, eq_ops in enumerate(equivalence_map):
        assert len(eq_ops) == 8
        for i_op, op in enumerate(eq_ops):
            _site_rep = xtal.IntegralSiteCoordinateRep(op, xtal_prim)
            generated = (_site_rep * prototype).sorted()
            assert isinstance(op, xtal.SymOp)
            assert generated == orbit[i_equiv]

    ## Check equivalence map indices ##
    equivalence_map_indices = clust.make_periodic_equivalence_map_indices(
        orbit=orbit,
        integral_site_coordinate_symgroup_rep=factor_group_site_rep,
    )

    assert len(equivalence_map_indices) == 6
    for i_equiv, eq_indices in enumerate(equivalence_map):
        assert len(eq_indices) == 8

    # fmt: off
    expected = [
      [0, 16, 18, 23, 25, 27, 32, 47],
      [2, 7, 12, 19, 28, 33, 38, 42],
      [1, 11, 14, 17, 26, 37, 40, 41],
      [3, 6, 21, 22, 30, 31, 43, 46],
      [4, 8, 9, 20, 29, 34, 35, 44],
      [5, 10, 13, 15, 24, 36, 39, 45]
    ]
    # fmt: on

    assert equivalence_map_indices == expected


def test_make_local_orbit():
    xtal_prim = xtal_prims.FCC(r=1.0, occ_dof=["A", "B", "Va"])
    prim_factor_group = sym_info.make_factor_group(xtal_prim)
    factor_group_site_rep = clust.make_integral_site_coordinate_symgroup_rep(
        prim_factor_group.elements, xtal_prim
    )
    phenomenal_cluster = clust.Cluster.from_list(
        [
            [0, 0, 0, 0],
            [0, 1, 0, 0],
        ]
    )
    phenomenal_group = clust.make_cluster_group(
        cluster=phenomenal_cluster,
        group=prim_factor_group,
        lattice=xtal_prim.lattice(),
        integral_site_coordinate_symgroup_rep=factor_group_site_rep,
    )
    assert len(phenomenal_group.elements) == 8
    phenomenal_group_site_rep = clust.make_integral_site_coordinate_symgroup_rep(
        phenomenal_group.elements, xtal_prim
    )
    assert len(phenomenal_group_site_rep) == 8

    cluster = clust.Cluster.from_list(
        [
            [0, 0, 1, 0],
        ]
    )

    ## Check orbit ##
    orbit = clust.make_local_orbit(
        orbit_element=cluster,
        integral_site_coordinate_symgroup_rep=phenomenal_group_site_rep,
    )
    assert len(orbit) == 4

    ## Check equivalence map ##
    equivalence_map = clust.make_local_equivalence_map(
        orbit=orbit,
        phenomenal_group=phenomenal_group,
        integral_site_coordinate_symgroup_rep=phenomenal_group_site_rep,
    )

    assert len(equivalence_map) == 4
    prototype = orbit[0]
    for i_equiv, eq_ops in enumerate(equivalence_map):
        assert len(eq_ops) == 2
        for i_op, op in enumerate(eq_ops):
            _site_rep = xtal.IntegralSiteCoordinateRep(op, xtal_prim)
            generated = (_site_rep * prototype).sorted()
            assert isinstance(op, xtal.SymOp)
            assert generated == orbit[i_equiv]

    ## Check equivalence map indices ##
    equivalence_map_indices = clust.make_local_equivalence_map_indices(
        orbit=orbit,
        integral_site_coordinate_symgroup_rep=phenomenal_group_site_rep,
    )

    assert len(equivalence_map_indices) == 4
    for i_equiv, eq_indices in enumerate(equivalence_map):
        assert len(eq_indices) == 2

    # fmt: off
    expected = [
        [0, 5],
        [3, 4],
        [1, 6],
        [2, 7],
    ]
    # fmt: on

    assert equivalence_map_indices == expected

    ## Check local cluster group
    cluster_group = clust.make_local_cluster_group(
        cluster=prototype,
        phenomenal_group=phenomenal_group,
        integral_site_coordinate_symgroup_rep=phenomenal_group_site_rep,
    )
    assert len(cluster_group.elements) == 2
    assert cluster_group.head_group_index == [
        phenomenal_group.head_group_index[0],
        phenomenal_group.head_group_index[5],
    ]
