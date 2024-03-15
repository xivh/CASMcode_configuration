import libcasm.clusterography as clust
import libcasm.occ_events as occ_events
import libcasm.sym_info as sym_info


def test_local_cluster_specs(fcc_1NN_A_Va_event):
    xtal_prim, phenomenal_occ_event = fcc_1NN_A_Va_event
    prim_factor_group = sym_info.make_factor_group(xtal_prim)
    symgroup_rep = occ_events.make_occevent_symgroup_rep(
        prim_factor_group.elements, xtal_prim
    )
    occevent_group = occ_events.make_occevent_group(
        occ_event=phenomenal_occ_event,
        group=prim_factor_group,
        lattice=xtal_prim.lattice(),
        occevent_symgroup_rep=symgroup_rep,
    )
    cluster_specs = clust.ClusterSpecs(
        xtal_prim=xtal_prim,
        generating_group=occevent_group,
        max_length=[0.0, 0.0],
        phenomenal=phenomenal_occ_event.cluster(),
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


def test_make_occevent_local_orbits(fcc_1NN_A_Va_event):
    xtal_prim, phenomenal_occ_event = fcc_1NN_A_Va_event
    cluster_specs = occ_events.make_occevent_cluster_specs(
        xtal_prim=xtal_prim,
        phenomenal_occ_event=phenomenal_occ_event,
        max_length=[0.0, 0.0],
        cutoff_radius=[0.0, 2.01],
    )
    orbits = cluster_specs.make_orbits()

    assert len(orbits) == 5
    assert len(orbits[0]) == 1
    assert len(orbits[1]) == 4
    assert len(orbits[2]) == 4
    assert len(orbits[3]) == 8
    assert len(orbits[4]) == 2


def test_make_phenomenal_occevent():
    assert True
