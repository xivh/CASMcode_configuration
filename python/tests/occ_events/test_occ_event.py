import libcasm.occ_events as occ_events
import libcasm.xtal as xtal
import libcasm.xtal.prims as xtal_prims


def test_OccEvent_construction_1():
    occ_event = occ_events.OccEvent()
    assert isinstance(occ_event, occ_events.OccEvent)
    assert occ_event.size() == 0


def test_OccEvent_construction_2(fcc_1NN_A_Va_event):
    prim, occ_event = fcc_1NN_A_Va_event

    assert isinstance(occ_event, occ_events.OccEvent)
    assert occ_event.size() == 2


def test_OccEvent_reverse(fcc_1NN_A_Va_event):
    prim, occ_event = fcc_1NN_A_Va_event

    assert occ_event.initial_occupation() == [0, 2]
    assert occ_event.final_occupation() == [2, 0]

    reversed_occ_event = occ_event.copy_reverse()
    assert reversed_occ_event.initial_occupation() == [2, 0]
    assert reversed_occ_event.final_occupation() == [0, 2]

    occ_event.reverse()
    assert occ_event.initial_occupation() == [2, 0]
    assert occ_event.final_occupation() == [0, 2]

    assert reversed_occ_event == occ_event


def test_OccEvent_to_from_dict():
    data = {
        "trajectories": [
            [
                {
                    # "chemical_name": "A",
                    "coordinate": [0, 0, 0, 0],
                    "occupant_index": 0,
                },
                {
                    # "chemical_name": "A",
                    "coordinate": [0, 0, 0, 1],
                    "occupant_index": 0,
                },
            ],
            [
                {
                    # "chemical_name": "Va",
                    "coordinate": [0, 0, 0, 1],
                    "occupant_index": 2,
                },
                {
                    # "chemical_name": "Va",
                    "coordinate": [0, 0, 0, 0],
                    "occupant_index": 2,
                },
            ],
        ]
    }

    prim = xtal_prims.FCC(a=1.0, occ_dof=["A", "B", "Va"])
    system = occ_events.OccSystem(prim)
    occ_event = occ_events.OccEvent.from_dict(data, system)

    assert isinstance(occ_event, occ_events.OccEvent)


def test_make_prim_periodic_orbit_1(fcc_1NN_A_Va_event):
    xtal_prim, occ_event = fcc_1NN_A_Va_event
    fg = xtal.make_factor_group(xtal_prim)
    occevent_symgroup_rep = occ_events.make_occevent_symgroup_rep(fg, xtal_prim)
    orbit = occ_events.make_prim_periodic_orbit(occ_event, occevent_symgroup_rep)

    assert isinstance(orbit, list)
    assert len(orbit) == 6
    for x in orbit:
        assert orbit[0] <= x


def test_make_prototype_occevent(fcc_1NN_A_Va_event):
    xtal_prim, occ_event = fcc_1NN_A_Va_event
    prototype = occ_events.make_canonical_occevent(xtal_prim, occ_event)
    assert isinstance(prototype, occ_events.OccEvent)


def test_save_load_occevent(fcc_1NN_A_Va_event, tmpdir):
    xtal_prim, occ_event = fcc_1NN_A_Va_event
    system = occ_events.OccSystem(xtal_prim)

    occ_events.save_occevent(tmpdir, "name", occ_event, system)
    loaded_occ_event = occ_events.load_occevent(tmpdir, "name", system)
    assert loaded_occ_event == occ_event


def test_get_occevent_coordinate(shared_datadir):
    import json

    import numpy as np

    import libcasm.clusterography as casmclust
    import libcasm.configuration as casmconfig

    with open(shared_datadir / "ZrO" / "equivalents_info.json", "r") as f:
        equivalents_info = json.load(f)
    prim = casmconfig.Prim.from_dict(equivalents_info["prototype"]["prim"])
    xtal_prim = prim.xtal_prim

    occ_system = occ_events.OccSystem(xtal_prim)
    with open(shared_datadir / "ZrO" / "event.json", "r") as f:
        event = occ_events.OccEvent.from_dict(json.load(f), system=occ_system)

    (
        phenomenal_clusters,
        equivalent_generating_op_indices,
    ) = casmclust.equivalents_info_from_dict(
        data=equivalents_info,
        xtal_prim=xtal_prim,
    )

    assert len(phenomenal_clusters) == 2
    assert equivalent_generating_op_indices == [0, 2]

    ### Generate phenomenal OccEvent consistent with local clexulator ###
    occevent_symgroup_rep = occ_events.make_occevent_symgroup_rep(
        prim.factor_group.elements, xtal_prim
    )
    occevent_orbit = []
    for i, generating_op_index in enumerate(equivalent_generating_op_indices):
        tmp = (occevent_symgroup_rep[generating_op_index] * event).standardize()
        trans = (
            phenomenal_clusters[i].sorted()[0].unitcell()
            - tmp.cluster().sorted()[0].unitcell()
        )
        occevent_orbit.append(tmp + trans)

    # Check that the clusters are consistent with equivalents_info, up to a translation
    assert len(occevent_orbit) == 2
    for i, equiv_occevent in enumerate(occevent_orbit):
        assert equiv_occevent.cluster().sorted() == phenomenal_clusters[i].sorted()

    # Make a supercell, get OccEvent coordinates within the supercell
    supercell = casmconfig.Supercell(
        prim=prim,
        transformation_matrix_to_super=np.eye(3, dtype="int") * 5,
    )

    # Check OccEvent in orbit
    unitcell_index, equivalent_index = occ_events.get_occevent_coordinate(
        occ_event=occevent_orbit[0],
        phenomenal_occevent=occevent_orbit,
        unitcell_index_converter=supercell.unitcell_index_converter,
    )
    assert unitcell_index == 0
    assert equivalent_index == 0

    # Check OccEvent in orbit
    unitcell_index, equivalent_index = occ_events.get_occevent_coordinate(
        occ_event=occevent_orbit[1],
        phenomenal_occevent=occevent_orbit,
        unitcell_index_converter=supercell.unitcell_index_converter,
    )
    assert unitcell_index == 0
    assert equivalent_index == 1

    for k in range(5):
        for j in range(5):
            for i in range(5):
                trans = np.array([i, j, k], dtype="int")
                unitcell_index = (
                    supercell.unitcell_index_converter.linear_unitcell_index(trans)
                )
                (
                    e_unitcell_index,
                    e_equivalent_index,
                ) = occ_events.get_occevent_coordinate(
                    occ_event=occevent_orbit[0] + trans,
                    phenomenal_occevent=occevent_orbit,
                    unitcell_index_converter=supercell.unitcell_index_converter,
                )
                print([i, j, k], unitcell_index, e_unitcell_index, e_equivalent_index)
                assert unitcell_index == e_unitcell_index

    # Check a translated OccEvent
    unitcell_index, equivalent_index = occ_events.get_occevent_coordinate(
        occ_event=occevent_orbit[0] + np.array([1, 0, 0], dtype="int"),
        phenomenal_occevent=occevent_orbit,
        unitcell_index_converter=supercell.unitcell_index_converter,
    )
    assert unitcell_index == supercell.unitcell_index_converter.linear_unitcell_index(
        np.array([1, 0, 0], dtype="int")
    )
    assert equivalent_index == 0

    # Check a translated OccEvent equivalent
    unitcell_index, equivalent_index = occ_events.get_occevent_coordinate(
        occ_event=occevent_orbit[1] + np.array([2, 2, 2]),
        phenomenal_occevent=occevent_orbit,
        unitcell_index_converter=supercell.unitcell_index_converter,
    )
    assert unitcell_index == supercell.unitcell_index_converter.linear_unitcell_index(
        np.array([2, 2, 2], dtype="int")
    )
    assert equivalent_index == 1

    # Check a translated OccEvent, where periodicity is used
    unitcell_index, equivalent_index = occ_events.get_occevent_coordinate(
        occ_event=occevent_orbit[0] + np.array([7, 2, 2]),
        phenomenal_occevent=occevent_orbit,
        unitcell_index_converter=supercell.unitcell_index_converter,
    )
    assert unitcell_index == supercell.unitcell_index_converter.linear_unitcell_index(
        np.array([2, 2, 2], dtype="int")
    )
    assert equivalent_index == 0

    # Check a translated OccEvent, where periodicity is used
    unitcell_index, equivalent_index = occ_events.get_occevent_coordinate(
        occ_event=occevent_orbit[1] + np.array([7, 2, 2]),
        phenomenal_occevent=occevent_orbit,
        unitcell_index_converter=supercell.unitcell_index_converter,
    )
    assert unitcell_index == supercell.unitcell_index_converter.linear_unitcell_index(
        np.array([2, 2, 2], dtype="int")
    )
    assert equivalent_index == 1

    # Check a transformed OccEvent: equivalent 0 -> equivalent 1
    unitcell_index, equivalent_index = occ_events.get_occevent_coordinate(
        occ_event=(occevent_symgroup_rep[2] * occevent_orbit[0]).standardize(),
        phenomenal_occevent=occevent_orbit,
        unitcell_index_converter=supercell.unitcell_index_converter,
    )
    assert equivalent_index == 1


def test_OccEvent_repr(fcc_1NN_A_Va_event):
    prim, occ_event = fcc_1NN_A_Va_event

    import io
    from contextlib import redirect_stdout

    f = io.StringIO()
    with redirect_stdout(f):
        print(occ_event)
    out = f.getvalue()
    assert "trajectories" in out

    print(occ_event)

    # Also test OccPosition.__repr__
    for traj in occ_event.trajectories():
        for pos in traj:
            f = io.StringIO()
            with redirect_stdout(f):
                print(pos)
            out = f.getvalue()
            assert "coordinate" in out
