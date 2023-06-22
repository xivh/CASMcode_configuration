import libcasm.xtal as xtal
import libcasm.xtal.prims as xtal_prims
import libcasm.occ_events as occ_events
from libcasm.occ_events import (
    OccEvent,
    OccSystem,
    make_occevent_symgroup_rep,
    make_prim_periodic_orbit,
)


def test_OccEvent_construction_1():
    occ_event = OccEvent()
    assert isinstance(occ_event, occ_events.OccEvent)
    assert occ_event.size() == 0


def test_OccEvent_construction_2(fcc_1NN_A_Va_event):
    prim, occ_event = fcc_1NN_A_Va_event

    assert isinstance(occ_event, occ_events.OccEvent)
    assert occ_event.size() == 2


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
