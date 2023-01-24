import libcasm.xtal as xtal
import libcasm.xtal.prims as xtal_prims
from libcasm.occ_events import OccEvent, OccSystem, make_occevent_symgroup_rep, make_prim_periodic_orbit


def test_OccEvent_construction_1():
    occ_event = OccEvent()
    assert isinstance(occ_event, OccEvent)
    assert occ_event.size() == 0


def test_OccEvent_construction_2(fcc_1NN_A_Va_event):
    prim, occ_event = fcc_1NN_A_Va_event

    assert isinstance(occ_event, OccEvent)
    assert occ_event.size() == 2


def test_OccEvent_to_from_dict():
    data = {
        "trajectories": [
            [
                {
                    #"chemical_name": "A",
                    "coordinate": [0, 0, 0, 0],
                    "occupant_index": 0
                },
                {
                    #"chemical_name": "A",
                    "coordinate": [0, 0, 0, 1],
                    "occupant_index": 0
                }
            ],
            [
                {
                    #"chemical_name": "Va",
                    "coordinate": [0, 0, 0, 1],
                    "occupant_index": 2
                },
                {
                    #"chemical_name": "Va",
                    "coordinate": [0, 0, 0, 0],
                    "occupant_index": 2
                }
            ]
        ]
    }

    prim = xtal_prims.FCC(a=1.0, occ_dof=["A", "B", "Va"])
    system = OccSystem(prim)
    occ_event = OccEvent.from_dict(data, system)

    assert isinstance(occ_event, OccEvent)


def test_make_prim_periodic_orbit_1(fcc_1NN_A_Va_event):
    xtal_prim, occ_event = fcc_1NN_A_Va_event
    fg = xtal.make_factor_group(xtal_prim)
    occevent_symgroup_rep = make_occevent_symgroup_rep(fg, xtal_prim)
    orbit = make_prim_periodic_orbit(occ_event, occevent_symgroup_rep)

    assert isinstance(orbit, list)
    assert len(orbit) == 6
