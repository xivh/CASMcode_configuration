import numpy as np

import libcasm.configuration as config
import libcasm.enumerate as enum
import libcasm.local_configuration as casmlocal
import libcasm.occ_events as occ_events
import libcasm.xtal as xtal


def test_make_occevent_suborbits(fcc_1NN_A_Va_event):
    # setup: FCC prim, 1NN A-Va exchange event
    xtal_prim, occ_event = fcc_1NN_A_Va_event
    prim = config.Prim(xtal_prim)
    occ_events.OccSystem(xtal_prim)

    def make_occevent_orbit(prim, occ_event):
        prim_factor_group = prim.factor_group
        prim_rep = occ_events.make_occevent_symgroup_rep(
            prim_factor_group.elements, prim.xtal_prim
        )
        return occ_events.make_prim_periodic_orbit(occ_event, prim_rep)

    orbit = make_occevent_orbit(prim, occ_event)
    assert len(orbit) == 6

    # conventional 4-site FCC
    T_motif = np.array(
        [
            [-1, 1, 1],
            [1, -1, 1],
            [1, 1, -1],
        ]
    )
    motif_supercell = config.Supercell(prim, T_motif)
    config.Configuration(motif_supercell)

    # supercell to fill: 3x3x3 of the conventional FCC
    supercell_1 = config.Supercell(prim, T_motif * 3)
    assert len(supercell_1.factor_group.elements) == 48

    suborbits_1 = casmlocal.make_occevent_suborbits(supercell_1, occ_event)
    assert len(suborbits_1) == 1
    assert len(suborbits_1[0]) == 6

    # supercell to fill: 2x2x3 of the conventional FCC
    S = np.array(
        [
            [2, 0, 0],
            [0, 2, 0],
            [0, 0, 3],
        ]
    )
    supercell_2 = config.Supercell(prim, T_motif @ S)
    assert len(supercell_2.factor_group.elements) == 16

    suborbits_2 = casmlocal.make_occevent_suborbits(supercell_2, occ_event)
    assert len(suborbits_2) == 2
    assert len(suborbits_2[0]) == 2
    assert len(suborbits_2[1]) == 4


def test_make_all_distinct_local_perturbations(fcc_1NN_A_Va_event):
    # setup: FCC prim, 1NN A-Va exchange event
    xtal_prim, phenomenal_occ_event = fcc_1NN_A_Va_event
    prim = config.Prim(xtal_prim)
    system = occ_events.OccSystem(xtal_prim)

    # make local-clusters to perturb to generate local environments:
    # null and point clusters only, 1NN sites to the event sites
    # expect null orbit and 4 point cluster orbits
    cluster_specs = occ_events.make_occevent_cluster_specs(
        xtal_prim=xtal_prim,
        phenomenal_occ_event=phenomenal_occ_event,
        max_length=[0.0, 0.0],
        cutoff_radius=[0.0, 2.01],
    )
    orbits = cluster_specs.make_orbits()
    local_clusters = [orbit[0] for orbit in orbits]

    # conventional 4-site FCC
    T_motif = np.array(
        [
            [-1, 1, 1],
            [1, -1, 1],
            [1, 1, -1],
        ]
    )
    motif_supercell = config.Supercell(prim, T_motif)
    motif = config.Configuration(motif_supercell)

    # supercell to fill: 3x3x3 of the conventional FCC
    supercell = config.Supercell(prim, T_motif * 3)

    # generate unique perturbation configurations
    configurations = enum.make_all_distinct_local_perturbations(
        supercell, phenomenal_occ_event, motif, local_clusters
    )

    # check the results:
    # expect 9 results:
    # - 1 is the unperturbed background
    # - + 4*2: for B or Va on the point cluster sites
    for i, c in enumerate(configurations):
        print(i)

        # show occupation count
        ssum = {0: 0, 1: 0, 2: 0}
        for s in c.occupation.tolist():
            ssum[s] += 1
        print("ssum:", ssum)

        # get interpolated structures
        structures = enum.make_occevent_simple_structures(
            c, phenomenal_occ_event, [0.0, 0.5, 1.0], system, False
        )
        print(xtal.pretty_json(structures[0].to_dict()))
        print()

    assert len(configurations) == 9
