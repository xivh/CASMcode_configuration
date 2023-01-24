import numpy as np
import libcasm.configuration as config
import libcasm.occ_events as occ_events
import libcasm.enumerate as enum
import libcasm.xtal as xtal


def test_make_all_distinct_local_perturbations(fcc_1NN_A_Va_event):

    xtal_prim, phenomenal_occ_event = fcc_1NN_A_Va_event
    prim = config.Prim(xtal_prim)
    system = occ_events.OccSystem(xtal_prim)
    cluster_specs = occ_events.make_occevent_cluster_specs(
        xtal_prim=xtal_prim,
        phenomenal_occ_event=phenomenal_occ_event,
        max_length=[0.0, 0.0],
        cutoff_radius=[0.0, 2.01])
    orbits = cluster_specs.make_orbits()
    local_clusters = [orbit[0] for orbit in orbits]

    T_motif = np.array([
        [-1, 1, 1],
        [1, -1, 1],
        [1, 1, -1],
    ])
    motif_supercell = config.Supercell(prim, T_motif)
    S_motif = motif_supercell.superlattice().column_vector_matrix()
    print(S_motif)
    motif = config.Configuration(motif_supercell)

    T = T_motif * 3
    supercell = config.Supercell(prim, T)
    S = supercell.superlattice().column_vector_matrix()
    print(S)

    configurations = enum.make_all_distinct_local_perturbations(
        supercell, phenomenal_occ_event, motif, local_clusters)

    for i, c in enumerate(configurations):
        print(i)

        # show occupation count
        ssum = {0: 0, 1: 0, 2: 0}
        for s in c.occupation().tolist():
            ssum[s] += 1
        print("ssum:", ssum)

        # get interpolated structures
        structures = enum.make_occ_event_simple_structures(
            c, phenomenal_occ_event, [0.0, 0.5, 1.0], system, False)
        print(xtal.pretty_json(structures[0].to_dict()))
        print()

    assert len(configurations) == 9
