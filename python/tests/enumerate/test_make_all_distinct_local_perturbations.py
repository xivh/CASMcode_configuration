import numpy as np
import libcasm.configuration as config
import libcasm.occ_events as occ_events
import libcasm.enumerate as enum
import libcasm.xtal as xtal


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
        cutoff_radius=[0.0, 2.01])
    orbits = cluster_specs.make_orbits()
    local_clusters = [orbit[0] for orbit in orbits]

    # conventional 4-site FCC
    T_motif = np.array([
        [-1, 1, 1],
        [1, -1, 1],
        [1, 1, -1],
    ])
    motif_supercell = config.Supercell(prim, T_motif)
    motif = config.Configuration(motif_supercell)

    # supercell to fill: 3x3x3 of the conventional FCC
    supercell = config.Supercell(prim, T_motif * 3)

    # generate unique perturbation configurations
    configurations = enum.make_all_distinct_local_perturbations(
        supercell, phenomenal_occ_event, motif, local_clusters)

    # check the results:
    # expect 9 results:
    # - 1 is the unperturbed background
    # - + 4*2: for B or Va on the point cluster sites
    for i, c in enumerate(configurations):
        print(i)

        # show occupation count
        ssum = {0: 0, 1: 0, 2: 0}
        for s in c.occupation().tolist():
            ssum[s] += 1
        print("ssum:", ssum)

        # get interpolated structures
        structures = enum.make_occevent_simple_structures(
            c, phenomenal_occ_event, [0.0, 0.5, 1.0], system, False)
        print(xtal.pretty_json(structures[0].to_dict()))
        print()

    assert len(configurations) == 9
