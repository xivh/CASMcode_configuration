import math

import libcasm.configuration as casmconfig
import libcasm.enumerate as casmenum
import libcasm.local_configuration as casmlocal
import libcasm.occ_events as occ_events


def expected_neb_supercells():
    return {
        "SCEL1_1_1_1_0_0_0/0": {
            "supercell_name": "SCEL108_6_6_3_3_3_0",
            "transformation_matrix_to_supercell": [[-3, 3, 3], [3, -3, 3], [3, 3, -3]],
        },
        "SCEL2_2_1_1_0_0_1/0": {
            "supercell_name": "SCEL108_6_6_3_3_3_0",
            "transformation_matrix_to_supercell": [[-3, 3, 3], [3, -3, 3], [3, 3, -3]],
        },
        "SCEL2_2_1_1_0_1_1/0": {
            "supercell_name": "SCEL128_8_4_4_0_4_4",
            "transformation_matrix_to_supercell": [[4, 0, 4], [0, 4, 4], [-4, -4, 0]],
        },
        "SCEL3_3_1_1_0_0_2/0": {
            "supercell_name": "SCEL108_6_6_3_3_3_0",
            "transformation_matrix_to_supercell": [[-3, 3, 3], [3, -3, 3], [3, 3, -3]],
        },
        "SCEL3_3_1_1_0_2_1/0": {
            "supercell_name": "SCEL108_6_6_3_3_3_0",
            "transformation_matrix_to_supercell": [[-3, 3, 3], [3, -3, 3], [3, 3, -3]],
        },
        "SCEL3_3_1_1_0_2_2/0": {
            "supercell_name": "SCEL108_6_6_3_3_3_0",
            "transformation_matrix_to_supercell": [[-3, 3, 3], [3, -3, 3], [3, 3, -3]],
        },
        "SCEL4_2_2_1_0_1_0/0": {
            "supercell_name": "SCEL128_8_4_4_0_4_0",
            "transformation_matrix_to_supercell": [[0, -4, 4], [4, 0, -4], [0, 4, 4]],
        },
        "SCEL4_2_2_1_1_1_0/0": {
            "supercell_name": "SCEL108_6_6_3_3_3_0",
            "transformation_matrix_to_supercell": [[-3, 3, 3], [3, -3, 3], [3, 3, -3]],
        },
        "SCEL4_4_1_1_0_0_0/0": {
            "supercell_name": "SCEL128_8_4_4_0_0_0",
            "transformation_matrix_to_supercell": [[0, 0, 8], [0, -4, -4], [4, 4, 0]],
        },
        "SCEL4_4_1_1_0_0_2/0": {
            "supercell_name": "SCEL128_8_4_4_0_0_4",
            "transformation_matrix_to_supercell": [[0, -4, -4], [0, 4, -4], [4, 0, 4]],
        },
        "SCEL4_4_1_1_0_0_3/0": {
            "supercell_name": "SCEL128_8_4_4_0_0_4",
            "transformation_matrix_to_supercell": [[0, -4, -4], [0, 4, -4], [4, 0, 4]],
        },
        "SCEL4_4_1_1_0_1_0/0": {
            "supercell_name": "SCEL128_8_4_4_0_4_0",
            "transformation_matrix_to_supercell": [[0, -4, 4], [4, 0, -4], [0, 4, 4]],
        },
        "SCEL4_4_1_1_0_2_1/0": {
            "supercell_name": "SCEL108_12_3_3_0_6_3",
            "transformation_matrix_to_supercell": [[-3, 3, 3], [3, -3, 3], [3, 3, -6]],
        },
        "SCEL5_1_1_5_0_0_0/0": {
            "supercell_name": "SCEL125_5_5_5_0_0_0",
            "transformation_matrix_to_supercell": [[5, 0, 0], [0, 5, 0], [0, 0, 5]],
        },
        "SCEL5_5_1_1_0_0_3/0": {
            "supercell_name": "SCEL125_5_5_5_0_0_0",
            "transformation_matrix_to_supercell": [[5, 0, 0], [0, 5, 0], [0, 0, 5]],
        },
        "SCEL5_5_1_1_0_0_4/0": {
            "supercell_name": "SCEL125_5_5_5_0_0_0",
            "transformation_matrix_to_supercell": [[5, 0, 0], [0, 5, 0], [0, 0, 5]],
        },
        "SCEL5_5_1_1_0_1_3/0": {
            "supercell_name": "SCEL125_5_5_5_0_0_0",
            "transformation_matrix_to_supercell": [[5, 0, 0], [0, 5, 0], [0, 0, 5]],
        },
        "SCEL5_5_1_1_0_4_3/0": {
            "supercell_name": "SCEL125_5_5_5_0_0_0",
            "transformation_matrix_to_supercell": [[5, 0, 0], [0, 5, 0], [0, 0, 5]],
        },
    }


def test_make_supercells_for_point_defects_1(fcc_1NN_A_Va_event):
    xtal_prim, event = fcc_1NN_A_Va_event
    r = 1.0
    a = 4 * r / math.sqrt(2)
    prim = casmconfig.Prim(xtal_prim)

    # Make event prim info
    system = occ_events.OccSystem(xtal_prim=prim.xtal_prim)
    event_prim_info = casmlocal.OccEventPrimSymInfo(
        prim=prim,
        system=system,
        prototype_event=event,
    )
    # print("# equivalent events:", len(event_prim_info.events))
    assert len(event_prim_info.events) == 6
    local_orbits = event_prim_info.make_local_orbits(
        max_length=[0, 0, 3.01 * a],
        cutoff_radius=[0, 1.01 * a],
    )
    # print("# orbits:", len(local_orbits[0]))
    assert len(local_orbits[0]) == 7

    cluster_sizes = [len(orbit[0]) for orbit in local_orbits[0]]
    # print(cluster_sizes)
    assert cluster_sizes == [0, 1, 1, 1, 1, 1, 1]

    orbit_sizes = [len(orbit) for orbit in local_orbits[0]]
    assert orbit_sizes == [1, 4, 4, 8, 4, 2, 4]

    # print("Orbit sizes:")
    # for i_orbit, orbit in enumerate(local_orbits[0]):
    #     orbit_size = len(orbit)
    #     cluster_size = 0
    #     if orbit_size > 0:
    #         cluster_size = len(orbit[0])
    #     print(
    #         f"- orbit: {i_orbit}, "
    #         f"cluster size: {cluster_size}, "
    #         f"orbit size: {orbit_size}"
    #     )
    # print()

    required_sites = casmenum.make_required_sites(
        phenomenal_clusters=event_prim_info.phenomenal_clusters,
        local_orbits=local_orbits,
    )
    # print("# of required sites:", len(required_sites[0]))
    assert len(required_sites[0]) == 28

    # Make configurations
    supercells = casmconfig.SupercellSet(prim=prim)
    configs = casmconfig.ConfigurationSet()
    config_enum = casmenum.ConfigEnumAllOccupations(
        prim=prim,
        supercell_set=supercells,
    )
    for i, configuration in enumerate(config_enum.by_supercell(max=5)):
        configs.add(configuration)
    assert len(supercells) == 18
    assert len(configs) == 291

    # Make supercells for point defect calculations
    neb_supercells = {}
    last_supercell = None
    for record in configs:
        if record.configuration.supercell == last_supercell:
            continue
        # print("~~~Configuration: ~~~")
        # print(record)
        # print()
        candidate_supercells = casmenum.make_supercells_for_point_defects(
            motif=record.configuration,
            base_max_volume=10,
            min_volume=1,
            max_volume=8**3,
            base_min_factor_group_size=None,
            required_sites=required_sites,
        )

        # casmenum.plot_point_defect_supercell_scores(
        #     candidate_supercells,
        #     record.configuration_name,
        #     min_volume=100,
        #     max_volume=216,
        #     require_has_required_sites=True,
        #     require_has_all_motif_operations=True,
        #     min_factor_group_size=None,
        #     min_voronoi_inner_radius=None,
        # )

        optimal_supercells = casmenum.find_optimal_point_defect_supercells(
            candidate_supercells=candidate_supercells,
            min_volume=100,
            max_volume=216,
            require_has_required_sites=True,
            require_has_all_motif_operations=True,
            min_factor_group_size=None,
            min_voronoi_inner_radius=None,
        )
        # print_optimal_point_defect_supercells(optimal_supercells)
        neb_supercells[record.configuration_name] = optimal_supercells[0][0].to_dict()

        last_supercell = record.configuration.supercell

        if len(optimal_supercells) == 0:
            raise ValueError(
                "No optimal supercells that match all criteria were found "
                "for point defect calculations "
                f"in configuration {record.configuration_name}."
            )

        # print(f"~~~ NEB supercells for {record.configuration_name} ~~~")
        # print()
        #
        # for i, (supercell, data) in enumerate(optimal_supercells):
        #     print(f"~~~ Supercell choice: {i} ~~~")
        #     casmenum.print_point_defect_supercell_info(supercell, data)
        #     print()

    import libcasm.xtal as xtal

    print("NEB supercells:")
    print(xtal.pretty_json(neb_supercells))
    for key, value in neb_supercells.items():
        print(f"configuration name: {key}")
        config = configs.get_by_name(key)
        print(f"configuration: {config}")
        print(f"supercell: {value}")
        print()
    assert neb_supercells == expected_neb_supercells()
