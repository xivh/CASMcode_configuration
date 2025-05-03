import sys

import numpy as np

import libcasm.configuration as casmconfig
import libcasm.enumerate as casmenum
import libcasm.local_configuration as casmlocal
import libcasm.occ_events as occ_events


def make_test_1_config(supercells):
    data = {
        "basis": "standard",
        "dof": {"occ": [0]},
        "supercell_name": "SCEL1_1_1_1_0_0_0",
        "transformation_matrix_to_supercell": [[1, 0, 0], [0, 1, 0], [0, 0, 1]],
    }
    # data = {
    #     "basis": "standard",
    #     "dof": {"occ": [1, 0]},
    #     "supercell_name": "SCEL2_2_1_1_0_1_1",
    #     "transformation_matrix_to_supercell": [[1, 0, 1], [0, 1, 1], [-1, -1, 0]],
    # }
    # data = {
    #     "basis": "standard",
    #     "dof": {"occ": [1, 0, 1, 0, 0, 0, 0, 0]},
    #     "transformation_matrix_to_supercell": [
    #         [-2, 1, 1],
    #         [2, -1, 1],
    #         [2, 1, -1],
    #     ],
    # }
    return casmconfig.Configuration.from_dict(data=data, supercells=supercells)


def make_test_1_neb_supercell(supercells):
    # data = {
    #     "supercell_name": "SCEL125_5_5_5_0_0_0",
    #     "transformation_matrix_to_supercell": [[5, 0, 0], [0, 5, 0], [0, 0, 5]],
    # }
    # data = {
    #     "supercell_name": "SCEL125_5_5_5_0_0_0",
    #     "transformation_matrix_to_supercell": [[3, 0, 0], [0, 4, 0], [0, 0, 4]],
    # }
    # data = {
    #     "supercell_name": "SCEL128_8_4_4_0_0_4",
    #     "transformation_matrix_to_supercell": [[0, -4, -4], [0, 4, -4], [4, 0, 4]],
    # }

    x = 4
    data = {
        "transformation_matrix_to_supercell": [
            [-x, x, x],
            [x, -x, x],
            [x, x, -x],
        ],
    }
    return casmconfig.Supercell.from_dict(data=data, supercells=supercells)


def test_make_distinct_local_configurations_L12(fcc_1NN_A_Va_event):
    xtal_prim, event = fcc_1NN_A_Va_event
    # r = 1.0
    # a = 4 * r / math.sqrt(2)
    prim = casmconfig.Prim(xtal_prim)

    T_fcc_conventional = np.array(
        [
            [-1, 1, 1],
            [1, -1, 1],
            [1, 1, -1],
        ],
        dtype="int",
    )
    # L12 ordering (A3, B1)
    config = casmconfig.Configuration(
        casmconfig.Supercell(
            prim=prim,
            transformation_matrix_to_super=T_fcc_conventional,
        )
    )
    config.set_occ(0, 1)

    system = occ_events.OccSystem(xtal_prim=prim.xtal_prim)
    event_info = casmlocal.OccEventSymInfo.init(
        prim=prim,
        system=system,
        prototype_event=event,
    )

    local_configurations = casmenum.make_distinct_local_configurations(
        background=config,
        event_info=event_info,
        fix="event",
    )

    print("Fix event:")
    print("Total # Local configurations:", len(local_configurations))
    for x in local_configurations:
        scelname = casmconfig.SupercellRecord(x.configuration.supercell).supercell_name
        print("-", scelname, x.configuration.occupation.tolist(), x.pos)
    assert len(local_configurations) == 2

    local_configurations = casmenum.make_distinct_local_configurations(
        background=config,
        event_info=event_info,
        fix="background",
    )

    print("Fix background:")
    print("Total # Local configurations:", len(local_configurations))
    for x in local_configurations:
        scelname = casmconfig.SupercellRecord(x.configuration.supercell).supercell_name
        print("-", scelname, x.configuration.occupation.tolist(), x.pos)
    assert len(local_configurations) == 2


def test_make_distinct_local_configurations_fcc(fcc_1NN_A_Va_event):
    xtal_prim, event = fcc_1NN_A_Va_event
    # r = 1.0
    # a = 4 * r / math.sqrt(2)
    prim = casmconfig.Prim(xtal_prim)

    system = occ_events.OccSystem(xtal_prim=prim.xtal_prim)
    event_info = casmlocal.OccEventSymInfo.init(
        prim=prim,
        system=system,
        prototype_event=event,
    )

    supercell_set = casmconfig.SupercellSet(prim=prim)

    print()
    config_enum = casmenum.ConfigEnumAllOccupations(
        prim=prim,
        supercell_set=supercell_set,
    )
    event_supercell_info = dict()
    for config in config_enum.by_supercell(max=5):
        if np.sum(config.occupation == 2) > 0:
            continue

        local_configurations_fix_background = (
            casmenum.make_distinct_local_configurations(
                background=config,
                event_info=event_info,
                fix="background",
            )
        )

        local_configurations_fix_event = casmenum.make_distinct_local_configurations(
            background=config,
            event_info=event_info,
            fix="event",
        )

        assert len(local_configurations_fix_background) == len(
            local_configurations_fix_event
        )
        n_local_configs = len(local_configurations_fix_background)

        scelname = casmconfig.SupercellRecord(config.supercell).supercell_name
        print(f"{scelname}, {config.occupation}: {n_local_configs}")

        def _make_canonical(x):
            if x.configuration.supercell not in event_supercell_info:
                event_supercell_info[x.configuration.supercell] = (
                    casmlocal.OccEventSupercellSymInfo(
                        event_prim_info=event_info.event_prim_info,
                        supercell=x.configuration.supercell,
                    )
                )
            curr = event_supercell_info[x.configuration.supercell]
            (
                config,
                pos,
                final_event_supercell_info,
            ) = curr.make_canonical_local_configuration(
                configuration=x.configuration,
                pos=x.pos,
                in_canonical_pos=True,
                in_canonical_supercell=True,
                apply_event_occupation=False,
            )
            return casmlocal.LocalConfiguration(
                pos=pos,
                configuration=config,
                event_info=x.event_info,
            )

        canonical_fix_background = []
        for x in local_configurations_fix_background:
            canonical_fix_background.append(_make_canonical(x))
        canonical_fix_background.sort()

        canonical_fix_event = []
        for x in local_configurations_fix_event:
            canonical_fix_event.append(_make_canonical(x))
        canonical_fix_event.sort()

        def _print(_local):
            for i, x in enumerate(_local):
                config = x.configuration
                scelname = casmconfig.SupercellRecord(config.supercell).supercell_name
                if config.supercell not in event_supercell_info:
                    event_supercell_info[config.supercell] = (
                        casmlocal.OccEventSupercellSymInfo(
                            event_prim_info=event_info.event_prim_info,
                            supercell=config.supercell,
                        )
                    )
                _curr = event_supercell_info[config.supercell]
                print(scelname, config.occupation.tolist(), x.pos)
                if x.pos[0] == 0:
                    event_group = _curr.event_group_rep(x.pos)
                    for op in event_group:
                        print(
                            f"- ({op.prim_factor_group_index()}):",
                            (op * config).occupation.tolist(),
                        )
            print()

        # print("local_configurations_fix_background:")
        # _print(local_configurations_fix_background)
        #
        # print("canonical_fix_background:")
        # _print(canonical_fix_background)
        #
        # print("local_configurations_fix_event:")
        # _print(local_configurations_fix_event)
        #
        # print("canonical_fix_event:")
        # _print(canonical_fix_event)

        assert len(canonical_fix_background) == len(local_configurations_fix_background)
        assert len(canonical_fix_event) == len(local_configurations_fix_event)

        assert canonical_fix_background == canonical_fix_event


def _print_result(result, existing: bool):
    print(
        f"- initial={result.i_initial}, "
        f"orbit={result.i_local_orbit}, "
        f"sites={result.sites}, "
        f"sublattices={result.sublattices}, "
        f"asymmetric_units={result.asymmetric_units}, "
        f"initial={result.initial_occupation}, ",
        f"final={result.final_occupation}",
        end="",
    )
    if existing:
        print(" (existing)")
    else:
        print()
    sys.stdout.flush()


def test_ConfigEnumLocalOccupations_L12(fcc_1NN_A_Va_event_L12):
    prim, event, event_info, config = fcc_1NN_A_Va_event_L12
    r = 1.0
    # a = 4 * r / math.sqrt(2)

    supercells = casmconfig.SupercellSet(prim=prim)
    local_config_list = casmlocal.LocalConfigurationList(event_info=event_info)

    # NEB supercell
    neb_supercell = make_test_1_neb_supercell(supercells)

    # Local-cluster specs
    cluster_specs = casmenum.make_first_n_orbits_cluster_specs(
        prim=prim,
        phenomenal=event,
        cutoff_radius=[0, 2.01 * r, 2.01 * r],
        # cutoff_radius=[0],
        make_all_possible_orbits=True,
    )

    verbose = False

    config_enum = casmenum.ConfigEnumLocalOccupations(
        event_info=event_info,
        supercell_set=supercells,
        verbose=verbose,
    )

    for result in config_enum.by_cluster_specs(
        background=config,
        supercell=neb_supercell,
        cluster_specs=cluster_specs,
        fix="background",
        # fix="both",
        # pos=(0, 0),
        neighborhood_from_orbits=[1],
    ):
        # Skip configurations with vacancies
        if 2 in result.final_occupation:
            continue

        #
        if result.canonical_local_configuration not in local_config_list:
            if verbose:
                _print_result(result, existing=False)
                # print(result)
            local_config_list.append(result.canonical_local_configuration)
        elif verbose:
            _print_result(result, existing=True)
            # print(result)

    if verbose:
        print()
        print("Reference:")
        print(config_enum.reference)
        print()
        print("Total # Local configurations:", len(local_config_list))

    assert len(local_config_list) == 17


# Skip configurations with vacancies
def skip_results_with_vacancies(result):
    """Valid for the FCC A-B-Va prim"""
    return 2 in result.final_occupation


def set_conventional_fcc_supercell(test, length: int):
    x = length
    test.set_supercell(
        transformation_matrix_to_super=np.array(
            [
                [-x, x, x],
                [x, -x, x],
                [x, x, -x],
            ],
            dtype="int",
        )
    )


def test_ConfigEnumLocalOccupations_L12_v2(L12_1NN_A_Va_local_enum):
    test = L12_1NN_A_Va_local_enum

    # Set supercell
    set_conventional_fcc_supercell(test, length=4)

    # Set local cluster specs
    r = 1.0
    test.make_local_cluster_specs(
        cutoff_radius=[0, 2.01 * r],
        make_all_possible_orbits=True,
    )

    # Run test - fix background
    test.check_enum(
        verbose=True,
        fix="background",
        neighborhood_from_orbits=[1],
        skip_result_f=skip_results_with_vacancies,
        n_expected=17,
    )
    first = test.fully_canonical_local_config_list.copy()

    # Run test - fix event
    test.check_enum(
        verbose=True,
        fix="event",
        neighborhood_from_orbits=[1],
        skip_result_f=skip_results_with_vacancies,
        n_expected=17,
    )
    second = test.fully_canonical_local_config_list.copy()

    for x in first:
        assert x in second
