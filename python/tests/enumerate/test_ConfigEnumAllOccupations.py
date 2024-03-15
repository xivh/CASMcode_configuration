import numpy as np

import libcasm.configuration as casmconfig
import libcasm.enumerate as casmenum
import libcasm.xtal as xtal
import libcasm.xtal.prims as xtal_prims


def test_ConfigEnumAllOccupations_by_supercell_FCC_1():
    xtal_prim = xtal_prims.FCC(
        r=0.5,
        occ_dof=["A", "B"],
    )
    prim = casmconfig.Prim(xtal_prim)
    supercell_set = casmconfig.SupercellSet(prim=prim)
    configuration_set = casmconfig.ConfigurationSet()
    config_enum = casmenum.ConfigEnumAllOccupations(
        prim=prim,
        supercell_set=supercell_set,
    )
    for i, configuration in enumerate(
        config_enum.by_supercell(
            supercells={
                "max": 4,
            },
        )
    ):
        configuration_set.add(configuration)
        assert isinstance(configuration, casmconfig.Configuration)

    # for i, record in enumerate(configuration_set):
    #     print(xtal.pretty_json(record.configuration.to_dict()))

    assert len(configuration_set) == 29


# TODO: test with fixed strain state and unique occupations


def test_ConfigEnumAllOccupations_by_supercell_FCC_2():
    xtal_prim = xtal_prims.FCC(
        r=0.5,
        occ_dof=["A", "B", "C"],
    )
    prim = casmconfig.Prim(xtal_prim)
    supercell_set = casmconfig.SupercellSet(prim=prim)
    configuration_set = casmconfig.ConfigurationSet()
    config_enum = casmenum.ConfigEnumAllOccupations(
        prim=prim,
        supercell_set=supercell_set,
    )
    for i, configuration in enumerate(
        config_enum.by_supercell(
            supercells={
                "max": 4,
            }
        )
    ):
        configuration_set.add(configuration)
        assert isinstance(configuration, casmconfig.Configuration)

    # for i, record in enumerate(configuration_set):
    #     print(xtal.pretty_json(record.configuration.to_dict()))

    assert len(configuration_set) == 126


def test_ConfigEnumAllOccupations_by_supercell_list_FCC_1():
    xtal_prim = xtal_prims.FCC(
        r=0.5,
        occ_dof=["A", "B"],
    )
    prim = casmconfig.Prim(xtal_prim)
    supercell_set = casmconfig.SupercellSet(prim=prim)
    scel_enum = casmenum.ScelEnum(
        prim=prim,
        supercell_set=supercell_set,
    )
    supercell_list = [x for x in scel_enum.by_volume(max=4)]

    configuration_set = casmconfig.ConfigurationSet()
    config_enum = casmenum.ConfigEnumAllOccupations(
        prim=prim,
        supercell_set=supercell_set,
    )
    for i, configuration in enumerate(
        config_enum.by_supercell_list(supercells=supercell_list)
    ):
        configuration_set.add(configuration)
        assert isinstance(configuration, casmconfig.Configuration)

    # for i, record in enumerate(configuration_set):
    #     print(xtal.pretty_json(record.configuration.to_dict()))

    assert len(configuration_set) == 29


def test_ConfigEnumAllOccupations_by_supercell_list_FCC_2():
    xtal_prim = xtal_prims.FCC(
        r=0.5,
        occ_dof=["A", "B", "C"],
    )
    prim = casmconfig.Prim(xtal_prim)
    supercell_set = casmconfig.SupercellSet(prim=prim)
    scel_enum = casmenum.ScelEnum(
        prim=prim,
        supercell_set=supercell_set,
    )
    supercell_list = [x for x in scel_enum.by_volume(max=4)]

    configuration_set = casmconfig.ConfigurationSet()
    config_enum = casmenum.ConfigEnumAllOccupations(
        prim=prim,
        supercell_set=supercell_set,
    )
    for i, configuration in enumerate(
        config_enum.by_supercell_list(supercells=supercell_list)
    ):
        configuration_set.add(configuration)
        assert isinstance(configuration, casmconfig.Configuration)

    # for i, record in enumerate(configuration_set):
    #     print(xtal.pretty_json(record.configuration.to_dict()))

    assert len(configuration_set) == 126


def test_ConfigEnumAllOccupations_by_sites_FCC_1():
    xtal_prim = xtal_prims.FCC(
        r=0.5,
        occ_dof=["A", "B"],
    )
    prim = casmconfig.Prim(xtal_prim)
    supercell_set = casmconfig.SupercellSet(prim=prim)
    supercell = casmconfig.Supercell(
        prim=prim,
        transformation_matrix_to_super=np.eye(3, dtype=int) * 2,
    )
    converter = supercell.site_index_converter
    background = casmconfig.Configuration(supercell)
    config_enum = casmenum.ConfigEnumAllOccupations(
        prim=prim,
        supercell_set=supercell_set,
    )

    # Test by_integral_site_coordinates:
    n_configs = []
    for i in range(supercell.n_sites):
        sites = []
        for j in range(i):
            sites.append(converter.integral_site_coordinate(j))

        generator = config_enum.by_integral_site_coordinates(
            background=background, sites=sites
        )

        n_configs.append(sum(1 for config in generator))

    assert n_configs == [0, 2, 4, 8, 16, 32, 64, 128]

    # Test by_linear_site_indices:
    n_configs = []
    for i in range(supercell.n_sites):
        sites = set(range(i))
        generator = config_enum.by_linear_site_indices(
            background=background, sites=sites
        )
        n_configs.append(sum(1 for config in generator))

    assert n_configs == [0, 2, 4, 8, 16, 32, 64, 128]


def test_ConfigEnumAllOccupations_by_cluster_FCC_1():
    xtal_prim = xtal_prims.FCC(
        a=1.0,
        occ_dof=["A", "B"],
    )
    print(xtal.pretty_json(xtal_prim.to_dict()))
    prim = casmconfig.Prim(xtal_prim)
    supercell_set = casmconfig.SupercellSet(prim=prim)
    configuration_set = casmconfig.ConfigurationSet()

    # conventional FCC cell
    supercell = casmconfig.Supercell(
        prim=prim,
        transformation_matrix_to_super=np.array(
            [
                [-1, 1, 1],
                [1, -1, 1],
                [1, 1, -1],
            ]
        ),
    )

    # L12 ordering (A3, B1)
    background = casmconfig.Configuration(supercell)
    background.set_occ(0, 1)

    config_enum = casmenum.ConfigEnumAllOccupations(
        prim=prim,
        supercell_set=supercell_set,
    )

    # 1x1x1, 2x2x2, 3x3x3 supercells
    supercells = {
        "max": 27,
        "diagonal_only": True,
        "fixed_shape": True,
        "unit_cell": supercell.transformation_matrix_to_super.tolist(),
    }

    # point clusters, 1NN & 2NN pair clusters
    cluster_specs = {
        "orbit_branch_specs": {
            "2": {"max_length": 1.01},
        },
    }

    info = casmenum.ConfigEnumInfo(config_enum, configuration_set)
    for configuration in config_enum.by_cluster(
        background=background,
        cluster_specs=cluster_specs,
        supercells=supercells,
    ):
        info.check()
        info.n_config_total += 1
        configuration_set.add(configuration)
    info.finish()

    assert len(configuration_set)
