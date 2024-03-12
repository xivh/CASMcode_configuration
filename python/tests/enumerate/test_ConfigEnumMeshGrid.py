import numpy as np

import libcasm.clexulator as casmclex
import libcasm.configuration as casmconfig
import libcasm.enumerate as casmenum
import libcasm.xtal as xtal
import libcasm.xtal.prims as xtal_prims


def test_ConfigEnumMeshGrid_by_range_FCC_1():
    """Test ConfigEnumMeshGrid.by_range

    Using: FCC, prim unit cell, default background, full basis, num=5,
    """
    xtal_prim = xtal_prims.FCC(
        r=0.5,
        occ_dof=["A", "B"],
        global_dof=[xtal.DoFSetBasis("Hstrain")],
    )
    prim = casmconfig.Prim(xtal_prim)
    supercell_set = casmconfig.SupercellSet(prim=prim)
    configs_as_enumerated = casmconfig.ConfigurationSet()
    canonical_configs = casmconfig.ConfigurationSet()

    dof_space = casmclex.DoFSpace(
        dof_key="Hstrain",
        xtal_prim=xtal_prim,
    )

    config_enum = casmenum.ConfigEnumMeshGrid(
        prim=prim,
        supercell_set=supercell_set,
    )

    supercell = supercell_set.add(
        np.array(
            [
                [1, 0, 0],
                [0, 1, 0],
                [0, 0, 1],
            ]
        )
    ).supercell
    background = casmconfig.Configuration(supercell)

    for i, configuration in enumerate(
        config_enum.by_range(
            background=background,
            dof_space=dof_space,
            start=-0.1,
            stop=0.1,
            num=5,
        )
    ):
        assert isinstance(configuration, casmconfig.Configuration)
        configs_as_enumerated.add(configuration)
        canonical_configs.add(casmconfig.make_canonical_configuration(configuration))

    # for i, record in enumerate(canonical_configs):
    #     print(xtal.pretty_json(record.configuration.to_dict()))

    assert len(configs_as_enumerated) == 5**6
    assert len(canonical_configs) == 900


def test_ConfigEnumMeshGrid_by_range_FCC_2():
    """Test ConfigEnumMeshGrid.by_range

    Using:
    - FCC,
    - prim unit cell, default background,
    - full basis,
    - num=5, skip_equivalents=True
    """
    xtal_prim = xtal_prims.FCC(
        r=0.5,
        occ_dof=["A", "B"],
        global_dof=[xtal.DoFSetBasis("Hstrain")],
    )
    prim = casmconfig.Prim(xtal_prim)
    supercell_set = casmconfig.SupercellSet(prim=prim)
    configs_as_enumerated = casmconfig.ConfigurationSet()
    canonical_configs = casmconfig.ConfigurationSet()

    dof_space = casmclex.DoFSpace(
        dof_key="Hstrain",
        xtal_prim=xtal_prim,
    )

    config_enum = casmenum.ConfigEnumMeshGrid(
        prim=prim,
        supercell_set=supercell_set,
    )

    supercell = supercell_set.add(
        np.array(
            [
                [1, 0, 0],
                [0, 1, 0],
                [0, 0, 1],
            ]
        )
    ).supercell
    background = casmconfig.Configuration(supercell)

    for i, configuration in enumerate(
        config_enum.by_range(
            background=background,
            dof_space=dof_space,
            start=-0.1,
            stop=0.1,
            num=5,
            skip_equivalents=True,
        )
    ):
        assert isinstance(configuration, casmconfig.Configuration)
        configs_as_enumerated.add(configuration)
        canonical_configs.add(casmconfig.make_canonical_configuration(configuration))

    # for i, record in enumerate(canonical_configs):
    #     print(xtal.pretty_json(record.configuration.to_dict()))

    assert len(configs_as_enumerated) == 900
    assert len(canonical_configs) == 900


def test_ConfigEnumMeshGrid_by_range_FCC_3():
    """Test ConfigEnumMeshGrid.by_range

    Using:
    - FCC,
    - prim unit cell, default background,
    - full basis,
    - step=0.05, skip_equivalents=True
    """
    xtal_prim = xtal_prims.FCC(
        r=0.5,
        occ_dof=["A", "B"],
        global_dof=[xtal.DoFSetBasis("Hstrain")],
    )
    prim = casmconfig.Prim(xtal_prim)
    supercell_set = casmconfig.SupercellSet(prim=prim)
    configs_as_enumerated = casmconfig.ConfigurationSet()
    canonical_configs = casmconfig.ConfigurationSet()

    dof_space = casmclex.DoFSpace(
        dof_key="Hstrain",
        xtal_prim=xtal_prim,
    )

    config_enum = casmenum.ConfigEnumMeshGrid(
        prim=prim,
        supercell_set=supercell_set,
    )

    supercell = supercell_set.add(
        np.array(
            [
                [1, 0, 0],
                [0, 1, 0],
                [0, 0, 1],
            ]
        )
    ).supercell
    background = casmconfig.Configuration(supercell)

    for i, configuration in enumerate(
        config_enum.by_range(
            background=background,
            dof_space=dof_space,
            start=-0.1,
            stop=0.101,
            step=0.05,
            skip_equivalents=True,
        )
    ):
        assert isinstance(configuration, casmconfig.Configuration)
        configs_as_enumerated.add(configuration)
        canonical_configs.add(casmconfig.make_canonical_configuration(configuration))

    # for i, record in enumerate(canonical_configs):
    #     print(xtal.pretty_json(record.configuration.to_dict()))

    assert len(configs_as_enumerated) == 900
    assert len(canonical_configs) == 900


def test_ConfigEnumMeshGrid_by_grid_coordinates_FCC_1():
    """Test ConfigEnumMeshGrid.by_grid_coordinates

    Using:
    - FCC,
    - prim unit cell, default background,
    - full basis,
    - step=0.05, skip_equivalents=True
    """
    xtal_prim = xtal_prims.FCC(
        r=0.5,
        occ_dof=["A", "B"],
        global_dof=[xtal.DoFSetBasis("Hstrain")],
    )
    prim = casmconfig.Prim(xtal_prim)
    supercell_set = casmconfig.SupercellSet(prim=prim)
    configs_as_enumerated = casmconfig.ConfigurationSet()
    canonical_configs = casmconfig.ConfigurationSet()

    dof_space = casmclex.DoFSpace(
        dof_key="Hstrain",
        xtal_prim=xtal_prim,
    )

    config_enum = casmenum.ConfigEnumMeshGrid(
        prim=prim,
        supercell_set=supercell_set,
    )

    supercell = supercell_set.add(
        np.array(
            [
                [1, 0, 0],
                [0, 1, 0],
                [0, 0, 1],
            ]
        )
    ).supercell
    background = casmconfig.Configuration(supercell)

    xi = []
    for d in range(dof_space.subspace_dim):
        xi.append(np.linspace(-0.1, 0.1, 5))

    for i, configuration in enumerate(
        config_enum.by_grid_coordinates(
            background=background,
            dof_space=dof_space,
            xi=xi,
            skip_equivalents=True,
        )
    ):
        assert isinstance(configuration, casmconfig.Configuration)
        configs_as_enumerated.add(configuration)
        canonical_configs.add(casmconfig.make_canonical_configuration(configuration))

    # for i, record in enumerate(canonical_configs):
    #     print(xtal.pretty_json(record.configuration.to_dict()))

    assert len(configs_as_enumerated) == 900
    assert len(canonical_configs) == 900


def test_ConfigEnumMeshGrid_by_irreducible_wedge_FCC_1():
    """Test ConfigEnumMeshGrid.by_irreducible_wedge

    Using:
    - FCC,
    - prim unit cell, default background,
    - full basis,
    - stop=0.1, num=5, skip_equivalents=default (False)
    """
    xtal_prim = xtal_prims.FCC(
        r=0.5,
        occ_dof=["A", "B"],
        global_dof=[xtal.DoFSetBasis("Hstrain")],
    )
    prim = casmconfig.Prim(xtal_prim)
    supercell_set = casmconfig.SupercellSet(prim=prim)
    configs_as_enumerated = casmconfig.ConfigurationSet()
    canonical_configs = casmconfig.ConfigurationSet()

    dof_space = casmclex.DoFSpace(
        dof_key="Hstrain",
        xtal_prim=xtal_prim,
    )
    dof_space_analysis_results = casmconfig.dof_space_analysis(
        dof_space=dof_space,
        prim=prim,
        calc_wedges=True,
    )
    symmetry_report = dof_space_analysis_results.symmetry_report

    config_enum = casmenum.ConfigEnumMeshGrid(
        prim=prim,
        supercell_set=supercell_set,
    )

    supercell = supercell_set.add(
        np.array(
            [
                [1, 0, 0],
                [0, 1, 0],
                [0, 0, 1],
            ]
        )
    ).supercell
    background = casmconfig.Configuration(supercell)

    for i, configuration in enumerate(
        config_enum.by_irreducible_wedge(
            background=background,
            dof_space=dof_space,
            irreducible_wedge=symmetry_report.irreducible_wedge,
            stop=0.1,
            num=5,
        )
    ):
        assert isinstance(configuration, casmconfig.Configuration)
        configs_as_enumerated.add(configuration)
        canonical_configs.add(casmconfig.make_canonical_configuration(configuration))

    # for i, record in enumerate(canonical_configs):
    #     print(xtal.pretty_json(record.configuration.to_dict()))

    assert len(configs_as_enumerated) == 4175
    assert len(canonical_configs) == 2157


def test_ConfigEnumMeshGrid_by_irreducible_wedge_FCC_2():
    """Test ConfigEnumMeshGrid.by_irreducible_wedge

    Using:
    - FCC,
    - prim unit cell, default background,
    - full basis,
    - stop=0.1, num=5, skip_equivalents=True
    """
    xtal_prim = xtal_prims.FCC(
        r=0.5,
        occ_dof=["A", "B"],
        global_dof=[xtal.DoFSetBasis("Hstrain")],
    )
    prim = casmconfig.Prim(xtal_prim)
    supercell_set = casmconfig.SupercellSet(prim=prim)
    configs_as_enumerated = casmconfig.ConfigurationSet()
    canonical_configs = casmconfig.ConfigurationSet()

    dof_space = casmclex.DoFSpace(
        dof_key="Hstrain",
        xtal_prim=xtal_prim,
    )
    dof_space_analysis_results = casmconfig.dof_space_analysis(
        dof_space=dof_space,
        prim=prim,
        calc_wedges=True,
    )
    symmetry_report = dof_space_analysis_results.symmetry_report

    config_enum = casmenum.ConfigEnumMeshGrid(
        prim=prim,
        supercell_set=supercell_set,
    )

    supercell = supercell_set.add(
        np.array(
            [
                [1, 0, 0],
                [0, 1, 0],
                [0, 0, 1],
            ]
        )
    ).supercell
    background = casmconfig.Configuration(supercell)

    for i, configuration in enumerate(
        config_enum.by_irreducible_wedge(
            background=background,
            dof_space=dof_space,
            irreducible_wedge=symmetry_report.irreducible_wedge,
            stop=0.1,
            num=5,
            skip_equivalents=True,
        )
    ):
        assert isinstance(configuration, casmconfig.Configuration)
        configs_as_enumerated.add(configuration)
        canonical_configs.add(casmconfig.make_canonical_configuration(configuration))

    # for i, record in enumerate(canonical_configs):
    #     print(xtal.pretty_json(record.configuration.to_dict()))

    assert len(configs_as_enumerated) == 2157
    assert len(canonical_configs) == 2157
