from math import sqrt

import numpy as np

import libcasm.configuration as casmconfig


def make_configuration(
    supercell: casmconfig.Supercell,
    occupation: list[int],
):
    config = casmconfig.Configuration(supercell)
    config.set_occupation(occupation)
    return config


def build_configurations_1(fcc_binary_prim: casmconfig.Prim):
    prim = fcc_binary_prim
    T_prim = np.eye(3, dtype=int)
    T_conventional = np.array(
        [  # conventional FCC cubic cell
            [-1, 1, 1],
            [1, -1, 1],
            [1, 1, -1],
        ],
        dtype=int,
    )

    # construct configurations, the resulting space will span the space spanned by
    # these configurations and the orbits of equivalent configurations
    configurations = {}

    supercell = casmconfig.make_canonical_supercell(casmconfig.Supercell(prim, T_prim))
    configurations["A1"] = make_configuration(supercell, [0])
    configurations["B1"] = make_configuration(supercell, [1])

    supercell = casmconfig.make_canonical_supercell(
        casmconfig.Supercell(prim, T_conventional)
    )
    configurations["A3B1"] = make_configuration(supercell, [1, 0, 0, 0])
    configurations["A1B3"] = make_configuration(supercell, [1, 1, 1, 0])

    return configurations


def test_config_space_analysis_1(FCC_binary_prim):
    prim = casmconfig.Prim(FCC_binary_prim)

    # Perform DoF space analysis
    results = casmconfig.config_space_analysis(
        configurations=build_configurations_1(prim),
        # dofs=None,
        # exclude_homogeneous_modes=None,
        # include_default_occ_modes=False,
        # sublattice_index_to_default_occ=None,
        # site_index_to_default_occ=None,
        # tol=casmglobal.TOL,
    )
    assert len(results) == 1

    symmetry_adapted_dof_space = results["occ"].symmetry_adapted_dof_space
    print("basis:\n", symmetry_adapted_dof_space.basis)

    assert symmetry_adapted_dof_space.basis.shape[0] == 8
    assert symmetry_adapted_dof_space.basis.shape[1] == 4

    ### equivalent configurations
    equivalent_configurations = results["occ"].equivalent_configurations
    assert len(equivalent_configurations["A1"]) == 1
    assert len(equivalent_configurations["B1"]) == 1
    assert len(equivalent_configurations["A3B1"]) == 4
    assert len(equivalent_configurations["A1B3"]) == 4

    ### equivalent dof values
    equivalent_dof_values = results["occ"].equivalent_dof_values
    assert len(equivalent_dof_values["A1"]) == 1
    assert len(equivalent_dof_values["B1"]) == 1
    assert len(equivalent_dof_values["A3B1"]) == 4
    assert len(equivalent_dof_values["A1B3"]) == 4

    expected = np.array(
        [
            [
                0.0,
                0.0,
                0.0,
                sqrt(6.0) / 3.0,
                0.0,
                -sqrt(6.0) / 6.0,
                0.0,
                -sqrt(6.0) / 6.0,
            ],
            [0.0, 0.0, 0.0, 0.0, 0.0, -sqrt(2.0) / 2.0, 0.0, sqrt(2.0) / 2.0],
            [
                0.0,
                -sqrt(3.0) / 2.0,
                0.0,
                sqrt(3.0) / 6.0,
                0.0,
                sqrt(3.0) / 6.0,
                0.0,
                sqrt(3.0) / 6.0,
            ],
            [0.0, 0.5, 0.0, 0.5, 0.0, 0.5, 0.0, 0.5],
        ]
    ).transpose()

    assert np.allclose(symmetry_adapted_dof_space.basis, expected)


def test_config_space_analysis_2(FCC_binary_prim):
    prim = casmconfig.Prim(FCC_binary_prim)

    sublattice_index_to_default_occ = {0: 1}

    # Perform DoF space analysis
    results = casmconfig.config_space_analysis(
        configurations=build_configurations_1(prim),
        # dofs=None,
        # exclude_homogeneous_modes=None,
        # include_default_occ_modes=False,
        sublattice_index_to_default_occ=sublattice_index_to_default_occ,
        # site_index_to_default_occ=None,
        # tol=casmglobal.TOL,
    )
    assert len(results) == 1

    symmetry_adapted_dof_space = results["occ"].symmetry_adapted_dof_space
    print("basis:\n", symmetry_adapted_dof_space.basis)

    assert symmetry_adapted_dof_space.basis.shape[0] == 8
    assert symmetry_adapted_dof_space.basis.shape[1] == 4

    ### equivalent configurations
    equivalent_configurations = results["occ"].equivalent_configurations
    assert len(equivalent_configurations["A1"]) == 1
    assert len(equivalent_configurations["B1"]) == 1
    assert len(equivalent_configurations["A3B1"]) == 4
    assert len(equivalent_configurations["A1B3"]) == 4

    ### equivalent dof values
    equivalent_dof_values = results["occ"].equivalent_dof_values
    assert len(equivalent_dof_values["A1"]) == 1
    assert len(equivalent_dof_values["B1"]) == 1
    assert len(equivalent_dof_values["A3B1"]) == 4
    assert len(equivalent_dof_values["A1B3"]) == 4

    expected = np.array(
        [
            [
                0.0,
                0.0,
                sqrt(6.0) / 3.0,
                0.0,
                -sqrt(6.0) / 6.0,
                0.0,
                -sqrt(6.0) / 6.0,
                0.0,
            ],
            [0.0, 0.0, 0.0, 0.0, -sqrt(2.0) / 2.0, 0.0, sqrt(2.0) / 2.0, 0.0],
            [
                -sqrt(3.0) / 2.0,
                0.0,
                sqrt(3.0) / 6.0,
                0.0,
                sqrt(3.0) / 6.0,
                0.0,
                sqrt(3.0) / 6.0,
                0.0,
            ],
            [0.5, 0.0, 0.5, 0.0, 0.5, 0.0, 0.5, 0.0],
        ]
    ).transpose()

    assert np.allclose(symmetry_adapted_dof_space.basis, expected)


def test_config_space_analysis_3(FCC_binary_prim):
    prim = casmconfig.Prim(FCC_binary_prim)

    site_index_to_default_occ = {0: 1, 1: 1, 2: 0, 3: 0}

    # Perform DoF space analysis
    results = casmconfig.config_space_analysis(
        configurations=build_configurations_1(prim),
        # dofs=None,
        # exclude_homogeneous_modes=None,
        # include_default_occ_modes=False,
        # sublattice_index_to_default_occ=None,
        site_index_to_default_occ=site_index_to_default_occ,
        # tol=casmglobal.TOL,
    )
    assert len(results) == 1

    symmetry_adapted_dof_space = results["occ"].symmetry_adapted_dof_space
    print("basis:\n", symmetry_adapted_dof_space.basis)

    assert symmetry_adapted_dof_space.basis.shape[0] == 8
    assert symmetry_adapted_dof_space.basis.shape[1] == 4

    ### equivalent configurations
    equivalent_configurations = results["occ"].equivalent_configurations
    assert len(equivalent_configurations["A1"]) == 1
    assert len(equivalent_configurations["B1"]) == 1
    assert len(equivalent_configurations["A3B1"]) == 4
    assert len(equivalent_configurations["A1B3"]) == 4

    ### equivalent dof values
    equivalent_dof_values = results["occ"].equivalent_dof_values
    assert len(equivalent_dof_values["A1"]) == 1
    assert len(equivalent_dof_values["B1"]) == 1
    assert len(equivalent_dof_values["A3B1"]) == 4
    assert len(equivalent_dof_values["A1B3"]) == 4

    expected = np.array(
        [
            [0.0, 0.0, 0.0, 0.0, 0.0, -sqrt(2.0) / 2.0, 0.0, sqrt(2.0) / 2.0],
            [-sqrt(2.0) / 2.0, 0.0, sqrt(2.0) / 2.0, 0.0, 0.0, 0.0, 0.0, 0.0],
            [-0.5, 0.0, -0.5, 0.0, 0.0, 0.5, 0.0, 0.5],
            [0.5, 0.0, 0.5, 0.0, 0.0, 0.5, 0.0, 0.5],
        ]
    ).transpose()

    assert np.allclose(symmetry_adapted_dof_space.basis, expected)


def test_config_space_analysis_4(FCC_binary_prim):
    prim = casmconfig.Prim(FCC_binary_prim)

    site_index_to_default_occ = {0: 0, 1: 0, 2: 1, 3: 1}

    # Perform DoF space analysis
    results = casmconfig.config_space_analysis(
        configurations=build_configurations_1(prim),
        # dofs=None,
        # exclude_homogeneous_modes=None,
        # include_default_occ_modes=False,
        # sublattice_index_to_default_occ=None,
        site_index_to_default_occ=site_index_to_default_occ,
        # tol=casmglobal.TOL,
    )
    assert len(results) == 1

    symmetry_adapted_dof_space = results["occ"].symmetry_adapted_dof_space
    print("basis:\n", symmetry_adapted_dof_space.basis)

    assert symmetry_adapted_dof_space.basis.shape[0] == 8
    assert symmetry_adapted_dof_space.basis.shape[1] == 4

    ### equivalent configurations
    equivalent_configurations = results["occ"].equivalent_configurations
    assert len(equivalent_configurations["A1"]) == 1
    assert len(equivalent_configurations["B1"]) == 1
    assert len(equivalent_configurations["A3B1"]) == 4
    assert len(equivalent_configurations["A1B3"]) == 4

    ### equivalent dof values
    equivalent_dof_values = results["occ"].equivalent_dof_values
    assert len(equivalent_dof_values["A1"]) == 1
    assert len(equivalent_dof_values["B1"]) == 1
    assert len(equivalent_dof_values["A3B1"]) == 4
    assert len(equivalent_dof_values["A1B3"]) == 4

    expected = np.array(
        [
            [0.0, 0.0, 0.0, 0.0, -sqrt(2.0) / 2.0, 0.0, sqrt(2.0) / 2.0, 0.0],
            [0.0, -sqrt(2.0) / 2.0, 0.0, sqrt(2.0) / 2.0, 0.0, 0.0, 0.0, 0.0],
            [0.0, -0.5, 0.0, -0.5, 0.5, 0.0, 0.5, 0.0],
            [0.0, 0.5, 0.0, 0.5, 0.5, 0.0, 0.5, 0.0],
        ]
    ).transpose()

    assert np.allclose(symmetry_adapted_dof_space.basis, expected)
