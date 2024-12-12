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


def is_same_space(found, expected):
    rank_found = np.linalg.matrix_rank(found)
    rank_expected = np.linalg.matrix_rank(expected)
    if rank_found != rank_expected:
        return False

    augmented = np.hstack((found, expected))
    rank_augmented = np.linalg.matrix_rank(augmented)
    return rank_augmented == rank_expected


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

    ### eigenvalues
    expected = np.array([2.0, 2.0, 2.0, 14.0])
    assert np.allclose(results["occ"].eigenvalues, expected)

    ### basis
    expected_1 = np.array(
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
        ]
    ).transpose()

    expected_2 = np.array(
        [
            [0.0, 0.5, 0.0, 0.5, 0.0, 0.5, 0.0, 0.5],
        ]
    ).transpose()

    found_1 = symmetry_adapted_dof_space.basis[:, :3]
    found_2 = symmetry_adapted_dof_space.basis[:, 3:]

    assert is_same_space(found_1, expected_1)
    assert is_same_space(found_2, expected_2)


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

    ### eigenvalues
    expected = np.array([2.0, 2.0, 2.0, 14.0])
    assert np.allclose(results["occ"].eigenvalues, expected)

    ### basis
    expected_1 = np.array(
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
        ]
    ).transpose()

    expected_2 = np.array(
        [
            [0.5, 0.0, 0.5, 0.0, 0.5, 0.0, 0.5, 0.0],
        ]
    ).transpose()

    found_1 = symmetry_adapted_dof_space.basis[:, :3]
    found_2 = symmetry_adapted_dof_space.basis[:, 3:]

    assert is_same_space(found_1, expected_1)
    assert is_same_space(found_2, expected_2)


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

    ### eigenvalues
    expected = np.array([2.0, 2.0, 4.0, 12.0])
    assert np.allclose(results["occ"].eigenvalues, expected)

    ### basis
    expected_1 = np.array(
        [
            [0.0, 0.0, 0.0, 0.0, 0.0, -sqrt(2.0) / 2.0, 0.0, sqrt(2.0) / 2.0],
            [-sqrt(2.0) / 2.0, 0.0, sqrt(2.0) / 2.0, 0.0, 0.0, 0.0, 0.0, 0.0],
        ]
    ).transpose()

    expected_2 = np.array(
        [
            [-0.5, 0.0, -0.5, 0.0, 0.0, 0.5, 0.0, 0.5],
        ]
    ).transpose()

    expected_3 = np.array(
        [
            [0.5, 0.0, 0.5, 0.0, 0.0, 0.5, 0.0, 0.5],
        ]
    ).transpose()

    found_1 = symmetry_adapted_dof_space.basis[:, :2]
    found_2 = symmetry_adapted_dof_space.basis[:, 2:3]
    found_3 = symmetry_adapted_dof_space.basis[:, 3:]

    assert is_same_space(found_1, expected_1)
    assert is_same_space(found_2, expected_2)
    assert is_same_space(found_3, expected_3)


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

    ### eigenvalues
    expected = np.array([2.0, 2.0, 4.0, 12.0])
    assert np.allclose(results["occ"].eigenvalues, expected)

    ### basis
    expected_1 = np.array(
        [
            [0.0, 0.0, 0.0, 0.0, -sqrt(2.0) / 2.0, 0.0, sqrt(2.0) / 2.0, 0.0],
            [0.0, -sqrt(2.0) / 2.0, 0.0, sqrt(2.0) / 2.0, 0.0, 0.0, 0.0, 0.0],
        ]
    ).transpose()

    expected_2 = np.array(
        [
            [0.0, -0.5, 0.0, -0.5, 0.5, 0.0, 0.5, 0.0],
        ]
    ).transpose()

    expected_3 = np.array(
        [
            [0.0, 0.5, 0.0, 0.5, 0.5, 0.0, 0.5, 0.0],
        ]
    ).transpose()

    found_1 = symmetry_adapted_dof_space.basis[:, 0:2]
    found_2 = symmetry_adapted_dof_space.basis[:, 2:3]
    found_3 = symmetry_adapted_dof_space.basis[:, 3:]

    assert is_same_space(found_1, expected_1)
    assert is_same_space(found_2, expected_2)
    assert is_same_space(found_3, expected_3)


def test_config_space_analysis_1_to_dict(FCC_binary_prim):
    prim = casmconfig.Prim(FCC_binary_prim)

    # Perform DoF space analysis
    results = casmconfig.config_space_analysis(
        configurations=build_configurations_1(prim),
    )
    assert len(results) == 1

    data = {key: value.to_dict() for key, value in results.items()}
    assert len(data) == 1
    assert "occ" in data
    assert "eigenvalues" in data["occ"]
    assert "equivalent_configurations" in data["occ"]
    assert "equivalent_dof_values" in data["occ"]
    assert "projector" in data["occ"]
    assert "standard_dof_space" in data["occ"]
    assert "symmetry_adapted_dof_space" in data["occ"]
