import numpy as np
import pytest

import libcasm.configuration as config
import libcasm.xtal as xtal


@pytest.fixture
def tetragonal_lattice():
    # Lattice vectors
    lattice_column_vector_matrix = np.array(
        [
            [1.0, 0.0, 0.0],  # a
            [0.0, 1.0, 0.0],  # a
            [0.0, 0.0, 2.0],  # c
        ]
    ).transpose()
    return xtal.Lattice(lattice_column_vector_matrix)


@pytest.fixture
def simple_cubic_binary_prim():
    # Lattice vectors
    lattice_column_vector_matrix = np.array(
        [
            [1.0, 0.0, 0.0],  # a
            [0.0, 1.0, 0.0],  # a
            [0.0, 0.0, 1.0],  # a
        ]
    ).transpose()
    lattice = xtal.Lattice(lattice_column_vector_matrix)

    # Basis sites positions, as columns of a matrix,
    # in fractional coordinates with respect to the lattice vectors
    coordinate_frac = np.array(
        [
            [0.0, 0.0, 0.0],
        ]
    ).transpose()

    # Occupation degrees of freedom (DoF)
    occupants = {}
    occ_dof = [["A", "B"]]

    # Local continuous degrees of freedom (DoF)
    local_dof = []

    # Global continuous degrees of freedom (DoF)
    global_dof = []

    return xtal.Prim(
        lattice=lattice,
        coordinate_frac=coordinate_frac,
        occ_dof=occ_dof,
        local_dof=local_dof,
        global_dof=global_dof,
        occupants=occupants,
    )


@pytest.fixture
def FCC_binary_GLstrain_disp_prim():
    # Lattice vectors
    lattice_column_vector_matrix = np.array(
        [
            [0.0, 1.0 / 2.0, 1.0 / 2.0],  # a
            [1.0 / 2.0, 0.0, 1.0 / 2.0],  # a
            [1.0 / 2.0, 1.0 / 2.0, 0.0],  # a
        ]
    ).transpose()
    lattice = xtal.Lattice(lattice_column_vector_matrix)

    # Basis sites positions, as columns of a matrix,
    # in fractional coordinates with respect to the lattice vectors
    coordinate_frac = np.array(
        [
            [0.0, 0.0, 0.0],
        ]
    ).transpose()

    # Occupation degrees of freedom (DoF)
    occupants = {}
    occ_dof = [["A", "B"]]

    # Local continuous degrees of freedom (DoF)
    disp_dof = xtal.DoFSetBasis("disp")  # Atomic displacement
    local_dof = [
        [disp_dof],  # local DoF, basis site b=0
    ]

    # Global continuous degrees of freedom (DoF)
    GLstrain_dof = xtal.DoFSetBasis("GLstrain")  # Green-Lagrange strain metric
    global_dof = [GLstrain_dof]

    return xtal.Prim(
        lattice=lattice,
        coordinate_frac=coordinate_frac,
        occ_dof=occ_dof,
        local_dof=local_dof,
        global_dof=global_dof,
        occupants=occupants,
    )


@pytest.fixture
def BCC_binary_GLstrain_disp_prim():
    # Lattice vectors
    lattice_column_vector_matrix = np.array(
        [
            [-1.0 / 2.0, 1.0 / 2.0, 1.0 / 2.0],  # a
            [1.0 / 2.0, -1.0 / 2.0, 1.0 / 2.0],  # a
            [1.0 / 2.0, 1.0 / 2.0, -1.0 / 2.0],  # a
        ]
    ).transpose()
    lattice = xtal.Lattice(lattice_column_vector_matrix)

    # Basis sites positions, as columns of a matrix,
    # in fractional coordinates with respect to the lattice vectors
    coordinate_frac = np.array(
        [
            [0.0, 0.0, 0.0],
        ]
    ).transpose()

    # Occupation degrees of freedom (DoF)
    occupants = {}
    occ_dof = [["A", "B"]]

    # Local continuous degrees of freedom (DoF)
    disp_dof = xtal.DoFSetBasis("disp")  # Atomic displacement
    local_dof = [
        [disp_dof],  # local DoF, basis site b=0
    ]

    # Global continuous degrees of freedom (DoF)
    GLstrain_dof = xtal.DoFSetBasis("GLstrain")  # Green-Lagrange strain metric
    global_dof = [GLstrain_dof]

    return xtal.Prim(
        lattice=lattice,
        coordinate_frac=coordinate_frac,
        occ_dof=occ_dof,
        local_dof=local_dof,
        global_dof=global_dof,
        occupants=occupants,
    )


@pytest.fixture
def simple_cubic_binary_SupercellSet_1_canonical():
    prim = config.Prim(simple_cubic_binary_prim)
    supercells = config.SupercellSet(prim)
    assert supercells.empty()

    T = np.array(
        [
            [2, 1, 0],
            [0, 1, 0],
            [0, 0, 1],
        ]
    )
    supercell = config.Supercell(prim, T)
    supercell = config.make_canonical_supercell(supercell)
    supercells.add(supercell)
    return supercells


@pytest.fixture
def simple_cubic_binary_SupercellSet_1_non_canonical():
    prim = config.Prim(simple_cubic_binary_prim)
    supercells = config.SupercellSet(prim)
    assert supercells.empty()

    T = np.array(
        [
            [2, 1, 0],
            [0, 1, 0],
            [0, 0, 1],
        ]
    )
    supercell = config.Supercell(prim, T)
    supercells.add(supercell)
    return supercells
