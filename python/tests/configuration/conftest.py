import builtins
from math import sqrt

import numpy as np
import pytest

import libcasm.configuration as config
import libcasm.xtal as xtal

from .functions import (
    make_discrete_magnetic_atom,
)


@pytest.fixture
def hide_available_pkg(monkeypatch):
    import_orig = builtins.__import__

    def mocked_import(name, *args, **kwargs):
        if name == "pkg":
            raise ImportError()
        return import_orig(name, *args, **kwargs)

    monkeypatch.setattr(builtins, "__import__", mocked_import)


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
def FCC_binary_prim():
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

    return xtal.Prim(
        lattice=lattice,
        coordinate_frac=coordinate_frac,
        occ_dof=occ_dof,
        occupants=occupants,
    )


@pytest.fixture
def FCC_binary_occ_fix_corner_prim():
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
            [0.5, 0.5, 0.0],
            [0.0, 0.5, 0.5],
            [0.5, 0.0, 0.5],
        ]
    ).transpose()

    # Occupation degrees of freedom (DoF)
    occupants = {}
    occ_dof = [
        ["A"],
        ["A", "B"],
        ["A", "B"],
        ["A", "B"],
    ]

    return xtal.Prim(
        lattice=lattice,
        coordinate_frac=coordinate_frac,
        occ_dof=occ_dof,
        occupants=occupants,
    )


@pytest.fixture
def FCC_binary_disp_fix_corner_prim():
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
            [0.5, 0.5, 0.0],
            [0.0, 0.5, 0.5],
            [0.5, 0.0, 0.5],
        ]
    ).transpose()

    # Occupation degrees of freedom (DoF)
    occupants = {}
    occ_dof = [
        ["A", "B"],
        ["A", "B"],
        ["A", "B"],
        ["A", "B"],
    ]

    # Local continuous degrees of freedom (DoF)
    disp_dof = xtal.DoFSetBasis("disp")  # Atomic displacement
    local_dof = [
        [],
        [disp_dof],
        [disp_dof],
        [disp_dof],
    ]

    return xtal.Prim(
        lattice=lattice,
        coordinate_frac=coordinate_frac,
        occ_dof=occ_dof,
        local_dof=local_dof,
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
def FCC_binary_Hstrain_noshear_prim():
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

    # Global continuous degrees of freedom (DoF)
    Hstrain_dof = xtal.DoFSetBasis(
        dofname="Hstrain",
        axis_names=["e_{1}", "e_{2}", "e_{3}"],
        basis=np.array(
            [
                [1.0 / sqrt(3), 1.0 / sqrt(3), 1.0 / sqrt(3), 0.0, 0.0, 0.0],
                [1.0 / sqrt(2), -1.0 / sqrt(2), 0.0, 0.0, 0.0, 0.0],
                [-1.0 / sqrt(6), -1.0 / sqrt(6), 2.0 / sqrt(6), 0.0, 0.0, 0.0],
            ]
        ).transpose(),
    )
    global_dof = [Hstrain_dof]

    return xtal.Prim(
        lattice=lattice,
        coordinate_frac=coordinate_frac,
        occ_dof=occ_dof,
        global_dof=global_dof,
        occupants=occupants,
    )


@pytest.fixture
def FCC_binary_Hstrain_disp_prim():
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
    Hstrain_dof = xtal.DoFSetBasis("Hstrain")  # Hencky strain metric
    global_dof = [Hstrain_dof]

    return xtal.Prim(
        lattice=lattice,
        coordinate_frac=coordinate_frac,
        occ_dof=occ_dof,
        local_dof=local_dof,
        global_dof=global_dof,
        occupants=occupants,
    )


@pytest.fixture
def FCC_binary_Hstrain_noshear_disp_nodz_prim():
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
    disp_dof = xtal.DoFSetBasis(
        dofname="disp",  # Atomic displacement
        axis_names=["d_{x}", "d_{y}"],
        basis=np.array(
            [
                [1.0, 0.0, 0.0],
                [0.0, 1.0, 0.0],
            ]
        ).transpose(),
    )
    local_dof = [
        [disp_dof],  # local DoF, basis site b=0
    ]

    # Global continuous degrees of freedom (DoF)
    Hstrain_dof = xtal.DoFSetBasis(
        dofname="Hstrain",
        axis_names=["e_{1}", "e_{2}", "e_{3}"],
        basis=np.array(
            [
                [1.0 / sqrt(3), 1.0 / sqrt(3), 1.0 / sqrt(3), 0.0, 0.0, 0.0],
                [1.0 / sqrt(2), -1.0 / sqrt(2), 0.0, 0.0, 0.0, 0.0],
                [-1.0 / sqrt(6), -1.0 / sqrt(6), 2.0 / sqrt(6), 0.0, 0.0, 0.0],
            ]
        ).transpose(),
    )
    global_dof = [Hstrain_dof]

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
def prim_ABC2():
    # Lattice vectors
    lattice_column_vector_matrix = np.array(
        [
            [0.0, 0.0, 3.30],  # a
            [4.87, 0.0, 0.0],  # a
            [0.0, 4.87, 0.0],  # a
        ]
    ).transpose()
    lattice = xtal.Lattice(lattice_column_vector_matrix)

    # Basis sites positions, as columns of a matrix,
    # in fractional coordinates with respect to the lattice vectors
    coordinate_frac = np.array(
        [
            [0.0, 0.0, 0.0],
            [0.500000000000, 0.500000000000, 0.500000000000],
            [0.304640000000, 0.000000000000, 0.304640000000],
            [0.804640000000, 0.500000000000, 0.195360010000],
            [0.195360010000, 0.500000000000, 0.804640000000],
            [0.695360000000, 0.000000000000, 0.695360000000],
        ]
    ).transpose()

    # Occupation degrees of freedom (DoF)
    occupants = {}
    occ_dof = [["A"], ["B"], ["C"], ["C"], ["C"], ["C"]]

    # Local continuous degrees of freedom (DoF)
    disp_dof = xtal.DoFSetBasis("disp")  # Atomic displacement
    local_dof = [
        [disp_dof],
        [disp_dof],
        [disp_dof],
        [disp_dof],
        [disp_dof],
        [disp_dof],
    ]

    return xtal.Prim(
        lattice=lattice,
        coordinate_frac=coordinate_frac,
        occ_dof=occ_dof,
        local_dof=local_dof,
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


@pytest.fixture
def ZrO_prim():
    # Lattice vectors
    lattice_column_vector_matrix = np.array(
        [
            [3.233986856383, 0.000000000000, 0.000000000000],  # a
            [-1.616993428191, 2.800714773133, 0.000000000000],  # a
            [0.000000000000, 0.000000000000, 5.168678340000],  # c
        ]
    ).transpose()
    lattice = xtal.Lattice(lattice_column_vector_matrix)

    # Basis sites positions, as columns of a matrix,
    # in fractional coordinates with respect to the lattice vectors
    coordinate_frac = np.array(
        [
            [0.0, 0.0, 0.0],
            [2.0 / 3.0, 1.0 / 3.0, 1.0 / 2.0],
            [1.0 / 3.0, 2.0 / 3.0, 1.0 / 4.0],
            [1.0 / 3.0, 2.0 / 3.0, 3.0 / 4.0],
        ]
    ).transpose()

    # Occupation degrees of freedom (DoF)
    occupants = {}
    occ_dof = [["Zr"], ["Zr"], ["O", "Va"], ["O", "Va"]]

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
def FCC_binary_discrete_Cmagspin_prim():
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
    occupants = {
        "A.up": make_discrete_magnetic_atom(name="A", value=1, flavor="C"),
        "A.down": make_discrete_magnetic_atom(name="A", value=-1, flavor="C"),
        "B.up": make_discrete_magnetic_atom(name="B", value=1, flavor="C"),
        "B.down": make_discrete_magnetic_atom(name="B", value=-1, flavor="C"),
    }
    occ_dof = [["A.up", "A.down", "B.up", "B.down"]]

    return xtal.Prim(
        lattice=lattice,
        coordinate_frac=coordinate_frac,
        occ_dof=occ_dof,
        occupants=occupants,
    )
