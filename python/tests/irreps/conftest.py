import numpy as np
import pytest

import libcasm.xtal as xtal


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
