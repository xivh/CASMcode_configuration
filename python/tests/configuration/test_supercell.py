import math
import numpy as np
import pytest
import libcasm.xtal as xtal
import libcasm.configuration as config


def test_simple_cubic_binary_supercell(simple_cubic_binary_prim):
    prim = config.Prim(simple_cubic_binary_prim)
    T = np.array(
        [
            [2, 1, 0],
            [0, 1, 0],
            [0, 0, 1],
        ]
    )
    supercell = config.Supercell(prim, T)
    S = supercell.superlattice().column_vector_matrix()
    assert math.isclose(S[0, 0], 2.0)
    assert math.isclose(S[0, 1], 1.0)
    assert math.isclose(S[1, 1], 1.0)
    assert math.isclose(S[2, 2], 1.0)
    assert supercell.n_sites() == 2
    assert supercell.n_unitcells() == 2
    assert supercell.n_occupants() == [2, 2]
    assert supercell.occ_dof() == [["A", "B"], ["A", "B"]]
    assert np.allclose(
        supercell.coordinate_frac(),
        np.array([[0.0, 0.0, 0.0], [1.0, 0.0, 0.0]]).transpose(),
    )
    assert np.allclose(
        supercell.coordinate_cart(),
        np.array([[0.0, 0.0, 0.0], [1.0, 0.0, 0.0]]).transpose(),
    )
    assert supercell.sublattice_indices() == [0, 0]
    assert (
        supercell.unitcell_indices() == np.array([[0, 0, 0], [1, 0, 0]]).transpose()
    ).all()


def test_simple_cubic_binary_supercell_compare(simple_cubic_binary_prim):
    prim = config.Prim(simple_cubic_binary_prim)
    T1 = np.array(
        [
            [2, 0, 0],
            [0, 1, 0],
            [0, 0, 1],
        ]
    )
    supercell1 = config.Supercell(prim, T1)

    T1b = np.array(
        [
            [2, 0, 0],
            [0, 1, 0],
            [0, 0, 1],
        ]
    )
    supercell1b = config.Supercell(prim, T1b)

    T2 = np.array(
        [
            [1, 0, 0],
            [0, 2, 0],
            [0, 0, 1],
        ]
    )
    supercell2 = config.Supercell(prim, T2)

    T3 = np.array(
        [
            [1, 0, 0],
            [0, 1, 0],
            [0, 0, 2],
        ]
    )
    supercell3 = config.Supercell(prim, T3)

    assert supercell2 > supercell1
    assert supercell3 > supercell1
    assert supercell3 > supercell2

    assert supercell1 < supercell2
    assert supercell1 < supercell3
    assert supercell2 < supercell3

    assert supercell1 != supercell2
    assert supercell1 == supercell1b

    expected_S_canonical = np.array(
        [
            [1.0, 0.0, 0.0],
            [0.0, 1.0, 0.0],
            [0.0, 0.0, 2.0],
        ]
    )
    canonical_supercell = config.make_canonical_supercell(supercell1)
    assert canonical_supercell == supercell3
    assert config.is_canonical_supercell(supercell1) == False
    assert config.is_canonical_supercell(supercell2) == False
    assert config.is_canonical_supercell(supercell3) == True
    assert np.allclose(
        canonical_supercell.superlattice().column_vector_matrix(), expected_S_canonical
    )

    equivalent_supercells = config.make_equivalent_supercells(supercell1)
    for scel in equivalent_supercells:
        print(scel.superlattice().column_vector_matrix())
    assert len(equivalent_supercells) == 3
