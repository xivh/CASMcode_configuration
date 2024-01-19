import numpy as np

import libcasm.configuration as casmconfig
import libcasm.enumerate as casmenum
import libcasm.xtal.prims as xtal_prims


def test_ScelEnum_FCC_1():
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
    for i, supercell in enumerate(scel_enum.by_volume(max=4)):
        assert isinstance(supercell, casmconfig.Supercell)

    assert len(supercell_set) == 13


def test_ScelEnum_FCC_2():
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
    for i, supercell in enumerate(scel_enum.by_volume(max=4, diagonal_only=True)):
        assert isinstance(supercell, casmconfig.Supercell)

    assert len(supercell_set) == 5


def test_ScelEnum_FCC_3():
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
    for i, supercell in enumerate(
        scel_enum.by_volume(max=10, diagonal_only=True, fixed_shape=True)
    ):
        assert isinstance(supercell, casmconfig.Supercell)

    assert len(supercell_set) == 2


def test_ScelEnum_FCC_4():
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
    unit_cell = np.array([[2, 0, 0], [0, 1, 0], [0, 0, 1]])
    for i, supercell in enumerate(scel_enum.by_volume(max=4, unit_cell=unit_cell)):
        assert isinstance(supercell, casmconfig.Supercell)
    assert len(supercell_set) == 19


def test_ScelEnum_HCP_1():
    xtal_prim = xtal_prims.HCP(
        r=0.5,
        occ_dof=["A", "B"],
    )
    prim = casmconfig.Prim(xtal_prim)
    supercell_set = casmconfig.SupercellSet(prim=prim)
    scel_enum = casmenum.ScelEnum(
        prim=prim,
        supercell_set=supercell_set,
    )
    for i, supercell in enumerate(scel_enum.by_volume(max=4)):
        assert isinstance(supercell, casmconfig.Supercell)

    assert len(supercell_set) == 20


def test_ScelEnum_lowsym_1(lowsym_occ_prim):
    xtal_prim = lowsym_occ_prim
    prim = casmconfig.Prim(xtal_prim)
    supercell_set = casmconfig.SupercellSet(prim=prim)
    scel_enum = casmenum.ScelEnum(
        prim=prim,
        supercell_set=supercell_set,
    )
    for i, supercell in enumerate(scel_enum.by_volume(max=4)):
        assert isinstance(supercell, casmconfig.Supercell)

    assert len(supercell_set) == 56
