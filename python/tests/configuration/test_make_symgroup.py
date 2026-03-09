import numpy as np

import libcasm.configuration as casmconfig
import libcasm.xtal as xtal
import libcasm.xtal.prims as xtal_prims


def test_make_dof_space_symmetry_1():
    xtal_prim = xtal_prims.BCC(
        a=3.55,
        occ_dof=["A", "B"],
        local_dof=[xtal.DoFSetBasis("disp")],
        global_dof=[xtal.DoFSetBasis("GLstrain")],
    )
    prim = casmconfig.Prim(xtal_prim=xtal_prim)
    T = np.eye(3, dtype=int) * 2
    supercell = casmconfig.Supercell(prim, T)
    group = supercell.symgroup_rep()
    zeros = np.zeros(3)

    # Displacement / Occ DoF space symmetry - default `point_group`
    symgroup = casmconfig.make_symgroup(
        group=group,
    )
    assert len(symgroup.elements) == 48 * 8
    count = 0
    for op in symgroup.elements:
        if not np.allclose(op.translation(), zeros):
            count += 1
    assert count == 48 * 7

    # Displacement / Occ DoF space symmetry - explicit `point_group=False`
    symgroup = casmconfig.make_symgroup(
        group=group,
        point_group=False,
    )
    assert len(symgroup.elements) == 48 * 8
    count = 0
    for op in symgroup.elements:
        if not np.allclose(op.translation(), zeros):
            count += 1
    assert count == 48 * 7

    # GLstrain DoF space symmetry - explicit point_group=True
    symgroup = casmconfig.make_symgroup(
        group=group,
        point_group=True,
    )
    assert len(symgroup.elements) == 48
    count = 0
    for op in symgroup.elements:
        if not np.allclose(op.translation(), zeros):
            count += 1
    assert count == 0


def test_make_dof_space_symmetry_2(ZrO_prim):
    xtal_prim = ZrO_prim
    prim = casmconfig.Prim(xtal_prim=xtal_prim)
    T = np.eye(3, dtype=int) * 2
    supercell = casmconfig.Supercell(prim, T)
    group = supercell.symgroup_rep()
    zeros = np.zeros(3)

    # Displacement / Occ DoF space symmetry - default `point_group`
    symgroup = casmconfig.make_symgroup(
        group=group,
    )
    assert len(symgroup.elements) == 24 * 8
    count = 0
    for op in symgroup.elements:
        if not np.allclose(op.translation(), zeros):
            count += 1
    assert count == 24 * 7 + 12

    # Displacement / Occ DoF space symmetry - explicit `point_group=False`
    symgroup = casmconfig.make_symgroup(
        group=group,
        point_group=False,
    )
    assert len(symgroup.elements) == 24 * 8
    count = 0
    for op in symgroup.elements:
        if not np.allclose(op.translation(), zeros):
            count += 1
    assert count == 24 * 7 + 12

    # GLstrain DoF space symmetry - explicit point_group=True
    symgroup = casmconfig.make_symgroup(
        group=group,
        point_group=True,
    )
    assert len(symgroup.elements) == 24
    count = 0
    for op in symgroup.elements:
        if not np.allclose(op.translation(), zeros):
            count += 1
    assert count == 0
