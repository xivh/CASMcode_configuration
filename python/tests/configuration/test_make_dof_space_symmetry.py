import numpy as np

import libcasm.clexulator as casmclex
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

    # Occupation DoF space symmetry
    dof_space = casmclex.DoFSpace(
        dof_key="occ",
        xtal_prim=prim.xtal_prim,
        transformation_matrix_to_super=supercell.transformation_matrix_to_super,
    )
    matrix_rep, symgroup = casmconfig.make_dof_space_symmetry(
        group=group,
        dof_space=dof_space,
    )
    assert len(matrix_rep) == 48 * 8
    assert len(symgroup.elements) == 48 * 8
    count = 0
    for op in symgroup.elements:
        if not np.allclose(op.translation(), zeros):
            count += 1
    assert count == 48 * 7

    # Displacement DoF space symmetry
    dof_space = casmclex.DoFSpace(
        dof_key="disp",
        xtal_prim=prim.xtal_prim,
        transformation_matrix_to_super=supercell.transformation_matrix_to_super,
    )
    matrix_rep, symgroup = casmconfig.make_dof_space_symmetry(
        group=group,
        dof_space=dof_space,
    )
    assert len(matrix_rep) == 48 * 8
    assert len(symgroup.elements) == 48 * 8
    count = 0
    for op in symgroup.elements:
        if not np.allclose(op.translation(), zeros):
            count += 1
    assert count == 48 * 7

    # GLstrain DoF space symmetry - result is a point group
    dof_space = casmclex.DoFSpace(
        dof_key="GLstrain",
        xtal_prim=prim.xtal_prim,
    )
    matrix_rep, symgroup = casmconfig.make_dof_space_symmetry(
        group=group,
        dof_space=dof_space,
    )
    assert len(matrix_rep) == 48
    assert len(symgroup.elements) == 48
    count = 0
    for op in symgroup.elements:
        if not np.allclose(op.translation(), zeros):
            count += 1
    assert count == 0


def test_make_dof_space_symmetry_2(ZrO_prim_GLstrain_disp):
    xtal_prim = ZrO_prim_GLstrain_disp
    prim = casmconfig.Prim(xtal_prim=xtal_prim)
    T = np.eye(3, dtype=int) * 2
    supercell = casmconfig.Supercell(prim, T)
    group = supercell.symgroup_rep()
    zeros = np.zeros(3)

    # Occupation DoF space symmetry
    dof_space = casmclex.DoFSpace(
        dof_key="occ",
        xtal_prim=prim.xtal_prim,
        transformation_matrix_to_super=supercell.transformation_matrix_to_super,
    )
    matrix_rep, symgroup = casmconfig.make_dof_space_symmetry(
        group=group,
        dof_space=dof_space,
    )
    assert len(matrix_rep) == 24 * 8
    assert len(symgroup.elements) == 24 * 8
    count = 0
    for op in symgroup.elements:
        if not np.allclose(op.translation(), zeros):
            count += 1
    assert count == 24 * 7 + 12

    # Displacement DoF space symmetry
    dof_space = casmclex.DoFSpace(
        dof_key="disp",
        xtal_prim=prim.xtal_prim,
        transformation_matrix_to_super=supercell.transformation_matrix_to_super,
    )
    matrix_rep, symgroup = casmconfig.make_dof_space_symmetry(
        group=group,
        dof_space=dof_space,
    )
    assert len(matrix_rep) == 24 * 8
    assert len(symgroup.elements) == 24 * 8
    count = 0
    for op in symgroup.elements:
        if not np.allclose(op.translation(), zeros):
            count += 1
    assert count == 24 * 7 + 12

    # GLstrain DoF space symmetry - result is a point group
    dof_space = casmclex.DoFSpace(
        dof_key="GLstrain",
        xtal_prim=prim.xtal_prim,
    )
    matrix_rep, symgroup = casmconfig.make_dof_space_symmetry(
        group=group,
        dof_space=dof_space,
    )
    assert len(matrix_rep) == 24
    assert len(symgroup.elements) == 24
    count = 0
    for op in symgroup.elements:
        if not np.allclose(op.translation(), zeros):
            count += 1
    assert count == 0
