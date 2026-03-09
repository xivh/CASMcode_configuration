import json
import pathlib

import numpy as np

import libcasm.clexulator as casmclex
import libcasm.configuration as casmconfig
import libcasm.group as casmgroup
import libcasm.irreps as casmirreps
import libcasm.xtal as xtal
import libcasm.xtal.prims as xtal_prims


def pretty(x):
    y = x.copy()
    y[np.abs(y) < 1e-8] = 0.0
    return y


def conventional_FCC_occ_symmetry_adapted_basis():
    # fmt: off
    return np.array([
        [ 0.,   0.,   0.,   0., ],
        [ 0.5,  0.5,  0.5,  0.5,],
        [ 0.,   0.,   0.,   0., ],
        [ 0.5,  0.5, -0.5, -0.5,],
        [ 0.,   0.,   0.,   0., ],
        [ 0.5, -0.5,  0.5, -0.5,],
        [ 0.,   0.,   0.,   0., ],
        [ 0.5, -0.5, -0.5,  0.5,],
    ])
    # fmt: on


def test_dof_space_analysis_2_generic(FCC_binary_prim):
    prim = casmconfig.Prim(FCC_binary_prim)
    T_dof_space = np.array(
        [  # conventional FCC cubic cell
            [-1, 1, 1],
            [1, -1, 1],
            [1, 1, -1],
        ],
        dtype=int,
    )
    supercell = casmconfig.Supercell(prim, T_dof_space)
    configuration = casmconfig.Configuration(
        supercell=supercell,
    )
    supercell_factor_group = casmconfig.make_invariant_subgroup(
        configuration=configuration,
    )
    symgroup = casmconfig.make_symgroup(supercell_factor_group)
    subset = casmgroup.Subset(group=symgroup)
    subset.all_subgroups(progress="none")
    subgroup_orbits = subset.all_subgroup_orbits()

    # construct occ DoFSpace with default basis
    dof_space = casmclex.DoFSpace(
        dof_key="occ",
        xtal_prim=FCC_binary_prim,
        transformation_matrix_to_super=T_dof_space,
    )

    matrix_rep = casmconfig.make_dof_space_rep(
        group=supercell_factor_group,
        dof_space=dof_space,
    )

    # Perform DoF space analysis
    irrep_decomposition = casmirreps.IrrepDecomposition(
        matrix_rep=matrix_rep,
        init_subspace=np.array(
            [
                [0.0, 0.0, 0.0, 0.0],
                [1.0, 0.0, 0.0, 0.0],
                [0.0, 0.0, 0.0, 0.0],
                [0.0, 1.0, 0.0, 0.0],
                [0.0, 0.0, 0.0, 0.0],
                [0.0, 0.0, 1.0, 0.0],
                [0.0, 0.0, 0.0, 0.0],
                [0.0, 0.0, 0.0, 1.0],
            ]
        ),
        subgroup_orbits=subgroup_orbits,
    )
    assert isinstance(irrep_decomposition, casmirreps.IrrepDecomposition)

    assert len(irrep_decomposition.irreps) == 2
    assert irrep_decomposition.symmetry_adapted_subspace.shape[0] == 8
    assert irrep_decomposition.symmetry_adapted_subspace.shape[1] == 4

    # print("Symmetry-adapted subspace:")
    # B = pretty(irrep_decomposition.symmetry_adapted_subspace)
    # for i in range(B.shape[1]):
    #     print(f"- {i}:", B[:, i].T)
    #     print()

    assert np.allclose(
        irrep_decomposition.symmetry_adapted_subspace,
        conventional_FCC_occ_symmetry_adapted_basis(),
    )

    sym_report = irrep_decomposition.make_symmetry_report(
        calc_wedges=True,
    )
    assert isinstance(sym_report, casmirreps.VectorSpaceSymReport)

    data = sym_report.to_dict()
    assert isinstance(data, dict)

    # print(xtal.pretty_json(data))


def test_BCC_2x2x2_GLstrain_disp_IrrepDecomposition():
    xtal_prim = xtal_prims.BCC(
        a=3.55,
        occ_dof=["A"],
        local_dof=[xtal.DoFSetBasis("disp")],
        global_dof=[xtal.DoFSetBasis("GLstrain")],
    )
    prim = casmconfig.Prim(xtal_prim=xtal_prim)
    T = np.eye(3, dtype=int) * 2
    supercell = casmconfig.Supercell(prim, T)

    ### Subgroup orbits ###
    symgroup = supercell.factor_group
    head_group = symgroup.head_group
    if head_group is None:
        head_group = symgroup
    indices = set(symgroup.head_group_index)
    subset = casmgroup.Subset(group=head_group, indices=indices)
    subset.all_subgroups(progress="none")
    subgroup_orbits = subset.all_subgroup_orbits()

    ### Irrep decomposition ###
    dof_key = "GLstrain"
    matrix_rep = prim.global_dof_matrix_rep(
        key=dof_key,
    )

    irrep_decomposition = casmirreps.IrrepDecomposition(
        matrix_rep=matrix_rep,
        subgroup_orbits=subgroup_orbits,
        # verbosity="standard",
    )

    irreps = irrep_decomposition.irreps
    assert len(irreps) == 3
    assert [irrep.irrep_dim for irrep in irreps] == [1, 2, 3]
    assert [irrep.irrep_type for irrep in irreps] == [0, 1, 2]

    B = irrep_decomposition.symmetry_adapted_subspace
    B = pretty(B)
    expected_B = np.array(
        [
            [0.57735027, 0.81649658, 0.0, 0.0, 0.0, 0.0],
            [0.57735027, -0.40824829, 0.70710678, 0.0, 0.0, 0.0],
            [0.57735027, -0.40824829, -0.70710678, 0.0, 0.0, 0.0],
            [0.0, 0.0, 0.0, 1.0, 0.0, 0.0],
            [0.0, 0.0, 0.0, 0.0, 1.0, 0.0],
            [0.0, 0.0, 0.0, 0.0, 0.0, 1.0],
        ]
    )
    assert np.allclose(B, expected_B)

    ### Irrep decomposition ###
    configuration = casmconfig.Configuration(
        supercell=supercell,
    )
    supercell_factor_group = casmconfig.make_invariant_subgroup(
        configuration=configuration,
    )
    symgroup = casmconfig.make_symgroup(supercell_factor_group)

    ### Subgroup orbits ###
    head_group = symgroup.head_group
    if head_group is None:
        head_group = symgroup
    indices = set(symgroup.head_group_index)
    subset = casmgroup.Subset(group=head_group, indices=indices)
    subset.all_subgroups(progress="none")
    subgroup_orbits = subset.all_subgroup_orbits()

    # construct disp DoFSpace with default basis
    dof_space = casmclex.DoFSpace(
        dof_key="disp",
        xtal_prim=prim.xtal_prim,
        transformation_matrix_to_super=supercell.transformation_matrix_to_super,
    )

    matrix_rep = casmconfig.make_dof_space_rep(
        group=supercell_factor_group,
        dof_space=dof_space,
    )

    # Perform DoF space analysis
    irrep_decomposition = casmirreps.IrrepDecomposition(
        matrix_rep=matrix_rep,
        subgroup_orbits=subgroup_orbits,
        # verbosity="standard",
    )

    irreps = irrep_decomposition.irreps
    assert len(irreps) == 5
    assert [irrep.irrep_dim for irrep in irreps] == [3, 3, 6, 6, 6]
    assert [irrep.irrep_type for irrep in irreps] == [0, 1, 2, 3, 4]

    B = irrep_decomposition.symmetry_adapted_subspace
    B = pretty(B)

    file = (
        pathlib.Path(__file__).parent
        / "data"
        / "test_BCC_2x2x2_disp_IrrepDecomposition.json"
    )
    # with open(file, "w") as f:
    #     f.write(xtal.pretty_json(B.T.tolist()))
    with open(file, "r") as f:
        expected_B = np.array(json.load(f)).T

    assert np.allclose(B, expected_B)
