import numpy as np

import libcasm.clexulator as casmclex
import libcasm.configuration as casmconfig
import libcasm.irreps as casmirreps


def conventional_FCC_occ_symmetry_adapted_basis():
    # fmt: off
    return np.array([
        [ 0.,   0.,   0.,   0., ],
        [-0.5, -0.5, -0.5, -0.5,],
        [ 0.,   0.,   0.,   0., ],
        [-0.5, -0.5,  0.5,  0.5,],
        [ 0.,   0.,   0.,   0., ],
        [-0.5,  0.5, -0.5,  0.5,],
        [ 0.,   0.,   0.,   0., ],
        [-0.5,  0.5,  0.5, -0.5,],
    ])
    # fmt: on


def test_MatrixRepGroup_1(FCC_binary_prim):
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

    group = casmirreps.MatrixRepGroup(
        elements=matrix_rep,
    )

    assert isinstance(group, casmirreps.MatrixRepGroup)
    assert len(group.elements) == 192
    assert len(group.multiplication_table) == 192

    # Perform DoF space analysis
    irrep_decomposition = casmirreps.IrrepDecomposition(
        matrix_rep=group.elements,
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
    )
    assert isinstance(irrep_decomposition, casmirreps.IrrepDecomposition)

    assert len(irrep_decomposition.irreps) == 2
    assert irrep_decomposition.symmetry_adapted_subspace.shape[0] == 8
    assert irrep_decomposition.symmetry_adapted_subspace.shape[1] == 4

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
