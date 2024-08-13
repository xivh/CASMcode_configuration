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


def conventional_FCC_occ_irrep_1_wedge_axes():
    # fmt: off
    return np.array([
        [ 0., ],
        [-0.5,],
        [ 0., ],
        [-0.5,],
        [ 0., ],
        [-0.5,],
        [ 0., ],
        [-0.5,],
    ])
    # fmt: off


def conventional_FCC_occ_irrep_2_wedge_axes():
    # fmt: off
    return np.array([
        [ 0.,         0.,         0.,        ],
        [-0.8660254, -0.5,       -0.28867513,],
        [ 0.,         0.,         0.        ,],
        [ 0.28867513, 0.5,       -0.28867513,],
        [ 0.,         0.,         0.        ,],
        [ 0.28867513, 0.5,        0.8660254 ,],
        [ 0.,         0.,         0.        ,],
        [ 0.28867513,-0.5,       -0.28867513,],
    ])
    # fmt: on


def test_dof_space_analysis_1(FCC_binary_prim):
    prim = casmconfig.Prim(FCC_binary_prim)
    T_dof_space = np.array([[1, 0, 0], [0, 1, 0], [0, 0, 1]], dtype=int)

    # construct occ DoFSpace with default basis
    dof_space = casmclex.DoFSpace(
        dof_key="occ",
        xtal_prim=FCC_binary_prim,
        transformation_matrix_to_super=T_dof_space,
    )

    # Perform DoF space analysis
    results = casmconfig.dof_space_analysis(
        dof_space=dof_space,
        prim=prim,
        # configuration=None,
        # exclude_homogeneous_modes=None,
        # include_default_occ_modes=False,
        # calc_wedges=False,
    )

    symmetry_adapted_dof_space = results.symmetry_adapted_dof_space
    assert isinstance(symmetry_adapted_dof_space, casmclex.DoFSpace)
    sym_report = results.symmetry_report
    assert isinstance(sym_report, casmirreps.VectorSpaceSymReport)

    assert len(sym_report.irreps) == 1
    assert sym_report.symmetry_adapted_subspace.shape[0] == 2
    assert sym_report.symmetry_adapted_subspace.shape[1] == 1

    assert np.allclose(
        symmetry_adapted_dof_space.basis, sym_report.symmetry_adapted_subspace
    )

    data = sym_report.to_dict()
    assert isinstance(data, dict)


def test_dof_space_analysis_1_generic(FCC_binary_prim):
    prim = casmconfig.Prim(FCC_binary_prim)
    T_dof_space = np.array([[1, 0, 0], [0, 1, 0], [0, 0, 1]], dtype=int)
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

    # Perform DoF space analysis
    irrep_decomposition = casmirreps.IrrepDecomposition(
        matrix_rep=matrix_rep,
        init_subspace=np.array(
            [
                [0.0],
                [1.0],
            ]
        ),
    )
    assert isinstance(irrep_decomposition, casmirreps.IrrepDecomposition)

    assert len(irrep_decomposition.irreps) == 1
    assert irrep_decomposition.symmetry_adapted_subspace.shape[0] == 2
    assert irrep_decomposition.symmetry_adapted_subspace.shape[1] == 1

    sym_report = irrep_decomposition.make_symmetry_report()
    assert isinstance(sym_report, casmirreps.VectorSpaceSymReport)

    data = sym_report.to_dict()
    assert isinstance(data, dict)


def test_dof_space_analysis_2(FCC_binary_prim):
    prim = casmconfig.Prim(FCC_binary_prim)
    T_dof_space = np.array(
        [  # conventional FCC cubic cell
            [-1, 1, 1],
            [1, -1, 1],
            [1, 1, -1],
        ],
        dtype=int,
    )

    # construct occ DoFSpace with default basis
    dof_space = casmclex.DoFSpace(
        dof_key="occ",
        xtal_prim=FCC_binary_prim,
        transformation_matrix_to_super=T_dof_space,
    )

    # Perform DoF space analysis
    results = casmconfig.dof_space_analysis(
        dof_space=dof_space,
        prim=prim,
        # configuration=None,
        # exclude_homogeneous_modes=None,
        # include_default_occ_modes=False,
        # sublattice_index_to_default_occ=None,
        # site_index_to_default_occ=None,
        calc_wedges=True,
    )

    symmetry_adapted_dof_space = results.symmetry_adapted_dof_space
    sym_report = results.symmetry_report

    assert len(sym_report.irreps) == 2
    assert sym_report.symmetry_adapted_subspace.shape[0] == 8
    assert sym_report.symmetry_adapted_subspace.shape[1] == 4

    assert np.allclose(
        symmetry_adapted_dof_space.basis, sym_report.symmetry_adapted_subspace
    )

    assert np.allclose(
        symmetry_adapted_dof_space.basis, conventional_FCC_occ_symmetry_adapted_basis()
    )

    data = sym_report.to_dict()
    assert isinstance(data, dict)

    data = results.to_dict()
    assert isinstance(data, dict)


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


def test_dof_space_analysis_2a(FCC_binary_prim):
    prim = casmconfig.Prim(FCC_binary_prim)
    T_dof_space = np.array(
        [  # conventional FCC cubic cell
            [-1, 1, 1],
            [1, -1, 1],
            [1, 1, -1],
        ],
        dtype=int,
    )

    # construct occ DoFSpace with default basis
    dof_space = casmclex.DoFSpace(
        dof_key="occ",
        xtal_prim=FCC_binary_prim,
        transformation_matrix_to_super=T_dof_space,
    )

    site_index_to_default_occ = {0: 0, 1: 0, 2: 0, 3: 1}

    # Perform DoF space analysis
    results = casmconfig.dof_space_analysis(
        dof_space=dof_space,
        prim=prim,
        # configuration=None,
        # exclude_homogeneous_modes=None,
        # include_default_occ_modes=False,
        # sublattice_index_to_default_occ=None,
        site_index_to_default_occ=site_index_to_default_occ,
        # calc_wedges=False,
    )

    symmetry_adapted_dof_space = results.symmetry_adapted_dof_space
    sym_report = results.symmetry_report
    # print("basis:\n", symmetry_adapted_dof_space.basis)

    assert len(sym_report.irreps) == 4
    assert sym_report.symmetry_adapted_subspace.shape[0] == 8
    assert sym_report.symmetry_adapted_subspace.shape[1] == 8

    assert np.allclose(
        symmetry_adapted_dof_space.basis, sym_report.symmetry_adapted_subspace
    )


def test_dof_space_analysis_2b(FCC_binary_prim):
    prim = casmconfig.Prim(FCC_binary_prim)
    T_dof_space = np.array(
        [  # conventional FCC cubic cell
            [-1, 1, 1],
            [1, -1, 1],
            [1, 1, -1],
        ],
        dtype=int,
    )

    # construct occ DoFSpace with default basis
    dof_space = casmclex.DoFSpace(
        dof_key="occ",
        xtal_prim=FCC_binary_prim,
        transformation_matrix_to_super=T_dof_space,
    )

    sublattice_index_to_default_occ = {0: 1}

    # Perform DoF space analysis
    results = casmconfig.dof_space_analysis(
        dof_space=dof_space,
        prim=prim,
        # configuration=None,
        # exclude_homogeneous_modes=None,
        # include_default_occ_modes=False,
        sublattice_index_to_default_occ=sublattice_index_to_default_occ,
        # site_index_to_default_occ=None,
        # calc_wedges=False,
    )

    symmetry_adapted_dof_space = results.symmetry_adapted_dof_space
    sym_report = results.symmetry_report
    # print("basis:\n", symmetry_adapted_dof_space.basis)

    assert len(sym_report.irreps) == 2
    assert sym_report.symmetry_adapted_subspace.shape[0] == 8
    assert sym_report.symmetry_adapted_subspace.shape[1] == 4

    assert np.allclose(
        symmetry_adapted_dof_space.basis, sym_report.symmetry_adapted_subspace
    )


def test_dof_space_analysis_2c(FCC_binary_prim):
    prim = casmconfig.Prim(FCC_binary_prim)
    T_dof_space = np.array(
        [  # conventional FCC cubic cell
            [-1, 1, 1],
            [1, -1, 1],
            [1, 1, -1],
        ],
        dtype=int,
    )

    # construct occ DoFSpace with default basis
    dof_space = casmclex.DoFSpace(
        dof_key="occ",
        xtal_prim=FCC_binary_prim,
        transformation_matrix_to_super=T_dof_space,
    )

    # Perform DoF space analysis
    results = casmconfig.dof_space_analysis(
        dof_space=dof_space,
        prim=prim,
        # configuration=None,
        # exclude_homogeneous_modes=None,
        # include_default_occ_modes=False,
        # sublattice_index_to_default_occ=None,
        # site_index_to_default_occ=None,
        calc_wedges=True,
    )

    symmetry_adapted_dof_space = results.symmetry_adapted_dof_space
    sym_report = results.symmetry_report

    assert len(sym_report.irreps) == 2
    assert sym_report.symmetry_adapted_subspace.shape[0] == 8
    assert sym_report.symmetry_adapted_subspace.shape[1] == 4

    assert np.allclose(
        symmetry_adapted_dof_space.basis, sym_report.symmetry_adapted_subspace
    )

    assert np.allclose(
        symmetry_adapted_dof_space.basis, conventional_FCC_occ_symmetry_adapted_basis()
    )

    assert len(sym_report.irreducible_wedge) == 1

    assert len(sym_report.irrep_names) == 2
    assert sym_report.irrep_names == ["irrep_1_1", "irrep_2_1"]

    assert len(sym_report.irrep_axes_indices) == 2
    assert sym_report.irrep_axes_indices == [[0], [1, 2, 3]]

    assert len(sym_report.irrep_wedge_axes) == 2
    assert np.allclose(
        sym_report.irrep_wedge_axes[0], conventional_FCC_occ_irrep_1_wedge_axes()
    )
    assert np.allclose(
        sym_report.irrep_wedge_axes[1], conventional_FCC_occ_irrep_2_wedge_axes()
    )


def test_dof_space_analysis_3(FCC_binary_occ_fix_corner_prim):
    prim = casmconfig.Prim(FCC_binary_occ_fix_corner_prim)
    T_dof_space = np.array([[1, 0, 0], [0, 1, 0], [0, 0, 1]], dtype=int)

    # construct occ DoFSpace with default basis
    dof_space = casmclex.DoFSpace(
        dof_key="occ",
        xtal_prim=FCC_binary_occ_fix_corner_prim,
        transformation_matrix_to_super=T_dof_space,
    )

    # Perform DoF space analysis
    results = casmconfig.dof_space_analysis(
        dof_space=dof_space,
        prim=prim,
        # configuration=None,
        # exclude_homogeneous_modes=None,
        # include_default_occ_modes=False,
        # calc_wedges=False,
    )

    symmetry_adapted_dof_space = results.symmetry_adapted_dof_space
    sym_report = results.symmetry_report

    assert len(sym_report.irreps) == 2
    assert sym_report.symmetry_adapted_subspace.shape[0] == 7
    assert sym_report.symmetry_adapted_subspace.shape[1] == 3

    assert np.allclose(
        symmetry_adapted_dof_space.basis, sym_report.symmetry_adapted_subspace
    )


def test_dof_space_analysis_4(FCC_binary_GLstrain_disp_prim):
    prim = casmconfig.Prim(FCC_binary_GLstrain_disp_prim)
    T_dof_space = np.array(
        [  # conventional FCC cubic cell
            [-1, 1, 1],
            [1, -1, 1],
            [1, 1, -1],
        ],
        dtype=int,
    )

    # construct disp DoFSpace with default basis
    dof_space = casmclex.DoFSpace(
        dof_key="disp",
        xtal_prim=FCC_binary_GLstrain_disp_prim,
        transformation_matrix_to_super=T_dof_space,
    )

    # Perform DoF space analysis
    results = casmconfig.dof_space_analysis(
        dof_space=dof_space,
        prim=prim,
        # configuration=None,
        # exclude_homogeneous_modes=None,
        # include_default_occ_modes=False,
        # calc_wedges=False,
    )

    symmetry_adapted_dof_space = results.symmetry_adapted_dof_space
    sym_report = results.symmetry_report

    assert len(sym_report.irreps) == 2
    assert sym_report.symmetry_adapted_subspace.shape[0] == 12
    assert sym_report.symmetry_adapted_subspace.shape[1] == 9

    assert np.allclose(
        symmetry_adapted_dof_space.basis, sym_report.symmetry_adapted_subspace
    )


def test_dof_space_analysis_5(prim_ABC2):
    prim = casmconfig.Prim(prim_ABC2)
    T_dof_space = np.array([[1, 0, 0], [0, 1, 0], [0, 0, 1]], dtype=int)

    # construct disp DoFSpace with default basis
    dof_space = casmclex.DoFSpace(
        dof_key="disp", xtal_prim=prim_ABC2, transformation_matrix_to_super=T_dof_space
    )

    # Perform DoF space analysis
    results = casmconfig.dof_space_analysis(
        dof_space=dof_space,
        prim=prim,
        # configuration=None,
        # exclude_homogeneous_modes=None,
        # include_default_occ_modes=False,
        # calc_wedges=False,
    )

    symmetry_adapted_dof_space = results.symmetry_adapted_dof_space
    sym_report = results.symmetry_report

    assert len(sym_report.irreps) == 15
    assert sym_report.symmetry_adapted_subspace.shape[0] == 18
    assert sym_report.symmetry_adapted_subspace.shape[1] == 15

    assert np.allclose(
        symmetry_adapted_dof_space.basis, sym_report.symmetry_adapted_subspace
    )


def test_dof_space_analysis_6(FCC_binary_disp_fix_corner_prim):
    prim = casmconfig.Prim(FCC_binary_disp_fix_corner_prim)
    T_dof_space = np.array([[1, 0, 0], [0, 1, 0], [0, 0, 1]], dtype=int)

    # construct disp DoFSpace with default basis
    dof_space = casmclex.DoFSpace(
        dof_key="disp",
        xtal_prim=FCC_binary_disp_fix_corner_prim,
        transformation_matrix_to_super=T_dof_space,
    )

    # Perform DoF space analysis
    results = casmconfig.dof_space_analysis(
        dof_space=dof_space,
        prim=prim,
        # configuration=None,
        # exclude_homogeneous_modes=None,
        # include_default_occ_modes=False,
        # calc_wedges=False,
    )

    symmetry_adapted_dof_space = results.symmetry_adapted_dof_space
    sym_report = results.symmetry_report

    assert len(sym_report.irreps) == 3
    assert sym_report.symmetry_adapted_subspace.shape[0] == 9
    assert sym_report.symmetry_adapted_subspace.shape[1] == 9

    assert np.allclose(
        symmetry_adapted_dof_space.basis, sym_report.symmetry_adapted_subspace
    )
