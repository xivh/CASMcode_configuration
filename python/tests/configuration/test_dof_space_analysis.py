import numpy as np

import libcasm.clexulator as casmclex
import libcasm.configuration as casmconfig
import libcasm.irreps as casmirreps


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
        # calc_wedges=False,
    )

    symmetry_adapted_dof_space = results.symmetry_adapted_dof_space
    sym_report = results.symmetry_report

    assert len(sym_report.irreps) == 2
    assert sym_report.symmetry_adapted_subspace.shape[0] == 8
    assert sym_report.symmetry_adapted_subspace.shape[1] == 4

    assert np.allclose(
        symmetry_adapted_dof_space.basis, sym_report.symmetry_adapted_subspace
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
