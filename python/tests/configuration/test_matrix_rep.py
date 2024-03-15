import numpy as np

import libcasm.clexulator as casmclex
import libcasm.configuration as config


def test_conventional_fcc_matrix_rep(FCC_binary_GLstrain_disp_prim):
    xtal_prim = FCC_binary_GLstrain_disp_prim
    prim = config.Prim(xtal_prim)
    T = np.array(
        [
            [-1, 1, 1],
            [1, -1, 1],
            [1, 1, -1],
        ]
    )
    supercell = config.make_canonical_supercell(config.Supercell(prim, T))
    assert supercell.n_sites == 4
    configuration = config.Configuration(supercell)
    invariant_subgroup = config.make_invariant_subgroup(configuration)
    site_indices = set(range(0, supercell.n_sites))

    occ_rep = config.make_local_dof_matrix_rep(invariant_subgroup, "occ", site_indices)
    assert len(occ_rep) == 4 * 48

    disp_rep = config.make_local_dof_matrix_rep(
        invariant_subgroup, "disp", site_indices
    )
    # print(disp_rep)
    assert len(disp_rep) == 4 * 48

    GLstrain_rep = config.make_global_dof_matrix_rep(invariant_subgroup, "GLstrain")
    # print(GLstrain_rep)
    assert len(GLstrain_rep) == 48


def test_conventional_bcc_matrix_rep(BCC_binary_GLstrain_disp_prim):
    xtal_prim = BCC_binary_GLstrain_disp_prim
    prim = config.Prim(xtal_prim)
    T = np.array(
        [
            [0, 1, 1],
            [1, 0, 1],
            [1, 1, 0],
        ]
    )
    supercell = config.make_canonical_supercell(config.Supercell(prim, T))
    assert supercell.n_sites == 2
    configuration = config.Configuration(supercell)
    invariant_subgroup = config.make_invariant_subgroup(configuration)
    site_indices = set(range(0, supercell.n_sites))

    occ_rep = config.make_local_dof_matrix_rep(invariant_subgroup, "occ", site_indices)
    # print(occ_rep)
    assert len(occ_rep) == 2 * 48

    disp_rep = config.make_local_dof_matrix_rep(
        invariant_subgroup, "disp", site_indices
    )
    # print(disp_rep)
    assert len(disp_rep) == 2 * 48

    GLstrain_rep = config.make_global_dof_matrix_rep(invariant_subgroup, "GLstrain")
    # print(GLstrain_rep)
    assert len(GLstrain_rep) == 48


def test_Hstrain_dof_space_rep_1(FCC_binary_Hstrain_disp_prim):
    xtal_prim = FCC_binary_Hstrain_disp_prim
    prim = config.Prim(xtal_prim)

    T_prim = np.eye(3, dtype=int)
    supercell = config.make_canonical_supercell(config.Supercell(prim, T_prim))
    configuration = config.Configuration(
        supercell=supercell,
    )
    supercell_factor_group = config.make_invariant_subgroup(
        configuration=configuration,
    )

    # construct Hstrain DoFSpace with default basis
    dof_space = casmclex.DoFSpace(
        dof_key="Hstrain",
        xtal_prim=xtal_prim,
    )
    matrix_rep = config.make_dof_space_rep(
        group=supercell_factor_group,
        dof_space=dof_space,
    )
    assert len(matrix_rep) == len(supercell_factor_group)

    Hstrain = np.array([0.01, 0.02, 0.03, 0.04, 0.05, 0.06])
    configuration.set_global_dof_values("Hstrain", Hstrain)
    assert np.allclose(configuration.global_dof_values("Hstrain"), Hstrain)
    assert np.isclose(configuration.global_dof_values("Hstrain")[0], Hstrain[0])

    # construct a OrderParameter calculator
    order_parameter = casmclex.OrderParameter(
        dof_space=dof_space,
    )
    dof_values = configuration.dof_values
    order_parameter.update(
        transformation_matrix_to_super=supercell.transformation_matrix_to_super,
        site_index_converter=supercell.site_index_converter,
        config_dof_values=dof_values,
    )

    # get values
    eta = order_parameter.value()
    print(eta)
    assert np.allclose(eta, Hstrain)

    # check transforming values, and compare against dof_space matrix_rep
    for i, M in enumerate(matrix_rep):
        transformed_eta = M @ eta
        transformed_config = supercell_factor_group[i] * configuration
        assert np.allclose(
            transformed_eta, transformed_config.global_dof_values("Hstrain")
        )


def test_disp_dof_space_rep_1(FCC_binary_Hstrain_disp_prim):
    """Apply SupercellSymOp"""
    xtal_prim = FCC_binary_Hstrain_disp_prim
    prim = config.Prim(xtal_prim)

    T_conventional = np.array(
        [  # conventional FCC cubic cell
            [-1, 1, 1],
            [1, -1, 1],
            [1, 1, -1],
        ],
        dtype=int,
    )

    supercell = config.make_canonical_supercell(config.Supercell(prim, T_conventional))
    configuration = config.Configuration(
        supercell=supercell,
    )
    supercell_factor_group = config.make_invariant_subgroup(
        configuration=configuration,
    )

    # construct disp DoFSpace with default basis
    dof_space = casmclex.DoFSpace(
        dof_key="disp",
        xtal_prim=xtal_prim,
        transformation_matrix_to_super=T_conventional,
    )
    matrix_rep = config.make_dof_space_rep(
        group=supercell_factor_group,
        dof_space=dof_space,
    )
    assert len(matrix_rep) == len(supercell_factor_group)

    # occupation = np.array([1, 1, 0, 0], dtype=int)
    # configuration.set_occupation(occupation)

    disp = np.array(
        [
            [0.01, 0.05, 0.09],
            [0.02, 0.06, 0.10],
            [0.03, 0.07, 0.11],
            [0.04, 0.08, 0.12],
        ]
    ).transpose()
    configuration.set_local_dof_values("disp", disp)
    assert np.allclose(configuration.local_dof_values("disp"), disp)
    assert np.isclose(configuration.local_dof_values("disp")[0, 0], disp[0, 0])

    # construct a OrderParameter calculator
    order_parameter = casmclex.OrderParameter(
        dof_space=dof_space,
    )
    dof_values = configuration.dof_values
    order_parameter.update(
        transformation_matrix_to_super=supercell.transformation_matrix_to_super,
        site_index_converter=supercell.site_index_converter,
        config_dof_values=dof_values,
    )

    # get values
    eta = order_parameter.value()
    print(eta)
    assert np.allclose(eta, disp.reshape((12,), order="F"))

    # check transforming values, and compare against dof_space matrix_rep
    for i, M in enumerate(matrix_rep):
        transformed_eta = M @ eta
        transformed_config = supercell_factor_group[i] * configuration
        transformed_disp = transformed_config.local_dof_values("disp")

        # print(f"~~~ i: {i} ~~~")
        # print(M)
        # print(eta)
        # print(transformed_disp)
        # print(transformed_disp.reshape((12,), order="F"))
        # print(transformed_eta)
        # print()
        assert np.allclose(transformed_eta, transformed_disp.reshape((12,), order="F"))

    # assert False
