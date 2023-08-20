import numpy as np

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
    assert supercell.n_sites() == 4
    configuration = config.Configuration(supercell)
    invariant_subgroup = config.make_invariant_subgroup(configuration)
    site_indices = set(range(0, supercell.n_sites()))

    occ_rep = config.make_local_dof_matrix_rep(invariant_subgroup, "occ", site_indices)
    # print(occ_rep)
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
    assert supercell.n_sites() == 2
    configuration = config.Configuration(supercell)
    invariant_subgroup = config.make_invariant_subgroup(configuration)
    site_indices = set(range(0, supercell.n_sites()))

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
