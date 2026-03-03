import json
import pathlib

import numpy as np

import libcasm.clexulator as casmclex
import libcasm.configuration as casmconfig
import libcasm.configuration.io.spglib as spglib_io
import libcasm.group as casmgroup
import libcasm.irreps as casmirreps
import libcasm.xtal as xtal


def conventional_FCC_occ_symmetry_adapted_basis():
    # fmt: off
    return np.array([
        [ 0.,  0.,   0.,   0., ],
        [0.5,  0.5,  0.5,  0.5,],
        [ 0.,  0.,   0.,   0., ],
        [0.5,  0.5, -0.5, -0.5,],
        [ 0.,  0.,   0.,   0., ],
        [0.5, -0.5,  0.5, -0.5,],
        [ 0.,  0.,   0.,   0., ],
        [0.5, -0.5, -0.5,  0.5,],
    ])
    # fmt: on


def conventional_FCC_occ_irrep_1_wedge_axes():
    # fmt: off
    return np.array([
        [ 0., ],
        [ 0.5,],
        [ 0., ],
        [ 0.5,],
        [ 0., ],
        [ 0.5,],
        [ 0., ],
        [ 0.5,],
    ])
    # fmt: off


def conventional_FCC_occ_irrep_2_wedge_axes():
    # fmt: off
    return np.array([
        [ 0.,          0.,   0.,        ],
        [0.8660254,    0.5,  0.28867513,],
        [ 0.,          0.,   0.        ,],
        [-0.28867513, -0.5,  0.28867513,],
        [ 0.,          0.,   0.        ,],
        [-0.28867513, -0.5, -0.8660254 ,],
        [ 0.,          0.,   0.        ,],
        [-0.28867513,  0.5,  0.28867513,],
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
        verbosity=None,
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

    # Check if basis.T @ basis is close to identity
    basis = symmetry_adapted_dof_space.basis
    identity_approx = basis.T @ basis
    # print("Diagonal of Basis.T @ Basis:\n", np.diag(identity_approx))
    # print("Basis.T @ Basis:\n", clean(identity_approx))
    assert np.allclose(identity_approx, np.eye(basis.shape[1]), atol=1e-5)


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
        verbosity="verbose",
    )
    assert isinstance(irrep_decomposition, casmirreps.IrrepDecomposition)

    assert len(irrep_decomposition.irreps) == 1
    assert irrep_decomposition.symmetry_adapted_subspace.shape[0] == 2
    assert irrep_decomposition.symmetry_adapted_subspace.shape[1] == 1

    sym_report = irrep_decomposition.make_symmetry_report()
    assert isinstance(sym_report, casmirreps.VectorSpaceSymReport)

    data = sym_report.to_dict()
    assert isinstance(data, dict)

    # Check if basis.T @ basis is close to identity
    basis = irrep_decomposition.symmetry_adapted_subspace
    identity_approx = basis.T @ basis
    # print("Diagonal of Basis.T @ Basis:\n", np.diag(identity_approx))
    # print("Basis.T @ Basis:\n", clean(identity_approx))
    assert np.allclose(identity_approx, np.eye(basis.shape[1]), atol=1e-5)


def plot_irrep_axes(irrep: casmirreps.IrrepInfo, index: int):
    from bokeh.layouts import gridplot
    from bokeh.models import Title
    from bokeh.plotting import figure, show

    i = index

    title = f"Irrep {i} Basis"
    B = irrep.trans_mat.real

    plots = []
    for j in range(B.shape[0]):
        # p = figure(width=800, height=300, title=f"Basis Vector {j}")
        # p.add_layout(Title(text=title, align="center"), "above")
        # x = np.arange(B.shape[1])
        # y = B[j, :]
        # p.vbar(x=x, top=y, width=0.5)
        # plots.append(p)

        # Each row of B corresponds to (x1, y1, z1, x2, y2, z2, ..., xn, yn, zn)
        # where (xi, yi, zi) are the components of the displacement vector for atom i.
        # Let's plot the x, y, z components separately for clarity:
        p_xyz = figure(width=400, height=200, title=f"Basis Vector {j} Components")
        p_xyz.add_layout(Title(text=title, align="center"), "above")
        x_atoms = np.arange(B.shape[1] // 3)
        y_x = B[j, 0::3]
        y_y = B[j, 1::3]
        y_z = B[j, 2::3]
        p_xyz.vbar(x=x_atoms - 0.2, top=y_x, width=0.2, color="red", legend_label="x")
        p_xyz.vbar(x=x_atoms, top=y_y, width=0.2, color="green", legend_label="y")
        p_xyz.vbar(x=x_atoms + 0.2, top=y_z, width=0.2, color="blue", legend_label="z")
        # p_xyz.legend.location = "top_right"
        # Move the auto-created legend outside the plot area (to the right)
        if p_xyz.legend:
            legend = p_xyz.legend[0]
            p_xyz.add_layout(legend, "right")
            legend.orientation = "vertical"
            legend.label_text_font_size = "10pt"

        plots.append(p_xyz)

    grid = gridplot(plots, ncols=3)
    show(grid)


def plot_irrep_directions_orbit(
    directions: list[np.array],
    irrep_index: int,
    orbit_index: int,
):
    from bokeh.layouts import gridplot
    from bokeh.models import Title
    from bokeh.plotting import figure, show

    i = irrep_index
    k = orbit_index

    title = f"Irrep {i}, Direction orbit {k}"

    plots = []
    for i_dir, d in enumerate(directions):
        p_dir = figure(
            width=400,
            height=200,
            title=f"{title}, Direction {i_dir}",
        )
        p_dir.add_layout(Title(text=title, align="center"), "above")
        x_atoms = np.arange(d.shape[0] // 3)
        y_x = d[0::3]
        y_y = d[1::3]
        y_z = d[2::3]
        p_dir.vbar(x=x_atoms - 0.2, top=y_x, width=0.2, color="red", legend_label="x")
        p_dir.vbar(x=x_atoms, top=y_y, width=0.2, color="green", legend_label="y")
        p_dir.vbar(x=x_atoms + 0.2, top=y_z, width=0.2, color="blue", legend_label="z")
        # p_dir.legend.location = "top_right"
        # Move the auto-created legend outside the plot area (to the right)
        if p_dir.legend:
            legend = p_dir.legend[0]
            p_dir.add_layout(legend, "right")
            legend.orientation = "vertical"
            legend.label_text_font_size = "10pt"

        plots.append(p_dir)

    grid = gridplot(plots, ncols=3)
    show(grid)


def print_report(report):
    for i, irrep in enumerate(report.irreps):
        with open(f"irrep.{i}.json", "w") as f:
            f.write(xtal.pretty_json(irrep.to_dict()))

    irreps_in = []
    i = 0
    path = pathlib.Path(f"irrep.{i}.json")
    while path.exists():
        with open(path, "r") as f:
            data = json.load(f)
            irrep = casmirreps.IrrepInfo.from_dict(data)
            irreps_in.append(irrep)
        i += 1
        path = pathlib.Path(f"irrep.{i}.json")

    irreps = irreps_in

    for i, irrep in enumerate(irreps):
        dirs_mult = [len(directions) for directions in irrep.directions]
        print(f"Irrep {i}: dim={irrep.irrep_dim}, dirs_mult={dirs_mult}")
    print()

    # selected_irrep = 1
    #
    # for i, irrep in enumerate(irreps):
    #     # if i != selected_irrep:
    #     #     continue
    #
    #     plot_irrep_axes(irrep=irrep, index=i)
    #
    #     # for k, directions in enumerate(irrep.directions):
    #     #     plot_irrep_directions_orbit(
    #     #         directions=directions,
    #     #         irrep_index=i,
    #     #         orbit_index=k,
    #     #     )

    # for i, irrep in enumerate(report.irreps):
    #     dirs_mult = [len(directions) for directions in irrep.directions]
    #     print(f"-- Irrep {i}: dim={irrep.irrep_dim}, dirs_mult={dirs_mult} --")
    #     for j, directions in enumerate(irrep.directions):
    #         print(f"- Direction Orbit {j}:")
    #         for k, d in enumerate(directions):
    #             print(f"- - {k}: {d.tolist()}")
    #     print()

    print("-- Summary --")
    print("Symmetry adapted space shape=", report.symmetry_adapted_subspace.shape)

    print()


def make_maximal_subgroups(supercell_symops: list[casmconfig.SupercellSymOp]):
    """Make maximal subgroups from supercell symmetry operations."""

    no_lattice_translation_group = []
    only_lattice_translation_group = []
    for op in supercell_symops:
        if op.translation_index() == 0:
            no_lattice_translation_group.append(op)
        if op.prim_factor_group_index() == 0:
            only_lattice_translation_group.append(op)
    print(f"# no_lattice_translation_group: {len(no_lattice_translation_group)}")
    print(f"# only_lattice_translation_group: {len(only_lattice_translation_group)}")


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
        # symmetrization="fast",
        calc_wedges=False,
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

    # Check if basis.T @ basis is close to identity
    basis = symmetry_adapted_dof_space.basis
    identity_approx = basis.T @ basis
    # print("Diagonal of Basis.T @ Basis:\n", np.diag(identity_approx))
    # print("Basis.T @ Basis:\n", clean(identity_approx))
    assert np.allclose(identity_approx, np.eye(basis.shape[1]), atol=1e-5)


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
    subset.all_subgroups()
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
        # symmetrization="fast",
        subgroup_orbits=subgroup_orbits,
        verbosity="verbose",
    )
    assert isinstance(irrep_decomposition, casmirreps.IrrepDecomposition)

    assert len(irrep_decomposition.irreps) == 2
    assert irrep_decomposition.symmetry_adapted_subspace.shape[0] == 8
    assert irrep_decomposition.symmetry_adapted_subspace.shape[1] == 4

    print("Symmetry adapted basis:\n", irrep_decomposition.symmetry_adapted_subspace)
    print(
        "Conventional FCC symmetry adapted basis:\n",
        conventional_FCC_occ_symmetry_adapted_basis(),
    )

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

    # Check if basis.T @ basis is close to identity
    basis = irrep_decomposition.symmetry_adapted_subspace
    identity_approx = basis.T @ basis
    # print("Diagonal of Basis.T @ Basis:\n", np.diag(identity_approx))
    # print("Basis.T @ Basis:\n", clean(identity_approx))
    assert np.allclose(identity_approx, np.eye(basis.shape[1]), atol=1e-5)


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
        symmetrization="fast",
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

    # Check if basis.T @ basis is close to identity
    basis = symmetry_adapted_dof_space.basis
    identity_approx = basis.T @ basis
    # print("Diagonal of Basis.T @ Basis:\n", np.diag(identity_approx))
    # print("Basis.T @ Basis:\n", clean(identity_approx))
    assert np.allclose(identity_approx, np.eye(basis.shape[1]), atol=1e-5)


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
        symmetrization="fast",
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

    # Check if basis.T @ basis is close to identity
    basis = symmetry_adapted_dof_space.basis
    identity_approx = basis.T @ basis
    # print("Diagonal of Basis.T @ Basis:\n", np.diag(identity_approx))
    # print("Basis.T @ Basis:\n", clean(identity_approx))
    assert np.allclose(identity_approx, np.eye(basis.shape[1]), atol=1e-5)


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
        symmetrization="fast",
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

    # Check if basis.T @ basis is close to identity
    basis = symmetry_adapted_dof_space.basis
    identity_approx = basis.T @ basis
    # print("Diagonal of Basis.T @ Basis:\n", np.diag(identity_approx))
    # print("Basis.T @ Basis:\n", clean(identity_approx))
    assert np.allclose(identity_approx, np.eye(basis.shape[1]), atol=1e-5)


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

    # Check if basis.T @ basis is close to identity
    basis = symmetry_adapted_dof_space.basis
    identity_approx = basis.T @ basis
    # print("Diagonal of Basis.T @ Basis:\n", np.diag(identity_approx))
    # print("Basis.T @ Basis:\n", clean(identity_approx))
    assert np.allclose(identity_approx, np.eye(basis.shape[1]), atol=1e-5)


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
        symmetrization="fast",
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

    # Check if basis.T @ basis is close to identity
    basis = symmetry_adapted_dof_space.basis
    identity_approx = basis.T @ basis
    # print("Diagonal of Basis.T @ Basis:\n", np.diag(identity_approx))
    # print("Basis.T @ Basis:\n", clean(identity_approx))
    assert np.allclose(identity_approx, np.eye(basis.shape[1]), atol=1e-5)


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
        # exclude_homogeneous_modes=False,
        # include_default_occ_modes=False,
        symmetrization="fast",
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

    # Check if basis.T @ basis is close to identity
    basis = symmetry_adapted_dof_space.basis
    identity_approx = basis.T @ basis
    # print("Diagonal of Basis.T @ Basis:\n", np.diag(identity_approx))
    # print("Basis.T @ Basis:\n", clean(identity_approx))
    assert np.allclose(identity_approx, np.eye(basis.shape[1]), atol=1e-5)

    # irreps = sym_report.irreps
    # for i, irrep in enumerate(irreps):
    #     print(f"Irrep {i}: dim={irrep.irrep_dim}, index={irrep.index}")
    #     # print characters
    #     print("- Characters:", irrep.characters)
    #     print()
    # assert False


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
        symmetrization="fast",
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

    # Check if basis.T @ basis is close to identity
    basis = symmetry_adapted_dof_space.basis
    identity_approx = basis.T @ basis
    # print("Diagonal of Basis.T @ Basis:\n", np.diag(identity_approx))
    # print("Basis.T @ Basis:\n", clean(identity_approx))
    assert np.allclose(identity_approx, np.eye(basis.shape[1]), atol=1e-5)


def test_dof_space_analysis_TlZn2Sb2_disp(TlZn2Sb2_disp_prim):
    xtal_prim = TlZn2Sb2_disp_prim

    prim = casmconfig.Prim(xtal_prim)
    spacegroup_type = spglib_io.get_spacegroup_type_from_symmetry(
        elements=prim.factor_group.elements,
        lattice=xtal_prim.lattice(),
    )
    # print("Factor group:")
    # print(prim.factor_group.brief_cart(lattice=xtal_prim.lattice()))
    # print(f"Space group: #{spacegroup_type.number}")
    assert spacegroup_type.number == 79

    T = (
        np.array(
            [
                [1, 0, 0],
                [0, 1, 0],
                [0, 0, 1],
            ],
            dtype=int,
        )
        * 2
    )
    supercell = casmconfig.Supercell(prim, T)
    configuration = casmconfig.Configuration(
        supercell=supercell,
    )

    dof_space = casmclex.DoFSpace(
        dof_key="disp",
        xtal_prim=xtal_prim,
        transformation_matrix_to_super=T,
    )

    results = casmconfig.dof_space_analysis(
        dof_space=dof_space,
        prim=prim,
        configuration=configuration,
        calc_wedges=False,
        exclude_homogeneous_modes=False,  # TODO: why is this necessary?
        verbosity="standard",
        commuter_method="deterministic",
    )
    assert isinstance(results, casmconfig.DoFSpaceAnalysisResults)
    basis = results.symmetry_adapted_dof_space.basis
    # symmetry_report = results.symmetry_report

    # print("# Irreps:", len(symmetry_report.irreps))
    # assert len(symmetry_report.irreps) == 30

    # print("Basis shape:", basis.shape)
    assert basis.shape == (supercell.n_sites * 3, supercell.n_sites * 3)
    # with np.printoptions(precision=6, suppress=True, linewidth=2000):
    #     print("Basis:\n", basis)

    # Check if basis.T @ basis is close to identity
    identity_approx = basis.T @ basis
    # print("Diagonal of Basis.T @ Basis:\n", np.diag(identity_approx))
    # print("Basis.T @ Basis:\n", clean(identity_approx))
    assert np.allclose(identity_approx, np.eye(basis.shape[1]), atol=1e-5)
