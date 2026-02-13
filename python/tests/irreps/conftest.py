import numpy as np
import pytest

import libcasm.clexulator as casmclex
import libcasm.configuration as casmconfig
import libcasm.group as casmgroup
import libcasm.irreps as casmirreps
import libcasm.xtal as xtal


@pytest.fixture(scope="session")
def FCC_binary_prim():
    # Lattice vectors
    lattice_column_vector_matrix = np.array(
        [
            [0.0, 1.0 / 2.0, 1.0 / 2.0],  # a
            [1.0 / 2.0, 0.0, 1.0 / 2.0],  # a
            [1.0 / 2.0, 1.0 / 2.0, 0.0],  # a
        ]
    ).transpose()
    lattice = xtal.Lattice(lattice_column_vector_matrix)

    # Basis sites positions, as columns of a matrix,
    # in fractional coordinates with respect to the lattice vectors
    coordinate_frac = np.array(
        [
            [0.0, 0.0, 0.0],
        ]
    ).transpose()

    # Occupation degrees of freedom (DoF)
    occupants = {}
    occ_dof = [["A", "B"]]

    return xtal.Prim(
        lattice=lattice,
        coordinate_frac=coordinate_frac,
        occ_dof=occ_dof,
        occupants=occupants,
    )


@pytest.fixture(scope="session")
def FCC_binary_irrep_decomposition(FCC_binary_prim):
    """IrrepDecomposition for FCC binary occ DoF in the conventional cell."""
    prim = casmconfig.Prim(FCC_binary_prim)
    T_dof_space = np.array(
        [
            [-1, 1, 1],
            [1, -1, 1],
            [1, 1, -1],
        ],
        dtype=int,
    )
    supercell = casmconfig.Supercell(prim, T_dof_space)
    configuration = casmconfig.Configuration(supercell=supercell)
    supercell_factor_group = casmconfig.make_invariant_subgroup(
        configuration=configuration,
    )
    symgroup = casmconfig.make_symgroup(supercell_factor_group)
    subset = casmgroup.Subset(group=symgroup)
    subset.all_subgroups()
    subgroup_orbits = subset.all_subgroup_orbits()

    dof_space = casmclex.DoFSpace(
        dof_key="occ",
        xtal_prim=FCC_binary_prim,
        transformation_matrix_to_super=T_dof_space,
    )
    matrix_rep = casmconfig.make_dof_space_rep(
        group=supercell_factor_group,
        dof_space=dof_space,
    )

    return casmirreps.IrrepDecomposition(
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


@pytest.fixture(scope="session")
def ABC2_disp_prim():
    # Lattice vectors
    lattice_column_vector_matrix = np.array(
        [
            [0.0, 0.0, 3.30],  # a
            [4.87, 0.0, 0.0],  # a
            [0.0, 4.87, 0.0],  # a
        ]
    ).transpose()
    lattice = xtal.Lattice(lattice_column_vector_matrix)

    # Basis sites positions, as columns of a matrix,
    # in fractional coordinates with respect to the lattice vectors
    coordinate_frac = np.array(
        [
            [0.0, 0.0, 0.0],
            [0.500000000000, 0.500000000000, 0.500000000000],
            [0.304640000000, 0.000000000000, 0.304640000000],
            [0.804640000000, 0.500000000000, 0.195360010000],
            [0.195360010000, 0.500000000000, 0.804640000000],
            [0.695360000000, 0.000000000000, 0.695360000000],
        ]
    ).transpose()

    # Occupation degrees of freedom (DoF)
    occupants = {}
    occ_dof = [["A"], ["B"], ["C"], ["C"], ["C"], ["C"]]

    # Local continuous degrees of freedom (DoF)
    disp_dof = xtal.DoFSetBasis("disp")  # Atomic displacement
    local_dof = [
        [disp_dof],
        [disp_dof],
        [disp_dof],
        [disp_dof],
        [disp_dof],
        [disp_dof],
    ]

    return xtal.Prim(
        lattice=lattice,
        coordinate_frac=coordinate_frac,
        occ_dof=occ_dof,
        local_dof=local_dof,
        occupants=occupants,
    )


@pytest.fixture(scope="session")
def ABC2_disp_irrep_decomposition(ABC2_disp_prim):
    """IrrepDecomposition for ABC2 prim with disp DoF."""
    xtal_prim = ABC2_disp_prim
    prim = casmconfig.Prim(xtal_prim)
    T_dof_space = np.eye(3, dtype=int) * 2  # 2x supercell in each direction
    supercell = casmconfig.Supercell(prim, T_dof_space)
    configuration = casmconfig.Configuration(supercell=supercell)
    supercell_factor_group = casmconfig.make_invariant_subgroup(
        configuration=configuration,
    )
    symgroup = casmconfig.make_symgroup(supercell_factor_group)
    subset = casmgroup.Subset(group=symgroup)
    subset.all_subgroups()
    subgroup_orbits = subset.all_subgroup_orbits()

    dof_space = casmclex.DoFSpace(
        dof_key="disp",
        xtal_prim=xtal_prim,
        transformation_matrix_to_super=T_dof_space,
    )
    matrix_rep = casmconfig.make_dof_space_rep(
        group=supercell_factor_group,
        dof_space=dof_space,
    )

    return casmirreps.IrrepDecomposition(
        matrix_rep=matrix_rep,
        subgroup_orbits=subgroup_orbits,
    )


@pytest.fixture(scope="session")
def TlZn2Sb2_disp_prim():
    L = np.array(
        [
            [-4.32450000000000, 4.32450000000000, 3.64349936250000],
            [4.32450000000000, -4.32450000000000, 3.64349936250000],
            [4.32450000000000, 4.32450000000000, -3.64349936250000],
        ]
    ).transpose()

    atom_type = ["Sb", "Sb", "Sb", "Sb", "Tl", "Tl", "Zn", "Zn", "Zn", "Zn"]

    atom_coordinate_frac = np.array(
        [
            [0.4434, 0.9434, 0.8260],  # asym 1
            [0.1174, 0.6174, 0.1740],  # asym 1
            [0.9434, 0.1174, 0.5000],  # asym 1
            [0.6174, 0.4434, 0.5000],  # asym 1
            [0.0000, 0.0000, 0.0000],  # asym 2
            [0.5104, 0.5104, 0.0000],  # asym 3
            [0.2794, 0.0474, 0.4980],  # asym 4
            [0.5494, 0.7814, 0.5020],  # asym 4
            [0.0474, 0.5494, 0.7680],  # asym 4
            [0.7814, 0.2794, 0.2320],  # asym 4
        ],
    ).transpose()
    return xtal.Prim(
        lattice=xtal.Lattice(L),
        coordinate_frac=atom_coordinate_frac,
        occ_dof=[[x] for x in atom_type],
        local_dof=[[xtal.DoFSetBasis("disp")] for x in atom_type],
    )


@pytest.fixture(scope="session")
def TlZn2Sb2_disp_irrep_decomposition(TlZn2Sb2_disp_prim):
    """IrrepDecomposition for ABC2 prim with disp DoF."""
    xtal_prim = TlZn2Sb2_disp_prim
    prim = casmconfig.Prim(xtal_prim)
    T_dof_space = np.eye(3, dtype=int) * 1  # 2x supercell in each direction
    supercell = casmconfig.Supercell(prim, T_dof_space)
    configuration = casmconfig.Configuration(supercell=supercell)
    supercell_factor_group = casmconfig.make_invariant_subgroup(
        configuration=configuration,
    )
    symgroup = casmconfig.make_symgroup(supercell_factor_group)
    subset = casmgroup.Subset(group=symgroup)
    subset.all_subgroups()
    subgroup_orbits = subset.all_subgroup_orbits()

    dof_space = casmclex.DoFSpace(
        dof_key="disp",
        xtal_prim=xtal_prim,
        transformation_matrix_to_super=T_dof_space,
    )
    matrix_rep = casmconfig.make_dof_space_rep(
        group=supercell_factor_group,
        dof_space=dof_space,
    )

    return casmirreps.IrrepDecomposition(
        matrix_rep=matrix_rep,
        subgroup_orbits=subgroup_orbits,
    )
