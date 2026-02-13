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
