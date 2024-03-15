from typing import Any

import numpy as np

import libcasm.configuration as config
import libcasm.xtal as xtal


def make_discrete_magnetic_atom(
    name: str,
    value: Any,
    flavor: str = "C",
) -> xtal.Occupant:
    """Construct a discrete magnetic atomic occupant

    Parameters
    ----------
    name: str
        A "chemical name", which must be identical for atoms to be found symmetrically
        equivalent. The names are case sensitive, and “Va” is reserved for vacancies.
    value: Any
        The discrete value of the magnetic spin to associate with the constructed
        Occupant. If the type is `int` or `float` the value is converted to a
        size 1 array of float, other types are converted using `numpy.asarray`.
    flavor: str
        The magnetic spin "flavor", which must be one of varieties supported by CASM:
        `C`, `NC`, `SO`.

    Returns
    -------
    discrete_magnetic_atom: xtal.Occupant
        An :class:`~libcasm.xtal.Occupant` consisting of a single atom with the
        specified magnetic spin flavor and value.
    """
    if isinstance(value, (int, float)):
        value = np.array([value], dtype=np.float64)
    else:
        value = np.asarray(value, dtype=np.float64)

    return xtal.Occupant(
        name=name,
        atoms=[
            xtal.AtomComponent(
                name=name,
                coordinate=[0.0, 0.0, 0.0],
                properties={flavor + "magspin": value},
            )
        ],
    )


def FCC_binary_discrete_Cmagspin_prim():
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
    occupants = {
        "A.up": make_discrete_magnetic_atom(name="A", value=1, flavor="C"),
        "A.down": make_discrete_magnetic_atom(name="A", value=-1, flavor="C"),
        "B.up": make_discrete_magnetic_atom(name="B", value=1, flavor="C"),
        "B.down": make_discrete_magnetic_atom(name="B", value=-1, flavor="C"),
    }
    occ_dof = [["A.up", "A.down", "B.up", "B.down"]]

    return xtal.Prim(
        lattice=lattice,
        coordinate_frac=coordinate_frac,
        occ_dof=occ_dof,
        occupants=occupants,
    )


def test_FCC_binary_discrete_Cmagspin_prim_1():
    xtal_prim = FCC_binary_discrete_Cmagspin_prim()
    factor_group = xtal.make_factor_group(xtal_prim)
    assert len(factor_group) == 96
    prim = config.Prim(xtal_prim)
    assert prim.xtal_prim.coordinate_frac().shape == (3, 1)
    assert len(prim.factor_group.elements) == 96

    assert prim.is_atomic
    assert prim.continuous_magspin_key is None
    assert prim.continuous_magspin_flavor is None
    assert prim.discrete_atomic_magspin_key == "Cmagspin"
    assert prim.discrete_atomic_magspin_flavor == "C"

    lattice = prim.xtal_prim.lattice()
    for op in prim.factor_group.elements:
        syminfo = xtal.SymInfo(op, lattice)
        print(xtal.pretty_json(syminfo.to_dict()))
