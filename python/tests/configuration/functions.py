import typing

import numpy as np
import spglib

import libcasm.xtal as xtal


def check_symmetry_dataset(
    symmetry_dataset: typing.Any,
    number: int,
    n_rotations: int,
    hall_number: typing.Optional[int] = None,
):
    if isinstance(symmetry_dataset, dict):
        # for spglib<2.5.0
        assert symmetry_dataset["number"] == number
        assert len(symmetry_dataset["rotations"]) == n_rotations
        if hall_number is not None:
            assert symmetry_dataset["hall_number"] == hall_number

    else:
        # for spglib>=2.5.0
        assert isinstance(symmetry_dataset, spglib.SpglibDataset)
        assert symmetry_dataset.number == number
        assert len(symmetry_dataset.rotations) == n_rotations

        if hall_number is not None:
            assert symmetry_dataset.hall_number == hall_number


def check_magnetic_symmetry_dataset(
    mag_symmetry_dataset: typing.Any,
    uni_number: int,
    n_operations: int,
    hall_number: typing.Optional[int] = None,
):
    if isinstance(mag_symmetry_dataset, dict):
        # for spglib<2.5.0
        assert mag_symmetry_dataset["uni_number"] == uni_number
        assert mag_symmetry_dataset["n_operations"] == n_operations
        if hall_number is not None:
            assert mag_symmetry_dataset["hall_number"] == hall_number

    else:
        # for spglib>=2.5.0
        assert isinstance(mag_symmetry_dataset, spglib.SpglibMagneticDataset)
        assert mag_symmetry_dataset.uni_number == uni_number
        assert mag_symmetry_dataset.n_operations == n_operations

        if hall_number is not None:
            assert mag_symmetry_dataset.hall_number == hall_number


def check_spacegroup_type(
    spacegroup_type: typing.Any,
    number: int,
    hall_number: typing.Optional[int] = None,
):
    if isinstance(spacegroup_type, dict):
        # for spglib<2.5.0
        assert spacegroup_type["number"] == number
        if hall_number is not None:
            assert spacegroup_type["hall_number"] == hall_number

    else:
        # for spglib>=2.5.0
        assert isinstance(spacegroup_type, spglib.SpaceGroupType)

        assert spacegroup_type.number == number
        if hall_number is not None:
            assert spacegroup_type.hall_number == hall_number


def make_discrete_magnetic_atom(
    name: str,
    value: typing.Any,
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
