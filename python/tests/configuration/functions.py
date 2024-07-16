import typing

import spglib


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
