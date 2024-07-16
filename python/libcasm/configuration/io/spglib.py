from __future__ import annotations

import copy
import dataclasses
import typing
import warnings

import numpy as np
import spglib

import libcasm.configuration as casmconfig
import libcasm.configuration.io.tools as io_tools
import libcasm.xtal as xtal


def asdict(
    obj: typing.Any,
):
    """If obj is a dataclass instance, convert to dict, converting numpy arrays to
    lists"""
    if obj is None:
        return None
    elif isinstance(obj, dict):
        data = copy.deepcopy(obj)
    else:
        data = dataclasses.asdict(obj)
    for key, value in data.items():
        if isinstance(value, np.ndarray):
            data[key] = value.tolist()
    return data


def as_lattice(obj: io_tools.CasmObjectWithLattice):
    """Get spglib lattice vectors from a lattice, structure, configuration, or prim

    Parameters
    ----------
    obj: :py:data:`CasmObjectWithLattice`
        A lattice, structure, configuration, or prim.

    Returns
    -------
    lattice: list[list[float]]
        The lattice vectors as rows (i.e. ``a == lattice[0]``, ``b == lattice[1]``,
        ``c == lattice[2]``).
    """
    return io_tools.get_lattice(obj).transpose().tolist()


def as_positions(obj: io_tools.CasmObjectWithPositions):
    """Get spglib positions from a structure, configuration, or prim

    Parameters
    ----------
    obj: :py:data:`CasmObjectWithPositions`
        An atomic structure, configuration, or prim.

    Returns
    -------
    positions: list[list[float]]
        The positions, as rows, in fractional coordinates.
    """
    return io_tools.get_coordinate_frac(obj).transpose().tolist()


def as_numbers(obj: io_tools.CasmObjectWithPositions):
    """Get spglib numbers from a structure, configuration, or prim

    Parameters
    ----------
    obj: :py:data:`CasmObjectWithPositions`
        An atomic structure, configuration, or prim.

    Returns
    -------
    numbers: list[int]
        For structures or configurations, a unique int for each atom type. For prim,
        the asymmetric unit indices.
    """
    return io_tools.get_unique_numbers(obj)


def as_magmoms(obj: io_tools.CasmObjectWithPositions):
    """Get spglib magmoms from a structure, configuration, or prim

    Parameters
    ----------
    obj: :py:data:`CasmObjectWithPositions`
        An atomic structure, configuration, or prim. Configurations are converted to
        structures using ``obj.to_structure()``.

    Returns
    -------
    magmoms: Union[list[float], list[list[float]], None]
        Magnetic moments associated with each position.

        For collinear magnetic spin, the return type is ``list[float]``.

        For non-collinear spins, the return type is ``list[list[float]]``.

        For prim with magnetic DoF, the values are set to value ``0.`` (collinear) or
        ``[0., 0., 0.]`` (non-collinear).

        Returns None if no magnetic spin properties or degrees of freedom (DoF) are
        present.

    """
    return io_tools.get_magmoms(obj)


def as_cell(obj: io_tools.CasmObjectWithPositions):
    """Convert a structure, configuration, or prim to a spglib cell tuple

    Parameters
    ----------
    obj: :py:data:`CasmObjectWithPositions`
        An atomic structure, configuration, or prim.

    Returns
    -------
    cell: tuple
        A spglib cell tuple, either ``(lattice, positions, numbers)`` or
        ``(lattice, positions, numbers, magmoms)``.

        Magmoms are included if magnetic spin properties or degrees of freedom (DoF)
        are present. For prim, magmoms are all set to value ``0.`` (collinear) or
        ``[0., 0., 0.]`` (non-collinear).
    """
    obj = io_tools.as_structure_if_config(obj)
    if isinstance(obj, (xtal.Prim, casmconfig.Prim)):
        if isinstance(obj, xtal.Prim):
            prim = casmconfig.Prim(obj)
        else:
            prim = obj
        if not prim.is_atomic:
            raise Exception("as_cell: not yet implemented for non-atomic prim")

    lattice = as_lattice(obj)
    positions = as_positions(obj)
    numbers = as_numbers(obj)
    magmoms = as_magmoms(obj)

    if magmoms:
        cell = (lattice, positions, numbers, magmoms)
    else:
        cell = (lattice, positions, numbers)
    return cell


def as_rotations(elements: list[xtal.SymOp], lattice: xtal.Lattice):
    """Get spglib rotations from a list of SymOp

    Parameters
    ----------
    elements: list[libcasm.xtal.SymOp]
        A list of symmetry group operations
    lattice: libcasm.xtal.Lattice
        The associated lattice

    Returns
    -------
    rotations: list[list[list[float]]], shape=(n_elements, 3, 3)
        The matrices from the symmetry group operations, in fractional coordinates.
    """
    rotations = []
    L = lattice.column_vector_matrix()
    Linv = np.linalg.inv(L)
    for op in elements:
        rotations.append(np.rint(Linv @ op.matrix() @ L).astype(int).tolist())
    return rotations


def as_translations(elements: list[xtal.SymOp], lattice: xtal.Lattice):
    """Get spglib translations from a list of SymOp

    Parameters
    ----------
    elements: list[libcasm.xtal.SymOp]
        A list of symmetry group operations
    lattice: libcasm.xtal.Lattice
        The associated lattice

    Returns
    -------
    translations: list[list[float]], shape=(n_elements, 3)
        The translation vectors from the symmetry group operations, in fractional
        coordinates.
    """
    Linv = np.linalg.inv(lattice.column_vector_matrix())
    translations = []
    for op in elements:
        translations.append((Linv @ op.translation()).tolist())
    return translations


def as_time_reversals(elements: list[xtal.SymOp]):
    """Get time reversal values from a list of SymOp

    Parameters
    ----------
    elements: list[libcasm.xtal.SymOp]
        A list of symmetry group operations

    Returns
    -------
    time_reversals: list[bool], shape=(n_elements,)
        The time_reversal values the symmetry group operations.
    """
    time_reversals = []
    for op in elements:
        time_reversals.append(op.time_reversal())
    return time_reversals


def get_symmetry_dataset(obj: io_tools.CasmObjectWithPositions):
    """Get the spglib symmetry dataset for a structure, configuration, or prim

    Notes
    -----

    - Magnetic spin properties or degrees of freedom (DoF) are not included in the
      symmetry group determination.

    Parameters
    ----------
    obj: :py:data:`CasmObjectWithPositions`
        An atomic structure, configuration, or prim.

    Returns
    -------
    dataset: Optional[spglib.SpglibDataset]
        A spglib symmetry dataset, as documented
        `here <https://spglib.readthedocs.io/en/latest/api/python-api/spglib/spglib.html#spglib.get_symmetry_dataset>`_.
        Returns None if getting the symmetry dataset fails. Note: for
        ``spglib<=2.4.0``, the return type is a dict.
    """
    cell = as_cell(obj)

    if len(cell) == 4:
        cell = (cell[0], cell[1], cell[2])
    dataset = spglib.get_symmetry_dataset(cell)
    return dataset


def get_magnetic_symmetry_dataset(
    obj: io_tools.CasmObjectWithPositions,
):
    """Get the spglib magnetic symmetry dataset for a structure, configuration, or prim

    Parameters
    ----------
    obj: :py:data:`CasmObjectWithPositions`
        An atomic structure, configuration, or prim. For prim, magmoms are all set to
        value ``0.`` (collinear) or ``[0., 0., 0.]`` (non-collinear).

    Returns
    -------
    dataset: Optional[spglib.SpglibMagneticDataset]
        A spglib magnetic symmetry dataset, as documented
        `here <https://spglib.readthedocs.io/en/latest/api/python-api/spglib/spglib.html#spglib.get_magnetic_symmetry_dataset>`_.
        Returns None if getting the magnetic symmetry dataset fails. Note: for
        ``spglib<=2.4.0``, the return type is a dict.
    """
    cell = as_cell(obj)

    if len(cell) == 3:
        cell = (cell[0], cell[1], cell[2], [0] * len(cell[2]))
    dataset = spglib.get_magnetic_symmetry_dataset(cell)
    return dataset


def get_spacegroup_type_from_symmetry(
    elements: list[xtal.SymOp],
    lattice: xtal.Lattice,
):
    """Get the spglib spacegroup type info from a list of SymOp

    Notes
    -----

    - To prevent the underlying spglib.get_spacegroup_type_from_symmetry method
      from segfaulting (as of v2.4.0), a ValueError is raised if `elements` contains
      operations with time reversal symmetry.

    Parameters
    ----------
    elements: list[libcasm.xtal.SymOp]
        A list of symmetry group operations.
    lattice: libcasm.xtal.Lattice
        The associated lattice

    Returns
    -------
    spacegroup_type: Optional[spglib.SpaceGroupType]
        A spglib.SpaceGroupType instance, as documented
        `here <https://spglib.readthedocs.io/en/latest/api/python-api/spglib/spglib.html#spglib.get_spacegroup_type_from_symmetry>`_,
        or None if getting the spacegroup type info fails. Note: for
        ``spglib<=2.4.0``, the return type is a dict.

    Raises
    ------
    ValueError
        If a symmetry operation has time reversal symmetry.
    """
    for op in elements:
        if op.time_reversal():
            raise ValueError(
                "Error in libcasm.configuration.io.spglib."
                "get_spacegroup_type_from_symmetry: "
                "not valid for operations with time reversal symmetry"
            )

    return spglib.get_spacegroup_type_from_symmetry(
        rotations=as_rotations(elements, lattice),
        translations=as_translations(elements, lattice),
        lattice=as_lattice(lattice),
    )


def get_magnetic_spacegroup_type_from_symmetry(
    elements: list[xtal.SymOp],
    lattice: xtal.Lattice,
):
    """Get the spglib magnetic spacegroup type info from a list of SymOp

    Parameters
    ----------
    elements: list[libcasm.xtal.SymOp]
        A list of symmetry group operations.
    lattice: libcasm.xtal.Lattice
        The associated lattice

    Returns
    -------
    magnetic_spacegroup_type: Optional[spglib.MagneticSpaceGroupType]
        A spglib.MagneticSpaceGroupType instance, as documented
        `here <https://spglib.readthedocs.io/en/latest/api/python-api/spglib/spglib.html#spglib.get_magnetic_spacegroup_type_from_symmetry>`_,
        or None if getting the magnetic spacegroup type info fails. Note: for
        ``spglib<=2.4.0``, the method is not available and this always returns None.

    """
    try:
        from spglib import get_magnetic_spacegroup_type_from_symmetry
    except ImportError:
        warnings.warn(
            "The function get_magnetic_spacegroup_type_from_symmetry is only available "
            "in spglib>=2.5.0."
        )
        return None

    return get_magnetic_spacegroup_type_from_symmetry(
        rotations=as_rotations(elements, lattice),
        translations=as_translations(elements, lattice),
        time_reversals=as_time_reversals(elements),
        lattice=as_lattice(lattice),
    )
