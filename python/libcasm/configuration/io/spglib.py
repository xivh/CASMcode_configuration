from __future__ import annotations

import typing

import numpy as np
import spglib

import libcasm.configuration as casmconfig
import libcasm.xtal
import libcasm.xtal as xtal

CasmObjectWithLattice: typing.TypeAlias = typing.Union[
    libcasm.xtal.Lattice,
    libcasm.xtal.Structure,
    libcasm.xtal.Prim,
    libcasm.configuration.Prim,
    libcasm.configuration.Configuration,
    libcasm.configuration.ConfigurationWithProperties,
]
"""TypeAlias for a CASM object with a lattice

May be:

- a :class:`libcasm.xtal.Lattice`,
- a :class:`libcasm.xtal.Structure`,
- a :class:`libcasm.xtal.Prim`,
- a :class:`libcasm.configuration.Prim`,
- a :class:`libcasm.configuration.Configuration`, or
- a :class:`libcasm.configuration.ConfigurationWithProperties`.

Note that :class:`~libcasm.configuration.Configuration` and 
:class:`~libcasm.configuration.ConfigurationWithProperties` are converted to 
:class:`~libcasm.xtal.Structure` using the equivalent of:

.. code-block:: Python

    if isinstance(obj, libcasm.configuration.Configuration):
        obj = obj.to_structure()
    if isinstance(obj, libcasm.configuration.ConfigurationWithProperties):
        obj = obj.to_structure()
"""

CasmObjectWithPositions: typing.TypeAlias = typing.Union[
    libcasm.xtal.Structure,
    libcasm.xtal.Prim,
    libcasm.configuration.Prim,
    libcasm.configuration.Configuration,
    libcasm.configuration.ConfigurationWithProperties,
]
"""TypeAlias for a CASM object with positions

May be:

- a :class:`libcasm.xtal.Structure`,
- a :class:`libcasm.xtal.Prim`,
- a :class:`libcasm.configuration.Prim`,
- a :class:`libcasm.configuration.Configuration`, or
- a :class:`libcasm.configuration.ConfigurationWithProperties`.

Note that :class:`~libcasm.configuration.Configuration` and 
:class:`~libcasm.configuration.ConfigurationWithProperties` are converted to 
:class:`~libcasm.xtal.Structure` using the equivalent of:

.. code-block:: Python

    if isinstance(obj, libcasm.configuration.Configuration):
        obj = obj.to_structure()
    if isinstance(obj, libcasm.configuration.ConfigurationWithProperties):
        obj = obj.to_structure()
"""


def _to_structure_if_config(
    obj: typing.Union[casmconfig.Configuration, casmconfig.ConfigurationWithProperties]
):
    if isinstance(obj, casmconfig.Configuration):
        return obj.to_structure()
    if isinstance(obj, casmconfig.ConfigurationWithProperties):
        return obj.to_structure()
    return obj


def as_lattice(obj: CasmObjectWithLattice):
    """Get spglib lattice vectors from a lattice, structure, configuration, or prim

    Parameters
    ----------
    obj: :py:data:`CasmObjectWithLattice`
        A lattice, structure, configuration, or prim.

    Returns
    -------
    lattice: list[list[float]]
        The lattice vectors as rows.
    """
    obj = _to_structure_if_config(obj)
    if isinstance(obj, xtal.Lattice):
        lattice = obj
    elif isinstance(obj, xtal.Prim):
        lattice = obj.lattice()
    elif isinstance(obj, casmconfig.Prim):
        lattice = obj.xtal_prim.lattice()
    elif isinstance(obj, xtal.Structure):
        lattice = obj.lattice()
    else:
        raise Exception(
            "Error in get_spglib_Lattice: "
            "not a lattice, structure, configuration, or prim"
        )
    return lattice.column_vector_matrix().transpose().tolist()


def as_positions(obj: CasmObjectWithPositions):
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
    obj = _to_structure_if_config(obj)
    if isinstance(obj, xtal.Prim):
        coordinate_frac = obj.coordinate_frac()
    elif isinstance(obj, casmconfig.Prim):
        coordinate_frac = obj.xtal_prim.coordinate_frac()
    elif isinstance(obj, xtal.Structure):
        coordinate_frac = obj.atom_coordinate_frac()
    else:
        raise Exception(
            "Error in get_spglib_Lattice: "
            "not a libcasm.xtal.Structure, "
            "libcasm.xtal.Prim, or libcasm.configuration.Prim"
        )

    return coordinate_frac.transpose().tolist()


def as_numbers(obj: CasmObjectWithPositions):
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
    obj = _to_structure_if_config(obj)
    if isinstance(obj, (xtal.Prim, casmconfig.Prim)):
        if isinstance(obj, xtal.Prim):
            xtal_prim = obj
        else:
            xtal_prim = obj.xtal_prim

        numbers = [None] * xtal_prim.coordinate_frac().shape[1]
        asym_indices = xtal.asymmetric_unit_indices(xtal_prim)
        for a, asym_unit in enumerate(asym_indices):
            for site_index in asym_unit:
                numbers[site_index] = a

    elif isinstance(obj, xtal.Structure):
        atom_type = obj.atom_type()
        unique_atom_types = list(set(atom_type))
        numbers = [unique_atom_types.index(atom) for atom in atom_type]

    else:
        raise Exception(
            "Error in get_spglib_Numbers: "
            "not a libcasm.xtal.Structure, "
            "libcasm.xtal.Prim, or libcasm.configuration.Prim"
        )

    return numbers


def as_magmoms(obj: CasmObjectWithPositions):
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
    obj = _to_structure_if_config(obj)
    if isinstance(obj, (xtal.Prim, casmconfig.Prim)):
        if isinstance(obj, xtal.Prim):
            prim = casmconfig.Prim(obj)
        else:
            prim = obj

        n_sites = prim.xtal_prim.coordinate_frac().shape[1]
        if prim.continuous_magspin_key:
            if prim.continuous_magspin_key[:1] == "C":
                magmoms = [0] * n_sites
            elif prim.continuous_magspin_key[:2] in ["NC", "SO"]:
                magmoms = [[0] * 3] * n_sites
            else:
                raise Exception(
                    "Error in as_magmoms: "
                    f"Unexpected continuous_magspin_key={prim.continuous_magspin_key}"
                )

        elif prim.discrete_atomic_magspin_key:
            if prim.discrete_atomic_magspin_key[:1] == "C":
                magmoms = [0] * n_sites
            elif prim.discrete_atomic_magspin_key[:2] in ["NC", "SO"]:
                magmoms = [[0] * 3] * n_sites
            else:
                raise Exception(
                    "Error in as_magmoms: Unexpected "
                    f"discrete_atomic_magspin_key={prim.continuous_magspin_key}"
                )

        else:
            magmoms = None

    elif isinstance(obj, xtal.Structure):
        atom_properties = obj.atom_properties()

        magmoms = None
        for key in atom_properties:
            if "magspin" in key:
                value = atom_properties[key]
                if value.shape[0] == 1:
                    magmoms = value[0, :].tolist()
                elif value.shape[0] == 3:
                    magmoms = value.transpose().tolist()
                else:
                    raise Exception(
                        "Error in as_magmoms: "
                        f"Invalid row dimension for atom_properties[{key}], "
                        f"must be 1 or 3, found {value.shape[0]}."
                    )
                break

    else:
        raise Exception(
            "Error in get_spglib_Magmoms: "
            "not a libcasm.xtal.Structure, "
            "libcasm.xtal.Prim, or libcasm.configuration.Prim"
        )

    return magmoms


def as_cell(obj: CasmObjectWithPositions):
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
    obj = _to_structure_if_config(obj)
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


def get_symmetry_dataset(obj: CasmObjectWithPositions):
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
    dataset: Optional[dict]
        A spglib symmetry dataset, as documented
        `here <https://spglib.readthedocs.io/en/latest/api/python-api/spglib/spglib.html#spglib.get_symmetry_dataset>`_.
        Returns None if getting the symmetry dataset fails.
    """
    cell = as_cell(obj)

    if len(cell) == 4:
        cell = (cell[0], cell[1], cell[2])
    dataset = spglib.get_symmetry_dataset(cell)
    return dataset


def get_magnetic_symmetry_dataset(
    obj: CasmObjectWithPositions,
):
    """Get the spglib magnetic symmetry dataset for a structure, configuration, or prim

    Parameters
    ----------
    obj: :py:data:`CasmObjectWithPositions`
        An atomic structure, configuration, or prim. For prim, magmoms are all set to
        value ``0.`` (collinear) or ``[0., 0., 0.]`` (non-collinear).

    Returns
    -------
    dataset: Optional[dict]
        A spglib magnetic symmetry dataset, as documented
        `here <https://spglib.readthedocs.io/en/latest/api/python-api/spglib/spglib.html#spglib.get_magnetic_symmetry_dataset>`_.
        Returns None if getting the magnetic symmetry dataset fails.
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
    spacegroup_type: Optional[dict]
        A spglib `spacegroup_type` dict, as documented
        `here <https://spglib.readthedocs.io/en/latest/api/python-api/spglib/spglib.html#spglib.get_spacegroup_type_from_symmetry>`_,
        or None if getting the spacegroup type info fails.

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
