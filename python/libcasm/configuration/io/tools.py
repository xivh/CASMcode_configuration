"""Methods for getting properties from CASM objects"""

import typing

import libcasm.configuration
import libcasm.configuration as casmconfig
import libcasm.xtal
import libcasm.xtal as xtal

try:
    # python 3.10+
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
    
    .. note::
        typing.TypeAlias was introduced in Python3.10; for Python3.9 
        CasmObjectWithLattice is typing.Any
    
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
    
    .. note::
        typing.TypeAlias was introduced in Python3.10; for Python3.9 
        CasmObjectWithPositions is typing.Any
    
    """
except AttributeError:
    # python 3.9
    CasmObjectWithLattice = typing.Any
    CasmObjectWithPositions = typing.Any


def as_structure_if_config(obj: typing.Any):
    """If a Configuration or ConfigurationWithProperties, convert to a Structure

    Parameters
    ----------
    obj: Any
        Any object

    Returns
    -------
    structure_if_config: Any
        If `obj` is a Configuration or ConfigurationWithProperties, call
        `obj.to_structure()` and return the resulting Structure. Otherwise, returns
        `obj`.
    """
    if isinstance(obj, casmconfig.Configuration):
        return obj.to_structure()
    if isinstance(obj, casmconfig.ConfigurationWithProperties):
        return obj.to_structure()
    return obj


def get_lattice(obj: CasmObjectWithLattice):
    """Get lattice vectors from a lattice, structure, configuration, or prim

    Parameters
    ----------
    obj: :py:data:`CasmObjectWithLattice`
        A lattice, structure, configuration, or prim.

    Returns
    -------
    lattice: np.ndarray[np.float64[3,3]]
        The lattice vectors as columns (i.e. ``a == lattice[:,0]``,
        ``b == lattice[:,1]``, ``c == lattice[:,2]``).
    """
    obj = as_structure_if_config(obj)
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
            "Error in get_lattice: " "not a lattice, structure, configuration, or prim"
        )
    return lattice.column_vector_matrix()


def get_coordinate_frac(obj: CasmObjectWithPositions):
    """Get fractional coordinates from a structure, configuration, or prim

    Parameters
    ----------
    obj: :py:data:`CasmObjectWithPositions`
        An atomic structure, configuration, or prim.

    Returns
    -------
    positions: np.ndarray[np.float64[3, n_sites]]
        The positions of atoms or sites, as columns, in fractional coordinates.
    """
    obj = as_structure_if_config(obj)
    if isinstance(obj, xtal.Prim):
        coordinate_frac = obj.coordinate_frac()
    elif isinstance(obj, casmconfig.Prim):
        coordinate_frac = obj.xtal_prim.coordinate_frac()
    elif isinstance(obj, xtal.Structure):
        coordinate_frac = obj.atom_coordinate_frac()
    else:
        raise Exception(
            "Error in get_positions: "
            "not a libcasm.xtal.Structure, "
            "libcasm.xtal.Prim, or libcasm.configuration.Prim"
        )

    return coordinate_frac


def get_coordinate_cart(obj: CasmObjectWithPositions):
    """Get Cartesian coordinates from a structure, configuration, or prim

    Parameters
    ----------
    obj: :py:data:`CasmObjectWithPositions`
        An atomic structure, configuration, or prim.

    Returns
    -------
    positions: np.ndarray[np.float64[3, n_sites]]
        The positions of atoms or sites, as columns, in Cartesian coordinates.
    """
    obj = as_structure_if_config(obj)
    if isinstance(obj, xtal.Prim):
        coordinate_cart = obj.coordinate_cart()
    elif isinstance(obj, casmconfig.Prim):
        coordinate_cart = obj.xtal_prim.coordinate_cart()
    elif isinstance(obj, xtal.Structure):
        coordinate_cart = obj.atom_coordinate_cart()
    else:
        raise Exception(
            "Error in get_coordinate_cart: "
            "not a libcasm.xtal.Structure, "
            "libcasm.xtal.Prim, or libcasm.configuration.Prim"
        )

    return coordinate_cart


def get_unique_numbers(
    obj: CasmObjectWithPositions,
):
    """Get unique occupant numbers from a structure, configuration, or prim

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
    obj = as_structure_if_config(obj)
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
            "Error in get_unique_numbers: "
            "not a libcasm.xtal.Structure, "
            "libcasm.xtal.Prim, or libcasm.configuration.Prim"
        )

    return numbers


def get_magmoms(obj: CasmObjectWithPositions):
    """Get magmoms from a structure, configuration, or prim

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
    obj = as_structure_if_config(obj)
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
                    "Error in get_magmoms: "
                    f"Unexpected continuous_magspin_key={prim.continuous_magspin_key}"
                )

        elif prim.discrete_atomic_magspin_key:
            if prim.discrete_atomic_magspin_key[:1] == "C":
                magmoms = [0] * n_sites
            elif prim.discrete_atomic_magspin_key[:2] in ["NC", "SO"]:
                magmoms = [[0] * 3] * n_sites
            else:
                raise Exception(
                    "Error in get_magmoms: Unexpected "
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
                        "Error in get_magmoms: "
                        f"Invalid row dimension for atom_properties[{key}], "
                        f"must be 1 or 3, found {value.shape[0]}."
                    )
                break

    else:
        raise Exception(
            "Error in get_magmoms: "
            "not a libcasm.xtal.Structure, "
            "libcasm.xtal.Prim, or libcasm.configuration.Prim"
        )

    return magmoms
