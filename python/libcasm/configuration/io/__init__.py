"""Additional methods for IO"""

from typing import Dict, List, Optional

import libcasm.configuration._configuration as _config

from ._symgroup import (
    symgroup_to_dict_with_group_classification,
)


def supercell_list_to_data(
    supercell_list: List[_config.Supercell],
) -> List[Dict]:
    """Represent a List[:class:`~libcasm.configuration.Supercell`] as a List[Dict]

    This method is equivalent to:

    .. code-block:: Python

        data_list = [x.to_dict() for x in supercell_list]


    Parameters
    ----------
    supercell_list: List[:class:`~libcasm.configuration.Supercell`]
        A list of :class:`~libcasm.configuration.Supercell`.

    Returns
    -------
    data_list: List[Dict]
        The representation of the supercell list as List[Dict].
    """
    return [x.to_dict() for x in supercell_list]


def supercell_list_from_data(
    data_list: List[Dict],
    prim: Optional[_config.Prim] = None,
    supercells: Optional[_config.SupercellSet] = None,
):
    """Construct a List[:class:`~libcasm.configuration.Supercell`] from a List[Dict]

    This method is equivalent to:

    .. code-block:: Python

        from libcasm.configuration import Supercell
        supercell_list = [Supercell.from_dict(data, supercells) for data in data_list]


    Parameters
    ----------
    data_list: List[Dict]
        A representation of List[:class:`~libcasm.configuration.Supercell`].
    prim: :class:`~libcasm.configuration.Prim`
        A :class:`~libcasm.configuration.Prim`, which is required if `supercells` is
        not provided.
    supercells: :class:`~libcasm.configuration.SupercellSet`
        A :class:`~libcasm.configuration.SupercellSet`, which may be provided to hold
        shared supercells in order to avoid duplicates.

    Returns
    -------
    supercell_list: List[:class:`~libcasm.configuration.Supercell`]
        The resulting list of supercells.
    """
    if prim is None and supercells is None:
        raise Exception(
            "Error in supercell_list_from_data: One of prim or supercells is required"
        )
    if supercells is None:
        supercells = _config.SupercellSet(prim)
    return [_config.Supercell.from_dict(data, supercells) for data in data_list]


def configuration_list_to_data(
    configuration_list: List[_config.Configuration],
) -> List[Dict]:
    """Represent a List[:class:`~libcasm.configuration.Configuration`] as a List[Dict]

    This method is equivalent to:

    .. code-block:: Python

        data_list = [x.to_dict() for x in configuration_list]


    Parameters
    ----------
    configuration_list: List[:class:`~libcasm.configuration.Configuration`]
        A list of :class:`~libcasm.configuration.Configuration`.

    Returns
    -------
    data_list: List[Dict]
        The representation of the configuration list as List[Dict].
    """
    return [x.to_dict() for x in configuration_list]


def configuration_list_from_data(
    data_list: List[Dict],
    prim: Optional[_config.Prim] = None,
    supercells: Optional[_config.SupercellSet] = None,
):
    """Construct a List[:class:`~libcasm.configuration.Configuration`] from a List[Dict]

    This method is equivalent to:

    .. code-block:: Python

        from libcasm.configuration import Configuration
        supercell_list = [
            Configuration.from_dict(data, supercells) for data in data_list
        ]


    Parameters
    ----------
    data_list: List[Dict]
        A representation of List[:class:`~libcasm.configuration.Configuration`].
    prim: :class:`~libcasm.configuration.Prim`
        A :class:`~libcasm.configuration.Prim`, which is required if `supercells` is
        not provided.
    supercells: :class:`~libcasm.configuration.SupercellSet`
        A :class:`~libcasm.configuration.SupercellSet`, which may be provided to hold
        shared supercells in order to avoid duplicates.

    Returns
    -------
    configuration_list: List[:class:`~libcasm.configuration.Configuration`]
        The resulting list of configurations.
    """
    if prim is None and supercells is None:
        raise Exception(
            "Error in configuration_list_from_data: "
            "One of prim or supercells is required"
        )
    if supercells is None:
        supercells = _config.SupercellSet(prim)
    return [_config.Configuration.from_dict(data, supercells) for data in data_list]
