from functools import total_ordering
from typing import Optional

import libcasm.clusterography
import libcasm.configuration
import libcasm.configuration as casmconfig
import libcasm.occ_events
import libcasm.xtal

from ._OccEventPrimSymInfo import OccEventPrimSymInfo
from ._OccEventSupercellSymInfo import OccEventSupercellSymInfo


@total_ordering
class LocalConfiguration:
    def __init__(
        self,
        pos: tuple[int, int],
        configuration: libcasm.configuration.Configuration,
    ):
        self.pos = pos
        """tuple[int, int]: The position of the phenomenal cluster or event in the 
        supercell, as `(unitcell_index, equivalent_index)`."""

        self.configuration = configuration
        """libcasm.configuration.Configuration: The configuration."""

    def copy(self):
        return LocalConfiguration(
            pos=self.pos,
            configuration=self.configuration.copy(),
        )

    def __copy__(self):
        return self.copy()

    def __deepcopy__(self, memodict={}):
        return self.copy()

    def __eq__(self, other):
        return self.pos == other.pos and self.configuration == other.configuration

    def __lt__(self, other):
        # Compare pos, then configuration
        if self.pos < other.pos:
            return True
        if self.pos > other.pos:
            return False
        return self.configuration < other.configuration

    def __repr__(self):
        return libcasm.xtal.pretty_json(self.to_dict())

    def to_dict(
        self,
        write_prim_basis: bool = False,
    ):
        """Represent the LocalConfiguration as a Python dict

        Parameters
        ----------
        write_prim_basis: bool = False
            If True, write DoF values using the prim basis. Default (False) is to
            write DoF values in the standard basis.

        Returns
        -------
        data: dict
            A Python dict representation of the LocalConfiguration.
        """
        return {
            "configuration": self.configuration.to_dict(
                write_prim_basis=write_prim_basis
            ),
            "pos": [self.pos[0], self.pos[1]],
        }

    @staticmethod
    def from_dict(
        data: dict,
        supercells: libcasm.configuration.SupercellSet,
    ):
        """Construct the LocalConfiguration from a Python dict

        Parameters
        ----------
        data : dict
            A :class:`~libcasm.enumerate.LocalConfiguration` as a dict.
        supercells : libcasm.configuration.SupercellSet
            A :class:`~libcasm.configuration.SupercellSet`, which holds shared
            supercells in order to avoid duplicates.

        Returns
        -------
        local_configuration : libcasm.enumerate.LocalConfiguration
            The :class:`~libcasm.configuration.Configuration` constructed from the dict.
        """
        return LocalConfiguration(
            configuration=libcasm.configuration.Configuration.from_dict(
                data=data["configuration"],
                supercells=supercells,
            ),
            pos=(data["pos"][0], data["pos"][1]),
        )


class LocalConfigurationList:
    def __init__(
        self,
        prim_sym_info: OccEventPrimSymInfo,
        local_configurations: list[LocalConfiguration] = [],
    ):
        self.prim_sym_info = prim_sym_info
        """OccEventPrimSymInfo: Information about the OccEvent with respect to the prim,
        which defines the meaning of the LocalConfiguration.pos attribute."""

        self.scel_sym_info = dict()
        """dict[libcasm.configuration.Supercell, OccEventSupercellSymInfo]: Information 
        about the OccEvent with respect to the supercell."""

        self._local_configurations = local_configurations
        """list[LocalConfiguration]: The list of local configurations."""

    def _check_value(self, value):
        e = ValueError(
            "Error: value must be a LocalConfiguration or "
            "tuple[OccEvent, Configuration]"
        )
        if isinstance(value, tuple):
            if (
                len(value) != 2
                or not isinstance(value[0], libcasm.occ_events.OccEvent)
                or not isinstance(value[1], libcasm.configuration.Configuration)
            ):
                raise e
            event = value[0]
            config = value[1]
            supercell = config.supercell
            if supercell not in self.scel_sym_info:
                self.scel_sym_info[supercell] = OccEventSupercellSymInfo(
                    prim_sym_info=self.prim_sym_info,
                    supercell=supercell,
                )
            pos = self.scel_sym_info[supercell].coordinate(event)
            return LocalConfiguration(pos=pos, configuration=config)
        elif not isinstance(value, LocalConfiguration):
            raise e
        return value

    def as_local_configuration(
        self,
        event: libcasm.occ_events.OccEvent,
        configuration: libcasm.configuration.Configuration,
    ) -> LocalConfiguration:
        supercell = configuration.supercell
        if supercell not in self.scel_sym_info:
            self.scel_sym_info[supercell] = OccEventSupercellSymInfo(
                prim_sym_info=self.prim_sym_info,
                supercell=supercell,
            )
        pos = self.scel_sym_info[supercell].coordinate(event)
        return LocalConfiguration(pos=pos, configuration=configuration)

    def as_tuple(
        self,
        local_configuration: LocalConfiguration,
    ) -> tuple[libcasm.occ_events.OccEvent, libcasm.configuration.Configuration]:
        """Convert a LocalConfiguration to a tuple of OccEvent and Configuration.

        Parameters
        ----------
        local_configuration : LocalConfiguration
            The LocalConfiguration to convert.

        Returns
        -------
        value: tuple[libcasm.occ_events.OccEvent, libcasm.configuration.Configuration]
            The equivalent tuple of OccEvent and Configuration.

        """
        pos = local_configuration.pos
        config = local_configuration.configuration
        unitcell = config.supercell.unitcell_index_converter.unitcell(pos[0])
        event = self.prim_sym_info.events[pos[1]] + unitcell
        return (event, config)

    def __contains__(self, value):
        value = self._check_value(value)
        return value in self._local_configurations

    def __len__(self):
        return len(self._local_configurations)

    def __getitem__(self, i):
        return self._local_configurations[i]

    def __setitem__(self, i, value):
        self._local_configurations[i] = value

    def __delitem__(self, i):
        del self._local_configurations[i]

    def __iter__(self):
        return iter(self._local_configurations)

    def append(self, value):
        self._local_configurations.append(value)

    def sort(self):
        self._local_configurations.sort()

    def index(self, value):
        return self._local_configurations.index(value)

    def clear(self):
        self._local_configurations.clear()

    def copy(self):
        return LocalConfigurationList(
            prim_sym_info=self.prim_sym_info.copy(),
            local_configurations=[lc.copy() for lc in self._local_configurations],
        )

    def to_dict(
        self,
        system: libcasm.occ_events.OccSystem,
    ):
        (
            prototype_event_data,
            equivalents_info_data,
        ) = self.prim_sym_info.to_data(system=system)
        return {
            "event": prototype_event_data,
            "equivalents_info": equivalents_info_data,
            "local_configurations": [lc.to_dict() for lc in self.local_configurations],
        }

    @staticmethod
    def from_dict(
        data: dict,
        prim: libcasm.configuration.Prim,
        system: libcasm.occ_events.OccSystem,
        supercells: libcasm.configuration.SupercellSet,
    ):
        return LocalConfigurationList(
            prim_sym_info=OccEventPrimSymInfo.from_data(
                prototype_event_data=data["prim_sym_info"],
                equivalents_info_data=data["equivalents_info"],
                prim=prim,
                system=system,
            ),
            local_configurations=[
                LocalConfiguration.from_dict(
                    data=lc,
                    supercells=supercells,
                )
                for lc in data["local_configurations"]
            ],
        )


class ConfigEnumLocalOccupations:
    def __init__(
        self,
        prim_sym_info: OccEventPrimSymInfo,
        local_configurations: LocalConfigurationList,
    ):
        self.prim_sym_info = prim_sym_info
        """OccEventPrimSymInfo: Information about the OccEvent with respect to the prim,
        which defines the meaning of the LocalConfiguration.pos attribute."""

    def by_cluster(
        self,
        background: casmconfig.Configuration,
        cluster_specs: dict,
        supercells: Optional[dict] = None,
        skip_non_primitive: bool = False,
        skip_equivalents: bool = True,
    ):
        pass
