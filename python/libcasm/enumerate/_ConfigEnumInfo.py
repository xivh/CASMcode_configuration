import sys
from typing import Any

import libcasm.configuration as casmconfig
import libcasm.xtal as xtal


class ConfigEnumInfo:
    """Helper for printing information during configuration enumeration

    .. rubric:: Example usage

    The following demonstrates constructing an occupation configuration enumerator,
    and using it to enumerate distinct occupations in supercells up to volume 4
    in ternary FCC prim. The ConfigEnumInfo object allows printing the current
    background configuration and sites on which the enumeration is occurring, along
    with information about the number of configurations generated.

    .. code-block:: Python

        import libcasm.configuration as casmconfig
        import libcasm.enumeration as casmenum


        def filter(configuration: casmconfig.Configuration):
            # A custom filter function
            # return True to keep configuration, False to skip
            return True


        xtal_prim = xtal_prims.FCC(
            r=0.5,
            occ_dof=["A", "B", "C"],
        )
        prim = casmconfig.Prim(xtal_prim)
        supercell_set = casmconfig.SupercellSet(prim=prim)
        configuration_set = casmconfig.ConfigurationSet()
        config_enum = casmenum.ConfigEnumAllOccupations(
            prim=prim,
            supercell_set=supercell_set,
        )
        info = ConfigEnumInfo(config_enum, configuration_set)
        for i, configuration in enumerate(
            config_enum.by_supercell(
                supercells={
                    "max": 12,
                }
            )
        ):
            info.check()
            info.n_config_total += 1
            if not filter(configuration):
                info.n_config_excluded += 1
            else:
                configuration_set.add(configuration)
        info.finish()

    Example output:

    ::

        ~~~~~~
        Background:
         {
          "dof": {
            "occ": [0]
          },
          "supercell_name": "SCEL1_1_1_1_0_0_0",
          "transformation_matrix_to_supercell": [
            [1, 0, 0],
            [0, 1, 0],
            [0, 0, 1]
          ]
        }
        Sites: {0}
        3 configurations (3 new, 0 excluded by filter)

        ~~~~~~
        Background:
         {
          "dof": {
            "occ": [0, 0]
          },
          "supercell_name": "SCEL2_2_1_1_0_1_1",
          "transformation_matrix_to_supercell": [
            [1, 0, 1],
            [0, 1, 1],
            [-1, -1, 0]
          ]
        }
        Sites: {0, 1}
        3 configurations (3 new, 0 excluded by filter)

        ...


    """

    def __init__(
        self, config_enum: Any, configuration_set: casmconfig.ConfigurationSet
    ):
        """
        .. rubric:: Constructor

        Parameters
        ----------
        config_enum: Any
            The current configuration enumerator. Expected to have the following
            attributes:

            - background: :class:`~libcasm.configuration.Configuration`, The current
              background configuration in which enumeration is being performed
            - sites: ``set[int]``, The linear site indices of the current set of sites
              on which enumeration is being performed
            - enum_index: ``int``, The current index over combinations of background and
              sites

        configuration_set: libcasm.configuration.ConfigurationSet
            A ConfigurationSet in which configurations are being saved.
        """
        self.config_enum = config_enum
        """Any: The current configuration enumerator. 
        
        Expected to have the following attributes:

        - background: :class:`~libcasm.configuration.Configuration`, The current
          background configuration in which enumeration is being performed
        - sites: ``set[int]``, The linear site indices of the current set of sites on 
          which enumeration is being performed
        - enum_index: ``int``, The current index over combinations of background and 
          sites
        """

        self.configuration_set = configuration_set
        """libcasm.configuration.ConfigurationSet: A ConfigurationSet in which 
        configurations are being saved."""

        self.curr_enum_index = None
        """Optional[int]: During enumeration, `curr_enum_index` is set to the index of 
        the current combination of background and sites"""

        self.n_config_before = len(self.configuration_set)
        """int: The size of `configuration_set` before the current combination of 
        background and sites"""

        self.n_config_total = 0
        """int: The number of configurations generated with the current combination of
        background and sites"""

        self.n_config_excluded = 0
        """int: The number of configurations generated with the current combination of
        background and sites, and excluded by a user's filter"""

    def _begin_background(self):
        self.n_config_total = 0
        self.n_config_excluded = 0

        print("~~~~~~")
        print(
            "Background:\n",
            xtal.pretty_json(self.config_enum.background.to_dict()),
            end="",
        )
        if hasattr(self.config_enum, "dof_space"):
            print(
                "DoFSpace:",
                xtal.pretty_json(self.config_enum.dof_space.to_dict()),
                end="",
            )
        elif hasattr(self.config_enum, "sites"):
            print("Sites:", self.config_enum.sites)
        sys.stdout.flush()

    def _finish_background(self):
        n_config_after = len(self.configuration_set)
        n_config_new = n_config_after - self.n_config_before
        print(
            f"{self.n_config_total} configurations "
            f"({n_config_new} new, {self.n_config_excluded} excluded by filter)"
        )
        print()
        sys.stdout.flush()
        self.n_config_before = n_config_after

    def check(self):
        """Call with each enumerated configuration"""
        if self.config_enum.enum_index != self.curr_enum_index:
            if self.curr_enum_index is not None:
                self._finish_background()
            self.curr_enum_index = self.config_enum.enum_index
            self._begin_background()

    def finish(self):
        """Call when enumeration is complete"""
        if self.curr_enum_index is not None:
            self._finish_background()
