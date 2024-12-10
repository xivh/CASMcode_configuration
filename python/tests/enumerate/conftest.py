import sys
from typing import Callable, Optional, Union

import numpy as np
import pytest

import libcasm.configuration as casmconfig
import libcasm.enumerate as casmenum
import libcasm.local_configuration as casmlocal
import libcasm.occ_events as occ_events
import libcasm.xtal as xtal
import libcasm.xtal.prims as xtal_prims


@pytest.fixture
def fcc_1NN_A_Va_event():
    from libcasm.occ_events import OccEvent, OccPosition

    xtal_prim = xtal_prims.FCC(r=1.0, occ_dof=["A", "B", "Va"])

    site1 = xtal.IntegralSiteCoordinate(sublattice=0, unitcell=[0, 0, 0])
    site2 = xtal.IntegralSiteCoordinate(sublattice=0, unitcell=[1, 0, 0])

    A_occ_index = 0
    Va_occ_index = 2

    A_initial_pos = OccPosition.molecule(site1, A_occ_index)
    A_final_pos = OccPosition.molecule(site2, A_occ_index)
    Va_initial_pos = OccPosition.molecule(site2, Va_occ_index)
    Va_final_pos = OccPosition.molecule(site1, Va_occ_index)

    return (
        xtal_prim,
        OccEvent([[A_initial_pos, A_final_pos], [Va_initial_pos, Va_final_pos]]),
    )


@pytest.fixture
def fcc_1NN_A_Va_event_L12(fcc_1NN_A_Va_event):
    xtal_prim, event = fcc_1NN_A_Va_event

    # r = 1.0
    # a = 4 * r / math.sqrt(2)
    prim = casmconfig.Prim(xtal_prim)
    system = occ_events.OccSystem(xtal_prim=prim.xtal_prim)
    event_info = casmlocal.OccEventSymInfo.init(
        prim=prim,
        system=system,
        prototype_event=event,
    )

    # L12 ordering (A3, B1)
    supercell = casmconfig.Supercell(
        prim=prim,
        transformation_matrix_to_super=np.array(
            [
                [-1, 1, 1],
                [1, -1, 1],
                [1, 1, -1],
            ],
            dtype="int",
        ),
    )
    config = casmconfig.Configuration(supercell)
    config.set_occ(0, 1)

    return (prim, event, event_info, config)


class ConfigEnumLocalTestData:
    def __init__(
        self,
        prim: casmconfig.Prim,
        event: occ_events.OccEvent,
        event_info: casmlocal.OccEventSymInfo,
    ):
        self.prim = prim
        self.event = event
        self.event_info = event_info

    def set_background(
        self,
        background: casmconfig.Configuration,
    ):
        self.background = background

    def set_supercell(
        self,
        supercell: Optional[casmconfig.Supercell] = None,
        transformation_matrix_to_super: Optional[np.ndarray] = None,
    ):
        if supercell is None:
            supercell = casmconfig.Supercell(
                prim=self.prim,
                transformation_matrix_to_super=transformation_matrix_to_super,
            )
        self.supercell = supercell

    def make_local_cluster_specs(
        self,
        n_orbits: Optional[list[int]] = None,
        cutoff_radius: Optional[list[float]] = None,
        make_all_possible_orbits: bool = False,
    ):
        self.local_cluster_specs = casmenum.make_first_n_orbits_cluster_specs(
            prim=self.prim,
            phenomenal=self.event,
            n_orbits=n_orbits,
            cutoff_radius=cutoff_radius,
            make_all_possible_orbits=make_all_possible_orbits,
        )

    def check_enum(
        self,
        verbose: bool = False,
        fix: str = "background",
        pos: Union[tuple[int, int], list[tuple[int, int]], None] = None,
        orbits: Optional[list[int]] = None,
        neighborhood_from_orbits: Optional[list[int]] = None,
        skip_result_f: Callable[
            [casmenum.ConfigEnumLocalOccupationsResult], bool
        ] = lambda result: False,
        n_expected: Optional[int] = None,
    ):
        """Setup and run test enumeration of local environments.

        Notes
        -----

        Will set:

        - `self.supercells` to a `SupercellSet` object.
        - `self.local_config_list` to a `LocalConfigurationList` object.
        - `self.canonical_local_config_list` to a `LocalConfigurationList` object.
        - `self.config_enum` to a `ConfigEnumLocalOccupations` object.

        Will assert:

        - The number of local configurations is equal to `n_expected` (if `n_expected`
          is not None).

        Parameters
        ----------
        verbose : bool, optional
            Print verbose output, by default False
        fix : Optional[str] = "background"
            See casmenum.ConfigEnumLocalOccupations.by_cluster_specs.
        pos : Union[tuple[int, int], list[tuple[int, int]], None] = None
            See casmenum.ConfigEnumLocalOccupations.by_cluster_specs.
        orbits : Optional[list[int]] = None
            See casmenum.ConfigEnumLocalOccupations.by_cluster_specs.
        neighborhood_from_orbits : Optional[list[int]] = None
            See casmenum.ConfigEnumLocalOccupations.by_cluster_specs.
        skip_result_f : Callable[[casmenum.ConfigEnumLocalOccupationsResult], bool]
            Optional filter function. If True, skip the result; If False, add it to
            `self.local_config_list` and `self.canonical_local_config_list` if not
            already present.
        n_expected : Optional[int] = None
            Expected number of local configurations. If not None, assert that the
            number of local configurations is equal to `n_expected`.

        """
        self.supercells = casmconfig.SupercellSet(prim=self.prim)
        """SupercellSet: Set of supercells for enumeration."""

        self.local_config_list = casmlocal.LocalConfigurationList(
            event_info=self.event_info,
        )
        """libcasm.local_configuration.LocalConfigurationList: List of local 
        configurations, as enumerated."""

        self.canonical_local_config_list = casmlocal.LocalConfigurationList(
            event_info=self.event_info,
        )
        """libcasm.local_configuration.LocalConfigurationList: List of local 
        configurations, made canonical with respect to the event invariant group for 
        the initial unperturbed configuration."""

        self.fully_canonical_local_config_list = casmlocal.LocalConfigurationList(
            event_info=self.event_info
        )
        """libcasm.local_configuration.LocalConfigurationList: List of local 
        configurations, made canonical with respect to the prim factor group."""

        self.config_enum = casmenum.ConfigEnumLocalOccupations(
            event_info=self.event_info,
            supercell_set=self.supercells,
            verbose=verbose,
        )

        for result in self.config_enum.by_cluster_specs(
            background=self.background,
            supercell=self.supercell,
            cluster_specs=self.local_cluster_specs,
            fix=fix,
            pos=pos,
            orbits=orbits,
            neighborhood_from_orbits=neighborhood_from_orbits,
        ):
            if skip_result_f(result):
                continue
            if (
                result.canonical_local_configuration
                not in self.canonical_local_config_list
            ):
                if verbose:
                    self.print_result_summary(result=result, existing=False)
                self.local_config_list.append(result.local_configuration)
                self.canonical_local_config_list.append(
                    result.canonical_local_configuration
                )

                # Make fully canonical LocalConfiguration
                x = result.local_configuration
                _list = self.fully_canonical_local_config_list
                x_fully_canonical = casmlocal.make_canonical_local_configuration(
                    x,
                    in_canonical_pos=True,
                    in_canonical_supercell=True,
                    apply_event_occupation=True,
                )
                self.fully_canonical_local_config_list.append(x_fully_canonical)
            elif verbose:
                self.print_result_summary(result=result, existing=True)

        if verbose:
            print("Total # Local configurations:", len(self.local_config_list))
            print()
            sys.stdout.flush()

        if n_expected is not None:
            assert len(self.local_config_list) == n_expected

    def print_result_summary(
        self,
        result: casmenum.ConfigEnumLocalOccupationsResult,
        existing: bool,
    ):
        print(
            f"- initial={result.i_initial}, "
            f"orbit={result.i_local_orbit}, "
            f"sites={result.sites}, "
            f"sublattices={result.sublattices}, "
            f"asymmetric_units={result.asymmetric_units}, "
            f"initial={result.initial_occupation}, ",
            f"final={result.final_occupation}",
            end="",
        )
        if existing:
            print(" (existing)")
        else:
            print()
        sys.stdout.flush()


@pytest.fixture
def L12_1NN_A_Va_local_enum(fcc_1NN_A_Va_event_L12):
    prim, event, event_info, config = fcc_1NN_A_Va_event_L12
    local_enum = ConfigEnumLocalTestData(prim=prim, event=event, event_info=event_info)
    local_enum.set_background(config)
    return local_enum


@pytest.fixture
def lowsym_Hstrain_prim():
    return xtal.Prim(
        lattice=xtal.Lattice(
            np.array(
                [
                    [1.0, 0.3, 0.4],  # a
                    [0.0, 1.2, 0.5],  # b
                    [0.0, 0.0, 1.4],  # c
                ]
            ).transpose()
        ),
        coordinate_frac=np.array(
            [
                [0.0, 0.0, 0.0],
                [0.4, 0.5, 0.6],
                [0.24, 0.25, 0.23],
            ]
        ).transpose(),
        occ_dof=[["A"], ["A"], ["A"]],
        global_dof=[xtal.DoFSetBasis("Hstrain")],
    )


@pytest.fixture
def lowsym_disp_prim():
    return xtal.Prim(
        lattice=xtal.Lattice(
            np.array(
                [
                    [1.0, 0.3, 0.4],  # a
                    [0.0, 1.2, 0.5],  # b
                    [0.0, 0.0, 1.4],  # c
                ]
            ).transpose()
        ),
        coordinate_frac=np.array(
            [
                [0.0, 0.0, 0.0],
                [0.4, 0.5, 0.6],
                [0.24, 0.25, 0.23],
            ]
        ).transpose(),
        occ_dof=[["A"], ["A"], ["A"]],
        local_dof=[
            [xtal.DoFSetBasis("disp")],
            [xtal.DoFSetBasis("disp")],
            [xtal.DoFSetBasis("disp")],
        ],
    )


@pytest.fixture
def lowsym_occ_prim():
    return xtal.Prim(
        lattice=xtal.Lattice(
            np.array(
                [
                    [1.0, 0.3, 0.4],  # a
                    [0.0, 1.2, 0.5],  # b
                    [0.0, 0.0, 1.4],  # c
                ]
            ).transpose()
        ),
        coordinate_frac=np.array(
            [
                [0.0, 0.0, 0.0],
                [0.4, 0.5, 0.6],
                [0.24, 0.25, 0.23],
            ]
        ).transpose(),
        occ_dof=[["A", "B"], ["A", "B"], ["A", "B"]],
    )


@pytest.fixture
def lowsym_Hstrain_disp_prim():
    return xtal.Prim(
        lattice=xtal.Lattice(
            np.array(
                [
                    [1.0, 0.3, 0.4],  # a
                    [0.0, 1.2, 0.5],  # b
                    [0.0, 0.0, 1.4],  # c
                ]
            ).transpose()
        ),
        coordinate_frac=np.array(
            [
                [0.0, 0.0, 0.0],
                [0.4, 0.5, 0.6],
                [0.24, 0.25, 0.23],
            ]
        ).transpose(),
        occ_dof=[["A"], ["A"], ["A"]],
        local_dof=[
            [xtal.DoFSetBasis("disp")],
            [xtal.DoFSetBasis("disp")],
            [xtal.DoFSetBasis("disp")],
        ],
        global_dof=[xtal.DoFSetBasis("Hstrain")],
    )
