import numpy as np

import libcasm.configuration as casmconfig
import libcasm.occ_events as occ_events

from ._local_configuration import (
    _make_canonical_local_configuration_about_event,
)
from ._methods import make_occevent_suborbits
from ._OccEventPrimSymInfo import OccEventPrimSymInfo


def _copy_apply_prim_factor_group_index(
    prim_factor_group_index: int,
    prim: casmconfig.Prim,
    occ_event: occ_events.OccEvent,
):
    """Copy and apply a supercell symmetry operation to an OccEvent, standardize it, and
    bring it within the supercell

    Parameters
    ----------
    prim_factor_group_index : int
        The index of the factor group operation to apply
    prim : libcasm.configuration.Prim
        The Prim.
    occ_event : occ_events.OccEvent
        The OccEvent to transform.

    Returns
    -------
    transformed_occ_event : occ_events.OccEvent
        The transformed OccEvent.
    """
    i_fg = prim_factor_group_index
    rep = occ_events.OccEventRep(
        integral_site_coordinate_rep=prim.integral_site_coordinate_symgroup_rep[i_fg],
        occupant_rep=prim.occ_symgroup_rep[i_fg],
        atom_position_rep=prim.atom_position_symgroup_rep[i_fg],
    )
    return rep * occ_event


def _make_consistent_suborbits(suborbits, event_prim_info):
    """Ensure `suborbits` events are consistent with `event_prim_info.events`"""
    default_supercell = casmconfig.Supercell(
        prim=event_prim_info.prim,
        transformation_matrix_to_super=np.eye(3, dtype="int"),
    )

    suborbit_event_indices = []
    for suborbit in suborbits:
        indices = []
        for suborbit_event in suborbit:
            pos = event_prim_info.coordinate(
                occ_event=suborbit_event,
                supercell=default_supercell,
            )
            indices.append(pos[1])
        indices.sort()
        suborbit_event_indices.append(indices)
    suborbit_event_indices.sort(key=lambda x: x[0])

    suborbits_out = []
    for indices in suborbit_event_indices:
        suborbit_out = [event_prim_info.events[i] for i in indices]
        suborbits_out.append(suborbit_out)
    return suborbits_out


def _get_suborbit_index(event, suborbits):
    for i_suborbit, suborbit in enumerate(suborbits):
        for suborbit_event in suborbit:
            if event == suborbit_event:
                return i_suborbit
        if event in suborbit:
            return i_suborbit
    return None


class OccEventSupercellSymInfo:
    """Information about an OccEvent with respect to a supercell factor group"""

    def __init__(
        self,
        event_prim_info: OccEventPrimSymInfo,
        supercell: casmconfig.Supercell,
    ):
        self.event_prim_info = event_prim_info
        """OccEventPrimSymInfo: OccEvent symmetry information with respect to the prim
        factor group."""

        self.supercell = supercell
        """libcasm.configuration.Supercell: The supercell."""

        suborbits = make_occevent_suborbits(
            supercell=supercell,
            occ_event=self.event_prim_info.events[0],
        )
        suborbits = _make_consistent_suborbits(suborbits, self.event_prim_info)
        self.suborbits = suborbits
        """list[list[libcasm.occ_events.OccEvent]]: The sub-orbits of OccEvent that are
        equivalent with respect to the supercell factor group."""

        self.n_suborbits = len(self.suborbits)
        """int: The number of sub-orbits of OccEvent in the supercell."""

        equivalent_index_to_suborbit_index = []
        suborbit_index_to_equivalent_index = [list() for _ in range(self.n_suborbits)]
        for i_equiv, event in enumerate(self.event_prim_info.events):
            i_suborbit = _get_suborbit_index(event, self.suborbits)
            if i_suborbit is None:
                raise Exception("Failed to find suborbit")
            equivalent_index_to_suborbit_index.append(i_suborbit)
            suborbit_index_to_equivalent_index[i_suborbit].append(i_equiv)

        self.equivalent_index_to_suborbit_index = equivalent_index_to_suborbit_index
        """list[int]: Maps equivalent index to suborbit index, according to
        `suborbit_index = equivalent_index_to_suborbit_index[equivalent_index]`"""

        self.suborbit_index_to_equivalent_index = suborbit_index_to_equivalent_index
        """list[list[int]]: Maps suborbit index to equivalent indices, where
        `equivalent_indices = suborbit_index_to_equivalent_index[suborbit_index]`
        is a list[int] of the equivalent indices that make up a suborbit of OccEvent
        that are symmetrically equivalent with respect to the supercell factor group."""

    def suborbit_prototype(self, suborbit_index: int):
        """Get the prototype event for the suborbit at the given index

        Parameters
        ----------
        suborbit_index: int
            The index of the suborbit.

        Returns
        -------
        prototype: libcasm.occ_events.OccEvent
            The prototype event for the suborbit.
        """
        return self.suborbits[suborbit_index][0]

    def suborbit_prototype_equivalent_index(self, suborbit_index: int):
        """Get the equivalent index of the prototype event for the suborbit at the
        given index

        Parameters
        ----------
        suborbit_index: int
            The index of the suborbit.

        Returns
        -------
        equivalent_index: int
            The equivalent index of the prototype event for the suborbit.
        """
        return self.suborbit_index_to_equivalent_index[suborbit_index][0]

    def event_group_rep(self, pos: tuple[int, int]):
        """Get the SupercellSymOp representation for the invariant group of the
        specified equivalent event.

        Parameters
        ----------
        pot: tuple[int, int]
            The event position as a `(unitcell_index, equivalent_index)` pair.

        Returns
        -------
        event_group_rep: list[libcasm.configuration.SupercellSymOp]
            The SupercellSymOp representation for the event invariant group.
        """
        t1 = casmconfig.SupercellSymOp(
            supercell=self.supercell,
            supercell_factor_group_index=0,
            translation_index=pos[0],
        )
        t1_inv = t1.inverse()
        invariant_group = self.event_prim_info.invariant_groups[pos[1]]
        g = self.supercell.local_symgroup_rep(invariant_group)
        return [t1 * g_i * t1_inv for g_i in g]

    def supercell_rep(
        self,
        initial: tuple[int, int],
        final: tuple[int, int],
    ):
        """Get the supercell representation for the transformation from initial event
        indices to final event indices

        Parameters
        ----------
        initial: tuple[int, int]
            The initial `(unitcell_index, equivalent_index)` pair.
        final: tuple[int, int]
            The final `(unitcell_index, equivalent_index)` pair.

        Returns
        -------
        supercell_rep: libcasm.configuration.SupercellSymOp
            The supercell representation for transforming from initial event to final
            event.

        Raises
        ------
        ValueError
            If the initial and final events are not equivalent by supercell factor group
            symmetry.
        """
        #   t2 * r_p2f * r_i2p * t1 * e_i
        # r_initial

        # transformation from {0, initial[1]} to {0, final[1]}
        initial_suborbit_index = self.equivalent_index_to_suborbit_index[initial[1]]
        final_suborbit_index = self.equivalent_index_to_suborbit_index[final[1]]
        if initial_suborbit_index != final_suborbit_index:
            raise ValueError(
                "Initial and final events are not equivalent by supercell factor group "
                "symmetry, so the transformation cannot be made in this supercell."
            )

        # Determine a prim fg operation from {0, initial[1]} to {0, final[1]}
        prim = self.event_prim_info.prim
        fg = prim.factor_group
        gen_indices = self.event_prim_info.equivalent_generating_op_indices
        scel_fg = self.supercell.factor_group

        # Find one prim fg operation from {0, initial[1]} to {0, final[1]}
        # (might not be in supercell factor group)
        i_prototype_to_initial = gen_indices[initial[1]]
        i_initial_to_prototype = fg.inv(i_prototype_to_initial)
        i_prototype_to_final = gen_indices[final[1]]
        i_initial_to_final_0 = fg.mult(i_prototype_to_final, i_initial_to_prototype)

        # Find an operation from {0, initial[1]} to {0, final[1]}
        # in the supercell factor group
        group = self.event_prim_info.invariant_groups[initial[1]]
        i_initial_to_final = None
        for g_i in group.head_group_index:
            i_test = fg.mult(i_initial_to_final_0, g_i)
            if i_test in scel_fg.head_group_index:
                i_initial_to_final = i_test
                break
        if i_initial_to_final is None:
            raise ValueError("Failed to find transformation in supercell factor group.")

        # Get the event OccEventRep
        event_rep = self.event_prim_info.occevent_symgroup_rep[i_initial_to_final]

        # Transform event by OccEventRep & get its position
        event_2 = event_rep * self.event(initial)
        pos_2 = self.coordinate(event_2)

        # Find translation (fractional wrt prim) to the final position
        frac_2 = self.supercell.unitcell_index_converter.unitcell(pos_2[0])
        frac_final = self.supercell.unitcell_index_converter.unitcell(final[0])

        # Combine the operations and return the SupercellSymOp
        r = casmconfig.SupercellSymOp(
            supercell=self.supercell,
            supercell_factor_group_index=scel_fg.head_group_index.index(
                i_initial_to_final
            ),
            translation_index=0,
        )
        t = casmconfig.SupercellSymOp(
            supercell=self.supercell,
            supercell_factor_group_index=0,
            translation_frac=frac_final - frac_2,
        )
        return t * r

    def canonical_pos(self, pos: tuple[int, int]):
        """Return the canonical equivalent event position

        Parameters
        ----------
        pos: tuple[int,int]
            The position of the event in the supercell, as
            `(unitcell_index, equivalent_index)`.

        Returns
        -------
        canonical_pos: tuple[int,int]
            The canonical position of the event in the supercell, as
            `(0, equivalent_index)`.
        """
        i_suborbit = self.equivalent_index_to_suborbit_index[pos[1]]
        return (0, self.suborbit_index_to_equivalent_index[i_suborbit][0])

    def coordinate(
        self,
        occ_event: occ_events.OccEvent,
    ):
        """
        Determine the coordinates `(unitcell_index, equivalent_index)` of a OccEvent,
        with `equivalent_index` referring to `self.events`, in `self.supercell`.

        Parameters
        ----------
        occ_event : libcasm.occ_events.OccEvent
            Input OccEvent, to find the coordinates of

        Returns
        -------
        (unitcell_index, equivalent_index) : tuple[int, int]
            The coordinates (unitcell_index, equivalent_index) of the
            input OccEvent, with `equivalent_index` referring to `self.events`,
            in `self.supercell`.

        Raises
        ------
        Exception
            If no match can be found, indicating the input OccEvent is not equivalent
            to any of the OccEvent in `self.events` up to a translation.
        """
        return self.event_prim_info.coordinate(
            occ_event=occ_event, supercell=self.supercell
        )

    def standardized_and_within(
        self,
        occ_event: occ_events.OccEvent,
    ):
        """

        Parameters
        ----------
        occ_event : libcasm.occ_events.OccEvent
            The OccEvent.

        Returns
        -------
        standardized_occ_event : libcasm.sym_info.SymGroup
            A copy of the input OccEvent, standardized and brought within the supercell.
        """
        pos = self.coordinate(occ_event)
        trans = self.supercell.unitcell_index_converter.unitcell(pos[0])
        return self.event_prim_info.events[pos[1]] + trans

    def copy_apply_supercell_symop(
        self,
        op: casmconfig.SupercellSymOp,
        occ_event: occ_events.OccEvent,
    ):
        """Copy and apply a supercell symmetry operation to an OccEvent, standardize it,
        and bring it within the supercell

        Parameters
        ----------
        op : libcasm.configuration.SupercellSymOp
            The supercell symmetry operation.
        occ_event : occ_events.OccEvent
            The OccEvent to transform.

        Returns
        -------
        transformed_occ_event : occ_events.OccEvent
            The transformed OccEvent.
        """
        prim = self.supercell.prim
        i_fg = op.prim_factor_group_index()
        site_rep = prim.integral_site_coordinate_symgroup_rep[i_fg]
        rep = occ_events.OccEventRep(
            integral_site_coordinate_rep=site_rep,
            occupant_rep=prim.occ_symgroup_rep[i_fg],
            atom_position_rep=prim.atom_position_symgroup_rep[i_fg],
        )
        return self.standardized_and_within(rep * occ_event + op.translation_frac())

    def event(
        self,
        pos: tuple[int, int],
    ):
        """Return the event at the specified position in the supercell

        Parameters
        ----------
        pos: tuple[int,int]
            The position of the event in the supercell, as
            `(unitcell_index, equivalent_index)`.

        Returns
        -------
        event: libcasm.occ_events.OccEvent
            The event at the specified position in the supercell.
        """
        unitcell = self.supercell.unitcell_index_converter.unitcell(pos[0])
        return self.event_prim_info.events[pos[1]] + unitcell

    def make_canonical_local_configuration(
        self,
        configuration: casmconfig.Configuration,
        pos: tuple[int, int],
        in_canonical_pos: bool = True,
        in_canonical_supercell: bool = False,
        apply_event_occupation: bool = True,
    ):
        """Make a canonical local configuration

        Parameters
        ----------
        configuration: libcasm.configuration.Configuration
            The configuration to transform.
        pos: tuple[int,int]
            The position of the event in the supercell of `configuration`, as
            `(unitcell_index, equivalent_index)`.
        in_canonical_pos: bool = True
            If True, transform `configuration` and `pos` to put `pos` in the
            canonical position in the supercell. Else, keep `pos` in its current
            position and only transform `configuration`.
        in_canonical_supercell: bool = False
            If True, first transform the `configuration` and `pos` to put them into the
            canonical equivalent supercell. If True, `in_canonical_pos` must also be
            True.
        apply_event_occupation: bool = True
            If True, apply the occupation of the event to the configuration. If False,
            maintain the current configuration occupation.

        Returns
        -------
        final_configuration: libcasm.configuration.Configuration
            The canonical local configuration for the requested event position.
        final_pos: tuple[int,int]
            The canonical position of the event in the supercell, as
            `(0, equivalent_index)`.
        final_event_supercell_info: OccEventSupercellSymInfo
            The canonical supercell event symmetry info.
        """
        if self.supercell != configuration.supercell:
            raise ValueError(
                "Error in "
                "OccEventSupercellSymInfo.make_canonical_local_configuration: "
                "self.supercell is not configuration.supercell"
            )

        this_is_canonical_supercell = casmconfig.is_canonical_supercell(
            supercell=self.supercell
        )

        if in_canonical_supercell is True and not this_is_canonical_supercell:
            if in_canonical_pos is False:
                raise ValueError(
                    "Error in "
                    "OccEventSupercellSymInfo.make_canonical_local_configuration: "
                    "in_canonical_supercell=True requires in_canonical_pos=True"
                )

            # Find i_fg, the index of the prim factor group operation that transforms
            # the current configuration's supercell to the canonical supercell
            supercell = configuration.supercell
            prim = supercell.prim
            S1 = supercell.superlattice
            canonical_supercell = casmconfig.make_canonical_supercell(supercell)
            S2 = canonical_supercell.superlattice
            is_superlat, T, i_fg = S2.is_equivalent_superlattice_of(
                S1, prim.factor_group.elements
            )
            if not is_superlat:
                raise ValueError(
                    "Error in "
                    "OccEventSupercellSymInfo.make_canonical_local_configuration: "
                    "failed to find transformation to canonical supercell"
                )

            # Apply i_fg-th operation to configuration and event
            tmp_configuration = casmconfig.copy_transformed_configuration(
                prim_factor_group_index=i_fg,
                translation=np.array([0, 0, 0]),
                motif=configuration,
                supercell=canonical_supercell,
            )
            tmp_event = _copy_apply_prim_factor_group_index(
                prim_factor_group_index=i_fg,
                prim=prim,
                occ_event=self.event(pos),
            )

            # Get the position of the transformed event in the canonical supercell
            _final_event_supercell_info = OccEventSupercellSymInfo(
                event_prim_info=self.event_prim_info,
                supercell=canonical_supercell,
            )
            # tmp_event = _final_event_supercell_info.standardized_within(tmp_event)
            tmp_pos = _final_event_supercell_info.coordinate(tmp_event)

            # Make the canonical local configuration in the canonical supercell
            (
                _final_config,
                _final_pos,
                _result_event_supercell_info,
            ) = _final_event_supercell_info.make_canonical_local_configuration(
                configuration=tmp_configuration,
                pos=tmp_pos,
                in_canonical_pos=True,
                in_canonical_supercell=True,
                apply_event_occupation=apply_event_occupation,
            )
            if _result_event_supercell_info != _final_event_supercell_info:
                raise ValueError(
                    "Error in "
                    "OccEventSupercellSymInfo.make_canonical_local_configuration: "
                    "result_event_supercell_info != final_event_supercell_info"
                )

            # Return result
            return (_final_config, _final_pos, _final_event_supercell_info)

        if in_canonical_supercell is True and not this_is_canonical_supercell:
            raise ValueError(
                "Error in "
                "OccEventSupercellSymInfo.make_canonical_local_configuration: "
                "in_canonical_supercell=True and self.supercell is not canonical"
            )

        # Transform to put event in canonical position
        if in_canonical_pos:
            final_pos = self.canonical_pos(pos)
            tmp_configuration = self.supercell_rep(pos, final_pos) * configuration
        else:
            final_pos = pos
            tmp_configuration = configuration

        if apply_event_occupation:
            final_configuration = _make_canonical_local_configuration_about_event(
                configuration=tmp_configuration,
                event=self.event(final_pos),
                event_group=self.event_group_rep(final_pos),
            )
        else:
            final_configuration = casmconfig.make_canonical_configuration(
                configuration=tmp_configuration,
                in_canonical_supercell=False,
                subgroup=self.event_group_rep(final_pos),
            )
        return (final_configuration, final_pos, self)
