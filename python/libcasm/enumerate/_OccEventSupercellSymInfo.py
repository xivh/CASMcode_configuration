import libcasm.clusterography
import libcasm.configuration
import libcasm.occ_events

from ._enumerate import (
    make_canonical_local_configuration,
)
from ._methods import make_occevent_suborbits
from ._OccEventPrimSymInfo import OccEventPrimSymInfo


def _get_suborbit_index(event, suborbits):
    for i_suborbit, suborbit in enumerate(suborbits):
        if event in suborbit:
            return i_suborbit
    return None


class OccEventSupercellSymInfo:
    """Information about an OccEvent with respect to a supercell factor group"""

    def __init__(
        self,
        prim_sym_info: OccEventPrimSymInfo,
        supercell: libcasm.configuration.Supercell,
    ):
        self.prim_sym_info = prim_sym_info
        """OccEventPrimSymInfo: OccEvent symmetry information with respect to the prim
        factor group."""

        self.supercell = supercell
        """libcasm.configuration.Supercell: The supercell."""

        suborbits = make_occevent_suborbits(
            supercell=supercell,
            occ_event=self.prim_sym_info.events[0],
        )
        self.suborbits = suborbits
        """list[list[libcasm.occ_events.OccEvent]]: The sub-orbits of OccEvent that are
        equivalent with respect to the supercell factor group."""

        self.n_suborbits = len(self.suborbits)
        """int: The number of sub-orbits of OccEvent in the supercell."""

        equivalent_index_to_suborbit_index = []
        suborbit_index_to_equivalent_index = [list() for _ in range(self.n_suborbits)]
        for i_equiv, event in enumerate(self.prim_sym_info.events):
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

    def event_group_rep(self, equivalent_index: int):
        """Get the SupercellSymOp representation for the invariant group of the
        specified equivalent event

        Parameters
        ----------
        equivalent_index: int
            The equivalent index of the event.

        Returns
        -------
        event_group_rep: list[libcasm.configuration.SupercellSymOp]
            The SupercellSymOp representation for the event invariant group.
        """
        invariant_group = self.prim_sym_info.invariant_groups[equivalent_index]
        return self.supercell.local_symgroup_rep(invariant_group)

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

        ## 3 steps:
        # - r1: Translate initial -> (0, initial[0])
        # - r2: Transform from (0, initial[1]) -> (0, final[1])
        # - r3: Translate from (0, final[1]) -> final
        # r_total = r3 * r2 * r1

        ##  Make r1: Translate initial -> (0, initial[0])
        r1 = libcasm.configuration.SupercellSymOp(
            supercell=self.supercell,
            supercell_factor_group_index=0,
            translation_index=initial[0],
        ).inverse()

        ## Make r2: Transform from (0, initial[1]) -> (0, final[1])

        # eq_i = op_i * eq_0 + trans_i
        # eq_f = op_f * eq_0 + trans_f

        # eq_f = op_f * (op_i.inv * (eq_i - trans_i) ) + trans_f
        # eq_f = r_total * eq_i

        # r_total = (t_f) * (op_f) * (op_i.inv) * (t_i.inv)
        # r2 = r2_c * r2_b * r2_a
        # r2_c = t_f, r2_b = op_f * (op_i.inv), r2_a = t_i.inv

        # transformation from {0, initial[1]} to {0, final[1]}
        initial_suborbit_index = self.equivalent_index_to_suborbit_index[initial[1]]
        final_suborbit_index = self.equivalent_index_to_suborbit_index[final[1]]
        if initial_suborbit_index != final_suborbit_index:
            raise ValueError(
                "Initial and final events are not equivalent by supercell factor group "
                "symmetry, so the transformation cannot be made in this supercell."
            )

        # Determine the prim fg operation for initial -> final
        prim = self.prim_sym_info.prim
        fg = prim.factor_group
        gen_indices = self.prim_sym_info.equivalent_generating_op_indices
        translations = self.prim_sym_info.translations
        scel_fg = self.supercell.factor_group

        # Prim fg indices for transforming initial -> final
        i_prototype_to_initial = gen_indices[initial[1]]
        i_initial_to_prototype = fg.inv(i_prototype_to_initial)
        i_prototype_to_final = gen_indices[final[1]]
        i_initial_to_final = fg.mul(i_prototype_to_final, i_initial_to_prototype)

        # The supercell_factor_group_index for transforming initial -> final
        # - Will raise if not found (should not occur, should be in same orbit)
        i_scel_initial_to_final = scel_fg.head_group_index.index(i_initial_to_final)

        r2_a = libcasm.configuration.SupercellSymOp(
            supercell=self.supercell,
            supercell_factor_group_index=0,
            translation_frac=-translations[initial[1]],
        ).inverse()
        r2_b = libcasm.configuration.SupercellSymOp(
            supercell=self.supercell,
            supercell_factor_group_index=i_scel_initial_to_final,
            translation_index=0,
        )
        r2_c = libcasm.configuration.SupercellSymOp(
            supercell=self.supercell,
            supercell_factor_group_index=0,
            translation_frac=translations[final[1]],
        )
        r2 = r2_c * r2_b * r2_a

        ##  Make r3: Translate from (0, final[1]) -> final
        r3 = libcasm.configuration.SupercellSymOp(
            supercell=self.supercell,
            supercell_factor_group_index=0,
            translation_index=final[0],
        )

        r_total = r3 * r2 * r1
        return r_total

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
        occ_event: libcasm.occ_events.OccEvent,
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
        return self.coordinate(occ_event=occ_event, supercell=self.supercell)

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
        return self.prim_sym_info.events[pos[1]] + unitcell

    def make_canonical_local_configuration(
        self,
        configuration: libcasm.configuration.Configuration,
        pos: tuple[int, int],
    ):
        """Make a canonical local configuration

        Parameters
        ----------
        configuration: libcasm.configuration.Configuration
            The configuration to transform.
        pos: tuple[int,int]
            The position of the event in the supercell of `configuration`, as
            `(unitcell_index, equivalent_index)`.

        Returns
        -------
        canonical_configuration: libcasm.configuration.Configuration
            The canonical local configuration for the specified event position.
        canonical_pos: tuple[int,int]
            The canonical position of the event in the supercell, as
            `(0, equivalent_index)`.
        """
        canonical_pos = self.canonical_pos(pos)
        canonical_configuration = make_canonical_local_configuration(
            configuration=self.supercell_rep(pos, canonical_pos) * configuration,
            event=self.prim_sym_info.events[canonical_pos[1]],
            event_group=self.event_group_rep(canonical_pos[1]),
        )
        return (canonical_configuration, canonical_pos)
