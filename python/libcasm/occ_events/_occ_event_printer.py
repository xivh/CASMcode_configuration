from typing import TextIO

from libcasm.occ_events._occ_events import (
    OccEvent,
    OccPosition,
    OccSystem,
)
from libcasm.xtal import (
    IntegralSiteCoordinate,
    cartesian_to_fractional,
)


class OccEventPrinter:
    """Pretty print :class:`~libcasm.occ_events.OccEvent`

    Example output:

    .. code-block:: Python

        Site Occupation:
            [0, 0, 0, 0]:  1 == B  ->  2 == Va
            [0, 1, 0, 0]:  2 == Va  ->  1 == B
        Trajectories:
            [[0, 0, 0, 0], 1] == B  ->  [[0, 1, 0, 0], 1] == B
            [[0, 1, 0, 0], 2] == Va  ->  [[0, 0, 0, 0], 2] == Va

    The "Site Occupation" format is:

    .. code-block:: Python

        r: p_i == n_i -> p_f == n_f

    Using the variables:

    - r: Site coordinate (includes atom component offset)
    - p: Occupation index
    - n: The "orientation name" of the occupant
    - _i, _f: Indicate initial and final, respectively


    For molecular occupant trajectories, the "Trajectories" format is:

    .. code-block:: Python

        [r_i, p_i] == n_i  ->  [r_f, p_f] == n_f

    For atom component trajectories, the "Trajectories" format is:

    .. code-block:: Python

        [r_i, p_i, a_i] == n_i, atom[a_i]=m_i  ->  [r_f, p_f, a_f] == n_f, atom[a_f]=m_f

    With the additional variables:

    - a: Atom position index
    - m: The name of the atom component

    Notes
    -----

    Including a note here.

    Parameters
    ----------
    f: TextIO
        Where to write text
    system: OccSystem
        Occupation system, for index conversions
    single_atom_occupant_as_mol: bool = True
        If True, print single atom occupant trajectories using the molecule trajectory
        format (i.e. do not show atom position index and atom name separately).
    coordinate_mode: str = 'integral'
        Mode for printing coordinates. Options are:

        - 'integral': Use :class:`~libcasm.xtal.IntegralSiteCoordinate` ([b, i, j, k])
        - 'cart': Use Cartesian coordinates
        - 'frac': Use fractional coordinates, with respect to the Prim lattice vectors


    .. rubric: Attributes

    f: TextIO
        Where to write text
    system: OccSystem
        Occupation system, for index conversions
    single_atom_occupant_as_mol: bool = True
        If True, print single atom occupant trajectories using the molecule trajectory
        format (i.e. do not show atom position index and atom name separately).
    coordinate_mode: str = 'integral'
        Mode for printing coordinates. Options are:

        - 'integral': Use :class:`~libcasm.xtal.IntegralSiteCoordinate` ([b, i, j, k])
        - 'cart': Use Cartesian coordinates
        - 'frac': Use fractional coordinates, with respect to the Prim lattice vectors

    """

    def __init__(
        self,
        f: TextIO,
        system: OccSystem,
        coordinate_mode: str = "integral",
        single_atom_occupant_as_mol: bool = True,
    ):
        """

        Parameters
        ----------
        f: TextIO
            Where to write text
        system: OccSystem
            Occupation system, for index conversions
        coordinate_mode: str = 'integral'
            Mode for printing coordinates. Options are:

            - 'integral': Use :class:`~libcasm.xtal.IntegralSiteCoordinate`
              ([b, i, j, k])
            - 'cart': Use Cartesian coordinates
            - 'frac': Use fractional coordinates, with respect to the Prim lattice
              vectors

        single_atom_occupant_as_mol: bool = True
            If True, print single atom occupant trajectories using the molecule
            trajectory format (i.e. do not show atom position index and atom name
            separately).
        """
        self.f = f
        self.system = system
        self.coordinate_mode = coordinate_mode.lower()
        self.single_atom_occupant_as_mol = single_atom_occupant_as_mol

        if self.coordinate_mode not in ["integral", "cart", "frac"]:
            raise Exception("Invalid coordinate_mode")

    def _coordinate(self, pos: OccPosition):
        """Return True if the molecule trajectory format should be used"""
        if self.coordinate_mode == "integral":
            return pos.integral_site_coordinate().to_list()
        else:
            r_cart = self.system.get_cartesian_coordinate(pos)
            if self.coordinate_mode == "cart":
                return r_cart.tolist()
            else:
                xtal_prim = self.system.xtal_prim()
                r_frac = cartesian_to_fractional(xtal_prim.lattice(), r_cart)
                return r_frac[:, 0].tolist()
        raise Exception("Invalid coordinate_mode")

    def _site_coordinate(self, site: IntegralSiteCoordinate):
        """Return True if the molecule trajectory format should be used"""
        if self.coordinate_mode == "integral":
            return site.to_list()
        else:
            xtal_prim = self.system.xtal_prim()
            if self.coordinate_mode == "cart":
                return site.coordinate_cart(xtal_prim).tolist()
            else:
                return site.coordinate_frac(xtal_prim).tolist()
        raise Exception("Invalid coordinate_mode")

    def _as_mol(self, pos: OccPosition):
        """Return True if the molecule trajectory format should be used"""
        if pos.is_atom():
            occupant = self.system.get_occupant(pos)
            n_atoms = len(occupant.atoms())
            if n_atoms == 1 and self.single_atom_occupant_as_mol:
                return True
            return False
        return True

    def _write_pos(self, pos: OccPosition):
        """Print cluster occupation (indices, then orientation names)

        Example output, if in reservoir:

            {chemical_name} (in reservoir)

        Example output, if molecule:

            {poslist} == {orientation_name}

        Example output, if atom component:

            {poslist} == {orientation_name}, atom[{atom_position_index}]={atom_name}

        Where 'poslist' is:

            [[b, i, j, k], {occupant_index}]

        Or:

            [[b, i, j, k], {occupant_index}, {atom_position_index}]

        Parameters
        ----------
        pos: ~libcasm.occ_events.OccPosition
            The occupant position

        """
        if pos.is_in_reservoir():
            chemical_name = self.system.get_chemical_name(pos)
            return self.f.write(f"{chemical_name} (in reservoir)")
        else:
            orientation_name = self.system.get_orientation_name(pos)
            poslist = [
                self._coordinate(pos),
                pos.occupant_index(),
            ]

            if self._as_mol(pos):
                return self.f.write(f"{poslist} == {orientation_name}")
            else:
                poslist.append(pos.atom_position_index())
                atom_position_index = pos.atom_position_index()
                atom_name = self.system.get_atom_name(pos)
                return self.f.write(
                    f"{poslist} == {orientation_name}, "
                    f"atom[{atom_position_index}]={atom_name}"
                )

    def _write_traj(self, pos0: OccPosition, pos1: OccPosition):
        """Print an occupant trajectory

        Example output, written to self.f:

            {pos0} -> {pos1}

        Or:

            {pos0} (no change)


        Parameters
        ----------
        cluster: ~libcasm.configuration.Cluster
            The cluster of sites

        occ: list[int]
            The occupation indices on the sites
        """
        self._write_pos(pos0)
        if pos0 == pos1:
            self.f.write(" (no change)")
        else:
            self.f.write("  ->  ")
            self._write_pos(pos1)

    def __call__(self, occ_event: OccEvent):
        """Write OccEvent to text output"""
        cluster = occ_event.cluster()
        occ_init = occ_event.initial_occupation()
        occ_final = occ_event.final_occupation()
        trajectories = occ_event.trajectories()
        occ_dof = self.system.xtal_prim().occ_dof()

        self.f.write("Site Occupation: \n")
        for i, site in enumerate(cluster):
            b = site.sublattice()
            name_init = occ_dof[b][occ_init[i]]
            name_final = occ_dof[b][occ_final[i]]
            OccPosition.occupant(site, occ_final[i])
            self.f.write("    ")
            self.f.write(str(self._site_coordinate(site)))
            self.f.write(":  " + str(occ_init[i]) + " == " + name_init)
            self.f.write("  ->  " + str(occ_final[i]) + " == " + name_final)
            self.f.write("\n")

        self.f.write("Trajectories:\n")
        for pos0, pos1 in trajectories:
            self.f.write("    ")
            self._write_traj(pos0, pos1)
            self.f.write("\n")
