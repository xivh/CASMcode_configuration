import sys
from typing import Optional, Union

import libcasm.clusterography as casmclust
import libcasm.configuration as casmconfig
import libcasm.local_configuration as casmlocal
import libcasm.xtal as xtal

from ._enumerate import (
    make_distinct_local_cluster_sites,
    make_distinct_local_perturbations,
)
from ._make_distinct_super_configurations import (
    make_distinct_super_configurations,
)
from ._point_defect_supercell_methods import (
    has_required_sites,
    make_required_sites,
)


def make_distinct_local_configurations(
    background: casmconfig.Configuration,
    event_info: casmlocal.OccEventSymInfo,
    fix: str = "background",
):
    """Return the distinct LocalConfiguration that can be generated as a combination of
    a background configuration and event.

    Parameters
    ----------
    background : casmconfig.Configuration
        A primitive background configuration.
    event_info : libcasm.local_configuration.OccEventSymInfo
        The OccEventSymInfo for the event.
    fix : str = "background"
        The object to be fixed during the method. Can be "background" or "event". If
        "background", the background is fixed and the event is varied. If "event", then
        the event is fixed and the background is varied. Note that if the event is
        fixed the resulting LocalConfiguration may have different supercells.

    Returns
    -------
    local_configurations : list[libcasm.local_configuration.LocalConfiguration]
        The symmetrically distinct LocalConfigurations.
    """
    if not casmconfig.is_primitive_configuration(background):
        raise ValueError(
            "Error in make_distinct_local_configurations: "
            "`background` must be a primitive configuration."
        )

    if fix == "background":
        background_group = casmconfig.make_invariant_subgroup(configuration=background)
        event_supercell_info = event_info.get_event_supercell_info(background.supercell)

        local_configurations = []
        for i_suborbit in range(event_supercell_info.n_suborbits):
            i_event = event_supercell_info.suborbit_index_to_equivalent_index[
                i_suborbit
            ][0]
            pos = (0, i_event)
            event_init = event_supercell_info.event(pos)
            event_group_rep = event_supercell_info.event_group_rep(pos)

            def _is_suborbit_generator_op(rep_init):
                for event_op in event_group_rep:
                    for background_op in background_group:
                        if background_op * rep_init * event_op < rep_init:
                            return False
                return True

            rep_it = casmconfig.SupercellSymOp.begin(background.supercell)
            rep_end = casmconfig.SupercellSymOp.end(background.supercell)

            while rep_it != rep_end:
                if _is_suborbit_generator_op(rep_it):
                    event_final = event_supercell_info.copy_apply_supercell_symop(
                        op=rep_it,
                        occ_event=event_init,
                    )
                    local_configurations.append(
                        casmlocal.LocalConfiguration(
                            pos=event_supercell_info.coordinate(event_final),
                            configuration=background,
                            event_info=event_info,
                        )
                    )
                rep_it.next()
        return local_configurations

    elif fix == "event":
        prim = background.supercell.prim

        superduperlattice = xtal.make_superduperlattice(
            lattices=[background.supercell.superlattice],
            mode="fully_commensurate",
            point_group=prim.crystal_point_group.elements,
        )
        superdupercell = casmconfig.Supercell(
            prim=prim,
            transformation_matrix_to_super=xtal.make_transformation_matrix_to_super(
                superlattice=superduperlattice,
                unit_lattice=prim.xtal_prim.lattice(),
            ),
        )
        super_backgrounds = casmconfig.make_all_super_configurations(
            motif=background,
            supercell=superdupercell,
        )
        event_supercell_info = event_info.get_event_supercell_info(superdupercell)
        event_group = event_supercell_info.event_group_rep((0, 0))
        distinct_backgrounds = []
        for x in super_backgrounds:
            if casmconfig.is_canonical_configuration(
                configuration=x,
                subgroup=event_group,
            ):
                prim_config = casmconfig.make_primitive_configuration(configuration=x)
                distinct_backgrounds.append(prim_config.copy())

        local_configurations = []
        for _background in distinct_backgrounds:
            local_configurations.append(
                casmlocal.LocalConfiguration(
                    pos=(0, 0),
                    configuration=_background,
                    event_info=event_info,
                )
            )

        return local_configurations


def _make_local_orbits_data(
    local_orbits: list[list[list[casmclust.Cluster]]],
    xtal_prim: xtal.Prim,
    phenomenal_clusters: list[casmclust.Cluster],
):
    """Represent local-cluster orbits for a list of events as a `list[list[list[dict]]]`

    Parameters
    ----------
    local_orbits : list[list[list[libcasm.clusterography.Cluster]]]
        The local-cluster orbits, for each event, where
        `local_orbits[i_event][i_local_orbit][i_cluster]` is the
        `i_cluster`-th cluster in the `i_local_orbit`-th local-cluster orbit
        about the `i_event`-th event.
    xtal_prim : libcasm.xtal.Prim
        The prim, used to calculate site distances.
    phenomenal_clusters : list[libcasm.clusterography.Cluster]
        The phenomenal clusters, for each event, where
        `phenomenal_clusters[i_event]` is the phenomenal cluster for the
        `i_event`-th event.

    Returns
    -------
    local_orbits_data : list[list[list[dict]]]
        The local-cluster orbits for each event, where
        `local_orbits_data[i_event][i_local_orbit][i_cluster]` is the
        `i_cluster`-th cluster in the `i_local_orbit`-th local-cluster orbit
        about the `i_event`-th event.
    """
    if local_orbits is not None:
        local_orbits_data = [
            [
                [
                    cluster.to_dict(
                        xtal_prim=xtal_prim,
                        phenomenal=phenomenal_clusters[i_event],
                    )
                    for cluster in local_orbit
                ]
                for local_orbit in _local_orbits
            ]
            for i_event, _local_orbits in enumerate(local_orbits)
        ]
    else:
        local_orbits_data = None
    return local_orbits_data


class ConfigEnumLocalOccupationsReference:
    """Data structure for storing contextual information common for results from
    :class:`~libcasm.local_configuration.ConfigEnumLocalOccupations`.

    Reference is used to store information about the background configurations
    that is common to all results from the enumeration of local perturbations.
    Individual results are returned as instances of :class:`Result`.
    """

    def __init__(
        self,
        event_info: casmlocal.OccEventSymInfo,
        event_supercell_info: casmlocal.OccEventSupercellSymInfo,
        initial: list[casmlocal.LocalConfiguration],
        prim_local_orbits: list[list[list[casmclust.Cluster]]] = None,
    ):
        if not isinstance(event_info, casmlocal.OccEventSymInfo):
            raise ValueError(
                "Error in ConfigEnumLocalOccupationsReference: "
                "`event_info` must be an instance of `OccEventSymInfo`."
            )
        if not isinstance(event_supercell_info, casmlocal.OccEventSupercellSymInfo):
            raise ValueError(
                "Error in ConfigEnumLocalOccupationsReference: "
                "`event_supercell_info` must be an instance of "
                "`OccEventSupercellSymInfo`."
            )

        self.event_info = event_info
        """libcasm.local_configuration.OccEventSymInfo: Shared information about the 
        OccEvent."""

        self.event_supercell_info = event_supercell_info
        """libcasm.local_configuration.OccEventSupercellSymInfo: Information about the 
        OccEvent with respect to the supercell."""

        self.initial = initial
        """list[libcasm.local_configuration.LocalConfiguration]: The list of initial 
        LocalConfiguration before perturbations. The configurations are expected
        to be equivalent and primitive, with events in distinct positions. Before
        generating perturbations, these configurations are filled into the 
        supercell without re-orientation and the events are placed at the same 
        coordinates relative to the origin unit cell."""

        # Generate event positions in the supercell for each initial configuration
        _trans_frac = []
        _pos = []
        _event = []
        for i_initial, local_config in enumerate(initial):
            # Validate that the initial local configuration can tile the supercell
            small_supercell = local_config.configuration.supercell
            large_supercell = event_supercell_info.supercell
            is_supercell, T = large_supercell.superlattice.is_superlattice_of(
                small_supercell.superlattice
            )
            if not is_supercell:
                print("~~~")
                print("Initial LocalConfiguration:")
                print(local_config)
                print("Supercell:")
                print(large_supercell)
                raise ValueError(
                    "Error in ConfigEnumLocalOccupations: "
                    "The initial background configuration "
                    "cannot tile the supercell."
                )

            # Determine `pos` in the supercell
            # from `pos` in the initial local config
            f_small = small_supercell.unitcell_index_converter
            f_large = large_supercell.unitcell_index_converter
            trans_frac = f_small.unitcell(local_config.pos[0])
            _trans_frac.append(trans_frac)
            pos = (f_large.linear_unitcell_index(trans_frac), local_config.pos[1])
            _pos.append(pos)

            # Get the event in the supercell at `pos`
            _event.append(event_supercell_info.event(pos))

        self.super_background = [
            casmconfig.copy_configuration(
                motif=x.configuration,
                supercell=event_supercell_info.supercell,
            )
            for x in initial
        ]
        """list[libcasm.configuration.Configuration]: The super configurations 
        formed by filling the `initial` local configurations into the supercell, 
        where the events are placed at the same coordinates relative to the origin 
        unit cell."""

        self.event = _event
        """list[libcasm.occ_events.OccEvent]: The events about which perturbations
        are enumerated, for each initial configuration."""

        self.event_trans_frac = _trans_frac
        """list[np.ndarray]: The fractional translation from the origin unit cell
        to the position of the phenomenal cluster or event, for each initial
        configuration."""

        self.event_pos = _pos
        """list[tuple[int, int]]: The position of the phenomenal cluster or event in
        the supercell, as `(unitcell_index, equivalent_index)`, for each initial
        configuration."""

        self.event_phenomenal_clusters = [x.cluster() for x in self.event]
        """list[libcasm.clusterography.Cluster]: The phenomenal clusters, for each
        event in 
        :py:attr:`ConfigEnumLocalOccupationsReference.event 
        <libcasm.enumerate.ConfigEnumLocalOccupationsReference.event>`."""

        f_large = event_supercell_info.supercell.site_index_converter
        self.event_sites = [
            [f_large.linear_site_index(site) for site in cluster]
            for cluster in self.event_phenomenal_clusters
        ]
        """list[list[int]]: The linear site indices for the event sites, for each
        event in 
        :py:attr:`ConfigEnumLocalOccupationsReference.event 
        <libcasm.enumerate.ConfigEnumLocalOccupationsReference.event>`."""

        # Generate asymmetric unit indices for each local configuration
        _asym_indices = []
        asymmetric_unit_indices_ref = None
        configuration_ref = None
        for local_config in initial:
            small_supercell = local_config.configuration.supercell
            asymmetric_unit_indices = casmconfig.asymmetric_unit_indices(
                configuration=local_config.configuration,
            )
            if asymmetric_unit_indices_ref is None:
                asymmetric_unit_indices_ref = asymmetric_unit_indices
                configuration_ref = local_config.configuration

            else:
                asymmetric_unit_indices = (
                    casmconfig.make_consistent_asymmetric_unit_indices(
                        initial=asymmetric_unit_indices,
                        configuration_init=local_config.configuration,
                        reference=asymmetric_unit_indices_ref,
                        configuration_ref=configuration_ref,
                    )
                )

                if asymmetric_unit_indices is None:
                    raise ValueError(
                        "Error in ConfigEnumLocalOccupationsReference: "
                        "Failed mapping asymmetric unit indices "
                        "between initial configurations."
                    )

            asym_by_small_site_index = [None] * small_supercell.n_sites
            for i_asym_unit, asym_unit in enumerate(asymmetric_unit_indices):
                for site in asym_unit:
                    asym_by_small_site_index[site] = i_asym_unit
            _asym_indices.append(asym_by_small_site_index)

        self.sublattice_indices = event_supercell_info.supercell.sublattice_indices()
        """list[int]: The prim sublattice indices for the sites in the supercell."""

        self.asymmetric_unit_indices = _asym_indices
        """list[list[int]]: The asymmetric unit indices for the sites in each
        `initial` configuration, where `asymmetric_unit_indices[i_initial][l]`
        gives the same integer value for all equivalent sites in the `i_initial`-th
        initial configuration. Asymmetric unit indices are arbitrarily determined
        for the first initial configuration and then made such that they are
        consistent for the other initial configurations."""

        self.occupation_indices = [
            casmconfig.copy_configuration(
                motif=x.configuration,
                supercell=self.event_supercell_info.supercell,
            ).occupation.tolist()
            for x in initial
        ]
        """list[list[int]]: The occupation indices for the sites in each super
        configuration formed by filling the `initial` configurations into the 
        supercell, where `occupation_indices[i_initial][l]` gives the occupation 
        index for the `l`-th site in the `i_initial`-th super configuration."""

        self.prim_local_orbits = prim_local_orbits
        """list[list[list[libcasm.clusterography.Cluster]]]: The 
        local-cluster orbits about each equivalent event in the origin unit cell.
        
        The local orbits for the equivalent events, where
        `prim_local_orbits[i_event][i_local_orbit][i_cluster]` is the 
        `i_cluster`-th cluster in the `i_local_orbit`-th local-cluster orbit about 
        `event_info.event_prim_info.events[i_event]`."""

        # Translate the local orbits to the position in the supercell
        # where the event is located, for each initial configuration
        _event_local_orbits = []
        for i_initial, local_config in enumerate(initial):
            # Get the position of the event in the supercell
            pos = self.event_pos[i_initial]

            # For each local-cluster orbit around the event...
            _curr = []
            for i_local_orbit, local_orbit in enumerate(prim_local_orbits[pos[1]]):
                # Make symmetrically distinct-local cluster sites,
                # as list[set[int]], appropriately translated to the correct `pos`
                # in the supercell
                _curr.append(
                    [
                        cluster + self.event_trans_frac[i_initial]
                        for cluster in local_orbit
                    ]
                )
            _event_local_orbits.append(_curr)

        self.event_local_orbits = _event_local_orbits
        """list[list[list[libcasm.clusterography.Cluster]]]: The local-cluster 
        orbits for each event around which perturbations are made,
        where `event_local_orbits[i_initial][i_local_orbit][i_cluster]` is the
        `i_cluster`-th cluster in the `i_local_orbit`-th local-cluster orbit 
        about `event[i_initial]`."""

    def to_dict(
        self,
        write_prim_basis: bool = False,
    ):
        """Represent the ConfigEnumLocalOccupationsReference as a Python dict.

        Parameters
        ----------
        write_prim_basis : bool, default=False
            If True, write DoF values using the prim basis. Default (False)
            is to write DoF values in the standard basis.

        Returns
        -------
        data : dict
            The ConfigEnumLocalOccupationsReference as a Python dict.
        """
        event_prim_info = self.event_info.event_prim_info
        system = event_prim_info.system
        (
            prototype_event_data,
            equivalents_info_data,
        ) = event_prim_info.to_data()

        xtal_prim = event_prim_info.prim.xtal_prim

        return dict(
            equivalents_info=equivalents_info_data,
            prototype_event=prototype_event_data,
            events=[event.to_dict(system=system) for event in event_prim_info.events],
            supercell=self.event_supercell_info.supercell.to_dict(),
            initial=[
                x.to_dict(write_prim_basis=write_prim_basis) for x in self.initial
            ],
            asymmetric_unit_indices=self.asymmetric_unit_indices,
            sublattice_indices=self.sublattice_indices,
            occupation_indices=self.occupation_indices,
            prim_local_orbits=_make_local_orbits_data(
                local_orbits=self.prim_local_orbits,
                xtal_prim=xtal_prim,
                phenomenal_clusters=event_prim_info.phenomenal_clusters,
            ),
            event_local_orbits=_make_local_orbits_data(
                local_orbits=self.event_local_orbits,
                xtal_prim=xtal_prim,
                phenomenal_clusters=self.event_phenomenal_clusters,
            ),
        )

    def __repr__(self):
        return xtal.pretty_json(self.to_dict())


class ConfigEnumLocalOccupationsResult:
    """Data structure for returning enumerated LocalConfiguration along with
    contextual information

    :class:`Reference` is used to store information about the background
    configurations that is common to all results from the enumeration of local
    perturbations. Individual results are returned as instances of Result.
    """

    def __init__(
        self,
        reference: ConfigEnumLocalOccupationsReference,
    ):
        self.reference = reference
        """libcasm.enumerate.ConfigEnumLocalOccupationsReference: The reference 
        information for the enumeration results."""

        self.local_configuration = None
        """libcasm.local_configuration.LocalConfiguration: The current 
        LocalConfiguration, as constructed from the background configuration and 
        perturbation."""

        self.canonical_local_configuration = None
        """libcasm.local_configuration.LocalConfiguration: The canonical equivalent 
        LocalConfiguration, after applying the event invariant group in the supercell 
        (leaving the event in the same position as the initial, unperturbed 
        configuration)."""

        self.i_initial = None
        """int: The index of the initial configuration in `reference.initial` 
        that was perturbed to create the current result."""

        self.pos = None
        """tuple[int, int]: The position of the phenomenal cluster or event in the 
        supercell, as `(unitcell_index, equivalent_index)`."""

        self.i_local_orbit = None
        """int: The index of the local orbit in `event_local_orbits[i_initial]` 
        that was perturbed."""

        self.sites = None
        """list[int]: The linear site indices for the cluster of sites on which the
        perturbation was applied."""

        self.sublattices = None
        """list[int]: The sublattice indices for the cluster of sites on which the
        perturbation was applied."""

        self.asymmetric_units = None
        """list[int]: The initial configuration asymmetric unit indices for the 
        cluster of sites on which the perturbation was applied, where the indices 
        are defined as in 
        :py:attr:`ConfigEnumLocalOccupationsReference.asymmetric_unit_indices 
        <libcasm.enumerate.ConfigEnumLocalOccupationsReference.asymmetric_unit_indices>`.
        """

        self.initial_occupation = None
        """list[int]: The initial occupation on the sites in `sites`."""

        self.final_occupation = None
        """list[int]: The final occupation on the sites in `sites`."""

    def to_dict(self):
        """Represent the ConfigEnumLocalOccupationsResult as a Python dict.

        Returns
        -------
        data : dict
            The Result as a Python dict.
        """
        return dict(
            local_configuration=self.local_configuration.to_dict(),
            canonical_local_configuration=self.canonical_local_configuration.to_dict(),
            i_initial=self.i_initial,
            pos=[self.pos[0], self.pos[1]],
            i_local_orbit=self.i_local_orbit,
            sites=self.sites,
            sublattices=self.sublattices,
            asymmetric_units=self.asymmetric_units,
            initial_occupation=self.initial_occupation,
            final_occupation=self.final_occupation,
        )

    def __repr__(self):
        return xtal.pretty_json(self.to_dict())


def _print_local_orbits_summary(
    local_orbits: list[list[casmclust.Cluster]],
    xtal_prim: xtal.Prim,
    phenomenal: casmclust.Cluster,
):
    for i_orbit, orbit in enumerate(local_orbits):
        cdist = orbit[0].distances(
            xtal_prim=xtal_prim,
        )
        if len(cdist) > 0:
            cdist_str = "[" + ", ".join([f"{d:.2f}" for d in cdist]) + "]"
        else:
            cdist_str = "None"

        pdist = orbit[0].phenomenal_distances(
            xtal_prim=xtal_prim,
            phenomenal=phenomenal,
        )
        if len(pdist) > 0:
            pdist_str = "[" + ", ".join([f"{d:.2f}" for d in pdist]) + "]"
        else:
            pdist_str = "None"

        print(
            f"- Orbit: {i_orbit}\n"
            f"  - Cluster size: {len(orbit[0])}\n"
            f"  - Orbit size: {len(orbit)}\n"
            f"  - Cluster site distances: {cdist_str}\n"
            f"  - Phenomenal site distances: {pdist_str}"
        )
    print()


class ConfigEnumLocalOccupations:
    """Enumerate occupations in the local environment of an event."""

    def __init__(
        self,
        event_info: casmlocal.OccEventSymInfo,
        supercell_set: Optional[casmconfig.SupercellSet] = None,
        verbose: bool = False,
    ):
        """
        .. rubric:: Constructor

        Parameters
        ----------
        event_info: libcasm.local_configuration.OccEventSymInfo
            The OccEventSymInfo for the event.
        supercell_set: Optional[casmconfig.SupercellSet] = None
            If not None, generated :class:`~casmconfig.Supercell` are constructed by
            adding in the :class:`~casmconfig.SupercellSet`.
        verbose: bool = False
            If True, print additional information about the enumeration process.

        """
        if not isinstance(event_info, casmlocal.OccEventSymInfo):
            raise ValueError(
                "Error in ConfigEnumLocalOccupations: "
                "`event_info` must be an instance of `OccEventSymInfo`."
            )

        self.event_info = event_info
        """libcasm.local_configuration.OccEventSymInfo: The OccEventSymInfo for the 
        event."""

        self.supercell_set = supercell_set
        """Optional[casmconfig.SupercellSet]: If not None, generated 
        :class:`~casmconfig.Supercell` are constructed by adding in the 
        :class:`~casmconfig.SupercellSet`."""

        self.reference = None
        """Optional[libcasm.enumerate.ConfigEnumLocalOccupationsReference]: The 
        reference information for the more recent enumeration results."""

        self.verbose = verbose
        """bool: If True, print additional information about the enumeration process."""

    def _make_prim_local_orbits(
        self,
        cluster_specs: casmclust.ClusterSpecs,
        supercell: casmconfig.Supercell,
    ):
        # Make local orbits for each equivalent event
        event_prim_info = self.event_info.event_prim_info
        prim_local_orbits = event_prim_info.make_local_orbits_from_cluster_specs(
            cluster_specs=cluster_specs,
        )

        if self.verbose:
            _print_local_orbits_summary(
                local_orbits=prim_local_orbits[0],
                xtal_prim=cluster_specs.xtal_prim(),
                phenomenal=event_prim_info.phenomenal_clusters[0],
            )

        # Check if orbits fit in supercell
        required_sites = make_required_sites(
            phenomenal_clusters=event_prim_info.phenomenal_clusters,
            local_orbits=prim_local_orbits,
        )
        if not has_required_sites(
            required_sites=required_sites,
            supercell=supercell,
        ):
            print(
                """
** WARNING: Small supercell / large neighborhood ******
**                                                   **
** Some sites in the neighborhood of the event map   **
** onto each other.                                  **
**                                                   **
** The methods `make_supercells_for_point_defects`,  **
** `find_optimal_point_defect_supercells`, and       **
** `plot_supercells_for_point_defects` may be useful **
** for choosing a supercell that avoids this.        **
**                                                   **
*******************************************************
"""
            )
        return prim_local_orbits

    def _check_supercell_symmetry(
        self,
        initial_configuration: casmconfig.Configuration,
        supercell: casmconfig.Supercell,
        fix: str,
    ):
        """Check that the supercell has higher symmetry than the initial configuration.

        Parameters
        ----------
        initial_configuration : casmconfig.Configuration
            The initial configuration; expected to be primitive for efficiency.
        supercell : casmconfig.Supercell
            The supercell in which to enumerate local perturbations.
        fix : str
            If not "both", then the supercell must not lower the symmetry of the initial
            configuration else a ValueError is raised.
        """
        group = casmconfig.make_invariant_subgroup(
            configuration=initial_configuration,
        )
        background_fg_indices = sorted(
            list(set([g_i.prim_factor_group_index() for g_i in group]))
        )

        scel_fg_indices = supercell.factor_group.head_group_index

        # Check that all background_fg_indices are also in scel_fg_indices
        if not all(i in scel_fg_indices for i in background_fg_indices):
            if fix == "both":
                print(
                    """
** WARNING: Symmetry reduction due to supercell *******
**                                                   **
** The choice of supercell has lower symmetry than   **
** the unperturbed background configuration. This    **
** is usually undesirable because the calculated     **
** energy will depend on the orientation of the      **
** event in the supercell.                           **
**                                                   **
** The methods `make_supercells_for_point_defects`,  **
** `find_optimal_point_defect_supercells`, and       **
** `plot_supercells_for_point_defects` may be useful **
** for choosing a supercell that avoids this.        **
**                                                   **
*******************************************************
"""
                )
            else:
                print(
                    """
** ERROR: Symmetry reduction due to supercell *********
**                                                   **
** The choice of supercell has lower symmetry than   **
** the unperturbed background configuration. This    **
** is usually undesirable because the calculated     **
** energy will depend on the orientation of the      **
** event in the supercell.                           **
**                                                   **
** The methods `make_supercells_for_point_defects`,  **
** `find_optimal_point_defect_supercells`, and       **
** `plot_supercells_for_point_defects` may be useful **
** for choosing a supercell that avoids this.        **
**
** Other options include using `fix="background"` to **
** fix the background configuration orientation, or  **
** using `pos` to specify one or more event          ** 
** positions with `fix="both"` to indicate that the  **
** orientations are specifically requested.          **
**                                                   **
*******************************************************
"""
                )

            if fix != "both":
                raise ValueError(
                    "Error in ConfigEnumLocalOccupations: "
                    "The supercell has lower symmetry than the initial background"
                    'and `fix` != "both".'
                )

    def _yield_results(
        self,
        result,
        i_initial: int,
        orbits: Optional[list[int]] = None,
        neighborhood_from_orbits: Optional[list[int]] = None,
    ):
        ref = result.reference
        small_supercell = ref.initial[i_initial].configuration.supercell
        large_supercell = ref.event_supercell_info.supercell

        # Get the background configuration in the supercell
        super_background = ref.super_background[i_initial]

        # Get the event position, event, and event group in the supercell
        pos = ref.event_pos[i_initial]
        event = ref.event[i_initial]
        event_group_rep = ref.event_supercell_info.event_group_rep(pos)

        # Get supercell site sublattice indices for result output
        sublattice_indices = large_supercell.sublattice_indices()

        # Set result attributes that don't change by cluster
        result.i_initial = i_initial
        result.pos = pos

        # Get distinct local-cluster sites (and i_local_orbit) to perturb
        sites_to_perturb = []
        i_local_orbit_for_sites = []
        allow_subcluster_perturbations = False
        if neighborhood_from_orbits is not None:
            # If `neighborhood_from_orbits` is True, then all local-cluster orbits
            # around the event are combined into a single set of sites to perturb.
            neighborhood_sites = set()
            _i_local_orbit = []
            for i_local_orbit, event_local_orbit in enumerate(
                ref.event_local_orbits[i_initial]
            ):
                if i_local_orbit not in neighborhood_from_orbits:
                    continue
                # Make symmetrically distinct-local cluster sites, as list[set[int]],
                # appropriately translated to the correct `pos` in the supercell
                _i_local_orbit.append(i_local_orbit)
                for cluster in event_local_orbit:
                    for site in cluster:
                        f_large = large_supercell.site_index_converter
                        site_index = f_large.linear_site_index(site)
                        neighborhood_sites.add(site_index)
            sites_to_perturb.append([neighborhood_sites])
            i_local_orbit_for_sites.append(_i_local_orbit)
            allow_subcluster_perturbations = True
        else:
            # For each local-cluster orbit around the event...
            for i_local_orbit, event_local_orbit in enumerate(
                ref.event_local_orbits[i_initial]
            ):
                if orbits is not None:
                    if i_local_orbit not in orbits:
                        continue
                # Make symmetrically distinct-local cluster sites, as list[set[int]],
                # appropriately translated to the correct `pos` in the supercell
                i_local_orbit_for_sites.append(i_local_orbit)
                distinct_local_cluster_sites = make_distinct_local_cluster_sites(
                    configuration=super_background,
                    event=event,
                    event_group=event_group_rep,
                    local_orbits=[event_local_orbit],
                )
                sites_to_perturb.append(distinct_local_cluster_sites)

        # For all distinct local-cluster sites to perturb...
        for i, sites in enumerate(sites_to_perturb):
            result.i_local_orbit = i_local_orbit_for_sites[i]

            # Make distinct local perturbations on the distinct local-cluster sites
            # Note: for a large neighborhood enumeration (i.e., using
            # `neighborhood_from_orbits`) this can take a long time... Could consider
            # implementing to yield each perturbation as it is generated, or allowing
            # for filters.
            perturbations = make_distinct_local_perturbations(
                configuration=super_background,
                event=event,
                event_group=event_group_rep,
                distinct_local_cluster_sites=sites,
                allow_subcluster_perturbations=allow_subcluster_perturbations,
            )

            # For each individual perturbation, store the results in `result` and yield
            for cluster_sites, config, canonical_config in perturbations:
                # Set result local configuration
                result.local_configuration = casmlocal.LocalConfiguration(
                    pos=pos,
                    configuration=config,
                    event_info=ref.event_info,
                )

                # Set result canonical local configuration
                result.canonical_local_configuration = casmlocal.LocalConfiguration(
                    pos=pos,
                    configuration=canonical_config,
                    event_info=ref.event_info,
                )

                # Set site indices for the cluster sites
                result.sites = cluster_sites

                # Set sublattice indices for the cluster sites
                result.sublattices = [
                    sublattice_indices[site_index] for site_index in cluster_sites
                ]

                # Set asymmetric unit indices for the cluster sites
                result.asymmetric_units = []
                f_small = small_supercell.site_index_converter
                f_large = large_supercell.site_index_converter
                for large_site_index in cluster_sites:
                    site = f_large.integral_site_coordinate(large_site_index)
                    small_site_index = f_small.linear_site_index(site)
                    result.asymmetric_units.append(
                        ref.asymmetric_unit_indices[i_initial][small_site_index]
                    )

                # Set initial and final occupations for the cluster sites
                result.initial_occupation = [
                    super_background.occ(s) for s in cluster_sites
                ]
                result.final_occupation = [config.occ(s) for s in cluster_sites]

                # Yield the result
                yield result

    def by_cluster(
        self,
        background: casmconfig.Configuration,
        supercell: casmconfig.Supercell,
        max_length: list[float],
        cutoff_radius: list[float],
        custom_generators: list[casmclust.ClusterOrbitGenerator] = [],
        fix: str = "background",
        pos: Union[tuple[int, int], list[tuple[int, int]], None] = None,
        orbits: Optional[list[int]] = None,
        neighborhood_from_orbits: Optional[list[int]] = None,
    ):
        """Enumerate occupation perturbations on local-clusters in the neighborhood of
        an event in a background configuration

        Parameters
        ----------
        background: casmconfig.Configuration,
            The background configuration on which enumeration takes place. This allows
            fixing other DoF and enumerating all occupations and/or specifying the
            starting occupation which is perturbed.
        supercell: casmconfig.Supercell
            The supercell in which to enumerate local perturbations.

            If the background configuration cannot tile the supercell, a message is
            printed and no results are yielded, but no exception is raised.

            If the supercell has lower symmetry than the prim and this results in
            symmetrically inequivalent ways of generating local configurations in the
            supercell which will have different calculated energies, then a warning
            message is printed. In this case, an exception is raised unless
            `fix="both"`.

            If the supercell is too small or the local neighborhood specified by
            `cluster_specs` is too large, such that periodic images of sites in the
            neighborhood overlap, then a warning message is printed, but no exception is
            raised.

        max_length : list[float]
            The maximum lengths for the local orbits.
        cutoff_radius : list[float]
            The cutoff radii for the local orbits.
        custom_generators : list[libcasm.clusterography.ClusterOrbitGenerator]
            Custom generators for the local orbits, specified about the prototype
            event, `event_info.event_prim_info.prototype_event`.
        fix: str = "background"
            The object to be fixed during the method. Can be "background", "event", or
            "both". If "background", the background is fixed and the event is varied to
            place it in all inequivalent positions before perturbing the local
            environment. If "event", then the event is fixed and the background is
            varied to generate inequivalent positions. If "both", then both the
            background and event are both fixed and perturbations enumerated
            for the event in this specific position only.
        orbits: Optional[list[int]] = None
            If not None, then only the listed local-cluster orbits are perturbed.
        neighborhood_from_orbits: Optional[list[int]] = None
            If not None, then all sites in the listed orbits are combined into a
            single neighborhood of sites to perturb.

        Yields
        ------
        result: libcasm.enumerate.ConfigEnumLocalOccupationsResult
            Contains the local configuration and additional information about the
            enumeration result.
        """
        print("~~~ ConfigEnumLocalOccupations.by_cluster ~~~")
        event_prim_info = self.event_info.event_prim_info
        cluster_specs = casmclust.ClusterSpecs(
            xtal_prim=event_prim_info.prim.xtal_prim,
            generating_group=event_prim_info.prototype_invariant_group,
            max_length=max_length,
            custom_generators=custom_generators,
            phenomenal=event_prim_info.prototype_event.cluster(),
            cutoff_radius=cutoff_radius,
        )
        for result in self.by_cluster_specs(
            background=background,
            supercell=supercell,
            cluster_specs=cluster_specs,
            fix=fix,
            pos=pos,
            orbits=orbits,
            neighborhood_from_orbits=neighborhood_from_orbits,
        ):
            yield result

    def _check_for_distinct_super_configurations(
        self,
        background: casmconfig.Configuration,
        supercell: casmconfig.Supercell,
        fix: str,
    ):
        if self.verbose:
            print("Generating distinct tilings of the background... ", end="")
            sys.stdout.flush()

        super_backgrounds = make_distinct_super_configurations(
            motif=background,
            supercell=supercell,
            fix="motif",
            supercell_set=self.supercell_set,
        )

        if self.verbose:
            print("DONE")
            print("Distinct background configurations: ", len(super_backgrounds))
            print()
            sys.stdout.flush()

        if len(super_backgrounds) == 0:
            print("** WARNING: Background configuration cannot tile the supercell. **")
            print()
            sys.stdout.flush()

        if len(super_backgrounds) > 1:
            print(
                """
** WARNING: Symmetry reduction due to supercell *******
**                                                   **
** The choice of supercell results in symmetrically  **
** distinct ways of tiling the background into the   **
** supercell. This is usually undesirable because    **
** the calculated energy will depend on the          **
** orientation of the event and environment in the   **
** supercell.                                        **
**                                                   **
** The methods `make_supercells_for_point_defects`,  **
** `find_optimal_point_defect_supercells`, and       **
** `plot_supercells_for_point_defects` may be useful **
** for choosing a supercell that avoids this.        **
**                                                   **
*******************************************************
"""
            )

        if len(super_backgrounds) > 1 and fix != "both":
            raise ValueError(
                "Error in ConfigEnumLocalOccupations.by_cluster_specs: "
                "The choice of supercell results in symmetrically distinct ways "
                "of tiling the background into the supercell and `fix` is not 'both'."
            )

        return super_backgrounds

    def _make_event_supercell_info(
        self,
        supercell: casmconfig.Supercell,
    ):
        if self.verbose:
            print("Generating supercell event info... ", end="")
            sys.stdout.flush()

        # Generate event_supercell_info for each super_background supercell
        event_supercell_info = self.event_info.get_event_supercell_info(supercell)

        if self.verbose:
            print("DONE")
            sys.stdout.flush()

        return event_supercell_info

    def by_cluster_specs(
        self,
        background: casmconfig.Configuration,
        supercell: casmconfig.Supercell,
        cluster_specs: casmclust.ClusterSpecs,
        fix: str = "background",
        pos: Union[tuple[int, int], list[tuple[int, int]], None] = None,
        orbits: Optional[list[int]] = None,
        neighborhood_from_orbits: Optional[list[int]] = None,
    ):
        """Enumerate occupation perturbations on local-clusters in the neighborhood of
        an event in a background configuration

        .. rubric:: Method

        1. Generate initial unperturbed local configurations of the
           event in the background using
           :func:`~libcasm.configuration.make_distinct_local_configurations` if `fix`
           is `background` or `event`, or by using the provided `background` and `pos`
           directly if `fix` is "both".
        2. Generate local-cluster orbits about the event in each initial local
           configuration according to the prim factor group symmetry using the provided
           `cluster_specs`.
        3. Generate the distinct local-cluster sites, for each local-cluster orbit,
           for each initial local configuration, using
           :func:`~libcasm.configuration.make_distinct_local_cluster_sites`.
        4. Generate distinct local perturbations on the distinct local-cluster sites
           using :func:`~libcasm.configuration.make_distinct_local_perturbations`.

        .. rubric:: Variations

        - The `orbits` argument allows selecting only some of the local-cluster orbits
          to perturb.
        - The `neighborhood_from_orbits` argument allows combining all sites in
          specified orbits into a single neighborhood of sites to perturb.

        .. rubric:: Results

        Each perturbation is yielded as a
        :class:`~libcasm.enumerate.ConfigEnumLocalOccupationsResult` object which
        contains the local configuration, both in its position as enumerated, and in
        canonical form to allow checking for duplicates, along with additional
        information about the enumeration result such as which orbit was perturbed and
        what the initial and final occupation was on the perturbed sites. This
        information can be used to filter results or just help to understand and check
        the results.

        The :class:`~libcasm.enumerate.ConfigEnumLocalOccupationsResult`
        object is re-used for each yield, but its members are new objects, except for
        the
        :py:attr:`~libcasm.enumerate.ConfigEnumLocalOccupationsResult.reference`
        object which remains constant for all results.

        Resulting :class:`~libcasm.local_configuration.LocalConfiguration` objects can
        be stored in a :class:`~libcasm.local_configuration.LocalConfigurationList`,
        for instance using:

        .. code-block:: Python

            local_config_list = LocalConfigurationList(...)
            config_enum = ConfigEnumLocalOccupations(...)
            for result in config_enum.by_cluster_specs(...):
                if result.local_configuration not in local_config_list:
                    local_config_list.append(result.local_configuration)


        Parameters
        ----------
        background: casmconfig.Configuration,
            The background configuration on which enumeration takes place. This allows
            fixing other DoF and enumerating all occupations and/or specifying the
            starting occupation which is perturbed.
        supercell: casmconfig.Supercell
            The supercell in which to enumerate local perturbations.

            If the background configuration cannot tile the supercell, a message is
            printed and no results are yielded, but no exception is raised.

            If the supercell has lower symmetry than the background configuration there
            will be symmetrically inequivalent ways of generating local configurations
            in the supercell for the same local configuration in the infinite crystal
            and a warning message will be printed. In this case an exception is raised
            unless `fix` is set to `"both"` to indicate that specified orientation in
            the supercell is specifically requested.

            If the supercell is too small or the local neighborhood specified by
            `cluster_specs` is too large, such that periodic images of sites in the
            neighborhood overlap, then a warning message is printed, but no exception is
            raised.

        cluster_specs : libcasm.clusterography.ClusterSpecs
            ClusterSpecs for generating the local-clusters which will be perturbed to
            enumerate distinct environments.
        fix: str = "background"
            The object to be fixed during the method. Can be "background", "event", or
            "both". If "background", the background is fixed and the event is varied to
            place it in all inequivalent positions before perturbing the local
            environment. If "event", then the event is fixed and the background is
            varied to generate inequivalent positions. If "both", then both the
            background and event are both fixed and perturbations enumerated
            for the event in this specific position only.
        pos: Union[tuple[int, int], list[tuple[int,int]], None] = None
            If `fix` is "both", then this specifies the position of one or more events
            in the supercell as `(unitcell_index, equivalent_index)` to use in
            combination with `background` to generate initial configurations.
        orbits: Optional[list[int]] = None
            An optional list of the indices of the local-cluster orbits to perturb. If
            not None, then only the listed local-cluster orbits are perturbed. By
            default, all local-cluster orbits are perturbed.
        neighborhood_from_orbits: Optional[list[int]] = None
            An optional list of the indices of the local-cluster orbits to combine into
            a single neighborhood to enumerate environments in. If not None, then all
            sites in the listed local-cluster orbits are combined into a
            single neighborhood of sites to perturb. By default, perturbations are
            enumerated on individual clusters.


        Yields
        ------
        result: libcasm.enumerate.ConfigEnumLocalOccupationsResult
            Contains the local configuration and additional information about the
            enumeration result. The same `result` object is re-used for each yield,
            but its members are new objects.
        """

        verbose = self.verbose

        if fix not in ["background", "event", "both"]:
            raise ValueError(
                "Error in ConfigEnumLocalOccupations.by_cluster_specs: "
                f"Invalid value for `fix`= {fix}. May be 'background', 'event', or "
                f"'both'."
            )
        if fix == "both" and pos is None:
            raise ValueError(
                "Error in ConfigEnumLocalOccupations.by_cluster_specs: "
                "If `fix` is 'both', then `pos` must be specified."
            )
        if not casmconfig.is_primitive_configuration(background):
            raise ValueError(
                "Error in ConfigEnumLocalOccupations.by_cluster_specs: "
                "The background configuration must be primitive."
            )

        if verbose:
            print("~~~ ConfigEnumLocalOccupations.by_cluster_specs ~~~")
            print()

        if verbose:
            event_prim_info = self.event_info.event_prim_info

            prim_fg_size = len(event_prim_info.prim.factor_group.elements)
            print("Prim:")
            print(f"- factor group size: {prim_fg_size}")
            print()
            sys.stdout.flush()

            is_primitive = casmconfig.is_primitive_configuration(background)
            group = casmconfig.make_invariant_subgroup(
                configuration=background,
            )
            background_fg_indices = sorted(
                list(set([g_i.prim_factor_group_index() for g_i in group]))
            )
            print("Background configuration:")
            print(f"- is primitive: {is_primitive}")
            print(f"- invariant group size: {len(group)}")
            print(f"- invariant group operations: {background_fg_indices}")
            print()
            sys.stdout.flush()

            print("Event orbit:")
            for i_event, event in enumerate(event_prim_info.events):
                group = event_prim_info.invariant_groups[i_event]
                cluster = event.cluster()
                print(f"- equivalent_index: {i_event}")
                print(f"  invariant group operations: {group.head_group_index}")
                print("  sites:")
                for site in cluster:
                    print(f"  - {site}")

            print()
            sys.stdout.flush()

            scel_fg_size = len(supercell.factor_group.elements)
            scel_fg_indices = supercell.factor_group.head_group_index
            print("Supercell:")
            print(f"- n_unitcells: {supercell.n_unitcells}")
            print(f"- factor group size: {scel_fg_size}")
            print(f"- factor group operations: {scel_fg_indices}")
            print()
            sys.stdout.flush()

        super_configurations = self._check_for_distinct_super_configurations(
            background=background,
            supercell=supercell,
            fix=fix,
        )

        if len(super_configurations) == 0:
            return

        event_supercell_info = self._make_event_supercell_info(
            supercell=supercell,
        )

        prim_local_orbits = self._make_prim_local_orbits(
            cluster_specs=cluster_specs,
            supercell=supercell,
        )

        if self.verbose:
            print("Generating distinct local configurations... ", end="")
            sys.stdout.flush()

        if fix == "both":
            if not isinstance(pos, list):
                pos = [pos]
            initial = [
                casmlocal.LocalConfiguration(
                    pos=_pos,
                    configuration=casmconfig.copy_configuration(
                        motif=background,
                        supercell=supercell,
                    ),
                    event_info=self.event_info,
                )
                for _pos in pos
            ]
        else:
            initial = make_distinct_local_configurations(
                background=background,
                event_info=self.event_info,
                fix=fix,
            )

        if self.verbose:
            print("DONE")
            sys.stdout.flush()

            print("Initial local configurations: ", len(initial))
            print()
            for i, x in enumerate(initial):
                print(f"Initial configuration {i}:")
                print(x)

        # If the supercell lowers the symmetry of the initial
        # configuration print a warning. Unless fix="both", also raise an error.
        for x in initial:
            self._check_supercell_symmetry(
                initial_configuration=x.configuration,
                supercell=supercell,
                fix=fix,
            )

        self.reference = ConfigEnumLocalOccupationsReference(
            event_info=self.event_info,
            event_supercell_info=event_supercell_info,
            initial=initial,
            prim_local_orbits=prim_local_orbits,
        )

        result = ConfigEnumLocalOccupationsResult(reference=self.reference)

        for i_initial, x in enumerate(self.reference.initial):
            for result in self._yield_results(
                result=result,
                i_initial=i_initial,
                orbits=orbits,
                neighborhood_from_orbits=neighborhood_from_orbits,
            ):
                yield result
