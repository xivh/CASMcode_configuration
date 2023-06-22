import libcasm.configuration
import libcasm.occ_events
import libcasm.clusterography
import libcasm.enumerate._enumerate as _enumerate


def make_occevent_suborbits(
    supercell: libcasm.configuration.Supercell, occ_event: libcasm.occ_events.OccEvent
) -> list[list[libcasm.occ_events.OccEvent]]:
    r"""Make the sub-orbits of OccEvent in a supercell

    This method generates the sub-orbits that occur due to a supercell that has less symmetry than the prim. Each sub-orbit is made up of OccEvent that are equivalent with respect to the supercell factor group. OccEvent in different sub-orbits are equivalent with respect to the prim factor group, but not the supercell factor group.

    Parameters
    ----------
    supercell: ~libcasm.configuration.Supercell
        The :class:`~libcasm.configuration.Supercell`.

    occ_event: ~libcasm.occ_events.OccEvent
        The :class:`~libcasm.occ_events.OccEvent`.

    Returns
    -------
    suborbits: list[list[~libcasm.occ_events.OccEvent]]
        The sub-orbits, where `suborbits[i][j]` is the j-th :class:`~libcasm.occ_events.OccEvent` in the i-th sub-orbit.
    """
    prim = supercell.prim()
    prim_factor_group = prim.factor_group()
    prim_rep = libcasm.occ_events.make_occevent_symgroup_rep(
        prim_factor_group.elements(), prim.xtal_prim()
    )
    orbit = libcasm.occ_events.make_prim_periodic_orbit(occ_event, prim_rep)

    def in_any_suborbit(x, suborbits):
        for suborbit in suborbits:
            if x in suborbit:
                return True
        return False

    suborbits = []
    scel_factor_group = supercell.factor_group()
    scel_rep = libcasm.occ_events.make_occevent_symgroup_rep(
        scel_factor_group.elements(), prim.xtal_prim()
    )
    for x in orbit:
        if in_any_suborbit(x, suborbits):
            continue
        new_suborbit = libcasm.occ_events.make_prim_periodic_orbit(x, scel_rep)
        if new_suborbit not in suborbits:
            suborbits.append(new_suborbit)
    return suborbits


def make_all_distinct_local_perturbations(
    supercell: libcasm.configuration.Supercell,
    occ_event: libcasm.occ_events.OccEvent,
    motif: libcasm.configuration.Configuration,
    local_clusters: list[list[libcasm.clusterography.Cluster]],
) -> list[libcasm.configuration.Configuration]:
    r"""
    Construct distinct local perturbations of a configuration

    This method constructs symmetrically distinct perturbations of a
    "motif" configuration, in a given supercell, on specified
    local-clusters around an :class:`~libcasm.occ_events.OccEvent`. The
    method works as follows:

    - Each local cluster is used to generated local-cluster orbits around the :class:`~libcasm.occ_events.OccEvent`, using the symmetry of the :class:`~libcasm.xtal.Prim` and :class:`~libcasm.occ_events.OccEvent`, without regard to the background :class:`~libcasm.configuration.Configuration` or :class:`~libcasm.configuration.Supercell`. Duplicate local-cluster orbits are discarded.
    - The motif :class:`~libcasm.configuration.Configuration` is tiled into the :class:`~libcasm.configuration.Supercell` to generate all possible background :class:`~libcasm.configuration.Configuration`, using :class:`~libcasm.configuration.make_all_super_configurations`.
    - The subgroup of the supercell factor group that leaves the :class:`~libcasm.occ_events.OccEvent` invariant is used to find symmetrically distinct background :class:`~libcasm.configuration.Configuration`.
    - For each distinct background :class:`~libcasm.configuration.Configuration`, the local-cluster orbits are used to find the distinct local-clusters, now taking into account the background :class:`~libcasm.configuration.Configuration` and :class:`~libcasm.configuration.Supercell`.
    - For each distinct local-cluster in a distinct background configuration, all possible changes in occupation are performed to generate local perturbation configurations.
    - Each local perturbation configuration is put into a canonical form, using the subgroup of the supercell factor group that leaves the :class:`~libcasm.occ_events.OccEvent` invariant, in order to identify the distinct local perturbation :class:`~libcasm.configuration.Configuration`.

    Parameters
    ----------
    supercell : ~libcasm.configuration.Supercell
      The supercell in which local environment configurations will
      be generated.

    occ_event: ~libcasm.occ_events.OccEvent
      The occupation event.

    motif: ~libcasm.configuration.Configuration
      The motif configuration is tiled into the supercell to generate
      background configurations. The occ_event is kept in the same position
      while prim factor group symmetry operations are applied to the motif
      configuration to generate all symmetrically distinct combinations of
      the background configuration and occ_event before generating
      perturbations. Only perfect tilings into the supercell are kept.

    local_clusters: list[list[~libcasm.clusterography.Cluster]]
      Local clusters, on which the occupation variables will be enumerated
      in order to generate local perturbation configurations. The initial
      cluster in each local-cluster orbit generated using
      :func:`~libcasm.occ_events.make_occevent_cluster_specs` is an
      appropriate input for this parameter. Each local cluster is first used
      to generated local-cluster orbits without regard to the background
      congifuration or supercell. Then, the local-clusters that are distinct
      taking the background configuration and supercell into account are
      perturbed with each possible occupation.

    Returns
    -------
    configurations : list[~libcasm.configuration.Configuration]
      The symmetrically distinct perturbations around the
      event on the specified local-clusters, in all distinct
      combinations of the event and motif configuration in
      the chosen supercell.

    """
    return _enumerate.make_all_distinct_local_perturbations(
        supercell, occ_event, motif, local_clusters
    )
