from typing import Optional, Union

import libcasm.clusterography
import libcasm.configuration
import libcasm.enumerate._enumerate as _enumerate
import libcasm.occ_events
import libcasm.xtal


def make_all_distinct_periodic_perturbations(
    supercell: libcasm.configuration.Supercell,
    motif: libcasm.configuration.Configuration,
    clusters: list[libcasm.clusterography.Cluster],
) -> list[libcasm.configuration.Configuration]:
    r"""
    Construct distinct local perturbations of a configuration

    This method constructs symmetrically distinct perturbations of a "motif"
    configuration, in a given supercell, on specified clusters. The method works as
    follows:

    - Each cluster is used to generate cluster orbits using the symmetry of the
      :class:`~libcasm.xtal.Prim`, without regard to the motif
      :class:`~libcasm.configuration.Configuration` or
      :class:`~libcasm.configuration.Supercell`. Duplicate cluster orbits are
      discarded.
    - The motif :class:`~libcasm.configuration.Configuration` is tiled into the
      :class:`~libcasm.configuration.Supercell` to generate all possible background
      :class:`~libcasm.configuration.Configuration`, using
      :class:`~libcasm.configuration.make_all_super_configurations`.
    - For each distinct background :class:`~libcasm.configuration.Configuration`, the
      cluster orbits are used to find the distinct clusters, now taking
      into account the background :class:`~libcasm.configuration.Configuration` and
      :class:`~libcasm.configuration.Supercell`.
    - For each distinct cluster in a distinct background configuration, all
      possible changes in occupation are performed to generate perturbation
      configurations.
    - Each perturbation configuration is put into a canonical form, using the
      supercell factor group in order to identify the distinct perturbation
      :class:`~libcasm.configuration.Configuration`.

    Parameters
    ----------
    supercell : ~libcasm.configuration.Supercell
        The supercell in which perturbation configurations will be generated.

    motif: ~libcasm.configuration.Configuration
        The motif configuration is tiled into the supercell to generate
        background configurations that are perturbed. Only perfect tilings into the
        supercell are kept.

    clusters: list[~libcasm.clusterography.Cluster]
        Clusters, on which the occupation variables will be enumerated in order
        to generate perturbation configurations. Each cluster is first used
        to generate orbits without regard to the motif or supercell. Then,
        the clusters that are distinct taking the motif and supercell into
        account are perturbed with each possible occupation.

    Returns
    -------
    configurations : list[~libcasm.configuration.Configuration]
        The symmetrically distinct perturbations of the motif configuration,
        in the supercell, on the specified clusters.
    """
    return _enumerate.make_all_distinct_periodic_perturbations(
        supercell, motif, clusters
    )


def make_occevent_suborbits(
    supercell: libcasm.configuration.Supercell, occ_event: libcasm.occ_events.OccEvent
) -> list[list[libcasm.occ_events.OccEvent]]:
    r"""Make the sub-orbits of OccEvent in a supercell

    This method generates the sub-orbits that occur due to a supercell that has less
    symmetry than the prim. Each sub-orbit is made up of OccEvent that are equivalent
    with respect to the supercell factor group. OccEvent in different sub-orbits are
    equivalent with respect to the prim factor group, but not the supercell factor
    group.

    Parameters
    ----------
    supercell: ~libcasm.configuration.Supercell
        The :class:`~libcasm.configuration.Supercell`.

    occ_event: ~libcasm.occ_events.OccEvent
        The :class:`~libcasm.occ_events.OccEvent`.

    Returns
    -------
    suborbits: list[list[~libcasm.occ_events.OccEvent]]
        The sub-orbits, where `suborbits[i][j]` is the j-th
        :class:`~libcasm.occ_events.OccEvent` in the i-th sub-orbit.
    """
    prim = supercell.prim
    prim_factor_group = prim.factor_group
    prim_rep = libcasm.occ_events.make_occevent_symgroup_rep(
        prim_factor_group.elements, prim.xtal_prim
    )
    orbit = libcasm.occ_events.make_prim_periodic_orbit(occ_event, prim_rep)

    def in_any_suborbit(x, suborbits):
        for suborbit in suborbits:
            if x in suborbit:
                return True
        return False

    suborbits = []
    scel_factor_group = supercell.factor_group
    scel_rep = libcasm.occ_events.make_occevent_symgroup_rep(
        scel_factor_group.elements, prim.xtal_prim
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
    local_clusters: list[libcasm.clusterography.Cluster],
) -> list[libcasm.configuration.Configuration]:
    r"""
    Construct distinct local perturbations of a configuration

    This method constructs symmetrically distinct perturbations of a
    "motif" configuration, in a given supercell, on specified
    local-clusters around an :class:`~libcasm.occ_events.OccEvent`. The
    method works as follows:

    - Each local cluster is used to generate local-cluster orbits around the
      :class:`~libcasm.occ_events.OccEvent`, using the symmetry of the
      :class:`~libcasm.xtal.Prim` and :class:`~libcasm.occ_events.OccEvent`, without
      regard to the background :class:`~libcasm.configuration.Configuration` or
      :class:`~libcasm.configuration.Supercell`. Duplicate local-cluster orbits are
      discarded.
    - The motif :class:`~libcasm.configuration.Configuration` is tiled into the
      :class:`~libcasm.configuration.Supercell` to generate all possible background
      :class:`~libcasm.configuration.Configuration`, using
      :class:`~libcasm.configuration.make_all_super_configurations`.
    - The subgroup of the supercell factor group that leaves the
      :class:`~libcasm.occ_events.OccEvent` invariant is used to find symmetrically
      distinct background :class:`~libcasm.configuration.Configuration`.
    - For each distinct background :class:`~libcasm.configuration.Configuration`, the
      local-cluster orbits are used to find the distinct local-clusters, now taking
      into account the background :class:`~libcasm.configuration.Configuration` and
      :class:`~libcasm.configuration.Supercell`.
    - For each distinct local-cluster in a distinct background configuration, all
      possible changes in occupation are performed to generate local perturbation
      configurations.
    - Each local perturbation configuration is put into a canonical form, using the
      subgroup of the supercell factor group that leaves the
      :class:`~libcasm.occ_events.OccEvent` invariant, in order to identify the
      distinct local perturbation :class:`~libcasm.configuration.Configuration`.

    Parameters
    ----------
    supercell : ~libcasm.configuration.Supercell
        The supercell in which local environment configurations will be generated.

    occ_event: ~libcasm.occ_events.OccEvent
        The occupation event.

    motif: ~libcasm.configuration.Configuration
        The motif configuration is tiled into the supercell to generate
        background configurations. The occ_event is kept in the same position
        while prim factor group symmetry operations are applied to the motif
        configuration to generate all symmetrically distinct combinations of
        the background configuration and occ_event before generating
        perturbations. Only perfect tilings into the supercell are kept.

    local_clusters: list[~libcasm.clusterography.Cluster]
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


def make_symop_inverse(
    op: libcasm.xtal.SymOp,
) -> libcasm.xtal.SymOp:
    """Make the inverse SymOp

    Parameters
    ----------
    op: libcasm.xtal.SymOp
        A SymOp

    Returns
    -------
    op_inv: libcasm.xtal.SymOp
        The inverse SymOp
    """
    return libcasm.xtal.SymOp(
        matrix=op.matrix().T,
        translation=-(op.matrix().T @ op.translation()),
        time_reversal=op.time_reversal(),
    )


def make_equivalents_generators(
    phenomenal_prototype: libcasm.clusterography.Cluster,
    generating_group: libcasm.sym_info.SymGroup,
    prim: libcasm.configuration.Prim,
) -> tuple[
    list[libcasm.xtal.SymOp], list[int], list[libcasm.xtal.IntegralSiteCoordinateRep]
]:
    """Make symmetry operations that generate all the equivalent local orbits in \
    the primitive cell.

    Notes
    -----

    - If the `cluster_specs.generating_group()` is a subgroup of the phenomenal
      cluster group, then there will be >1 distinct sets of local orbits around the
      phenomenal cluster.
    - This method finds the cluster group operations that generate the distinct set of
      local orbits, then it uses the phenomenal cluster orbit equivalence map to
      find operations that transform the local orbits sets to the equivalent
      phenomenal clusters in the primitive cell.

    Parameters
    ----------
    phenomenal_prototype: casmclust.Cluster
        The prototype phenomenal cluster. The prototype phenomenal cluster must be
        chosen from one of the equivalents that is generated by
        :func:`~libcasm.clusterography.make_periodic_orbit` using the prim factor
        group.
    generating_group: sym_info.SymGroup
        The local orbits generating group.
    prim: libcasm.configuration.Prim
        The prim, with symmetry info

    Returns
    -------
    ops: list[libcasm.xtal.SymOp]
        Symmetry operations that generate the equivalent local orbits in
        the primitive unit cell.

    indices: list[int]
        Indices of prim factor group operations corresponding to `ops`.

    site_reps: list[libcasm.xtal.IntegralSiteCoordinateRep]
        Symmetry group representation for transforming IntegralSiteCoordinate
        and Cluster corresponding to `ops`.
    """
    import libcasm.clusterography as casmclust

    # collect needed sym groups
    factor_grp = prim.factor_group
    generating_grp = generating_group
    site_rep = prim.integral_site_coordinate_symgroup_rep

    phenomenal_cluster_grp = casmclust.make_cluster_group(
        cluster=phenomenal_prototype,
        group=factor_grp,
        lattice=prim.xtal_prim.lattice(),
        integral_site_coordinate_symgroup_rep=site_rep,
    )

    # find indices of factor group operations
    # that generate equivalent (but rotated) local basis sets
    # about the phenomenal cluster.
    i_factor_grp_equiv_on_phenomenal = set()
    for i_cluster_grp in range(len(phenomenal_cluster_grp.elements)):
        i_factor_grp = phenomenal_cluster_grp.head_group_index[i_cluster_grp]

        i_factor_grp_min = i_factor_grp
        for j_generating_grp in range(len(generating_grp.elements)):
            j_factor_grp = generating_grp.head_group_index[j_generating_grp]

            k_product = factor_grp.mult(j_factor_grp, i_factor_grp)
            if k_product < i_factor_grp_min:
                i_factor_grp_min = k_product

        i_factor_grp_equiv_on_phenomenal.add(i_factor_grp_min)

    # generate the phenomenal cluster orbit and equivalence map
    phenomenal_orbit = casmclust.make_periodic_orbit(
        orbit_element=phenomenal_prototype,
        integral_site_coordinate_symgroup_rep=site_rep,
    )
    equivalence_map = casmclust.make_periodic_equivalence_map(
        orbit=phenomenal_orbit,
        symgroup=factor_grp,
        lattice=prim.xtal_prim.lattice(),
        integral_site_coordinate_symgroup_rep=site_rep,
    )
    equivalence_map_indices = casmclust.make_periodic_equivalence_map_indices(
        orbit=phenomenal_orbit,
        integral_site_coordinate_symgroup_rep=site_rep,
    )

    # find op and factor group index that transform phenomenal to orbit prototype
    to_prototype_op = None
    i_to_prototype = None
    sorted_phenomenal = phenomenal_prototype.sorted()
    for i_clust, clust in enumerate(phenomenal_orbit):
        sorted_clust = clust.sorted()
        if sorted_clust == sorted_phenomenal:
            to_prototype_op = make_symop_inverse(equivalence_map[i_clust][0])
            i_to_prototype = factor_grp.inv(equivalence_map_indices[i_clust][0])
            break
    if to_prototype_op is None:
        raise Exception(
            "Error in make_equivalents_generators: Failed to find to_prototype_op."
            "The phenomenal cluster must be chosen from one of the"
            "equivalents generated by `libcasm.clusterography.make_periodic_orbit`."
        )

    # equiv bset = to_equiv_op * (to_prototype_op * (cluster_grp_op * phenomenal bset))
    #                                              ^ distinct bset on phenomenal
    #                             ^ bset on prototype
    #              ^ bset on other phenomenal
    generating_indices = []
    for i_factor_grp in i_factor_grp_equiv_on_phenomenal:
        for i_equiv in range(len(equivalence_map)):
            i_to_equiv = equivalence_map_indices[i_equiv][0]
            i_generating_op = factor_grp.mult(
                i_to_equiv, factor_grp.mult(i_to_prototype, i_factor_grp)
            )

            generating_indices.append(i_generating_op)
    generating_indices.sort()

    generating_ops = []
    for i in generating_indices:
        generating_ops.append(factor_grp.elements[i])

    generating_site_reps = casmclust.make_integral_site_coordinate_symgroup_rep(
        group_elements=generating_ops,
        xtal_prim=prim.xtal_prim,
    )

    return (generating_ops, generating_indices, generating_site_reps)


def make_occevent_equivalents_generators(
    prototype_event: libcasm.occ_events.OccEvent,
    prim: libcasm.configuration.Prim,
) -> tuple[
    list[libcasm.xtal.SymOp], list[int], list[libcasm.xtal.IntegralSiteCoordinateRep]
]:
    """Make symmetry operations that generate all the equivalent OccEvent in \
    the primitive cell.

    Parameters
    ----------
    prototype_event: libcasm.occ_events.OccEvent
        The prototype OccEvent. The underlying cluster `prototype_event.cluster` must be
        chosen from one of the equivalents that is generated by
        :func:`~libcasm.clusterography.make_periodic_orbit` using the prim factor
        group.
    prim: libcasm.configuration.Prim
        The prim, with symmetry info

    Returns
    -------
    ops: list[libcasm.xtal.SymOp]
        Symmetry operations that generate the equivalent local orbits in
        the primitive unit cell.

    indices: list[int]
        Indices of prim factor group operations corresponding to `ops`.

    site_reps: list[libcasm.xtal.IntegralSiteCoordinateRep]
        Symmetry group representation for transforming IntegralSiteCoordinate
        and Cluster corresponding to `ops`.
    """
    occevent_symgroup_rep = libcasm.occ_events.make_occevent_symgroup_rep(
        prim.factor_group.elements, prim.xtal_prim
    )
    generating_group = libcasm.occ_events.make_occevent_group(
        occ_event=prototype_event,
        group=prim.factor_group,
        lattice=prim.xtal_prim.lattice(),
        occevent_symgroup_rep=occevent_symgroup_rep,
    )
    return make_equivalents_generators(
        phenomenal_prototype=prototype_event.cluster(),
        generating_group=generating_group,
        prim=prim,
    )


def make_first_n_orbits(
    prim: libcasm.configuration.Prim,
    phenomenal: Union[
        libcasm.clusterography.Cluster, libcasm.occ_events.OccEvent, None
    ] = None,
    n_orbits: Optional[list[int]] = [],
    cutoff_radius: Optional[list[float]] = [],
):
    """Construct the first `n` local cluster orbits of each cluster size.

    Parameters
    ----------
    prim: libcasm.configuration.Prim
        The Prim.
    phenomenal: Union[libcasm.clusterography.Cluster, libcasm.occ_events.OccEvent, \
    None] = None
        If provided, generate local-cluster orbits using the invariant group of the
        phenomenal cluster or event. By default, periodic cluster orbits are
        generated.
    n_orbits: Optional[list[int]] = []
        The number of orbits to include, by number of sites in the cluster. The null
        cluster is always included, and all point cluster orbits are always included
        for periodic orbits. Example: for periodic orbits, `[0, 0, 6, 4]` specifies
        that the null cluster, all point clusters, the first 6 pair cluster orbits, and
        the first 4 triplet cluster orbits should be included.
    cutoff_radius: Optional[list[float]] = []
        For local clusters, the maximum distance of sites from any phenomenal cluster
        site to include in the local environment, by number of sites in the cluster.
        The null cluster value (element 0) is arbitrary.

    Returns
    -------
    orbits: list[list[libcasm.clusterography.Cluster]]
        The resulting orbits, including the null cluster orbit and all point cluster
        orbits.
    """

    from libcasm.clusterography import (
        Cluster,
        make_cluster_group,
        make_integral_site_coordinate_symgroup_rep,
    )
    from libcasm.occ_events import (
        OccEvent,
        make_occevent_group,
        make_occevent_symgroup_rep,
    )

    xtal_prim = prim.xtal_prim

    if phenomenal is None:
        generating_group = prim.factor_group
        phenomenal_cluster = None
        while len(n_orbits) < 2:
            n_orbits.append(0)
    elif isinstance(phenomenal, Cluster):
        symgroup_rep = make_integral_site_coordinate_symgroup_rep(
            prim.factor_group.elements, prim.xtal_prim
        )
        generating_group = make_cluster_group(
            cluster=phenomenal,
            group=prim.factor_group,
            lattice=prim.xtal_prim.lattice(),
            integral_site_coordinate_symgroup_rep=symgroup_rep,
        )
        phenomenal_cluster = phenomenal
        while len(n_orbits) < 1:
            n_orbits.append(0)
    elif isinstance(phenomenal, OccEvent):
        symgroup_rep = make_occevent_symgroup_rep(
            prim.factor_group.elements, prim.xtal_prim
        )
        generating_group = make_occevent_group(
            occ_event=phenomenal,
            group=prim.factor_group,
            lattice=prim.xtal_prim.lattice(),
            occevent_symgroup_rep=symgroup_rep,
        )
        phenomenal_cluster = phenomenal.cluster()
        while len(n_orbits) < 1:
            n_orbits.append(0)
    else:
        raise ValueError(
            "Error in make_first_n_orbits:"
            "`phenomenal` must be a Cluster, OccEvent, or None"
        )

    l_1 = min(xtal_prim.lattice().lengths_and_angles()[:3])
    l_current = l_1

    total_n_orbits = 0
    while True:
        max_length = [0.0, 0.0]
        while len(max_length) < len(n_orbits):
            max_length.append(l_current)

        cluster_specs = libcasm.clusterography.ClusterSpecs(
            xtal_prim=xtal_prim,
            phenomenal=phenomenal_cluster,
            include_phenomenal_sites=False,
            generating_group=generating_group,
            max_length=max_length,
            cutoff_radius=cutoff_radius,
        )
        orbits = cluster_specs.make_orbits()
        if len(orbits) <= total_n_orbits:
            if phenomenal is None:
                raise ValueError(
                    "Error in make_first_n_orbits:"
                    "Failed to generate new orbits. Unknown error."
                )
            else:
                raise ValueError(
                    "Error in make_first_n_orbits:"
                    "Failed to generate new orbits. Try increasing the cutoff radius."
                )
        total_n_orbits = len(orbits)

        branches = [list() for _ in range(len(cluster_specs.max_length()))]
        for orbit in orbits:
            i_branch = len(orbit[0])
            branches[i_branch].append(orbit)

        complete = True
        for i_branch, branch in enumerate(branches):
            if phenomenal is None and i_branch < 2:
                continue
            elif i_branch < 1:
                continue
            if len(branch) < n_orbits[i_branch]:
                complete = False
                break

        if complete:
            orbits = []
            for i_branch, branch in enumerate(branches):
                if i_branch == 0:
                    orbits += branch
                elif phenomenal is None and i_branch == 1:
                    orbits += branch
                else:
                    orbits += branch[: n_orbits[i_branch]]

            return orbits
        else:
            l_current += l_1

    return


def make_first_n_orbits_cluster_specs(
    prim: libcasm.configuration.Prim,
    phenomenal: Union[
        libcasm.clusterography.Cluster, libcasm.occ_events.OccEvent, None
    ] = None,
    n_orbits: Optional[list[int]] = [],
    cutoff_radius: Optional[list[float]] = [],
):
    """Construct ClusterSpecs for generating the first `n` local cluster orbits of each
    cluster size using `custom_generators` for all orbits.

    Parameters
    ----------
    prim: libcasm.configuration.Prim
        The Prim.
    phenomenal: Union[libcasm.clusterography.Cluster, libcasm.occ_events.OccEvent, \
    None] = None
        If provided, generate local-cluster orbits using the invariant group of the
        phenomenal cluster or event. By default, periodic cluster orbits are
        generated.
    n_orbits: Optional[list[int]] = []
        The number of orbits to include, by number of sites in the cluster. The null
        cluster is always included, and all point cluster orbits are always included
        for periodic orbits. Example: for periodic orbits, `[0, 0, 6, 4]` specifies
        that the null cluster, all point clusters, the first 6 pair cluster orbits, and
        the first 4 triplet cluster orbits should be included.
    cutoff_radius: Optional[list[float]] = []
        For local clusters, the maximum distance of sites from any phenomenal cluster
        site to include in the local environment, by number of sites in the cluster.
        The null cluster value (element 0) is arbitrary.

    Returns
    -------
    cluster_specs: libcasm.clusterography.ClusterSpecs
        ClusterSpecs for generating the specified orbits using `custom_generators` for
        all orbits.
    """

    orbits = make_first_n_orbits(
        prim=prim,
        phenomenal=phenomenal,
        n_orbits=n_orbits,
        cutoff_radius=cutoff_radius,
    )
    custom_generators = [
        libcasm.clusterography.ClusterOrbitGenerator(
            prototype=orbit[0],
            include_subclusters=False,
        )
        for orbit in orbits
    ]

    from libcasm.clusterography import (
        Cluster,
        make_cluster_group,
        make_integral_site_coordinate_symgroup_rep,
    )
    from libcasm.occ_events import (
        OccEvent,
        make_occevent_group,
        make_occevent_symgroup_rep,
    )

    if phenomenal is None:
        generating_group = prim.factor_group
        phenomenal_cluster = None
    elif isinstance(phenomenal, Cluster):
        symgroup_rep = make_integral_site_coordinate_symgroup_rep(
            prim.factor_group.elements, prim.xtal_prim
        )
        generating_group = make_cluster_group(
            cluster=phenomenal,
            group=prim.factor_group,
            lattice=prim.xtal_prim.lattice(),
            integral_site_coordinate_symgroup_rep=symgroup_rep,
        )
        phenomenal_cluster = phenomenal
    elif isinstance(phenomenal, OccEvent):
        symgroup_rep = make_occevent_symgroup_rep(
            prim.factor_group.elements, prim.xtal_prim
        )
        generating_group = make_occevent_group(
            occ_event=phenomenal,
            group=prim.factor_group,
            lattice=prim.xtal_prim.lattice(),
            occevent_symgroup_rep=symgroup_rep,
        )
        phenomenal_cluster = phenomenal.cluster()
    else:
        raise ValueError(
            "Error in make_first_n_orbits_cluster_specs:"
            "`phenomenal` must be a Cluster, OccEvent, or None"
        )

    return libcasm.clusterography.ClusterSpecs(
        xtal_prim=prim.xtal_prim,
        phenomenal=phenomenal_cluster,
        generating_group=generating_group,
        custom_generators=custom_generators,
        include_phenomenal_sites=False,
    )
