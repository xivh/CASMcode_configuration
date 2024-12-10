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


def make_first_n_orbits(
    prim: libcasm.configuration.Prim,
    phenomenal: Union[
        libcasm.clusterography.Cluster, libcasm.occ_events.OccEvent, None
    ] = None,
    n_orbits: Optional[list[int]] = None,
    cutoff_radius: Optional[list[float]] = None,
    make_all_possible_orbits: bool = False,
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
    n_orbits: Optional[list[int]] = None
        The number of orbits to include, by number of sites in the cluster. The null
        cluster is always included, and all point cluster orbits are always included
        for periodic orbits. Example: for periodic orbits, `[0, 0, 6, 4]` specifies
        that the null cluster, all point clusters, the first 6 pair cluster orbits, and
        the first 4 triplet cluster orbits should be included.
    cutoff_radius: Optional[list[float]] = None
        For local clusters, the maximum distance of sites from any phenomenal cluster
        site to include in the local environment, by number of sites in the cluster.
        The null cluster value (element 0) is arbitrary.
    make_all_possible_orbits: bool = False
        If True, ignore `n_orbits` and generate all possible local-cluster orbits for
        the specified `cutoff_radius`. If False, raise ValueError if the specified
        `n_orbits` cannot be generated. Default is False. Requires `phenomenal` and
        `cutoff_radius` to be specified.

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

    if n_orbits is None:
        n_orbits = []
    if cutoff_radius is None:
        cutoff_radius = []

    # Set `generating_group`, `phenomenal_cluster`, and minimum `n_orbits`
    if phenomenal is None:
        generating_group = prim.factor_group
        phenomenal_cluster = None
        while len(n_orbits) < 2:
            n_orbits.append(0)
        if make_all_possible_orbits:
            raise ValueError(
                "Error in make_first_n_orbits:"
                "`make_all_possible_orbits` requires `phenomenal` and `cutoff_radius`."
            )
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

    # If `make_all_possible_orbits`, generate all possible orbits
    if make_all_possible_orbits:
        max_cluster_distance = phenomenal_cluster.distances(xtal_prim=xtal_prim)[-1]
        max_length = []
        for i_branch, r in enumerate(cutoff_radius):
            if i_branch < 2:
                max_length.append(0.0)
            else:
                max_length.append(r * 2.0 + max_cluster_distance * 2.0)

        cluster_specs = libcasm.clusterography.ClusterSpecs(
            xtal_prim=xtal_prim,
            phenomenal=phenomenal_cluster,
            include_phenomenal_sites=False,
            generating_group=generating_group,
            max_length=max_length,
            cutoff_radius=cutoff_radius,
        )
        return cluster_specs.make_orbits()

    # Otherwise, prepapre to make first `n` orbits in each branch...

    # Check `n_orbits` and `cuttoff_radius` consistency
    if phenomenal is not None:
        if len(n_orbits) != len(cutoff_radius):
            raise ValueError(
                "Error in make_first_n_orbits:"
                "`n_orbits` and `cutoff_radius` must have the same length "
                "for local-cluster orbits."
            )

    # Prepare to initial `max_length` for each branch...

    # Use the minimum lattice vector length
    # as the initial `max_length` for branches 2, 3, ...
    l_min = min(xtal_prim.lattice().lengths_and_angles()[:3])
    max_length = [0.0, 0.0]
    while len(max_length) < len(n_orbits):
        max_length.append(l_min * 2.0)

    # Method:
    # - Generate orbits
    # - Check if `n_orbits` are generated for each branch
    # - If not complete:
    #   - If no change from previous step, raise ValueError
    #   - Otherwise increase `max_length` for the incomplete branches
    #   - Go to the next iteration
    # - If complete:
    #   - Collect the first `n` orbits from each branch and return them
    total_n_orbits = 0
    while True:
        # Generate orbits
        cluster_specs = libcasm.clusterography.ClusterSpecs(
            xtal_prim=xtal_prim,
            phenomenal=phenomenal_cluster,
            include_phenomenal_sites=False,
            generating_group=generating_group,
            max_length=max_length,
            cutoff_radius=cutoff_radius,
        )
        orbits = cluster_specs.make_orbits()

        # If no change from previous step, raise ValueError
        all_possible_orbits_found = len(orbits) <= total_n_orbits
        if all_possible_orbits_found:
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

        # Organize orbits by branch
        branches = [list() for _ in range(len(cluster_specs.max_length()))]
        for orbit in orbits:
            i_branch = len(orbit[0])
            branches[i_branch].append(orbit)

        # Check if `n_orbits` are generated for each branch
        complete = True
        for i_branch, branch in enumerate(branches):
            if phenomenal is None and i_branch < 2:
                continue
            elif i_branch < 1:
                continue
            if len(branch) < n_orbits[i_branch]:
                complete = False
                max_length[i_branch] += l_min

        # If complete, collect the first `n` orbits from each branch and return them
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

        # Otherwise, continue to the next iteration

    return


def make_first_n_orbits_cluster_specs(
    prim: libcasm.configuration.Prim,
    phenomenal: Union[
        libcasm.clusterography.Cluster, libcasm.occ_events.OccEvent, None
    ] = None,
    n_orbits: Optional[list[int]] = None,
    cutoff_radius: Optional[list[float]] = None,
    make_all_possible_orbits: bool = False,
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
    n_orbits: Optional[list[int]] = None
        The number of orbits to include, by number of sites in the cluster. The null
        cluster is always included, and all point cluster orbits are always included
        for periodic orbits. Example: for periodic orbits, `[0, 0, 6, 4]` specifies
        that the null cluster, all point clusters, the first 6 pair cluster orbits, and
        the first 4 triplet cluster orbits should be included.
    cutoff_radius: Optional[list[float]] = None
        For local clusters, the maximum distance of sites from any phenomenal cluster
        site to include in the local environment, by number of sites in the cluster.
        The null cluster value (element 0) is arbitrary.
    make_all_possible_orbits: bool = False
        If True, ignore `n_orbits` and generate all possible local-cluster orbits for
        the specified `cutoff_radius`. If False, raise ValueError if the specified
        `n_orbits` cannot be generated. Default is False. Requires `phenomenal` and
        `cutoff_radius` to be specified.

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
        make_all_possible_orbits=make_all_possible_orbits,
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
