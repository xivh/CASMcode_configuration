import libcasm.clusterography._clusterography as _clust
import libcasm.sym_info as sym_info
import libcasm.xtal as xtal


def make_periodic_cluster_specs(
    xtal_prim: xtal.Prim,
    max_length: list[float],
    custom_generators: list[_clust.ClusterOrbitGenerator] = [],
) -> _clust.ClusterSpecs:
    """Construct ClusterSpecs for periodic orbits using the prim factor group

    Parameters
    ----------
    xtal_prim: libcasm.xtal.Prim
        The Prim structure
    max_length: list[float]
        The maximum site-to-site distance to allow in clusters, by number
        of sites in the cluster. Example: `[0.0, 0.0, 5.0, 4.0]` specifies
        that pair clusters up to distance 5.0 and triplet clusters up to
        distance 4.0 should be included. The null cluster and point
        cluster values (elements 0 and 1) are arbitrary.
    custom_generators: list[libcasm.clusterography.ClusterOrbitGenerator]=[]
          Specifies clusters that should be uses to construct orbits
          regardless of the max_length or cutoff_radius parameters

    Returns
    -------
    cluster_specs: libcasm.clusterography.ClusterSpecs
        The resulting ClusterSpecs
    """
    prim_factor_group = sym_info.make_factor_group(xtal_prim)
    return _clust.ClusterSpecs(
        xtal_prim=xtal_prim,
        generating_group=prim_factor_group,
        max_length=max_length,
        custom_generators=custom_generators,
    )


def make_local_cluster_specs(
    xtal_prim: xtal.Prim,
    phenomenal_cluster: _clust.Cluster,
    max_length: list[float],
    cutoff_radius: list[float],
    custom_generators: list[_clust.ClusterOrbitGenerator] = [],
) -> _clust.ClusterSpecs:
    """Construct ClusterSpecs for local-cluster orbits around a cluster

    Parameters
    ----------
    xtal_prim: libcasm.xtal.Prim
        The Prim structure
    phenomenal_cluster: libcasm.clusterography.Cluster
        The orbit generating group is the subgroup of the prim factor
        group that leaves `phenomenal_cluster` invariant.
    max_length: list[float]
        The maximum site-to-site distance to allow in clusters, by number
        of sites in the cluster. Example: `[0.0, 0.0, 5.0, 4.0]` specifies
        that pair clusters up to distance 5.0 and triplet clusters up to
        distance 4.0 should be included. The null cluster and point
        cluster values (elements 0 and 1) are arbitrary.
    cutoff_radius: list[float]
        For local clusters, the maximum distance of sites from any
        phenomenal cluster site to include in the local environment, by
        number of sites in the cluster. The null cluster value
        (element 0) is arbitrary.
    custom_generators: list[libcasm.clusterography.ClusterOrbitGenerator]=[]
          Specifies clusters that should be uses to construct orbits
          regardless of the max_length or cutoff_radius parameters

    Returns
    -------
    cluster_specs: libcasm.clusterography.ClusterSpecs
        The resulting ClusterSpecs
    """
    prim_factor_group = sym_info.make_factor_group(xtal_prim)
    symgroup_rep = _clust.make_integral_site_coordinate_symgroup_rep(
        prim_factor_group.elements, xtal_prim
    )
    occevent_group = _clust.make_cluster_group(
        cluster=phenomenal_cluster,
        group=prim_factor_group,
        lattice=xtal_prim.lattice(),
        integral_site_coordinate_symgroup_rep=symgroup_rep,
    )
    return _clust.ClusterSpecs(
        xtal_prim=xtal_prim,
        generating_group=occevent_group,
        max_length=max_length,
        phenomenal=phenomenal_cluster,
        cutoff_radius=cutoff_radius,
        custom_generators=custom_generators,
    )
