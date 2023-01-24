import libcasm.clusterography as clust
import libcasm.occ_events._occ_events as _occ_events
import libcasm.sym_info as sym_info
import libcasm.xtal as xtal

def make_occevent_cluster_specs(
    xtal_prim: xtal.Prim,
    phenomenal_occ_event: _occ_events.OccEvent,
    max_length: list[float],
    cutoff_radius: list[float],
    custom_generators: list[clust.ClusterOrbitGenerator] = []) -> clust.ClusterSpecs:
    """Construct ClusterSpecs for local-cluster orbits around an OccEvent

    Parameters
    ----------
    xtal_prim: libcasm.xtal.Prim
        The Prim structure
    phenomenal_occ_event: libcasm.occ_events.OccEvent
        The orbit generating group is the subgroup of the prim factor
        group that leaves `phenomenal_occ_event` invariant.
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
    symgroup_rep = _occ_events.make_occevent_symgroup_rep(
        prim_factor_group.elements(),
        xtal_prim)
    occevent_group = _occ_events.make_occevent_group(
        occ_event=phenomenal_occ_event,
        group=prim_factor_group,
        lattice=xtal_prim.lattice(),
        occevent_symgroup_rep=symgroup_rep)
    return clust.ClusterSpecs(
        xtal_prim=xtal_prim,
        generating_group=occevent_group,
        max_length=max_length,
        phenomenal=phenomenal_occ_event.cluster(),
        cutoff_radius=cutoff_radius,
        custom_generators=custom_generators)
