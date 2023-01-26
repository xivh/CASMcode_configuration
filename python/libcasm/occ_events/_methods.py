import json
import pathlib
import libcasm.clusterography as clust
import libcasm.occ_events._occ_events as _occ_events
import libcasm.sym_info as sym_info
import libcasm.xtal as xtal

def save_occevent(
    root: pathlib.Path,
    name: str,
    occ_event: _occ_events.OccEvent,
    system: _occ_events.OccSystem):
    """Save an OccEvent

    Saves an :class:`~libcasm.occ_events.OccEvent` to:

    - root / "events" / "event.<name>" / "event.json"

    Parameters
    ----------
    root: pathlib.Path
        Path to the parent of the "events" directory.
    name: str
        A name for the event, such as "1NN_A_Va".
    occ_event: ~libcasm.occ_events.OccEvent
        The event to save
    system: ~libcasm.occ_events.OccSystem
        An :class:`~libcasm.occ_events.OccSystem`, used for indexing

    """
    root = pathlib.Path(root)
    path = root / "events" / ("event." + name) / "event.json"
    path.parent.mkdir(parents=True, exist_ok=True)
    with open(path, 'w') as f:
        data = occ_event.to_dict(system)
        f.write(xtal.pretty_json(data))

def load_occevent(
    root: pathlib.Path,
    name: str,
    system: _occ_events.OccSystem) -> _occ_events.OccEvent:
    """Load a saved OccEvent

    Loads an :class:`~libcasm.occ_events.OccEvent` from:

    - root / "events" / "event.<name>" / "event.json"

    Parameters
    ----------
    root: pathlib.Path
        Path to the parent of the "events" directory.
    name: str
        A unique name for the event, such as "1NN_A_Va".
    system: ~libcasm.occ_events.OccSystem
        An :class:`~libcasm.occ_events.OccSystem`, used for indexing

    Returns
    -------
    occ_event: ~libcasm.occ_events.OccEvent
        The saved :class:`~libcasm.occ_events.OccEvent`.
    """
    root = pathlib.Path(root)
    path = root / "events" / ("event." + name) / "event.json"
    with open(path, 'r') as f:
        data = json.load(f)
        occ_event = _occ_events.OccEvent.from_dict(data, system)
    return occ_event

def make_prototype_occevent(
    xtal_prim: xtal.Prim,
    occ_event: _occ_events.OccEvent):
    """Construct the prototype equivalent OccEvent

    Generates an orbit of :class:`~libcasm.occ_events.OccEvent` and
    returns the first element of the orbit.

    Parameters
    ----------
    xtal_prim: libcasm.xtal.Prim
        The Prim structure
    occ_event: ~libcasm.occ_events.OccEvent
        An :class:`~libcasm.occ_events.OccEvent`.

    Returns
    -------
    prototype occ_event: ~libcasm.occ_events.OccEvent
        The prototype equivalent :class:`~libcasm.occ_events.OccEvent`.
    """
    prim_factor_group = sym_info.make_factor_group(xtal_prim)
    occevent_symgroup_rep = _occ_events.make_occevent_symgroup_rep(
        prim_factor_group.elements(), xtal_prim)
    occevent_orbit = _occ_events.make_prim_periodic_orbit(
        occ_event,
        occevent_symgroup_rep)
    return occevent_orbit[0]

def make_occevent_cluster_specs(
    xtal_prim: xtal.Prim,
    phenomenal_occ_event: _occ_events.OccEvent,
    max_length: list[float],
    cutoff_radius: list[float],
    custom_generators: list[clust.ClusterOrbitGenerator] = []
) -> clust.ClusterSpecs:
    """Construct ClusterSpecs for local-cluster orbits around an OccEvent

    Parameters
    ----------
    xtal_prim: libcasm.xtal.Prim
        The Prim structure
    phenomenal_occ_event: ~libcasm.occ_events.OccEvent
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
    custom_generators: list[~libcasm.clusterography.ClusterOrbitGenerator]=[]
          Specifies clusters that should be uses to construct orbits
          regardless of the max_length or cutoff_radius parameters

    Returns
    -------
    cluster_specs: ~libcasm.clusterography.ClusterSpecs
        The resulting ClusterSpecs
    """
    prim_factor_group = sym_info.make_factor_group(xtal_prim)
    symgroup_rep = _occ_events.make_occevent_symgroup_rep(
        prim_factor_group.elements(), xtal_prim)
    occevent_group = _occ_events.make_occevent_group(
        occ_event=phenomenal_occ_event,
        group=prim_factor_group,
        lattice=xtal_prim.lattice(),
        occevent_symgroup_rep=symgroup_rep)
    return clust.ClusterSpecs(xtal_prim=xtal_prim,
                              generating_group=occevent_group,
                              max_length=max_length,
                              phenomenal=phenomenal_occ_event.cluster(),
                              cutoff_radius=cutoff_radius,
                              custom_generators=custom_generators)
