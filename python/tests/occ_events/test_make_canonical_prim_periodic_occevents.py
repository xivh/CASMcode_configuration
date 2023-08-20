import math
import sys

import libcasm.clusterography as clust
import libcasm.occ_events as occ_events
import libcasm.sym_info as sym_info
import libcasm.xtal.prims as xtal_prims


def test_make_canonical_prim_periodic_occevents():
    r = 1.0  # ideal atom radius
    a = math.sqrt(((4 * r) ** 2) / 2.0)  # conventional FCC lattice parameter
    tol = 1e-5
    xtal_prim = xtal_prims.FCC(r=r, occ_dof=["A", "B", "Va"])

    # The OccSystem provides index conversions
    system = occ_events.OccSystem(xtal_prim)

    # The maximum site-to-site distance to allow in clusters,
    # by number of sites in the cluster. The null cluster and
    # point cluster values (elements 0 and 1) are arbitrary
    # for periodic clusters.
    max_length = [
        0.0,  # null-cluster orbit
        0.0,  # point-cluster orbits
        a + tol,  # pair-cluster orbits, including 2NN sites
        a + tol,  # triplet-cluster orbits, including 2NN sites
    ]

    # Custom generators is a list[clust.ClusterOrbitGenerator]
    # that allows specifying custom clusters to include,
    # independent of the max_length cutoff,
    # and optionally also including subclusters
    custom_generators = []

    # Construct ClusterSpecs, with generating group equal to
    # the invariant group of prototype_occ_event
    cluster_specs = clust.ClusterSpecs(
        xtal_prim=xtal_prim,
        generating_group=sym_info.make_factor_group(xtal_prim),
        max_length=max_length,
        custom_generators=custom_generators,
    )

    # null, point, 2 pair, 2 triplet
    assert len(cluster_specs.make_orbits()) == 6

    # `occevent_counter_params` is a dict that sets filters
    # See the `make_canonical_prim_periodic_occevents` documentation
    # for the list of options
    occevent_counter_params = {
        "max_cluster_size": 3,
        "print_state_info": False,
    }

    # `custom_occevents` is a list[occ_events.OccEvent]
    # that allows specifying custom OccEvent to include,
    # independent of the cluster_specs,
    # and not subject to filtering,
    # but still subject to removing duplicates
    custom_occevents = []

    canonical_occevents = occ_events.make_canonical_prim_periodic_occevents(
        system, cluster_specs, occevent_counter_params, custom_occevents
    )

    print_event = occ_events.OccEventPrinter(
        f=sys.stdout, system=system, coordinate_mode="cart"
    )

    for i, x in enumerate(canonical_occevents):
        print(i)
        print_event(x)
        print()

    assert len(canonical_occevents) == 24
