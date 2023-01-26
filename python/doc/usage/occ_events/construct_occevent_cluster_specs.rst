Constructing local-cluster orbits
=================================


Construct a prototype OccEvent
------------------------------

To generate local-cluster orbits or local basis functions, one prototype :class:`~libcasm.occ_events.OccEvent` from the orbit of symmetrically equivalent :class:`~libcasm.occ_events.OccEvent` must be provided. For convenience in checking for unique :class:`~libcasm.occ_events.OccEvent`, it is useful to construct an orbit of :class:`~libcasm.occ_events.OccEvent` and use the :class:`~libcasm.occ_events.OccEvent` that is sorted into the initial position. This can be done as follows:

.. code-block:: Python

    import libcasm.occ_events as occ_events
    import libcasm.sym_info as sym_info

    # Construct the Prim
    xtal_prim = xtal.Prim(...)

    # Construct one example of the OccEvent
    occ_event = occ_events.OccEvent(...)

    # Construct an orbit of OccEvent, and choose the prototype
    # (the first element) as the phenomenal OccEvent about which
    # to construct a local basis set
    #
    # Note the use of sym_info.make_factor_group:
    prim_factor_group = sym_info.make_factor_group(xtal_prim)
    occevent_symgroup_rep = occ_events.make_occevent_symgroup_rep(
        prim_factor_group.elements(), xtal_prim)
    occevent_orbit = occ_events.make_prim_periodic_orbit(
        occ_event,
        occevent_symgroup_rep)
    prototype_occ_event = occevent_orbit[0]


Save/load OccEvent
------------------

A standard location to save a prototype event for future use is in:

- <CASM project directory> / events / event.<name_of_event> / event.json

The functions :func:`~libcasm.occ_events.save_occevent` and :func:`~libcasm.occ_events.load_occevent` methods can be used to save :class:`~libcasm.occ_events.OccEvent` to the standard location and later load them:

.. code-block:: Python

    # root: pathlib.path
    # prototype_occ_event: Prototype OccEvent

    # The OccSystem provides index conversions
    system = occ_events.OccSystem(xtal_prim)

    # Save an OccEvent:
    occ_events.save_occevent(root, "1NN_A_Va", prototype_occ_event, system)

    # Load an OccEvent:
    loaded_occ_event = occ_events.load_occevent(root, "1NN_A_Va", system)

    assert loaded_occ_event == prototype_occ_event


ClusterSpecs for local-cluster orbits
-------------------------------------

The subgroup of the prim factor group that leaves an :class:`~libcasm.occ_events.OccEvent` invariant is the generating group for local-cluster orbits and local basis functions for properties of that event.

The :class:`~libcasm.clusterography.ClusterSpecs` class encapsulates all the parameters needed for constructing cluster orbits. A :class:`~libcasm.clusterography.ClusterSpecs` object with the correct generating group can be constructed using :func:`~libcasm.occ_events.make_occevent_cluster_specs`:

.. code-block:: Python

    # xtal_prim: xtal.Prim

    # Construct ClusterSpecs, with generating group equal to
    # the invariant group of prototype_occ_event
    cluster_specs = occ_events.make_occevent_cluster_specs(
        xtal_prim=xtal_prim,
        phenomenal_occ_event=prototype_occ_event,
        max_length=[0.0, 0.0],
        cutoff_radius=[0.0, 2.01])


Local-cluster orbits
--------------------

Once the :class:`~libcasm.clusterography.ClusterSpecs` instance is constructed, local-cluster orbits can be generated using :func:`~libcasm.clusterography.ClusterSpecs.make_orbits`:

.. code-block:: Python

    # Construct local cluster orbits
    local_cluster_orbits = cluster_specs.make_orbits()

Local basis sets
----------------

The :class:`~libcasm.clusterography.ClusterSpecs` instance can be output to JSON for use as input for constructing local basis sets using :func:`~libcasm.clusterography.ClusterSpecs.to_dict`:

.. code-block:: Python

    # Output cluster specs JSON for local basis set construction
    cluster_specs_json = cluster_specs.to_dict()


OccEvent invariant group
------------------------

The subgroup of the prim factor group that leaves an :class:`~libcasm.occ_events.OccEvent` invariant is the generating group for local basis functions of properties of that event. It can be constructed explicitly using :func:`~libcasm.occ_events.make_occevent_group`:

.. code-block:: Python

    import libcasm.occ_events as occ_events
    import libcasm.sym_info as sym_info

    # xtal_prim: xtal.Prim
    # occ_event: occ_events.OccEvent

    # Note the use of sym_info.make_factor_group:
    prim_factor_group = sym_info.make_factor_group(xtal_prim)
    occevent_symgroup_rep = occ_events.make_occevent_symgroup_rep(
        prim_factor_group.elements(), xtal_prim)

    # Construct the group which leaves the phenomenal OccEvent invariant
    invariant_group = occ_events.make_occevent_group(
        occ_event=prototype_occ_event,
        group=prim_factor_group,
        lattice=xtal_prim.lattice(),
        occevent_symgroup_rep=occevent_symgroup_rep)

The objects ``prim_factor_group`` and ``invariant_group`` are instances of :class:`~libcasm.sym_info.SymGroup`, with the relationship that ``invariant_group`` is a sub-group of ``prim_factor_group``, which is called the "head group". The class :class:`~libcasm.sym_info.SymGroup` provides more information than a simple ``list[libcasm.xtal.SymOp]``, including the multiplication table and the head group indices of the sub-group operations.
