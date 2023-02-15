"""Cluster orbits

The :py:mod:`libcasm.clusterography` module supports construction of orbits of clusters

Primarily, this purpose of this module is to provide:

- :class:`~libcasm.clusterography.Cluster`, a class which represents a cluster of :class:`~libcasm.xtal.IntegralSiteCoordinate`.
- :func:`~libcasm.clusterography.make_prim_periodic_orbit`, a function that  constructs an orbit of :class:`~libcasm.clusterography.Cluster` obeying the periodic translation symmetry of a primitive crystal structure and allowed degrees of freedom (DoF) (:class:`~libcasm.xtal.Prim`).
- :func:`~libcasm.clusterography.make_local_orbit`, a function that constructs an orbit of :class:`~libcasm.clusterography.Cluster` without translation symmetry
- :class:`~libcasm.clusterography.ClusterSpecs`, a class which collects parameters controlling the generation of orbits, for either periodic or local-clusters, and generates all symmetrically distinct orbits

The :py:mod:`libcasm.clusterography` module has dependencies on:

- :py:mod:`libcasm.xtal`
- :py:mod:`libcasm.sym_info`

"""
from ._clusterography import *
from ._methods import *

__all__ = ["Cluster"]
