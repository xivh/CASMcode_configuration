"""Crystallographic symmetry

The libcasm.sym_info package includes:

- Methods to construct the prim factor group and point group
- The :class:`~libcasm.sym_info.SymGroup` class to store group elements, multiplication
  and inverse index tables, and group-subgroup relationships.

Currently, symmetry representations are not accessed in Python directly through
:py:mod:`libcasm.sym_info`, but can be obtained when necessary through, for example,
:class:`~libcasm.configuration.Prim` and
:func:`~libcasm.occ_events.make_occevent_symgroup_rep`.


The :py:mod:`libcasm.sym_info` module only has a dependency on
:py:mod:`libcasm.xtal`.
"""

from ._sym_info import (
    SymGroup,
    make_factor_group,
    make_point_group,
)
