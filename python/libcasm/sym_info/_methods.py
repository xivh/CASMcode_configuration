import libcasm.sym_info._sym_info as _sym_info
import libcasm.xtal as xtal
from libcasm.sym_info._sym_group import SymGroup


def make_factor_group(xtal_prim: xtal.Prim) -> SymGroup:
    """
    Construct the prim factor group as a SymGroup

    Parameters
    ----------
    xtal_prim: ~libcasm.xtal.Prim
        The prim

    Returns
    -------
    factor_group : SymGroup
        The group which leaves `xtal_prim` invariant
    """
    return _sym_info.make_factor_group(xtal_prim)


def make_point_group(xtal_prim: xtal.Prim, factor_group: SymGroup) -> SymGroup:
    """
    Construct the prim point group as a SymGroup

    Parameters
    ----------
    xtal_prim: ~libcasm.xtal.Prim
        The prim

    factor_group: SymGroup
        The factor group of `prim`, as calculated by
        :func:`~libcasm.sym_info.make_factor_group`.

    Returns
    -------
    point_group : SymGroup
        The prim point group operations, constructed by removing the translation from
        prim factor group operations. Degenerate symmetry operations are not added. The
        resulting point_group is its own head group, not a subgroup of the input
        factor_group.
    """
    return _sym_info.make_point_group(xtal_prim, factor_group)
