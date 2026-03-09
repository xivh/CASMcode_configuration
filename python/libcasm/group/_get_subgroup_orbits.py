import json

import libcasm.group as casmgroup
from libcasm.group.sqlite_cache import UserCache


def get_cyclic_subgroup_orbits(
    symgroup: casmgroup.Group,
) -> list[list[int]]:
    """Get cyclic subgroup orbits for a given symmetry group.

    Parameters
    ----------
    symgroup: casmgroup.Group
        The symmetry group for which to compute cyclic subgroup orbits.

    Returns
    -------
    cyclic_subgroup_orbits: list[list[int]]
        All cyclic subgroup orbits.

    """
    head_group = symgroup.head_group
    if head_group is None:
        head_group = symgroup
    indices = set(symgroup.head_group_index)

    subset = casmgroup.Subset(group=head_group, indices=indices)
    subset.cyclic_subgroups()
    cyclic_subgroup_orbits = subset.cyclic_subgroup_orbits()

    return cyclic_subgroup_orbits


def get_all_subgroup_orbits(
    symgroup: casmgroup.Group,
) -> list[list[list[int]]]:
    """Get subgroup orbits for a given symmetry group, with caching.

    Notes
    -----

    A :class:`~libcasm.group.Subset` object is created from the head group and indices
    of the given symmetry group, and the subgroup orbits are computed from this Subset
    object. The result is completely determined by the multiplication table of the head
    group and the indices of the elements in the subgroup. Therefore, the result
    is cached in ``~/.config/casm/sqlite_cache/subgroup_orbits.db`` using
    :class:`~libcasm.group.sqlite_cache.UserCache` to avoid redundant calculations.

    The cache key is constructed using:

    .. code-block:: Python

        head_group = symgroup.head_group
        if head_group is None:
            head_group = symgroup
        indices = set(symgroup.head_group_index)
        data = dict(
            multiplication_table=head_group.multiplication_table,
            indices=list(indices),
        )
        key = json.dumps(data, sort_keys=True)

    The stored value is:

    .. code-block:: Python

        subset = casmgroup.Subset(group=head_group, indices=indices)
        subset.all_subgroups()
        subgroup_orbits = subset.all_subgroup_orbits()
        value = json.dumps(subgroup_orbits, sort_keys=True)


    Parameters
    ----------
    symgroup: casmgroup.Group
        The symmetry group for which to compute subgroup orbits.

    Returns
    -------
    subgroup_orbits: list[list[list[int]]]
        All subgroup orbits.

    """
    ucache = UserCache()
    head_group = symgroup.head_group
    if head_group is None:
        head_group = symgroup
    indices = set(symgroup.head_group_index)

    data = dict(
        multiplication_table=head_group.multiplication_table,
        indices=list(indices),
    )
    key = json.dumps(data, sort_keys=True)
    value = ucache.get(cache="subgroup_orbits", key=key)
    if value is None:
        subset = casmgroup.Subset(group=head_group, indices=indices)
        subset.all_subgroups()
        subgroup_orbits = subset.all_subgroup_orbits()
        value = json.dumps(subgroup_orbits, sort_keys=True)
        ucache.store(cache="subgroup_orbits", key=key, value=value)
    else:
        subgroup_orbits = json.loads(value)
    return subgroup_orbits
