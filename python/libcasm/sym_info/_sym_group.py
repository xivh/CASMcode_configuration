import libcasm.xtal as xtal
import libcasm.sym_info._sym_info as _sym_info
from typing import Optional, NewType

SymGroupType = NewType("SymGroup", None)


class SymGroup(_sym_info.SymGroup):
    r"""Data structure holding group elements and other group info, such as group-subgroup relationships.

    The :class:`~libcasm.sym_info.SymGroup` class goes beyond storing a ``list[SymOp]`` to include the multiplication table and inverse element table. The SymGroup class may represent a head group, or a subgroup.

    When representing a subgroup, the head group can be obtained as a shared pointer from :func:`~libcasm.sym_info.SymGroup.head_group`. Pointers to subgroups are not stored by the head group. For subgroups, there is a list of indices indicating which element in the head group each subgroup element corresponds to (:func:`~libcasm.sym_info.SymGroup.head_group_index`).


    The :class:`libcasm.sym_info` package provides factory functions for the common use cases of constructing the prim factor group and prim point group:

    .. rubric:: Factory Functions

    +------------------------------------------------+------------------------------------------+
    | :func:`~libcasm.sym_info.make_factor_group`    | Make the prim factor group               |
    +------------------------------------------------+------------------------------------------+
    | :func:`~libcasm.sym_info.make_point_group`     | Make the prim point group                |
    +------------------------------------------------+------------------------------------------+


    Parameters
    ----------
    element : list[~libcasm.xtal.SymOp]
        List of elements of the symmetry group.
    multiplication_table : list[list[int]]
        The group multiplication table, where ``k == multiplication_table[i][j]`` indicates that ``element[k] == element[i] * element[j]``.

    """
    def __init__(self, element: list[xtal.SymOp],
                 multiplication_table: list[list[int]]):
        super().__init__(element, multiplication_table)

    def conjugacy_classes(self) -> list[list[int]]:
        r"""
        Return conjugacy classes

        Returns
        -------
        conjugacy_classes : list[list[int]]
            The value ``conjugacy_classes[i][j]`` is the index of the j-th element in conjugacy class i.

        """
        return super().conjugacy_classes()

    def elements(self) -> list[xtal.SymOp]:
        r"""Returns the elements list."""
        return super().elements()

    def head_group(self) -> Optional[SymGroupType]:
        r"""Return the head group of a subgroup

        Returns
        -------
        head_group: Optional[SymGroup]
            If this is a subgroup, returns the head group. Otherwise, returns None."""
        return super().head_group()

    def head_group_index(self) -> list[int]:
        r"""Return the list of head group indices

        Returns
        -------
        head_group_index: list[int]
            Specifies the index into the head group elements list corresponding to each element in this group.
        """
        return super().head_group_index()

    def inv(self, i: int) -> int:
        r"""Return the index of the inverse element.

        Parameters
        ----------
        i : int
            Index of an element in this group

        Returns
        -------
        k: int
            The index in this group's elements list of the inverse of the i-th element in this group.
        """
        return super().inv(i)

    def inverse_index(self) -> list[int]:
        r"""Return the inverse element indices

        Returns
        -------
        inverse_index: int
            The index in this group's elements list of each element in this group.
        """
        return super().inverse_index()

    def is_subgroup(self) -> bool:
        r"""Return True if this is a subgroup"""
        return super().is_subgroup()

    def make_subgroup(
            self,
            head_group_index: set[int],
            element: Optional[list[xtal.SymOp]] = None) -> SymGroupType:
        r"""Construct a subgroup

        Parameters
        ----------
        head_group_index : set[int]
            Indices of elements to include in the subgroup.
        element : Optional[list[~libcasm.xtal.SymOp]]
            List of elements in the subgroup. This is optional, and if not provided the head group elements corresponding to the head_group_index will be used. It may be provided to specify particular elements for a subgroup of a factor group, such as the elements of a cluster group. Must be listed in ascending head group index order.
        """
        return super().is_subgroup()

    def mult(self, i: int, j: int) -> int:
        r"""Return the element product index

        Parameters
        ----------
        i : int
            Index of an element in this group

        j : int
            Index of an element in this group

        Returns
        -------
        k: int
            The element product index form the group multiplication table, where ``k == multiplication_table[i][j]`` indicates that ``element[k] == element[i] * element[j]``.
        """
        return super().inv(i)

    def multiplication_table(self) -> list[list[int]]:
        r"""Return the group multiplication table

        Returns
        -------
        multiplication_table: list[list[int]]
            The group multiplication table, where ``k == multiplication_table[i][j]`` indicates that ``element[k] == element[i] * element[j]``.
        """
        return super().inv(i)
