from typing import Optional

import numpy as np

import libcasm.configuration as casmconfig

from ._enumerate import (
    SuperlatticeEnumBase,
)


class SuperlatticeEnum:
    """Enumerate symmetrically distinct supercell lattices

    Supercell lattices satisfy:

    .. code-block:: Python

        S = L @ T,

    where `S` and `L` are, respectively, the superlattice and unit lattice vectors
    as columns of shape=(3, 3) matrices, and `T` is an integer shape=(3,3)
    transformation matrix.

    Superlattices `S1` and `S2` are symmetrically equivalent if there exists `p` and
    `A` such that:

    .. code-block:: Python

      S2 = p.matrix() @ S1 @ A,

    where `p` is an element in the crystal point group of the prim, and `A` is a
    unimodular matrix (integer matrix, with abs(det(A))==1).

    SuperlatticeEnum generates symmetrically distinct supercells by iterating over
    possible transformation matrices `T` using the hermite normal form and applying
    symmetry operations to check for uniqueness. Resulting superlattices are put into
    the CASM canonical lattice form.

    See also `Lattice Canonical Form`_.

    .. _`Lattice Canonical Form`: https://prisms-center.github.io/CASMcode_docs/formats/lattice_canonical_form/

    """

    def __init__(
        self,
        prim: casmconfig.Prim,
    ):
        """
        .. rubric:: Constructor

        Parameters
        ----------
        prim: casmconfig.Prim
            The Prim
        """
        self._prim = prim

    @property
    def prim(self) -> casmconfig.Prim:
        """The Prim"""
        return self._prim

    def by_volume(
        self,
        max: int,
        min: int = 1,
        unit_cell: Optional[np.ndarray] = None,
        dirs: str = "abc",
        diagonal_only: bool = False,
        fixed_shape: bool = False,
    ):
        """Yields symmetrically distinct superlattices for a range of volumes

        Parameters
        ----------
        max : int
            The maximum volume superlattice to enumerate. The volume is measured
            relative the unit cell being used to generate supercells.
        min : int, default=1
            The minimum volume superlattice to enumerate. The volume is measured
            relative the unit cell being used to generate supercells.
        dirs : str, default="abc"
            A string indicating which lattice vectors to enumerate over. Some
            combination of 'a', 'b', and 'c', where 'a' indicates the first lattice
            vector of the unit cell, 'b' the second, and 'c' the third.
        unit_cell: Optional[np.ndarray] = None,
            An integer shape=(3,3) transformation matrix `U` allows specifying an
            alternative unit cell that can be used to generate superlattices of the
            form `S = (L @ U) @ T`. If None, `U` is set to the identity matrix.
        diagonal_only: bool = False
            If true, restrict :math:`T` to diagonal matrices.
        fixed_shape: bool = False
            If true, restrict :math:`T` to diagonal matrices with diagonal coefficients
            :math:`[m, 1, 1]` (1d), :math:`[m, m, 1]` (2d), or :math:`[m, m, m]` (3d),
            where the dimension is determined from `len(dirs)`.

        Yields
        ------
        superlattice: libcasm.xtal.Lattice
            A :class:`~libcasm.xtal.Lattice`, guaranteed to be in canonical
            form.
        """
        superlat_enum = SuperlatticeEnumBase(
            unit_lattice=self.prim.xtal_prim.lattice(),
            point_group=self.prim.crystal_point_group.elements,
            max_volume=max,
            min_volume=min,
            dirs=dirs,
            unit_cell=unit_cell,
            diagonal_only=diagonal_only,
            fixed_shape=fixed_shape,
        )

        while superlat_enum.is_valid():
            yield superlat_enum.value()
            superlat_enum.advance()
