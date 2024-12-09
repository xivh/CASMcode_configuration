from typing import Optional

import numpy as np

import libcasm.configuration as casmconfig
import libcasm.xtal as xtal

from ._SuperlatticeEnum import SuperlatticeEnum


class ScelEnum:
    """Enumerate symmetrically distinct supercells

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

    ScelEnum generates symmetrically distinct supercells by iterating over possible
    transformation matrices `T` using the hermite normal form and applying symmetry
    operations to check for uniqueness. Resulting superlattices are put into the CASM
    canonical lattice form and used to construct a
    :class:`~libcasm.configuration.Supercell`.

    See also `Lattice Canonical Form`_.

    .. _`Lattice Canonical Form`: https://prisms-center.github.io/CASMcode_docs/formats/lattice_canonical_form/


    .. rubric:: Example usage

    The following example uses ScelEnum to enumerate symmetrically distinct supercells
    with volume up to 4 times the prim lattice volume, for an FCC prim.

    .. code-block:: Python

        import libcasm.configuration as casmconfig
        import libcasm.enumerate as casmenum
        import libcasm.xtal.prims as xtal_prims

        xtal_prim = xtal_prims.FCC(
            r=0.5,
            occ_dof=["A", "B"],
        )
        prim = casmconfig.Prim(xtal_prim)
        scel_enum = casmenum.ScelEnum(
            prim=prim,
        )
        for i, supercell in enumerate(scel_enum.by_volume(max=4)):
            print(f"Supercell {i}:")
            print("Transformation matrix to supercell:")
            print(supercell.transformation_matrix_to_super)
            print("Superlattice (column_vector_matrix):")
            print(supercell.superlattice.column_vector_matrix())
            print()

    ::

        Supercell 0:
        Transformation matrix to supercell:
        [[1 0 0]
         [0 1 0]
         [0 0 1]]
        Superlattice (column_vector_matrix):
        [[0.         0.70710678 0.70710678]
         [0.70710678 0.         0.70710678]
         [0.70710678 0.70710678 0.        ]]

        Supercell 1:
        Transformation matrix to supercell:
        [[ 1  0  1]
         [ 0  1  1]
         [-1 -1  0]]
        Superlattice (column_vector_matrix):
        [[-0.70710678  0.          0.70710678]
         [ 0.         -0.70710678  0.70710678]
         [ 0.70710678  0.70710678  1.41421356]]
        ...


    Often, ScelEnum is used as part of a configuration enumeration procedure and
    it is convenient that generated supercells are also inserted into a
    :class:`~libcasm.configuration.SupercellSet`. This is done automatically if
    the :class:`~libcasm.configuration.SupercellSet` is provided to the ScelEnum
    constructor.

    .. code-block:: Python

        import libcasm.configuration as casmconfig
        import libcasm.enumerate as casmenum
        import libcasm.xtal.prims as xtal_prims

        xtal_prim = xtal_prims.FCC(
            r=0.5,
            occ_dof=["A", "B"],
        )
        prim = casmconfig.Prim(xtal_prim)
        supercell_set = casmconfig.SupercellSet(prim)
        scel_enum = casmenum.ScelEnum(
            prim=prim,
            supercell_set=supercell_set,
        )
        for i, supercell in enumerate(scel_enum.by_volume(max=4)):
            # ... do something ...
            pass
        print(f"Number of Supercells: {len(supercell_set)}")

    ::

        Number of Supercells: 13

    """

    def __init__(
        self,
        prim: casmconfig.Prim,
        supercell_set: Optional[casmconfig.SupercellSet] = None,
    ):
        """
        .. rubric:: Constructor

        Parameters
        ----------
        prim: casmconfig.Prim
            The Prim
        supercell_set: Optional[casmconfig.SupercellSet] = None
            If not None, generated :class:`~casmconfig.Supercell` are
            constructed by adding in the :class:`~casmconfig.SupercellSet`.
        """
        self._prim = prim
        self._supercell_set = supercell_set

    @property
    def prim(self) -> casmconfig.Prim:
        """The Prim"""
        return self._prim

    @property
    def supercell_set(self) -> Optional[casmconfig.SupercellSet]:
        """If not None, generated :class:`~casmconfig.Supercell` are constructed by
        adding in the :class:`~casmconfig.SupercellSet`.
        """
        return self._supercell_set

    def by_volume(
        self,
        max: int,
        min: int = 1,
        unit_cell: Optional[np.ndarray] = None,
        dirs: str = "abc",
        diagonal_only: bool = False,
        fixed_shape: bool = False,
    ):
        """Yields symmetrically distinct supercells for a range of volumes

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
        supercell: casmconfig.Supercell
            A :class:`~casmconfig.Supercell`, guaranteed to be in canonical
            form.
        """
        prim_lattice = self.prim.xtal_prim.lattice()
        superlat_enum = SuperlatticeEnum(
            prim=self.prim,
        )
        for superlattice in superlat_enum.by_volume(
            max=max,
            min=min,
            dirs=dirs,
            unit_cell=unit_cell,
            diagonal_only=diagonal_only,
            fixed_shape=fixed_shape,
        ):
            T = xtal.make_transformation_matrix_to_super(
                superlattice=superlattice,
                unit_lattice=prim_lattice,
            )
            if self.supercell_set is None:
                yield casmconfig.Supercell(
                    self.prim,
                    transformation_matrix_to_super=T,
                )
            else:
                record = self.supercell_set.add_by_transformation_matrix_to_super(
                    transformation_matrix_to_super=T,
                )
                yield record.supercell
