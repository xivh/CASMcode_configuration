from typing import Optional

import numpy as np

from libcasm.configuration import (
    Configuration,
    Prim,
    SupercellSet,
    copy_configuration,
    make_canonical_configuration,
    make_equivalent_supercells,
)

from ._ScelEnum import ScelEnum


class SuperConfigEnum:
    def __init__(
        self,
        prim: Prim,
        supercell_set: Optional[SupercellSet] = None,
    ):
        """
        .. rubric:: Constructor

        Parameters
        ----------
        prim: libcasm.configuration.Prim
            The Prim
        supercell_set: Optional[libcasm.configuration.SupercellSet] = None
            If not None, generated :class:`Supercell` are constructed by
            adding in the :class:`~SupercellSet`.
        """
        self._prim = prim
        self._supercell_set = supercell_set

    @property
    def prim(self) -> Prim:
        """The Prim"""
        return self._prim

    @property
    def supercell_set(self) -> Optional[SupercellSet]:
        """If not None, generated :class:`~casmconfig.Supercell` are constructed by
        adding in the :class:`~casmconfig.SupercellSet`."""
        return self._supercell_set

    def by_supercell(
        self,
        motif: Configuration,
        max: int,
        min: int = 1,
        unit_cell: Optional[np.ndarray] = None,
        dirs: str = "abc",
        diagonal_only: bool = False,
        fixed_shape: bool = False,
    ):
        """Make super configurations of the motif, without changing the orientation of
        the motif

        Method:

        - For each supercell generated by `ScelEnum.by_volume`, equivalent supercells
          with respect to the prim factor group are generated.
        - If the motif tiles a supercells exactly, a super configuration is
          generated in that supercell.
        - If the super configuration is unique, as determined by checking a running
          list of canonical super configurations in the canonical supercell, then it
          is yielded.

        Parameters
        ----------
        motif: libcasm.configuraiton.Configuration
            The configuration to generate super configurations from.
        max : Optional[int] = None
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
        config: casmconfig.Configuration
            A :class:`~casmconfig.Configuration`. Configurations might not be in the
            canonical supercell.
        """
        if motif.supercell.prim is not self.prim:
            raise ValueError(
                "Error in SuperConfigEnum.with_fixed_orientation: "
                "motif.supercell.prim must be the same as prim used at construction"
            )

        supercell_set = self.supercell_set

        scel_enum = ScelEnum(
            prim=motif.supercell.prim,
            supercell_set=supercell_set,
        )
        for supercell in scel_enum.by_volume(
            max=max,
            min=min,
            unit_cell=unit_cell,
            dirs=dirs,
            diagonal_only=diagonal_only,
            fixed_shape=fixed_shape,
        ):
            _super_configurations = []
            _canonical_super_configurations = []
            for equiv in make_equivalent_supercells(supercell):
                is_superlat, T = equiv.superlattice.is_superlattice_of(
                    motif.supercell.superlattice
                )
                if not is_superlat:
                    continue
                super_configuration = copy_configuration(motif, equiv)
                canonical_super_configuration = make_canonical_configuration(
                    configuration=super_configuration,
                    in_canonical_supercell=True,
                )
                if canonical_super_configuration not in _canonical_super_configurations:
                    _super_configurations.append(super_configuration)
                    _canonical_super_configurations.append(
                        canonical_super_configuration
                    )
                    yield super_configuration.copy()
