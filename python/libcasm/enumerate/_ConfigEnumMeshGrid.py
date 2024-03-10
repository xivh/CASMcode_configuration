import copy
from typing import Optional

import numpy as np

import libcasm.counter
import libcasm.clexulator as casmclex
import libcasm.clusterography as casmclust
import libcasm.configuration as casmconfig
import libcasm.xtal as xtal


class ConfigEnumMeshGrid:
    """Enumerate configuration with continuous DoF values in a mesh grid"""

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
            If not None, generated :class:`~casmconfig.Supercell` are constructed by
            adding in the :class:`~casmconfig.SupercellSet`.
        """
        self._prim = prim
        self._supercell_set = supercell_set

        # Set and updated during an enumeration
        self._background = None
        self._sites = None
        self._enum_index = None

    @property
    def prim(self) -> casmconfig.Prim:
        """The Prim"""
        return self._prim

    @property
    def supercell_set(self) -> Optional[casmconfig.SupercellSet]:
        """If not None, generated :class:`~casmconfig.Supercell` are constructed by
        adding in the :class:`~casmconfig.SupercellSet`."""
        return self._supercell_set

    @property
    def background(self) -> Optional[casmconfig.Configuration]:
        """During enumeration, `background` is set to the current background
        configuration in which occupations are enumerated."""
        return self._background

    def by_dof_space(
        self,
        *xi,
        background: casmconfig.Configuration,
        dof_space: casmclex.DoFSpace,
    ):
        """Enumerate on a meshgrid specified with respect to a DoFSpace basis

        Parameters
        ----------
        x1, x2,..., xn : array_like
            1-D arrays representing the coordinates of a grid. x1 gives the coordinates
            along the first `dof_space` basis vector, x2 along the second `dof_space`
            vector, etc. There must be one array per `dof_space` basis vector.
        background: casmconfig.Configuration
            The background configuration on which enumeration takes place.
        dof_space: casmclex.DoFSpace
            Specifies the DoF being enumerated and the basis of the DoFSpace are the
            axes on which the meshgrid is constructed and enumeration takes place.
            For local DoF, the dof_space supercell must tile the background supercell.

        Yields
        ------
        config: casmconfig.Configuration
            A :class:`~casmconfig.Configuration` on a point in the DoF values meshgrid.
        """
        self._background = background

        if dof_space.dof_key == "occ":
            raise Exception('Invalid dof_space: dof_space.dof_key == "occ"')

        # Make the supercell for the DoFSpace
        if dof_space.transformation_matrix_to_super is None:
            dof_space_supercell = background.supercell()
        else:
            record = self.supercell_set.add_by_transformation_matrix_to_super(
                dof_space.transformation_matrix_to_super
            )
            dof_space_supercell = record.supercell

        # For
        # If DoF are global, or the background supercell is same as DoFSpace supercell,
        # then we can just copy DoF values onto the background.
        # Otherwise, we need to make a super-configuration and then copy DoF values
        make_superconfig = False
        if not dof_space.is_global:
            if background.supercell() != dof_space_supercell:
                if (
                    not background.supercell()
                    .superlattice()
                    .is_superlattice_of(dof_space_supercell.superlattice())
                ):
                    raise Exception(
                        "Error in ConfigEnumMeshGrid.by_dof_space: background supercell is not tiled by the dof_space supercell"
                    )
                make_superconfig = True

        dim = len(xi)
        counter = libcasm.counter.IntCounter(
            initial=[0] * dim,
            final=[len(x) - 1 for x in xi],
            increment=[1] * dim,
        )

        for indices in counter:
            # copy background
            generated_config = copy.copy(background)

            # generate dof space coordinate
            eta = np.array([xi[i] for i in indices])

            # make config with generated dof values in dof_space supercell
            # (all other DoF values have default value)
            dof_space_config = casmconfig.Configuration(dof_space_supercell)
            dof_space_config.set_dof_space_values(
                dof_space=dof_space,
                dof_space_coordinate=eta,
            )

            # make config with generated dof values in background supercell
            if make_superconfig:
                dof_space_config = casmconfig.copy_configuration(
                    dof_space_config,
                    background.supercell(),
                )

            # copy generated dof values into background supercell
            if dof_space.is_global:
                generated_config.set_global_dof_values(
                    key=dof_space.dof_key,
                    dof_values=dof_space_config.global_dof_values(dof_space.dof_key),
                )
            else:
                # set values site-by-site to respect dof_space.site_indices
                value = dof_space_config.local_dof_values(dof_space.dof_key)
                for l in dof_space.site_indices:
                    generated_config.set_local_dof_site_value(
                        key=dof_space.dof_key,
                        l=l,
                        dof_values=value[:, l],
                    )
            yield generated_config

        return None
