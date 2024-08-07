from typing import Optional

import numpy as np

import libcasm.clusterography as casmclust
import libcasm.configuration as casmconfig
import libcasm.xtal as xtal

from ._enumerate import (
    ConfigEnumAllOccupationsBase,
    make_distinct_cluster_sites,
)
from ._ScelEnum import ScelEnum
from ._SuperConfigEnum import SuperConfigEnum


def _make_sublat_sites(supercell: casmconfig.Supercell, sublats: set[int]):
    converter = supercell.site_index_converter
    sites = set()
    for i in range(supercell.n_sites):
        integral_site_coordinate = converter.integral_site_coordinate(i)
        if integral_site_coordinate.sublattice() in sublats:
            sites.add(i)
    return sites


class ConfigEnumAllOccupations:
    """Enumerate configuration occupations"""

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

    @property
    def sites(self) -> Optional[set[int]]:
        """During enumeration, `sites` is set to the linear site indices on which
        occupations are enumerated in the current background configuration."""
        return self._sites

    @property
    def enum_index(self) -> Optional[int]:
        """During enumeration, `enum_index` is incremented to count over the
        combinations of background configuration and sites on which enumeration is
        performed, starting from 0."""
        return self._enum_index

    def _begin(self):
        """Initialize values"""
        self._background = None
        self._sites = None
        self._enum_index = None

    def _set_motif(self, motif: Optional[casmconfig.Configuration] = None):
        """Check motif and use volume 1 default configuration if not provided"""
        if motif is None:
            T = np.eye(3, dtype=int)
            if self.supercell_set is not None:
                record = self.supercell_set.add_by_transformation_matrix_to_super(
                    transformation_matrix_to_super=T,
                )
                supercell = record.supercell
            else:
                supercell = casmconfig.Supercell(
                    prim=self.prim,
                    transformation_matrix_to_super=T,
                )
            motif = casmconfig.Configuration(supercell)
        return motif

    def _make_supercell_list(
        self,
        background: casmconfig.Configuration,
        supercells: Optional[dict] = None,
    ):
        """Make supercell list from `supercells` or just background supercell"""
        if supercells is None:
            if self.supercell_set is not None:
                self.supercell_set.add_supercell(background.supercell)
            supercell_list = [background.supercell]
        else:
            scel_enum = ScelEnum(
                prim=self.prim,
                supercell_set=self.supercell_set,
            )
            supercell_list = [
                supercell for supercell in scel_enum.by_volume(**supercells)
            ]
        return supercell_list

    def _make_super_backgrounds(
        self,
        background: casmconfig.Configuration,
        supercells: Optional[dict] = None,
    ):
        """Fill the background configuration into enumerated supercells, maintaining
        the orientation of the background configuration.

        Parameters
        ----------
        background: casmconfig.Configuration
            The background configuration to fill into enumerated supercells.
        supercells: Optional[dict] = None
            Parameters to forward to
            :func:`ScelEnum.by_volume <libcasm.enumerate.ScelEnum.by_volume>`
            to specify the supercells that the background configuration will be filled
            into. If None, only the exact background configuration is returned.

        Returns
        -------
        super_backgrounds: list[casmconfig.Configuration]
            A list of :class:`~casmconfig.Configuration` objects, each a distinct
            super configuration of the background.
        """
        if supercells is None:
            return [background.copy()]
        super_config_enum = SuperConfigEnum(
            prim=self.prim,
            supercell_set=self.supercell_set,
        )
        super_backgrounds = []
        for config in super_config_enum.by_supercell(
            motif=background,
            **supercells,
        ):
            super_backgrounds.append(config.copy())
        return super_backgrounds

    def _by_site(
        self,
        background: casmconfig.Configuration,
        sites: set[int],
        skip_non_primitive: bool,
        skip_equivalents: bool,
        use_background_invariant_group: bool,
        which_dofs: Optional[set[str]] = None,
    ):
        """Run the inner loop of enumerating occupations on sites in a background

        This method will set and update self.background, self.sites, and
        self.enum_index as it yields configurations.

        Parameters
        ----------
        background: casmconfig.Configuration
            During enumeration, `background` is set to the current background
            configuration in which occupations are enumerated.
        sites: set[int]
            During enumeration, `sites` is set to the linear site indices on which
            occupations are enumerated in the current background configuration.
        skip_non_primitive: bool
            If True, enumeration skips non-primitive configurations.
        skip_equivalents: bool
            If True, enumeration skips non-unique configurations with respect
            to either (i) the complete set of symmetry operations that leave the
            supercell lattice vectors invariant (these are canonical configurations),
            or (ii) the subgroup that leaves the background configuration invariant and
            does not mix the given sites and other sites.
        use_background_invariant_group: bool
            If True, check for uniqueness using the subgroup that leaves the background
            configuration invariant and does not mix the given sites and other sites;
            otherwise use the complete set of symmetry operations that leave the
            supercell lattice vectors invariant.
        which_dofs: Optional[set[str]] = None
            When ``skip_equivalents is True`` and
            ``use_background_invariant_group is True``, the names of the degrees of
            freedom (DoF) in the background that must be invariant. The default is that
            all DoF must be invariant.

        Yields
        ------
        config: casmconfig.Configuration
            A :class:`~casmconfig.Configuration`.
        """
        if which_dofs is None:
            which_dofs = set(["all"])
        self._background = background
        self._sites = sites
        if self._enum_index is None:
            self._enum_index = 0
        else:
            self._enum_index += 1
        config_enum = ConfigEnumAllOccupationsBase(
            background=background,
            sites=sites,
        )
        if skip_equivalents:
            if use_background_invariant_group:
                canonicalization_group = casmconfig.make_invariant_subgroup(
                    configuration=background,
                    site_indices=sites,
                    which_dofs=which_dofs,
                )
            else:
                canonicalization_group = None
        while config_enum.is_valid():
            if skip_non_primitive and not casmconfig.is_primitive_configuration(
                configuration=config_enum.value()
            ):
                config_enum.advance()
                continue
            if skip_equivalents and not casmconfig.is_canonical_configuration(
                configuration=config_enum.value(),
                subgroup=canonicalization_group,
            ):
                config_enum.advance()
                continue
            yield config_enum.value()
            config_enum.advance()

    def by_supercell(
        self,
        max: Optional[int] = None,
        min: int = 1,
        unit_cell: Optional[np.ndarray] = None,
        dirs: str = "abc",
        diagonal_only: bool = False,
        fixed_shape: bool = False,
        supercells: Optional[dict] = None,
        motif: Optional[casmconfig.Configuration] = None,
        skip_non_primitive: bool = True,
        skip_non_canonical: bool = True,
    ):
        """Enumerate all occupations in a series of enumerated supercells

        Parameters
        ----------
        max: Optional[int] = None
            The maximum volume superlattice to enumerate. The volume is measured
            relative the unit cell being used to generate supercells (i.e. the
            determinant of the `unit_cell` parameter). Is required if `supercells` is
            None.
        min: int = 1
            The minimum volume superlattice to enumerate. The volume is measured
            relative the unit cell being used to generate supercells.
        dirs: str = "abc"
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
        supercells: Optional[dict] = None
            Parameters to forward to
            :func:`ScelEnum.by_volume <libcasm.enumerate.ScelEnum.by_volume>` to
            specify the supercells that the motif configuration will be filled into.

            .. deprecated:: 2.0a5
                Give the
                :func:`ScelEnum.by_volume <libcasm.enumerate.ScelEnum.by_volume>`
                parameters directly, or by using `**supercells`, instead of using this
                parameter.

        motif: Optional[casmconfig.Configuration] = None
            The background configuration on which enumeration takes place. The motif is
            filled into each supercell using
            :func:`~libcasm.configuration.make_distinct_super_configurations`. If the
            motif does not tile exactly into a supercell that supercell is skipped.
            Providing a motif allows enumerating all occupations with other degrees of
            freedom (DoF) fixed. If None, the default configuration in the volume 1
            supercell is used.

            .. deprecated:: 2.0a5
                Use the method
                :func:`~libcasm.enumerate.ConfigEnumAllOccupations.by_supercell_with_continuous_dof`
                instead.

        skip_non_primitive: bool = True
            If True, enumeration skips non-primitive configurations. All DoF are
            included in the check for primitive configurations.
        skip_non_canonical: bool = True
            If True, enumeration skips non-canonical configurations with respect
            to the symmetry operations that leave the supercell lattice vectors
            invariant.

        Yields
        ------
        config: casmconfig.Configuration
            A :class:`~casmconfig.Configuration`. All configurations are in the
            canonical supercell.
        """
        import warnings

        if max is None:
            if supercells is None or "max" not in supercells:
                raise ValueError(
                    "Error in ConfigEnumAllOccupations.by_supercell: "
                    "`max` is required. Could not be obtained from supercells."
                )
            else:
                max = supercells["max"]
                warnings.warn(
                    "The `max` parameter is required. "
                    "Obtaining max from supercells, but the 'supercells' parameter is "
                    "deprecated. Give the ScelEnum.by_volume parameters directly, or "
                    "by using `**supercells`, instead of using 'supercells'.",
                    DeprecationWarning,
                    stacklevel=2,
                )

        _supercells = dict(
            max=max,
            min=min,
            unit_cell=unit_cell,
            dirs=dirs,
            diagonal_only=diagonal_only,
            fixed_shape=fixed_shape,
        )

        if supercells is not None:
            warnings.warn(
                "The 'supercells' parameter is deprecated. Give the "
                "ScelEnum.by_volume parameters directly, or by using "
                "`**supercells`, instead of using 'supercells'.",
                DeprecationWarning,
                stacklevel=2,
            )
            _supercells.update(supercells)

        if motif is not None:
            import warnings

            warnings.warn(
                "The 'motif' parameter is deprecated. Use the method "
                "'ConfigEnumAllOccupations.by_supercell_with_continuous_dof' instead.",
                DeprecationWarning,
                stacklevel=2,
            )
            for config in self.by_supercell_with_continuous_dof(
                source=motif,
                **_supercells,
                skip_non_primitive=skip_non_primitive,
                skip_equivalents=skip_non_canonical,
            ):
                yield config
            return

        self._begin()
        scel_enum = ScelEnum(
            prim=self.prim,
            supercell_set=self.supercell_set,
        )
        for supercell in scel_enum.by_volume(
            **_supercells,
        ):
            sites = set(range(supercell.n_sites))
            default_config = casmconfig.Configuration(supercell)
            for config in self._by_site(
                background=default_config,
                sites=sites,
                skip_non_primitive=skip_non_primitive,
                skip_equivalents=skip_non_canonical,
                use_background_invariant_group=False,
            ):
                yield config

    def by_supercell_list(
        self,
        supercells: list[casmconfig.Supercell],
        skip_non_primitive: bool = True,
        skip_non_canonical: bool = True,
    ):
        """Enumerate all occupations in a list of supercells explicitly provided

        Parameters
        ----------
        supercells: list[casmconfig.Supercell]
            An explicit list of supercells in which to perform enumeration.
        skip_non_primitive: bool = True
            If True, enumeration skips non-primitive configurations. All DoF are
            included in the check for primitive configurations.
        skip_non_canonical: bool = True
            If True, enumeration skips non-canonical configurations with respect
            to the symmetry operations that leave the supercell lattice vectors
            invariant.

        Yields
        ------
        config: casmconfig.Configuration
            A :class:`~casmconfig.Configuration`.
        """
        self._begin()
        for supercell in supercells:
            if self.supercell_set is not None:
                self.supercell_set.add_supercell(supercell)
            sites = set(range(supercell.n_sites))
            default_config = casmconfig.Configuration(supercell)
            for config in self._by_site(
                background=default_config,
                sites=sites,
                skip_non_primitive=skip_non_primitive,
                skip_equivalents=skip_non_canonical,
                use_background_invariant_group=False,
            ):
                yield config

    def by_supercell_with_continuous_dof(
        self,
        source: casmconfig.Configuration,
        max: int,
        min: int = 1,
        unit_cell: Optional[np.ndarray] = None,
        dirs: str = "abc",
        diagonal_only: bool = False,
        fixed_shape: bool = False,
        skip_non_primitive: bool = True,
        skip_equivalents: bool = True,
    ):
        """Enumerate all occupations in a series of enumerated supercells, with
        non-default continuous DoF

        Parameters
        ----------
        source: casmconfig.Configuration
            Any global continuous DoF from `source` are set to the same value in all
            generated supercells. Any local continuous DoF of `source` are filled into
            each supercell using :func:`~libcasm.configuration.copy_configuration`. If
            `source` does not tile exactly into a supercell in its canonical form, or
            any equivalent supercell with respect to the prim factor group, then that
            supercell is skipped.
        max: Optional[int] = None
            The maximum volume superlattice to enumerate. The volume is measured
            relative the unit cell being used to generate supercells (i.e. the
            determinant of the `unit_cell` parameter). Is required if `supercells` is
            None.
        min: int = 1
            The minimum volume superlattice to enumerate. The volume is measured
            relative the unit cell being used to generate supercells.
        dirs: str = "abc"
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
        skip_non_primitive: bool = True
            If True, enumeration skips non-primitive configurations. All DoF are
            included in the check for primitive configurations.
        skip_equivalents: bool = True
            If True, enumeration skips non-unique configurations with respect
            to the subgroup that leaves the supercell lattice vectors
            and continuous DoF of the background configuration invariant.

        Yields
        ------
        config: casmconfig.Configuration
            A :class:`~casmconfig.Configuration`.
        """
        self._begin()
        supercells = dict(
            max=max,
            min=min,
            unit_cell=unit_cell,
            dirs=dirs,
            diagonal_only=diagonal_only,
            fixed_shape=fixed_shape,
        )
        super_backgrounds = self._make_super_backgrounds(
            background=source,
            supercells=supercells,
        )

        # Get continuous DoF
        which_dofs = set()
        for _global_dof in self.prim.xtal_prim.global_dof():
            which_dofs.add(_global_dof.dofname())
        for _sublattice_dof in self.prim.xtal_prim.local_dof():
            for _local_dof in _sublattice_dof:
                which_dofs.add(_local_dof.dofname())

        for super_background in super_backgrounds:
            sites = set(range(super_background.supercell.n_sites))
            for config in self._by_site(
                background=super_background,
                sites=sites,
                skip_non_primitive=skip_non_primitive,
                skip_equivalents=skip_equivalents,
                use_background_invariant_group=True,
                which_dofs=which_dofs,
            ):
                yield config

    def by_linear_site_indices(
        self,
        background: casmconfig.Configuration,
        sites: set[int],
        skip_non_primitive: bool = False,
        skip_equivalents: bool = True,
    ):
        """Enumerate occupation perturbations of a background configuration on
        specified sites

        Parameters
        ----------
        background: casmconfig.Configuration
            The background configuration on which enumeration takes place.
        sites: set[int]
            The sites, as linear site indices, on which occupations are enumerated,
            while all other degrees of freedom (DoF) are fixed.
        skip_non_primitive: bool = False
            If True, enumeration skips non-primitive configurations. All DoF are
            included in the check for primitive configurations.
        skip_equivalents: bool = True
            If True, enumeration skips equivalent configurations with respect
            to the subgroup that leaves the background configuration invariant
            and does not mix the given sites and other sites.

        Yields
        ------
        config: casmconfig.Configuration
            A :class:`~casmconfig.Configuration`.
        """
        self._begin()
        for config in self._by_site(
            background=background,
            sites=sites,
            skip_non_primitive=skip_non_primitive,
            skip_equivalents=skip_equivalents,
            use_background_invariant_group=True,
        ):
            yield config

    def by_integral_site_coordinates(
        self,
        background: casmconfig.Configuration,
        sites: list[xtal.IntegralSiteCoordinate],
        skip_non_primitive: bool = False,
        skip_equivalents: bool = True,
    ):
        """Enumerate occupation perturbations of a background configuration on
        specified sites

        Parameters
        ----------
        background: casmconfig.Configuration
            The background configuration on which enumeration takes place.
        sites: list[xtal.IntegralSiteCoordinate]
            The sites on which occupations are enumerated, while all other degrees of
            freedom (DoF) are fixed.
        skip_non_primitive: bool = False
            If True, enumeration skips non-primitive configurations. All DoF are
            included in the check for primitive configurations.
        skip_equivalents: bool = True
            If True, enumeration skips equivalent configurations with respect
            to the subgroup that leaves the background configuration invariant
            and does not mix the given sites and other sites.

        Yields
        ------
        config: casmconfig.Configuration
            A :class:`~casmconfig.Configuration`.
        """
        self._begin()
        converter = background.supercell.site_index_converter
        site_indices = set([converter.linear_site_index(site) for site in sites])
        for config in self._by_site(
            background=background,
            sites=site_indices,
            skip_non_primitive=skip_non_primitive,
            skip_equivalents=skip_equivalents,
            use_background_invariant_group=True,
        ):
            yield config

    def by_sublattice(
        self,
        background: casmconfig.Configuration,
        sublats: set[int],
        supercells: Optional[dict] = None,
        skip_non_primitive: bool = True,
        skip_equivalents: bool = True,
    ):
        """Enumerate occupation perturbations of a background configuration on
        specified sublattices


        Parameters
        ----------
        background: casmconfig.Configuration,
            The background configuration on which enumeration takes place. This allows
            fixing other DoF and enumerating all occupations and/or specifying the
            starting occupation which is perturbed.
        sublats: set[int]
            If not None, occupation enumeration only occurs on the specified
            sublattices.
        supercells: Optional[dict] = None
            Parameters to forward to
            :func:`ScelEnum.by_volume <libcasm.enumerate.ScelEnum.by_volume>`
            to specify the supercells that the background configuration will be filled
            into. If None, only the exact background configuration is perturbed.
        skip_non_primitive: bool = False
            If True, enumeration skips non-primitive configurations. All DoF are
            included in the check for primitive configurations.
        skip_equivalents: bool = True
            If True, enumeration skips equivalent configurations with respect
            to the subgroup that leaves the background configuration invariant
            and does not mix the selected sublattice sites and other sites.

        Yields
        ------
        config: casmconfig.Configuration
            A :class:`~casmconfig.Configuration`.
        """
        self._begin()
        super_backgrounds = self._make_super_backgrounds(
            background=background,
            supercells=supercells,
        )
        for super_background in super_backgrounds:
            sublat_sites = _make_sublat_sites(super_background.supercell, sublats)
            for config in self._by_site(
                background=super_background,
                sites=sublat_sites,
                skip_non_primitive=skip_non_primitive,
                skip_equivalents=skip_equivalents,
                use_background_invariant_group=True,
            ):
                yield config

    def by_cluster(
        self,
        background: casmconfig.Configuration,
        cluster_specs: dict,
        supercells: Optional[dict] = None,
        skip_non_primitive: bool = False,
        skip_equivalents: bool = True,
    ):
        """Enumerate occupation perturbations of a background configuration on
        specified clusters

        Parameters
        ----------
        background: casmconfig.Configuration,
            The background configuration on which enumeration takes place. This allows
            fixing other DoF and enumerating all occupations and/or specifying the
            starting occupation which is perturbed.
        cluster_specs: dict
            Parameters used to construct a
            :class:`~libcasm.clusterography.ClusterSpecs` instance using
            :func:`~libcasm.clusterography.ClusterSpecs.from_dict`. Specifies orbits of
            clusters in the infinite crystal. Enumeration is done on all symmetrically
            distinct clusters in those orbits within the background configurations
            generated by filling `background` into the supercells specified by the
            `supercells` parameter.
        supercells: Optional[dict] = None
            Parameters to forward to
            :func:`ScelEnum.by_volume <libcasm.enumerate.ScelEnum.by_volume>`
            to specify the supercells that the background configuration will be filled
            into. If None, only the exact background configuration is perturbed.
        skip_non_primitive: bool = False
            If True, enumeration skips non-primitive configurations. All DoF are
            included in the check for primitive configurations.
        skip_equivalents: bool = True
            If True, enumeration skips equivalent configurations with respect
            to the subgroup that leaves the background configuration invariant
            and does not mix the cluster sites and other sites.

        Yields
        ------
        config: casmconfig.Configuration
            A :class:`~casmconfig.Configuration`.
        """
        self._begin()
        super_backgrounds = self._make_super_backgrounds(
            background=background,
            supercells=supercells,
        )

        cspecs = casmclust.ClusterSpecs.from_dict(
            data=cluster_specs,
            xtal_prim=self.prim.xtal_prim,
            prim_factor_group=self.prim.factor_group,
            integral_site_coordinate_symgroup_rep=self.prim.integral_site_coordinate_symgroup_rep,
        )
        orbits = cspecs.make_orbits()

        for super_background in super_backgrounds:
            distinct_cluster_sites = make_distinct_cluster_sites(
                configuration=super_background,
                orbits=orbits,
            )
            for cluster_sites in distinct_cluster_sites:
                for config in self._by_site(
                    background=super_background,
                    sites=cluster_sites,
                    skip_non_primitive=skip_non_primitive,
                    skip_equivalents=skip_equivalents,
                    use_background_invariant_group=True,
                ):
                    yield config

    def by_cluster_list(
        self,
        background: casmconfig.Configuration,
        clusters: list[casmclust.Cluster],
        supercells: Optional[dict] = None,
        skip_non_primitive: bool = False,
        skip_equivalents: bool = True,
    ):
        """Enumerate occupation perturbations of a background configuration on
        an explicitly provided list of clusters

        Parameters
        ----------
        background: casmconfig.Configuration,
            The background configuration on which enumeration takes place. This allows
            fixing other DoF and enumerating all occupations and/or specifying the
            starting occupation which is perturbed.
        clusters: list[Clusters]
            Each cluster is used to generate an orbit of clusters in the infinite
            crystal. Enumeration is done on all symmetrically distinct clusters in
            those orbits within the background configurations generated by filling
            the supercells specified by the `supercells` parameter.
        supercells: Optional[dict] = None
            Parameters to forward to
            :func:`ScelEnum.by_volume <libcasm.enumerate.ScelEnum.by_volume>`
            to specify the supercells that the background configuration will be filled
            into. If None, only the exact background configuration is perturbed.
        skip_non_primitive: bool = False
            If True, enumeration skips non-primitive configurations. All DoF are
            included in the check for primitive configurations.
        skip_equivalents: bool = True
            If True, enumeration skips equivalent configurations with respect
            to the subgroup that leaves the background configuration invariant
            and does not mix the cluster sites and other sites.

        Yields
        ------
        config: casmconfig.Configuration
            A :class:`~casmconfig.Configuration`.
        """
        self._begin()
        super_backgrounds = self._make_super_backgrounds(
            background=background,
            supercells=supercells,
        )

        orbits = []
        orbit_prototypes = []
        for cluster in clusters:
            orbit = casmclust.make_periodic_orbit(
                orbit_element=cluster,
                integral_site_coordinate_symgroup_rep=self.prim.integral_site_coordinate_symgroup_rep,
            )
            if orbit[0] not in orbit_prototypes:
                orbit_prototypes.append(orbit[0])
                orbits.append(orbit)

        for super_background in super_backgrounds:
            distinct_cluster_sites = make_distinct_cluster_sites(
                configuration=super_background,
                orbits=orbits,
            )
            for cluster_sites in distinct_cluster_sites:
                for config in self._by_site(
                    background=super_background,
                    sites=cluster_sites,
                    skip_non_primitive=skip_non_primitive,
                    skip_equivalents=skip_equivalents,
                    use_background_invariant_group=True,
                ):
                    yield config
