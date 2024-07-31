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

    def _by_site(
        self,
        background: casmconfig.Configuration,
        sites: set[int],
        skip_non_primitive: bool,
        skip_non_canonical: bool,
        use_background_invariant_group: bool,
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
        skip_non_canonical: bool
            If True, enumeration skips non-canonical configurations with respect
            to either (i) the complete set of symmetry operations that leave the
            supercell lattice vectors invariant, or (ii) the subgroup that leaves the
            background configuration invariant and does not mix the given sites and
            other sites.
        use_background_invariant_group: bool
            If True, use the subgroup that leaves the background configuration invariant
            and does not mix the given sites and other sites; otherwise use the complete
            set of symmetry operations that leave the supercell lattice vectors
            invariant.

        Yields
        ------
        config: casmconfig.Configuration
            A :class:`~casmconfig.Configuration`.
        """
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
        if skip_non_canonical:
            if use_background_invariant_group:
                canonicalization_group = casmconfig.make_invariant_subgroup(
                    configuration=background,
                    site_indices=sites,
                )
            else:
                canonicalization_group = None
        while config_enum.is_valid():
            if skip_non_primitive and not casmconfig.is_primitive_configuration(
                configuration=config_enum.value()
            ):
                config_enum.advance()
                continue
            if skip_non_canonical and not casmconfig.is_canonical_configuration(
                configuration=config_enum.value(),
                subgroup=canonicalization_group,
            ):
                config_enum.advance()
                continue
            yield config_enum.value()
            config_enum.advance()

    def by_supercell(
        self,
        supercells: dict,
        motif: Optional[casmconfig.Configuration] = None,
        skip_non_primitive: bool = True,
        skip_non_canonical: bool = True,
    ):
        """Enumerate all occupations in a series of enumerated supercells

        Parameters
        ----------
        supercells: dict
            Parameters to forward to :class:`~libcasm.configuration.ScelEnum`, to
            specify the supercells that the motif configruation will be filled into.
        motif: Optional[casmconfig.Configuration] = None
            The background configuration on which enumeration takes place. The motif is
            filled into each supercell using
            :func:`~libcasm.configuration.make_distinct_super_configurations`. If the
            motif does not tile exactly into a supercell that supercell is skipped.
            Providing a motif allows enumerating all occupations with other degrees of
            freedom (DoF) fixed. If None, the default configuration in the volume 1
            supercell is used.
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
        motif = self._set_motif(motif)
        scel_enum = ScelEnum(
            prim=self.prim,
            supercell_set=self.supercell_set,
        )
        for supercell in scel_enum.by_volume(**supercells):
            sites = set(range(supercell.n_sites))
            super_configurations = casmconfig.make_distinct_super_configurations(
                motif=motif, supercell=supercell
            )
            for background in super_configurations:
                for config in self._by_site(
                    background=background,
                    sites=sites,
                    skip_non_primitive=skip_non_primitive,
                    skip_non_canonical=skip_non_canonical,
                    use_background_invariant_group=False,
                ):
                    yield config

    def by_supercell_list(
        self,
        supercells: list[casmconfig.Supercell],
        motif: Optional[casmconfig.Configuration] = None,
        skip_non_primitive: bool = True,
        skip_non_canonical: bool = True,
    ):
        """Enumerate all occupations in a list of supercells explicitly provided

        Parameters
        ----------
        supercells: list[casmconfig.Supercell]
            An explicit list of supercells in which to perform enumeration.
        motif: Optional[casmconfig.Configuration] = None
            The background configuration on which enumeration takes place. The motif is
            filled into each supercell using
            :func:`~libcasm.configuration.make_distinct_super_configurations`. If the
            motif does not tile exactly into a supercell that supercell is skipped.
            Providing a motif allows enumerating all occupations with other degrees of
            freedom (DoF) fixed. If None, the default configuration in the volume 1
            supercell is used.
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
        motif = self._set_motif(motif)
        for supercell in supercells:
            sites = set(range(supercell.n_sites))
            super_configurations = casmconfig.make_distinct_super_configurations(
                motif=motif, supercell=supercell
            )
            for background in super_configurations:
                for config in self._by_site(
                    background=background,
                    sites=sites,
                    skip_non_primitive=skip_non_primitive,
                    skip_non_canonical=skip_non_canonical,
                    use_background_invariant_group=False,
                ):
                    yield config

    def by_linear_site_indices(
        self,
        background: casmconfig.Configuration,
        sites: set[int],
        skip_non_primitive: bool = False,
        skip_non_canonical: bool = False,
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
        skip_non_canonical: bool = False
            If True, enumeration skips non-canonical configurations with respect
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
            skip_non_canonical=skip_non_canonical,
            use_background_invariant_group=True,
        ):
            yield config

    def by_integral_site_coordinates(
        self,
        background: casmconfig.Configuration,
        sites: list[xtal.IntegralSiteCoordinate],
        skip_non_primitive: bool = False,
        skip_non_canonical: bool = False,
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
        skip_non_canonical: bool = False
            If True, enumeration skips non-canonical configurations with respect
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
            skip_non_canonical=skip_non_canonical,
            use_background_invariant_group=True,
        ):
            yield config

    def by_sublattice(
        self,
        background: casmconfig.Configuration,
        sublats: set[int],
        supercells: Optional[dict] = None,
        skip_non_primitive: bool = True,
        skip_non_canonical: bool = True,
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
            Parameters to forward to ScelEnum, to specify the supercells that the
            background configuration will be filled into. If None, only the
            exact background configuration is perturbed.
        skip_non_primitive: bool = False
            If True, enumeration skips non-primitive configurations. All DoF are
            included in the check for primitive configurations.
        skip_non_canonical: bool = False
            If True, enumeration skips non-canonical configurations with respect
            to the subgroup that leaves the background configuration invariant
            and does not mix the selected sublattice sites and other sites.

        Yields
        ------
        config: casmconfig.Configuration
            A :class:`~casmconfig.Configuration`.
        """
        self._begin()
        supercell_list = self._make_supercell_list(
            background=background,
            supercells=supercells,
        )

        for supercell in supercell_list:
            sublat_sites = _make_sublat_sites(supercell, sublats)
            super_backgrounds = casmconfig.make_distinct_super_configurations(
                motif=background, supercell=supercell
            )
            for super_background in super_backgrounds:
                for config in self._by_site(
                    background=super_background,
                    sites=sublat_sites,
                    skip_non_primitive=skip_non_primitive,
                    skip_non_canonical=skip_non_canonical,
                    use_background_invariant_group=True,
                ):
                    yield config

    def by_cluster(
        self,
        background: casmconfig.Configuration,
        cluster_specs: dict,
        supercells: Optional[dict] = None,
        skip_non_primitive: bool = False,
        skip_non_canonical: bool = False,
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
            Parameters to forward to ScelEnum, to specify the supercells that the
            background configuration will be filled into. If None, only the exact
            background configuration is perturbed.
        skip_non_primitive: bool = False
            If True, enumeration skips non-primitive configurations. All DoF are
            included in the check for primitive configurations.
        skip_non_canonical: bool = False
            If True, enumeration skips non-canonical configurations with respect
            to the subgroup that leaves the background configuration invariant
            and does not mix the cluster sites and other sites.

        Yields
        ------
        config: casmconfig.Configuration
            A :class:`~casmconfig.Configuration`.
        """
        self._begin()
        supercell_list = self._make_supercell_list(
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

        for supercell in supercell_list:
            super_backgrounds = casmconfig.make_distinct_super_configurations(
                motif=background, supercell=supercell
            )
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
                        skip_non_canonical=skip_non_canonical,
                        use_background_invariant_group=True,
                    ):
                        yield config

    def by_cluster_list(
        self,
        background: casmconfig.Configuration,
        clusters: list[casmclust.Cluster],
        supercells: Optional[dict] = None,
        skip_non_primitive: bool = False,
        skip_non_canonical: bool = False,
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
            Parameters to forward to ScelEnum, to specify the supercells that the
            background configuration will be filled into. If None, only the exact
            background configuration is perturbed.
        skip_non_primitive: bool = False
            If True, enumeration skips non-primitive configurations. All DoF are
            included in the check for primitive configurations.
        skip_non_canonical: bool = False
            If True, enumeration skips non-canonical configurations with respect
            to the subgroup that leaves the background configuration invariant
            and does not mix the cluster sites and other sites.

        Yields
        ------
        config: casmconfig.Configuration
            A :class:`~casmconfig.Configuration`.
        """
        self._begin()
        supercell_list = self._make_supercell_list(
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

        for supercell in supercell_list:
            super_backgrounds = casmconfig.make_distinct_super_configurations(
                motif=background, supercell=supercell
            )
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
                        skip_non_canonical=skip_non_canonical,
                        use_background_invariant_group=True,
                    ):
                        yield config
