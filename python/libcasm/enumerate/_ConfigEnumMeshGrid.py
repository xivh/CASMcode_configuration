import copy
from typing import Optional, TypeVar, Union

import numpy as np

import libcasm.casmglobal
import libcasm.clexulator as casmclex
import libcasm.configuration as casmconfig
import libcasm.counter
from libcasm.configuration._misc import (
    equivalent_order_parameters_index,
    make_canonical_order_parameters,
)
from libcasm.irreps import (
    SubWedge,
)

array_like = TypeVar("array_like")


def _is_corner_point(
    x: np.ndarray,
    shape_factor: np.ndarray,
    abs_tol: float = libcasm.casmglobal.TOL,
) -> bool:
    if x.transpose() @ shape_factor @ x > 1.0 + abs_tol:
        return True
    return False


def _subwedge_counter(
    subwedge: SubWedge,
    stop: float,
    num: int,
) -> tuple[list[np.ndarray], libcasm.counter.IntCounter]:
    """Generate a counter for integer coordinates of subwedge points"""
    dim = subwedge.trans_mat.shape[1]
    start = [0.0] * dim
    stop = [stop] * dim
    num = [num] * dim

    # For axes with multiplicity==1, include both positive and negative coordinates;
    # otherwise, only include positive
    subwedge_axis_index = 0
    for irrep_wedge in subwedge.irrep_wedges:
        for i, m in enumerate(irrep_wedge.mult):
            if m == 1:
                start[subwedge_axis_index] = -stop[subwedge_axis_index]

                # keep point at 0.0 * axis
                num[subwedge_axis_index] = 2 * num[subwedge_axis_index] - 1
        subwedge_axis_index += 1

    xi = [
        np.linspace(_start, _stop, num=_num)
        for _start, _stop, _num in zip(start, stop, num)
    ]

    counter = libcasm.counter.IntCounter(
        initial=[0] * dim,
        final=[len(x) - 1 for x in xi],
        increment=[1] * dim,
    )

    return (xi, counter)


def meshgrid_points(
    background: casmconfig.Configuration,
    dof_space: casmclex.DoFSpace,
    xi: list[array_like],
    skip_equivalents: bool = False,
    abs_tol: float = libcasm.casmglobal.TOL,
):
    """Generate points in a meshgrid

    Parameters
    ----------
    background: casmconfig.Configuration
        The background configuration on which enumeration takes place.
    dof_space: casmclex.DoFSpace
        Specifies the DoF being enumerated and the basis of the DoFSpace are the
        axes on which the meshgrid is constructed and enumeration takes place.
        For local DoF, the dof_space supercell must tile the background supercell.
        Not supported for ``dof_space.dof_key == "occ"``, which raises.
    xi : list[array_like]
        1-D arrays, `[x1, x2, ...]`, representing the coordinates of a grid
        generated as if by `np.meshgrid`. x1 gives the coordinates along the first
        `dof_space` basis vector, x2 along the second `dof_space` vector, etc.
        There must be one array per `dof_space` basis vector.
    skip_equivalents: bool = False
        If True, skip symmetrically equivalent points, using the symmetry of
        the `background` configuration.
    abs_tol: float = :data:`~libcasm.casmglobal.TOL`
        The absolute tolerance used for checking equivalence.

    Yields
    ------
    eta: np.ndarray
        A point in the meshgrid.
    """

    if skip_equivalents:
        fg = casmconfig.make_invariant_subgroup(background)
        dof_space_rep = casmconfig.make_dof_space_rep(group=fg, dof_space=dof_space)
        canonical_eta_list = []

    dim = len(xi)
    counter = libcasm.counter.IntCounter(
        initial=[0] * dim,
        final=[len(x) - 1 for x in xi],
        increment=[1] * dim,
    )

    for indices in counter:
        eta = np.array([x[i] for x, i in zip(xi, indices)])
        if skip_equivalents:
            if (
                equivalent_order_parameters_index(
                    canonical_eta_list=canonical_eta_list,
                    eta=eta,
                    dof_space_rep=dof_space_rep,
                    abs_tol=abs_tol,
                )
                is not None
            ):
                continue
            else:
                canonical_eta_list.append(
                    make_canonical_order_parameters(
                        eta=eta, dof_space_rep=dof_space_rep, abs_tol=abs_tol
                    )
                )
        yield eta


def irreducible_wedge_points(
    background: casmconfig.Configuration,
    dof_space: casmclex.DoFSpace,
    irreducible_wedge: list[SubWedge],
    stop: float,
    num: int,
    trim_corners: bool = True,
    skip_equivalents: bool = False,
    abs_tol: float = libcasm.casmglobal.TOL,
):
    """Generate points in the irreducible wedge

    Parameters
    ----------
    background: casmconfig.Configuration
        The background configuration on which enumeration takes place.
    dof_space: casmclex.DoFSpace
        Specifies the DoF being enumerated and the basis of the DoFSpace are the
        axes on which the meshgrid is constructed and enumeration takes place.
        For local DoF, the dof_space supercell must tile the background supercell.
        Not supported for ``dof_space.dof_key == "occ"``, which raises.
    irreducible_wedge: list[SubWedge]
        The irreducible wedge, from a :class:`VectorSpaceSymReport` calculated
        using :func:`casmconfig.dof_space_analysis` of `background`.
    stop: float
        The ending value of the sequence of grid points along each wedge edge
        vector. The same value is used along each edge vector. Must be positive.

        The start value is 0.0, unless the corresponding vector has a symmetric
        multiplicity of 1, in which case the start value is ``-stop``. This treats
        cases such as :math:`e_1` strain for which positive and negative values are
        symmetrically distinct.
    num: int
        Number of grid points to generate along each wedge edge vector, using
        `numpy.linspace`. The origin and `stop` value are included in the sequence.
        Must be positive.

        If the corresponding vector has symmetric multiplicity of 1, then
        ``num*2-1`` is used along that vector to include positive and negative
        coordinates while keeping the same spacing.
    trim_corners: bool = True
        If True, skip grid points that lie outside the ellipsoid inscribed
        within the extrema of the grid.
    skip_equivalents: bool = False
        If True, skip symmetrically equivalent points, using the symmetry of
        the `background` configuration.
    abs_tol: float = :data:`~libcasm.casmglobal.TOL`
        The absolute tolerance used for checking equivalence.

    Yields
    ------
    (subwedge_index, eta): Tuple[int, np.ndarray]

        subwedge_index: int
            The index of the current SubWedge in the irreducible wedge.

        eta: np.ndarray
            A point in the irreducible wedge, as coordinates with respect to
            `dof_space.basis`.
    """

    if skip_equivalents:
        fg = casmconfig.make_invariant_subgroup(background)
        dof_space_rep = casmconfig.make_dof_space_rep(group=fg, dof_space=dof_space)
        canonical_eta_list = []
    if trim_corners:
        dim = dof_space.subspace_dim
        shape_factor = np.eye(dim)
        for i in range(dim):
            shape_factor[i, i] = 1.0 / (abs(stop) ** 2)

    for subwedge_index, subwedge in enumerate(irreducible_wedge):
        xi, counter = _subwedge_counter(subwedge, stop, num)

        for indices in counter:
            eta_subwedge = np.array([x[i] for x, i in zip(xi, indices)])
            if trim_corners:
                if _is_corner_point(
                    x=eta_subwedge, shape_factor=shape_factor, abs_tol=abs_tol
                ):
                    continue
            eta = dof_space.basis_inv @ subwedge.trans_mat @ eta_subwedge

            if skip_equivalents:
                if (
                    equivalent_order_parameters_index(
                        canonical_eta_list=canonical_eta_list,
                        eta=eta,
                        dof_space_rep=dof_space_rep,
                        abs_tol=abs_tol,
                    )
                    is not None
                ):
                    continue
                else:
                    canonical_eta_list.append(
                        make_canonical_order_parameters(
                            eta=eta, dof_space_rep=dof_space_rep, abs_tol=abs_tol
                        )
                    )
            yield (subwedge_index, eta)


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
        self._dof_space = None
        self._enum_index = None
        self._order_parameters = None
        self._subwedge_index = None

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
        configuration in which enumeration is performed."""
        return self._background

    @property
    def dof_space(self) -> Optional[casmclex.DoFSpace]:
        """During enumeration, `dof_space` is set to the current DoFspace in which
        enumeration is performed."""
        return self._dof_space

    @property
    def enum_index(self) -> Optional[int]:
        """During enumeration, `enum_index` is incremented to count over the
        combinations of background configuration and dof space on which enumeration is
        performed, starting from 0."""
        return self._enum_index

    @property
    def order_parameters(self) -> Optional[np.ndarray]:
        """During enumeration, `order_parameters` is set to the current order parameter
        values use to construct the generated configuration."""
        return self._order_parameters

    @property
    def subwedge_index(self) -> Optional[int]:
        """During enumeration by irreducible wedge, `subwedge_index` is set to the
        index of the current SubWedge of the irreducible wedge."""
        return self._subwedge_index

    def _by_grid_coordinates(
        self,
        background: casmconfig.Configuration,
        dof_space: casmclex.DoFSpace,
        xi: list[array_like],
        skip_equivalents: bool = False,
        abs_tol: float = libcasm.casmglobal.TOL,
    ):
        """Enumerate on a meshgrid, specified with coordinates along each axis of a \
        DoFSpace basis

        Parameters
        ----------
        background: casmconfig.Configuration
            The background configuration on which enumeration takes place.
        dof_space: casmclex.DoFSpace
            Specifies the DoF being enumerated and the basis of the DoFSpace are the
            axes on which the meshgrid is constructed and enumeration takes place.
            For local DoF, the dof_space supercell must tile the background supercell.
            Not supported for ``dof_space.dof_key == "occ"``, which raises.
        xi : list[array_like]
            1-D arrays, `[x1, x2, ...]`, representing the coordinates of a grid
            generated as if by `np.meshgrid`. x1 gives the coordinates along the first
            `dof_space` basis vector, x2 along the second `dof_space` vector, etc.
            There must be one array per `dof_space` basis vector.
        skip_equivalents: bool = False
            If True, skip symmetrically equivalent configurations.
        abs_tol: float = :data:`~libcasm.casmglobal.TOL`
            The absolute tolerance used for checking equivalence.

        Yields
        ------
        config: casmconfig.Configuration
            A :class:`~casmconfig.Configuration` on a point in the DoF values meshgrid.
        """
        self._background = background
        self._dof_space = dof_space
        if self._enum_index is None:
            self._enum_index = 0
        else:
            self._enum_index += 1

        if self.supercell_set is not None:
            self.supercell_set.add_supercell(self.background.supercell)

        if dof_space.dof_key == "occ":
            raise Exception(
                "Error in ConfigEnumMeshGrid._by_grid_coordinates: "
                'Invalid dof_space, dof_space.dof_key == "occ" is not supported'
            )

        # Note: Do skip_equivalents work here, using
        # casmconfig.make_canonical_configuration, instead of generating points without
        # equivalents using meshgrid_points(..., skip_equivalents=True) because:
        # 1) This is exact: does not depend on background choice
        # 2) This seems to be faster: meshgrid_points currently uses comparisons of the
        #    order parameter vectors in Python, vs Configuration comparisons in C++

        if skip_equivalents:
            is_canonical_background_supercell = casmconfig.is_canonical_supercell(
                background.supercell
            )
            if is_canonical_background_supercell:
                canonical_configs = casmconfig.ConfigurationSet()
            else:
                canonical_configs = []

        for eta in meshgrid_points(
            background=background,
            dof_space=dof_space,
            xi=xi,
            skip_equivalents=False,  # see note above
            abs_tol=abs_tol,
        ):
            config = copy.copy(self.background)
            config.set_order_parameters(
                dof_space=dof_space,
                order_parameters=eta,
            )
            if skip_equivalents:
                canonical_config = casmconfig.make_canonical_configuration(config)
                if canonical_config in canonical_configs:
                    continue
                else:
                    if is_canonical_background_supercell:
                        canonical_configs.add(canonical_config)
                    else:
                        canonical_configs.append(canonical_config)
            self._order_parameters = eta
            yield config

    def by_grid_coordinates(
        self,
        background: casmconfig.Configuration,
        dof_space: casmclex.DoFSpace,
        xi: list[array_like],
        skip_equivalents: bool = False,
        abs_tol: float = libcasm.casmglobal.TOL,
    ):
        """Enumerate on a meshgrid, specified with coordinates along each axis of a \
        DoFSpace basis

        Parameters
        ----------
        background: casmconfig.Configuration
            The background configuration on which enumeration takes place.
        dof_space: casmclex.DoFSpace
            Specifies the DoF being enumerated and the basis of the DoFSpace are the
            axes on which the meshgrid is constructed and enumeration takes place.
            For local DoF, the dof_space supercell must tile the background supercell.
            Not supported for ``dof_space.dof_key == "occ"``, which raises.
        xi : list[array_like]
            1-D arrays, `[x1, x2, ...]`, representing the coordinates of a grid
            generated as if by `np.meshgrid`. x1 gives the coordinates along the first
            `dof_space` basis vector, x2 along the second `dof_space` vector, etc.
            There must be one array per `dof_space` basis vector.
        skip_equivalents: bool = False
            If True, skip symmetrically equivalent configurations.
        abs_tol: float = :data:`~libcasm.casmglobal.TOL`
            The absolute tolerance used for checking equivalence.

        Yields
        ------
        config: casmconfig.Configuration
            A :class:`~casmconfig.Configuration` on a point in the DoF values meshgrid.
        """
        # these will be set by self._by_grid_coordinates
        self._background = None
        self._dof_space = None
        self._enum_index = None

        if dof_space.dof_key == "occ":
            raise Exception(
                "Error in ConfigEnumMeshGrid.by_grid_coordinates: "
                'Invalid dof_space, dof_space.dof_key == "occ" is not supported'
            )

        for config in self._by_grid_coordinates(
            background=background,
            dof_space=dof_space,
            xi=xi,
            skip_equivalents=skip_equivalents,
            abs_tol=abs_tol,
        ):
            yield config

    def by_range(
        self,
        background: casmconfig.Configuration,
        dof_space: casmclex.DoFSpace,
        start: Union[float, list[float]],
        stop: Union[float, list[float]],
        num: Union[int, list[int], None] = None,
        step: Union[float, list[float], None] = None,
        skip_equivalents: bool = False,
        abs_tol: float = libcasm.casmglobal.TOL,
    ):
        """Enumerate on a meshgrid, with coordinates specified by a range along each \
        axis of a DoFSpace basis

        Parameters
        ----------
        background: casmconfig.Configuration
            The background configuration on which enumeration takes place.
        dof_space: casmclex.DoFSpace
            Specifies the DoF being enumerated and the basis of the DoFSpace are the
            axes on which the meshgrid is constructed and enumeration takes place.
            For local DoF, the dof_space supercell must tile the background supercell.
            Not supported for ``dof_space.dof_key == "occ"``, which raises.
        start: Union[float, list[float]]
            The starting value of the sequence of grid points along each `dof_space`
            basis vector. If scalar-valued, the same value is used along each
            `dof_space` basis vector.
        stop: Union[float, list[float]]
            The ending value of the sequence of grid points along each `dof_space`
            basis vector. If scalar-valued, the same value is used along each
            `dof_space` basis vector.
        num: Union[int, list[int], None] = None
            Number of grid points to generate along each `dof_space`
            basis vector, using `numpy.linspace`. If scalar-valued, the same value is
            used along each `dof_space` basis vector. The `stop` value is included in
            the sequence. Must be non-negative. One and only one of `num` or `step`
            must be provided.
        step: Union[float, list[float], None]=None
            Spacing between grid points along each `dof_space` basis vector, in a
            sequence generated using `numpy.arange`. If scalar-valued, the same value is
            used along each `dof_space` basis vector. The sequence does not include the
            `stop` value, except in some cases where `step` is not an integer and
            floating point round-off affects the length of the sequence. Must be
            non-negative. One and only one of `num` or `step` must be provided.
        skip_equivalents: bool = False
            If True, skip symmetrically equivalent points.
        abs_tol: float = :data:`~libcasm.casmglobal.TOL`
            The absolute tolerance used for checking equivalence.

        Yields
        ------
        config: casmconfig.Configuration
            A :class:`~casmconfig.Configuration` on a point in the DoF values meshgrid.
        """
        # these will be set by self._by_grid_coordinates
        self._background = None
        self._dof_space = None
        self._enum_index = None

        if dof_space.dof_key == "occ":
            raise Exception(
                "Error in ConfigEnumMeshGrid.by_range: "
                'Invalid dof_space, dof_space.dof_key == "occ" is not supported'
            )

        dim = dof_space.subspace_dim
        if isinstance(start, float):
            start = [start] * dim
        if isinstance(stop, float):
            stop = [stop] * dim
        if step is not None and num is None:
            if isinstance(step, float):
                step = [step] * dim
            xi = [
                np.arange(_start, _stop, _step)
                for _start, _stop, _step in zip(start, stop, step)
            ]
        elif num is not None and step is None:
            if isinstance(num, int):
                num = [num] * dim
            xi = [
                np.linspace(_start, _stop, num=_num)
                for _start, _stop, _num in zip(start, stop, num)
            ]
        else:
            raise Exception(
                "Error in ConfigEnumMeshGrid.by_range: "
                "One and only one of 'step' or 'num' is required"
            )

        for config in self._by_grid_coordinates(
            background=background,
            dof_space=dof_space,
            xi=xi,
            skip_equivalents=skip_equivalents,
            abs_tol=abs_tol,
        ):
            yield config

    def by_irreducible_wedge(
        self,
        background: casmconfig.Configuration,
        dof_space: casmclex.DoFSpace,
        irreducible_wedge: list[SubWedge],
        stop: float,
        num: int,
        trim_corners: bool = True,
        skip_equivalents: bool = False,
        abs_tol: float = libcasm.casmglobal.TOL,
    ):
        """Enumerate on a meshgrid in each subwedge of the irreducible wedge of a \
        DoFSpace

        Parameters
        ----------
        background: casmconfig.Configuration
            The background configuration on which enumeration takes place.
        dof_space: casmclex.DoFSpace
            Specifies the DoF being enumerated and the basis of the DoFSpace are the
            axes on which the meshgrid is constructed and enumeration takes place.
            For local DoF, the dof_space supercell must tile the background supercell.
            Not supported for ``dof_space.dof_key == "occ"``, which raises.
        irreducible_wedge: list[SubWedge]
            The irreducible wedge, from a :class:`VectorSpaceSymReport` calculated
            using :func:`casmconfig.dof_space_analysis` of `background`.
        stop: float
            The ending value of the sequence of grid points along each wedge edge
            vector. The same value is used along each edge vector. Must be positive.

            The start value is 0.0, unless the corresponding vector has a symmetric
            multiplicity of 1, in which case the start value is ``-stop``. This treats
            cases such as :math:`e_1` strain for which positive and negative values are
            symmetrically distinct.
        num: int
            Number of grid points to generate along each wedge edge vector, using
            `numpy.linspace`. The origin and `stop` value are included in the sequence.
            Must be positive.

            If the corresponding vector has symmetric multiplicity of 1, then
            ``num*2-1`` is used along that vector to include positive and negative
            coordinates while keeping the same spacing.
        skip_equivalents: bool = False
            If True, skip symmetrically equivalent configurations.
        trim_corners: bool = True
            If True, skip grid points that lie outside the ellipsoid inscribed
            within the extrema of the grid.

        Yields
        ------
        config: casmconfig.Configuration
            A :class:`~casmconfig.Configuration` on a point in the irreducible wedge.
        """
        self._background = background
        self._dof_space = dof_space
        self._enum_index = 0

        if self.supercell_set is not None:
            self.supercell_set.add_supercell(self.background.supercell)

        if dof_space.dof_key == "occ":
            raise Exception(
                "Error in ConfigEnumMeshGrid.grid_coordinates: "
                'Invalid dof_space, dof_space.dof_key == "occ" is not supported'
            )

        # Note: Do skip_equivalents work here, using
        # casmconfig.make_canonical_configuration, instead of generating points without
        # equivalents using meshgrid_points(..., skip_equivalents=True) because:
        # 1) This is exact: does not depend on background choice
        # 2) This seems to be faster: meshgrid_points currently uses comparisons of the
        #    order parameter vectors in Python, vs Configuration comparisons in C++

        if skip_equivalents:
            is_canonical_background_supercell = casmconfig.is_canonical_supercell(
                background.supercell
            )
            if is_canonical_background_supercell:
                canonical_configs = casmconfig.ConfigurationSet()
            else:
                canonical_configs = []

        for subwedge_index, eta in irreducible_wedge_points(
            background=background,
            dof_space=dof_space,
            irreducible_wedge=irreducible_wedge,
            stop=stop,
            num=num,
            trim_corners=trim_corners,
            skip_equivalents=False,  # see note above
            abs_tol=abs_tol,
        ):
            config = copy.copy(self.background)
            config.set_order_parameters(
                dof_space=dof_space,
                order_parameters=eta,
            )
            if skip_equivalents:
                canonical_config = casmconfig.make_canonical_configuration(config)
                if canonical_config in canonical_configs:
                    continue
                else:
                    if is_canonical_background_supercell:
                        canonical_configs.add(canonical_config)
                    else:
                        canonical_configs.append(canonical_config)
            self._subwedge_index = subwedge_index
            self._order_parameters = eta
            yield config
