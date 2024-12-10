import math
import sys
from typing import Optional

import libcasm.clusterography as casmclust
import libcasm.configuration as casmconfig
import libcasm.enumerate as casmenum
import libcasm.xtal as xtal


def _format_figure(
    p,
    font_size_1="18pt",
    font_size_2="14pt",
    font_name="times",
):
    """Format the bokeh figure title, x-axis, and y-axis.

    Parameters
    ----------
    font_size_1: str
        The font size for the title and axis labels.
    font_size_2: str
        The font size for the axis major tick labels.
    font_name: str
        The font name.

    """

    p.title.text_font = font_name
    p.title.text_font_size = font_size_1

    p.xaxis.axis_label_text_font = font_name
    p.xaxis.axis_label_text_font_size = font_size_1
    p.xaxis.major_label_text_font = font_name
    p.xaxis.major_label_text_font_size = font_size_2

    p.yaxis.axis_label_text_font = font_name
    p.yaxis.axis_label_text_font_size = font_size_1
    p.yaxis.major_label_text_font = font_name
    p.yaxis.major_label_text_font_size = font_size_2


def _make_pareto(
    n_unitcells: list[int],
    voronoi_inner_radius: list[float],
    has_required_sites: list[bool],
    has_all_motif_operations: list[bool],
):
    # Find pareto (n_unitcells, dist),
    # restricted to has_required_sites==True and has_all_motif_operations==True
    pareto = []
    n_points = len(n_unitcells)
    tol = 1e-5
    for i in range(n_points):
        if not has_required_sites[i]:
            continue
        if not has_all_motif_operations[i]:
            continue

        point_set = []
        for j in range(n_points):
            if n_unitcells[i] >= n_unitcells[j]:
                point_set.append(j)
        is_max = True
        for j in point_set:
            if j == i:
                continue
            if n_unitcells[i] == n_unitcells[j] and math.isclose(
                voronoi_inner_radius[i], voronoi_inner_radius[j], abs_tol=tol
            ):
                continue
            if voronoi_inner_radius[i] < voronoi_inner_radius[j] + tol:
                is_max = False
                break
        if is_max:
            pareto.append(i)
    return pareto


def _score_point_defect_supercell(supercell: casmconfig.Supercell):
    """Return a score for ranking supercells for point defect calculations.

    The score is `(supercell_factor_group_size, -n_unitcells, voronoi_inner_radius).`
    This metric tends to prefer high symmetry, small, compact supercells as the "base"
    supercell, which can then be multiplied by a scaling factor to make a larger
    supercell that point defect calculations can be performed in.

    Parameters
    ----------
    supercell: casmconfig.Supercell
        The supercell being scored.

    Returns
    -------
    supercell_factor_group_size: int
        The supercell factor group size
    negative_n_unitcells: int
        The negative of the number of unit cells.
    voronoi_inner_radius: float
        The superlattice voronoi inner radius distance, which is the minimum
        defect periodic image distance.

    """
    scel_fg_size = len(supercell.factor_group.elements)
    n_unitcells = supercell.n_unitcells
    voronoi_inner_radius = supercell.superlattice.voronoi_inner_radius()
    return (scel_fg_size, -n_unitcells, voronoi_inner_radius)


def make_required_sites(
    phenomenal_clusters: list[casmclust.Cluster],
    local_orbits: list[list[list[casmclust.Cluster]]],
):
    """Make the required sites for a point defect calculation.

    Parameters
    ----------
    phenomenal_clusters: list[casmclust.Cluster]
        Phenomenal clusters.

    local_orbits: list[list[list[casmclust.Cluster]]]
        Local orbits, where ``local_orbits[i_equiv]`` is a list of orbits for the
        equivalent event ``i_equiv``.

    Returns
    -------
    required_sites: list[list[xtal.IntegralSiteCoordinate]]
        Sets of required sites for the point defect calculation, where
        `required_sites[i]` is the `i`-th set. For each set, it is required that the
        periodic images of the sites in the supercell are present without overlap.

    """
    required_sites = []
    if len(local_orbits) != len(phenomenal_clusters):
        raise ValueError(
            "Error in make_required_sites: "
            "len(local_orbits) != len(phenomenal_clusters)."
        )
    required_sites = []
    for i_equiv, _local_orbits in enumerate(local_orbits):
        _required_sites = list()
        for site in phenomenal_clusters[i_equiv]:
            if site not in _required_sites:
                _required_sites.append(site)
        for orbit in _local_orbits:
            for cluster in orbit:
                for site in cluster:
                    if site not in _required_sites:
                        _required_sites.append(site)
        _required_sites.sort()
        required_sites.append(_required_sites)
    return required_sites


def has_required_sites(
    required_sites: list[list[xtal.IntegralSiteCoordinate]],
    supercell: casmconfig.Supercell,
):
    """Check if the supercell has the required sites for a point defect calculation.

    Parameters
    ----------
    required_sites: list[list[xtal.IntegralSiteCoordinate]]
        Sets of required sites for the point defect calculation, where
        `required_sites[i]` is the `i`-th set. For each set, it is required that the
        periodic images of the sites in the supercell are present without overlap.
    supercell: casmconfig.Supercell
        The supercell to check.

    Returns
    -------
    result: bool
        True if the supercell has each set of required sites without overlap; else
        False.

    """
    for _required_sites in required_sites:
        site_indices = set()
        for site in _required_sites:
            site_indices.add(supercell.linear_site_index(site))
        if len(site_indices) != len(_required_sites):
            return False
    return True


def make_equivalent_supercells_with_required_operations(
    supercell: casmconfig.Supercell,
    required_operations: list[int],
    supercell_set: Optional[casmconfig.SupercellSet] = None,
) -> list[casmconfig.Supercell]:
    """Make equivalent supercells with required operations.

    Parameters
    ----------
    supercell: casmconfig.Supercell
        The initial supercell.
    required_operations: set[int]
        The indices of prim factor group operations that are required in the
        in supercell factor group.
    supercell_set: Optional[casmconfig.SupercellSet] = None
            If not None, generated :class:`~casmconfig.Supercell` are constructed by
            adding in the :class:`~casmconfig.SupercellSet`.

    Returns
    -------
    matching: list[casmconfig.Supercell]
        Equivalent supercells to `supercell` with the required operations.
    """
    candidates = casmconfig.make_equivalent_supercells(supercell)
    matching = []
    for candidate in candidates:
        operations = set(candidate.factor_group.head_group_index)
        if required_operations.issubset(operations):
            if supercell_set is not None:
                supercell_set.add(candidate)
            matching.append(candidate)
    return matching


def make_supercells_for_point_defects(
    motif: casmconfig.Configuration,
    base_max_volume: int,
    min_volume: int,
    max_volume: int,
    base_min_factor_group_size: Optional[int] = None,
    required_sites: Optional[list[list[xtal.IntegralSiteCoordinate]]] = None,
    supercell_set: Optional[casmconfig.SupercellSet] = None,
):
    """Find the best supercells for point defect calculations which can be filled by a
    particular motif configuration.

    Supercells are scored by the (1) size of their supercell factor group (larger is
    better) and (2) the minimum periodic distance between point defects (larger is
    better).

    Parameters
    ----------
    motif: casmconfig.Configuration
        The motif configuration.
    base_max_volume: int
        Score unit cells with volume up to this value.
    min_volume: int
        Return supercells for point defect calculations with at least this volume
    max_volume: int
        Return supercells for point defect calculations with at most this volume
    base_min_factor_group_size: Optional[int] = None
        If not None, ignore unit cells with factor group size less than this value.
    required_sites: Optional[list[list[xtal.IntegralSiteCoordinate]]] = None
        Sets of required sites for the point defect calculation, where
        `required_sites[i]` is the `i`-th set. If provided, it is checked if the
        periodic images of the sites in the supercell are present without overlap,
        for each set. Supercells are included in results in either case.
    supercell_set: Optional[casmconfig.SupercellSet] = None
            If not None, generated :class:`~casmconfig.Supercell` are constructed by
            adding in the :class:`~casmconfig.SupercellSet`.

    Returns
    -------
    candidate_supercells: Optional[dict]
        The best supercells, if any were found which meets the criteria; else None.
        The dictionary keys are the candidate supercells (using an equivalent which
        has all motif invariant group operations if possible) and the values hold data
        about the supercells:

        - "unitcell": The unit supercell used to generate the final supercell.
        - "unitcell_name": The name of the unit supercell.
        - "has_required_sites": True if the supercell has the required sites without
          overlap (or if there are no required sites).
        - "has_all_motif_operations": True if at least one equivalent supercell
          has all motif invariant group operations in the supercell factor group.

    """
    prim = motif.supercell.prim
    scel_enum = casmenum.ScelEnum(prim=prim)

    # Get all prim factor group operations in the motif invariant group
    required_operations = set()
    for op in casmconfig.make_invariant_subgroup(motif):
        required_operations.add(op.prim_factor_group_index())

    # Generate candidate unit cells
    candidate_unit_cells = []
    for supercell in scel_enum.by_volume(min=1, max=base_max_volume):
        score = _score_point_defect_supercell(supercell)
        if (
            base_min_factor_group_size is not None
            and score[0] < base_min_factor_group_size
        ):
            continue
        candidate_unit_cells.append((supercell, score))

    # Sort candidate unit cells by score:
    # (factor_group_size, -n_unitcells, voronoi_inner_radius)
    candidate_unit_cells.sort(key=lambda x: x[1], reverse=True)

    # Generate the candidate supercells
    candidate_supercells = {}
    for candidate, score in candidate_unit_cells:
        m = 0
        T1 = candidate.transformation_matrix_to_super
        while True:
            m += 1
            supercell = casmconfig.Supercell(prim, T1 * m)

            if supercell.n_unitcells > max_volume:
                break
            if supercell.n_unitcells < min_volume or supercell in candidate_supercells:
                continue

            # Check if motif can tile the supercell
            S = supercell.superlattice
            is_superlat, _, _ = S.is_equivalent_superlattice_of(
                lattice2=motif.supercell.superlattice,
                point_group=prim.factor_group.elements,
            )
            if not is_superlat:
                continue

            # Check if at least one equivalent supercell has all required operations
            matching = make_equivalent_supercells_with_required_operations(
                supercell=supercell,
                required_operations=required_operations,
                supercell_set=supercell_set,
            )
            has_all_motif_operations = len(matching) > 0

            if has_all_motif_operations:
                key = matching[0]
            else:
                key = supercell

            # All filters passed: include supercell
            if supercell_set is not None:
                supercell_set.add(key)
            candidate_supercells[key] = dict(
                unitcell=candidate,
                unitcell_name=casmconfig.SupercellRecord(candidate).supercell_name,
                has_required_sites=(
                    True
                    if required_sites is None
                    else has_required_sites(required_sites, supercell)
                ),
                has_all_motif_operations=has_all_motif_operations,
            )

    if len(candidate_supercells) == 0:
        candidate_supercells = None
    return candidate_supercells


def filter_point_defect_supercells(
    candidate_supercells: dict,
    min_volume: Optional[int] = None,
    max_volume: Optional[int] = None,
    require_has_all_motif_operations: bool = True,
    require_has_required_sites: bool = True,
    min_factor_group_size: Optional[int] = None,
    min_voronoi_inner_radius: Optional[float] = None,
):
    """Filter candidate supercells for point defect calculations.

    Parameters
    ----------
    candidate_supercells: dict
        The candidate supercells, as returned by `make_supercells_for_point_defects`.
    min_volume: Optional[int] = None
        If not None, ignore supercells with fewer unit cells than this value.
    max_volume: Optional[int] = None
        If not None, ignore supercells with more unit cells than this value.
    require_has_all_motif_operations: bool = True
        If not True, ignore supercells that do not have at least
        one equivalent supercell which includes all required operations in the
        supercell factor group.
    require_has_required_sites: bool = True
        If not True, ignore supercells that do not have all required sites.
    min_factor_group_size: Optional[int] = None
        If not None, ignore supercells with factor group size less than this value.
    min_voronoi_inner_radius: Optional[float] = None
        If not None, ignore supercells with minimum periodic image distance less than
        this value.

    Returns
    -------
    filtered_supercells: dict
        The supercells are for which `has_required_sites` is True and the optional
        criteria hold. The dictionary values hold data about the supercells:

        - "unitcell": The unit supercell used to generate the final supercell.
        - "unitcell_name": The name of the unit supercell.
        - "has_required_sites": True if the supercell has the required sites without
          overlap.
        - "has_all_motif_operations": True if at least one equivalent supercell
          has all motif invariant group operations in the supercell factor group.

    """
    x = {}
    for supercell, info in candidate_supercells.items():
        if require_has_required_sites and not info["has_required_sites"]:
            continue
        if require_has_all_motif_operations and not info["has_all_motif_operations"]:
            continue
        if min_volume is not None and supercell.n_unitcells < min_volume:
            continue
        if max_volume is not None and supercell.n_unitcells > max_volume:
            continue
        if (
            min_factor_group_size is not None
            and len(supercell.factor_group.elements) < min_factor_group_size
        ):
            continue
        if (
            min_voronoi_inner_radius is not None
            and supercell.superlattice.voronoi_inner_radius() < min_voronoi_inner_radius
        ):
            continue

        x[supercell] = info
    return x


def find_optimal_point_defect_supercells(
    candidate_supercells: dict,
    min_volume: Optional[int] = None,
    max_volume: Optional[int] = None,
    require_has_all_motif_operations: bool = True,
    require_has_required_sites: bool = True,
    min_factor_group_size: Optional[int] = None,
    min_voronoi_inner_radius: Optional[float] = None,
):
    """Find the Pareto-optimal supercells for point defect calculations.

    Parameters
    ----------
    candidate_supercells: dict
        The candidate supercells, as returned by `make_supercells_for_point_defects`.
    min_volume: Optional[int] = None
        If not None, ignore supercells with fewer unit cells than this value.
    max_volume: Optional[int] = None
        If not None, ignore supercells with more unit cells than this value.
    require_has_all_motif_operations: bool = True
        If not True, ignore supercells that do not have at least
        one equivalent supercell which includes all required operations in the
        supercell factor group.
    require_has_required_sites: bool = True
        If not True, ignore supercells that do not have all required sites.
    min_factor_group_size: Optional[int] = None
        If not None, ignore supercells with factor group size less than this value.
    min_voronoi_inner_radius: Optional[float] = None
        If not None, ignore supercells with minimum periodic image distance less than
        this value.

    Returns
    -------
    optimal_supercells: list[Supercell, dict]
        The Pareto-optimal supercells and info, sorted by number of unit cells. The
        Pareto-optimal supercells are found for the set of supercells for which
        `has_required_sites` is True and the optional criteria hold - they may not be
        Pareto-optimal considering all supercells. The dictionary values hold data
        about the supercells:

        - "unitcell": The unit supercell used to generate the final supercell.
        - "unitcell_name": The name of the unit supercell.
        - "has_required_sites": True if the supercell has the required sites without
          overlap.
        - "has_all_motif_operations": True if at least one equivalent supercell
          has all motif invariant group operations in the supercell factor group.


    """
    x = filter_point_defect_supercells(
        candidate_supercells=candidate_supercells,
        min_volume=min_volume,
        max_volume=max_volume,
        require_has_all_motif_operations=require_has_all_motif_operations,
        require_has_required_sites=require_has_required_sites,
        min_factor_group_size=min_factor_group_size,
        min_voronoi_inner_radius=min_voronoi_inner_radius,
    )

    if len(x) == 0:
        return []

    supercells = [supercell for supercell in x]
    values = [x[supercell] for supercell in supercells]
    n_unitcells = [supercell.n_unitcells for supercell in supercells]
    voronoi_inner_radius = [
        supercell.superlattice.voronoi_inner_radius() for supercell in supercells
    ]
    has_required_sites = [value["has_required_sites"] for value in values]
    has_all_motif_operations = [value["has_all_motif_operations"] for value in values]

    # Find pareto (n_unitcells, dist)
    pareto = _make_pareto(
        n_unitcells, voronoi_inner_radius, has_required_sites, has_all_motif_operations
    )

    return sorted(
        [(supercells[i], values[i]) for i in pareto],
        key=lambda x: -x[0].n_unitcells,
        reverse=True,
    )


def print_point_defect_supercell_info(
    supercell: casmconfig.Supercell,
    data: dict,
):
    record = casmconfig.SupercellRecord(supercell)
    base_record = casmconfig.SupercellRecord(data["unitcell"])
    print(f"Supercell name: {record.supercell_name}")
    print("Number of unit cells:", supercell.n_unitcells)
    print(
        "Minimum periodic image distance:",
        supercell.superlattice.voronoi_inner_radius(),
    )
    print("Factor group size:", len(supercell.factor_group.elements))
    print("Has all motif operations:", data["has_all_motif_operations"])
    print("Has required sites:", data["has_required_sites"])
    print("Superlattice vectors: (column vector matrix)")
    print(supercell.superlattice.column_vector_matrix())
    print("Transformation matrix to supercell:")
    print(supercell.transformation_matrix_to_super)
    print(f"Base supercell: {base_record.supercell_name}")
    print("Base supercell transformation matrix:")
    print(data["unitcell"].transformation_matrix_to_super)
    sys.stdout.flush()
    print()


def plot_point_defect_supercell_scores(
    candidate_supercells: dict,
    motif_name: str,
    min_volume: Optional[int] = None,
    max_volume: Optional[int] = None,
    require_has_all_motif_operations=True,
    require_has_required_sites=True,
    min_factor_group_size: Optional[int] = None,
    min_voronoi_inner_radius: Optional[float] = None,
    color_optimal: str = "darkorange",
    color_non_optimal: str = "steelblue",
    size_matching: int = 10,
    size_non_matching: int = 5,
    marker_optimal: str = "star",
    marker_matching: str = "circle",
    marker_non_matching: str = "x",
):
    """Use bokeh to plot the scores of supercells for point defect calculations.

    Plots the scores of supercells for point defect calculations, where:

    - x-axis: number of unit cells in the supercell
    - y-axis: minimum periodic image distance between point defects
    - size: supercells that have the required sites and match the optional criteria
      (`min_volume`, `max_volume`, `min_factor_group_size`, and
      `min_voronoi_inner_radius`) are plotted with a larger marker size than those
      which are non-matching
    - color: pareto-optimal supercells for which `has_required_sites` and
      `has_all_motif_operations` are True are highlighted

    Parameters
    ----------
    candidate_supercells: dict
        The candidate supercells, as returned by `make_supercells_for_point_defects`.
    motif_name: str
        The name of the motif configuration to use in the plot title.
    min_volume: Optional[int] = None
        If not None, ignore supercells with fewer unit cells than this value.
    max_volume: Optional[int] = None
        If not None, ignore supercells with more unit cells than this value.
    require_has_all_motif_operations: bool = True
        If not True, ignore supercells that do not have at least
        one equivalent supercell which includes all required operations in the
        supercell factor group.
    require_has_required_sites: bool = True
        If not True, ignore supercells that do not have all required sites.
    min_factor_group_size: Optional[int] = None
        If not None, ignore supercells with factor group size less than this value.
    min_voronoi_inner_radius: Optional[float] = None
        If not None, ignore supercells with minimum periodic image distance less than
        this value.
    color_optimal: str = "darkorange"
        The color for pareto-optimal supercells.
    color_non_optimal: str = "steelblue"
        The color for non-pareto-optimal supercells.
    size_matching: int = 10
        The size of the markers for supercells that match the criteria.
    size_non_matching: int = 5
        The size of the markers for supercells that do not match the criteria.
    marker_optimal: str = "star"
        The marker for pareto-optimal supercells in the matching set.
    marker_matching: str = "circle"
        The marker for supercells that are in the matching set, but not optimal.
    marker_non_matching: str = "x"
        The marker for supercells that are not in the matching set.
    """
    from bokeh.models import ColumnDataSource
    from bokeh.plotting import figure, show

    # Make a ColumnDataSource
    x = candidate_supercells
    supercells = [supercell for supercell in candidate_supercells]
    values = [x[supercell] for supercell in supercells]
    factor_group_size = [
        len(supercell.factor_group.elements) for supercell in supercells
    ]
    n_unitcells = [supercell.n_unitcells for supercell in supercells]
    voronoi_inner_radius = [
        supercell.superlattice.voronoi_inner_radius() for supercell in supercells
    ]
    has_required_sites = [value["has_required_sites"] for value in values]
    has_all_motif_operations = [value["has_all_motif_operations"] for value in values]
    supercell_name = [
        casmconfig.SupercellRecord(supercell).supercell_name
        for supercell in candidate_supercells
    ]
    T = [supercell.transformation_matrix_to_super for supercell in supercells]
    v1_frac = [_T[:, 0] for _T in T]
    v2_frac = [_T[:, 1] for _T in T]
    v3_frac = [_T[:, 2] for _T in T]
    S = [supercell.superlattice.column_vector_matrix() for supercell in supercells]
    v1_cart = [[round(x, 3) for x in _S[:, 0]] for _S in S]
    v2_cart = [[round(x, 3) for x in _S[:, 1]] for _S in S]
    v3_cart = [[round(x, 3) for x in _S[:, 2]] for _S in S]

    filtered_supercells = filter_point_defect_supercells(
        candidate_supercells=candidate_supercells,
        min_volume=min_volume,
        max_volume=max_volume,
        require_has_all_motif_operations=require_has_all_motif_operations,
        require_has_required_sites=require_has_required_sites,
        min_factor_group_size=min_factor_group_size,
        min_voronoi_inner_radius=min_voronoi_inner_radius,
    )
    # _sizes = [
    #     size_matching if supercells[i] in filtered_supercells else size_non_matching
    #     for i in range(len(supercells))
    # ]

    # Find optimal supercells: list[tuple[Supercell, dict]]
    optimal_supercells = find_optimal_point_defect_supercells(
        candidate_supercells=filtered_supercells,
    )
    # Just the supercells:
    _optimal_supercells = [x[0] for x in optimal_supercells]

    # _colors = [
    #     color_optimal if supercells[i] in _optimal_supercells else color_non_optimal
    #     for i in range(len(supercells))
    # ]

    _label = []
    _color = []
    _marker = []
    _sizes = []
    for i in range(len(supercells)):
        if supercells[i] in _optimal_supercells:
            _label.append("Optimal")
            _color.append(color_optimal)
            _marker.append(marker_optimal)
            _sizes.append(size_matching)
        elif supercells[i] in filtered_supercells:
            _label.append("Matching")
            _color.append(color_non_optimal)
            _marker.append(marker_matching)
            _sizes.append(size_matching)
        else:
            _label.append("Non-matching")
            _color.append(color_non_optimal)
            _marker.append(marker_non_matching)
            _sizes.append(size_non_matching)

    source = ColumnDataSource(
        data=dict(
            supercell_name=supercell_name,
            factor_group_size=factor_group_size,
            n_unitcells=n_unitcells,
            voronoi_inner_radius=voronoi_inner_radius,
            has_required_sites=has_required_sites,
            has_all_motif_operations=has_all_motif_operations,
            v1_frac=v1_frac,
            v2_frac=v2_frac,
            v3_frac=v3_frac,
            v1_cart=v1_cart,
            v2_cart=v2_cart,
            v3_cart=v3_cart,
            color=_color,
            sizes=_sizes,
            label=_label,
            marker=_marker,
        )
    )

    # from bokeh.io import output_notebook
    # output_notebook()

    figure_params = {
        "tools": "pan,box_zoom,wheel_zoom,save,hover,lasso_select,box_select,tap,reset"
    }
    if motif_name is not None:
        figure_params["title"] = "Supercells for " + motif_name
    figure_params["tooltips"] = [
        ("supercell", "@supercell_name"),
        ("factor group size", "@factor_group_size"),
        ("n_unitcells", "@n_unitcells"),
        ("voronoi_inner_radius", "@voronoi_inner_radius"),
        ("has_required_sites", "@has_required_sites"),
        ("has_all_motif_operations", "@has_all_motif_operations"),
        ("v1_frac", "@v1_frac"),
        ("v2_frac", "@v2_frac"),
        ("v3_frac", "@v3_frac"),
        ("v1_cart", "@v1_cart"),
        ("v2_cart", "@v2_cart"),
        ("v3_cart", "@v3_cart"),
    ]

    p = figure(**figure_params)

    p.scatter(
        x="n_unitcells",
        y="voronoi_inner_radius",
        size="sizes",
        color="color",
        marker="marker",
        source=source,
        legend_group="label",
    )

    # label axes
    p.xaxis.axis_label = "Number of unitcells"
    p.yaxis.axis_label = "Minimum periodic image distance"

    p.legend.location = "bottom_right"

    _format_figure(p)
    show(p)
