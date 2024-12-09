import typing

from libcasm.xtal import IntegralSiteCoordinate

from ._configuration import (
    Configuration,
    ConfigurationWithProperties,
    SupercellSymOp,
    apply_to_configuration,
    apply_to_configuration_with_properties,
    apply_to_integral_site_coordinate,
    copy_apply_to_configuration,
    copy_apply_to_configuration_with_properties,
    copy_apply_to_integral_site_coordinate,
    copy_transformed_configuration,
)


def apply(
    op: SupercellSymOp,
    obj: typing.Union[
        Configuration, ConfigurationWithProperties, IntegralSiteCoordinate
    ],
) -> typing.Union[Configuration, ConfigurationWithProperties, IntegralSiteCoordinate]:
    """
    Apply a SupercellSymOp to a Configuration, ConfigurationWithProperties, or
    IntegralSiteCoordinate

    Parameters
    ----------
    op: SupercellSymOp
        A symmetry operation to apply
    obj: typing.Union[Configuration, ConfigurationWithProperties, \
    libcasm.xtal.IntegralSiteCoordinate]
        The object to which the symmetry operation should be applied

    Returns
    -------
    obj: typing.Union[Configuration, ConfigurationWithProperties, \
    libcasm.xtal.IntegralSiteCoordinate]
        The object after the symmetry operation has been applied
    """
    if isinstance(obj, Configuration):
        return apply_to_configuration(op, obj)
    elif isinstance(obj, ConfigurationWithProperties):
        return apply_to_configuration_with_properties(op, obj)
    elif isinstance(obj, IntegralSiteCoordinate):
        return apply_to_integral_site_coordinate(op, obj)
    else:
        raise TypeError(
            "Error in libcasm.configuration.apply: "
            "`obj` must be a Configuration, ConfigurationWithProperties, "
            f"or IntegralSiteCoordinate, not {type(obj)}"
        )


def copy_apply(
    op: SupercellSymOp,
    obj: typing.Union[
        Configuration, ConfigurationWithProperties, IntegralSiteCoordinate
    ],
) -> typing.Union[Configuration, ConfigurationWithProperties, IntegralSiteCoordinate]:
    """
    Copy a Configuration, ConfigurationWithProperties, or IntegralSiteCoordinate,
    then apply a SupercellSymOp

    Parameters
    ----------
    op: SupercellSymOp
        A symmetry operation to apply
    obj: typing.Union[Configuration, ConfigurationWithProperties, \
    libcasm.xtal.IntegralSiteCoordinate]
        The object to be copied and then transformed

    Returns
    -------
    obj: typing.Union[Configuration, ConfigurationWithProperties, \
    libcasm.xtal.IntegralSiteCoordinate]
        The new object, after the symmetry operation has been applied
    """
    if isinstance(obj, Configuration):
        return copy_apply_to_configuration(op, obj)
    elif isinstance(obj, ConfigurationWithProperties):
        return copy_apply_to_configuration_with_properties(op, obj)
    elif isinstance(obj, IntegralSiteCoordinate):
        return copy_apply_to_integral_site_coordinate(op, obj)
    else:
        raise TypeError(
            "Error in libcasm.configuration.copy_apply: "
            "`obj` must be a Configuration, ConfigurationWithProperties, "
            f"or IntegralSiteCoordinate, not {type(obj)}"
        )


def find_mapping_operation(
    configuration: Configuration,
    configuration_ref: Configuration,
):
    """Find a symmetry operation that maps the configuration to an equivalent reference
    configuration which may be in a different supercell

    Parameters
    ----------
    configuration : libcasm.configuration.Configuration
        The configuration to map.
    configuration_ref : libcasm.configuration.Configuration
        The reference configuration to map to.

    Returns
    -------
    mapping_op : libcasm.xtal.SymOp
        The symmetry operation that maps `configuration` to `configuration_ref`.
    """
    # Find a symmetry operation that maps
    # the configuration to the reference configuration
    fg = configuration.supercell.prim.factor_group
    for i, op in enumerate(fg.elements):
        test_lat = op * configuration.supercell.superlattice
        if not test_lat.is_equivalent_to(configuration_ref.supercell.superlattice):
            continue
        config_init = copy_transformed_configuration(
            prim_factor_group_index=i,
            translation=[0, 0, 0],
            motif=configuration,
            supercell=configuration_ref.supercell,
        )
        rep_it = SupercellSymOp.begin(config_init.supercell)
        rep_end = SupercellSymOp.end(config_init.supercell)
        while rep_it != rep_end:
            test_config = rep_it * config_init
            if test_config == configuration_ref:
                return rep_it.to_symop() * op
            rep_it.next()
    return None


def make_consistent_asymmetric_unit_indices(
    initial: list[list[int]],
    configuration_init: Configuration,
    reference: list[list[int]],
    configuration_ref: Configuration,
):
    """Update asymmetric unit indices which are consistent with a reference set
    generated for an equivalent reference configuration.

    Parameters
    ----------
    initial : list[list[int]]
        The asymmetric unit indices for `configuration`.
    configuration_init : libcasm.configuration.Configuration
        The configuration for which the asymmetric unit indices are given.
    reference : list[list[int]]
        The asymmetric unit indices for `configuration_ref`.
    configuration_ref : libcasm.configuration.Configuration
        The reference configuration for which the reference asymmetric unit indices
        are given.

    Returns
    -------
    final : Optional[list[list[int]]]
        The asymmetric unit indices for `configuration`, made consistent with the
        `asymmetric_unit_indices_ref` indices for `configuration_ref`. If the
        configurations are not equivalent, then None is returned.

    """
    mapping_op = find_mapping_operation(
        configuration=configuration_init,
        configuration_ref=configuration_ref,
    )
    if mapping_op is None:
        return None

    # Find how the asymmetric units map
    xtal_prim = configuration_init.supercell.prim.xtal_prim
    f = configuration_init.supercell.site_index_converter
    f_ref = configuration_ref.supercell.site_index_converter
    final = [None] * len(reference)
    for i_asym, asym_unit in enumerate(initial):
        site = f.integral_site_coordinate(asym_unit[0])
        coord_cart = site.coordinate_cart(prim=xtal_prim)
        transformed_coord_cart = mapping_op * coord_cart
        transformed_site = IntegralSiteCoordinate.from_coordinate_cart(
            coordinate_cart=transformed_coord_cart,
            prim=xtal_prim,
        )
        transformed_site_index = f_ref.linear_site_index(transformed_site)
        i_asym_final = None
        for i_asym_ref, asym_unit_ref in enumerate(reference):
            if transformed_site_index in asym_unit_ref:
                i_asym_final = i_asym_ref
                break
        if i_asym_final is None:
            raise ValueError("Error: Failed making consistent asymmetric unit indices.")
        final[i_asym_final] = list(asym_unit)
    return final
