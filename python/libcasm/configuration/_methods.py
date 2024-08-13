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
