from typing import Optional

import libcasm.configuration as casmconfig

from ._SuperConfigEnum import SuperConfigEnum


def make_distinct_super_configurations(
    motif: casmconfig.Configuration,
    supercell: casmconfig.Supercell,
    fix: str = "supercell",
    supercell_set: Optional[casmconfig.SupercellSet] = None,
) -> list[casmconfig.Configuration]:
    """
    Make configurations that fill a supercell and are equivalent with respect to the
    prim factor group, but distinct by supercell factor group operations

    Parameters
    ----------
    motif: libcasm.configuration.Configuration
        The initial configuration, with DoF values to be filled into `supercell` or
        equivalents of `supercell`.
    supercell: libcasm.configuration.Supercell
        The supercell to be filled by the motif configuration.
    fix: str = "supercell"
        The object to be fixed during the method. Can be "supercell" or "motif". If
        "supercell", the supercell is fixed and the motif is re-orientated using
        prim factor group operations. If "motif", then the motif is fixed and the
        supercell is re-orientated.
    supercell_set: Optional[libcasm.configuration.SupercellSet] = None
        If not None, generated :class:`Supercell` are constructed by
        adding in the :class:`~SupercellSet`.

    Returns
    -------
    super_configurations: list[Configuration]
        A list of distinct super configurations of the motif in the supercell or
        equivalent supercells.
    """
    if fix == "supercell":
        return casmconfig.make_distinct_super_configurations(
            motif=motif,
            supercell=supercell,
        )
    elif fix == "motif":
        super_config_enum = SuperConfigEnum(
            prim=motif.supercell.prim,
            supercell_set=supercell_set,
        )
        super_backgrounds = []
        for config in super_config_enum.by_supercell_list(
            motif=motif,
            supercells=[supercell],
        ):
            super_backgrounds.append(config)
        return super_backgrounds
    else:
        raise ValueError(
            f"Error in make_distinct_super_configurations: "
            f"`fix` must be 'supercell' or 'motif', not '{fix}'"
        )
