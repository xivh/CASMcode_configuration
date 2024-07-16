from __future__ import annotations

import numpy as np

import libcasm.configuration as casmconfig
import libcasm.configuration.io.tools as io_tools
import libcasm.sym_info as sym_info
import libcasm.xtal as xtal

from .spglib import asdict as spg_asdict


def symgroup_to_dict_with_group_classification(
    obj: io_tools.CasmObjectWithLattice,
    symgroup: sym_info.SymGroup,
) -> dict:
    """Use spglib to add group classification data for a lattice, structure,
    configuration, or prim to SymGroup data

    Notes
    -----

    - For symmetry determination of structures or configurations, the spglib `numbers`
      are set to a unique int for each atom type.
    - For symmetry determination of prim, the spglib `numbers` are set to the
      asymmetric unit indices for each site, which do consider magnetic spin degrees
      of freedom (DoF).
    - For prim with magnetic DoF, the magmoms are set to value ``0.`` (collinear)
      or ``[0., 0., 0.]`` (non-collinear).

    Parameters
    ----------
    obj: :py:data:`~libcasm.configuration.io.spglib.CasmObjectWithPositions`
        A lattice, or an atomic structure, configuration, or prim. Configurations are
        converted to structures using ``obj.to_structure()``. Lattice are converted
        to a `prim` with a single basis site at the origin.
    symgroup: libcasm.sym_info.SymGroup
        The SymGroup

    Returns
    -------
    data: dict
        The ``SymGroup.to_dict()`` output with the following added to the
        ``data["group_classification"]`` attribute:

        - ``"spacegroup_type"``: Space group type information from spglib, based on the
          lattice, positions, and numbers.

        - ``"spacegroup_type_from_casm_symmetry"``: Space group type information
          from spglib, based on the symmetry operations found by CASM. Only added if no
          magnetic spin properties or DoF and it differs from `"spacegroup_type"`.

        - ``"magnetic_spacegroup_type"``: Space group type information
          from spglib, based on the lattice, positions, numbers, and magmoms. Only
          added if there are magnetic spin properties or DoF.

        - ``"magnetic_spacegroup_type_from_casm_symmetry"``: Space group type
          information from spglib, based on the symmetry operations found by CASM. Only
          added if there are magnetic spin properties or DoF and it differs from
          `"magnetic_spacegroup_type"`. Set to `None` for ``spglib<2.5.0``, because the
          method used is not available.

        If spglib is not available, ``data["group_classification"]`` is set to `None`.

    """

    obj = io_tools.as_structure_if_config(obj)

    if isinstance(obj, xtal.Lattice):
        lattice = obj
        obj = xtal.Prim(
            lattice=lattice,
            coordinate_frac=np.zeros((3, 1)),
            occ_dof=[["A"]],
        )
    elif isinstance(obj, (xtal.Prim, xtal.Structure)):
        lattice = obj.lattice()
    elif isinstance(obj, casmconfig.Prim):
        lattice = obj.xtal_prim.lattice()
    data = symgroup.to_dict(lattice=lattice)

    try:
        import spglib

    except ImportError:
        data["group_classification"] = None
        return data

    from .spglib import (
        as_cell,
        get_magnetic_spacegroup_type_from_symmetry,
        get_magnetic_symmetry_dataset,
        get_spacegroup_type_from_symmetry,
        get_symmetry_dataset,
    )

    cell = as_cell(obj)
    is_magnetic = False
    if len(cell) == 4:
        is_magnetic = True

    symmetry_dataset = spg_asdict(get_symmetry_dataset(obj))
    spacegroup_type = None
    if symmetry_dataset is not None:
        spacegroup_type = spg_asdict(
            spglib.get_spacegroup_type(symmetry_dataset["hall_number"])
        )
    grpcls = data["group_classification"]
    grpcls["spacegroup_type"] = spacegroup_type

    if not is_magnetic:
        from_casm = spg_asdict(
            get_spacegroup_type_from_symmetry(symgroup.elements, lattice)
        )
        if from_casm != spacegroup_type:
            grpcls["spacegroup_type_from_casm_symmetry"] = from_casm

    if is_magnetic:
        mag_symmetry_dataset = spg_asdict(get_magnetic_symmetry_dataset(obj))
        mag_spacegroup_type = None
        if mag_symmetry_dataset is not None:
            mag_spacegroup_type = spg_asdict(
                spglib.get_magnetic_spacegroup_type(mag_symmetry_dataset["uni_number"])
            )
        grpcls["magnetic_spacegroup_type"] = mag_spacegroup_type

        from_casm = spg_asdict(
            get_magnetic_spacegroup_type_from_symmetry(symgroup.elements, lattice)
        )
        if from_casm != mag_spacegroup_type:
            grpcls["magnetic_spacegroup_type_from_casm_symmetry"] = from_casm

    return data
