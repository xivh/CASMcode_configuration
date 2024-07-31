import numpy as np
import pytest

import libcasm.configuration as config
import libcasm.configuration.io as config_io
import libcasm.configuration.io.spglib as spglib_io
import libcasm.xtal as xtal

from .functions import (
    check_magnetic_symmetry_dataset,
    check_symmetry_dataset,
    make_discrete_magnetic_atom,
)


def test_FCC_binary_discrete_Cmagspin_prim_1(FCC_binary_discrete_Cmagspin_prim):
    xtal_prim = FCC_binary_discrete_Cmagspin_prim
    factor_group = xtal.make_factor_group(xtal_prim)
    assert len(factor_group) == 96
    prim = config.Prim(xtal_prim)
    assert prim.xtal_prim.coordinate_frac().shape == (3, 1)
    assert len(prim.factor_group.elements) == 96

    assert prim.is_atomic
    assert prim.continuous_magspin_key is None
    assert prim.continuous_magspin_flavor is None
    assert prim.discrete_atomic_magspin_key == "Cmagspin"
    assert prim.discrete_atomic_magspin_flavor == "C"

    lattice = prim.xtal_prim.lattice()
    for op in prim.factor_group.elements:
        syminfo = xtal.SymInfo(op, lattice)
        print(xtal.pretty_json(syminfo.to_dict()))


def test_prim_magspin_spglib_io(FCC_binary_discrete_Cmagspin_prim):
    xtal_prim = FCC_binary_discrete_Cmagspin_prim
    prim = config.Prim(xtal_prim)
    cell = spglib_io.as_cell(prim)

    assert isinstance(cell, tuple)
    assert len(cell) == 4

    symmetry_dataset = spglib_io.get_symmetry_dataset(prim)
    check_symmetry_dataset(
        symmetry_dataset,
        number=225,
        n_rotations=48,
        hall_number=523,
    )

    mag_symmetry_dataset = spglib_io.get_magnetic_symmetry_dataset(prim)
    check_magnetic_symmetry_dataset(
        mag_symmetry_dataset,
        uni_number=1619,
        n_operations=96,
        hall_number=523,
    )

    # print(prim.factor_group.brief_cart(xtal_prim.lattice()))

    with pytest.raises(
        ValueError, match="not valid for operations with time reversal symmetry"
    ):
        _ = spglib_io.get_spacegroup_type_from_symmetry(
            elements=prim.factor_group.elements,
            lattice=xtal_prim.lattice(),
        )


def test_symgroup_to_dict_with_group_classification(FCC_binary_discrete_Cmagspin_prim):
    xtal_prim = FCC_binary_discrete_Cmagspin_prim
    prim = config.Prim(xtal_prim)

    data = config_io.symgroup_to_dict_with_group_classification(prim, prim.factor_group)

    assert isinstance(data, dict)
    assert "spacegroup_type" in data["group_classification"]
    assert data["group_classification"]["spacegroup_type"]["number"] == 225
    assert "magnetic_spacegroup_type" in data["group_classification"]
    assert (
        data["group_classification"]["magnetic_spacegroup_type"]["uni_number"] == 1619
    )
    assert "spacegroup_type_from_casm_symmetry" not in data["group_classification"]

    # for spglib<2.5.0 this returns None because the method is not available;
    # for spglib>=2.5.0 this should not exist (because the same magnetic spacegroup
    # type should be found by spglib using its own symmetry operations or the ones
    # found by CASM and therefore both are not included)
    if "magnetic_spacegroup_type_from_casm_symmetry" in data["group_classification"]:
        assert (
            data["group_classification"]["magnetic_spacegroup_type_from_casm_symmetry"]
            is None
        )


def test_magspin_occ_symgroup_rep():
    # Lattice vectors
    lattice = xtal.Lattice(np.eye(3))

    # Basis sites positions, as columns of a matrix,
    # in fractional coordinates with respect to the lattice vectors
    coordinate_frac = np.array(
        [
            [0.0, 0.0, 0.0],
            [0.0, 0.5, 0.5],
            [0.5, 0.0, 0.5],
            [0.5, 0.5, 0.0],
        ]
    ).transpose()

    # Occupation degrees of freedom (DoF)
    occupants = {
        "A.up": make_discrete_magnetic_atom(name="A", value=1, flavor="C"),
        "A.down": make_discrete_magnetic_atom(name="A", value=-1, flavor="C"),
        "B.up": make_discrete_magnetic_atom(name="B", value=1, flavor="C"),
        "B.down": make_discrete_magnetic_atom(name="B", value=-1, flavor="C"),
    }
    occ_dof = [
        ["A.up", "A.down"],
        ["A.up", "A.down", "B.up", "B.down"],
        ["A.up", "B.up", "B.down", "A.down"],
        ["A.up", "A.down", "B.up", "B.down"],
    ]

    xtal_prim = xtal.Prim(
        lattice=lattice,
        coordinate_frac=coordinate_frac,
        occ_dof=occ_dof,
        occupants=occupants,
    )
    prim = config.Prim(xtal_prim)

    assert len(prim.factor_group.elements) == 96

    occ_symgroup_rep = prim.occ_symgroup_rep

    up_down_mappings = 0

    for i_factor_group, occ_op_rep in enumerate(occ_symgroup_rep):
        site_rep = prim.integral_site_coordinate_symgroup_rep[i_factor_group]
        for i_sublat_before, occ_sublat_rep in enumerate(occ_op_rep):
            site_before = xtal.IntegralSiteCoordinate(i_sublat_before, [0, 0, 0])
            site_after = site_rep * site_before
            i_sublat_after = site_after.sublattice()
            for i_occ_before in range(len(occ_sublat_rep)):
                i_occ_after = occ_sublat_rep[i_occ_before]

                orientation_name_before = occ_dof[i_sublat_before][i_occ_before]
                orientation_name_after = occ_dof[i_sublat_after][i_occ_after]

                # assert occupants map (chemical name match)
                assert (
                    occupants[orientation_name_before].name()
                    == occupants[orientation_name_after].name()
                )

                if orientation_name_before != orientation_name_after:
                    up_down_mappings += 1

    assert up_down_mappings != 0
