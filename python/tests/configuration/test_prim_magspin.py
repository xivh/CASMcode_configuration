from typing import Any

import numpy as np
import pytest

import libcasm.configuration as config
import libcasm.configuration.io as config_io
import libcasm.configuration.io.spglib as spglib_io
import libcasm.xtal as xtal


def make_discrete_magnetic_atom(
    name: str,
    value: Any,
    flavor: str = "C",
) -> xtal.Occupant:
    """Construct a discrete magnetic atomic occupant

    Parameters
    ----------
    name: str
        A "chemical name", which must be identical for atoms to be found symmetrically
        equivalent. The names are case sensitive, and “Va” is reserved for vacancies.
    value: Any
        The discrete value of the magnetic spin to associate with the constructed
        Occupant. If the type is `int` or `float` the value is converted to a
        size 1 array of float, other types are converted using `numpy.asarray`.
    flavor: str
        The magnetic spin "flavor", which must be one of varieties supported by CASM:
        `C`, `NC`, `SO`.

    Returns
    -------
    discrete_magnetic_atom: xtal.Occupant
        An :class:`~libcasm.xtal.Occupant` consisting of a single atom with the
        specified magnetic spin flavor and value.
    """
    if isinstance(value, (int, float)):
        value = np.array([value], dtype=np.float64)
    else:
        value = np.asarray(value, dtype=np.float64)

    return xtal.Occupant(
        name=name,
        atoms=[
            xtal.AtomComponent(
                name=name,
                coordinate=[0.0, 0.0, 0.0],
                properties={flavor + "magspin": value},
            )
        ],
    )


def FCC_binary_discrete_Cmagspin_prim():
    # Lattice vectors
    lattice_column_vector_matrix = np.array(
        [
            [0.0, 1.0 / 2.0, 1.0 / 2.0],  # a
            [1.0 / 2.0, 0.0, 1.0 / 2.0],  # a
            [1.0 / 2.0, 1.0 / 2.0, 0.0],  # a
        ]
    ).transpose()
    lattice = xtal.Lattice(lattice_column_vector_matrix)

    # Basis sites positions, as columns of a matrix,
    # in fractional coordinates with respect to the lattice vectors
    coordinate_frac = np.array(
        [
            [0.0, 0.0, 0.0],
        ]
    ).transpose()

    # Occupation degrees of freedom (DoF)
    occupants = {}
    occupants = {
        "A.up": make_discrete_magnetic_atom(name="A", value=1, flavor="C"),
        "A.down": make_discrete_magnetic_atom(name="A", value=-1, flavor="C"),
        "B.up": make_discrete_magnetic_atom(name="B", value=1, flavor="C"),
        "B.down": make_discrete_magnetic_atom(name="B", value=-1, flavor="C"),
    }
    occ_dof = [["A.up", "A.down", "B.up", "B.down"]]

    return xtal.Prim(
        lattice=lattice,
        coordinate_frac=coordinate_frac,
        occ_dof=occ_dof,
        occupants=occupants,
    )


def test_FCC_binary_discrete_Cmagspin_prim_1():
    xtal_prim = FCC_binary_discrete_Cmagspin_prim()
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


def test_prim_magspin_spglib_io():
    xtal_prim = FCC_binary_discrete_Cmagspin_prim()
    prim = config.Prim(xtal_prim)
    cell = spglib_io.as_cell(prim)

    assert isinstance(cell, tuple)
    assert len(cell) == 4

    symmetry_dataset = spglib_io.get_symmetry_dataset(prim)

    # print(symmetry_dataset.keys())
    # for key, value in symmetry_dataset.items():
    #     print("symmetry", key)
    #     print(value)
    #     print()

    assert isinstance(symmetry_dataset, dict)
    assert symmetry_dataset["number"] == 225
    assert symmetry_dataset["hall_number"] == 523
    assert len(symmetry_dataset["rotations"]) == 48

    mag_symmetry_dataset = spglib_io.get_magnetic_symmetry_dataset(prim)

    # print(mag_symmetry_dataset.keys())
    # for key, value in mag_symmetry_dataset.items():
    #     print("mag_symmetry", key)
    #     print(value)
    #     print()

    assert isinstance(mag_symmetry_dataset, dict)
    assert mag_symmetry_dataset["uni_number"] == 1619
    assert mag_symmetry_dataset["hall_number"] == 523
    assert mag_symmetry_dataset["n_operations"] == 96

    # print(prim.factor_group.brief_cart(xtal_prim.lattice()))

    with pytest.raises(
        ValueError, match="not valid for operations with time reversal symmetry"
    ):
        _ = spglib_io.get_spacegroup_type_from_symmetry(
            elements=prim.factor_group.elements,
            lattice=xtal_prim.lattice(),
        )


def test_symgroup_to_dict_with_group_classification():
    xtal_prim = FCC_binary_discrete_Cmagspin_prim()
    prim = config.Prim(xtal_prim)

    data = config_io.symgroup_to_dict_with_group_classification(prim, prim.factor_group)

    assert isinstance(data, dict)
    assert "spacegroup_type" in data["group_classification"]
    assert data["group_classification"]["spacegroup_type"]["number"] == 225
    assert "magnetic_spacegroup_type" in data["group_classification"]
    assert (
        data["group_classification"]["magnetic_spacegroup_type"]["uni_number"] == 1619
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
