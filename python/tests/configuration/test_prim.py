import json

import numpy as np
import pytest

import libcasm.configuration as config
import libcasm.configuration.io as config_io
import libcasm.configuration.io.spglib as spglib_io
import libcasm.sym_info as sym_info
import libcasm.xtal as xtal

from .functions import (
    check_magnetic_symmetry_dataset,
    check_spacegroup_type,
    check_symmetry_dataset,
)


def test_simple_cubic_binary_factor_group(simple_cubic_binary_prim):
    xtal_prim = simple_cubic_binary_prim
    factor_group = xtal.make_factor_group(xtal_prim)
    assert len(factor_group) == 48
    prim = config.Prim(xtal_prim)
    assert prim.xtal_prim.coordinate_frac().shape == (3, 1)
    assert len(prim.factor_group.elements) == 48

    assert prim.is_atomic
    assert prim.continuous_magspin_key is None
    assert prim.continuous_magspin_flavor is None
    assert prim.discrete_atomic_magspin_key is None
    assert prim.discrete_atomic_magspin_flavor is None

    lattice = prim.xtal_prim.lattice()
    for op in prim.factor_group.elements:
        syminfo = xtal.SymInfo(op, lattice)
        print(xtal.pretty_json(syminfo.to_dict()))


@pytest.mark.xfail(reason="known tolerance issue")
def test_symop_prec():
    # Relates to libcasm-xtal commit 28f1140

    # This tests whether symmetry operation matrices
    # have exact zeros instead of very small numbers.
    # It is not a failure, just normal floating point error
    # - but it would be easier to read output if the elements
    # that should be exactly zero are exactly zero,
    # and may be addressed in a future change.

    import libcasm.xtal.lattices as xtal_lattices
    import libcasm.xtal.prims as xtal_prims

    tetragonal_prim = xtal.Prim(
        lattice=xtal_lattices.tetragonal(a=1.0, c=1.5),
        coordinate_frac=np.zeros((3, 1)),
        occ_dof=[["A", "B"]],
        title="tetragonal",
    )

    rhombohedral_prim = xtal.Prim(
        lattice=xtal_lattices.rhombohedral(a=1.0, alpha=70.0),
        coordinate_frac=np.zeros((3, 1)),
        occ_dof=[["A", "B"]],
        title="tetragonal",
    )

    monoclinic_prim = xtal.Prim(
        lattice=xtal_lattices.monoclinic(
            a=1.0,
            b=1.2,
            c=1.5,
            beta=82.0,
        ),
        coordinate_frac=np.zeros((3, 1)),
        occ_dof=[["A", "B"]],
        title="monoclinic",
    )

    triclinic_prim = xtal.Prim(
        lattice=xtal_lattices.triclinic(
            a=1.0, b=1.2, c=1.5, alpha=73.0, beta=82.0, gamma=105.0
        ),
        coordinate_frac=np.zeros((3, 1)),
        occ_dof=[["A", "B"]],
        title="triclinic",
    )

    test_prims = [
        ("cubic", xtal_prims.cubic(a=1.0)),
        ("FCC-a", xtal_prims.FCC(a=1.0)),
        ("FCC-b", xtal_prims.FCC(r=0.5)),  # -> fails with libcasm-xtal 2.0a11
        ("FCC-c", xtal_prims.FCC(r=1.0)),  # -> fails
        ("BCC-a", xtal_prims.BCC(a=1.0)),
        ("BCC-b", xtal_prims.BCC(r=0.5)),  # -> fails
        ("BCC-c", xtal_prims.BCC(r=1.0)),  # -> fails
        ("HCP-a", xtal_prims.HCP(r=0.6)),
        ("HCP-a", xtal_prims.HCP(r=1.0)),
        ("HCP-b", xtal_prims.HCP(a=1.0, c=1.64)),
        ("tetragonal", tetragonal_prim),
        ("rhombohedral", rhombohedral_prim),  # -> fails
        ("monoclinic", monoclinic_prim),  # -> fails
        ("triclinic", triclinic_prim),  # -> fails
    ]

    for prim_type, x in test_prims:
        # print(f"prim_type: {prim_type}")
        # print("lattice:\n", x.lattice().column_vector_matrix())
        # print(np.linalg.inv(x.lattice().column_vector_matrix()))
        # for i, op in enumerate(prim.factor_group.elements):
        #     print(f"operation {i}")
        #     print(op.matrix())
        #     print(op.translation())
        #     print()

        prim = config.Prim(x)
        for i, op in enumerate(prim.factor_group.elements):
            R = op.matrix()
            tau = op.translation()
            assert (
                (R == 0.0) | (np.abs(R) > 1e-10)
            ).all(), f"prim_type: {prim_type}, operation {i}, matrix:\n{R}"
            assert (
                (tau == 0.0) | (np.abs(tau) > 1e-10)
            ).all(), f"prim_type: {prim_type}, operation {i}, translation:\n{tau}"


def test_simple_cubic_binary_conjugacy_classes(simple_cubic_binary_prim):
    xtal_prim = simple_cubic_binary_prim
    prim = config.Prim(xtal_prim)
    conjugacy_classes = prim.factor_group.conjugacy_classes()
    assert len(conjugacy_classes) == 10


def test_simple_cubic_binary_crystal_point_group(simple_cubic_binary_prim):
    xtal_prim = simple_cubic_binary_prim
    point_group = xtal.make_prim_crystal_point_group(xtal_prim)
    assert len(point_group) == 48
    prim = config.Prim(xtal_prim)
    assert prim.xtal_prim.coordinate_frac().shape == (3, 1)
    assert len(prim.crystal_point_group.elements) == 48


def test_simple_cubic_binary_lattice_point_group(simple_cubic_binary_prim):
    xtal_prim = simple_cubic_binary_prim
    prim = config.Prim(xtal_prim)
    assert len(prim.lattice_point_group.elements) == 48


def test_simple_cubic_binary_to_json_deprecated(simple_cubic_binary_prim):
    xtal_prim = simple_cubic_binary_prim
    prim = config.Prim(xtal_prim)
    prim_json = json.loads(prim.to_json())
    assert "basis" in prim_json
    assert "coordinate_mode" in prim_json
    assert "lattice_vectors" in prim_json


def test_from_json_deprecated():
    prim_json_str = """{
        "basis": [{
            "coordinate": [0.0, 0.0, 0.0],
            "occupants": ["A", "B"]}
        ],
        "coordinate_mode": "Fractional",
        "lattice_vectors": [
            [1.0, 0.0, 0.0],
            [0.0, 1.0, 0.0],
            [0.0, 0.0, 1.0]
        ],
        "title": "prim"}"""
    prim = config.Prim.from_json(prim_json_str)
    assert prim.xtal_prim.coordinate_frac().shape == (3, 1)
    assert len(prim.factor_group.elements) == 48


def test_simple_cubic_binary_to_dict(simple_cubic_binary_prim):
    xtal_prim = simple_cubic_binary_prim
    prim = config.Prim(xtal_prim)
    prim_data = prim.to_dict()
    assert "basis" in prim_data
    assert "coordinate_mode" in prim_data
    assert "lattice_vectors" in prim_data


def test_from_dict():
    prim_data = {
        "basis": [{"coordinate": [0.0, 0.0, 0.0], "occupants": ["A", "B"]}],
        "coordinate_mode": "Fractional",
        "lattice_vectors": [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]],
        "title": "prim",
    }
    prim = config.Prim.from_dict(prim_data)
    assert prim.xtal_prim.coordinate_frac().shape == (3, 1)
    assert len(prim.factor_group.elements) == 48


def test_simple_cubic_binary_repr(simple_cubic_binary_prim):
    xtal_prim = simple_cubic_binary_prim
    prim = config.Prim(xtal_prim)

    # Test print
    import io
    from contextlib import redirect_stdout

    f = io.StringIO()
    with redirect_stdout(f):
        print(prim)
    out = f.getvalue()

    assert "basis" in out
    assert "coordinate_mode" in out
    assert "lattice_vectors" in out


def test_symgroup_brief(simple_cubic_binary_prim):
    xtal_prim = simple_cubic_binary_prim
    prim = config.Prim(xtal_prim)
    lattice = xtal_prim.lattice()

    brief_cart_str = prim.factor_group.brief_cart(lattice, index_from=1)
    assert isinstance(brief_cart_str, str)

    brief_frac_str = prim.factor_group.brief_frac(lattice, index_from=1)
    assert isinstance(brief_frac_str, str)


def test_symgroup_to_from_dict(simple_cubic_binary_prim):
    xtal_prim = simple_cubic_binary_prim
    prim = config.Prim(xtal_prim)
    lattice = xtal_prim.lattice()

    data = prim.factor_group.to_dict(lattice)
    assert isinstance(data, dict)

    # print(xtal.pretty_json(data))

    factor_group_in = sym_info.SymGroup.from_dict(data, lattice)
    assert isinstance(factor_group_in, sym_info.SymGroup)
    assert len(factor_group_in.elements) == 48

    data_again = factor_group_in.to_dict(lattice)
    assert data == data_again


def test_prim_spglib_io(ZrO_prim):
    xtal_prim = ZrO_prim
    prim = config.Prim(xtal_prim)
    cell = spglib_io.as_cell(prim)

    assert isinstance(cell, tuple)
    assert len(cell) == 3

    symmetry_dataset = spglib_io.get_symmetry_dataset(prim)
    check_symmetry_dataset(
        symmetry_dataset,
        number=194,
        n_rotations=24,
        hall_number=488,
    )

    mag_symmetry_dataset = spglib_io.get_magnetic_symmetry_dataset(prim)
    check_magnetic_symmetry_dataset(
        mag_symmetry_dataset,
        uni_number=1494,
        n_operations=48,
        hall_number=488,
    )

    spacegroup_type = spglib_io.get_spacegroup_type_from_symmetry(
        elements=prim.factor_group.elements,
        lattice=xtal_prim.lattice(),
    )
    # from libcasm.configuration.io.spglib import asdict as spg_asdict
    # print(spg_asdict(spacegroup_type))
    check_spacegroup_type(
        spacegroup_type,
        number=194,
        hall_number=488,
    )


def test_symgroup_to_dict_with_group_classification_1(simple_cubic_binary_prim):
    xtal_prim = simple_cubic_binary_prim
    prim = config.Prim(xtal_prim)

    data = config_io.symgroup_to_dict_with_group_classification(prim, prim.factor_group)

    assert isinstance(data, dict)
    assert "spacegroup_type" in data["group_classification"]
    assert isinstance(data["group_classification"]["spacegroup_type"], dict)
    assert data["group_classification"]["spacegroup_type"]["number"] == 221
    assert "spacegroup_type_from_casm_symmetry" not in data["group_classification"]
    assert "magnetic_spacegroup_type" not in data["group_classification"]
    assert (
        "magnetic_spacegroup_type_from_casm_symmetry"
        not in data["group_classification"]
    )


def test_symgroup_to_dict_with_group_classification_2():
    lattice = xtal.Lattice(
        column_vector_matrix=np.array(
            [
                [0.0, 0.5, 0.5],
                [0.5, 0.0, 0.5],
                [0.5, 0.5, 0.0],
            ]
        ).transpose(),
    )
    coordinate_cart = np.array(
        [
            [0.0, 0.0, 0.0],
            [0.25, 0.25, 0.25],
        ]
    ).transpose()
    xtal_prim = xtal.Prim(
        lattice=lattice,
        coordinate_frac=xtal.fractional_to_cartesian(lattice, coordinate_cart),
        occ_dof=[["A"], ["A"]],
    )
    prim = config.Prim(xtal_prim)

    # factor group
    data = config_io.symgroup_to_dict_with_group_classification(
        xtal_prim, prim.factor_group
    )

    assert isinstance(data, dict)
    # print(xtal.pretty_json(data["group_classification"]))
    assert "spacegroup_type" in data["group_classification"]
    assert data["group_classification"]["spacegroup_type"]["number"] == 227
    assert "spacegroup_type_from_casm_symmetry" not in data["group_classification"]
    assert "magnetic_spacegroup_type" not in data["group_classification"]
    assert (
        "magnetic_spacegroup_type_from_casm_symmetry"
        not in data["group_classification"]
    )

    # lattice point group
    data = config_io.symgroup_to_dict_with_group_classification(
        xtal_prim.lattice(), prim.lattice_point_group
    )

    assert isinstance(data, dict)
    # print(xtal.pretty_json(data["group_classification"]))
    assert "spacegroup_type" in data["group_classification"]
    assert data["group_classification"]["spacegroup_type"]["number"] == 225
    assert "spacegroup_type_from_casm_symmetry" not in data["group_classification"]
    assert "magnetic_spacegroup_type" not in data["group_classification"]
    assert (
        "magnetic_spacegroup_type_from_casm_symmetry"
        not in data["group_classification"]
    )


def test_prim_occ_symgroup_rep():
    occ_dof = [
        ["A", "B"],
        ["B", "C", "D"],
        ["B", "D", "C"],
        ["C", "D", "B"],
    ]
    xtal_prim = xtal.Prim(
        lattice=xtal.Lattice(
            column_vector_matrix=np.eye(3),
        ),
        coordinate_frac=np.array(
            [
                [0.0, 0.0, 0.0],
                [0.0, 0.5, 0.5],
                [0.5, 0.0, 0.5],
                [0.5, 0.5, 0.0],
            ]
        ).T,
        occ_dof=occ_dof,
    )
    prim = config.Prim(xtal_prim)
    occ_symgroup_rep = prim.occ_symgroup_rep

    for i_factor_group, occ_op_rep in enumerate(occ_symgroup_rep):
        site_rep = prim.integral_site_coordinate_symgroup_rep[i_factor_group]
        for i_sublat_before, occ_sublat_rep in enumerate(occ_op_rep):
            site_before = xtal.IntegralSiteCoordinate(i_sublat_before, [0, 0, 0])
            site_after = site_rep * site_before
            i_sublat_after = site_after.sublattice()
            for i_occ_before in range(len(occ_sublat_rep)):
                i_occ_after = occ_sublat_rep[i_occ_before]
                assert (
                    occ_dof[i_sublat_before][i_occ_before]
                    == occ_dof[i_sublat_after][i_occ_after]
                )


def test_prim_atom_position_symgroup_rep():
    mol_x = xtal.Occupant(
        name="mol",
        atoms=[
            xtal.AtomComponent(name="B", coordinate=[-0.1, 0.0, 0.0], properties={}),
            xtal.AtomComponent(name="B", coordinate=[0.1, 0.0, 0.0], properties={}),
        ],
    )
    mol_y = xtal.Occupant(
        name="mol",
        atoms=[
            xtal.AtomComponent(name="B", coordinate=[0.0, -0.1, 0.0], properties={}),
            xtal.AtomComponent(name="B", coordinate=[0.0, 0.1, 0.0], properties={}),
        ],
    )
    mol_z = xtal.Occupant(
        name="mol",
        atoms=[
            xtal.AtomComponent(name="B", coordinate=[0.0, 0.0, -0.1], properties={}),
            xtal.AtomComponent(name="B", coordinate=[0.0, 0.0, 0.1], properties={}),
        ],
    )
    atom_A = xtal.Occupant(
        name="A",
        atoms=[
            xtal.AtomComponent(name="A", coordinate=[0.0, 0.0, 0.0], properties={}),
        ],
    )
    occupants = {"mol.x": mol_x, "mol.y": mol_y, "mol.z": mol_z, "A": atom_A}
    occ_dof = [
        ["A"],
        ["mol.x", "mol.y", "mol.z", "A"],
        ["mol.x", "mol.z", "mol.y", "A"],
        ["mol.y", "mol.z", "A", "mol.x"],
    ]
    xtal_prim = xtal.Prim(
        lattice=xtal.Lattice(
            column_vector_matrix=np.eye(3),
        ),
        coordinate_frac=np.array(
            [
                [0.0, 0.0, 0.0],
                [0.0, 0.5, 0.5],
                [0.5, 0.0, 0.5],
                [0.5, 0.5, 0.0],
            ]
        ).T,
        occ_dof=occ_dof,
        occupants=occupants,
    )
    prim = config.Prim(xtal_prim)
    assert len(prim.factor_group.elements) == 48

    coordinate_cart = xtal_prim.coordinate_cart()
    occ_symgroup_rep = prim.occ_symgroup_rep
    atom_position_symgroup_rep = prim.atom_position_symgroup_rep

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

                atom_position_rep = atom_position_symgroup_rep[i_factor_group][
                    i_sublat_before
                ][i_occ_before]
                for i_atom_before in range(len(atom_position_rep)):
                    i_atom_after = atom_position_rep[i_atom_before]

                    occ_before = occupants[occ_dof[i_sublat_before][i_occ_before]]
                    atom_before = occ_before.atoms()[i_atom_before]
                    occ_after = occupants[occ_dof[i_sublat_after][i_occ_after]]
                    atom_after = occ_after.atoms()[i_atom_after]

                    # assert atom names map
                    assert atom_before.name() == atom_after.name()

                    # assert atom positions map
                    cart_before = (
                        coordinate_cart[:, i_sublat_before] + atom_before.coordinate()
                    )
                    sym_op = prim.factor_group.elements[i_factor_group]
                    cart_after = (
                        coordinate_cart[:, i_sublat_after] + atom_after.coordinate()
                    )

                    d = xtal.min_periodic_displacement(
                        lattice=xtal_prim.lattice(),
                        r1=sym_op.matrix() @ cart_before + sym_op.translation(),
                        r2=cart_after,
                    )
                    assert np.allclose(d, np.zeros((3,)))
