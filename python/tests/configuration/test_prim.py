import json

import libcasm.configuration as config
import libcasm.xtal as xtal


def test_simple_cubic_binary_factor_group(simple_cubic_binary_prim):
    xtal_prim = simple_cubic_binary_prim
    factor_group = xtal.make_factor_group(xtal_prim)
    assert len(factor_group) == 48
    prim = config.Prim(xtal_prim)
    assert prim.xtal_prim().coordinate_frac().shape == (3, 1)
    assert len(prim.factor_group().elements()) == 48

    lattice = prim.xtal_prim().lattice()
    for op in prim.factor_group().elements():
        syminfo = xtal.SymInfo(op, lattice)
        print(syminfo.to_json())


def test_simple_cubic_binary_conjugacy_classes(simple_cubic_binary_prim):
    xtal_prim = simple_cubic_binary_prim
    prim = config.Prim(xtal_prim)
    conjugacy_classes = prim.factor_group().conjugacy_classes()
    assert len(conjugacy_classes) == 10


def test_simple_cubic_binary_crystal_point_group(simple_cubic_binary_prim):
    xtal_prim = simple_cubic_binary_prim
    point_group = xtal.make_prim_crystal_point_group(xtal_prim)
    assert len(point_group) == 48
    prim = config.Prim(xtal_prim)
    assert prim.xtal_prim().coordinate_frac().shape == (3, 1)
    assert len(prim.crystal_point_group().elements()) == 48


def test_simple_cubic_binary_to_json(simple_cubic_binary_prim):
    xtal_prim = simple_cubic_binary_prim
    prim = config.Prim(xtal_prim)
    prim_json = json.loads(prim.to_json())
    assert "basis" in prim_json
    assert "coordinate_mode" in prim_json
    assert "lattice_vectors" in prim_json


def test_from_json():
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
    assert prim.xtal_prim().coordinate_frac().shape == (3, 1)
    assert len(prim.factor_group().elements()) == 48
