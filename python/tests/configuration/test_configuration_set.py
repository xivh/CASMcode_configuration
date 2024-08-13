import copy

import numpy as np
import pytest

import libcasm.configuration as config
import libcasm.xtal as xtal


def test_ConfigurationSet_constructor_1(simple_cubic_binary_prim):
    config.Prim(simple_cubic_binary_prim)
    configurations = config.ConfigurationSet()

    assert isinstance(configurations, config.ConfigurationSet)
    assert configurations.empty()


def test_ConfigurationSet_add_remove_discard_get_1(simple_cubic_binary_prim):
    prim = config.Prim(simple_cubic_binary_prim)
    configurations = config.ConfigurationSet()

    T = np.array(
        [
            [2, 1, 0],
            [0, 1, 0],
            [0, 0, 1],
        ]
    )
    supercell = config.make_canonical_supercell(config.Supercell(prim, T))
    configuration = config.Configuration(supercell)

    # add remove
    record = configurations.add(configuration)

    assert len(configurations) == 1
    assert isinstance(record, config.ConfigurationRecord)
    assert isinstance(record.configuration, config.Configuration)
    assert record.configuration_id == "0"
    assert record.supercell_name == "SCEL2_1_2_1_1_0_0"
    assert record.configuration_name == "SCEL2_1_2_1_1_0_0/0"
    assert record.configuration == configuration
    assert record.configuration is not configuration

    configurations.remove(configuration)
    assert len(configurations) == 0

    with pytest.raises(KeyError):
        configurations.remove(configuration)

    # add discard
    record = configurations.add(configuration)
    assert len(configurations) == 1

    configurations.discard(configuration)
    assert len(configurations) == 0

    configurations.discard(configuration)
    assert len(configurations) == 0

    # add get
    record_1 = configurations.add(configuration)
    assert len(configurations) == 1

    record_2 = configurations.get(configuration)
    assert record_1 == record_2

    record_3 = configurations.get(record_1.configuration_name)
    assert record_1 == record_3

    configuration_name = copy.copy(record_1.configuration_name)
    configurations.remove(configuration_name)

    record_4 = configurations.get(configuration)
    assert record_4 is None

    record_5 = configurations.get(configuration_name)
    assert record_5 is None


def test_ConfigurationSet_to_dict_1(simple_cubic_binary_prim):
    prim = config.Prim(simple_cubic_binary_prim)
    configurations = config.ConfigurationSet()

    T = np.array(
        [
            [2, 0, 0],
            [0, 2, 0],
            [0, 0, 1],
        ]
    )
    supercell = config.make_canonical_supercell(config.Supercell(prim, T))
    configuration = config.Configuration(supercell)

    # add configurations
    for i in range(supercell.n_unitcells):
        print(f"{i}:")
        x = copy.copy(configuration)
        x.set_occ(i, 1)
        print(xtal.pretty_json(x.to_dict(write_prim_basis=False)))
        print(xtal.pretty_json(x.to_dict(write_prim_basis=True)))
        configurations.add(x)
        print()
    assert len(configurations) == 4

    # output to JSON data
    data = configurations.to_dict()
    print(xtal.pretty_json(data))

    # clear
    configurations.clear()
    assert len(configurations) == 0

    # read back in
    supercells = config.SupercellSet(prim)
    configurations_in = config.ConfigurationSet.from_dict(data, supercells)
    assert len(configurations_in) == 4


def test_ConfigurationSet_to_dict_2(FCC_binary_Hstrain_noshear_prim):
    prim = config.Prim(FCC_binary_Hstrain_noshear_prim)
    configurations = config.ConfigurationSet()

    T = np.array(
        [
            [2, 0, 0],
            [0, 2, 0],
            [0, 0, 1],
        ]
    )
    supercell = config.make_canonical_supercell(config.Supercell(prim, T))
    configuration = config.Configuration(supercell)

    for i in range(supercell.n_unitcells):
        x = copy.copy(configuration)
        x.set_occ(i, 1)
        configurations.add(x)
    assert len(configurations) == 4

    configurations.clear()
    assert len(configurations) == 0

    for i in range(supercell.n_unitcells):
        x = copy.copy(configuration)
        x.set_occ(i, 1)
        configurations.add(config.make_canonical_configuration(x))
    assert len(configurations) == 1


def test_ConfigurationRecord_repr(simple_cubic_binary_prim):
    prim = config.Prim(simple_cubic_binary_prim)
    configurations = config.ConfigurationSet()

    T = np.array(
        [
            [2, 0, 0],
            [0, 2, 0],
            [0, 0, 1],
        ]
    )
    supercell = config.make_canonical_supercell(config.Supercell(prim, T))
    configuration = config.Configuration(supercell)

    # add configurations
    for i in range(supercell.n_unitcells):
        x = configuration.copy()
        x.set_occ(i, 1)
        configurations.add(x)
    assert len(configurations) == 4

    # Test print ConfigurationRecord
    import io
    from contextlib import redirect_stdout

    for record in configurations:
        assert isinstance(record, config.ConfigurationRecord)
        f = io.StringIO()
        with redirect_stdout(f):
            print(record)
        out = f.getvalue()
        assert "configuration" in out
        assert "supercell_name" in out
        assert "configuration_id" in out
        assert "configuration_name" in out
