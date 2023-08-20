import numpy as np

import libcasm.configuration as config
import libcasm.configuration.io as config_io


def test_supercell_list_io(simple_cubic_binary_prim):
    prim = config.Prim(simple_cubic_binary_prim)

    T = np.array(
        [
            [2, 1, 0],
            [0, 1, 0],
            [0, 0, 1],
        ]
    )
    supercell1 = config.Supercell(prim, T)
    supercell2 = config.make_canonical_supercell(supercell1)

    # add supercells to List[Supercell]
    supercell_list = [
        supercell1,
        supercell2,
    ]

    # ~~~
    data_list = config_io.supercell_list_to_data(supercell_list)
    for x in data_list:
        assert isinstance(x, dict)

    # using prim arg
    supercell_list_2 = config_io.supercell_list_from_data(data_list, prim=prim)
    assert isinstance(supercell_list_2, list)
    assert len(supercell_list_2) == 2

    # using supercells arg
    supercellset = config.SupercellSet(prim)
    supercell_list_3 = config_io.supercell_list_from_data(
        data_list, supercells=supercellset
    )
    assert isinstance(supercell_list_3, list)
    assert len(supercell_list_3) == 2
    assert len(supercellset) == 2


def test_configuration_list_io(simple_cubic_binary_prim):
    prim = config.Prim(simple_cubic_binary_prim)

    T = np.array(
        [
            [2, 1, 0],
            [0, 1, 0],
            [0, 0, 1],
        ]
    )
    supercell = config.Supercell(prim, T)
    canonical_supercell = config.make_canonical_supercell(supercell)

    configuration1 = config.Configuration(supercell)
    configuration2 = config.Configuration(canonical_supercell)

    # add configurations to List[Configuration]
    configuration_list = [
        configuration1,
        configuration2,
    ]

    # ~~~
    data_list = config_io.configuration_list_to_data(configuration_list)
    for x in data_list:
        assert isinstance(x, dict)

    # using prim arg
    configuration_list_2 = config_io.configuration_list_from_data(data_list, prim=prim)
    assert isinstance(configuration_list_2, list)
    assert len(configuration_list_2) == 2

    # using supercells arg
    supercellset = config.SupercellSet(prim)
    configuration_list_3 = config_io.configuration_list_from_data(
        data_list, supercells=supercellset
    )
    assert isinstance(configuration_list_3, list)
    assert len(configuration_list_3) == 2
    assert len(supercellset) == 2
