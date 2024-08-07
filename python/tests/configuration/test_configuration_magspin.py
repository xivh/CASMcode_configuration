import numpy as np

import libcasm.configuration as casmconfig
import libcasm.enumerate as casmenum
import libcasm.xtal as xtal


def test_Cmagspin_enum(FCC_binary_discrete_Cmagspin_prim):
    prim = casmconfig.Prim(FCC_binary_discrete_Cmagspin_prim)
    assert len(prim.factor_group.elements) == 96

    supercell_set = casmconfig.SupercellSet(prim)
    config_enum = casmenum.ConfigEnumAllOccupations(
        prim=prim, supercell_set=supercell_set
    )
    config_list = []
    for config in config_enum.by_supercell(max=1):
        config_list.append(config.copy())

    for config in config_list:
        print(xtal.pretty_json(config.to_dict()))
    assert len(config_list) == 2

    assert (config_list[0].occupation == [1]).all()
    assert (config_list[1].occupation == [3]).all()


def test_Cmagspin_to_structure(FCC_binary_discrete_Cmagspin_prim):
    prim = casmconfig.Prim(FCC_binary_discrete_Cmagspin_prim)

    supercell = casmconfig.Supercell(
        prim=prim,
        transformation_matrix_to_super=np.array(
            [
                [-1, 1, 1],
                [1, -1, 1],
                [1, 1, -1],
            ]
        ),
    )
    assert supercell.n_unitcells == 4

    configuration = casmconfig.Configuration(
        supercell=supercell,
    )
    assert (configuration.occupation == np.array([0, 0, 0, 0])).all()

    structure = configuration.to_structure()
    assert isinstance(structure, xtal.Structure)

    ### Default configuration ###
    data = structure.to_dict()
    assert isinstance(data, dict)

    assert np.allclose(
        data["atom_properties"]["Cmagspin"]["value"],
        np.array([[1.0, 1.0, 1.0, 1.0]]).T,
    )
    assert data["atom_type"] == ["A", "A", "A", "A"]

    ### Configuration [A.down, A.up, A.up, A.up] ###
    configuration.set_occ(0, 1)
    structure = configuration.to_structure()
    data = structure.to_dict()
    assert np.allclose(
        data["atom_properties"]["Cmagspin"]["value"],
        np.array([[-1.0, 1.0, 1.0, 1.0]]).T,
    )
    assert data["atom_type"] == ["A", "A", "A", "A"]

    ### Configuration [A.down, B.up, A.up, A.up] ###
    configuration.set_occ(1, 2)
    structure = configuration.to_structure()
    data = structure.to_dict()
    assert np.allclose(
        data["atom_properties"]["Cmagspin"]["value"],
        np.array([[-1.0, 1.0, 1.0, 1.0]]).T,
    )
    assert data["atom_type"] == ["A", "B", "A", "A"]

    ### Configuration [A.down, B.down, A.up, A.up] ###
    configuration.set_occ(1, 3)
    structure = configuration.to_structure()
    data = structure.to_dict()
    assert np.allclose(
        data["atom_properties"]["Cmagspin"]["value"],
        np.array([[-1.0, -1.0, 1.0, 1.0]]).T,
    )
    assert data["atom_type"] == ["A", "B", "A", "A"]
