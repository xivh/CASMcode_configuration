import copy
import math

import numpy as np
from sortedcontainers import SortedList

import libcasm.clexulator as casmclex
import libcasm.configuration as casmconfig
import libcasm.configuration.io.spglib as spglib_io

from .functions import check_symmetry_dataset


def test_configuration_constructor(simple_cubic_binary_prim):
    prim = casmconfig.Prim(simple_cubic_binary_prim)
    T = np.array(
        [
            [1, 0, 0],
            [0, 1, 0],
            [0, 0, 1],
        ]
    )
    supercell = casmconfig.make_canonical_supercell(casmconfig.Supercell(prim, T))
    configuration = casmconfig.Configuration(supercell)
    assert isinstance(configuration, casmconfig.Configuration)


def test_configuration_dof_values(simple_cubic_binary_prim):
    prim = casmconfig.Prim(simple_cubic_binary_prim)
    T = np.array(
        [
            [1, 0, 0],
            [0, 1, 0],
            [0, 0, 1],
        ]
    )
    supercell = casmconfig.make_canonical_supercell(casmconfig.Supercell(prim, T))
    configuration = casmconfig.Configuration(supercell)
    dof_values = configuration.dof_values
    assert isinstance(dof_values, casmclex.ConfigDoFValues)
    assert dof_values is configuration.dof_values

    assert configuration.occ(0) == 0
    configuration.dof_values.set_occupation([1])
    assert configuration.occ(0) == 1


def test_configuration_occupation(simple_cubic_binary_prim):
    prim = casmconfig.Prim(simple_cubic_binary_prim)
    T = np.array(
        [
            [4, 0, 0],
            [0, 4, 0],
            [0, 0, 4],
        ]
    )
    supercell = casmconfig.make_canonical_supercell(casmconfig.Supercell(prim, T))
    configuration = casmconfig.Configuration(supercell)

    occupation = configuration.occupation
    assert occupation.shape == (64,)
    assert configuration.occ(0) == 0
    assert (occupation == np.array([0] * 64)).all()

    configuration.set_occ(0, 1)
    occupation = configuration.occupation
    assert configuration.occ(0) == 1
    assert (occupation == np.array([1] + [0] * 63)).all()

    configuration.set_occupation([1, 1] + [0] * 62)
    assert (occupation == np.array([1, 1] + [0] * 62)).all()

    occupation = np.array([1] * 3 + [0] * 61, dtype=int)
    configuration.set_occupation(occupation)
    assert (configuration.occupation == np.array([1] * 3 + [0] * 61)).all()
    assert configuration.occupation.shape == (64,)

    supercell_occ_dof = configuration.supercell.occ_dof()
    occupation_dof = [
        supercell_occ_dof[site_index][s]
        for site_index, s in enumerate(configuration.occupation)
    ]
    assert occupation_dof == ["B"] * 3 + ["A"] * 61


def test_canonical_configuration_occupation(simple_cubic_binary_prim):
    prim = casmconfig.Prim(simple_cubic_binary_prim)
    T = np.array(
        [
            [4, 0, 0],
            [0, 4, 0],
            [0, 0, 4],
        ]
    )
    supercell = casmconfig.make_canonical_supercell(casmconfig.Supercell(prim, T))
    configuration = casmconfig.Configuration(supercell)

    configuration.set_occ(0, 1)
    assert casmconfig.is_canonical_configuration(configuration) is True

    configuration.set_occ(0, 0)
    configuration.set_occ(10, 1)
    canon_config = casmconfig.make_canonical_configuration(configuration)
    assert casmconfig.is_canonical_configuration(canon_config) is True
    assert (canon_config.occupation == np.array([1] + [0] * 63)).all()


def test_configuration_invariant_subgroup(simple_cubic_binary_prim):
    prim = casmconfig.Prim(simple_cubic_binary_prim)
    T = np.array(
        [
            [4, 0, 0],
            [0, 4, 0],
            [0, 0, 4],
        ]
    )
    supercell = casmconfig.make_canonical_supercell(casmconfig.Supercell(prim, T))
    configuration = casmconfig.Configuration(supercell)

    invariant_subgroup = casmconfig.make_invariant_subgroup(configuration)
    assert len(invariant_subgroup) == 64 * 48

    configuration.set_occ(0, 1)
    invariant_subgroup = casmconfig.make_invariant_subgroup(configuration)
    assert len(invariant_subgroup) == 48


def test_configuration_spglib_io(simple_cubic_binary_prim):
    prim = casmconfig.Prim(simple_cubic_binary_prim)
    T = np.array(
        [
            [4, 0, 0],
            [0, 4, 0],
            [0, 0, 4],
        ]
    )
    supercell = casmconfig.make_canonical_supercell(casmconfig.Supercell(prim, T))
    configuration = casmconfig.Configuration(supercell)

    cell = spglib_io.as_cell(configuration)

    assert isinstance(cell, tuple)
    assert len(cell) == 3

    symmetry_dataset = spglib_io.get_symmetry_dataset(configuration)
    check_symmetry_dataset(symmetry_dataset, number=221, n_rotations=64 * 48)

    configuration.set_occ(0, 1)
    symmetry_dataset = spglib_io.get_symmetry_dataset(configuration)
    check_symmetry_dataset(symmetry_dataset, number=221, n_rotations=48)


def test_configuration_apply_1(simple_cubic_binary_prim):
    prim = casmconfig.Prim(simple_cubic_binary_prim)
    T = np.array(
        [
            [4, 0, 0],
            [0, 4, 0],
            [0, 0, 4],
        ]
    )
    supercell = casmconfig.make_canonical_supercell(casmconfig.Supercell(prim, T))
    config = casmconfig.Configuration(supercell)
    config.set_occ(l=0, s=1)
    config.set_occ(l=1, s=1)

    for rep in config.supercell.symgroup_rep():
        assert isinstance(rep, casmconfig.SupercellSymOp)

        # option 1: use the multiplication operator
        transformed_config_1 = rep * config
        assert transformed_config_1 is not config

        # option 2: use the copy_apply method
        transformed_config_2 = casmconfig.copy_apply(rep, config)
        assert transformed_config_2 is not config

        # option 3: copy, then apply rep in-place
        copied_config = config.copy()
        assert copied_config is not config
        transformed_config_3 = casmconfig.apply(rep, copied_config)
        assert transformed_config_3 is copied_config


def test_configuration_apply_2(simple_cubic_binary_prim):
    prim = casmconfig.Prim(simple_cubic_binary_prim)
    T = np.array(
        [
            [4, 0, 0],
            [0, 4, 0],
            [0, 0, 4],
        ]
    )
    supercell = casmconfig.make_canonical_supercell(casmconfig.Supercell(prim, T))
    configuration = casmconfig.Configuration(supercell)
    configuration.set_occ(l=0, s=1)
    configuration.set_occ(l=1, s=1)

    # Test SupercellSymOp using begin & end
    rep = casmconfig.SupercellSymOp.begin(supercell)
    end = casmconfig.SupercellSymOp.end(supercell)
    i = 0
    equivs = SortedList()
    while rep != end:
        transformed = rep * configuration
        if transformed not in equivs:
            equivs.add(transformed)
        rep.next()
        i += 1
    assert i == 48 * 64
    assert len(equivs) == 192

    # Test iterating using Supercell.symgroup_rep directly
    i = 0
    equivs = SortedList()
    for rep in supercell.symgroup_rep():
        transformed = rep * configuration
        if transformed not in equivs:
            equivs.add(transformed)
        i += 1
    assert i == 48 * 64
    assert len(equivs) == 192

    # Test SupercellSymOp storing Supercell.symgroup_rep()
    i = 0
    equivs = SortedList()
    configuration_rep = supercell.symgroup_rep()
    assert len(configuration_rep) == 48 * 64
    for rep in configuration_rep:
        transformed = rep * configuration
        if transformed not in equivs:
            equivs.add(transformed)
        i += 1
    assert i == 48 * 64
    assert len(equivs) == 192

    # Test SupercellSymOp copy.deepcopy
    configuration_rep = []
    rep = casmconfig.SupercellSymOp.begin(supercell)
    end = casmconfig.SupercellSymOp.end(supercell)
    while rep != end:
        configuration_rep.append(copy.deepcopy(rep))
        rep.next()
    assert len(configuration_rep) == 48 * 64
    assert configuration_rep[0] is not configuration_rep[1]

    i = 0
    equivs = SortedList()
    for rep in configuration_rep:
        transformed = rep * configuration
        if transformed not in equivs:
            equivs.add(transformed)
        i += 1
    assert i == 48 * 64
    assert len(equivs) == 192

    # Test SupercellSymOp.copy
    configuration_rep = []
    rep = casmconfig.SupercellSymOp.begin(supercell)
    end = casmconfig.SupercellSymOp.end(supercell)
    while rep != end:
        configuration_rep.append(rep.copy())
        rep.next()
    assert len(configuration_rep) == 48 * 64
    assert configuration_rep[0] is not configuration_rep[1]

    i = 0
    equivs = SortedList()
    for rep in configuration_rep:
        transformed = rep * configuration
        if transformed not in equivs:
            equivs.add(transformed)
        i += 1
    assert i == 48 * 64
    assert len(equivs) == 192

    # Test SupercellSymOp.__rep__
    import io
    from contextlib import redirect_stdout

    for rep in supercell.symgroup_rep():
        f = io.StringIO()
        with redirect_stdout(f):
            print(rep)
        out = f.getvalue()
        assert "combined_op" in out
        assert "supercell_factor_group_index" in out
        assert "translation_index" in out


def test_copy_configuration(simple_cubic_binary_prim):
    import copy

    prim = casmconfig.Prim(simple_cubic_binary_prim)
    T = np.array(
        [
            [4, 0, 0],
            [0, 4, 0],
            [0, 0, 4],
        ]
    )
    supercell = casmconfig.make_canonical_supercell(casmconfig.Supercell(prim, T))
    configuration1 = casmconfig.Configuration(supercell)
    configuration2 = copy.copy(configuration1)
    configuration3 = copy.deepcopy(configuration1)

    assert isinstance(configuration1, casmconfig.Configuration)
    assert isinstance(configuration2, casmconfig.Configuration)
    assert isinstance(configuration3, casmconfig.Configuration)
    assert configuration2 is not configuration1
    assert configuration3 is not configuration1


def test_configuration_to_from_dict(FCC_binary_Hstrain_noshear_disp_nodz_prim):
    import io
    from contextlib import redirect_stdout

    prim = casmconfig.Prim(FCC_binary_Hstrain_noshear_disp_nodz_prim)
    T = np.array(
        [
            [-1, 1, 1],
            [1, -1, 1],
            [1, 1, -1],
        ],
        dtype="int",
    )
    supercell = casmconfig.make_canonical_supercell(casmconfig.Supercell(prim, T))

    disp_init_prim = np.array(
        [
            [0.01, 0.05],
            [0.02, 0.06],
            [0.03, 0.07],
            [0.04, 0.08],
        ]
    ).transpose()
    disp_init_standard = np.array(
        [
            [0.01, 0.05, 0.0],
            [0.02, 0.06, 0.0],
            [0.03, 0.07, 0.0],
            [0.04, 0.08, 0.0],
        ]
    ).transpose()
    Hstrain_init_prim = np.array([0.1, 0.0, 0.0])
    Hstrain_init_standard = np.array(
        [
            0.1 / math.sqrt(3.0),
            0.1 / math.sqrt(3.0),
            0.1 / math.sqrt(3.0),
            0.0,
            0.0,
            0.0,
        ]
    )

    configuration = casmconfig.Configuration(supercell)
    configuration.set_local_dof_values(
        key="disp",
        dof_values=disp_init_prim,
    )
    configuration.set_global_dof_values(key="Hstrain", dof_values=Hstrain_init_prim)

    # default: write DoF values in standard basis
    data = configuration.to_dict()
    assert isinstance(data, dict)
    assert "basis" in data and data["basis"] == "standard"
    assert np.allclose(
        data["dof"]["global_dofs"]["Hstrain"]["values"],
        Hstrain_init_standard,
    )
    assert np.allclose(
        data["dof"]["local_dofs"]["disp"]["values"],
        disp_init_standard.transpose(),
    )

    # read back in
    supercell_set = casmconfig.SupercellSet(prim)
    configuration_in_1 = casmconfig.Configuration.from_dict(
        data=data, supercells=supercell_set
    )
    assert isinstance(configuration_in_1, casmconfig.Configuration)
    assert (configuration_in_1.supercell.transformation_matrix_to_super == T).all()
    assert len(supercell_set) == 1
    assert np.allclose(
        configuration_in_1.global_dof_values("Hstrain"),
        Hstrain_init_prim,
    )
    assert np.allclose(
        configuration_in_1.local_dof_values("disp"),
        disp_init_prim,
    )
    f = io.StringIO()
    with redirect_stdout(f):
        print(configuration_in_1)
    out = f.getvalue()
    assert "dof" in out

    # write DoF values in prim basis
    data = configuration.to_dict(write_prim_basis=True)
    assert isinstance(data, dict)
    assert "basis" in data and data["basis"] == "prim"
    assert np.allclose(
        data["dof"]["global_dofs"]["Hstrain"]["values"],
        Hstrain_init_prim,
    )
    assert np.allclose(
        data["dof"]["local_dofs"]["disp"]["values"],
        disp_init_prim.transpose(),
    )

    # read back in
    supercell_set = casmconfig.SupercellSet(prim)
    configuration_in_2 = casmconfig.Configuration.from_dict(
        data=data, supercells=supercell_set
    )
    assert isinstance(configuration_in_2, casmconfig.Configuration)
    assert (configuration_in_2.supercell.transformation_matrix_to_super == T).all()
    assert len(supercell_set) == 1
    assert np.allclose(
        configuration_in_1.global_dof_values("Hstrain"),
        Hstrain_init_prim,
    )
    assert np.allclose(
        configuration_in_1.local_dof_values("disp"),
        disp_init_prim,
    )
    f = io.StringIO()
    with redirect_stdout(f):
        print(configuration_in_2)
    out = f.getvalue()
    assert "dof" in out
