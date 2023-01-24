import math
import numpy as np
import pytest
from sortedcontainers import SortedList
import libcasm.xtal as xtal
import libcasm.clexulator as casmclex
import libcasm.configuration as config


def test_configuration_constructor(simple_cubic_binary_prim):
    prim = config.Prim(simple_cubic_binary_prim)
    T = np.array([
        [1, 0, 0],
        [0, 1, 0],
        [0, 0, 1],
    ])
    supercell = config.make_canonical_supercell(config.Supercell(prim, T))
    configuration = config.Configuration(supercell)
    assert type(configuration) == config.Configuration

def test_configuration_dof_values(simple_cubic_binary_prim):
    prim = config.Prim(simple_cubic_binary_prim)
    T = np.array([
        [1, 0, 0],
        [0, 1, 0],
        [0, 0, 1],
    ])
    supercell = config.make_canonical_supercell(config.Supercell(prim, T))
    configuration = config.Configuration(supercell)
    dof_values = configuration.dof_values()
    assert type(dof_values) == casmclex.ConfigDoFValues

    # Configuration.dof_values() returns a *copy*
    assert configuration.occ(0) == 0
    configuration.dof_values().set_occupation([1])
    assert configuration.occ(0) == 0

def test_configuration_occupation(simple_cubic_binary_prim):
    prim = config.Prim(simple_cubic_binary_prim)
    T = np.array([
        [4, 0, 0],
        [0, 4, 0],
        [0, 0, 4],
    ])
    supercell = config.make_canonical_supercell(config.Supercell(prim, T))
    configuration = config.Configuration(supercell)

    occupation = configuration.occupation()
    assert configuration.occ(0) == 0
    assert (occupation == np.array([0] * 64)).all()

    configuration.set_occ(0, 1)
    occupation = configuration.occupation()
    assert configuration.occ(0) == 1
    assert (occupation == np.array([1] + [0] * 63)).all()


def test_canonical_configuration_occupation(simple_cubic_binary_prim):
    prim = config.Prim(simple_cubic_binary_prim)
    T = np.array([
        [4, 0, 0],
        [0, 4, 0],
        [0, 0, 4],
    ])
    supercell = config.make_canonical_supercell(config.Supercell(prim, T))
    configuration = config.Configuration(supercell)

    configuration.set_occ(0, 1)
    assert config.is_canonical_configuration(configuration) == True

    configuration.set_occ(0, 0)
    configuration.set_occ(10, 1)
    canon_config = config.make_canonical_configuration(configuration)
    assert config.is_canonical_configuration(canon_config) == True
    assert (canon_config.occupation() == np.array([1] + [0] * 63)).all()


def test_configuration_invariant_subgroup(simple_cubic_binary_prim):
    prim = config.Prim(simple_cubic_binary_prim)
    T = np.array([
        [4, 0, 0],
        [0, 4, 0],
        [0, 0, 4],
    ])
    supercell = config.make_canonical_supercell(config.Supercell(prim, T))
    configuration = config.Configuration(supercell)

    invariant_subgroup = config.make_invariant_subgroup(configuration)
    assert len(invariant_subgroup) == 64 * 48

    configuration.set_occ(0, 1)
    invariant_subgroup = config.make_invariant_subgroup(configuration)
    assert len(invariant_subgroup) == 48


def test_configuration_apply(simple_cubic_binary_prim):
    prim = config.Prim(simple_cubic_binary_prim)
    T = np.array([
        [4, 0, 0],
        [0, 4, 0],
        [0, 0, 4],
    ])
    supercell = config.make_canonical_supercell(config.Supercell(prim, T))
    configuration = config.Configuration(supercell)
    configuration.set_occ(l=0, s=1)
    configuration.set_occ(l=1, s=1)

    rep = config.SupercellSymOp.begin(supercell)
    end = config.SupercellSymOp.end(supercell)
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

def test_copy_configuration(simple_cubic_binary_prim):
    import copy
    prim = config.Prim(simple_cubic_binary_prim)
    T = np.array([
        [4, 0, 0],
        [0, 4, 0],
        [0, 0, 4],
    ])
    supercell = config.make_canonical_supercell(config.Supercell(prim, T))
    configuration1 = config.Configuration(supercell)
    configuration2 = copy.copy(configuration1)

    assert isinstance(configuration1, config.Configuration)
    assert isinstance(configuration2, config.Configuration)
    assert configuration2 is not configuration1
