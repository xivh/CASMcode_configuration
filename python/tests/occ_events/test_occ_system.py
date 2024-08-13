import libcasm.xtal.prims as xtal_prims
from libcasm.occ_events import OccSystem


def test_OccSystem_construction_1():
    prim = xtal_prims.FCC(a=1.0, occ_dof=["A", "B", "Va"])
    system = OccSystem(prim)
    assert isinstance(system, OccSystem)
    assert system.chemical_name_list() == ["A", "B", "Va"]
    assert system.is_vacancy_list() == [False, False, True]
    assert system.orientation_name_list() == ["A", "B", "Va"]


def test_OccSystem_construction_2():
    prim = xtal_prims.FCC(a=1.0, occ_dof=["A", "B", "Va"])
    system = OccSystem(prim, chemical_name_list=["Va", "B", "A"])
    assert isinstance(system, OccSystem)
    assert system.chemical_name_list() == ["Va", "B", "A"]
    assert system.is_vacancy_list() == [True, False, False]
    assert system.orientation_name_list() == ["A", "B", "Va"]


def test_OccSystem_construction_3():
    prim = xtal_prims.FCC(a=1.0, occ_dof=["A", "B", "Va"])
    system = OccSystem(prim, vacancy_name_list=[])
    assert isinstance(system, OccSystem)
    assert system.chemical_name_list() == ["A", "B", "Va"]
    assert system.is_vacancy_list() == [False, False, False]
    assert system.orientation_name_list() == ["A", "B", "Va"]


def test_OccSystem_to_from_dict():
    prim = xtal_prims.FCC(a=1.0, occ_dof=["A", "B", "Va"])
    system = OccSystem(prim)

    data = system.to_dict()
    assert isinstance(data, dict)
    assert data["chemical_name_list"] == ["A", "B", "Va"]
    assert data["is_vacancy_list"] == [False, False, True]
    assert data["orientation_name_list"] == ["A", "B", "Va"]

    system_in = OccSystem.from_dict(data, prim)
    assert isinstance(system_in, OccSystem)
    assert system_in.chemical_name_list() == ["A", "B", "Va"]
    assert system_in.is_vacancy_list() == [False, False, True]
    assert system_in.orientation_name_list() == ["A", "B", "Va"]


def test_OccSystem_repr():
    prim = xtal_prims.FCC(a=1.0, occ_dof=["A", "B", "Va"])
    system = OccSystem(prim)

    import io
    from contextlib import redirect_stdout

    f = io.StringIO()
    with redirect_stdout(f):
        print(system)
    out = f.getvalue()
    assert "atom_name_list" in out
