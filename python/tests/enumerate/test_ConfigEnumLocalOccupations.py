import libcasm.configuration as casmconfig
import libcasm.xtal.prims as xtal_prims


def test_ConfigEnumLocalOccupations_1():
    xtal_prim = xtal_prims.FCC(a=1.0, occ_dof=["A", "B", "Va"])
    prim = casmconfig.Prim(xtal_prim)

    assert prim is not None
