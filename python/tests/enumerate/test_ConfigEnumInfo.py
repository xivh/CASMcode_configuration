import libcasm.configuration as casmconfig
import libcasm.enumerate as casmenum
import libcasm.xtal.prims as xtal_prims


def test_ConfigEnumInfo_1():
    """Checks for running without error. Does not check output format."""

    def filter(configuration: casmconfig.Configuration):
        # A custom filter function
        # return True to keep configuration, False to skip
        return True

    xtal_prim = xtal_prims.FCC(
        r=0.5,
        occ_dof=["A", "B", "C"],
    )
    prim = casmconfig.Prim(xtal_prim)
    supercell_set = casmconfig.SupercellSet(prim=prim)
    configuration_set = casmconfig.ConfigurationSet()
    config_enum = casmenum.ConfigEnumAllOccupations(
        prim=prim,
        supercell_set=supercell_set,
    )
    info = casmenum.ConfigEnumInfo(config_enum, configuration_set)
    for i, configuration in enumerate(
        config_enum.by_supercell(
            supercells={
                "max": 4,
            }
        )
    ):
        info.check()
        info.n_config_total += 1
        if not filter(configuration):
            info.n_config_excluded += 1
        else:
            configuration_set.add(configuration)
        assert isinstance(configuration, casmconfig.Configuration)
    info.finish()

    assert len(configuration_set) == 126
