import numpy as np
import pytest

import libcasm.xtal as xtal
import libcasm.xtal.prims as xtal_prims
from libcasm.occ_events import OccEvent, OccPosition


@pytest.fixture
def fcc_1NN_A_Va_event():
    xtal_prim = xtal_prims.FCC(r=1.0, occ_dof=["A", "B", "Va"])

    site1 = xtal.IntegralSiteCoordinate(sublattice=0, unitcell=[0, 0, 0])
    site2 = xtal.IntegralSiteCoordinate(sublattice=0, unitcell=[1, 0, 0])

    A_occ_index = 0
    Va_occ_index = 2

    A_initial_pos = OccPosition.molecule(site1, A_occ_index)
    A_final_pos = OccPosition.molecule(site2, A_occ_index)
    Va_initial_pos = OccPosition.molecule(site2, Va_occ_index)
    Va_final_pos = OccPosition.molecule(site1, Va_occ_index)

    return (
        xtal_prim,
        OccEvent([[A_initial_pos, A_final_pos], [Va_initial_pos, Va_final_pos]]),
    )


@pytest.fixture
def lowsym_Hstrain_prim():
    return xtal.Prim(
        lattice=xtal.Lattice(
            np.array(
                [
                    [1.0, 0.3, 0.4],  # a
                    [0.0, 1.2, 0.5],  # b
                    [0.0, 0.0, 1.4],  # c
                ]
            ).transpose()
        ),
        coordinate_frac=np.array(
            [
                [0.0, 0.0, 0.0],
                [0.4, 0.5, 0.6],
                [0.24, 0.25, 0.23],
            ]
        ).transpose(),
        occ_dof=[["A"], ["A"], ["A"]],
        global_dof=[xtal.DoFSetBasis("Hstrain")],
    )


@pytest.fixture
def lowsym_disp_prim():
    return xtal.Prim(
        lattice=xtal.Lattice(
            np.array(
                [
                    [1.0, 0.3, 0.4],  # a
                    [0.0, 1.2, 0.5],  # b
                    [0.0, 0.0, 1.4],  # c
                ]
            ).transpose()
        ),
        coordinate_frac=np.array(
            [
                [0.0, 0.0, 0.0],
                [0.4, 0.5, 0.6],
                [0.24, 0.25, 0.23],
            ]
        ).transpose(),
        occ_dof=[["A"], ["A"], ["A"]],
        local_dof=[
            [xtal.DoFSetBasis("disp")],
            [xtal.DoFSetBasis("disp")],
            [xtal.DoFSetBasis("disp")],
        ],
    )


@pytest.fixture
def lowsym_occ_prim():
    return xtal.Prim(
        lattice=xtal.Lattice(
            np.array(
                [
                    [1.0, 0.3, 0.4],  # a
                    [0.0, 1.2, 0.5],  # b
                    [0.0, 0.0, 1.4],  # c
                ]
            ).transpose()
        ),
        coordinate_frac=np.array(
            [
                [0.0, 0.0, 0.0],
                [0.4, 0.5, 0.6],
                [0.24, 0.25, 0.23],
            ]
        ).transpose(),
        occ_dof=[["A", "B"], ["A", "B"], ["A", "B"]],
    )


@pytest.fixture
def lowsym_Hstrain_disp_prim():
    return xtal.Prim(
        lattice=xtal.Lattice(
            np.array(
                [
                    [1.0, 0.3, 0.4],  # a
                    [0.0, 1.2, 0.5],  # b
                    [0.0, 0.0, 1.4],  # c
                ]
            ).transpose()
        ),
        coordinate_frac=np.array(
            [
                [0.0, 0.0, 0.0],
                [0.4, 0.5, 0.6],
                [0.24, 0.25, 0.23],
            ]
        ).transpose(),
        occ_dof=[["A"], ["A"], ["A"]],
        local_dof=[
            [xtal.DoFSetBasis("disp")],
            [xtal.DoFSetBasis("disp")],
            [xtal.DoFSetBasis("disp")],
        ],
        global_dof=[xtal.DoFSetBasis("Hstrain")],
    )
