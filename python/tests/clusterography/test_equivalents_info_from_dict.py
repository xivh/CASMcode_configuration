import json

import libcasm.clusterography as casmclust
import libcasm.xtal as xtal


def test_equivalents_info_from_dict_ZrO(shared_datadir):
    with open(shared_datadir / "ZrO" / "equivalents_info.json", "r") as f:
        equivalents_info = json.load(f)

    xtal_prim = xtal.Prim.from_dict(equivalents_info["prototype"]["prim"])

    (
        phenomenal_clusters,
        equivalent_generating_op_indices,
    ) = casmclust.equivalents_info_from_dict(
        data=equivalents_info,
        xtal_prim=xtal_prim,
    )

    assert len(phenomenal_clusters) == 2
    assert equivalent_generating_op_indices == [0, 2]
