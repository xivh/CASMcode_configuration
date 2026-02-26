import numpy as np

import libcasm.irreps as casmirreps


def clean(arr):
    arr = np.where(np.abs(arr) < 1e-5, 0.0, arr)
    return arr


def test_IrrepInfo_constructor():
    """Test constructing an IrrepInfo directly."""
    trans_mat = np.array([[1.0 + 0j, 0.0 + 0j, 0.0 + 0j]])
    characters = np.array([1.0 + 0j, 1.0 + 0j, 1.0 + 0j])
    irrep_info = casmirreps.IrrepInfo(trans_mat=trans_mat, characters=characters)

    assert isinstance(irrep_info, casmirreps.IrrepInfo)
    assert irrep_info.irrep_dim == 1
    assert irrep_info.vector_dim == 3
    assert irrep_info.index is None
    assert irrep_info.directions is None


def test_IrrepInfo_attributes(FCC_binary_irrep_decomposition):
    """Test IrrepInfo attributes from an irrep decomposition."""
    irreps = FCC_binary_irrep_decomposition.irreps

    assert len(irreps) == 2

    for irrep_info in irreps:
        assert isinstance(irrep_info, casmirreps.IrrepInfo)

        # Check irrep_dim and vector_dim
        assert irrep_info.vector_dim == 8
        assert irrep_info.irrep_dim > 0

        # Check trans_mat shape
        assert irrep_info.trans_mat.shape == (
            irrep_info.irrep_dim,
            irrep_info.vector_dim,
        )

        # Check characters length matches number of group operations
        assert len(irrep_info.characters) == len(
            FCC_binary_irrep_decomposition.matrix_rep
        )


def test_IrrepInfo_identity_irrep(FCC_binary_irrep_decomposition):
    """Test the identity irrep properties."""
    irreps = FCC_binary_irrep_decomposition.irreps

    # One of the irreps should be the identity irrep (1-dimensional, all chars = 1)
    identity_irreps = [ir for ir in irreps if ir.is_identity]
    assert len(identity_irreps) == 1

    identity = identity_irreps[0]
    assert identity.irrep_dim == 1
    assert np.allclose(np.abs(identity.characters), np.ones(len(identity.characters)))


def frobenius_schur_indicator_assertions(irrep_info):
    """Helper function to check Frobenius-Schur indicator consistency."""
    # Frobenius-Schur indicator should be -1, 0, or 1
    assert irrep_info.frobenius_schur_indicator in (-1, 0, 1)

    # Check consistency with is_real, is_complex_irrep, is_pseudo_real
    if irrep_info.frobenius_schur_indicator == 1:
        assert irrep_info.is_real is True
        assert irrep_info.is_complex_irrep is False
        assert irrep_info.is_pseudo_real is False
    elif irrep_info.frobenius_schur_indicator == 0:
        assert irrep_info.is_real is False
        assert irrep_info.is_complex_irrep is True
        assert irrep_info.is_pseudo_real is False
    elif irrep_info.frobenius_schur_indicator == -1:
        assert irrep_info.is_real is False
        assert irrep_info.is_complex_irrep is False
        assert irrep_info.is_pseudo_real is True


def print_summary(irreps):
    for i, x in enumerate(irreps):
        print(f"Irrep {i}: dim={x.irrep_dim}, type={x.irrep_type}, index={x.index}")
    print()
    # print("Symmetry-adapted basis:")
    # B = clean(irrep_decomposition.symmetry_adapted_subspace)
    # for i in range(B.shape[1]):
    #     print(f"{i}: {B[:, i]}")
    # print()
    # assert False


def test_IrrepInfo_1(FCC_binary_irrep_decomposition):
    """Test the Frobenius-Schur indicator."""
    irrep_decomposition = FCC_binary_irrep_decomposition
    irreps = irrep_decomposition.irreps

    # print_summary(irreps)

    for irrep_info in irreps:
        frobenius_schur_indicator_assertions(irrep_info)


# def test_IrrepInfo_2(FCC_disp_vol32_irrep_decomposition):
#     """Test the Frobenius-Schur indicator."""
#     irrep_decomposition = FCC_disp_vol32_irrep_decomposition
#     irreps = irrep_decomposition.irreps
#
#     # print_summary(irreps)
#
#     for irrep_info in irreps:
#         frobenius_schur_indicator_assertions(irrep_info)


# ABC2_disp_irrep_decomposition
def test_IrrepInfo_3(ABC2_disp_irrep_decomposition):
    """Test the Frobenius-Schur indicator."""
    irrep_decomposition = ABC2_disp_irrep_decomposition
    irreps = irrep_decomposition.irreps

    # print_summary(irreps)

    for irrep_info in irreps:
        frobenius_schur_indicator_assertions(irrep_info)


def test_IrrepInfo_4(TlZn2Sb2_disp_irrep_decomposition):
    """Test the Frobenius-Schur indicator."""
    irrep_decomposition = TlZn2Sb2_disp_irrep_decomposition
    irreps = irrep_decomposition.irreps

    # print_summary(irreps)

    for irrep_info in irreps:
        frobenius_schur_indicator_assertions(irrep_info)


def test_IrrepInfo_is_gerade(FCC_binary_irrep_decomposition):
    """Test the is_gerade property."""
    irreps = FCC_binary_irrep_decomposition.irreps

    for irrep_info in irreps:
        assert isinstance(irrep_info.is_gerade, bool)


def test_IrrepInfo_directions(FCC_binary_irrep_decomposition):
    """Test the directions attribute."""
    irreps = FCC_binary_irrep_decomposition.irreps

    for irrep_info in irreps:
        directions = irrep_info.directions
        assert directions is not None
        assert isinstance(directions, list)
        for orbit in directions:
            assert isinstance(orbit, list)
            for direction in orbit:
                assert isinstance(direction, np.ndarray)
                assert len(direction) == irrep_info.vector_dim


def test_IrrepInfo_index(FCC_binary_irrep_decomposition):
    """Test the index attribute."""
    irreps = FCC_binary_irrep_decomposition.irreps

    for irrep_info in irreps:
        assert irrep_info.index is not None
        assert isinstance(irrep_info.index, int)
        assert irrep_info.index >= 0


def test_IrrepInfo_comparison(FCC_binary_irrep_decomposition):
    """Test comparison operators."""
    irreps = FCC_binary_irrep_decomposition.irreps

    assert len(irreps) == 2

    # Test equality with self
    assert irreps[0] == irreps[0]
    assert irreps[1] == irreps[1]

    # Test inequality between different irreps
    assert irreps[0] != irreps[1]

    # Test less-than (identity should sort first)
    identity_idx = 0 if irreps[0].is_identity else 1
    other_idx = 1 - identity_idx
    assert irreps[identity_idx] < irreps[other_idx]


def test_IrrepInfo_to_dict(FCC_binary_irrep_decomposition):
    """Test to_dict serialization."""
    irreps = FCC_binary_irrep_decomposition.irreps

    for irrep_info in irreps:
        data = irrep_info.to_dict()
        assert isinstance(data, dict)


def test_IrrepInfo_from_dict_roundtrip(FCC_binary_irrep_decomposition):
    """Test from_dict / to_dict roundtrip."""
    irreps = FCC_binary_irrep_decomposition.irreps

    for irrep_info in irreps:
        data = irrep_info.to_dict()
        reconstructed = casmirreps.IrrepInfo.from_dict(data)

        assert isinstance(reconstructed, casmirreps.IrrepInfo)
        assert reconstructed.irrep_dim == irrep_info.irrep_dim
        assert reconstructed.vector_dim == irrep_info.vector_dim
        assert np.allclose(reconstructed.trans_mat, irrep_info.trans_mat)
        assert np.allclose(reconstructed.characters, irrep_info.characters)
        assert (
            reconstructed.frobenius_schur_indicator
            == irrep_info.frobenius_schur_indicator
        )
        assert reconstructed.index == irrep_info.index
