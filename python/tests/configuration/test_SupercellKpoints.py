"""Tests for SupercellKpoints, SupercellDoF, DiscreteFourierTransform, and standalone
functions make_plane_wave_basis, make_kpoint_irreps, and make_unique_kpoint_irreps."""

import numpy as np
import pytest

import libcasm.configuration as casmconfig


def _check_kpoint_orbits(kpoints, kpoint_orbits):
    """Check that k-point orbits partition all k-points exactly once."""
    n_kpoints = kpoints.shape[1]
    covered = [False] * n_kpoints
    for orbit in kpoint_orbits:
        for idx in orbit:
            assert 0 <= idx < n_kpoints, f"k-point index {idx} out of range"
            assert not covered[idx], f"k-point {idx} appears in multiple orbits"
            covered[idx] = True
    assert all(covered), "Not all k-points are covered by orbits"


def _check_orthonormal(B, tol=1e-8):
    """Check that columns of B are orthonormal."""
    dim = B.shape[1]
    assert np.allclose(
        B.T @ B, np.eye(dim), atol=tol
    ), f"Basis is not orthonormal:\n{B.T @ B}"


class TestBCC2x2x2Disp:
    """Tests for SupercellKpoints and SupercellDoF with 2x2x2 BCC supercell and disp
    DoF.

    BCC primitive cell has 1 basis site with 3 disp components.
    The 2x2x2 supercell has 8 unit cells and 24 total disp DoF.
    """

    @pytest.fixture
    def sc_kpts(self, BCC_binary_GLstrain_disp_prim):
        prim = casmconfig.Prim(xtal_prim=BCC_binary_GLstrain_disp_prim)
        T = np.eye(3, dtype=int) * 2
        supercell = casmconfig.Supercell(prim, T)
        return casmconfig.SupercellKpoints(supercell=supercell)

    @pytest.fixture
    def sc_dof(self, BCC_binary_GLstrain_disp_prim):
        prim = casmconfig.Prim(xtal_prim=BCC_binary_GLstrain_disp_prim)
        T = np.eye(3, dtype=int) * 2
        supercell = casmconfig.Supercell(prim, T)
        return casmconfig.SupercellDoF(supercell=supercell, dof_key="disp")

    def test_construction(self, sc_kpts, sc_dof):
        assert sc_kpts.supercell is not None
        assert sc_kpts.coordinates is not None
        assert sc_kpts.orbits is not None
        assert sc_kpts.little_groups is not None
        assert sc_dof.dof_key == "disp"
        assert sc_dof.dof_space is not None
        assert sc_dof.dof_id is not None

    def test_kpoints_shape(self, sc_kpts):
        # 2x2x2 = 8 k-points
        assert sc_kpts.coordinates.shape == (3, 8)
        assert sc_kpts.indices.shape == (3, 8)

    def test_dof_space(self, sc_dof):
        # 1 prim site * 3 disp components * 8 unitcells = 24 total DoF
        assert sc_dof.dof_space.basis.shape == (24, 24)
        # 1 sublattice * 3 disp components = 3 unique dof IDs
        assert len(sc_dof.dof_id) == 3
        assert len(sc_dof.axis_dof_id) == 24

    def test_kpoint_orbits(self, sc_kpts):
        _check_kpoint_orbits(sc_kpts.coordinates, sc_kpts.orbits)
        # BCC Oh symmetry creates 3 orbits for 2x2x2:
        # Gamma (1), 6 zone-face centers, R-point (1)
        assert len(sc_kpts.orbits) == 3
        assert sc_kpts.orbits == [[0], [1, 2, 3, 4, 5, 6], [7]]

    def test_little_groups(self, sc_kpts):
        # Each k-point has a non-empty little group
        assert len(sc_kpts.little_groups) == 8
        for lg in sc_kpts.little_groups:
            assert len(lg) > 0

    def test_make_plane_wave_basis(self, sc_kpts, sc_dof):
        # Test make_plane_wave_basis for one k-point per orbit
        for orbit in sc_kpts.orbits:
            kpoint_index = orbit[0]
            basis = casmconfig.make_plane_wave_basis(
                supercell_kpoints=sc_kpts,
                supercell_dof=sc_dof,
                kpoint_index=kpoint_index,
                as_complex=False,
            )
            assert basis.ndim == 2
            assert basis.shape[0] == 24
            _check_orthonormal(basis)

    def test_make_plane_wave_basis_complex(self, sc_kpts, sc_dof):
        for orbit in sc_kpts.orbits:
            kpoint_index = orbit[0]
            basis = casmconfig.make_plane_wave_basis(
                supercell_kpoints=sc_kpts,
                supercell_dof=sc_dof,
                kpoint_index=kpoint_index,
                as_complex=True,
            )
            assert basis.ndim == 2
            assert basis.shape[0] == 24
            # Check orthonormality for complex basis
            dim = basis.shape[1]
            assert np.allclose(
                basis.conj().T @ basis, np.eye(dim), atol=1e-8
            ), "Complex Bloch basis is not orthonormal"

    def test_make_kpoint_irreps(self, sc_kpts, sc_dof):
        # Test make_kpoint_irreps for one k-point per orbit
        for orbit in sc_kpts.orbits:
            kpoint_index = orbit[0]
            irrep_decomp, dof_space = casmconfig.make_kpoint_irreps(
                supercell_kpoints=sc_kpts,
                supercell_dof=sc_dof,
                kpoint_index=kpoint_index,
            )
            assert irrep_decomp is not None
            assert dof_space is not None
            assert len(irrep_decomp.irreps) > 0
            # Symmetry-adapted subspace should be orthonormal
            B = irrep_decomp.symmetry_adapted_subspace
            _check_orthonormal(B)

    def test_make_unique_kpoint_irreps(self, sc_kpts, sc_dof):
        kpoint_irreps, dof_space, axis_irrep_info = (
            casmconfig.make_unique_kpoint_irreps(
                supercell_kpoints=sc_kpts,
                supercell_dof=sc_dof,
            )
        )
        B = dof_space.basis
        # Spans the full 24-dimensional disp space
        assert B.shape == (24, 24)
        _check_orthonormal(B)
        assert len(kpoint_irreps) > 0
        assert len(axis_irrep_info) == 24

    def test_axis_irrep_info(self, sc_kpts, sc_dof):
        """Each axis_irrep_info entry has valid orbit_index and kpoint_irreps_index."""
        kpoint_irreps, _, axis_irrep_info = casmconfig.make_unique_kpoint_irreps(
            supercell_kpoints=sc_kpts,
            supercell_dof=sc_dof,
        )
        n_orbits = len(sc_kpts.orbits)
        for info in axis_irrep_info:
            assert 0 <= info.orbit_index < n_orbits
            assert 0 <= info.kpoint_irreps_index < len(kpoint_irreps)
            assert info.irrep_char_index >= 0

    def test_discrete_fourier_transform_construction(self, sc_kpts, sc_dof):
        dft = casmconfig.DiscreteFourierTransform(
            supercell_kpoints=sc_kpts,
            supercell_dof=sc_dof,
        )
        n_unitcells = sc_kpts.supercell.n_unitcells
        n_dof_id = len(sc_dof.dof_id)
        assert dft.dft_M.shape == (n_unitcells, n_unitcells)
        assert dft.idft_M.shape == (n_unitcells, n_unitcells)
        assert dft.dft_phase.shape == (n_dof_id, n_unitcells)
        assert dft.idft_phase.shape == (n_dof_id, n_unitcells)


class TestZrO2x2x2Disp:
    """Tests for SupercellKpoints and SupercellDoF with 2x2x2 ZrO supercell and disp
    DoF.

    ZrO primitive cell has 4 basis sites each with 3 disp components.
    The 2x2x2 supercell has 8 unit cells and 96 total disp DoF.
    """

    @pytest.fixture
    def sc_kpts(self, ZrO_prim_GLstrain_disp):
        prim = casmconfig.Prim(xtal_prim=ZrO_prim_GLstrain_disp)
        T = np.eye(3, dtype=int) * 2
        supercell = casmconfig.Supercell(prim, T)
        return casmconfig.SupercellKpoints(supercell=supercell)

    @pytest.fixture
    def sc_dof(self, ZrO_prim_GLstrain_disp):
        prim = casmconfig.Prim(xtal_prim=ZrO_prim_GLstrain_disp)
        T = np.eye(3, dtype=int) * 2
        supercell = casmconfig.Supercell(prim, T)
        return casmconfig.SupercellDoF(supercell=supercell, dof_key="disp")

    def test_construction(self, sc_kpts, sc_dof):
        assert sc_kpts.supercell is not None
        assert sc_kpts.coordinates is not None
        assert sc_kpts.orbits is not None
        assert sc_kpts.little_groups is not None
        assert sc_dof.dof_key == "disp"
        assert sc_dof.dof_space is not None

    def test_kpoints_shape(self, sc_kpts):
        # 2x2x2 = 8 k-points
        assert sc_kpts.coordinates.shape == (3, 8)
        assert sc_kpts.indices.shape == (3, 8)

    def test_dof_space(self, sc_dof):
        # 4 prim sites * 3 disp components * 8 unitcells = 96 total DoF
        assert sc_dof.dof_space.basis.shape == (96, 96)
        # 4 sublattices * 3 disp components = 12 unique dof IDs
        assert len(sc_dof.dof_id) == 12
        assert len(sc_dof.axis_dof_id) == 96

    def test_kpoint_orbits(self, sc_kpts):
        # All k-points are covered exactly once by orbits
        _check_kpoint_orbits(sc_kpts.coordinates, sc_kpts.orbits)
        # Total k-points across orbits = 8
        total = sum(len(orbit) for orbit in sc_kpts.orbits)
        assert total == 8

    def test_little_groups(self, sc_kpts):
        assert len(sc_kpts.little_groups) == 8
        for lg in sc_kpts.little_groups:
            assert len(lg) > 0

    def test_make_plane_wave_basis(self, sc_kpts, sc_dof):
        for orbit in sc_kpts.orbits:
            kpoint_index = orbit[0]
            basis = casmconfig.make_plane_wave_basis(
                supercell_kpoints=sc_kpts,
                supercell_dof=sc_dof,
                kpoint_index=kpoint_index,
            )
            assert basis.ndim == 2
            assert basis.shape[0] == 96
            _check_orthonormal(basis)

    def test_make_kpoint_irreps(self, sc_kpts, sc_dof):
        for orbit in sc_kpts.orbits:
            kpoint_index = orbit[0]
            irrep_decomp, dof_space = casmconfig.make_kpoint_irreps(
                supercell_kpoints=sc_kpts,
                supercell_dof=sc_dof,
                kpoint_index=kpoint_index,
            )
            assert irrep_decomp is not None
            assert len(irrep_decomp.irreps) > 0
            B = irrep_decomp.symmetry_adapted_subspace
            _check_orthonormal(B)

    def test_make_unique_kpoint_irreps(self, sc_kpts, sc_dof):
        kpoint_irreps, dof_space, axis_irrep_info = (
            casmconfig.make_unique_kpoint_irreps(
                supercell_kpoints=sc_kpts,
                supercell_dof=sc_dof,
            )
        )
        B = dof_space.basis
        # Spans the full 96-dimensional disp space
        assert B.shape == (96, 96)
        _check_orthonormal(B)
        assert len(kpoint_irreps) > 0
        assert len(axis_irrep_info) == 96

    def test_discrete_fourier_transform_construction(self, sc_kpts, sc_dof):
        dft = casmconfig.DiscreteFourierTransform(
            supercell_kpoints=sc_kpts,
            supercell_dof=sc_dof,
        )
        n_unitcells = sc_kpts.supercell.n_unitcells
        n_dof_id = len(sc_dof.dof_id)
        assert dft.dft_M.shape == (n_unitcells, n_unitcells)
        assert dft.dft_phase.shape == (n_dof_id, n_unitcells)


class TestZrO1x2x3Disp:
    """Tests for SupercellKpoints and SupercellDoF with 1x2x3 ZrO supercell and disp
    DoF.

    ZrO primitive cell has 4 basis sites each with 3 disp components.
    The 1x2x3 supercell has 6 unit cells and 72 total disp DoF.
    """

    @pytest.fixture
    def sc_kpts(self, ZrO_prim_GLstrain_disp):
        prim = casmconfig.Prim(xtal_prim=ZrO_prim_GLstrain_disp)
        T = np.diag([1, 2, 3]).astype(int)
        supercell = casmconfig.Supercell(prim, T)
        return casmconfig.SupercellKpoints(supercell=supercell)

    @pytest.fixture
    def sc_dof(self, ZrO_prim_GLstrain_disp):
        prim = casmconfig.Prim(xtal_prim=ZrO_prim_GLstrain_disp)
        T = np.diag([1, 2, 3]).astype(int)
        supercell = casmconfig.Supercell(prim, T)
        return casmconfig.SupercellDoF(supercell=supercell, dof_key="disp")

    def test_construction(self, sc_kpts, sc_dof):
        assert sc_kpts.supercell is not None
        assert sc_kpts.coordinates is not None
        assert sc_kpts.orbits is not None
        assert sc_kpts.little_groups is not None
        assert sc_dof.dof_key == "disp"
        assert sc_dof.dof_space is not None

    def test_kpoints_shape(self, sc_kpts):
        # 1x2x3 = 6 k-points
        assert sc_kpts.coordinates.shape == (3, 6)
        assert sc_kpts.indices.shape == (3, 6)

    def test_dof_space(self, sc_dof):
        # 4 prim sites * 3 disp components * 6 unitcells = 72 total DoF
        assert sc_dof.dof_space.basis.shape == (72, 72)
        # 4 sublattices * 3 disp components = 12 unique dof IDs
        assert len(sc_dof.dof_id) == 12
        assert len(sc_dof.axis_dof_id) == 72

    def test_kpoint_orbits(self, sc_kpts):
        _check_kpoint_orbits(sc_kpts.coordinates, sc_kpts.orbits)
        total = sum(len(orbit) for orbit in sc_kpts.orbits)
        assert total == 6

    def test_little_groups(self, sc_kpts):
        assert len(sc_kpts.little_groups) == 6
        for lg in sc_kpts.little_groups:
            assert len(lg) > 0

    def test_make_plane_wave_basis(self, sc_kpts, sc_dof):
        for orbit in sc_kpts.orbits:
            kpoint_index = orbit[0]
            basis = casmconfig.make_plane_wave_basis(
                supercell_kpoints=sc_kpts,
                supercell_dof=sc_dof,
                kpoint_index=kpoint_index,
            )
            assert basis.ndim == 2
            assert basis.shape[0] == 72
            _check_orthonormal(basis)

    def test_make_kpoint_irreps(self, sc_kpts, sc_dof):
        for orbit in sc_kpts.orbits:
            kpoint_index = orbit[0]
            irrep_decomp, dof_space = casmconfig.make_kpoint_irreps(
                supercell_kpoints=sc_kpts,
                supercell_dof=sc_dof,
                kpoint_index=kpoint_index,
            )
            assert irrep_decomp is not None
            assert len(irrep_decomp.irreps) > 0
            B = irrep_decomp.symmetry_adapted_subspace
            _check_orthonormal(B)

    def test_make_unique_kpoint_irreps(self, sc_kpts, sc_dof):
        kpoint_irreps, dof_space, axis_irrep_info = (
            casmconfig.make_unique_kpoint_irreps(
                supercell_kpoints=sc_kpts,
                supercell_dof=sc_dof,
            )
        )
        B = dof_space.basis
        # Spans the full 72-dimensional disp space
        assert B.shape == (72, 72)
        _check_orthonormal(B)
        assert len(kpoint_irreps) > 0
        assert len(axis_irrep_info) == 72

    def test_discrete_fourier_transform_construction(self, sc_kpts, sc_dof):
        dft = casmconfig.DiscreteFourierTransform(
            supercell_kpoints=sc_kpts,
            supercell_dof=sc_dof,
        )
        n_unitcells = sc_kpts.supercell.n_unitcells
        n_dof_id = len(sc_dof.dof_id)
        assert dft.dft_M.shape == (n_unitcells, n_unitcells)
        assert dft.dft_phase.shape == (n_dof_id, n_unitcells)
