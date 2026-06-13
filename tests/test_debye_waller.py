"""Tests for the Debye-Waller factor implementation (added 2026-06-13).

Covers code/debye_waller.py (Mode 2: phonon-based BZ sum) and its
integration into SingleQAnalyzer via config.DEBYE_WALLER_MODE.

Run from project root: pytest tests/ -v
"""
import numpy as np
import pytest

from code import config
from code.constants import const
from code.force_constants import ForceConstants
from code.debye_waller import (compute_U_tensors_phonon, compute_U_tensors_cif,
                                debye_waller_factor, _conv_cell_vectors,
                                _reciprocal_magnitudes, _u_cif_to_cart)

AUTE2_FC = "data/AuTe_2_m.fc"


@pytest.fixture(scope="module")
def aute2():
    xtal = ForceConstants(AUTE2_FC)
    Phi = xtal.convert_to_eV_per_Angstrom2()
    masses = np.array([xtal.masses[xtal.atom_type_map[i] - 1]
                       for i in range(xtal.nat)])
    return xtal, Phi, masses


@pytest.fixture
def restore_dw_config():
    """Restore config.DEBYE_WALLER_MODE/QMESH after a test that mutates them."""
    saved_mode = config.DEBYE_WALLER_MODE
    saved_mesh = config.DEBYE_WALLER_QMESH
    try:
        yield
    finally:
        config.DEBYE_WALLER_MODE = saved_mode
        config.DEBYE_WALLER_QMESH = saved_mesh


def _kT_THz():
    return 207 * const.c * 100 / 1e12


# ---------------------------------------------------------------------------
# compute_U_tensors_phonon
# ---------------------------------------------------------------------------

def test_U_tensors_shape_and_sanity(aute2):
    xtal, Phi, masses = aute2
    kT = _kT_THz()
    U = compute_U_tensors_phonon(xtal, Phi, masses, kT, n_mesh=8)

    assert U.shape == (xtal.nat, 3, 3)
    for k in range(xtal.nat):
        assert np.allclose(U[k], U[k].T)
        evals = np.linalg.eigvalsh(U[k])
        assert np.all(evals > 0)
        U_iso = np.trace(U[k]) / 3
        assert 0.001 < U_iso < 0.1


# ---------------------------------------------------------------------------
# debye_waller_factor convention check
# ---------------------------------------------------------------------------

def test_debye_waller_factor_isotropic_convention():
    Q_cart = np.array([1.5, -0.7, 2.3])
    Q_mag = np.linalg.norm(Q_cart)
    Uiso = 0.02
    U = Uiso * np.eye(3)

    sin_theta_over_lambda = Q_mag / (4 * np.pi)
    expected = np.exp(-8 * np.pi**2 * Uiso * sin_theta_over_lambda**2)

    assert debye_waller_factor(Q_cart, U) == pytest.approx(expected)


# ---------------------------------------------------------------------------
# Integration with SingleQAnalyzer
# ---------------------------------------------------------------------------

def test_dw_mode_none_unchanged(aute2, restore_dw_config):
    from code.single_q_analysis import SingleQAnalyzer
    xtal, Phi, masses = aute2
    kT = _kT_THz()
    config.DEBYE_WALLER_MODE = 'none'

    a = SingleQAnalyzer(xtal, Phi, masses, kT)
    assert a.dw_U is None

    r = a.analyze((0, 0, 2), coords='conventional', print_results=False)
    assert np.allclose(r['debye_waller_factors'], 1.0)
    assert np.allclose(r['form_factors_eff'],
                        [r['form_factors'][s] for s in a.symbols])


def test_dw_factor_scales_ixs_quadratically(aute2):
    """Wiring check: f_eff = f * T (same T for every atom) must scale
    IXS_stokes by exactly T^2, since DW enters |F|^2 quadratically."""
    from code.single_q_analysis import SingleQAnalyzer
    xtal, Phi, masses = aute2
    kT = _kT_THz()

    a = SingleQAnalyzer(xtal, Phi, masses, kT)
    a.dw_U = None
    r_off = a.analyze((0, 0, 2), coords='conventional', print_results=False)

    Uiso = 0.02
    a.dw_U = np.array([Uiso * np.eye(3)] * xtal.nat)
    r_on = a.analyze((0, 0, 2), coords='conventional', print_results=False)

    T = debye_waller_factor(r_off['Q_cart'], Uiso * np.eye(3))
    assert 0 < T < 1
    assert np.allclose(r_on['debye_waller_factors'], T)

    nz = r_off['IXS_stokes'] > 1e-12
    assert np.allclose(r_on['IXS_stokes'][nz], r_off['IXS_stokes'][nz] * T**2)
    assert np.allclose(r_on['IXS_antistokes'][nz], r_off['IXS_antistokes'][nz] * T**2)


def test_dw_mode_phonon_reduces_ixs(aute2, restore_dw_config):
    from code.single_q_analysis import SingleQAnalyzer
    xtal, Phi, masses = aute2
    kT = _kT_THz()

    config.DEBYE_WALLER_MODE = 'none'
    a_off = SingleQAnalyzer(xtal, Phi, masses, kT)
    r_off = a_off.analyze((0, 0, 2), coords='conventional', print_results=False)

    config.DEBYE_WALLER_MODE = 'phonon'
    config.DEBYE_WALLER_QMESH = 8
    a_on = SingleQAnalyzer(xtal, Phi, masses, kT)
    assert a_on.dw_U is not None

    r_on = a_on.analyze((0, 0, 2), coords='conventional', print_results=False)
    dw = r_on['debye_waller_factors']
    assert np.all(dw > 0) and np.all(dw < 1)

    nz = r_off['IXS_stokes'] > 1e-12
    assert np.all(r_on['IXS_stokes'][nz] < r_off['IXS_stokes'][nz])

    # DW factors should be closer to 1 (less suppression) at small |Q| than
    # at a larger |Q| along the same direction.
    r_small = a_on.analyze((0, 0, 0.2), coords='conventional', print_results=False)
    assert np.all(r_small['debye_waller_factors'] > dw)


# ---------------------------------------------------------------------------
# compute_U_tensors_cif (Mode 1: fixed literature ADPs, Reithmayer et al.)
# ---------------------------------------------------------------------------

def test_u_cif_to_cart_orthorhombic_limit():
    """For beta=90 (orthorhombic), U_cart must equal U_cif exactly."""
    A = _conv_cell_vectors(a=7.189, b=4.407, c=5.069, beta=90.0)
    D = np.diag(_reciprocal_magnitudes(A))

    U_cif = np.array([
        [0.0167, 0.0,    0.0017],
        [0.0,    0.0103, 0.0],
        [0.0017, 0.0,    0.0104]])

    U_cart = _u_cif_to_cart(U_cif, A, D)
    assert np.allclose(U_cart, U_cif)


def test_U_tensors_cif_shape_and_values(aute2):
    xtal, _, _ = aute2
    U = compute_U_tensors_cif(xtal)

    assert U.shape == (3, 3, 3)
    for k in range(3):
        assert np.allclose(U[k], U[k].T)
        evals = np.linalg.eigvalsh(U[k])
        assert np.all(evals > 0)

    # Te1 and Te2 are symmetry-equivalent (m_y mirror => U12=U23=0 in the
    # CIF), so they must get the identical tensor.
    assert np.allclose(U[1], U[2])

    U_iso = np.trace(U, axis1=1, axis2=2) / 3
    assert U_iso[0] == pytest.approx(0.01653, abs=1e-4)   # Au
    assert U_iso[1] == pytest.approx(0.01247, abs=1e-4)   # Te


def test_dw_mode_cif_integration(aute2, restore_dw_config):
    from code.single_q_analysis import SingleQAnalyzer
    xtal, Phi, masses = aute2
    kT = _kT_THz()

    config.DEBYE_WALLER_MODE = 'cif'
    a = SingleQAnalyzer(xtal, Phi, masses, kT)
    assert a.dw_U is not None
    assert np.allclose(a.dw_U, compute_U_tensors_cif(xtal))

    r = a.analyze((0, 0, 2), coords='conventional', print_results=False)
    dw = r['debye_waller_factors']
    assert np.all(dw > 0) and np.all(dw < 1)
