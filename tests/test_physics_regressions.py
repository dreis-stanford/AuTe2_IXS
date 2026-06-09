"""Regression tests for the three physics bugs fixed 2026-06-09.

1. .fc mass unit conversion (Rydberg a.u. -> amu, divisor 911.444 not 931.494)
2. Eigenvector reshape: DOF layout is atom-major -> reshape(nat, 3)
3. kT cm^-1 -> THz conversion (the '* 100' that the separator sed turned into '* 80')

Plus general sanity: Gamma-point acoustic modes, Bose Stokes/anti-Stokes ratio,
and conv<->prim round trips.

Run from project root:  pytest tests/ -v
"""
import re
import numpy as np
import pytest

from code.constants import const
from code.force_constants import ForceConstants
from code.phonons import calc_Dq, calc_freq_eig, convert_frequencies
from code.aute2_structure import aute2_conv2prim_k, aute2_prim2conv_k

SI_FC = "data/Test__Silicon_dispersion/Qgrid_888/Cg.fc"
AUTE2_FC = "data/AuTe_2_m.fc"


@pytest.fixture(scope="module")
def si():
    xtal = ForceConstants(SI_FC)
    Phi = xtal.convert_to_eV_per_Angstrom2()
    masses = np.array([xtal.masses[xtal.atom_type_map[i] - 1]
                       for i in range(xtal.nat)])
    return xtal, Phi, masses


@pytest.fixture(scope="module")
def aute2():
    xtal = ForceConstants(AUTE2_FC)
    Phi = xtal.convert_to_eV_per_Angstrom2()
    masses = np.array([xtal.masses[xtal.atom_type_map[i] - 1]
                       for i in range(xtal.nat)])
    return xtal, Phi, masses


# ---------------------------------------------------------------------------
# Bug 3: mass unit conversion
# ---------------------------------------------------------------------------

def test_si_mass_is_amu(si):
    xtal, _, masses = si
    # Standard atomic weight of Si = 28.0855 amu; file stores it in Ry a.u.
    assert masses[0] == pytest.approx(28.0855, rel=1e-4)


def test_aute2_masses_are_amu(aute2):
    xtal, _, masses = aute2
    by_symbol = {xtal.symbols[xtal.atom_type_map[i] - 1]: masses[i]
                 for i in range(xtal.nat)}
    assert by_symbol['Au'] == pytest.approx(196.97, rel=1e-3)
    assert by_symbol['Te'] == pytest.approx(127.60, rel=1e-3)


# ---------------------------------------------------------------------------
# Sanity: acoustic modes at Gamma
# ---------------------------------------------------------------------------

def test_gamma_acoustic_modes_near_zero(si):
    xtal, Phi, masses = si
    D = calc_Dq(np.zeros((1, 3)), xtal.uvw, Phi, masses)
    w_raw, _ = calc_freq_eig(D, sort_modes=False, rotate_ev=False)
    w_THz = convert_frequencies(w_raw, 'THz')[0]
    # 3 acoustic modes ~0 (ASR enforced), optical modes well above
    assert np.all(np.abs(sorted(np.abs(w_THz))[:3]) < 0.05)
    assert np.all(np.array(sorted(np.abs(w_THz))[3:]) > 10.0)


# ---------------------------------------------------------------------------
# Bug 2: L-char eigenvector layout (the known Si [100] issue)
# ---------------------------------------------------------------------------

def _lchar(ev_mode, Q_hat, masses, nat):
    """L-char exactly as implemented in single_q_analysis*.py (post-fix)."""
    ev_reshaped = ev_mode.reshape(nat, 3)
    tot, lon = 0.0, 0.0
    for iat in range(nat):
        ephys = ev_reshaped[iat, :] / np.sqrt(masses[iat])
        tot += np.real(np.vdot(ephys, ephys))
        lon += np.abs(np.dot(Q_hat, ephys)) ** 2
    return lon / tot


def test_si_lchar_clean_split(si):
    """Along a pure direction, Si modes must be cleanly L (1.0) or T (0.0)."""
    xtal, Phi, masses = si
    q = np.array([[0.1, 0.0, 0.0]])  # primitive rlu
    D = calc_Dq(q, xtal.uvw, Phi, masses)
    w_raw, ev = calc_freq_eig(D, sort_modes=False, rotate_ev=True)

    Q_cart = (q @ xtal.b_l)[0]
    Q_hat = Q_cart / np.linalg.norm(Q_cart)

    L = sorted(_lchar(ev[:, m, 0], Q_hat, masses, xtal.nat)
               for m in range(3 * xtal.nat))
    # 4 transverse (TA x2, TO x2) ~ 0; 2 longitudinal (LA, LO) ~ 1
    assert np.allclose(L[:4], 0.0, atol=0.02), f"T modes not clean: {L}"
    assert np.allclose(L[4:], 1.0, atol=0.02), f"L modes not clean: {L}"


def test_lchar_synthetic_layout():
    """A constructed pure-longitudinal mode must give L-char = 1."""
    nat = 2
    masses = np.array([1.0, 2.0])
    Q_hat = np.array([1.0, 0.0, 0.0])
    # atom-major: [x1,y1,z1, x2,y2,z2], displacement along x only
    ev_mode = np.array([1.0, 0.0, 0.0, 0.5, 0.0, 0.0], dtype=complex)
    assert _lchar(ev_mode, Q_hat, masses, nat) == pytest.approx(1.0)
    # and a pure transverse mode gives 0
    ev_t = np.array([0.0, 1.0, 0.0, 0.0, 0.5, 0.0], dtype=complex)
    assert _lchar(ev_t, Q_hat, masses, nat) == pytest.approx(0.0)


def test_complex_amplitude_uses_vdot():
    """Guard against |dot(e,e)| regression for complex eigenvectors."""
    z = np.array([1 + 1j, 1 - 1j, 0.0])
    assert np.real(np.vdot(z, z)) == pytest.approx(4.0)
    # The old buggy expression gives a different number for this vector:
    assert abs(np.dot(z, z)) != pytest.approx(4.0)


# ---------------------------------------------------------------------------
# Bug 1: kT conversion constant (sed corruption guard)
# ---------------------------------------------------------------------------

def test_kt_conversion_in_source():
    """The cm^-1 -> THz conversion must use c * 100 (m/s -> cm/s), not * 80."""
    for fname in ("code/single_q_analysis.py", "code/single_q_analysis_si.py"):
        with open(fname) as f:
            src = f.read()
        for m in re.finditer(r"kT_THz\s*=\s*kT_cm\s*\*\s*const\.c\s*\*\s*(\d+)", src):
            assert m.group(1) == "100", f"{fname}: kT uses * {m.group(1)}, expected * 100"


def test_kt_value():
    kT_THz = 207 * const.c * 100 / 1e12
    # 207 cm^-1 = 6.205 THz = kT at ~297.8 K
    assert kT_THz == pytest.approx(6.2057, rel=1e-3)
    kT_K = kT_THz * const.THz2meV / 1000 / const.kb
    assert kT_K == pytest.approx(297.8, rel=2e-3)


# ---------------------------------------------------------------------------
# Sanity: Bose factor and Stokes/anti-Stokes ratio
# ---------------------------------------------------------------------------

def test_stokes_antistokes_ratio(si):
    from code.ixs import calc_ixs
    xtal, Phi, masses = si
    q = np.array([[2.1, 0.0, 0.0]])  # away from Gamma, nonzero structure factor
    D = calc_Dq(np.array([[0.1, 0.0, 0.0]]), xtal.uvw, Phi, masses)
    w_raw, ev = calc_freq_eig(D, sort_modes=False, rotate_ev=True)
    w_THz = convert_frequencies(w_raw, 'THz')

    kT_THz = 207 * const.c * 100 / 1e12
    fQ = np.full((1, xtal.nat), 10.0)
    Is, Ias, n, F, info = calc_ixs(w_THz, ev, q, xtal.b_l, xtal.xs, fQ,
                                   masses, kT_THz)
    nz = Is.flatten() > 1e-12
    ratio = Is.flatten()[nz] / Ias.flatten()[nz]
    expected = (n.flatten()[nz] + 1) / n.flatten()[nz]
    assert np.allclose(ratio, expected, rtol=1e-10)


# ---------------------------------------------------------------------------
# Eigenvector Q.e table (verified 2026-06-09)
# ---------------------------------------------------------------------------

def test_qe_projections_si(si):
    """Pin the verified per-atom Q.epsilon values for Si at on-grid q.

    At q=(0.25,0,0) primitive (on the 4x4x4 FC grid, interpolation exact),
    Si LA/LO modes have equal |Q.e| on both atoms (1/sqrt(2)) and relative
    phases of ~+48 deg (LA) / ~-132 deg (LO): the uniform per-(111)-plane
    phase advance, plus 180 for the optic mode. Transverse modes have
    |Q.e| ~ 0 and their phases must be blanked (NaN).
    """
    from code.single_q_analysis import SingleQAnalyzer
    xtal, Phi, masses = si
    q = np.array([[0.25, 0.0, 0.0]])
    D = calc_Dq(q, xtal.uvw, Phi, masses)
    w_raw, ev = calc_freq_eig(D, sort_modes=False, rotate_ev=True)
    w_THz = convert_frequencies(w_raw, 'THz')[0]

    Q_cart = (q @ xtal.b_l)[0]
    Q_hat = Q_cart / np.linalg.norm(Q_cart)

    order = np.argsort(w_THz)
    la, lo = order[2], order[3]  # LA and LO are the longitudinal pair here

    mags_la, ph_la = SingleQAnalyzer.qe_projections(ev[:, la, 0], Q_hat, xtal.nat)
    mags_lo, ph_lo = SingleQAnalyzer.qe_projections(ev[:, lo, 0], Q_hat, xtal.nat)

    # Equal amplitudes 1/sqrt(2) on both atoms
    assert np.allclose(mags_la, 1/np.sqrt(2), atol=0.02)
    assert np.allclose(mags_lo, 1/np.sqrt(2), atol=0.02)

    # Relative phase: reference atom = 0 deg; other = per-plane advance
    rel_la = ph_la[np.isfinite(ph_la) & (ph_la != 0)]
    rel_lo = ph_lo[np.isfinite(ph_lo) & (ph_lo != 0)]
    assert abs(abs(rel_la[0]) - 47.8) < 2.0
    assert abs(abs(rel_lo[0]) - 132.2) < 2.0
    # LA and LO differ by ~180
    assert abs(abs(rel_la[0]) + abs(rel_lo[0]) - 180.0) < 1.0

    # Transverse modes: near-zero |Q.e|, phases blanked as NaN
    for t in (order[0], order[1]):
        mags_t, ph_t = SingleQAnalyzer.qe_projections(ev[:, t, 0], Q_hat, xtal.nat)
        assert np.all(mags_t < 0.05)
        assert np.isnan(ph_t[1 - int(np.argmax(mags_t))])


def test_qe_projections_reconstruct_F(si):
    """f/sqrt(M) * e^{-2pi i Q.r} * (Q.e) summed over atoms must equal calc_ixs F."""
    from code.ixs import calc_ixs
    xtal, Phi, masses = si
    Q_full = np.array([2.1, 0.0, 0.0])
    q_red = np.array([[0.1, 0.0, 0.0]])
    D = calc_Dq(q_red, xtal.uvw, Phi, masses)
    w_raw, ev = calc_freq_eig(D, sort_modes=False, rotate_ev=True)
    w_THz = convert_frequencies(w_raw, 'THz')

    Q_cart = Q_full @ xtal.b_l
    fQ = np.full((1, xtal.nat), 10.0)
    kT = 207 * const.c * 100 / 1e12
    _, _, _, F, _ = calc_ixs(w_THz, ev, Q_full.reshape(1, 3), xtal.b_l,
                             xtal.xs, fQ, masses, kT)

    for m in range(3 * xtal.nat):
        e_m = ev[:, m, 0].reshape(xtal.nat, 3)
        F_rec = sum(fQ[0, s] / np.sqrt(masses[s])
                    * np.exp(-2j * np.pi * np.dot(Q_full, xtal.xs[s]))
                    * np.dot(Q_cart, e_m[s])
                    for s in range(xtal.nat))
        assert abs(F_rec - F[0, m]) < 1e-10 * max(1.0, abs(F[0, m]))


# ---------------------------------------------------------------------------
# Sanity: coordinate transforms
# ---------------------------------------------------------------------------

def test_conv_prim_round_trip():
    rng = np.random.default_rng(0)
    for _ in range(10):
        q = rng.uniform(-3, 3, 3)
        assert np.allclose(aute2_conv2prim_k(aute2_prim2conv_k(q)), q)
        assert np.allclose(aute2_prim2conv_k(aute2_conv2prim_k(q)), q)
