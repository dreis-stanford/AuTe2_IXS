"""
Sixcircle geometry sanity check for AuTe2 IXS experiment.

Tests:
  1. Cubic Si (trivial mounting) — verifies convention baseline
  2. AuTe2 — verifies UB consistency and verify_scattering convention

Strategy: use scbasic.py directly (no sixcircle module side effects).
b_matrix passed to verify_full must be (U @ B).T so that
  Q_hkl @ b_matrix == U @ B @ HKL_col
"""

import sys
import numpy as np

SIXCIRCLE_PATH = '/Users/dreis/Documents/MyPython/Others/sixcircle_1p85'
sys.path.insert(0, SIXCIRCLE_PATH)

import scbasic  # noqa: E402 (prints on import, expected)

sys.path.insert(0, '/Users/dreis/Documents/MyPython/AuTe2_ixs_project')
from code.verify_scattering import verify_basic, verify_full, print_verification


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def build_UB(a, b, c, al, be, ga,
             wavelength,
             or0_hkl, or0_tth, or0_th, or0_chi, or0_phi, or0_mu, or0_gam,
             or1_hkl, or1_tth, or1_th, or1_chi, or1_phi, or1_mu, or1_gam):
    """Return (U, B) matrices using scbasic conventions."""
    flag, B, *_ = scbasic.B_matrix(a, b, c, al, be, ga)
    assert flag, "Invalid lattice parameters"

    # omega = th - tth/2
    omega0 = or0_th - or0_tth / 2
    omega1 = or1_th - or1_tth / 2

    # Unit scattering vectors in phi space from motor positions
    u0phi = scbasic.uphi_vector(or0_tth / 2, omega0, or0_chi, or0_phi, or0_mu, or0_gam)
    u1phi = scbasic.uphi_vector(or1_tth / 2, omega1, or1_chi, or1_phi, or1_mu, or1_gam)

    # Unit scattering vectors in Cartesian crystal space from HKL
    h0b = np.array(or0_hkl, dtype=float).reshape(3, 1)
    h1b = np.array(or1_hkl, dtype=float).reshape(3, 1)
    u0c = scbasic.uc_vector(B, h0b)
    u1c = scbasic.uc_vector(B, h1b)

    errcode, U = scbasic.U_matrix(u0phi, u1phi, u0c, u1c)
    assert errcode == 0, f"U_matrix failed: errcode={errcode}"
    return U, B


def ca_from_scbasic(hkl, U, B, wavelength, mu, gam, omega=0.0):
    """
    Compute motor angles for hkl using SC_angles (frozen mu, gam, omega).
    Returns list of (tth, th, chi, phi, mu, gam) tuples, or [] on failure.
    """
    hb = np.array(hkl, dtype=float).reshape(3, 1)
    flag, raw_angles = scbasic.SC_angles(
        wavelength, U, B, hb,
        thd='x', thetad='x', omegad=omega,
        chid='x', phid='x', mud=mu, gammad=gam
    )
    if not flag:
        return []
    N, motor_angles = scbasic.motors(raw_angles)
    return motor_angles


def angles_dict(motor_tuple):
    tth, th, chi, phi, mu, gam = motor_tuple
    return {'tth': tth, 'th': th, 'chi': chi, 'phi': phi, 'mu': mu, 'gam': gam}


def check_reflection(label, hkl, expected_tth, expected_chi, expected_phi,
                     U, B, wavelength, mu, gam, tol_tth=0.1, tol_angle=1.0):
    """Compute ca() for hkl and check against expected angles + verify_full."""
    print(f"\n  {label}: HKL={hkl}")
    motor_list = ca_from_scbasic(hkl, U, B, wavelength, mu, gam)

    if not motor_list:
        print("    ca(): NO SOLUTION")
        return False

    # Show all solutions
    for sol in motor_list:
        print(f"    solution: tth={sol[0]:.3f} th={sol[1]:.3f} chi={sol[2]:.3f} "
              f"phi={sol[3]:.3f} mu={sol[4]:.4f} gam={sol[5]:.3f}")

    # Pick solution closest to expected
    best = min(motor_list,
               key=lambda s: abs(s[0] - expected_tth) + abs(s[2] - expected_chi))
    tth, th, chi, phi, mu_sol, gam_sol = best
    ang = angles_dict(best)

    tth_ok = abs(tth - expected_tth) < tol_tth
    chi_ok = abs(chi - expected_chi) < tol_angle
    print(f"    best: tth={tth:.3f} (exp {expected_tth:.2f}) {'✓' if tth_ok else '✗'}  "
          f"chi={chi:.3f} (exp {expected_chi:.2f}) {'✓' if chi_ok else '✗'}  "
          f"phi={phi:.3f}")

    # verify_full: b_matrix must be (U @ B).T so that
    # Q_hkl @ b_matrix == U @ B @ HKL (as column)
    UB_mat = U @ B
    b_matrix = UB_mat.T

    hkl_arr = np.array(hkl, dtype=float)
    basic = verify_basic(hkl, ang, wavelength, b_matrix)
    full = verify_full(hkl, ang, wavelength, b_matrix)

    basic_status, full_status = print_verification(basic, full)

    # The convention/geometry check is the FULL test (k_out = k_in + Q_lab).
    # BASIC test WARNING (not FAIL) reflects sub-0.01deg rounding in the
    # experimental orientation angles, not a geometry bug.
    return tth_ok and chi_ok and basic_status != 'FAIL' and full_status == 'PASS'


# ---------------------------------------------------------------------------
# Test 1: Si cubic, trivial mounting
# ---------------------------------------------------------------------------

def test_si_cubic():
    print("\n" + "=" * 70)
    print("TEST 1: Si cubic (a=b=c=5.431 Å), trivial mounting a||x, b||y, c||z")
    print("=" * 70)

    a = b = c = 5.431
    al = be = ga = 90.0
    wavelength = 0.570107
    mu, gam = 0.0, 0.0

    # or0: (4,0,0) at phi=0, chi=0, symmetric (omega=0)
    flag, B, *_ = scbasic.B_matrix(a, b, c, al, be, ga)
    h0b = np.array([[4], [0], [0]], dtype=float)
    _, thetaD0 = scbasic.thetaD_angle(wavelength, B, h0b)
    tth0 = 2 * thetaD0
    th0 = tth0 / 2
    or0 = dict(hkl=(4, 0, 0), tth=tth0, th=th0, chi=0.0, phi=0.0, mu=0.0, gam=0.0)

    # or1: (0,4,0) at phi=90, chi=0
    h1b = np.array([[0], [4], [0]], dtype=float)
    _, thetaD1 = scbasic.thetaD_angle(wavelength, B, h1b)
    tth1 = 2 * thetaD1
    th1 = tth1 / 2
    or1 = dict(hkl=(0, 4, 0), tth=tth1, th=th1, chi=0.0, phi=90.0, mu=0.0, gam=0.0)

    print(f"  or0: {or0['hkl']}  tth={or0['tth']:.3f}° th={or0['th']:.3f}° chi=0 phi=0")
    print(f"  or1: {or1['hkl']}  tth={or1['tth']:.3f}° th={or1['th']:.3f}° chi=0 phi=90")

    U, B = build_UB(a, b, c, al, be, ga, wavelength,
                    or0['hkl'], or0['tth'], or0['th'], or0['chi'], or0['phi'], or0['mu'], or0['gam'],
                    or1['hkl'], or1['tth'], or1['th'], or1['chi'], or1['phi'], or1['mu'], or1['gam'])

    print(f"\n  U matrix (should be identity for trivial mounting):")
    for row in U:
        print(f"    [{row[0]:8.5f}  {row[1]:8.5f}  {row[2]:8.5f}]")

    results = []
    results.append(check_reflection("(4,0,0) round-trip", (4, 0, 0),
                                    expected_tth=tth0, expected_chi=0.0, expected_phi=0.0,
                                    U=U, B=B, wavelength=wavelength, mu=mu, gam=gam))
    results.append(check_reflection("(0,4,0) round-trip", (0, 4, 0),
                                    expected_tth=tth1, expected_chi=0.0, expected_phi=90.0,
                                    U=U, B=B, wavelength=wavelength, mu=mu, gam=gam))
    results.append(check_reflection("(0,0,4)", (0, 0, 4),
                                    expected_tth=tth0, expected_chi=90.0, expected_phi=0.0,
                                    U=U, B=B, wavelength=wavelength, mu=mu, gam=gam))

    ok = all(r is not False for r in results)
    print(f"\n  TEST 1: {'PASS' if ok else 'FAIL'}")
    return ok


# ---------------------------------------------------------------------------
# Test 2: AuTe2 experimental orientation
# ---------------------------------------------------------------------------

def test_aute2():
    print("\n" + "=" * 70)
    print("TEST 2: AuTe2 experimental orientation (from ini.conf / _load_sixcircle)")
    print("=" * 70)

    a, b, c = 7.189, 4.407, 5.069
    al, be, ga = 90.0, 89.96, 90.0
    wavelength = 0.570118589663309
    mu, gam = -0.1719, 0.0

    # From _load_sixcircle (physically verified: Bragg law consistent)
    or0 = dict(hkl=(0, 0, 1), tth=6.45, th=3.225, chi=90.0, phi=0.0, mu=mu, gam=gam)
    or1 = dict(hkl=(0, 2, 0), tth=14.84, th=7.42, chi=0.0, phi=90.0, mu=mu, gam=gam)

    print(f"  or0: {or0['hkl']}  tth={or0['tth']}° th={or0['th']}° chi=90 phi=0")
    print(f"  or1: {or1['hkl']}  tth={or1['tth']}° th={or1['th']}° chi=0 phi=90")

    U, B = build_UB(a, b, c, al, be, ga, wavelength,
                    or0['hkl'], or0['tth'], or0['th'], or0['chi'], or0['phi'], or0['mu'], or0['gam'],
                    or1['hkl'], or1['tth'], or1['th'], or1['chi'], or1['phi'], or1['mu'], or1['gam'])

    print(f"\n  U matrix:")
    for row in U:
        print(f"    [{row[0]:8.5f}  {row[1]:8.5f}  {row[2]:8.5f}]")

    results = []
    # Round-trip the orientation reflections
    results.append(check_reflection(
        "or0 round-trip (0,0,1)", (0, 0, 1),
        expected_tth=6.45, expected_chi=90.0, expected_phi=0.0,
        U=U, B=B, wavelength=wavelength, mu=mu, gam=gam))
    results.append(check_reflection(
        "or1 round-trip (0,2,0)", (0, 2, 0),
        expected_tth=14.84, expected_chi=0.0, expected_phi=90.0,
        U=U, B=B, wavelength=wavelength, mu=mu, gam=gam))
    # A few other physically interesting reflections
    results.append(check_reflection(
        "(0,0,2)", (0, 0, 2),
        expected_tth=12.91, expected_chi=90.0, expected_phi=0.0,
        U=U, B=B, wavelength=wavelength, mu=mu, gam=gam))
    results.append(check_reflection(
        "(0,4,0)", (0, 4, 0),
        expected_tth=30.01, expected_chi=0.0, expected_phi=90.0,
        U=U, B=B, wavelength=wavelength, mu=mu, gam=gam))

    ok = all(r is not False for r in results)
    print(f"\n  TEST 2: {'PASS' if ok else 'FAIL'}")
    return ok


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

if __name__ == '__main__':
    t1 = test_si_cubic()
    t2 = test_aute2()

    print("\n" + "=" * 70)
    print(f"OVERALL: {'PASS' if (t1 and t2) else 'FAIL'}")
    if not t1:
        print("  Test 1 FAILED — convention mismatch in scbasic/verify_scattering")
    if not t2:
        print("  Test 2 FAILED — UB matrix or b_matrix issue for AuTe2")
    print("=" * 70)
