"""
Verify six-circle diffractometer scattering geometry

Lab frame (from sixcircle convention):
- x: horizontal (tth/detector rotates about z, but x is perpendicular to beam)
- y: beam direction (source at +y, beam travels along -y to sample)
- z: vertical (theta rotates about z)

Chi: Right-handed rotation about -y (incident beam direction when th=0, mu=0)

Initial orientation (all angles = 0):
- c-axis along +z (vertical)
- b-axis along -y (beam direction)
- a-axis along +x
"""

import numpy as np


def rotate_Q_to_lab_frame(Q_crystal, chi, phi, theta, mu):
    """
    Rotate Q from crystal frame to lab frame
    
    Initial crystal orientation (all angles = 0):
    - a-axis along +x
    - b-axis along -y (beam)
    - c-axis along +z (vertical/normal)
    
    Rotation sequence: R_mu @ R_theta @ R_chi @ R_phi
    
    Parameters:
    -----------
    Q_crystal : array (3,)
        Q in crystal Cartesian coordinates (a, b, c basis)
    chi, phi, theta, mu : float
        Angles in degrees
        
    Returns:
    --------
    Q_lab : array (3,)
        Q in lab frame (x, y, z)
    """
    
    # Convert to radians
    chi_r = np.radians(chi)
    phi_r = np.radians(phi)
    theta_r = np.radians(theta)
    mu_r = np.radians(mu)
    
    # phi: rotation about z-axis (when chi=0, phi and theta are degenerate)
    R_phi = np.array([
        [np.cos(phi_r), -np.sin(phi_r), 0],
        [np.sin(phi_r),  np.cos(phi_r), 0],
        [0, 0, 1]
    ])
    
    # chi: RIGHT-HANDED rotation about -y axis (beam direction)
    # This is equivalent to LEFT-HANDED rotation about +y
    # So we negate the angle
    R_chi = np.array([
        [np.cos(-chi_r), 0, np.sin(-chi_r)],
        [0, 1, 0],
        [-np.sin(-chi_r), 0, np.cos(-chi_r)]
    ])
    
    # theta: rotation about z-axis
    R_theta = np.array([
        [np.cos(theta_r), -np.sin(theta_r), 0],
        [np.sin(theta_r),  np.cos(theta_r), 0],
        [0, 0, 1]
    ])
    
    # mu: rotation about x-axis (beam tilt correction)
    R_mu = np.array([
        [1, 0, 0],
        [0, np.cos(mu_r), -np.sin(mu_r)],
        [0, np.sin(mu_r),  np.cos(mu_r)]
    ])
    
    # Combined rotation: mu, then theta, then chi, then phi
    R_total = R_mu @ R_theta @ R_chi @ R_phi
    
    Q_lab = R_total @ Q_crystal
    
    return Q_lab


def verify_basic(Q_hkl, angles, wavelength, b_matrix):
    """Basic verification: |Q| = 2|k|sin(θ)"""
    
    Q_hkl_arr = np.array(Q_hkl)
    Q_crystal = Q_hkl_arr @ b_matrix
    Q_mag = np.linalg.norm(Q_crystal)
    
    k_mag = 2 * np.pi / wavelength
    theta = np.radians(abs(angles['tth']) / 2)
    Q_mag_from_angle = 2 * k_mag * np.sin(theta)
    
    error = abs(Q_mag - Q_mag_from_angle)
    error_rel = error / Q_mag if Q_mag > 0 else 0
    
    return {
        'Q_hkl': Q_hkl,
        'Q_mag': Q_mag,
        'Q_mag_from_angle': Q_mag_from_angle,
        'error': error,
        'error_rel': error_rel,
        'k_mag': k_mag,
        'theta_deg': np.degrees(theta),
    }


def verify_full(Q_hkl, angles, wavelength, b_matrix):
    """
    Full verification: Rotate Q to lab frame and check k_out = k_in + Q_lab
    
    Beam geometry:
    - k_in along -y: [0, -1, 0]
    - k_out with gamma=0: [sin(tth), -cos(tth), 0]
    """
    
    # Q from lattice
    Q_hkl_arr = np.array(Q_hkl)
    Q_crystal = Q_hkl_arr @ b_matrix
    
    # Rotate Q to lab frame
    Q_lab = rotate_Q_to_lab_frame(
        Q_crystal,
        angles['chi'],
        angles['phi'],
        angles['th'],  # theta (not tth)
        angles['mu']
    )
    
    # Wave vectors
    k_mag = 2 * np.pi / wavelength
    mu_r = np.radians(angles.get('mu', 0.0))
    k_in = k_mag * np.array([0, -np.cos(mu_r), -np.sin(mu_r)])

    # Expected k_out from detector angles (with gamma=0)
    tth_r = np.radians(angles['tth'])
    gam_r = np.radians(angles.get('gam', 0.0))
    k_out_expected = k_mag * np.array([
        np.sin(tth_r) * np.cos(gam_r),
        -np.cos(tth_r) * np.cos(gam_r),
        -np.sin(gam_r)
    ])
    
    # k_out from momentum conservation
    k_out_calc = k_in + Q_lab
    k_out_mag = np.linalg.norm(k_out_calc)
    
    # Errors
    k_mag_error = abs(k_out_mag - k_mag)
    k_mag_error_rel = k_mag_error / k_mag
    
    k_dir_error = np.linalg.norm(k_out_calc - k_out_expected)
    
    return {
        'Q_hkl': Q_hkl,
        'Q_crystal': Q_crystal,
        'Q_lab': Q_lab,
        'k_in': k_in,
        'k_out_calc': k_out_calc,
        'k_out_expected': k_out_expected,
        'k_mag': k_mag,
        'k_out_mag': k_out_mag,
        'k_mag_error': k_mag_error,
        'k_mag_error_rel': k_mag_error_rel,
        'k_dir_error': k_dir_error,
        'angles': angles,
    }


def print_verification(basic_results, full_results=None):
    """Print verification results"""
    
    print("\n" + "="*80)
    print("Diffractometer Geometry Verification")
    print("="*80)
    
    # Basic test
    print(f"\nBASIC TEST: |Q| = 2|k|sin(θ)")
    print(f"  Q(hkl) = {basic_results['Q_hkl']}")
    print(f"  |Q| from lattice = {basic_results['Q_mag']:.6f} (2π/Å)")
    print(f"  |Q| from angle   = {basic_results['Q_mag_from_angle']:.6f} (2π/Å)")
    print(f"  Error = {basic_results['error']:.2e} (relative: {basic_results['error_rel']:.2e})")
    
    if basic_results['error_rel'] < 1e-4:
        print("  ✓ PASS")
        basic_status = "PASS"
    elif basic_results['error_rel'] < 1e-2:
        print("  ⚠ WARNING")
        basic_status = "WARNING"
    else:
        print("  ✗ FAIL")
        basic_status = "FAIL"
    
    # Full test
    if full_results:
        print(f"\nFULL TEST: k_out = k_in + Q_lab")
        print(f"  Q_lab = [{full_results['Q_lab'][0]:7.4f}, {full_results['Q_lab'][1]:7.4f}, {full_results['Q_lab'][2]:7.4f}]")
        print(f"  k_out (calculated) = [{full_results['k_out_calc'][0]:7.4f}, {full_results['k_out_calc'][1]:7.4f}, {full_results['k_out_calc'][2]:7.4f}]")
        print(f"  k_out (expected)   = [{full_results['k_out_expected'][0]:7.4f}, {full_results['k_out_expected'][1]:7.4f}, {full_results['k_out_expected'][2]:7.4f}]")
        print(f"  |k_out| = {full_results['k_out_mag']:.6f} (should be {full_results['k_mag']:.6f})")
        print(f"  Elastic error = {full_results['k_mag_error_rel']:.2e}")
        print(f"  Direction error = {full_results['k_dir_error']:.4f}")
        
        if full_results['k_mag_error_rel'] < 1e-5 and full_results['k_dir_error'] < 0.01:
            print("  ✓ PASS")
            full_status = "PASS"
        elif full_results['k_mag_error_rel'] < 1e-3 and full_results['k_dir_error'] < 0.1:
            print("  ⚠ WARNING")
            full_status = "WARNING"
        else:
            print("  ✗ FAIL")
            full_status = "FAIL"
    else:
        full_status = None
    
    print("="*80 + "\n")
    
    return basic_status, full_status
