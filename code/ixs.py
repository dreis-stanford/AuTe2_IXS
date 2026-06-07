"""
Inelastic X-ray Scattering (IXS) cross-section calculations
Converted from MATLAB CalcIXS.m
"""

import numpy as np
from .constants import const

def calc_ixs(w, ev, Q_rlu, b_l, xs, fQ, masses, kT, 
             units='arbitrary', per_steradian=True):
    """
    Calculate IXS cross-sections for phonons
    
    Matches MATLAB CalcIXS.m exactly
    
    Parameters:
    -----------
    w : ndarray (nQ, ndof)
        Phonon frequencies in THz
    ev : ndarray (ndof, ndof, nQ)
        Eigenvectors (3*ns x ndof x nQ)
        Mass-weighted eigenvectors from calc_freq_eig
    Q_rlu : ndarray (nQ, 3)
        Scattering vectors in reciprocal lattice units
    b_l : ndarray (3, 3)
        Reciprocal lattice vectors (rows: b1, b2, b3) in Cartesian (2π/Å)
    xs : ndarray (ns, 3)
        Atomic positions (fractional coordinates)
    fQ : ndarray (nQ, ns)
        Atomic form factors for each Q and atom
    masses : ndarray (ns,)
        Atomic masses in amu
    kT : float
        Temperature in THz (same units as w)
    units : str
        'arbitrary', 're_squared', or 'barn'
    per_steradian : bool
        If True, cross-section per steradian
    
    Returns:
    --------
    Is : ndarray (nQ, ndof)
        Stokes scattering intensity (phonon creation)
    Ias : ndarray (nQ, ndof)
        Anti-Stokes intensity (phonon annihilation)
    n : ndarray (nQ, ndof)
        Bose-Einstein occupation factor
    F : ndarray (nQ, ndof)
        Structure factor (complex)
    cross_section_info : dict
        Information about units and normalization
    """
    
    nQ, ndof = w.shape
    ns = xs.shape[0]
    
    # Ensure masses is 1D
    masses = np.asarray(masses).flatten()
    
    # 1. Frequency thresholding to prevent division by zero at Gamma
    min_w = 0.1  # THz
    w_safe = w.copy()
    w_safe[w_safe < min_w] = min_w
    
    # 2. Thermal Occupation (Bose-Einstein factor)
    # n(ω) = 1 / (exp(ω/kT) - 1)
    n = 1.0 / (np.exp(w_safe / kT) - 1)
    
    # 3. Geometric Phase Factors [nQ x ns]
    # Phase: exp(-i * 2π * Q_rlu · r)
    # Note: Using -2πi convention for scattering
    phases = np.exp(-2j * np.pi * (Q_rlu @ xs.T))  # [nQ x ns]
    
    # 4. Convert Q to Cartesian coordinates
    Q_cart = Q_rlu @ b_l  # [nQ x 3] in units of 2π/Å
    
    # 5. Calculate Dynamic Structure Factor F(Q, ω)
    # Initialize
    F = np.zeros((nQ, ndof), dtype=complex)
    
    # Loop over atoms in the basis
    for s in range(ns):
        # Extract displacement vectors for atom s: [3 x ndof x nQ]
        atom_ev = ev[3*s:3*s+3, :, :]
        
        # Calculate Q · e_s for each mode and q-point
        # Q_cart: [nQ x 3]
        # atom_ev: [3 x ndof x nQ]
        # We need Q · e for each (q, mode) pair
        
        # Reshape Q_cart for broadcasting: [nQ x 3 x 1]
        Q_reshaped = Q_cart[:, :, np.newaxis]
        
        # Dot product along Cartesian dimension: [nQ x ndof]
        # Sum over axis 1 (Cartesian components)
        Qdot_e = np.sum(Q_reshaped * atom_ev.transpose(2, 0, 1), axis=1)
        
        # Add contribution from atom s
        # fQ[:, s]: [nQ] form factors for atom s
        # masses[s]: scalar mass
        # phases[:, s]: [nQ] phase factors
        # Qdot_e: [nQ x ndof]
        
        f_over_sqrtM = fQ[:, s] / np.sqrt(masses[s])  # [nQ]
        phase_s = phases[:, s]  # [nQ]
        
        # Broadcast and accumulate
        F += (f_over_sqrtM * phase_s)[:, np.newaxis] * Qdot_e
    
    # 6. Calculate Intensities
    # Common factor: |F|² / ω
    common_factor = np.abs(F)**2 / w_safe
    
    # Stokes (phonon creation): (n + 1) factor
    Is = common_factor * (n + 1)
    
    # Anti-Stokes (phonon annihilation): n factor
    Ias = common_factor * n
    
    # 7. Apply Physical Normalization
    # Classical electron radius
    r_e_angstrom = 2.8179403262e-5  # Å
    r_e_cm = 2.8179403262e-13       # cm
    barn_to_cm2 = 1e-24             # 1 barn = 10^-24 cm²
    
    if units == 're_squared':
        # Cross-section in r_e² per unit cell
        normalization = r_e_angstrom**2
        
        if per_steradian:
            unit_string = 'r_e²/(unit cell·sr)'
        else:
            unit_string = 'r_e²/unit cell'
    
    elif units == 'barn':
        # Cross-section in barn per unit cell
        # r_e² = 7.94 × 10^-26 cm² = 79.4 millibarn
        normalization = (r_e_cm**2) / barn_to_cm2  # Convert r_e² (cm²) to barns
        
        if per_steradian:
            unit_string = 'barn/(unit cell·sr)'
        else:
            unit_string = 'barn/unit cell'
    
    elif units == 'arbitrary':
        normalization = 1.0
        unit_string = 'arbitrary units'
    
    else:
        raise ValueError(f"Unknown units: {units}. Use 'arbitrary', 're_squared', or 'barn'")
    
    # Apply normalization
    Is = Is * normalization
    Ias = Ias * normalization
    
    # 8. Store cross-section information
    cross_section_info = {
        'units': unit_string,
        'normalization_factor': normalization,
        'r_e_angstrom': r_e_angstrom,
        'r_e_cm': r_e_cm,
        'r_e_squared_in_barns': (r_e_cm**2) / barn_to_cm2,
        'per_steradian': per_steradian
    }
    
    return Is, Ias, n, F, cross_section_info


if __name__ == "__main__":
    import sys
    import os
    sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
    
    from force_constants import ForceConstants
    from phonons import calc_Dq, calc_freq_eig, convert_frequencies
    from form_factors import CalcAtomicfQ
    
    print("=" * 70)
    print("Testing IXS Cross-Section Calculator")
    print("=" * 70)
    
    # Load force constants
    fc_file = "data/AuTe_2_m.fc"
    xtal = ForceConstants(fc_file)
    Phi = xtal.convert_to_eV_per_Angstrom2()
    
    # Get masses for each atom
    masses = np.array([xtal.masses[xtal.atom_type_map[i]-1] 
                       for i in range(xtal.nat)])
    
    print(f"\nMasses: {masses} amu")
    print(f"Symbols: {[xtal.symbols[xtal.atom_type_map[i]-1] for i in range(xtal.nat)]}")
    
    # Test at a few Q points
    print("\n" + "=" * 70)
    print("Calculating IXS at Q = (0.5, 0, 0)")
    print("=" * 70)
    
    Q_test = np.array([[0.5, 0.0, 0.0]])
    
    # Calculate dynamical matrix
    D = calc_Dq(Q_test, xtal.uvw, Phi, masses)
    
    # Solve for modes
    w_raw, ev = calc_freq_eig(D, sort_modes=False, rotate_ev=True)
    
    # Convert to THz
    w_THz = convert_frequencies(w_raw, 'THz')
    w_cm = convert_frequencies(w_raw, 'cm-1')
    
    print(f"\nPhonon frequencies:")
    for i in range(w_THz.shape[1]):
        print(f"  Mode {i+1}: {w_THz[0,i]:8.3f} THz ({w_cm[0,i]:8.2f} cm⁻¹)")
    
    # Calculate Q in Cartesian (for form factors)
    Q_cart = Q_test @ xtal.b_l  # [1 x 3] in 2π/Å
    Q_mag = np.linalg.norm(Q_cart, axis=1)  # [1]
    
    print(f"\nQ magnitude: {Q_mag[0]:.4f} (2π/Å)")
    print(f"sin(θ)/λ = Q/(4π) = {Q_mag[0]/(4*np.pi):.4f} Å⁻¹")
    
    # Calculate form factors
    fAu = CalcAtomicfQ(Q_mag, 'Au', scale=4*np.pi, use_xraylib=False)
    fTe = CalcAtomicfQ(Q_mag, 'Te', scale=4*np.pi, use_xraylib=False)
    
    # Create fQ matrix [nQ x ns]
    fQ = np.column_stack([fAu, fTe, fTe])
    
    print(f"\nForm factors:")
    print(f"  f(Au) = {fAu[0]:.2f}")
    print(f"  f(Te) = {fTe[0]:.2f}")
    
    # Temperature
    kT_cm = 207  # cm⁻¹ (wavenumbers)
    kT_THz = kT_cm * const.c * 100 / 1e12  # Convert to THz
    
    print(f"\nTemperature: {kT_cm:.1f} cm⁻¹ = {kT_THz:.2f} THz")
    
    # Calculate IXS cross-sections
    print("\n" + "=" * 70)
    print("Calculating IXS cross-sections...")
    print("=" * 70)
    
    Is, Ias, n, F, info = calc_ixs(
        w_THz, ev, Q_test, xtal.b_l, xtal.xs, fQ, masses, kT_THz,
        units='barn', per_steradian=True
    )
    
    print(f"\nCross-section units: {info['units']}")
    print(f"Normalization factor: {info['normalization_factor']:.3e}")
    
    print(f"\n{'Mode':>6} {'ω (THz)':>10} {'ω (cm⁻¹)':>10} {'n(ω)':>10} {'IXS(S)':>12} {'IXS(AS)':>12}")
    print("-" * 70)
    
    for i in range(w_THz.shape[1]):
        ixs_s = Is[0, i]
        ixs_as = Ias[0, i]
        
        # Format IXS values
        if ixs_s < 1e-10:
            ixs_s_str = "~0"
        elif ixs_s < 0.01:
            ixs_s_str = f"{ixs_s:.2e}"
        else:
            ixs_s_str = f"{ixs_s:.4f}"
        
        if ixs_as < 1e-10:
            ixs_as_str = "~0"
        elif ixs_as < 0.01:
            ixs_as_str = f"{ixs_as:.2e}"
        else:
            ixs_as_str = f"{ixs_as:.4f}"
        
        print(f"{i+1:6d} {w_THz[0,i]:10.3f} {w_cm[0,i]:10.2f} "
              f"{n[0,i]:10.4f} {ixs_s_str:>12} {ixs_as_str:>12}")
    
    print("\n" + "=" * 70)
    print("IXS calculator ready!")
    print("=" * 70)
    print("\nNext: Create single_q_analysis.py for interactive analysis")
