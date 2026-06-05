#!/usr/bin/env python
"""
Example 3: Low-level API usage

This example shows how to use the individual modules directly
for custom analysis workflows.
"""

import numpy as np
import sys
import os

sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', 'code'))

from force_constants import ForceConstants
from phonons import calc_Dq, calc_freq_eig, convert_frequencies
from form_factors import CalcAtomicfQ
from ixs import calc_ixs
from constants import const

def main():
    print("="*70)
    print("Example 3: Low-level API Usage")
    print("="*70)
    
    # 1. Load force constants
    print("\n1. Loading force constants...")
    xtal = ForceConstants("../data/AuTe_2_m.fc")
    Phi = xtal.convert_to_eV_per_Angstrom2()
    
    masses = np.array([xtal.masses[xtal.atom_type_map[i]-1] 
                       for i in range(xtal.nat)])
    
    print(f"   Crystal: {xtal.nat} atoms in primitive cell")
    print(f"   Masses: {masses} amu")
    
    # 2. Calculate phonons at a Q-point
    print("\n2. Calculating phonons...")
    
    Q_full = np.array([[0.0, 0.0, 4.1]])  # Full Q (outside first BZ)
    Q_reduced = Q_full - np.round(Q_full)  # Reduce to first BZ
    
    print(f"   Q (full):    {Q_full[0]}")
    print(f"   Q (reduced): {Q_reduced[0]}")
    
    # Calculate dynamical matrix at reduced q
    Dq = calc_Dq(Q_reduced, xtal.uvw, Phi, masses)
    
    print(f"   Dynamical matrix shape: {Dq.shape}")
    
    # Solve for frequencies and eigenvectors
    w_raw, ev = calc_freq_eig(Dq, sort_modes=True, rotate_ev=True)
    
    # Convert to different units
    w_THz = convert_frequencies(w_raw, 'THz')
    w_cm = convert_frequencies(w_raw, 'cm-1')
    w_meV = convert_frequencies(w_raw, 'meV')
    
    print(f"   Number of modes: {w_cm.shape[1]}")
    print(f"   Frequency range: {w_cm[0, 0]:.2f} - {w_cm[0, -1]:.2f} cm⁻¹")
    
    # 3. Calculate form factors
    print("\n3. Calculating atomic form factors...")
    
    Q_cart_full = Q_full @ xtal.b_l
    Q_mag_full = np.linalg.norm(Q_cart_full)
    
    print(f"   |Q| = {Q_mag_full:.4f} (2π/Å)")
    print(f"   sin(θ)/λ = {Q_mag_full/(4*np.pi):.4f} Å⁻¹")
    
    # Calculate form factors using Cromer-Mann
    fAu_CM = CalcAtomicfQ(Q_mag_full, 'Au', scale=4*np.pi, use_xraylib=False)
    fTe_CM = CalcAtomicfQ(Q_mag_full, 'Te', scale=4*np.pi, use_xraylib=False)
    
    print(f"   f(Au) = {fAu_CM:.2f} (Cromer-Mann)")
    print(f"   f(Te) = {fTe_CM:.2f} (Cromer-Mann)")
    
    # Try with xraylib if available
    try:
        fAu_XR = CalcAtomicfQ(Q_mag_full, 'Au', scale=4*np.pi, use_xraylib=True)
        fTe_XR = CalcAtomicfQ(Q_mag_full, 'Te', scale=4*np.pi, use_xraylib=True)
        print(f"   f(Au) = {fAu_XR:.2f} (xraylib)")
        print(f"   f(Te) = {fTe_XR:.2f} (xraylib)")
    except:
        print("   (xraylib not available)")
    
    fQ_matrix = np.array([[fAu_CM, fTe_CM, fTe_CM]])
    
    # 4. Calculate IXS cross-sections
    print("\n4. Calculating IXS cross-sections...")
    
    kT_THz = 6.2  # Temperature
    
    Is, Ias, n, F, xs_info = calc_ixs(
        w_THz, ev, Q_full,  # Use full Q for phase factors
        xtal.b_l, xtal.xs, fQ_matrix, masses, kT_THz,
        units='barn', per_steradian=True
    )
    
    print(f"   Units: {xs_info['units']}")
    print(f"   Normalization: {xs_info['normalization_factor']:.3e}")
    
    # 5. Analyze results
    print("\n5. Results for each mode:")
    print("   " + "-"*66)
    print(f"   {'Mode':>4} {'Freq(cm⁻¹)':>12} {'Freq(meV)':>10} {'IXS(S)':>12} {'IXS(AS)':>12} {'n(ω)':>8}")
    print("   " + "-"*66)
    
    for i in range(w_cm.shape[1]):
        print(f"   {i+1:4d} {w_cm[0,i]:12.2f} {w_meV[0,i]:10.2f} "
              f"{Is[0,i]:12.3f} {Ias[0,i]:12.3f} {n[0,i]:8.3f}")
    
    print("   " + "-"*66)
    
    # 6. Calculate total scattering
    print("\n6. Integrated scattering:")
    
    total_stokes = np.sum(Is)
    total_antistokes = np.sum(Ias)
    total_ixs = total_stokes + total_antistokes
    
    print(f"   Total Stokes:      {total_stokes:.2f} barn")
    print(f"   Total Anti-Stokes: {total_antistokes:.2f} barn")
    print(f"   Total IXS:         {total_ixs:.2f} barn")
    
    # 7. Structure factor analysis
    print("\n7. Structure factor analysis:")
    
    print(f"   {'Mode':>4} {'|F|²':>12} {'ω(THz)':>10} {'|F|²/ω':>12}")
    print("   " + "-"*42)
    
    for i in range(min(5, w_THz.shape[1])):  # Show first 5 modes
        F_sq = np.abs(F[0, i])**2
        F_sq_over_w = F_sq / w_THz[0, i]
        print(f"   {i+1:4d} {F_sq:12.6f} {w_THz[0,i]:10.4f} {F_sq_over_w:12.6f}")
    
    print("   " + "...")
    
    # 8. Temperature dependence
    print("\n8. Temperature dependence (mode 3):")
    
    mode_idx = 2
    temperatures_K = [100, 200, 300, 400, 500]
    
    print(f"   Mode {mode_idx+1}: {w_cm[0, mode_idx]:.2f} cm⁻¹")
    print(f"   {'T (K)':>8} {'kT (THz)':>10} {'n(ω)':>10} {'IXS(S)':>12} {'IXS(AS)':>12}")
    print("   " + "-"*54)
    
    for T_K in temperatures_K:
        kT = T_K * const.kB / const.THz2meV * 1000
        
        # Recalculate IXS at this temperature
        Is_T, Ias_T, n_T, _, _ = calc_ixs(
            w_THz, ev, Q_full,
            xtal.b_l, xtal.xs, fQ_matrix, masses, kT,
            units='barn', per_steradian=True
        )
        
        print(f"   {T_K:8.0f} {kT:10.2f} {n_T[0, mode_idx]:10.3f} "
              f"{Is_T[0, mode_idx]:12.3f} {Ias_T[0, mode_idx]:12.3f}")
    
    print("\n" + "="*70)
    print("Example complete!")
    print("="*70)
    print("\nKey observations:")
    print("  - Phonons calculated at reduced q (first BZ)")
    print("  - Form factors use full |Q| (including G vector)")
    print("  - IXS intensity ∝ |F|²/ω × (n+1) for Stokes")
    print("  - Anti-Stokes suppressed at low T (small n)")

if __name__ == "__main__":
    main()
