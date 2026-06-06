#!/usr/bin/env python
"""
Interactive IXS analysis for Silicon
Supports primitive and conventional coordinates
"""

import sys
sys.path.insert(0, 'code')

from force_constants import ForceConstants
from phonons import calc_Dq, calc_freq_eig, convert_frequencies
from form_factors import CalcAtomicfQ
from ixs import calc_ixs
from fcc_structure import fcc_conv2prim_k
import numpy as np

def analyze_si_q(xtal, Phi, masses, kT_THz, Q_prim):
    """Analyze Si at single Q-point"""
    
    Q_prim = np.asarray(Q_prim).reshape(1, 3)
    Q_rounded = np.round(Q_prim)
    Q_reduced = Q_prim - Q_rounded
    
    Dq = calc_Dq(Q_reduced, xtal.uvw, Phi, masses)
    w_raw, ev = calc_freq_eig(Dq, sort_modes=True, rotate_ev=True)
    
    w_THz = convert_frequencies(w_raw, 'THz')
    w_cm = convert_frequencies(w_raw, 'cm-1').flatten()
    w_meV = convert_frequencies(w_raw, 'meV').flatten()
    
    Q_cart_full = Q_prim @ xtal.b_l
    Q_mag_full = np.linalg.norm(Q_cart_full)
    
    fSi = CalcAtomicfQ(Q_mag_full, 'Si', scale=4*np.pi, use_xraylib=False)
    fQ_matrix = np.array([[fSi, fSi]])
    
    Is, Ias, n, F, info = calc_ixs(
        w_THz, ev, Q_prim, xtal.b_l, xtal.xs,
        fQ_matrix, masses, kT_THz, units='barn', per_steradian=True
    )
    
    # Print
    print("\n╔" + "="*64 + "╗")
    print("║" + "  IXS Analysis - Silicon".center(64) + "║")
    print("╚" + "="*64 + "╝")
    
    print(f"\nQ (primitive):    [{Q_prim[0,0]:.4f}, {Q_prim[0,1]:.4f}, {Q_prim[0,2]:.4f}]")
    print(f"|Q| = {Q_mag_full:.4f} (2π/Å)")
    print(f"Reduced q:        [{Q_reduced[0,0]:.4f}, {Q_reduced[0,1]:.4f}, {Q_reduced[0,2]:.4f}]")
    if not np.allclose(Q_prim, Q_reduced):
        print(f"G vector:         [{Q_rounded[0,0]:.0f}, {Q_rounded[0,1]:.0f}, {Q_rounded[0,2]:.0f}]")
    
    print(f"\nT = {kT_THz:.2f} THz (~{kT_THz*47.99:.0f} K)")
    print(f"f(Si) = {fSi:.2f}")
    
    print("\n" + "="*80)
    print(f"{'Mode':>4} {'ω(cm⁻¹)':>10} {'ω(meV)':>10} {'IXS(S)':>12} {'IXS(AS)':>12}")
    print("-"*80)
    
    for i in range(len(w_cm)):
        print(f"{i+1:4d} {w_cm[i]:10.2f} {w_meV[i]:10.3f} {Is[0,i]:12.3e} {Ias[0,i]:12.3e}")
    
    print("="*80)

def main():
    fc_file = "data/Test__Silicon_dispersion/Qgrid_444/Cg.fc"
    
    print("="*70)
    print("Silicon IXS Interactive Analysis")
    print("="*70)
    
    choice = input("\nSupercell (1=4x4x4, 2=8x8x8) [1]: ").strip()
    if choice == '2':
        fc_file = "data/Test__Silicon_dispersion/Qgrid_888/Cg.fc"
    
    xtal = ForceConstants(fc_file)
    Phi = xtal.convert_to_eV_per_Angstrom2()
    masses = np.array([xtal.masses[0], xtal.masses[0]])
    kT_THz = 6.2
    
    print("\n" + "="*70)
    print("High-symmetry points (conventional): Γ=(0,0,0) X=(0,0.5,0.5)")
    print("                                      K=(3/8,3/8,3/4) L=(0.5,0.5,0.5)")
    print("="*70)
    
    coord_choice = input("\nCoordinates (1=primitive, 2=conventional) [1]: ").strip()
    use_conv = (coord_choice == '2')
    
    print(f"\nUsing {'conventional' if use_conv else 'primitive'} coordinates")
    print("Type 'switch' to change, Enter to quit\n")
    
    while True:
        try:
            coord_name = "conv" if use_conv else "prim"
            user_input = input(f"Q ({coord_name}): ").strip()
            
            if not user_input:
                break
            if user_input.lower() == 'switch':
                use_conv = not use_conv
                print(f"Switched to {'conventional' if use_conv else 'primitive'}")
                continue
            
            Q_parts = user_input.split()
            if len(Q_parts) != 3:
                print("Error: Need 3 numbers")
                continue
            
            Q = np.array([float(x) for x in Q_parts])
            
            # Convert if needed
            if use_conv:
                Q_prim = fcc_conv2prim_k(Q)
                print(f"Converted to primitive: {Q_prim}")
            else:
                Q_prim = Q
            
            analyze_si_q(xtal, Phi, masses, kT_THz, Q_prim)
            
        except KeyboardInterrupt:
            break
        except Exception as e:
            print(f"Error: {e}")

if __name__ == "__main__":
    main()
