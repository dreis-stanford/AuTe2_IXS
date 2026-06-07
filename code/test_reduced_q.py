"""
Test if we should use reduced q (first BZ) instead of full Q for IXS
"""

import numpy as np
import sys
sys.path.insert(0, 'code')

from .force_constants import ForceConstants
from .phonons import calc_Dq, calc_freq_eig, convert_frequencies
from .form_factors import CalcAtomicfQ
from .ixs import calc_ixs

# Load structure
fc_file = "data/AuTe_2_m.fc"
xtal = ForceConstants(fc_file)
Phi = xtal.convert_to_eV_per_Angstrom2()

masses = np.array([xtal.masses[xtal.atom_type_map[i]-1] 
                   for i in range(xtal.nat)])

# Full Q
Q_full = np.array([[0.0, 0.0, 4.1]])
kT_THz = 6.21

print("="*80)
print("Comparing IXS with full Q vs reduced q")
print("="*80)

# Find reduced q
Q_full_rounded = np.round(Q_full)
Q_reduced = Q_full - Q_full_rounded

print(f"\nFull Q:     {Q_full[0]}")
print(f"Nearest G:  {Q_full_rounded[0]}")
print(f"Reduced q:  {Q_reduced[0]}")

# Phonons are calculated at REDUCED q
print("\n" + "="*80)
print("Calculate phonons at REDUCED q")
print("="*80)

Dq = calc_Dq(Q_reduced, xtal.uvw, Phi, masses)
w_raw, ev = calc_freq_eig(Dq, sort_modes=False, rotate_ev=True)
w_cm = convert_frequencies(w_raw, 'cm-1').flatten()
w_THz = convert_frequencies(w_raw, 'THz')

print(f"\nMode 3 frequency: {w_cm[2]:.2f} cm⁻¹")

# Form factors at FULL |Q|
print("\n" + "="*80)
print("Option 1: Use FULL Q for form factors and structure factor")
print("="*80)

Q_cart_full = Q_full @ xtal.b_l
Q_mag_full = np.linalg.norm(Q_cart_full)

print(f"|Q_full| = {Q_mag_full:.4f} (2π/Å)")

fAu_full = CalcAtomicfQ(Q_mag_full, 'Au', scale=4*np.pi, use_xraylib=False)
fTe_full = CalcAtomicfQ(Q_mag_full, 'Te', scale=4*np.pi, use_xraylib=False)
fQ_full = np.array([[fAu_full, fTe_full, fTe_full]])

print(f"f(Au) = {fAu_full:.2f}")
print(f"f(Te) = {fTe_full:.2f}")

Is_full, _, _, F_full, _ = calc_ixs(
    w_THz, ev, Q_full, xtal.b_l, xtal.xs, fQ_full, masses, kT_THz,
    units='barn', per_steradian=True
)

print(f"\nMode 3:")
print(f"  |F|² = {np.abs(F_full[0,2])**2:.6f}")
print(f"  IXS  = {Is_full[0,2]:.3f} barn")

# Form factors at REDUCED |q|
print("\n" + "="*80)
print("Option 2: Use REDUCED q for form factors and structure factor")
print("="*80)

Q_cart_reduced = Q_reduced @ xtal.b_l
Q_mag_reduced = np.linalg.norm(Q_cart_reduced)

print(f"|q_reduced| = {Q_mag_reduced:.4f} (2π/Å)")

fAu_red = CalcAtomicfQ(Q_mag_reduced, 'Au', scale=4*np.pi, use_xraylib=False)
fTe_red = CalcAtomicfQ(Q_mag_reduced, 'Te', scale=4*np.pi, use_xraylib=False)
fQ_red = np.array([[fAu_red, fTe_red, fTe_red]])

print(f"f(Au) = {fAu_red:.2f}")
print(f"f(Te) = {fTe_red:.2f}")

Is_red, _, _, F_red, _ = calc_ixs(
    w_THz, ev, Q_reduced, xtal.b_l, xtal.xs, fQ_red, masses, kT_THz,
    units='barn', per_steradian=True
)

print(f"\nMode 3:")
print(f"  |F|² = {np.abs(F_red[0,2])**2:.6f}")
print(f"  IXS  = {Is_red[0,2]:.3f} barn")

print("\n" + "="*80)
print(f"MATLAB reports for Q=(0,0,4.1), Mode 3:")
print(f"  IXS = 655.4 barn")
print(f"  Uses |Q_full| = {Q_mag_full:.4f} for form factors")
print(f"  Uses q_reduced for phonons and phase factors")
print("="*80)

