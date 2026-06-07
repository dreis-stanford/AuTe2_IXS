"""
Test the difference between rotate_ev=True and False
"""

import numpy as np
import sys
sys.path.insert(0, 'code')

from .force_constants import ForceConstants
from .phonons import calc_Dq, calc_freq_eig, convert_frequencies
from form_factors import CalcAtomicfQ
from ixs import calc_ixs

# Load structure
fc_file = "data/AuTe_2_m.fc"
xtal = ForceConstants(fc_file)
Phi = xtal.convert_to_eV_per_Angstrom2()

masses = np.array([xtal.masses[xtal.atom_type_map[i]-1] 
                   for i in range(xtal.nat)])

Q_prim = np.array([[0.0, 0.0, 0.1]])
kT_THz = 6.21

Dq = calc_Dq(Q_prim, xtal.uvw, Phi, masses)

print("="*80)
print("Comparing rotate_ev=False vs rotate_ev=True")
print("="*80)

# Without rotation
print("\n" + "="*80)
print("WITHOUT rotate_ev (rotate_ev=False)")
print("="*80)

w_raw, ev_no_rot = calc_freq_eig(Dq, sort_modes=False, rotate_ev=False)
w_cm = convert_frequencies(w_raw, 'cm-1').flatten()

# Calculate IXS
Q_cart = Q_prim @ xtal.b_l
Q_mag = np.linalg.norm(Q_cart)
fAu = CalcAtomicfQ(Q_mag, 'Au', scale=4*np.pi, use_xraylib=False)
fTe = CalcAtomicfQ(Q_mag, 'Te', scale=4*np.pi, use_xraylib=False)
fQ = np.array([[fAu, fTe, fTe]])

w_THz = convert_frequencies(w_raw, 'THz')
Is_no_rot, _, _, F_no_rot, _ = calc_ixs(
    w_THz, ev_no_rot, Q_prim, xtal.b_l, xtal.xs, fQ, masses, kT_THz,
    units='barn', per_steradian=True
)

print(f"\nMode 3: {w_cm[2]:.2f} cm⁻¹")
print(f"  |F|² = {np.abs(F_no_rot[0,2])**2:.6f}")
print(f"  IXS = {Is_no_rot[0,2]:.3f} barn")

# Show eigenvector
ev_mode = ev_no_rot[:, 2, 0].reshape(3, 3)
print(f"\n  Eigenvector amplitudes:")
for i in range(3):
    ev_atom = ev_mode[:, i] / np.sqrt(masses[i])
    amp = np.linalg.norm(ev_atom)
    symbol = xtal.symbols[xtal.atom_type_map[i]-1]
    print(f"    {symbol}: |u| = {amp:.6f}")

# WITH rotation
print("\n" + "="*80)
print("WITH rotate_ev (rotate_ev=True)")
print("="*80)

w_raw2, ev_rot = calc_freq_eig(Dq, sort_modes=False, rotate_ev=True)

Is_rot, _, _, F_rot, _ = calc_ixs(
    w_THz, ev_rot, Q_prim, xtal.b_l, xtal.xs, fQ, masses, kT_THz,
    units='barn', per_steradian=True
)

print(f"\nMode 3: {w_cm[2]:.2f} cm⁻¹")
print(f"  |F|² = {np.abs(F_rot[0,2])**2:.6f}")
print(f"  IXS = {Is_rot[0,2]:.3f} barn")

# Show eigenvector
ev_mode = ev_rot[:, 2, 0].reshape(3, 3)
print(f"\n  Eigenvector amplitudes:")
for i in range(3):
    ev_atom = ev_mode[:, i] / np.sqrt(masses[i])
    amp = np.linalg.norm(ev_atom)
    symbol = xtal.symbols[xtal.atom_type_map[i]-1]
    print(f"    {symbol}: |u| = {amp:.6f}")

print("\n" + "="*80)
print(f"MATLAB reports for Mode 3:")
print(f"  IXS = 655.4 barn")
print("="*80)

print(f"\nPython with rotate_ev=False: {Is_no_rot[0,2]:.3f} barn")
print(f"Python with rotate_ev=True:  {Is_rot[0,2]:.3f} barn")

print("\nDifference in |F|²:")
print(f"  Without rotation: {np.abs(F_no_rot[0,2])**2:.6f}")
print(f"  With rotation:    {np.abs(F_rot[0,2])**2:.6f}")
print(f"  Ratio:            {np.abs(F_rot[0,2])**2 / np.abs(F_no_rot[0,2])**2:.1f}")

