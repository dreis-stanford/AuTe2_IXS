"""
Debug the IXS calculation step-by-step
Compare with MATLAB at Q=(0,0,0.1), Mode 3
"""

import numpy as np
import sys
sys.path.insert(0, 'code')

from .constants import const
from .force_constants import ForceConstants
from .phonons import calc_Dq, calc_freq_eig, convert_frequencies
from .form_factors import CalcAtomicfQ

# Load structure
fc_file = "data/AuTe_2_m.fc"
xtal = ForceConstants(fc_file)
Phi = xtal.convert_to_eV_per_Angstrom2()

masses = np.array([xtal.masses[xtal.atom_type_map[i]-1] 
                   for i in range(xtal.nat)])

Q_prim = np.array([[0.0, 0.0, 0.1]])

# Temperature
kT_cm = 207
kT_THz = kT_cm * const.c * 100 / 1e12

print("="*80)
print("Step-by-step IXS calculation for Mode 3")
print("="*80)

# Calculate modes
Dq = calc_Dq(Q_prim, xtal.uvw, Phi, masses)
w_raw, ev = calc_freq_eig(Dq, sort_modes=False, rotate_ev=False)
w_THz = convert_frequencies(w_raw, 'THz').flatten()
w_cm = convert_frequencies(w_raw, 'cm-1').flatten()

mode_idx = 2
print(f"\nMode {mode_idx+1}: {w_cm[mode_idx]:.2f} cm⁻¹ = {w_THz[mode_idx]:.3f} THz")

# Q vector
Q_cart = Q_prim @ xtal.b_l
Q_mag = np.linalg.norm(Q_cart)

print(f"\nQ_prim = {Q_prim[0]}")
print(f"Q_cart = {Q_cart[0]} (2π/Å)")
print(f"|Q| = {Q_mag:.4f} (2π/Å)")

# Form factors
fAu = CalcAtomicfQ(Q_mag, 'Au', scale=4*np.pi, use_xraylib=False)
fTe = CalcAtomicfQ(Q_mag, 'Te', scale=4*np.pi, use_xraylib=False)

print(f"\nForm factors:")
print(f"  f(Au) = {fAu:.2f}")
print(f"  f(Te) = {fTe:.2f}")

# Structure factor calculation (manual)
print("\n" + "="*80)
print("Structure Factor Calculation")
print("="*80)

ev_mode = ev[:, mode_idx, 0].reshape(3, 3)

# Geometric phases
phases = np.exp(-2j * np.pi * (Q_prim @ xtal.xs.T))
print(f"\nPhase factors:")
for i in range(3):
    print(f"  Atom {i}: exp(-2πi Q·r) = {phases[0,i]}")

# Calculate F for this mode
F = 0.0 + 0.0j

print(f"\nStructure factor contributions:")

for iat in range(3):
    ev_atom = ev_mode[:, iat]
    ev_phys = ev_atom / np.sqrt(masses[iat])
    
    # Q · e
    Q_dot_e = np.sum(Q_cart[0] * ev_phys)
    
    # Form factor
    f_atom = fAu if iat == 0 else fTe
    
    # Contribution
    contrib = f_atom * phases[0, iat] * Q_dot_e / np.sqrt(masses[iat])
    F += contrib
    
    symbol = xtal.symbols[xtal.atom_type_map[iat]-1]
    print(f"  {symbol}: f={f_atom:.2f}, phase={phases[0,iat]:.3f}, "
          f"Q·e={Q_dot_e:.6f}, contrib={contrib:.6f}")

print(f"\nTotal F = {F:.6f}")
print(f"|F|² = {np.abs(F)**2:.6f}")

# Bose factor
n = 1.0 / (np.exp(w_THz[mode_idx] / kT_THz) - 1)
print(f"\nBose factor n(ω) = {n:.4f}")

# IXS formula
common_factor = np.abs(F)**2 / w_THz[mode_idx]
Is_manual = common_factor * (n + 1)
Ias_manual = common_factor * n

print(f"\nCommon factor |F|²/ω = {common_factor:.6f}")

# Normalization
r_e_cm = 2.8179403262e-13  # cm
barn_to_cm2 = 1e-24
norm_barn = (r_e_cm**2) / barn_to_cm2

print(f"\nNormalization to barns:")
print(f"  r_e² = {r_e_cm**2:.6e} cm²")
print(f"  r_e² in barns = {norm_barn:.6f}")

Is_barn = Is_manual * norm_barn
Ias_barn = Ias_manual * norm_barn

print(f"\nFinal IXS (Python):")
print(f"  IXS(Stokes)      = {Is_barn:.3f} barn")
print(f"  IXS(Anti-Stokes) = {Ias_barn:.3f} barn")

print(f"\nMATLAB reports:")
print(f"  IXS(Stokes)      = 655.4 barn")
print(f"  IXS(Anti-Stokes) = 593.1 barn")

ratio = 655.4 / Is_barn if Is_barn > 0 else 0
print(f"\nRatio MATLAB/Python = {ratio:.1f}")

# Check if there's a factor of frequency difference
print(f"\nω (THz) = {w_THz[mode_idx]:.4f}")
print(f"ω (cm⁻¹) = {w_cm[mode_idx]:.4f}")
print(f"Ratio × ω(THz) = {ratio * w_THz[mode_idx]:.1f}")
print(f"Ratio × ω(cm⁻¹) = {ratio * w_cm[mode_idx]:.1f}")

