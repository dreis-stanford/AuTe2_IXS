"""
Debug eigenvector normalization and L-char calculation
Compare with MATLAB
"""

import numpy as np
import sys
sys.path.insert(0, 'code')

from .force_constants import ForceConstants
from .phonons import calc_Dq, calc_freq_eig, convert_frequencies
from aute2_structure import aute2_conv2prim_k, aute2_prim2conv_k

# Load structure
fc_file = "data/AuTe_2_m.fc"
xtal = ForceConstants(fc_file)
Phi = xtal.convert_to_eV_per_Angstrom2()

masses = np.array([xtal.masses[xtal.atom_type_map[i]-1] 
                   for i in range(xtal.nat)])

print("="*70)
print("Eigenvector Normalization Diagnostics")
print("="*70)

# Test Q point
Q_prim = np.array([[0.0, 0.0, 0.1]])
Q_cart = Q_prim @ xtal.b_l
Q_mag = np.linalg.norm(Q_cart)
Q_hat = Q_cart / Q_mag

print(f"\nQ = {Q_prim[0]}")
print(f"Q_cart = {Q_cart[0]}")
print(f"|Q| = {Q_mag:.4f}")
print(f"Q_hat = {Q_hat[0]}")

# Calculate modes
Dq = calc_Dq(Q_prim, xtal.uvw, Phi, masses)
w_raw, ev = calc_freq_eig(Dq, sort_modes=True, rotate_ev=True)
w_cm = convert_frequencies(w_raw, 'cm-1').flatten()

print(f"\nFrequencies: {w_cm[:3]}")

# Focus on mode 3 (the longitudinal one)
mode_idx = 2  # 0-indexed

print(f"\n{'='*70}")
print(f"Analyzing Mode {mode_idx+1}: {w_cm[mode_idx]:.2f} cm⁻¹")
print(f"{'='*70}")

ev_mode = ev[:, mode_idx, 0]
ev_reshaped = ev_mode.reshape(3, xtal.nat)

print("\nMass-weighted eigenvector (from calc_freq_eig):")
for iat in range(xtal.nat):
    ev_atom = ev_reshaped[:, iat]
    print(f"  Atom {iat} ({xtal.symbols[xtal.atom_type_map[iat]-1]}): "
          f"[{ev_atom[0]:8.5f}, {ev_atom[1]:8.5f}, {ev_atom[2]:8.5f}]  "
          f"|e| = {np.linalg.norm(ev_atom):.5f}")

# Check normalization
norm_total_sq = np.sum(np.abs(ev_mode)**2)
print(f"\nTotal norm squared (mass-weighted): {norm_total_sq:.6f} (should be ~1)")

print("\nPhysical displacements (undo mass weighting):")
for iat in range(xtal.nat):
    ev_atom = ev_reshaped[:, iat]
    ev_atom_physical = ev_atom / np.sqrt(masses[iat])
    print(f"  Atom {iat} ({xtal.symbols[xtal.atom_type_map[iat]-1]}, m={masses[iat]:.2f}): "
          f"[{ev_atom_physical[0]:8.5f}, {ev_atom_physical[1]:8.5f}, {ev_atom_physical[2]:8.5f}]  "
          f"|u| = {np.linalg.norm(ev_atom_physical):.5f}")

# Calculate L-char the current way
print("\n" + "="*70)
print("L-char Calculation (Current Method)")
print("="*70)

total_amp_sq = 0.0
long_amp_sq = 0.0

for iat in range(xtal.nat):
    ev_atom = ev_reshaped[:, iat]
    ev_atom_physical = ev_atom / np.sqrt(masses[iat])
    
    # Method 1: Real part of dot product
    atom_amp_sq = np.abs(np.dot(ev_atom_physical, ev_atom_physical))
    total_amp_sq += atom_amp_sq
    
    long_proj = np.dot(Q_hat[0], ev_atom_physical)
    long_amp_sq += np.abs(long_proj)**2
    
    print(f"Atom {iat}: amp²={atom_amp_sq:.6f}, long_proj={long_proj:.6f}, long²={np.abs(long_proj)**2:.6f}")

L_char = long_amp_sq / total_amp_sq
print(f"\nL-char = {long_amp_sq:.6f} / {total_amp_sq:.6f} = {L_char:.3f}")

# Alternative: use real parts explicitly
print("\n" + "="*70)
print("L-char Calculation (Alternative: Real Parts)")
print("="*70)

total_amp_sq_alt = 0.0
long_amp_sq_alt = 0.0

for iat in range(xtal.nat):
    ev_atom = ev_reshaped[:, iat]
    ev_atom_physical = ev_atom / np.sqrt(masses[iat])
    
    # Take real parts
    ev_real = np.real(ev_atom_physical)
    
    atom_amp_sq = np.sum(ev_real**2)
    total_amp_sq_alt += atom_amp_sq
    
    long_proj = np.dot(Q_hat[0], ev_real)
    long_amp_sq_alt += long_proj**2
    
    print(f"Atom {iat}: amp²={atom_amp_sq:.6f}, long_proj={long_proj:.6f}, long²={long_proj**2:.6f}")

L_char_alt = long_amp_sq_alt / total_amp_sq_alt
print(f"\nL-char (alt) = {long_amp_sq_alt:.6f} / {total_amp_sq_alt:.6f} = {L_char_alt:.3f}")

# MATLAB method: sum over all components, then divide
print("\n" + "="*70)
print("L-char Calculation (MATLAB-style)")
print("="*70)

# Collect all displacement components
u_all = []
for iat in range(xtal.nat):
    ev_atom = ev_reshaped[:, iat]
    ev_atom_physical = ev_atom / np.sqrt(masses[iat])
    u_all.extend(np.real(ev_atom_physical))

u_all = np.array(u_all)

# Calculate Q component for each atom
Q_components = []
for iat in range(xtal.nat):
    ev_atom = ev_reshaped[:, iat]
    ev_atom_physical = np.real(ev_atom / np.sqrt(masses[iat]))
    Q_dot = np.dot(Q_hat[0], ev_atom_physical)
    Q_components.append(Q_dot)

Q_components = np.array(Q_components)

print(f"Total displacement norm: {np.linalg.norm(u_all):.6f}")
print(f"Q components: {Q_components}")
print(f"Q component norm: {np.linalg.norm(Q_components):.6f}")

# From MATLAB: long_amp_sq = sum of |Q·u_atom|²
# total_amp_sq = sum of |u_atom|²

L_char_matlab = np.sum(Q_components**2) / np.linalg.norm(u_all)**2
print(f"\nL-char (MATLAB-style) = {L_char_matlab:.3f}")

print("\n" + "="*70)
print("Expected from MATLAB: 0.967")
print("="*70)

