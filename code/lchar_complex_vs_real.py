"""
Compare L-char calculation using complex vs real eigenvectors
MATLAB might use a different convention
"""

import numpy as np
import sys
sys.path.insert(0, 'code')

from .force_constants import ForceConstants
from .phonons import calc_Dq, calc_freq_eig, convert_frequencies

# Load structure
fc_file = "data/AuTe_2_m.fc"
xtal = ForceConstants(fc_file)
Phi = xtal.convert_to_eV_per_Angstrom2()

masses = np.array([xtal.masses[xtal.atom_type_map[i]-1] 
                   for i in range(xtal.nat)])

Q_prim = np.array([[0.0, 0.0, 0.1]])
Q_cart = Q_prim @ xtal.b_l
Q_mag = np.linalg.norm(Q_cart)
Q_hat = Q_cart[0] / Q_mag

Dq = calc_Dq(Q_prim, xtal.uvw, Phi, masses)
w_raw, ev = calc_freq_eig(Dq, sort_modes=False, rotate_ev=False)
w_cm = convert_frequencies(w_raw, 'cm-1').flatten()

print("="*70)
print("Comparing L-char calculations for Mode 3 (20.67 cm⁻¹)")
print("="*70)

mode_idx = 2
ev_mode = ev[:, mode_idx, 0].reshape(3, 3)

print("\nMethod 1: Using REAL parts of eigenvectors")
print("-" * 70)

total_amp_sq = 0.0
long_amp_sq = 0.0

for iat in range(3):
    ev_atom = ev_mode[:, iat] / np.sqrt(masses[iat])
    ev_real = np.real(ev_atom)
    
    atom_amp_sq = np.sum(ev_real**2)
    total_amp_sq += atom_amp_sq
    
    long_proj = np.dot(Q_hat, ev_real)
    long_amp_sq += long_proj**2
    
    symbol = xtal.symbols[xtal.atom_type_map[iat]-1]
    print(f"{symbol}: u_real = [{ev_real[0]:.6f}, {ev_real[1]:.6f}, {ev_real[2]:.6f}], "
          f"Q·u = {long_proj:.6f}")

L_char_real = long_amp_sq / total_amp_sq
print(f"\nL-char (real parts) = {L_char_real:.3f}")

print("\n" + "="*70)
print("Method 2: Using COMPLEX eigenvectors (absolute values)")
print("-" * 70)

total_amp_sq = 0.0
long_amp_sq = 0.0

for iat in range(3):
    ev_atom = ev_mode[:, iat] / np.sqrt(masses[iat])
    
    # Use absolute values
    atom_amp_sq = np.sum(np.abs(ev_atom)**2)
    total_amp_sq += atom_amp_sq
    
    # Complex dot product
    long_proj = np.dot(Q_hat, ev_atom)
    long_amp_sq += np.abs(long_proj)**2
    
    symbol = xtal.symbols[xtal.atom_type_map[iat]-1]
    print(f"{symbol}: |u| = {np.linalg.norm(ev_atom):.6f}, |Q·u| = {np.abs(long_proj):.6f}")

L_char_complex = long_amp_sq / total_amp_sq
print(f"\nL-char (complex, abs values) = {L_char_complex:.3f}")

print("\n" + "="*70)
print("Method 3: Phase-averaged (rotate eigenvector through all phases)")
print("-" * 70)

# Try rotating eigenvector by different phases
phases = np.linspace(0, 2*np.pi, 100)
L_chars = []

for phase in phases:
    ev_rotated = ev[:, mode_idx, 0] * np.exp(1j * phase)
    ev_rotated_reshaped = ev_rotated.reshape(3, 3)
    
    total_amp_sq = 0.0
    long_amp_sq = 0.0
    
    for iat in range(3):
        ev_atom = ev_rotated_reshaped[:, iat] / np.sqrt(masses[iat])
        ev_real = np.real(ev_atom)
        
        atom_amp_sq = np.sum(ev_real**2)
        total_amp_sq += atom_amp_sq
        
        long_proj = np.dot(Q_hat, ev_real)
        long_amp_sq += long_proj**2
    
    L_char = long_amp_sq / total_amp_sq
    L_chars.append(L_char)

print(f"L-char range over all phases: [{min(L_chars):.3f}, {max(L_chars):.3f}]")
print(f"Average L-char: {np.mean(L_chars):.3f}")

# Find phase that gives maximum L-char
max_idx = np.argmax(L_chars)
best_phase = phases[max_idx]
print(f"\nBest phase: {best_phase:.3f} rad ({np.degrees(best_phase):.1f}°)")
print(f"L-char at best phase: {L_chars[max_idx]:.3f}")

# Show eigenvector at best phase
ev_best = ev[:, mode_idx, 0] * np.exp(1j * best_phase)
ev_best_reshaped = ev_best.reshape(3, 3)

print("\nEigenvector at best phase (real parts):")
for iat in range(3):
    ev_atom = ev_best_reshaped[:, iat] / np.sqrt(masses[iat])
    ev_real = np.real(ev_atom)
    symbol = xtal.symbols[xtal.atom_type_map[iat]-1]
    print(f"  {symbol}: [{ev_real[0]:.6f}, {ev_real[1]:.6f}, {ev_real[2]:.6f}]")

print("\n" + "="*70)
print(f"MATLAB reports: L-char = 0.967")
print("="*70)

