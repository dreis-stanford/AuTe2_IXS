"""
Check for degeneracies and proper mode combinations
"""

import numpy as np
import sys
sys.path.insert(0, 'code')

from force_constants import ForceConstants
from phonons import calc_Dq, calc_freq_eig, convert_frequencies

# Load structure
fc_file = "data/AuTe_2_m.fc"
xtal = ForceConstants(fc_file)
Phi = xtal.convert_to_eV_per_Angstrom2()

masses = np.array([xtal.masses[xtal.atom_type_map[i]-1] 
                   for i in range(xtal.nat)])

# Test Q point
Q_prim = np.array([[0.0, 0.0, 0.1]])
Q_cart = Q_prim @ xtal.b_l
Q_mag = np.linalg.norm(Q_cart)
Q_hat = Q_cart[0] / Q_mag

print("="*70)
print("Checking for degenerate modes at Q=(0,0,0.1)")
print("="*70)

Dq = calc_Dq(Q_prim, xtal.uvw, Phi, masses)
w_raw, ev = calc_freq_eig(Dq, sort_modes=False, rotate_ev=False)
w_cm = convert_frequencies(w_raw, 'cm-1').flatten()

print("\nFrequencies:")
for i in range(9):
    print(f"  Mode {i+1}: {w_cm[i]:8.2f} cm⁻¹")

# Check frequency gaps
print("\nFrequency gaps:")
for i in range(8):
    gap = w_cm[i+1] - w_cm[i]
    print(f"  Mode {i+1} → {i+2}: {gap:6.2f} cm⁻¹")

# For Q along z, we expect:
# - Modes with motion in xy plane: degenerate (or nearly so)
# - Mode with motion along z: separate

print("\n" + "="*70)
print("Analyzing eigenvector directions (without rotation)")
print("="*70)

for imode in range(9):
    ev_mode = ev[:, imode, 0].reshape(3, 3)
    
    print(f"\nMode {imode+1} ({w_cm[imode]:6.2f} cm⁻¹):")
    
    # Check each atom's displacement direction
    for iat in range(3):
        ev_atom = ev_mode[:, iat] / np.sqrt(masses[iat])
        
        # Project onto z-axis
        z_component = ev_atom[2]
        xy_component = np.sqrt(np.abs(ev_atom[0])**2 + np.abs(ev_atom[1])**2)
        
        symbol = xtal.symbols[xtal.atom_type_map[iat]-1]
        print(f"  {symbol}: z={np.abs(z_component):.4f}, xy={np.abs(xy_component):.4f}")

# Now check mode 3 specifically
print("\n" + "="*70)
print("Mode 3 detailed analysis")
print("="*70)

mode3_idx = 2
ev_mode3 = ev[:, mode3_idx, 0].reshape(3, 3)

print(f"\nMode 3 eigenvector (mass-weighted):")
for iat in range(3):
    ev_atom = ev_mode3[:, iat]
    symbol = xtal.symbols[xtal.atom_type_map[iat]-1]
    print(f"  {symbol}: [{ev_atom[0]:.4f}, {ev_atom[1]:.4f}, {ev_atom[2]:.4f}]")

print(f"\nMode 3 eigenvector (physical):")
for iat in range(3):
    ev_atom = ev_mode3[:, iat] / np.sqrt(masses[iat])
    symbol = xtal.symbols[xtal.atom_type_map[iat]-1]
    amp = np.linalg.norm(ev_atom)
    print(f"  {symbol}: [{np.real(ev_atom[0]):.6f}, {np.real(ev_atom[1]):.6f}, {np.real(ev_atom[2]):.6f}]  |u|={np.abs(amp):.6f}")

# Calculate what L-char SHOULD be for a pure longitudinal mode
print(f"\nQ direction: z-axis (Q_hat ≈ [0, 0, 1])")
print(f"Actual Q_hat: [{Q_hat[0]:.6f}, {Q_hat[1]:.6f}, {Q_hat[2]:.6f}]")

# For pure longitudinal: all atoms move along Q (z-direction)
# For pure transverse: all atoms move perpendicular to Q (xy-plane)

