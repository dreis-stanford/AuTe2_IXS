"""
Check dynamical matrix symmetry
The two Te atoms should be related by symmetry
"""

import numpy as np
import sys
sys.path.insert(0, 'code')

from force_constants import ForceConstants
from phonons import calc_Dq

# Load structure
fc_file = "data/AuTe_2_m.fc"
xtal = ForceConstants(fc_file)
Phi = xtal.convert_to_eV_per_Angstrom2()

masses = np.array([xtal.masses[xtal.atom_type_map[i]-1] 
                   for i in range(xtal.nat)])

print("="*70)
print("Checking Dynamical Matrix Symmetry at Q=(0,0,0.1)")
print("="*70)

# Test Q point
Q_prim = np.array([[0.0, 0.0, 0.1]])
Dq = calc_Dq(Q_prim, xtal.uvw, Phi, masses)

D = Dq[:, :, 0]  # 9x9 matrix

print(f"\nDynamical matrix shape: {D.shape}")

# Check hermiticity
hermiticity_error = np.max(np.abs(D - D.conj().T))
print(f"Hermiticity error: {hermiticity_error:.3e}")

# Check blocks corresponding to Te atoms
# Atom 0: rows/cols 0:3 (Au)
# Atom 1: rows/cols 3:6 (Te1)
# Atom 2: rows/cols 6:9 (Te2)

print("\n" + "="*70)
print("Comparing Te1-Te1 vs Te2-Te2 blocks")
print("="*70)

Te1_block = D[3:6, 3:6]
Te2_block = D[6:9, 6:9]

print("\nTe1-Te1 block:")
print(Te1_block)

print("\nTe2-Te2 block:")
print(Te2_block)

diff = Te1_block - Te2_block
print("\nDifference (Te1-Te1) - (Te2-Te2):")
print(diff)
print(f"Max difference: {np.max(np.abs(diff)):.3e}")

print("\n" + "="*70)
print("Comparing Au-Te1 vs Au-Te2 blocks")
print("="*70)

Au_Te1_block = D[0:3, 3:6]
Au_Te2_block = D[0:3, 6:9]

print("\nAu-Te1 block:")
print(Au_Te1_block)

print("\nAu-Te2 block:")
print(Au_Te2_block)

diff2 = Au_Te1_block - Au_Te2_block
print("\nDifference Au-Te1 - Au-Te2:")
print(diff2)
print(f"Max difference: {np.max(np.abs(diff2)):.3e}")

# Atomic positions
print("\n" + "="*70)
print("Atomic positions (fractional)")
print("="*70)
for i in range(3):
    symbol = xtal.symbols[xtal.atom_type_map[i]-1]
    pos = xtal.xs[i]
    print(f"Atom {i} ({symbol}): [{pos[0]:.4f}, {pos[1]:.4f}, {pos[2]:.4f}]")

# Check if Te1 and Te2 are related by inversion through origin
Te1_pos = xtal.xs[1]
Te2_pos = xtal.xs[2]
Te2_inverted = -Te1_pos

print(f"\nTe1 position: {Te1_pos}")
print(f"Te2 position: {Te2_pos}")
print(f"-Te1 position: {Te2_inverted}")
print(f"Te2 - (-Te1): {Te2_pos - Te2_inverted}")

# At Q=(0,0,0.1), check phase factors
Q_vec = Q_prim[0]
phase_Te1 = np.exp(-2j * np.pi * np.dot(Q_vec, Te1_pos))
phase_Te2 = np.exp(-2j * np.pi * np.dot(Q_vec, Te2_pos))

print(f"\nPhase factor for Te1: {phase_Te1}")
print(f"Phase factor for Te2: {phase_Te2}")
print(f"Phase ratio Te2/Te1: {phase_Te2/phase_Te1}")

