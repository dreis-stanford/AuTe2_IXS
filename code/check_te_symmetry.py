"""
Check if Te atoms are truly equivalent at different Q points
"""

import numpy as np
import sys
sys.path.insert(0, 'code')

from force_constants import ForceConstants

fc_file = "data/AuTe_2_m.fc"
xtal = ForceConstants(fc_file)

print("="*70)
print("Checking Te atom symmetry")
print("="*70)

Te1_pos = xtal.xs[1]
Te2_pos = xtal.xs[2]

print(f"\nTe1 position: [{Te1_pos[0]:.4f}, {Te1_pos[1]:.4f}, {Te1_pos[2]:.4f}]")
print(f"Te2 position: [{Te2_pos[0]:.4f}, {Te2_pos[1]:.4f}, {Te2_pos[2]:.4f}]")

# Check inversion symmetry
print(f"\n-Te1:         [{-Te1_pos[0]:.4f}, {-Te1_pos[1]:.4f}, {-Te1_pos[2]:.4f}]")
print(f"Te2 - (-Te1): [{Te2_pos[0]+Te1_pos[0]:.4f}, {Te2_pos[1]+Te1_pos[1]:.4f}, {Te2_pos[2]+Te1_pos[2]:.4f}]")

# At different Q points, check phase factors
test_Q = [
    [0.5, 0.0, 0.0],
    [0.0, 0.0, 0.1],
    [0.0, 0.0, 0.0],  # Gamma
]

print("\n" + "="*70)
print("Phase factors at different Q:")
print("="*70)

for Q in test_Q:
    Q_arr = np.array(Q)
    phase1 = np.exp(-2j * np.pi * np.dot(Q_arr, Te1_pos))
    phase2 = np.exp(-2j * np.pi * np.dot(Q_arr, Te2_pos))
    
    print(f"\nQ = [{Q[0]:.1f}, {Q[1]:.1f}, {Q[2]:.1f}]:")
    print(f"  phase(Te1) = {phase1:.4f}")
    print(f"  phase(Te2) = {phase2:.4f}")
    print(f"  phase2/phase1 = {phase2/phase1:.4f}")
    print(f"  |phase2/phase1| = {np.abs(phase2/phase1):.4f}")
    
    # Are they equivalent?
    phase_ratio = phase2/phase1
    is_real = np.abs(np.imag(phase_ratio)) < 0.01
    is_plus_minus_1 = np.abs(np.abs(phase_ratio) - 1.0) < 0.01
    
    if is_plus_minus_1:
        if is_real and np.real(phase_ratio) > 0:
            print("  → Te atoms IN PHASE (should move together)")
        elif is_real and np.real(phase_ratio) < 0:
            print("  → Te atoms OUT OF PHASE (should move oppositely)")
        else:
            print("  → Phase ratio is ±1 but complex")
    else:
        print("  → Te atoms NOT equivalent at this Q (different phases)")

print("\n" + "="*70)
print("Conclusion:")
print("="*70)
print("The Te atoms are only equivalent at special Q points where")
print("the phase ratio is exactly ±1. At generic Q, they break symmetry.")
print("\nMATLAB must be doing additional eigenvector rotation that we're not.")

