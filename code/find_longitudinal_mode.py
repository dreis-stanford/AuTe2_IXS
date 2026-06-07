"""
Find the proper longitudinal mode combination
Modes 1, 2, 3 might need to be linearly combined
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

# Test Q point
Q_prim = np.array([[0.0, 0.0, 0.1]])
Q_cart = Q_prim @ xtal.b_l
Q_mag = np.linalg.norm(Q_cart)
Q_hat = Q_cart[0] / Q_mag

print("="*70)
print("Search for longitudinal mode combination")
print("="*70)

Dq = calc_Dq(Q_prim, xtal.uvw, Phi, masses)
w_raw, ev = calc_freq_eig(Dq, sort_modes=False, rotate_ev=False)
w_cm = convert_frequencies(w_raw, 'cm-1').flatten()

# Try all linear combinations of modes 1, 2, 3
print("\nTrying different linear combinations of modes 1, 2, 3:")
print("Looking for one where all 3 atoms participate and L-char is high\n")

best_L = 0
best_combo = None
best_participation = None

# Try many random combinations
np.random.seed(42)
for trial in range(1000):
    # Random coefficients
    if trial < 9:
        # First try pure modes and simple combinations
        if trial < 3:
            coeffs = np.zeros(3)
            coeffs[trial] = 1.0
        elif trial < 6:
            coeffs = np.zeros(3)
            coeffs[trial-3] = 1.0/np.sqrt(2)
            coeffs[(trial-3+1)%3] = 1.0/np.sqrt(2)
        else:
            coeffs = np.ones(3) / np.sqrt(3)
    else:
        coeffs = np.random.randn(3) + 1j * np.random.randn(3)
        coeffs = coeffs / np.linalg.norm(coeffs)
    
    # Build linear combination
    ev_combined = np.zeros(9, dtype=complex)
    for i in range(3):
        ev_combined += coeffs[i] * ev[:, i, 0]
    
    ev_combined = ev_combined / np.linalg.norm(ev_combined)
    ev_reshaped = ev_combined.reshape(3, 3)
    
    # Calculate participation
    participation = []
    for iat in range(3):
        ev_atom = ev_reshaped[:, iat] / np.sqrt(masses[iat])
        amp = np.linalg.norm(ev_atom)
        participation.append(amp)
    
    total = sum(participation)
    part_pct = [p/total*100 for p in participation]
    
    # Calculate L-char
    total_amp_sq = 0.0
    long_amp_sq = 0.0
    
    for iat in range(3):
        ev_atom = ev_reshaped[:, iat] / np.sqrt(masses[iat])
        atom_amp_sq = np.sum(np.abs(ev_atom)**2)
        total_amp_sq += atom_amp_sq
        
        long_proj = np.dot(Q_hat, ev_atom)
        long_amp_sq += np.abs(long_proj)**2
    
    L_char = long_amp_sq / total_amp_sq if total_amp_sq > 1e-12 else 0
    
    # Check if all atoms participate reasonably
    min_participation = min(part_pct)
    
    if min_participation > 20 and L_char > best_L:
        best_L = L_char
        best_combo = coeffs.copy()
        best_participation = part_pct
        
        if trial < 20 or L_char > 0.9:
            print(f"Trial {trial:4d}: L={L_char:.3f}, "
                  f"Au={part_pct[0]:4.1f}%, Te1={part_pct[1]:4.1f}%, Te2={part_pct[2]:4.1f}%")

print(f"\n" + "="*70)
print(f"Best combination found:")
print(f"  L-char = {best_L:.3f}")
print(f"  Participation: Au={best_participation[0]:.1f}%, "
      f"Te1={best_participation[1]:.1f}%, Te2={best_participation[2]:.1f}%")
print(f"  Coefficients: {best_combo}")
print(f"\n  Target (MATLAB): L=0.967, Au=33.8%, Te1=33.7%, Te2=32.5%")
print("="*70)

