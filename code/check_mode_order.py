"""
Check if mode ordering matches between Python and MATLAB
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

# Calculate modes WITHOUT sorting
print("="*70)
print("Modes WITHOUT sorting")
print("="*70)

Dq = calc_Dq(Q_prim, xtal.uvw, Phi, masses)
w_raw_unsorted, ev_unsorted = calc_freq_eig(Dq, sort_modes=False, rotate_ev=True)
w_cm_unsorted = convert_frequencies(w_raw_unsorted, 'cm-1').flatten()

print("\nFrequencies (unsorted):")
for i, w in enumerate(w_cm_unsorted):
    print(f"  Mode {i+1}: {w:8.2f} cm⁻¹")

# Calculate modes WITH sorting
print("\n" + "="*70)
print("Modes WITH sorting")
print("="*70)

w_raw_sorted, ev_sorted = calc_freq_eig(Dq, sort_modes=True, rotate_ev=True)
w_cm_sorted = convert_frequencies(w_raw_sorted, 'cm-1').flatten()

print("\nFrequencies (sorted):")
for i, w in enumerate(w_cm_sorted):
    print(f"  Mode {i+1}: {w:8.2f} cm⁻¹")

# Expected from MATLAB
matlab_freqs = [9.86, 11.77, 20.67, 32.28, 33.50, 110.77, 119.55, 123.77, 126.11]

print("\n" + "="*70)
print("MATLAB frequencies:")
for i, w in enumerate(matlab_freqs):
    print(f"  Mode {i+1}: {w:8.2f} cm⁻¹")

print("\n" + "="*70)
print("Checking longitudinal mode at ~20.67 cm⁻¹")
print("="*70)

# Find mode closest to 20.67
target_freq = 20.67
diff = np.abs(w_cm_sorted - target_freq)
mode_idx = np.argmin(diff)

print(f"\nClosest mode: {mode_idx+1} at {w_cm_sorted[mode_idx]:.2f} cm⁻¹")

# Check all modes near 20.67
for i in range(9):
    ev_mode = ev_sorted[:, i, 0].reshape(3, xtal.nat)
    
    # Physical displacements
    amps = []
    for iat in range(xtal.nat):
        ev_atom = ev_mode[:, iat]
        ev_phys = ev_atom / np.sqrt(masses[iat])
        amp = np.linalg.norm(ev_phys)
        amps.append(amp)
    
    total = sum(amps)
    participation = [a/total*100 for a in amps]
    
    print(f"\nMode {i+1} ({w_cm_sorted[i]:6.2f} cm⁻¹): "
          f"Au={participation[0]:5.1f}%, Te1={participation[1]:5.1f}%, Te2={participation[2]:5.1f}%")

print("\n" + "="*70)
print("MATLAB Mode 3: Au=33.8%, Te1=33.7%, Te2=32.5%")
print("="*70)

