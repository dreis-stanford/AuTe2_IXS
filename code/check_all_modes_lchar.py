"""
Calculate L-char for ALL modes to see which one is truly longitudinal
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
print("L-character for ALL modes (no eigenvector rotation)")
print("="*70)

Dq = calc_Dq(Q_prim, xtal.uvw, Phi, masses)
w_raw, ev = calc_freq_eig(Dq, sort_modes=False, rotate_ev=False)
w_cm = convert_frequencies(w_raw, 'cm-1').flatten()

print(f"\nQ_hat = [{Q_hat[0]:.6f}, {Q_hat[1]:.6f}, {Q_hat[2]:.6f}]")
print("(nearly pure z-direction)\n")

print("Mode  Freq(cm⁻¹)   L-char   Au%   Te1%  Te2%  Description")
print("-" * 70)

for imode in range(9):
    ev_mode = ev[:, imode, 0].reshape(3, 3)
    
    # Calculate participation
    participation = []
    for iat in range(3):
        ev_atom = ev_mode[:, iat] / np.sqrt(masses[iat])
        amp = np.linalg.norm(ev_atom)
        participation.append(amp)
    
    total = sum(participation)
    part_pct = [p/total*100 if total > 1e-12 else 0 for p in participation]
    
    # Calculate L-char (use REAL parts like MATLAB probably does)
    total_amp_sq = 0.0
    long_amp_sq = 0.0
    
    for iat in range(3):
        ev_atom = ev_mode[:, iat] / np.sqrt(masses[iat])
        ev_real = np.real(ev_atom)  # Take real part
        
        atom_amp_sq = np.sum(ev_real**2)
        total_amp_sq += atom_amp_sq
        
        long_proj = np.dot(Q_hat, ev_real)
        long_amp_sq += long_proj**2
    
    L_char = long_amp_sq / total_amp_sq if total_amp_sq > 1e-12 else 0
    
    # Classify
    if L_char < 0.3:
        desc = "Transverse"
    elif L_char > 0.7:
        desc = "Longitudinal"
    else:
        desc = "Mixed"
    
    print(f"{imode+1:2d}    {w_cm[imode]:8.2f}    {L_char:6.3f}  "
          f"{part_pct[0]:4.1f}  {part_pct[1]:4.1f}  {part_pct[2]:4.1f}  {desc}")

print("\n" + "="*70)
print("MATLAB Mode 3: 20.67 cm⁻¹, L=0.967, Au=33.8%, Te1=33.7%, Te2=32.5%")
print("="*70)

# Now check: which mode has BOTH high L-char AND equal participation?
print("\nLooking for mode with L>0.7 AND all atoms participating ~equally:")
for imode in range(9):
    ev_mode = ev[:, imode, 0].reshape(3, 3)
    
    participation = []
    for iat in range(3):
        ev_atom = ev_mode[:, iat] / np.sqrt(masses[iat])
        amp = np.linalg.norm(ev_atom)
        participation.append(amp)
    
    total = sum(participation)
    part_pct = [p/total*100 if total > 1e-12 else 0 for p in participation]
    
    total_amp_sq = 0.0
    long_amp_sq = 0.0
    
    for iat in range(3):
        ev_atom = ev_mode[:, iat] / np.sqrt(masses[iat])
        ev_real = np.real(ev_atom)
        
        atom_amp_sq = np.sum(ev_real**2)
        total_amp_sq += atom_amp_sq
        
        long_proj = np.dot(Q_hat, ev_real)
        long_amp_sq += long_proj**2
    
    L_char = long_amp_sq / total_amp_sq if total_amp_sq > 1e-12 else 0
    
    min_part = min(part_pct)
    max_part = max(part_pct)
    
    if L_char > 0.7 and min_part > 25:
        print(f"  → Mode {imode+1}: {w_cm[imode]:.2f} cm⁻¹, L={L_char:.3f}, "
              f"participation range {min_part:.1f}-{max_part:.1f}%")

