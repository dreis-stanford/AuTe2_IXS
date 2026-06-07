"""
Calculate IXS for ALL modes and compare with MATLAB
Find which Python mode corresponds to MATLAB's high-intensity mode 3
"""

import numpy as np
import sys
sys.path.insert(0, 'code')

from .force_constants import ForceConstants
from .phonons import calc_Dq, calc_freq_eig, convert_frequencies
from .form_factors import CalcAtomicfQ
from .ixs import calc_ixs

# Load structure
fc_file = "data/AuTe_2_m.fc"
xtal = ForceConstants(fc_file)
Phi = xtal.convert_to_eV_per_Angstrom2()

masses = np.array([xtal.masses[xtal.atom_type_map[i]-1] 
                   for i in range(xtal.nat)])

Q_prim = np.array([[0.0, 0.0, 0.1]])

# Temperature
kT_cm = 207
kT_THz = kT_cm * 2.99792458e10 / 1e12

# Calculate modes
Dq = calc_Dq(Q_prim, xtal.uvw, Phi, masses)
w_raw, ev = calc_freq_eig(Dq, sort_modes=False, rotate_ev=False)
w_THz = convert_frequencies(w_raw, 'THz').flatten()
w_cm = convert_frequencies(w_raw, 'cm-1').flatten()

# Form factors
Q_cart = Q_prim @ xtal.b_l
Q_mag = np.linalg.norm(Q_cart)
fAu = CalcAtomicfQ(Q_mag, 'Au', scale=4*np.pi, use_xraylib=False)
fTe = CalcAtomicfQ(Q_mag, 'Te', scale=4*np.pi, use_xraylib=False)
fQ_matrix = np.array([[fAu, fTe, fTe]])

# Calculate IXS
Is, Ias, n, F, info = calc_ixs(
    w_THz.reshape(1, -1), ev, Q_prim,
    xtal.b_l, xtal.xs, fQ_matrix, masses, kT_THz,
    units='barn', per_steradian=True
)

print("="*80)
print("Python IXS for all modes")
print("="*80)

print("\nMode  Freq(cm⁻¹)   IXS(S)      IXS(AS)     |F|²")
print("-" * 80)

for i in range(9):
    F_sq = np.abs(F[0, i])**2
    print(f"{i+1:2d}    {w_cm[i]:8.2f}   {Is[0,i]:10.3f}  {Ias[0,i]:10.3f}  {F_sq:10.3e}")

print("\n" + "="*80)
print("MATLAB IXS values")
print("="*80)

matlab_data = [
    (9.86,   "~0",     "~0"),
    (11.77,  8.572,    8.098),
    (20.67,  655.4,    593.1),
    (32.28,  "~0",     "~0"),
    (33.50,  8.590,    7.306),
    (110.77, 1.609,    0.942),
    (119.55, 1.765,    0.990),
    (123.77, "~0",     "~0"),
    (126.11, 7.407,    4.028),
]

print("\nMode  Freq(cm⁻¹)   IXS(S)      IXS(AS)")
print("-" * 60)
for i, (freq, ixs_s, ixs_as) in enumerate(matlab_data):
    if isinstance(ixs_s, str):
        print(f"{i+1:2d}    {freq:8.2f}   {ixs_s:>10s}  {ixs_as:>10s}")
    else:
        print(f"{i+1:2d}    {freq:8.2f}   {ixs_s:10.3f}  {ixs_as:10.3f}")

print("\n" + "="*80)
print("Searching for matches...")
print("="*80)

# Find which Python mode matches each MATLAB mode
matlab_freqs = [9.86, 11.77, 20.67, 32.28, 33.50, 110.77, 119.55, 123.77, 126.11]
matlab_ixs = [0, 8.572, 655.4, 0, 8.590, 1.609, 1.765, 0, 7.407]

for i, (mat_freq, mat_ixs) in enumerate(zip(matlab_freqs, matlab_ixs)):
    # Find Python mode with matching frequency
    py_idx = np.argmin(np.abs(w_cm - mat_freq))
    py_freq = w_cm[py_idx]
    py_ixs = Is[0, py_idx]
    
    match_freq = np.abs(py_freq - mat_freq) < 0.1
    match_ixs = np.abs(py_ixs - mat_ixs) < 1.0 if mat_ixs > 1 else True
    
    status = "✓" if (match_freq and match_ixs) else "✗"
    
    print(f"{status} MATLAB Mode {i+1} ({mat_freq:.2f} cm⁻¹, IXS={mat_ixs:.1f})")
    print(f"  → Python Mode {py_idx+1} ({py_freq:.2f} cm⁻¹, IXS={py_ixs:.1f})")
    
    if not match_ixs and mat_ixs > 1:
        print(f"  ⚠ IXS MISMATCH: {py_ixs:.1f} vs {mat_ixs:.1f}")

