#!/usr/bin/env python3
"""
Test that analyzer array calculations match individual Q-point calculations
"""

import numpy as np
import sys
import os

# Add code directory to path
sys.path.insert(0, os.path.join(os.path.dirname(__file__), 'code'))

# Now import with absolute imports
from code.force_constants import ForceConstants
from code.single_q_analysis import SingleQAnalyzer
from code.constants import const

print("="*80)
print("Testing Analyzer Array vs Individual Calculations")
print("="*80)

# Load force constants
fc_file = "data/AuTe_2_m.fc"
print(f"\nLoading {fc_file}...")
xtal = ForceConstants(fc_file)
Phi = xtal.convert_to_eV_per_Angstrom2()
masses = np.array([xtal.masses[xtal.atom_type_map[i]-1] for i in range(xtal.nat)])

# Temperature
kT_cm = 207
kT_THz = kT_cm * const.c * 100 / 1e12

# Create analyzer
analyzer = SingleQAnalyzer(xtal, Phi, masses, kT_THz)

# Test Q point
Q_center = np.array([1.0, 0.0, 0.0])

print(f"\nTest Q (center): [{Q_center[0]:.1f}, {Q_center[1]:.1f}, {Q_center[2]:.1f}]")
print("\nCalculating analyzer array...")

# Get array results
array_results = analyzer.analyze_array(Q_center, coords='conventional', 
                                       print_results=False)

print(f"✓ Got {len(array_results)} analyzer positions\n")

# Now test a few individual analyzers
test_analyzers = ['a06', 'a03', 'a20', 'a35']  # Center, corners

print("="*80)
print("Comparing Individual vs Array Calculations")
print("="*80)
print(f"\n{'Ana':4s}  {'Δfreq (meV)':>12s}  {'ΔQ_H':>8s}  {'ΔQ_K':>8s}  {'ΔQ_L':>8s}  {'ΔIXS':>8s}")
print("-"*80)

max_freq_diff = 0.0
max_ixs_diff = 0.0
max_q_diff = 0.0

for ana_name in test_analyzers:
    # Find this analyzer in array results
    array_result = next(r for r in array_results if r['analyzer'] == ana_name)
    
    # Calculate individually
    Q_ana = array_result['Q_conv']
    individual_result = analyzer.analyze(Q_ana, coords='conventional', 
                                        print_results=False)
    
    # Compare frequencies (use first 3 modes)
    freq_array = array_result['frequencies_meV'][:3]
    freq_indiv = individual_result['frequencies_meV'][:3]
    freq_diff = np.abs(freq_array - freq_indiv)
    max_freq_error = np.max(freq_diff)
    
    # Compare Q positions
    Q_diff = np.abs(array_result['Q_conv'] - individual_result['Q_conv'])
    max_q_error = np.max(Q_diff)
    
    # Compare IXS (first mode)
    ixs_array = array_result['IXS_stokes'][0]
    ixs_indiv = individual_result['IXS_stokes'][0]
    ixs_diff = np.abs(ixs_array - ixs_indiv)
    
    # Track maximums
    max_freq_diff = max(max_freq_diff, max_freq_error)
    max_ixs_diff = max(max_ixs_diff, ixs_diff)
    max_q_diff = max(max_q_diff, max_q_error)
    
    print(f"{ana_name:4s}  {max_freq_error:12.6f}  "
          f"{Q_diff[0]:8.6f}  {Q_diff[1]:8.6f}  {Q_diff[2]:8.6f}  "
          f"{ixs_diff:8.6f}")

print("="*80)
print("\nMaximum Differences:")
print(f"  Frequency: {max_freq_diff:.6f} meV")
print(f"  Q position: {max_q_diff:.6f} r.l.u.")
print(f"  IXS cross-section: {max_ixs_diff:.6f} barn/uc·sr")

# Check if differences are acceptable (should be essentially zero)
tolerance_freq = 1e-6  # meV
tolerance_q = 1e-10    # r.l.u.
tolerance_ixs = 1e-6   # barn

print("\nTest Results:")
if max_freq_diff < tolerance_freq:
    print(f"  ✓ Frequencies match (< {tolerance_freq} meV)")
else:
    print(f"  ✗ Frequency mismatch! ({max_freq_diff:.2e} meV)")

if max_q_diff < tolerance_q:
    print(f"  ✓ Q positions match (< {tolerance_q} r.l.u.)")
else:
    print(f"  ✗ Q position mismatch! ({max_q_diff:.2e} r.l.u.)")

if max_ixs_diff < tolerance_ixs:
    print(f"  ✓ IXS cross-sections match (< {tolerance_ixs} barn)")
else:
    print(f"  ✗ IXS mismatch! ({max_ixs_diff:.2e} barn)")

if max_freq_diff < tolerance_freq and max_q_diff < tolerance_q and max_ixs_diff < tolerance_ixs:
    print("\n" + "="*80)
    print("✓ ALL TESTS PASSED - Array and individual calculations agree!")
    print("="*80)
    sys.exit(0)
else:
    print("\n" + "="*80)
    print("✗ TESTS FAILED - Discrepancies found")
    print("="*80)
    sys.exit(1)
