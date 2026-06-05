#!/usr/bin/env python
"""
Example 1: Analyze IXS at a single Q-point

This example shows how to calculate IXS cross-sections at a specific
Q-vector and access the results programmatically.
"""

import numpy as np
import sys
import os

# Add code directory to path
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', 'code'))

from single_q_analysis import SingleQAnalyzer
from force_constants import ForceConstants
from constants import const

def main():
    print("="*70)
    print("Example 1: Single Q-point IXS Analysis")
    print("="*70)
    
    # Load force constants
    fc_file = "../data/AuTe_2_m.fc"
    print(f"\nLoading force constants from: {fc_file}")
    
    xtal = ForceConstants(fc_file)
    Phi = xtal.convert_to_eV_per_Angstrom2()
    
    # Get masses for each atom in the basis
    masses = np.array([xtal.masses[xtal.atom_type_map[i]-1] 
                       for i in range(xtal.nat)])
    
    # Set temperature
    T_kelvin = 300  # Room temperature
    kT_meV = const.kB * T_kelvin  # in eV
    kT_THz = kT_meV / const.THz2meV * 1000  # Convert to THz
    
    print(f"Temperature: {T_kelvin} K = {kT_THz:.2f} THz")
    
    # Create analyzer
    analyzer = SingleQAnalyzer(xtal, Phi, masses, kT_THz)
    
    # Analyze at Q = (0.5, 0, 0) in primitive coordinates
    Q_prim = [0.5, 0.0, 0.0]
    
    print(f"\nAnalyzing Q = {Q_prim} (primitive coordinates)")
    print("-"*70)
    
    # Run analysis (with printed output)
    result = analyzer.analyze(Q_prim, coords='primitive', print_results=True)
    
    # Access results programmatically
    print("\n" + "="*70)
    print("Accessing Results Programmatically")
    print("="*70)
    
    print(f"\nQ (primitive):    {result['Q_prim']}")
    print(f"Q (conventional): {result['Q_conv']}")
    print(f"Q (Cartesian):    {result['Q_cart']} (2π/Å)")
    print(f"|Q| = {result['Q_mag']:.4f} (2π/Å)")
    
    print(f"\nNumber of modes: {len(result['frequencies_cm'])}")
    
    # Find mode with highest IXS intensity
    max_ixs_idx = np.argmax(result['IXS_stokes'])
    
    print(f"\nStrongest IXS mode:")
    print(f"  Mode {max_ixs_idx + 1}")
    print(f"  Frequency: {result['frequencies_cm'][max_ixs_idx]:.2f} cm⁻¹")
    print(f"            {result['frequencies_meV'][max_ixs_idx]:.2f} meV")
    print(f"  IXS (Stokes):      {result['IXS_stokes'][max_ixs_idx]:.3f} barn")
    print(f"  IXS (Anti-Stokes): {result['IXS_antistokes'][max_ixs_idx]:.3f} barn")
    print(f"  L-character: {result['long_char'][max_ixs_idx]:.3f}")
    print(f"  Polarization: {result['pol_type'][max_ixs_idx]}")
    
    # List all modes with significant IXS intensity
    threshold = 0.1  # barn
    significant = result['IXS_stokes'] > threshold
    
    print(f"\nModes with IXS > {threshold} barn:")
    for i in np.where(significant)[0]:
        print(f"  Mode {i+1}: {result['frequencies_cm'][i]:6.2f} cm⁻¹, "
              f"IXS = {result['IXS_stokes'][i]:6.3f} barn")
    
    print("\n" + "="*70)
    print("Example complete!")
    print("="*70)

if __name__ == "__main__":
    main()
