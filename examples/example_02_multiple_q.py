#!/usr/bin/env python
"""
Example 2: Compare IXS at multiple Q-points

This example calculates IXS along a path in reciprocal space
and compares the intensities.
"""

import numpy as np
import matplotlib.pyplot as plt
import sys
import os

sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', 'code'))

from single_q_analysis import SingleQAnalyzer
from force_constants import ForceConstants

def main():
    print("="*70)
    print("Example 2: Multiple Q-points Analysis")
    print("="*70)
    
    # Load structure
    xtal = ForceConstants("../data/AuTe_2_m.fc")
    Phi = xtal.convert_to_eV_per_Angstrom2()
    masses = np.array([xtal.masses[xtal.atom_type_map[i]-1] 
                       for i in range(xtal.nat)])
    
    kT_THz = 6.2  # ~300 K
    analyzer = SingleQAnalyzer(xtal, Phi, masses, kT_THz)
    
    # Define path in reciprocal space: Γ → X
    # Γ = (0, 0, 0), X = (0.5, 0, 0) in primitive coordinates
    n_points = 11
    q_path = np.linspace(0, 0.5, n_points)
    
    print(f"\nCalculating along Γ → X ({n_points} points)")
    print("-"*70)
    
    # Store results
    all_freqs = []
    all_ixs = []
    q_magnitudes = []
    
    for q in q_path:
        Q = [q, 0.0, 0.0]
        
        # Analyze without printing
        result = analyzer.analyze(Q, coords='primitive', print_results=False)
        
        all_freqs.append(result['frequencies_cm'])
        all_ixs.append(result['IXS_stokes'])
        q_magnitudes.append(result['Q_mag'])
        
        print(f"  q = [{q:.3f}, 0, 0]: |Q| = {result['Q_mag']:.4f} (2π/Å)")
    
    # Convert to arrays
    all_freqs = np.array(all_freqs)  # Shape: (n_points, n_modes)
    all_ixs = np.array(all_ixs)
    
    # Plot dispersion with IXS intensity
    print("\nPlotting results...")
    
    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(10, 8), sharex=True)
    
    # Plot 1: Phonon dispersion
    n_modes = all_freqs.shape[1]
    for i in range(n_modes):
        ax1.plot(q_path, all_freqs[:, i], 'o-', markersize=4)
    
    ax1.set_ylabel('Frequency (cm⁻¹)', fontsize=12)
    ax1.set_title('Phonon Dispersion: Γ → X', fontsize=14)
    ax1.grid(True, alpha=0.3)
    ax1.set_xlim(0, 0.5)
    
    # Plot 2: IXS intensity
    for i in range(n_modes):
        # Size proportional to IXS intensity
        sizes = all_ixs[:, i] * 10  # Scale for visibility
        ax2.scatter(q_path, all_freqs[:, i], s=sizes, alpha=0.6)
    
    ax2.set_xlabel('q along [h, 0, 0] (r.l.u.)', fontsize=12)
    ax2.set_ylabel('Frequency (cm⁻¹)', fontsize=12)
    ax2.set_title('IXS Intensity (marker size ∝ intensity)', fontsize=14)
    ax2.grid(True, alpha=0.3)
    
    plt.tight_layout()
    
    output_file = 'ixs_dispersion_gamma_x.png'
    plt.savefig(output_file, dpi=150, bbox_inches='tight')
    print(f"Saved: {output_file}")
    
    # Print summary statistics
    print("\n" + "="*70)
    print("Summary Statistics")
    print("="*70)
    
    total_ixs = np.sum(all_ixs, axis=1)
    max_total_idx = np.argmax(total_ixs)
    
    print(f"\nTotal IXS intensity vs q:")
    for i, q in enumerate(q_path):
        print(f"  q = {q:.3f}: Total IXS = {total_ixs[i]:.2f} barn")
    
    print(f"\nMaximum total intensity at q = {q_path[max_total_idx]:.3f}")
    
    # Find which mode dominates
    for i, q in enumerate([0.0, 0.25, 0.5]):
        idx = np.argmin(np.abs(q_path - q))
        dominant_mode = np.argmax(all_ixs[idx, :])
        
        print(f"\nAt q = {q:.2f}:")
        print(f"  Dominant mode: {dominant_mode + 1}")
        print(f"  Frequency: {all_freqs[idx, dominant_mode]:.2f} cm⁻¹")
        print(f"  IXS: {all_ixs[idx, dominant_mode]:.3f} barn")
    
    print("\n" + "="*70)
    print("Example complete!")
    print("="*70)
    
    plt.show()

if __name__ == "__main__":
    main()
