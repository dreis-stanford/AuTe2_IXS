"""
Plot phonon dispersion for Silicon
Converted from test_silicon.py for dispersion plotting
"""

import numpy as np
import matplotlib.pyplot as plt
from .force_constants import ForceConstants
from .phonons import calc_Dq, calc_freq_eig, convert_frequencies

def plot_si_dispersion(fc_file='data/Test__Silicon_dispersion/Qgrid_888/Cg.fc', 
                       n_points=101, show=True, block=False):
    """
    Plot phonon dispersion along high-symmetry path for Silicon (FCC)
    
    Parameters:
    -----------
    fc_file : str
        Path to force constants file
    n_points : int
        Number of points per segment (default 101 to match QE)
    show : bool
        Show the plot
    block : bool
        Block execution until plot is closed
    """
    
    print("="*70)
    print("Silicon Phonon Dispersion Calculation")
    print("="*70)
    
    # Load force constants
    xtal = ForceConstants(fc_file)
    Phi = xtal.convert_to_eV_per_Angstrom2()
    
    # Silicon has 2 atoms of same type
    masses = np.array([xtal.masses[0], xtal.masses[0]])
    
    print(f"\nSilicon: {xtal.nat} atoms, mass = {masses[0]:.2f} amu")
    print(f"Using force constants from: {fc_file}")
    
    # High-symmetry path: G - K - X - G - L - X - W - L
    # Coordinates in PRIMITIVE reciprocal lattice units
    path = [
        ('G', np.array([0.0, 0.0, 0.0])),
        ('K', np.array([0.75, 0.375, 0.375])),
        ('X', np.array([1.0, 0.5, 0.5])),
        ('G', np.array([1.0, 1.0, 1.0])),      # Equivalent to (0,0,0)
        ('L', np.array([0.5, 0.5, 0.5])),
        ('X', np.array([0.5, 0.5, 0.0])),      # Different X point
        ('W', np.array([0.75, 0.5, 0.25])),
        ('L', np.array([0.5, 0.5, 0.5])),
    ]
    
    # Generate q-points along path
    q_path = []
    labels = [path[0][0]]
    ticks = [0]
    
    for i in range(len(path)-1):
        start = path[i][1]
        end = path[i+1][1]
        
        # Use n_points per segment, avoiding duplication
        n_use = n_points if i == len(path)-2 else n_points-1
        for j in range(n_use):
            t = j / (n_points - 1)
            q_path.append((1-t)*start + t*end)
        
        labels.append(path[i+1][0])
        ticks.append(len(q_path)-1)
    
    q_path = np.array(q_path)
    
    print(f"Path: {' - '.join([p[0] for p in path])}")
    print(f"Total q-points: {len(q_path)}")
    
    # Calculate dynamical matrices
    print("\nCalculating dynamical matrices...")
    D = calc_Dq(q_path, xtal.uvw, Phi, masses)
    
    # Solve for frequencies and eigenvectors
    print("Solving for phonon modes...")
    w_raw, ev = calc_freq_eig(D, sort_modes=True, rotate_ev=True)
    
    # Convert to meV
    w_meV = convert_frequencies(w_raw, 'meV')
    
    print(f"Frequency range: {w_meV.min():.1f} to {w_meV.max():.1f} meV")
    
    # Plot
    print("\nCreating plot...")
    fig, ax = plt.subplots(figsize=(10, 6))
    
    for i in range(6):
        ax.plot(w_meV[:, i], 'b-', linewidth=1.5)
    
    for tick in ticks:
        ax.axvline(tick, color='k', linestyle='--', linewidth=0.5, alpha=0.5)
    
    ax.set_xticks(ticks)
    ax.set_xticklabels(labels)
    ax.set_ylabel('Frequency [meV]', fontsize=12)
    ax.set_xlabel('q-(path)', fontsize=12)
    ax.set_title('Phonon dispersion (Si) [PhX+MatdynX]', fontsize=14)
    ax.set_ylim(0, 70)
    ax.set_xlim(0, len(q_path)-1)
    ax.grid(True, alpha=0.3)
    
    plt.tight_layout()
    plt.savefig('si_dispersion.png', dpi=300)
    print(f"\nSaved plot to: si_dispersion.png")
    
    if show:
        plt.show(block=block)
    
    return fig, ax, q_path, w_meV


if __name__ == "__main__":
    import sys
    import os
    sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
    
    fig, ax, q_path, w_meV = plot_si_dispersion(block=True)
    
    print("\n" + "="*70)
    print("Dispersion calculation complete!")
    print("="*70)
