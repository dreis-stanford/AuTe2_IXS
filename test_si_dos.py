"""
Quick test of Silicon phonon DOS
"""

import numpy as np
import matplotlib.pyplot as plt
import sys
sys.path.insert(0, 'code')

from force_constants import ForceConstants
from phonons import calc_Dq, calc_freq_eig, convert_frequencies, generate_bz_grid

print("="*70)
print("Silicon Phonon Density of States")
print("="*70)

# Load force constants
fc_file = "data/Test__Silicon_dispersion/Qgrid_888/Cg.fc"
xtal = ForceConstants(fc_file)
Phi = xtal.convert_to_eV_per_Angstrom2()
masses = np.array([xtal.masses[0], xtal.masses[0]])

print(f"\nCalculating DOS on BZ grid...")

# Generate BZ grid (try different sizes)
for N in [10, 20]:
    print(f"\n{N}x{N}x{N} grid ({N**3} q-points):")
    
    Q_grid = generate_bz_grid(N)
    print(f"  Calculating phonons...")
    D_grid = calc_Dq(Q_grid, xtal.uvw, Phi, masses)
    w_raw, _ = calc_freq_eig(D_grid, sort_modes=False)
    w_meV = convert_frequencies(w_raw, 'meV').flatten()
    
    print(f"  Frequency range: {w_meV.min():.2f} to {w_meV.max():.2f} meV")
    
    # Calculate histogram
    dos, bins = np.histogram(w_meV, bins=100, range=(0, 70))
    bin_centers = (bins[1:] + bins[:-1]) / 2
    
    # Plot
    fig, ax = plt.subplots(figsize=(8, 5))
    ax.plot(bin_centers, dos, '-', linewidth=2, label=f'{N}³ grid')
    ax.set_xlabel('Frequency [meV]', fontsize=12)
    ax.set_ylabel('DOS [states/meV]', fontsize=12)
    ax.set_title(f'Silicon Phonon DOS ({N}x{N}x{N} grid)', fontsize=14)
    ax.grid(True, alpha=0.3)
    ax.set_xlim(0, 70)
    ax.legend()
    
    plt.tight_layout()
    plt.savefig(f'si_dos_{N}x{N}x{N}.png', dpi=200)
    print(f"  Saved: si_dos_{N}x{N}x{N}.png")

print("\n" + "="*70)
print("DOS calculation complete!")
print("Expected features for Si:")
print("  - Gap between acoustic and optic branches")
print("  - Max frequency ~63 meV")
print("  - Van Hove singularities")
print("="*70)

plt.show()
