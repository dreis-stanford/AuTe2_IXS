"""
Complete example: AuTe2 IXS analysis with cross-sections
"""

from sixcircle_minimal import init_sixcircle
from .aute2_structure import AuTe2
from ixs_cross_section import create_ixs_calculator
import numpy as np
import matplotlib.pyplot as plt

print("=" * 60)
print("AuTe2 IXS Complete Analysis")
print("=" * 60)

# Setup
crystal = AuTe2()
sc = init_sixcircle(21.747)
sc.setup_crystal(crystal.a, crystal.b, crystal.c,
                 crystal.alpha, crystal.beta, crystal.gamma,
                 crystal.name)

ixs = create_ixs_calculator(crystal, verbose=True)

# Choose a reflection
h, k, l = 1, 0, 0
result = sc.ca(h, k, l)

# Calculate spectrum
print("\n" + "=" * 60)
print(f"IXS spectrum for ({h},{k},{l})")
print("=" * 60)

energies = np.linspace(-20, 120, 200)  # meV
energies, cross_sections = ixs.calculate_spectrum(h, k, l, result['Q'], energies)

# Plot
plt.figure(figsize=(10, 6))
plt.plot(energies, cross_sections, 'b-', linewidth=2)
plt.xlabel('Energy Transfer (meV)', fontsize=12)
plt.ylabel('Cross-section (arb. units)', fontsize=12)
plt.title(f'{crystal.name} IXS: Q=({h},{k},{l}), |Q|={result["Q"]:.3f} Å⁻¹', 
          fontsize=14)
plt.grid(True, alpha=0.3)
plt.axvline(x=50, color='r', linestyle='--', alpha=0.5, label='Phonon energy (placeholder)')
plt.legend()
plt.tight_layout()
plt.savefig('aute2_ixs_spectrum.png', dpi=150)

print("\n✓ Saved: aute2_ixs_spectrum.png")
print(f"✓ Peak cross-section: {np.max(cross_sections):.3e}")
print("\nReady to add your MATLAB physics!")
