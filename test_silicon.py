"""
Test loading Silicon force constants
"""

import numpy as np
import sys
sys.path.insert(0, 'code')

from force_constants import ForceConstants
from phonons import calc_Dq, calc_freq_eig, convert_frequencies

print("="*70)
print("Testing Silicon Force Constants")
print("="*70)

# Try loading the 4x4x4 supercell
fc_files = [
    "data/Test__Silicon_dispersion/Qgrid_444/Cg.fc",
    "data/Test__Silicon_dispersion/Qgrid_888/Cg.fc"
]

for fc_file in fc_files:
    print(f"\n{'='*70}")
    print(f"Loading: {fc_file}")
    print('='*70)
    
    try:
        xtal = ForceConstants(fc_file)
        
        print(f"\nCrystal structure:")
        print(f"  Atoms in primitive cell: {xtal.nat}")
        print(f"  Atom types: {xtal.ntyp}")
        print(f"  Symbols: {xtal.symbols}")
        print(f"  Masses: {xtal.masses} amu")
        
        print(f"\nLattice parameters:")
        print(f"  a = {xtal.a:.4f} Å")
        print(f"  b = {xtal.b:.4f} Å") 
        print(f"  c = {xtal.c:.4f} Å")
        print(f"  α = {xtal.alpha:.2f}°")
        print(f"  β = {xtal.beta:.2f}°")
        print(f"  γ = {xtal.gamma:.2f}°")
        
        print(f"\nAtomic positions (fractional):")
        for i in range(xtal.nat):
            symbol = xtal.symbols[xtal.atom_type_map[i]-1]
            pos = xtal.xs[i]
            print(f"  Atom {i+1} ({symbol}): [{pos[0]:7.4f}, {pos[1]:7.4f}, {pos[2]:7.4f}]")
        
        # Convert force constants
        Phi = xtal.convert_to_eV_per_Angstrom2()
        
        # Get masses
        masses = np.array([xtal.masses[xtal.atom_type_map[i]-1] 
                          for i in range(xtal.nat)])
        
        # Test phonons at Gamma point
        print(f"\nTesting phonons at Γ point (0,0,0):")
        Q_gamma = np.array([[0, 0, 0]])
        Dq = calc_Dq(Q_gamma, xtal.uvw, Phi, masses)
        w_raw, ev = calc_freq_eig(Dq, sort_modes=True)
        w_THz = convert_frequencies(w_raw, 'THz').flatten()
        w_cm = convert_frequencies(w_raw, 'cm-1').flatten()
        
        print(f"\n  {'Mode':>4} {'ω (THz)':>10} {'ω (cm⁻¹)':>12}")
        print(f"  {'-'*28}")
        for i in range(len(w_cm)):
            print(f"  {i+1:4d} {w_THz[i]:10.4f} {w_cm[i]:12.2f}")
        
        print(f"\n  ✓ First 3 modes should be ~0 (acoustic)")
        
        # Test at X point
        print(f"\nTesting phonons at X point (0.5,0,0):")
        Q_X = np.array([[0.5, 0, 0]])
        Dq = calc_Dq(Q_X, xtal.uvw, Phi, masses)
        w_raw, ev = calc_freq_eig(Dq, sort_modes=True)
        w_THz = convert_frequencies(w_raw, 'THz').flatten()
        w_cm = convert_frequencies(w_raw, 'cm-1').flatten()
        
        print(f"\n  {'Mode':>4} {'ω (THz)':>10} {'ω (cm⁻¹)':>12}")
        print(f"  {'-'*28}")
        for i in range(len(w_cm)):
            print(f"  {i+1:4d} {w_THz[i]:10.4f} {w_cm[i]:12.2f}")
        
        print(f"\n✓ Successfully loaded and tested {fc_file}")
        
    except Exception as e:
        print(f"\n✗ Error loading {fc_file}:")
        print(f"  {e}")
        import traceback
        traceback.print_exc()

print("\n" + "="*70)
print("Test complete!")
print("="*70)
