"""
Test symmetrizing eigenvectors for equivalent atoms
"""

import numpy as np
import sys
sys.path.insert(0, 'code')

from force_constants import ForceConstants
from phonons import calc_Dq, calc_freq_eig, convert_frequencies

def symmetrize_eigenvectors(ev, atom_types, masses):
    """
    Symmetrize eigenvectors so equivalent atoms have consistent phases
    
    For each mode, adjust phases so equivalent atoms move "together"
    
    Parameters:
    -----------
    ev : ndarray (3*nat, nmodes, nQ)
        Eigenvectors
    atom_types : list
        Atom type indices (e.g., [0, 1, 1] for Au, Te, Te)
    masses : ndarray
        Atomic masses
    
    Returns:
    --------
    ev_sym : ndarray
        Symmetrized eigenvectors
    """
    
    ev_sym = ev.copy()
    nat = len(atom_types)
    nmodes = ev.shape[1]
    nQ = ev.shape[2]
    
    for iq in range(nQ):
        for imode in range(nmodes):
            ev_mode = ev[:, imode, iq].reshape(3, nat)
            
            # Find groups of equivalent atoms
            unique_types = list(set(atom_types))
            
            for atom_type in unique_types:
                # Find all atoms of this type
                indices = [i for i, t in enumerate(atom_types) if t == atom_type]
                
                if len(indices) <= 1:
                    continue  # Only one atom of this type
                
                # Get eigenvector components for these atoms
                # Choose first atom as reference
                ref_idx = indices[0]
                ev_ref = ev_mode[:, ref_idx] / np.sqrt(masses[ref_idx])
                
                # For each other atom of same type, find best phase to match reference
                for idx in indices[1:]:
                    ev_atom = ev_mode[:, idx] / np.sqrt(masses[idx])
                    
                    # Try to maximize overlap with reference
                    # Best phase: arg(ev_ref^* · ev_atom)
                    overlap = np.dot(np.conj(ev_ref), ev_atom)
                    
                    if np.abs(overlap) > 1e-10:
                        # Rotate to make overlap real and positive
                        phase_correction = np.conj(overlap) / np.abs(overlap)
                        ev_mode[:, idx] *= phase_correction
            
            # Store back
            ev_sym[:, imode, iq] = ev_mode.flatten()
    
    return ev_sym


# Test it
fc_file = "data/AuTe_2_m.fc"
xtal = ForceConstants(fc_file)
Phi = xtal.convert_to_eV_per_Angstrom2()

masses = np.array([xtal.masses[xtal.atom_type_map[i]-1] 
                   for i in range(xtal.nat)])

atom_types = [xtal.atom_type_map[i]-1 for i in range(xtal.nat)]  # [0, 1, 1] for Au, Te, Te

print("="*70)
print("Testing eigenvector symmetrization")
print("="*70)
print(f"\nAtom types: {atom_types}")
print(f"Symbols: {[xtal.symbols[t] for t in atom_types]}")

# Test at Q = (0.5, 0, 0) where we saw issues
Q_prim = np.array([[0.5, 0.0, 0.0]])

# Calculate with reduced q
Q_reduced = Q_prim - np.round(Q_prim)
Dq = calc_Dq(Q_reduced, xtal.uvw, Phi, masses)
w_raw, ev_orig = calc_freq_eig(Dq, sort_modes=True, rotate_ev=False)
w_cm = convert_frequencies(w_raw, 'cm-1').flatten()

# Symmetrize
ev_sym = symmetrize_eigenvectors(ev_orig, atom_types, masses)

print("\n" + "="*70)
print("Comparing atomic participation before/after symmetrization")
print("="*70)

for imode in [0, 1, 5]:  # Check a few modes
    print(f"\nMode {imode+1}: {w_cm[imode]:.2f} cm⁻¹")
    
    # Original
    ev_mode_orig = ev_orig[:, imode, 0].reshape(3, xtal.nat)
    part_orig = []
    for iat in range(xtal.nat):
        ev_atom = ev_mode_orig[:, iat] / np.sqrt(masses[iat])
        part_orig.append(np.linalg.norm(ev_atom))
    total = sum(part_orig)
    part_orig_pct = [p/total*100 for p in part_orig]
    
    # Symmetrized
    ev_mode_sym = ev_sym[:, imode, 0].reshape(3, xtal.nat)
    part_sym = []
    for iat in range(xtal.nat):
        ev_atom = ev_mode_sym[:, iat] / np.sqrt(masses[iat])
        part_sym.append(np.linalg.norm(ev_atom))
    total = sum(part_sym)
    part_sym_pct = [p/total*100 for p in part_sym]
    
    print(f"  Original:    Au={part_orig_pct[0]:5.1f}%, Te1={part_orig_pct[1]:5.1f}%, Te2={part_orig_pct[2]:5.1f}%")
    print(f"  Symmetrized: Au={part_sym_pct[0]:5.1f}%, Te1={part_sym_pct[1]:5.1f}%, Te2={part_sym_pct[2]:5.1f}%")

print("\n" + "="*70)
print("Does symmetrization help?")
print("="*70)

