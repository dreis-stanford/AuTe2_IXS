"""
Parser for force constants file (AuTe_2_m.fc format)
Calculates phonon dispersions for AuTe2
"""

import numpy as np
from scipy import linalg

class ForceConstants:
    """Parse and use force constants for phonon calculations"""
    
    def __init__(self, filename):
        """
        Load force constants from file
        
        Parameters:
        -----------
        filename : str
            Path to .fc file
        """
        self.parse_file(filename)
        
    def parse_file(self, filename):
        """Parse the force constants file"""
        
        with open(filename, 'r') as f:
            lines = f.readlines()
        
        # First line: number of atoms and supercell
        parts = lines[0].split()
        self.n_atoms_prim = int(parts[0])  # Atoms in primitive cell
        self.n_atoms_sc = int(parts[1])     # Atoms in supercell
        
        # Second line: lattice vectors (Angstroms)
        parts = lines[1].split()
        self.a = float(parts[0])
        self.b = float(parts[1])
        self.c = float(parts[2])
        self.alpha = float(parts[3])
        self.beta = float(parts[4])
        self.gamma = float(parts[5])
        
        print(f"Lattice: a={self.a:.4f} b={self.b:.4f} c={self.c:.4f} Å")
        print(f"Angles: α={self.alpha:.2f}° β={self.beta:.2f}° γ={self.gamma:.2f}°")
        
        # Reciprocal lattice vectors (lines 2-4)
        self.recip_vecs = []
        for i in range(2, 5):
            vec = [float(x) for x in lines[i].split()]
            self.recip_vecs.append(np.array(vec))
        
        # Atomic species and masses
        idx = 5
        self.species = []
        self.masses = []  # in atomic mass units
        
        for i in range(self.n_atoms_prim):
            parts = lines[idx].split()
            species_name = parts[1].strip("'")
            mass = float(parts[2])  # This is mass * some factor
            self.species.append(species_name)
            self.masses.append(mass)
            idx += 1
        
        print(f"Species: {self.species}")
        print(f"Masses: {self.masses}")
        
        # Atomic positions in fractional coordinates
        self.positions = []
        for i in range(self.n_atoms_sc):
            parts = lines[idx].split()
            pos = [float(parts[1]), float(parts[2]), float(parts[3])]
            self.positions.append(np.array(pos))
            idx += 1
        
        print(f"Number of atoms in supercell: {len(self.positions)}")
        
        # Skip F line
        idx += 1
        
        # Force constants tensor
        # This is complex - store raw data for now
        self.fc_data = lines[idx:]
        
        print(f"✓ Force constants file parsed successfully")
        
    def get_structure_info(self):
        """Return structure information"""
        return {
            'a': self.a,
            'b': self.b, 
            'c': self.c,
            'alpha': self.alpha,
            'beta': self.beta,
            'gamma': self.gamma,
            'species': self.species,
            'masses': self.masses,
            'positions': self.positions
        }

if __name__ == "__main__":
    import sys
    
    # Test parsing
    fc_file = "../data/AuTe_2_m.fc"
    
    print("=" * 60)
    print("Testing Force Constants Parser")
    print("=" * 60)
    
    try:
        fc = ForceConstants(fc_file)
        info = fc.get_structure_info()
        
        print("\n" + "=" * 60)
        print("Structure Information")
        print("=" * 60)
        print(f"Lattice: {info['a']:.3f} × {info['b']:.3f} × {info['c']:.3f} Å")
        print(f"Species: {', '.join(info['species'])}")
        print(f"Atoms in primitive cell: {len(info['species'])}")
        
    except FileNotFoundError:
        print(f"\nFile not found. Please save your .fc file to:")
        print(f"  data/AuTe_2_m.fc")
        print("\nCreate data directory:")
        print("  mkdir -p data/")
