"""
Force constants parser for phonon calculations
Converted from readFC.m - reads .fc files from DFT calculations
"""

import numpy as np
from .constants import const
import re

class ForceConstants:
    """
    Parse and process force constants file
    Equivalent to MATLAB readFC.m
    """
    
    def __init__(self, filename):
        """
        Load force constants from .fc file
        
        Parameters:
        -----------
        filename : str
            Path to force constants file
        """
        self.filename = filename
        self.parse_file()
        
    def parse_file(self):
        """Parse the .fc file - matches readFC.m exactly"""
        
        with open(self.filename, 'r') as f:
            lines = [line.rstrip() for line in f.readlines()]
        
        idx = 0
        
        # 1. Header & Lattice
        line1 = [float(x) for x in lines[idx].split()]
        self.ntyp = int(line1[0])  # Number of atom types
        self.nat = int(line1[1])   # Number of atoms in primitive cell
        ascale = line1[3] * const.Bohr  # Scale factor: Bohr to Angstrom
        idx += 1
        
        print(f"Reading force constants file: {self.filename}")
        print(f"  Number of atom types: {self.ntyp}")
        print(f"  Atoms in primitive cell: {self.nat}")
        print(f"  Lattice scale factor: {ascale:.6f} Å")
        
        # Direct lattice vectors (Cartesian, in Angstrom)
        self.a_l = np.zeros((3, 3))
        for i in range(3):
            self.a_l[i, :] = np.array([float(x) for x in lines[idx].split()]) * ascale
            idx += 1
        
        self.aunits = 'Angstrom'
        
        # Lattice parameters from vectors
        self.a = np.linalg.norm(self.a_l[0, :])
        self.b = np.linalg.norm(self.a_l[1, :])
        self.c = np.linalg.norm(self.a_l[2, :])
        
        # Angles between lattice vectors
        self.alpha = np.degrees(np.arccos(
            np.dot(self.a_l[1, :], self.a_l[2, :]) / (self.b * self.c)))
        self.beta = np.degrees(np.arccos(
            np.dot(self.a_l[2, :], self.a_l[0, :]) / (self.c * self.a)))
        self.gamma = np.degrees(np.arccos(
            np.dot(self.a_l[0, :], self.a_l[1, :]) / (self.a * self.b)))
        
        # Unit cell volume
        self.Vcell = np.dot(self.a_l[0, :], 
                           np.cross(self.a_l[1, :], self.a_l[2, :]))
        
        # Reciprocal lattice vectors (Cartesian, 2π/Angstrom)
        self.b_l = np.zeros((3, 3))
        self.b_l[0, :] = 2*np.pi/self.Vcell * np.cross(self.a_l[1, :], self.a_l[2, :])
        self.b_l[1, :] = 2*np.pi/self.Vcell * np.cross(self.a_l[2, :], self.a_l[0, :])
        self.b_l[2, :] = 2*np.pi/self.Vcell * np.cross(self.a_l[0, :], self.a_l[1, :])
        self.bunits = '2pi/Angstrom'
        
        print(f"\nLattice parameters:")
        print(f"  a = {self.a:.4f} Å")
        print(f"  b = {self.b:.4f} Å")
        print(f"  c = {self.c:.4f} Å")
        print(f"  α = {self.alpha:.2f}°")
        print(f"  β = {self.beta:.2f}°")
        print(f"  γ = {self.gamma:.2f}°")
        print(f"  Volume = {self.Vcell:.4f} Å³")
        
        # 2. Species and Masses
        # Format: index 'symbol' mass_value
        # Use regex to handle varying whitespace
        self.masses = np.zeros(self.ntyp)
        self.symbols = []
        
        for i in range(self.ntyp):
            line = lines[idx]
            # Match pattern: number 'text' number
            match = re.search(r"(\d+)\s+'([^']+)'\s+([\d\.eE\+\-]+)", line)
            if match:
                type_idx = int(match.group(1)) - 1  # Convert to 0-based
                symbol = match.group(2).strip()
                m_val = float(match.group(3))
                
                # Convert from file units to amu
                self.masses[type_idx] = m_val / const.M_u * 1e6  # Now in amu
                self.symbols.append(symbol)
            else:
                raise ValueError(f"Could not parse mass line: {line}")
            idx += 1
        
        self.massunit = 'amu'
        
        print(f"\nAtomic species:")
        for i in range(self.ntyp):
            print(f"  {i+1}. {self.symbols[i]:3s}: {self.masses[i]:7.3f} amu")
        
        # 3. Atomic Positions
        self.rs = np.zeros((self.nat, 3))  # Cartesian positions
        self.atom_type_map = np.zeros(self.nat, dtype=int)
        
        for i in range(self.nat):
            row = [float(x) for x in lines[idx].split()]
            self.atom_type_map[i] = int(row[1])
            self.rs[i, :] = np.array(row[2:5]) * ascale
            idx += 1
        
        # Fractional coordinates
        self.xs = self.rs @ np.linalg.inv(self.a_l)
        
        print(f"\nAtomic positions (fractional):")
        for i in range(self.nat):
            atom_type = self.atom_type_map[i]
            symbol = self.symbols[atom_type - 1]
            print(f"  Atom {i+1} ({symbol}): [{self.xs[i,0]:7.4f}, {self.xs[i,1]:7.4f}, {self.xs[i,2]:7.4f}]")
        
        # 4. Dielectric flag
        has_dielectric = lines[idx].strip()
        idx += 1
        # If 'T', would read epsilon_infinity here (not implemented)
        
        # 5. Grid and Force Constants
        # Skip any lines that don't look like a grid (e.g., supercell vectors)
        while idx < len(lines):
            line_parts = lines[idx].split()
            try:
                # Try to parse as integers
                grid_line = [int(x) for x in line_parts]
                # Check if it looks like a grid (should be 3 integers)
                if len(grid_line) == 3 and all(g > 0 for g in grid_line):
                    self.grid = np.array(grid_line)
                    break
            except ValueError:
                # This line has floats or other non-integer data, skip it
                pass
            idx += 1
        else:
            raise ValueError("Could not find supercell grid in force constants file")
        n_cells = int(np.prod(self.grid))
        idx += 1
        
        print(f"\nForce constant supercell grid: {self.grid}")
        print(f"  Total supercell images: {n_cells}")
        
        # Initialize force constants tensor
        # phi[3*nat, 3*nat, n_cells]
        self.phi = np.zeros((3*self.nat, 3*self.nat, n_cells))
        self.uvw = np.zeros((n_cells, 3), dtype=int)
        
        first_in_loop = True
        
        print(f"\nReading force constants tensor...")
        
        # Read force constants
        # Loop order: mu, nu, r, s (Cartesian indices and atoms)
        for mu in range(1, 4):  # Cartesian direction 1
            for nu in range(1, 4):  # Cartesian direction 2
                for r in range(1, self.nat + 1):  # Atom 1
                    for s in range(1, self.nat + 1):  # Atom 2
                        # Skip index line
                        idx += 1
                        
                        # Read force constants for all supercell translations
                        for w in range(1, self.grid[2] + 1):
                            for v in range(1, self.grid[1] + 1):
                                for u in range(1, self.grid[0] + 1):
                                    # Cell index (0-based for Python)
                                    k = (u-1) + self.grid[0]*(v-1) + \
                                        self.grid[0]*self.grid[1]*(w-1)
                                    
                                    if first_in_loop:
                                        self.uvw[k, :] = [u, v, w]
                                    
                                    data_line = [float(x) for x in lines[idx].split()]
                                    idx += 1
                                    
                                    # Store force constant
                                    # Convert to 0-based indexing
                                    i_idx = (mu-1) + 3*(r-1)
                                    j_idx = (nu-1) + 3*(s-1)
                                    self.phi[i_idx, j_idx, k] = data_line[3]
                        
                        first_in_loop = False
        
        print(f"  ✓ Force constants loaded")
        
        # Shift uvw to center around origin
        # From 1-based [1,2,3,4] to 0-based [0,1,2,3] to centered [-2,-1,0,1]
        self.uvw = self.uvw - 1
        self.uvw = np.mod(self.uvw + 2, 4) - 2
        
        print(f"  ✓ Supercell translations centered around origin")
        
        # Check symmetry of phi at origin
        diff = self.phi[:, :, 0] - self.phi[:, :, 0].T
        max_asymmetry = np.max(np.abs(diff))
        if max_asymmetry > 1e-6:
            print(f"  ⚠ Warning: Force constant matrix at origin not symmetric")
            print(f"    Max asymmetry: {max_asymmetry:.3e}")
        
        # Enforce acoustic sum rule
        print(f"  Enforcing acoustic sum rule...")
        for d in range(3 * self.nat):
            current_sum = np.sum(self.phi[d, :, :])
            self.phi[d, d, 0] -= current_sum
        
        # Verify sum rule
        row_sums = np.sum(self.phi, axis=(1, 2))
        max_sum = np.max(np.abs(row_sums))
        print(f"  ✓ Acoustic sum rule enforced (max residual: {max_sum:.3e})")
        
        print(f"\n{'='*60}")
        print(f"Force constants loaded successfully!")
        print(f"{'='*60}\n")
    
    def convert_to_eV_per_Angstrom2(self):
        """
        Convert force constants to eV/Ų units
        Matches MATLAB: Phi = xtal.phi * 1/2 * const.alpha^2 * const.Me / const.Bohr^2
        
        Returns:
        --------
        Phi : ndarray
            Force constants in eV/Ų
        """
        conversion = 0.5 * const.alpha**2 * const.Me / const.Bohr**2
        Phi = self.phi * conversion
        
        print(f"Force constants converted to eV/Ų")
        print(f"  Conversion factor: {conversion:.6e}")
        
        return Phi
    
    def get_structure_info(self):
        """Return structure information as dictionary"""
        return {
            'a': self.a,
            'b': self.b,
            'c': self.c,
            'alpha': self.alpha,
            'beta': self.beta,
            'gamma': self.gamma,
            'Vcell': self.Vcell,
            'a_l': self.a_l,
            'b_l': self.b_l,
            'symbols': self.symbols,
            'masses': self.masses,
            'positions_cartesian': self.rs,
            'positions_fractional': self.xs,
            'atom_type_map': self.atom_type_map,
            'grid': self.grid,
            'uvw': self.uvw
        }
    
    def __repr__(self):
        return (f"ForceConstants('{self.filename}': "
                f"{self.nat} atoms, {' '.join(self.symbols)}, "
                f"grid={self.grid})")


if __name__ == "__main__":
    import os
    import sys
    
    # Add parent directory to path for imports
    sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
    
    # Test the parser
    fc_file = "data/AuTe_2_m.fc"
    
    if not os.path.exists(fc_file):
        print(f"File not found: {fc_file}")
        print(f"Current directory: {os.getcwd()}")
        print(f"Looking for: {os.path.abspath(fc_file)}")
    else:
        print("=" * 60)
        print("Testing Force Constants Parser")
        print("=" * 60)
        print()
        
        # Load force constants
        xtal = ForceConstants(fc_file)
        
        # Convert to eV/Angstrom^2
        Phi = xtal.convert_to_eV_per_Angstrom2()
        
        # Test results
        print(f"\nTest Results:")
        print(f"  Structure: {xtal}")
        print(f"  Force constants shape: {Phi.shape}")
        print(f"  Force constant range: [{np.min(Phi):.3e}, {np.max(Phi):.3e}] eV/Ų")
        
        # Get structure info
        info = xtal.get_structure_info()
        print(f"\nStructure info dictionary has {len(info)} fields")
