"""
Physical constants for phonon and IXS calculations
Converted from loadconstants.m
"""

import numpy as np

class PhysicalConstants:
    """Physical constants in convenient units for IXS calculations"""
    
    def __init__(self):
        # Fundamental constants
        self.hbarc = 1973.269718  # eV·Å
        self.hc = self.hbarc * 2 * np.pi  # eV·Å
        self.alpha = 7.2973525643e-3  # Fine structure constant
        
        # Energy conversions
        self.THz2eV = self.hc / 3e6  # 0.00414 eV/THz
        self.eV2THz = 1.0 / self.THz2eV
        self.meV2THz = self.eV2THz / 1000
        self.THz2meV = 1000 * self.THz2eV
        
        # Particle masses
        self.Me = 510998.95069  # Electron mass in eV
        self.M_u = 938.27208816 / 1.007276466621 * 1e6  # Atomic mass unit in eV
        self.eV2kg = 9.1093837139e-31 / self.Me
        
        # Charge and speed of light
        self.q_e = 1.602176634e-19  # Electron charge (C)
        self.c = 299792458  # Speed of light (m/s)
        
        # Unit conversions
        self.Ninvm2THz = self.c**2 / self.q_e * 1e-24 / self.THz2eV
        
        # Lengths
        self.re = self.hbarc / self.Me * self.alpha  # Classical electron radius (Å)
        self.Bohr = self.re / self.alpha**2  # Bohr radius (Å)
        
        # Energy scales
        self.eVtoHartree = 1 / (self.Me * self.alpha**2)  # 1/27.2114
        
        # Boltzmann constant (for reference)
        self.kB = 8.617333262e-5  # eV/K
        
        # Wavenumber conversions
        self.cm_inv_to_THz = self.c * 100 / 1e12  # cm^-1 to THz
        self.THz_to_cm_inv = 1 / self.cm_inv_to_THz
        self.cm_inv_to_meV = self.cm_inv_to_THz * self.THz2meV
        self.meV_to_cm_inv = 1 / self.cm_inv_to_meV
        
    def __repr__(self):
        return f"PhysicalConstants(ℏc={self.hbarc:.2f} eV·Å, α={self.alpha:.6f})"

# Create global instance
const = PhysicalConstants()

if __name__ == "__main__":
    print("Physical Constants for IXS Calculations")
    print("=" * 50)
    print(f"ℏc = {const.hbarc:.4f} eV·Å")
    print(f"hc = {const.hc:.4f} eV·Å")
    print(f"α (fine structure) = {const.alpha:.10f}")
    print(f"\nEnergy conversions:")
    print(f"  1 THz = {const.THz2eV:.6f} eV = {const.THz2meV:.4f} meV")
    print(f"  1 cm⁻¹ = {const.cm_inv_to_THz:.6f} THz = {const.cm_inv_to_meV:.4f} meV")
    print(f"\nMasses:")
    print(f"  Electron: {const.Me:.2f} eV/c²")
    print(f"  Atomic mass unit: {const.M_u:.2e} eV/c²")
    print(f"\nLengths:")
    print(f"  Bohr radius: {const.Bohr:.6f} Å")
    print(f"  Classical e⁻ radius: {const.re:.6e} Å")
    print(f"\nSpeed of light: {const.c:.0f} m/s")
