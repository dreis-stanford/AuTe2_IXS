"""
Physical constants for phonon calculations
Converted from loadconstants.m
"""

import numpy as np

class Constants:
    """Physical constants in convenient units"""
    
    # Fundamental constants
    hbarc = 1.973269718e3  # eV·Å
    hc = hbarc * 2 * np.pi  # eV·Å
    alpha = 7.2973525643e-3  # Fine structure constant
    
    # Energy/frequency conversions
    THz2eV = hc / 3e6  # ~0.0041 eV/THz
    eV2THz = 1.0 / THz2eV  # THz/eV
    
    # Mass conversions
    Me = 510998.95069  # Electron mass in eV/c^2
    M_u = 938.27208816 / 1.007276466621 * 1e6  # Atomic mass unit in eV/c^2 (~931.49 MeV)
    eV2kg = 9.1093837139e-31 / Me  # Convert eV to kg
    
    # SI constants
    q_e = 1.602176634e-19  # Electron charge in Coulombs
    c = 299792458  # Speed of light in m/s
    
    # Derived constants
    Ninvm2THz = c**2 / q_e * 1e-24 / THz2eV  # Convert N/m to THz
    
    # Atomic units
    re = hbarc / Me * alpha  # Classical electron radius
    Bohr = re / alpha**2  # Bohr radius in Angstroms (~0.529 Å)
    eVtoHartree = 1 / (Me * alpha**2)  # ~1/27.2114
    HartreeToEV = Me * alpha**2  # Inverse
    
    # Additional useful conversions
    Ry = Me * alpha**2 / 2  # Rydberg in eV (~13.6 eV)
    Ha = Me * alpha**2  # Hartree in eV (~27.2 eV)
    
    # Frequency unit conversions
    meV2THz = 1000.0 / THz2eV  # meV to THz
    THz2meV = THz2eV / 1000.0   # THz to meV
    invCm2meV = 0.12398419843   # cm^-1 to meV
    meV2invCm = 1.0 / invCm2meV # meV to cm^-1
    
    # Thermal
    kb = 8.617333262e-5  # Boltzmann constant in eV/K
    hbar = 1.054571817e-34  # Reduced Planck constant in J·s
    amu = 1.66053906660e-27  # Atomic mass unit in kg

# Create singleton instance
const = Constants()

# Print some key values for verification
if __name__ == '__main__':
    print("Physical Constants Loaded:")
    print(f"  ℏc = {const.hbarc:.6f} eV·Å")
    print(f"  Bohr radius = {const.Bohr:.9f} Å")
    print(f"  THz to eV = {const.THz2eV:.6f}")
    print(f"  THz to meV = {const.THz2meV:.6f}")
    print(f"  Hartree to eV = {const.HartreeToEV:.6f}")
    print(f"  Rydberg = {const.Ry:.6f} eV")
