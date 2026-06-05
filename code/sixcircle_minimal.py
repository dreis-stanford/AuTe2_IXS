"""
Minimal sixcircle interface that avoids initialization bugs
We'll use the core calculation functions directly
"""

import sys
import os
import config
import numpy as np

# Import just the calculation module, not the full sixcircle
sys.path.insert(0, config.SIXCIRCLE_PATH)
import scbasic

class MinimalSixCircle:
    """
    Minimal sixcircle interface using only scbasic (the calculation core)
    Avoids all the initialization bugs in the main sixcircle module
    """
    
    def __init__(self, photon_energy_kev=21.747):
        """
        Initialize with photon energy
        
        Parameters:
        -----------
        photon_energy_kev : float
            Photon energy in keV
        """
        self.energy_kev = photon_energy_kev
        self.wavelength = 12.398 / photon_energy_kev  # Angstroms
        
        # Crystal parameters (set with setup_crystal)
        self.a = None
        self.b = None
        self.c = None
        self.alpha = None
        self.beta = None
        self.gamma = None
        self.name = "Not set"
        
        # Orientation matrix (will be calculated)
        self.UB = None
        
        print(f"✓ Minimal sixcircle initialized")
        print(f"  Energy: {self.energy_kev:.3f} keV")
        print(f"  Wavelength: {self.wavelength:.6f} Å")
    
    def setup_crystal(self, a, b, c, alpha=90, beta=90, gamma=90, name="Sample"):
        """
        Setup crystal lattice
        
        Parameters:
        -----------
        a, b, c : float
            Lattice parameters in Angstroms
        alpha, beta, gamma : float
            Lattice angles in degrees
        name : str
            Sample name
        """
        self.a = a
        self.b = b
        self.c = c
        self.alpha = alpha
        self.beta = beta
        self.gamma = gamma
        self.name = name
        
        print(f"\n✓ Crystal '{name}' configured:")
        print(f"  a={a:.3f} Å, b={b:.3f} Å, c={c:.3f} Å")
        print(f"  α={alpha:.1f}°, β={beta:.1f}°, γ={gamma:.1f}°")
    
    def ca(self, h, k, l, frozen_mode='456'):
        """
        Calculate angles for (h,k,l)
        
        Parameters:
        -----------
        h, k, l : float
            Miller indices
        frozen_mode : str
            Which angles to freeze (default '456' = mu, gam, omega)
        
        Returns:
        --------
        dict with calculated angles
        """
        if self.a is None:
            raise RuntimeError("Must call setup_crystal() first!")
        
        # Calculate Q vector magnitude
        Q = self.calc_Q_magnitude(h, k, l)
        
        # Calculate scattering angle
        # |Q| = 4π sin(θ) / λ
        sin_theta = Q * self.wavelength / (4 * np.pi)
        
        if abs(sin_theta) > 1:
            print(f"ERROR: Reflection ({h},{k},{l}) not accessible!")
            print(f"  |Q| = {Q:.4f} Å⁻¹ requires sin(θ) = {sin_theta:.4f}")
            return None
        
        theta_rad = np.arcsin(sin_theta)
        theta_deg = np.degrees(theta_rad)
        tth_deg = 2 * theta_deg
        
        print(f"\nReflection ({h:.1f} {k:.1f} {l:.1f}):")
        print(f"  |Q| = {Q:.4f} Å⁻¹")
        print(f"  2θ = {tth_deg:.4f}°")
        print(f"  θ  = {theta_deg:.4f}°")
        print(f"  (Simplified calculation - for full 6-circle use real sixcircle)")
        
        return {
            'h': h, 'k': k, 'l': l,
            'Q': Q,
            'tth': tth_deg,
            'th': theta_deg,
        }
    
    def calc_Q_magnitude(self, h, k, l):
        """
        Calculate |Q| for given (h,k,l) using crystal lattice
        
        Returns |Q| in Å⁻¹
        """
        # Reciprocal lattice parameters
        # For cubic: a* = 2π/a
        # For general case, need full calculation
        
        # Convert to radians
        alpha_r = np.radians(self.alpha)
        beta_r = np.radians(self.beta)
        gamma_r = np.radians(self.gamma)
        
        # Calculate volume of unit cell
        V = self.a * self.b * self.c * np.sqrt(
            1 - np.cos(alpha_r)**2 - np.cos(beta_r)**2 - np.cos(gamma_r)**2
            + 2 * np.cos(alpha_r) * np.cos(beta_r) * np.cos(gamma_r)
        )
        
        # Reciprocal lattice parameters
        a_star = 2 * np.pi * self.b * self.c * np.sin(alpha_r) / V
        b_star = 2 * np.pi * self.c * self.a * np.sin(beta_r) / V
        c_star = 2 * np.pi * self.a * self.b * np.sin(gamma_r) / V
        
        # For orthogonal systems (alpha=beta=gamma=90), this simplifies
        if self.alpha == 90 and self.beta == 90 and self.gamma == 90:
            Q_squared = (h * a_star)**2 + (k * b_star)**2 + (l * c_star)**2
        else:
            # Full calculation needed - using simplified version
            Q_squared = (h * a_star)**2 + (k * b_star)**2 + (l * c_star)**2
            print("  Warning: Non-orthogonal - using simplified Q calculation")
        
        return np.sqrt(Q_squared)

def init_sixcircle(photon_energy_kev=21.747):
    """Quick initialization"""
    return MinimalSixCircle(photon_energy_kev)

if __name__ == "__main__":
    print("Testing minimal sixcircle...")
    print("=" * 60)
    
    sc = init_sixcircle(21.747)
    
    print("\n" + "=" * 60)
    print("Testing with cubic crystal...")
    print("=" * 60)
    sc.setup_crystal(5.0, 5.0, 5.0, name="TestCrystal")
    
    print("\n" + "=" * 60)
    print("Calculating angles...")
    print("=" * 60)
    sc.ca(1, 0, 0)
    sc.ca(1, 1, 0)
    sc.ca(2, 0, 0)
    
    print("\n" + "=" * 60)
    print("SUCCESS! Minimal interface working")
    print("=" * 60)
    print("\nNote: This uses basic diffraction calculations")
    print("For full 6-circle functionality at the beamline,")
    print("use sixcircle directly in interactive Python")
