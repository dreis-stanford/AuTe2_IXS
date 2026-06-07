"""
Handle modulated/superstructure reflections for AuTe2
Deals with satellite reflections from CDW modulation
"""

import numpy as np
from typing import List, Tuple, Dict, Optional

from . import config

class ModulatedStructure:
    """
    Handle satellite reflections for modulated AuTe2 structure
    """
    
    def __init__(self, 
                 q_mod: Tuple[float, float, float] = None,
                 lattice_params: Dict[str, float] = None):
        """
        Parameters:
        -----------
        q_mod : Modulation vector in r.l.u. (default from config)
        lattice_params : Lattice parameters (default from config)
        """
        self.q_mod = np.array(q_mod if q_mod else config.get_modulation_vector())
        self.lattice = lattice_params if lattice_params else config.get_lattice_parameters()
    
    def satellite_hkl(self, 
                      h: int, k: int, l: int, 
                      m: int) -> Tuple[float, float, float]:
        """
        Calculate m-th order satellite position
        
        Q_satellite = Q_main + m * q_mod
        
        Parameters:
        -----------
        h, k, l : Main reflection Miller indices
        m : Satellite order (±1, ±2, ...)
        
        Returns:
        --------
        (H, K, L) position in r.l.u.
        
        Example:
        --------
        >>> mod = ModulatedStructure(q_mod=(0.28, 0, 0))
        >>> mod.satellite_hkl(2, 0, 0, m=1)
        (2.28, 0.0, 0.0)
        >>> mod.satellite_hkl(2, 0, 0, m=-1)
        (1.72, 0.0, 0.0)
        """
        hkl_main = np.array([h, k, l])
        hkl_sat = hkl_main + m * self.q_mod
        return tuple(hkl_sat)
    
    def generate_reflection_list(self,
                                 h_range: Tuple[int, int] = (-3, 4),
                                 k_range: Tuple[int, int] = (-3, 4),
                                 l_range: Tuple[int, int] = (-3, 4),
                                 max_satellite_order: int = 2,
                                 max_q: float = None) -> List[Dict]:
        """
        Generate list of main and satellite reflections
        
        Parameters:
        -----------
        h_range, k_range, l_range : Range of main reflections
        max_satellite_order : Maximum |m| for satellites
        max_q : Maximum |Q| in nm^-1 (default from config)
        
        Returns:
        --------
        List of dicts with keys:
            'hkl' : (H, K, L) position
            'type' : 'main' or 'satellite'
            'order' : satellite order (0 for main)
            'parent' : parent (h,k,l) for satellites
            'q_magnitude' : |Q| in nm^-1
        """
        if max_q is None:
            max_q = config.MAX_Q
        
        reflections = []
        
        for h in range(*h_range):
            for k in range(*k_range):
                for l in range(*l_range):
                    if h == 0 and k == 0 and l == 0:
                        continue
                    
                    # Main reflection
                    hkl_main = (h, k, l)
                    q_mag_main = self.calculate_q_magnitude(hkl_main)
                    
                    if q_mag_main <= max_q:
                        reflections.append({
                            'hkl': hkl_main,
                            'type': 'main',
                            'order': 0,
                            'parent': None,
                            'q_magnitude': q_mag_main
                        })
                    
                    # Satellites
                    for m in range(-max_satellite_order, max_satellite_order + 1):
                        if m == 0:
                            continue
                        
                        hkl_sat = self.satellite_hkl(h, k, l, m)
                        q_mag_sat = self.calculate_q_magnitude(hkl_sat)
                        
                        if q_mag_sat <= max_q:
                            reflections.append({
                                'hkl': hkl_sat,
                                'type': 'satellite',
                                'order': m,
                                'parent': hkl_main,
                                'q_magnitude': q_mag_sat
                            })
        
        return reflections
    
    def calculate_q_magnitude(self, hkl: Tuple[float, float, float]) -> float:
        """
        Calculate |Q| for given HKL using metric tensor
        
        Returns |Q| in nm^-1
        """
        h, k, l = hkl
        
        a = self.lattice['a']
        b = self.lattice['b']
        c = self.lattice['c']
        alpha = np.radians(self.lattice['alpha'])
        beta = np.radians(self.lattice['beta'])
        gamma = np.radians(self.lattice['gamma'])
        
        # Reciprocal lattice parameters
        V = a * b * c * np.sqrt(1 + 2*np.cos(alpha)*np.cos(beta)*np.cos(gamma)
                                - np.cos(alpha)**2 - np.cos(beta)**2 - np.cos(gamma)**2)
        
        a_star = 2 * np.pi * b * c * np.sin(alpha) / V
        b_star = 2 * np.pi * a * c * np.sin(beta) / V
        c_star = 2 * np.pi * a * b * np.sin(gamma) / V
        
        # Metric tensor in reciprocal space
        cos_alpha_star = (np.cos(beta)*np.cos(gamma) - np.cos(alpha)) / (np.sin(beta)*np.sin(gamma))
        cos_beta_star = (np.cos(alpha)*np.cos(gamma) - np.cos(beta)) / (np.sin(alpha)*np.sin(gamma))
        cos_gamma_star = (np.cos(alpha)*np.cos(beta) - np.cos(gamma)) / (np.sin(alpha)*np.sin(beta))
        
        # Calculate |Q|^2
        q_squared = (h * a_star)**2 + (k * b_star)**2 + (l * c_star)**2 + \
                    2*h*k * a_star * b_star * cos_gamma_star + \
                    2*h*l * a_star * c_star * cos_beta_star + \
                    2*k*l * b_star * c_star * cos_alpha_star
        
        # Return in nm^-1
        return np.sqrt(q_squared) * 10.0
    
    def find_satellites_near_q(self,
                               q_target: Tuple[float, float, float],
                               tolerance: float = 0.1) -> List[Dict]:
        """
        Find satellite reflections near a target Q
        
        Useful for identifying which satellites contribute to phonon at given Q
        
        Parameters:
        -----------
        q_target : Target Q in r.l.u.
        tolerance : Distance tolerance in r.l.u.
        
        Returns:
        --------
        List of nearby reflections sorted by distance
        """
        reflections = self.generate_reflection_list()
        q_target = np.array(q_target)
        
        nearby = []
        for refl in reflections:
            q_refl = np.array(refl['hkl'])
            distance = np.linalg.norm(q_refl - q_target)
            
            if distance <= tolerance:
                refl['distance'] = distance
                nearby.append(refl)
        
        # Sort by distance
        nearby.sort(key=lambda x: x['distance'])
        
        return nearby
    
    def bragg_accessible_satellites(self,
                                   wavelength: float = 0.570107,
                                   max_2theta: float = 150.0) -> List[Dict]:
        """
        Find Bragg-accessible satellite reflections
        
        Parameters:
        -----------
        wavelength : X-ray wavelength in Angstroms
        max_2theta : Maximum scattering angle in degrees
        
        Returns:
        --------
        List of accessible reflections with scattering angles
        """
        reflections = self.generate_reflection_list()
        accessible = []
        
        for refl in reflections:
            q_mag_ang = refl['q_magnitude'] / 10.0  # nm^-1 to Angstrom^-1
            
            # Bragg condition: sin(theta) = lambda * |Q| / (4*pi)
            sin_theta = wavelength * q_mag_ang / (4 * np.pi)
            
            if 0 < sin_theta <= 1.0:
                theta = np.degrees(np.arcsin(sin_theta))
                two_theta = 2 * theta
                
                if two_theta <= max_2theta:
                    refl['theta'] = theta
                    refl['two_theta'] = two_theta
                    accessible.append(refl)
        
        # Sort by 2theta
        accessible.sort(key=lambda x: x['two_theta'])
        
        return accessible
    
    def print_satellite_info(self, max_order: int = 2):
        """Print information about satellite reflections"""
        print("\n" + "="*70)
        print(f"Modulated Structure - {config.SAMPLE_NAME}")
        print("="*70)
        print(f"Modulation vector: q_mod = {tuple(self.q_mod)}")
        print(f"  |q_mod| = {np.linalg.norm(self.q_mod):.4f} r.l.u.")
        print(f"  = {self.calculate_q_magnitude(tuple(self.q_mod)):.3f} nm^-1")
        
        print(f"\nSatellite positions for (2,0,0) main reflection:")
        for m in range(-max_order, max_order + 1):
            hkl_sat = self.satellite_hkl(2, 0, 0, m)
            q_mag = self.calculate_q_magnitude(hkl_sat)
            sat_type = "main" if m == 0 else f"m={m:+d} satellite"
            print(f"  {sat_type:15s}: ({hkl_sat[0]:6.3f}, {hkl_sat[1]:6.3f}, {hkl_sat[2]:6.3f})  "
                  f"|Q| = {q_mag:6.2f} nm^-1")
        
        print("="*70 + "\n")
