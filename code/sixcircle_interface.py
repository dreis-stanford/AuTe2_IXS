"""
Interface to sixcircle diffractometer control for AuTe2 IXS experiments
Integrates with Alfred Baron's sixcircle code from SPring-8
"""

import numpy as np
from typing import List, Tuple, Dict, Optional
import sys

from . import config

class SixCircleInterface:
    """
    Interface to sixcircle package for diffractometer control
    """
    
    def __init__(self, simulation_mode: bool = None):
        """
        Initialize sixcircle interface
        
        Parameters:
        -----------
        simulation_mode : Force simulation mode (default: auto-detect)
        """
        if simulation_mode is None:
            simulation_mode = not config.SIXCIRCLE_AVAILABLE
        
        self.simulation_mode = simulation_mode
        self.sixc = None
        self.sixc_rqd = None
        
        if not simulation_mode and config.SIXCIRCLE_PATH:
            self._load_sixcircle()
    
    def _load_sixcircle(self):
        """Load sixcircle modules"""
        try:
            sys.path.insert(0, config.SIXCIRCLE_PATH)
            import sixcircle
            import sixcircle_rqd
            
            self.sixc = sixcircle
            self.sixc_rqd = sixcircle_rqd
            
            print("✓ Sixcircle modules loaded successfully")
            self.simulation_mode = False
            
        except ImportError as e:
            print(f"✗ Could not load sixcircle: {e}")
            print("  Running in simulation mode")
            self.simulation_mode = True
    
    def setup_experiment(self):
        """
        Setup sixcircle for AuTe2 experiment using config parameters
        """
        lattice = config.get_lattice_parameters()
        
        if self.simulation_mode:
            print("\n" + "="*60)
            print("SIMULATION MODE - Sixcircle Setup")
            print("="*60)
            print(f"Sample: {config.SAMPLE_NAME}")
            print(f"Lattice: a={lattice['a']:.4f} b={lattice['b']:.4f} c={lattice['c']:.4f}")
            print(f"         α={lattice['alpha']:.1f}° β={lattice['beta']:.1f}° γ={lattice['gamma']:.1f}°")
            print(f"Si order: {config.SI_ORDER}")
            print(f"Beamline: BL{config.BEAMLINE}LXU")
            print("="*60 + "\n")
            return
        
        # Real sixcircle setup
        print("\nConfiguring sixcircle...")
        
        # Set lattice
        self.sixc.setlat(
            lattice['a'], lattice['b'], lattice['c'],
            lattice['alpha'], lattice['beta'], lattice['gamma']
        )
        
        # Set wavelength via Si order
        self.sixc_rqd.setorder(config.SI_ORDER)
        
        # Set beamline
        self.sixc_rqd.setbl(config.BEAMLINE)
        
        # Set incident beam
        self.sixc_rqd.setincident(config.INCIDENT_TYPE)
        
        # Set frozen angles
        self.sixc.setfrozen(config.FROZEN_ANGLES)
        
        print(f"✓ Sixcircle configured for {config.SAMPLE_NAME}")
    
    def set_orientation(self, 
                       or0_hkl: Tuple[float, float, float],
                       or0_angles: Dict[str, float],
                       or1_hkl: Tuple[float, float, float],
                       or1_angles: Dict[str, float]):
        """
        Set orientation matrix from two reflections
        
        Parameters:
        -----------
        or0_hkl : Primary reflection (H, K, L)
        or0_angles : Angles where or0 was found {'tth': ?, 'th': ?, 'chi': ?, 'phi': ?}
        or1_hkl : Secondary reflection (H, K, L)
        or1_angles : Angles where or1 was found
        """
        if self.simulation_mode:
            print(f"\nSIMULATION: Setting orientation")
            print(f"  or0: {or0_hkl} at {or0_angles}")
            print(f"  or1: {or1_hkl} at {or1_angles}")
            return
        
        # Move to or0 and set it
        self.sixc.mv(**or0_angles)
        self.sixc.or0(*or0_hkl)
        
        # Move to or1 and set it
        self.sixc.mv(**or1_angles)
        self.sixc.or1(*or1_hkl)
        
        # Check consistency
        print("\nOrientation matrix set. Checking consistency...")
        self.sixc.or_check()
    
    def calculate_ixs_positions(self, 
                               hkl: Tuple[float, float, float],
                               agap_v: float = None,
                               agap_h: float = None) -> Dict:
        """
        Calculate IXS analyzer positions for measurement at given Q
        """
        if agap_v is None:
            agap_v = config.DEFAULT_AGAP_V
        if agap_h is None:
            agap_h = config.DEFAULT_AGAP_H
        
        if self.simulation_mode:
            return self._simulate_ixs_positions(hkl, agap_v, agap_h)
        
        # Set analyzer gaps
        self.sixc_rqd.mvgap(agapv=agap_v, agaph=agap_h)
        
        # Calculate analyzer positions
        print(f"\nCalculating IXS geometry for Q = {hkl}")
        print(f"Analyzer gaps: V={agap_v}μm, H={agap_h}μm")
        
        self.sixc_rqd.ca6(*hkl)
        
        return {'center_hkl': hkl, 'status': 'calculated'}
    
    def _simulate_ixs_positions(self, hkl, agap_v, agap_h):
        """Simulate IXS analyzer positions"""
        print(f"\nSIMULATION: IXS positions for Q = {hkl}")
        print(f"  Analyzer gaps: V={agap_v}μm, H={agap_h}μm")
        
        # Simulate some analyzer positions
        q_mag = np.linalg.norm(hkl)
        dq = 0.04  # Typical Q resolution
        
        analyzers = {}
        for i in range(-3, 4):
            h_offset = i * dq
            analyzers[f'a{i+10:02d}'] = {
                'hkl': tuple(np.array(hkl) + np.array([h_offset, 0, 0])),
                'q_mag': q_mag + h_offset,
                'dq': dq
            }
        
        return {
            'center_hkl': hkl,
            'analyzers': analyzers,
            'gaps': {'v': agap_v, 'h': agap_h}
        }
    
    def move_to_hkl(self, hkl: Tuple[float, float, float], 
                    check_only: bool = False):
        """
        Move to (or calculate) angles for given HKL
        
        Parameters:
        -----------
        hkl : Target (H, K, L)
        check_only : If True, only calculate without moving
        """
        if self.simulation_mode:
            angles = self._simulate_angles(hkl)
            print(f"\nSIMULATION: {'Would move' if not check_only else 'Angles'} to {hkl}")
            for key, val in angles.items():
                print(f"  {key:5s} = {val:7.3f}°")
            return angles
        
        if check_only:
            self.sixc.ca(*hkl)
        else:
            self.sixc.br(*hkl)
        
        return self.get_current_angles()
    
    def get_current_angles(self) -> Dict[str, float]:
        """Get current diffractometer angles"""
        if self.simulation_mode:
            return {'tth': 0, 'th': 0, 'chi': 0, 'phi': 0, 'mu': -0.17, 'gam': 0}
        
        # Would extract from sixcircle state
        # Implementation depends on sixcircle version
        return {}
    
    def _simulate_angles(self, hkl):
        """Simulate angle calculation for testing"""
        lattice = config.get_lattice_parameters()
        
        # Simple calculation assuming orthorhombic
        a_star = 2 * np.pi / lattice['a']
        b_star = 2 * np.pi / lattice['b']
        c_star = 2 * np.pi / lattice['c']
        
        h, k, l = hkl
        Q_mag = np.sqrt((h*a_star)**2 + (k*b_star)**2 + (l*c_star)**2)
        
        wavelength = 0.570107  # Angstroms for Si(11,11,11)
        theta_bragg = np.degrees(np.arcsin(Q_mag * wavelength / (4 * np.pi)))
        
        return {
            'tth': 2 * theta_bragg,
            'th': theta_bragg,
            'chi': 0.0,
            'phi': 0.0,
            'mu': -0.1719,
            'gam': 0.0
        }
    
    def save_configuration(self, filename: str):
        """Save current sixcircle configuration"""
        if self.simulation_mode:
            print(f"SIMULATION: Would save to {filename}")
            return
        
        self.sixc.save(filename)
        print(f"Configuration saved to {filename}")
    
    def load_configuration(self, filename: str):
        """Load saved sixcircle configuration"""
        if self.simulation_mode:
            print(f"SIMULATION: Would load from {filename}")
            return
        
        self.sixc.load(filename)
        print(f"Configuration loaded from {filename}")
