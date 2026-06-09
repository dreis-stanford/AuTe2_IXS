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
        """Load sixcircle module from configured path"""
        import sys
        import os
        
        sixcircle_path = config.SIXCIRCLE_PATH
        
        if not sixcircle_path or not os.path.exists(sixcircle_path):
            raise ImportError(f"Sixcircle not found at {sixcircle_path}")
        
        sys.path.insert(0, sixcircle_path)
        
        # Temporarily rename config files to prevent auto-load
        config_files = ['ini.conf', 'sixcircle_last_UB']
        backups = []
        for conf_file in config_files:
            if os.path.exists(conf_file):
                backup = conf_file + '.bak'
                os.rename(conf_file, backup)
                backups.append((conf_file, backup))
        
        try:
            # Import sixcircle (won't auto-load config now)
            import sixcircle_rqd
            import sixcircle
            
            # Manually initialize all required variables
            sixcircle.LAMBDA = 0.570118589663309
            sixcircle.g_sample = 'AuTe2_exp'
            sixcircle.g_aa = 7.189
            sixcircle.g_bb = 4.407
            sixcircle.g_cc = 5.069
            sixcircle.g_al = 90.0
            sixcircle.g_be = 89.96
            sixcircle.g_ga = 90.0
            
            # Orientation reflection 0: (0,0,1) at chi=90
            sixcircle.g_h0, sixcircle.g_k0, sixcircle.g_l0 = 0.0, 0.0, 1.0
            sixcircle.g_u00, sixcircle.g_u01, sixcircle.g_u02 = 6.45, 3.225, 90.0
            sixcircle.g_u03, sixcircle.g_u04, sixcircle.g_u05 = 0.0, -0.1719, 0.0
            
            # Orientation reflection 1: (0,2,0) at chi=0, phi=90
            sixcircle.g_h1, sixcircle.g_k1, sixcircle.g_l1 = 0.0, 2.0, 0.0
            sixcircle.g_u10, sixcircle.g_u11, sixcircle.g_u12 = 14.84, 7.42, 0.0
            sixcircle.g_u13, sixcircle.g_u14, sixcircle.g_u15 = 90.0, -0.1719, 0.0
            
            # Azimuthal reference
            sixcircle.g_haz, sixcircle.g_kaz, sixcircle.g_laz = 0.0, 0.0, 1.0
            
            # Set angle limits (required for ca() function)
            sixcircle.L_TTH = -180.0
            sixcircle.U_TTH = 180.0
            sixcircle.L_TH = -180.0
            sixcircle.U_TH = 180.0
            sixcircle.L_CHI = -180.0
            sixcircle.U_CHI = 180.0
            sixcircle.L_PHI = -180.0
            sixcircle.U_PHI = 180.0
            sixcircle.F_TTH = 0
            sixcircle.F_TH = 0
            sixcircle.F_CHI = 0
            sixcircle.F_PHI = 0
            sixcircle.L_MU = -180.0
            sixcircle.U_MU = 180.0
            sixcircle.L_GAM = -180.0
            sixcircle.U_GAM = 180.0
            sixcircle.F_MU = 0
            sixcircle.F_GAM = 0
            sixcircle.L_ALPHA = -180.0
            sixcircle.U_ALPHA = 180.0
            sixcircle.L_BETA = -180.0
            sixcircle.U_BETA = 180.0
            sixcircle.F_ALPHA = 0
            sixcircle.F_BETA = 0
            
            # Calculate UB matrix
            sixcircle.UB()
            
            self.sixc = sixcircle
            self.sixc_rqd = sixcircle_rqd
            
        finally:
            # Restore config files
            for conf_file, backup in backups:
                if os.path.exists(backup):
                    os.rename(backup, conf_file)
        
        print("✓ Sixcircle modules loaded successfully")

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
            # Calculate angles quietly (caller will print if desired)
            angles = self._simulate_angles(hkl)
            return angles
        
        if check_only:
            # Capture output from ca() to parse angles
            import sys
            import io
            from contextlib import redirect_stdout
            
            f = io.StringIO()
            with redirect_stdout(f):
                self.sixc.ca(*hkl)
            output = f.getvalue()
            
            # Parse angles from output
            # Look for the line with angle values (after "tth        th ...")
            angles = {}
            lines = output.split('\n')
            for i, line in enumerate(lines):
                if 'tth' in line and 'th' in line and 'chi' in line:
                    # Next non-empty line has the values
                    for j in range(i+1, len(lines)):
                        if lines[j].strip():
                            values = lines[j].split()
                            if len(values) >= 6:
                                try:
                                    angles = {
                                        'tth': float(values[0]),
                                        'th': float(values[1]),
                                        'chi': float(values[2]),
                                        'phi': float(values[3]),
                                        'mu': float(values[4]),
                                        'gam': float(values[5])
                                    }
                                except (ValueError, IndexError):
                                    pass
                            break
                    break
            
            return angles
        else:
            self.sixc.br(*hkl)
            return self.get_current_angles()
    
    def get_current_angles(self) -> Dict[str, float]:
        """Get current diffractometer angles"""
        if self.simulation_mode:
            return {'tth': 0, 'th': 0, 'chi': 0, 'phi': 0, 'mu': -0.17, 'gam': 0}
        
        # Extract angles from sixcircle module variables
        return {
            'tth': self.sixc.TTH,
            'th': self.sixc.TH,
            'chi': self.sixc.CHI,
            'phi': self.sixc.PHI,
            'mu': self.sixc.MU,
            'gam': self.sixc.GAM
        }
    
    def _setup_simulation_UB(self):
        """
        Set up UB matrix for simulation mode
        Assumes: c||z, b||y, a||x at th=chi=phi=mu=gam=0
        """
        import numpy as np
        sys.path.insert(0, config.SIXCIRCLE_PATH)
        import sixcircle
        import sixcircle_rqd
        
        # Manually initialize required sixcircle globals (avoid ini() which loads files)
        sixcircle.TTH = 0.0
        sixcircle.TH = 0.0
        sixcircle.CHI = 0.0
        sixcircle.PHI = 0.0
        sixcircle.MU = 0.0
        sixcircle.GAM = 0.0
        sixcircle.LAMBDA = 0.570107  # Si(11,11,11)
        sixcircle.g_frozen = '456'  # Freeze mu, gam, omega
        sixcircle.PRE = 4
        
        # Set lattice parameters
        lattice = config.get_lattice_parameters()
        sixcircle.setlat(
            lattice['a'], lattice['b'], lattice['c'],
            lattice['alpha'], lattice['beta'], lattice['gamma']
        )
        
        # Set wavelength
        sixcircle_rqd.setorder(config.SI_ORDER)
        
        # Calculate (0,0,1) angles
        from .force_constants import ForceConstants
        from .aute2_structure import aute2_conv2prim_k
        
        fc_file = "data/AuTe_2_m.fc"
        xtal = ForceConstants(fc_file)
        
        # (0,0,1) in conventional -> primitive -> Cartesian
        Q_001_prim = aute2_conv2prim_k(np.array([0, 0, 1]))
        Q_001_cart = Q_001_prim @ xtal.b_l
        Q_001_mag = np.linalg.norm(Q_001_cart)
        
        wavelength = 0.570107  # Si(11,11,11)
        theta_001 = np.degrees(np.arcsin(Q_001_mag * wavelength / (4 * np.pi)))
        
        # For c||z: need chi=90deg to bring c into horizontal scattering plane
        or0_angles = {
            'tth': 2 * theta_001,
            'th': theta_001,
            'chi': 90.0,
            'phi': 0.0,
            'mu': 0.0,
            'gam': 0.0
        }
        
        # Calculate (1,0,0) angles  
        Q_100_prim = aute2_conv2prim_k(np.array([1, 0, 0]))
        Q_100_cart = Q_100_prim @ xtal.b_l
        Q_100_mag = np.linalg.norm(Q_100_cart)
        theta_100 = np.degrees(np.arcsin(Q_100_mag * wavelength / (4 * np.pi)))
        
        # For a||x (perpendicular to beam): chi=0, need to rotate about z
        or1_angles = {
            'tth': 2 * theta_100,
            'th': theta_100 - 90.0,
            'chi': 0.0,
            'phi': 0.0,
            'mu': 0.0,
            'gam': 0.0
        }
        
        # Set motor positions and orientation reflections
        sixcircle.TTH = or0_angles['tth']
        sixcircle.TH = or0_angles['th']
        sixcircle.CHI = or0_angles['chi']
        sixcircle.PHI = or0_angles['phi']
        sixcircle.MU = or0_angles['mu']
        sixcircle.GAM = or0_angles['gam']
        sixcircle.or0(0, 0, 1)
        
        sixcircle.TTH = or1_angles['tth']
        sixcircle.TH = or1_angles['th']
        sixcircle.CHI = or1_angles['chi']
        sixcircle.PHI = or1_angles['phi']
        sixcircle.MU = or1_angles['mu']
        sixcircle.GAM = or1_angles['gam']
        sixcircle.or1(1, 0, 0)
        
        print("UB matrix configured (simulation)")
        print(f"  or0: (0,0,1) at chi=90deg, th={theta_001:.2f}deg")
        print(f"  or1: (1,0,0) at chi=0deg, th={theta_100-90:.2f}deg")
        
        self._ub_initialized = True

    def _simulate_angles(self, hkl):
        """Calculate angles using sixcircle UB matrix in simulation mode"""
        import numpy as np
        
        if not hasattr(self, '_ub_initialized') or not self._ub_initialized:
            self._setup_simulation_UB()
        
        sys.path.insert(0, config.SIXCIRCLE_PATH)
        import sixcircle
        
        # Use sixcircle to calculate angles
        h, k, l = hkl
        try:
            sixcircle.ca(h, k, l)
        except Exception as e:
            import traceback
            print(f"\nError in sixcircle.ca({h}, {k}, {l}):")
            traceback.print_exc()
            raise
        
        # Return calculated angles
        return {
            'tth': sixcircle.A[0],
            'th': sixcircle.A[1],
            'chi': sixcircle.A[2],
            'phi': sixcircle.A[3],
            'mu': sixcircle.A[4],
            'gam': sixcircle.A[5]
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
