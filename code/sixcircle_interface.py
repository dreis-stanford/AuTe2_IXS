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

        # The sixcircle Python package is pure calculation (mv()/br() only
        # update module-level globals — no hardware I/O), so we load it
        # whenever it's available, even in simulation mode, to compute real
        # diffractometer angles. `simulation_mode` only gates the `move`
        # command in the interactive tool.
        if config.SIXCIRCLE_AVAILABLE:
            try:
                self._load_sixcircle()
            except Exception as e:
                if not simulation_mode:
                    raise
                print(f"⚠ Could not load sixcircle ({e}); "
                      f"angle calculation unavailable.")
    
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

            # Populate the BL43LXU analyzer-array geometry globals (x_spac,
            # z_spac, a_radi, x_n, z_n, x_cen, z_cen, ...) from
            # BL43XU_CONST.mac, used by analyze_array for per-analyzer
            # (tth, gam) offsets. Explicit and idempotent, independent of
            # whatever init_rqd() did at import time via previous.bl.conf.
            sixcircle_rqd.setbl(config.BEAMLINE)

            # Manually initialize all required variables. Use setorder()
            # (not a bare sixcircle.LAMBDA assignment) so sixcircle_rqd's
            # silicon_order global -- used by pa() and analyzer_mask.py --
            # and the persisted previous.bl.conf stay in sync with
            # config.SI_ORDER (previous.bl.conf otherwise keeps reporting a
            # stale Si order from before config.SI_ORDER was last changed).
            sixcircle_rqd.setorder(config.SI_ORDER)
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
            
            # Azimuthal reference / surface normal (used for alpha, beta)
            sn = config.SURFACE_NORMAL
            sixcircle.g_haz, sixcircle.g_kaz, sixcircle.g_laz = (
                float(sn[0]), float(sn[1]), float(sn[2]))

            sixcircle.PRE = 4
            self.sixc = sixcircle
            self.sixc_rqd = sixcircle_rqd

            # Frozen mode/values and angle limits from config (these are what
            # ca()/ca6() use to pick a solution; nothing here moves hardware).
            self.apply_frozen()
            self.apply_limits()

            # Calculate UB matrix
            sixcircle.UB()

        finally:
            # Restore config files
            for conf_file, backup in backups:
                if os.path.exists(backup):
                    os.rename(backup, conf_file)
        
        print("✓ Sixcircle modules loaded successfully")

    # ------------------------------------------------------------------
    # Diffractometer modes / constraints (SPEC-like, calculation-only)
    # ------------------------------------------------------------------
    # Digit codes used by sixcircle for frozen angles:
    _ANGLE_CODES = {'0': 'tth', '1': 'th', '2': 'chi', '3': 'phi', '4': 'mu',
                    '5': 'gam', '6': 'omega', '7': 'azimuth', '8': 'alpha',
                    '9': 'beta'}
    _FROZEN_GLOBALS = {'tth': 'F_TTH', 'th': 'F_TH', 'chi': 'F_CHI',
                       'phi': 'F_PHI', 'mu': 'F_MU', 'gam': 'F_GAM',
                       'omega': 'F_OMEGA', 'azimuth': 'F_AZIMUTH',
                       'alpha': 'F_ALPHA', 'beta': 'F_BETA'}
    _LIMIT_GLOBALS = {'tth': ('L_TTH', 'U_TTH'), 'th': ('L_TH', 'U_TH'),
                      'chi': ('L_CHI', 'U_CHI'), 'phi': ('L_PHI', 'U_PHI'),
                      'mu': ('L_MU', 'U_MU'), 'gam': ('L_GAM', 'U_GAM'),
                      'alpha': ('L_ALPHA', 'U_ALPHA'),
                      'beta': ('L_BETA', 'U_BETA')}

    def apply_frozen(self):
        """Push config's frozen mode (FROZEN_ANGLES) and held values
        (FROZEN_VALUES) into the loaded sixcircle module."""
        if self.sixc is None:
            return
        self.sixc.g_frozen = config.FROZEN_ANGLES
        for name, gvar in self._FROZEN_GLOBALS.items():
            setattr(self.sixc, gvar, float(config.FROZEN_VALUES.get(name, 0.0)))

    def apply_limits(self):
        """Push config's ANGLE_LIMITS into the loaded sixcircle module."""
        if self.sixc is None:
            return
        for name, (lvar, uvar) in self._LIMIT_GLOBALS.items():
            lo, hi = config.ANGLE_LIMITS[name]
            setattr(self.sixc, lvar, float(lo))
            setattr(self.sixc, uvar, float(hi))

    def set_frozen(self, code: str):
        """Set which three angles are frozen (e.g. '456') and re-apply values."""
        config.set_frozen_angles(code)
        self.apply_frozen()

    def set_frozen_values(self, **kwargs):
        """Set the held values of frozen angles, e.g. set_frozen_values(mu=-0.17)."""
        config.set_frozen_values(**kwargs)
        self.apply_frozen()

    def set_limit(self, angle: str, lower: float, upper: float):
        """Set the (lower, upper) degree limit for one angle and apply it."""
        config.set_angle_limit(angle, lower, upper)
        self.apply_limits()

    def reset_limits(self):
        """Reset all angle limits to ±180° and apply."""
        config.reset_angle_limits()
        self.apply_limits()

    def describe_constraints(self, include_limits: bool = True) -> str:
        """Return a human-readable summary of the current frozen mode (+ limits)."""
        frozen = ", ".join(
            f"{self._ANGLE_CODES[d]}={config.FROZEN_VALUES.get(self._ANGLE_CODES[d], 0.0):.4g}"
            for d in config.FROZEN_ANGLES)
        lines = [f"Frozen ({config.FROZEN_ANGLES}): {frozen}"]
        if include_limits:
            lines.append("Limits (deg):")
            for name, (lo, hi) in config.ANGLE_LIMITS.items():
                flag = "  (default)" if (lo, hi) == (-180.0, 180.0) else ""
                lines.append(f"    {name:>5}: [{lo:8.3f}, {hi:8.3f}]{flag}")
        return "\n".join(lines)

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
    
    def calculate_angles(self, hkl: Tuple[float, float, float],
                         ignore_limits: bool = False) -> Optional[Dict[str, float]]:
        """
        Calculate the six diffractometer angles for a reflection (no motion).

        Uses sixcircle's ca_s() solver under the current frozen mode and angle
        limits. Returns a dict {tth, th, chi, phi, mu, gam, omega, alpha, beta}
        for the first valid solution, or None if the reflection is
        inaccessible within the limits. Raises RuntimeError if the sixcircle
        library is not loaded.

        If ignore_limits is True, the configured angle limits are temporarily
        widened to ±180° for this calculation (then restored), so the
        returned solution may lie outside config.ANGLE_LIMITS.
        """
        if self.sixc is None:
            raise RuntimeError(
                "Angle calculation requires the sixcircle library, which is "
                "not loaded (check config.SIXCIRCLE_PATH).")
        h, k, l = (float(hkl[0]), float(hkl[1]), float(hkl[2]))
        if ignore_limits:
            saved = {}
            for name, (lvar, uvar) in self._LIMIT_GLOBALS.items():
                saved[name] = (getattr(self.sixc, lvar), getattr(self.sixc, uvar))
                setattr(self.sixc, lvar, -180.0)
                setattr(self.sixc, uvar, 180.0)
            try:
                flag, pos = self.sixc.ca_s(h, k, l)
            finally:
                for name, (lvar, uvar) in self._LIMIT_GLOBALS.items():
                    lo, hi = saved[name]
                    setattr(self.sixc, lvar, lo)
                    setattr(self.sixc, uvar, hi)
        else:
            flag, pos = self.sixc.ca_s(h, k, l)
        if not flag or not pos:
            return None
        tth, th, chi, phi, mu, gam = pos[0][:6]
        omega = pos[0][7]
        alpha, beta = pos[0][9], pos[0][10]
        return {'tth': tth, 'th': th, 'chi': chi,
                'phi': phi, 'mu': mu, 'gam': gam, 'omega': omega,
                'alpha': alpha, 'beta': beta}

    def analyzer_array_offsets(self, hkl: Tuple[float, float, float]
                                ) -> Optional[Dict[str, dict]]:
        """
        Compute (dH, dK, dL) and (dtth, dgam) offsets from hkl for each of the
        28 analyzers in the BL43LXU array, for the arm centered on hkl.

        For each analyzer, tth and gam are offset from the center position by
        the fixed angular steps set by the array geometry (x_spac, z_spac,
        a_radi, x_n, z_n, x_cen, z_cen from BL43XU_CONST.mac via setbl()),
        then converted to (H,K,L) via the real UB matrix -- the same
        (tth, gam) -> (H,K,L) forward-kinematics ca6() uses, without its
        printing/file output. Analyzer names follow ca6()'s convention
        (a03-a09, a14-a20, a25-a31, a35-a41).

        Returns a dict {analyzer_name: info}, or None if hkl is inaccessible
        under the current frozen angles/limits. Raises RuntimeError if
        sixcircle is not loaded. Each `info` dict has:
          - 'dQ':   np.array([dH, dK, dL]), offset from the center Q
          - 'dtth', 'dgam': angular offsets (degrees) from the center
            (tth, gam), i.e. the analyzer array's effective 2D layout.
            dtth is symmetric about 0 (x_cen is the middle column), but with
            BL43LXU's current geometry (z_cen=0, the first row), dgam is <= 0
            -- the array extends to one side (lower gamma) only, not
            symmetrically. This is the real array geometry, not a bug.
          - 'xi', 'zi': horizontal/vertical grid indices (0-based) into the
            x_n x z_n analyzer array
        """
        if self.sixc is None:
            raise RuntimeError(
                "Analyzer array calculation requires the sixcircle library, "
                "which is not loaded (check config.SIXCIRCLE_PATH).")
        angles = self.calculate_angles(hkl)
        if angles is None:
            return None

        sixc = self.sixc
        rqd = self.sixc_rqd
        ca_tth, ca_gam = angles['tth'], angles['gam']

        o_FLAG_WH = sixc.FLAG_WH
        sixc.FLAG_WH = False
        try:
            sixc.mv(tth=ca_tth, th=angles['th'], chi=angles['chi'],
                    phi=angles['phi'], mu=angles['mu'], gam=ca_gam)
            Q_center = np.array([sixc.H, sixc.K, sixc.L])

            offsets = {}
            for zi in range(rqd.z_n):
                for xi in range(rqd.x_n):
                    if zi <= 2:
                        name = 'a' + str(11 * zi + 6 + (xi - rqd.x_cen)).zfill(2)
                    else:
                        name = 'a' + str(11 * zi + 6 + (xi - rqd.x_cen) - 1).zfill(2)

                    A_x = rqd.x_off + rqd.x_spac * (xi - rqd.x_cen)
                    A_z = rqd.z_off + rqd.z_spac * (rqd.z_cen - zi)
                    # sixcircle's mv() requires plain Python floats (it
                    # rejects numpy scalars with "Invalid motor position")
                    A_tth = float(ca_tth + np.degrees(np.arctan2(A_x, rqd.a_radi)))
                    A_gam = float(ca_gam + np.degrees(
                        np.arctan2(A_z, np.sqrt(rqd.a_radi**2 + A_x**2))))

                    sixc.mv(tth=A_tth, gam=A_gam)
                    Q_analyzer = np.array([sixc.H, sixc.K, sixc.L])
                    offsets[name] = {
                        'dQ': Q_analyzer - Q_center,
                        'dtth': A_tth - ca_tth,
                        'dgam': A_gam - ca_gam,
                        'xi': xi,
                        'zi': zi,
                    }

            # Restore center position
            sixc.mv(tth=ca_tth, th=angles['th'], chi=angles['chi'],
                    phi=angles['phi'], mu=angles['mu'], gam=ca_gam)
        finally:
            sixc.FLAG_WH = o_FLAG_WH

        return offsets

    def get_UB(self):
        """Return the 3x3 UB matrix (U @ B) currently loaded in sixcircle, or
        None if unavailable. Columns of UB are a*, b*, c* in the phi frame;
        UB @ (h,k,l) is the reciprocal vector, inv(UB).T columns are real a,b,c."""
        if self.sixc is None:
            return None
        return np.array(self.sixc.M_U, float) @ np.array(self.sixc.M_B, float)

    def move_to_hkl(self, hkl: Tuple[float, float, float],
                    check_only: bool = False):
        """
        Calculate angles for a given HKL.

        This interface is calculation-only and never drives the diffractometer;
        `check_only` is accepted for backward compatibility but has no effect.
        Returns the calculated angles (see calculate_angles).
        """
        return self.calculate_angles(hkl)

    def get_current_angles(self) -> Dict[str, float]:
        """Get the diffractometer angles currently stored in the sixcircle module."""
        if self.sixc is None:
            return {'tth': 0, 'th': 0, 'chi': 0, 'phi': 0, 'mu': 0, 'gam': 0}

        # Extract angles from sixcircle module variables
        return {
            'tth': self.sixc.TTH,
            'th': self.sixc.TH,
            'chi': self.sixc.CHI,
            'phi': self.sixc.PHI,
            'mu': self.sixc.MU,
            'gam': self.sixc.GAM
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
