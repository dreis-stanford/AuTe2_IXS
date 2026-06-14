"""
Optimize Q-point selection for IXS measurements
Determine which Q-points give maximum information about phonon eigenvectors
"""

import numpy as np
from typing import List, Tuple, Dict, Optional

from . import config

class QPointOptimizer:
    """
    Optimize selection of Q-points for IXS measurements
    to maximize information about phonon modes
    """
    
    def __init__(self, structure_factor_calculator=None):
        """
        Parameters:
        -----------
        structure_factor_calculator : Object that can calculate S(Q,ω)
                                     (your existing DFT/phonon code)
        """
        self.sf_calc = structure_factor_calculator
    
    def rank_q_points_for_mode(self,
                               q_mode: Tuple[float, float, float],
                               branch_index: int,
                               candidate_q_points: List[Tuple[float, float, float]],
                               eigenvector: np.ndarray = None) -> List[Dict]:
        """
        Rank candidate Q-points by expected intensity for a phonon mode
        
        Parameters:
        -----------
        q_mode : Phonon wavevector in r.l.u.
        branch_index : Which phonon branch
        candidate_q_points : List of Bragg peaks to consider
        eigenvector : Phonon eigenvector (if known from DFT)
        
        Returns:
        --------
        List of dicts sorted by expected intensity, each containing:
            'Q_bragg' : Bragg peak position
            'Q_transfer' : Q_bragg - q_mode (momentum transfer)
            'expected_intensity' : Relative intensity
            'scattering_vector' : Q direction
        """
        ranked = []
        
        for Q_bragg in candidate_q_points:
            Q_bragg_vec = np.array(Q_bragg)
            q_mode_vec = np.array(q_mode)
            
            # Momentum transfer
            Q_transfer = Q_bragg_vec - q_mode_vec
            
            # Calculate expected intensity
            if eigenvector is not None and self.sf_calc is not None:
                # Use actual structure factor calculation
                intensity = self._calculate_structure_factor(
                    Q_transfer, q_mode_vec, branch_index, eigenvector
                )
            else:
                # Simple geometric estimate: intensity ~ |Q·e|^2
                # where e is displacement direction
                intensity = self._estimate_intensity(Q_transfer, eigenvector)
            
            ranked.append({
                'Q_bragg': Q_bragg,
                'Q_transfer': tuple(Q_transfer),
                'q_phonon': q_mode,
                'branch': branch_index,
                'expected_intensity': intensity,
                'Q_magnitude': np.linalg.norm(Q_transfer)
            })
        
        # Sort by expected intensity (highest first)
        ranked.sort(key=lambda x: x['expected_intensity'], reverse=True)
        
        return ranked
    
    def _estimate_intensity(self, Q: np.ndarray, eigenvector: np.ndarray = None) -> float:
        """
        Estimate relative intensity using simple model
        
        For one-phonon scattering: I ~ |Q·e|^2 * exp(-2W)
        """
        if eigenvector is None:
            # Without eigenvector info, assume isotropic
            return 1.0
        
        # Normalize Q
        Q_norm = Q / np.linalg.norm(Q) if np.linalg.norm(Q) > 0 else Q
        
        # Project eigenvector onto Q direction
        # eigenvector shape should be (n_atoms, 3)
        if eigenvector.ndim == 1:
            eigenvector = eigenvector.reshape(-1, 3)
        
        # Sum over all atoms
        intensity = 0.0
        for atom_displacement in eigenvector:
            projection = np.dot(Q_norm, atom_displacement)
            intensity += projection**2
        
        return intensity
    
    def _calculate_structure_factor(self, Q, q, branch, eigenvector):
        """
        Calculate structure factor using full formula
        This should interface with your existing IXS calculation code
        """
        if self.sf_calc is None:
            return self._estimate_intensity(Q, eigenvector)
        
        # Interface with your existing code
        # return self.sf_calc.calculate_one_phonon(Q, q, branch)
        return 1.0
    
    def optimal_q_coverage(self,
                          q_path: List[Tuple[float, float, float]],
                          available_bragg_peaks: List[Tuple[float, float, float]],
                          min_intensity_threshold: float = 0.1) -> Dict:
        """
        Find optimal set of Bragg peaks to cover a Q-path
        
        Parameters:
        -----------
        q_path : List of Q-points along dispersion to measure
        available_bragg_peaks : List of accessible Bragg reflections
        min_intensity_threshold : Minimum relative intensity to consider
        
        Returns:
        --------
        Dict mapping each q-point to best Bragg peak(s)
        """
        coverage = {}
        
        for q in q_path:
            # Rank all Bragg peaks for this q
            ranked = self.rank_q_points_for_mode(
                q, branch_index=0,  # Could iterate over branches
                candidate_q_points=available_bragg_peaks
            )
            
            # Filter by intensity threshold
            good_peaks = [r for r in ranked 
                         if r['expected_intensity'] >= min_intensity_threshold]
            
            if good_peaks:
                coverage[q] = {
                    'primary': good_peaks[0],
                    'alternatives': good_peaks[1:3],  # Next 2 best
                    'n_accessible': len(good_peaks)
                }
            else:
                coverage[q] = {
                    'primary': None,
                    'alternatives': [],
                    'n_accessible': 0
                }
        
        return coverage
    
    def suggest_measurement_plan(self,
                                q_points: List[Tuple[float, float, float]],
                                bragg_peaks: List[Tuple[float, float, float]],
                                max_measurements: int = 20) -> List[Dict]:
        """
        Suggest an optimal measurement plan
        
        Tries to maximize information while minimizing number of 
        sample reorientations
        
        Parameters:
        -----------
        q_points : Phonon Q-points to measure
        bragg_peaks : Available Bragg peaks
        max_measurements : Maximum number of measurements
        
        Returns:
        --------
        List of suggested measurements with priorities
        """
        measurements = []
        
        # For each Q, find best Bragg peak
        for q in q_points:
            ranked = self.rank_q_points_for_mode(
                q, 0, bragg_peaks
            )
            
            if ranked:
                best = ranked[0]
                measurements.append({
                    'q_phonon': q,
                    'Q_bragg': best['Q_bragg'],
                    'priority': best['expected_intensity'],
                    'Q_transfer': best['Q_transfer']
                })
        
        # Sort by priority
        measurements.sort(key=lambda x: x['priority'], reverse=True)
        
        # Limit to max_measurements
        measurements = measurements[:max_measurements]
        
        # Group by nearby Bragg peaks (to minimize reorientations)
        measurements = self._group_by_orientation(measurements)
        
        return measurements
    
    def _group_by_orientation(self, measurements: List[Dict]) -> List[Dict]:
        """
        Group measurements that can be done with similar sample orientations
        """
        # Simple implementation: sort by Bragg peak position
        # More sophisticated version would cluster by orientation matrix
        
        for i, meas in enumerate(measurements):
            meas['group'] = i  # Placeholder
            meas['reorientation_needed'] = True  # Placeholder
        
        return measurements
    
    def longitudinal_vs_transverse(self,
                                  q: Tuple[float, float, float],
                                  Q_bragg: Tuple[float, float, float]) -> Dict:
        """
        Determine if measurement is more sensitive to longitudinal
        or transverse phonons
        
        Longitudinal: polarization || q
        Transverse: polarization ⊥ q
        
        Parameters:
        -----------
        q : Phonon wavevector
        Q_bragg : Bragg peak for measurement
        
        Returns:
        --------
        Dict with sensitivity estimates
        """
        q_vec = np.array(q)
        Q_vec = np.array(Q_bragg)
        Q_transfer = Q_vec - q_vec
        
        # Normalize
        q_hat = q_vec / np.linalg.norm(q_vec) if np.linalg.norm(q_vec) > 0 else q_vec
        Q_hat = Q_transfer / np.linalg.norm(Q_transfer) if np.linalg.norm(Q_transfer) > 0 else Q_transfer
        
        # Sensitivity to longitudinal mode (polarization || q)
        # I_L ~ |Q · q_hat|^2
        long_sensitivity = np.abs(np.dot(Q_hat, q_hat))**2
        
        # Sensitivity to transverse modes (polarization ⊥ q)
        # I_T ~ 1 - |Q · q_hat|^2
        trans_sensitivity = 1 - long_sensitivity
        
        return {
            'longitudinal_sensitivity': long_sensitivity,
            'transverse_sensitivity': trans_sensitivity,
            'Q_angle_to_q': np.degrees(np.arccos(np.clip(np.dot(Q_hat, q_hat), -1, 1))),
            'mode_type': 'longitudinal' if long_sensitivity > trans_sensitivity else 'transverse'
        }
    
    def print_optimization_summary(self,
                                  q_points: List[Tuple[float, float, float]],
                                  bragg_peaks: List[Tuple[float, float, float]]):
        """
        Print summary of Q-point optimization
        """
        print("\n" + "="*70)
        print("Q-Point Optimization Summary")
        print("="*70)
        print(f"Phonon Q-points to measure: {len(q_points)}")
        print(f"Available Bragg peaks: {len(bragg_peaks)}")
        
        coverage = self.optimal_q_coverage(q_points, bragg_peaks)
        
        accessible = sum(1 for v in coverage.values() if v['primary'] is not None)
        print(f"Accessible with good intensity: {accessible}/{len(q_points)}")
        
        print("\nTop 5 recommended measurements:")
        plan = self.suggest_measurement_plan(q_points, bragg_peaks, max_measurements=5)
        
        for i, meas in enumerate(plan, 1):
            q = meas['q_phonon']
            Q = meas['Q_bragg']
            sens = self.longitudinal_vs_transverse(q, Q)
            
            print(f"\n  {i}. q = {q}")
            print(f"     Q_Bragg = {Q}")
            print(f"     Priority: {meas['priority']:.3f}")
            print(f"     Sensitivity: {sens['mode_type']} "
                  f"(L={sens['longitudinal_sensitivity']:.2f}, "
                  f"T={sens['transverse_sensitivity']:.2f})")
        
        print("="*70 + "\n")
