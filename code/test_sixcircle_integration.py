"""
Test and demonstrate sixcircle integration for AuTe2 IXS
"""

import numpy as np
from . import config
from .sixcircle_interface import SixCircleInterface
from .modulated_structure import ModulatedStructure
from .q_optimizer import QPointOptimizer

def test_configuration():
    """Test configuration loading"""
    print("\n" + "="*70)
    print("Testing Configuration")
    print("="*70)
    
    config.print_config()
    
    # Test updating parameters
    print("Updating lattice parameter a to 8.16...")
    config.update_lattice_parameters(a=8.16)
    print(f"New a = {config.LATTICE_PARAMS['a']}")
    
    # Reset
    config.update_lattice_parameters(a=8.15)

def test_sixcircle_interface():
    """Test sixcircle interface in simulation mode"""
    print("\n" + "="*70)
    print("Testing Sixcircle Interface")
    print("="*70)
    
    # Initialize (will auto-detect simulation mode)
    sixc = SixCircleInterface()
    
    # Setup experiment
    sixc.setup_experiment()
    
    # Test orientation setting
    print("\nSetting orientation matrix...")
    sixc.set_orientation(
        or0_hkl=(2, 0, 0),
        or0_angles={'tth': 25.3, 'th': 12.65, 'chi': 0, 'phi': 0},
        or1_hkl=(0, 2, 0),
        or1_angles={'tth': 26.1, 'th': 13.05, 'chi': 0, 'phi': 90}
    )
    
    # Test moving to HKL
    print("\nCalculating angles for (1, 1, 0)...")
    angles = sixc.move_to_hkl((1, 1, 0), check_only=True)
    
    # Test IXS positions
    print("\nCalculating IXS analyzer positions...")
    ixs_info = sixc.calculate_ixs_positions((1.1, 1.1, 0))
    
    if 'analyzers' in ixs_info:
        print(f"\nFound {len(ixs_info['analyzers'])} analyzer positions")
        for name, info in list(ixs_info['analyzers'].items())[:3]:
            print(f"  {name}: HKL = {info['hkl']}")

def test_modulated_structure():
    """Test modulated structure handling"""
    print("\n" + "="*70)
    print("Testing Modulated Structure")
    print("="*70)
    
    mod = ModulatedStructure()
    
    # Print info
    mod.print_satellite_info(max_order=2)
    
    # Generate satellite list
    print("Generating satellite reflections...")
    reflections = mod.generate_reflection_list(
        h_range=(0, 3),
        k_range=(0, 3),
        l_range=(0, 2),
        max_satellite_order=1,
        max_q=40.0
    )
    
    print(f"Found {len(reflections)} reflections (main + satellites)")
    
    # Count by type
    n_main = sum(1 for r in reflections if r['type'] == 'main')
    n_sat = sum(1 for r in reflections if r['type'] == 'satellite')
    print(f"  Main reflections: {n_main}")
    print(f"  Satellites: {n_sat}")
    
    # Show some examples
    print("\nFirst 5 reflections:")
    for r in reflections[:5]:
        print(f"  {r['type']:10s} {r['hkl']}, |Q| = {r['q_magnitude']:.2f} nm^-1, order = {r['order']}")
    
    # Find satellites near a phonon Q
    print("\nFinding satellites near q = (0.3, 0, 0)...")
    near_q = mod.find_satellites_near_q((0.3, 0, 0), tolerance=0.1)
    print(f"Found {len(near_q)} reflections within 0.1 r.l.u.")
    for r in near_q[:3]:
        print(f"  {r['hkl']}, distance = {r['distance']:.4f} r.l.u.")
    
    # Check Bragg accessibility
    print("\nChecking Bragg-accessible satellites...")
    accessible = mod.bragg_accessible_satellites(max_2theta=140.0)
    print(f"Found {len(accessible)} accessible reflections (2θ < 140°)")
    print("\nFirst 5 by scattering angle:")
    for r in accessible[:5]:
        print(f"  {r['hkl']}, 2θ = {r['two_theta']:.2f}°, {r['type']}")

def test_q_optimizer():
    """Test Q-point optimization"""
    print("\n" + "="*70)
    print("Testing Q-Point Optimizer")
    print("="*70)
    
    opt = QPointOptimizer()
    
    # Define some phonon Q-points to measure (along Γ-X)
    q_points = [
        (0.1, 0, 0),
        (0.2, 0, 0),
        (0.3, 0, 0),
        (0.4, 0, 0),
        (0.5, 0, 0),  # X point
    ]
    
    # Define available Bragg peaks
    bragg_peaks = [
        (1, 0, 0),
        (2, 0, 0),
        (3, 0, 0),
        (1, 1, 0),
        (2, 1, 0),
        (0, 1, 0),
        (1, 0, 1),
    ]
    
    print(f"\nPhonon Q-points: {len(q_points)}")
    print(f"Available Bragg peaks: {len(bragg_peaks)}")
    
    # Rank Bragg peaks for one Q-point
    print("\nRanking Bragg peaks for q = (0.3, 0, 0)...")
    ranked = opt.rank_q_points_for_mode(
        q_mode=(0.3, 0, 0),
        branch_index=0,
        candidate_q_points=bragg_peaks
    )
    
    print("\nTop 3 Bragg peaks:")
    for i, r in enumerate(ranked[:3], 1):
        sens = opt.longitudinal_vs_transverse(r['q_phonon'], r['Q_bragg'])
        print(f"  {i}. Q_Bragg = {r['Q_bragg']}")
        print(f"     Q_transfer = {r['Q_transfer']}")
        print(f"     Intensity = {r['expected_intensity']:.3f}")
        print(f"     Type: {sens['mode_type']} (angle = {sens['Q_angle_to_q']:.1f}°)")
    
    # Get measurement plan
    print("\nGenerating optimal measurement plan...")
    opt.print_optimization_summary(q_points, bragg_peaks)

def run_all_tests():
    """Run all integration tests"""
    print("\n" + "#"*70)
    print("# AuTe2 IXS - Sixcircle Integration Tests")
    print("#"*70)
    
    test_configuration()
    test_sixcircle_interface()
    test_modulated_structure()
    test_q_optimizer()
    
    print("\n" + "#"*70)
    print("# All tests completed!")
    print("#"*70 + "\n")

if __name__ == '__main__':
    run_all_tests()
