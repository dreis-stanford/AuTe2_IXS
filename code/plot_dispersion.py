"""
Plot phonon dispersion for AuTe2
Converted from MATLAB dispersion plotting script
"""

import numpy as np
import matplotlib.pyplot as plt
from force_constants import ForceConstants
from phonons import calc_Dq, calc_freq_eig, convert_frequencies, generate_bz_grid

def plot_aute2_dispersion(fc_file='data/AuTe_2_m.fc', n_points=201, show=True, block=False):
    """
    Plot phonon dispersion along high-symmetry path for AuTe2
    
    Parameters:
    -----------
    fc_file : str
        Path to force constants file
    n_points : int
        Number of points between high-symmetry points
    show : bool
        Show the plot (default True)
    block : bool
        Block execution until plot is closed (default False)
    """
    
    print("=" * 70)
    print("AuTe2 Phonon Dispersion Calculation")
    print("=" * 70)
    
    # Load force constants
    xtal = ForceConstants(fc_file)
    Phi = xtal.convert_to_eV_per_Angstrom2()
    
    # Get masses for each atom
    masses = np.array([xtal.masses[xtal.atom_type_map[i]-1] 
                       for i in range(xtal.nat)])
    
    print(f"\nCalculating dispersion with {n_points} points per segment...")
    
    # Define high-symmetry points in primitive coordinates
    QG = np.array([0.0, 0.0, 0.0])  # Γ
    QA = np.array([0.5, 0.0, 0.0])  # A
    QM = np.array([0.5, 0.5, 0.5])  # M
    QY = np.array([0.0, 0.5, 0.5])  # Y
    QV = np.array([0.0, 0.5, 0.0])  # V
    QL = np.array([0.5, 0.5, 0.0])  # L
    QZ = np.array([0.0, 0.0, 0.5])  # Z
    
    # Create path: Z→Γ→A→M→Y→Γ→V→L→A
    xi = np.linspace(0, 1, n_points, endpoint=False)
    
    path_segments = [
        (QZ, QG, 'Z', 'Γ'),
        (QG, QA, 'Γ', 'A'),
        (QA, QM, 'A', 'M'),
        (QM, QY, 'M', 'Y'),
        (QY, QG, 'Y', 'Γ'),
        (QG, QV, 'Γ', 'V'),
        (QV, QL, 'V', 'L'),
        (QL, QA, 'L', 'A'),
    ]
    
    q_path = []
    labels = []
    segment_indices = [0]
    
    for q_start, q_end, label_start, label_end in path_segments:
        # Generate points along segment
        for i in range(n_points):
            q_path.append(q_start + xi[i] * (q_end - q_start))
        
        if len(labels) == 0:
            labels.append(label_start)
        labels.append(label_end)
        segment_indices.append(len(q_path))
    
    q_path = np.array(q_path)
    
    print(f"  Total path points: {len(q_path)}")
    print(f"  Path: {' → '.join(labels)}")
    
    # Calculate dynamical matrices along path
    print("\nCalculating dynamical matrices...")
    D_path = calc_Dq(q_path, xtal.uvw, Phi, masses)
    
    # Solve for frequencies and eigenvectors
    print("Solving for phonon modes...")
    w_raw, ev = calc_freq_eig(D_path, sort_modes=True, rotate_ev=True)
    
    # Convert to cm^-1
    w_cm = convert_frequencies(w_raw, 'cm-1')
    
    print(f"  Frequency range: {w_cm.min():.2f} to {w_cm.max():.2f} cm⁻¹")
    
    # Calculate DOS
    print("\nCalculating density of states...")
    Q_BZ = generate_bz_grid(15)
    D_BZ = calc_Dq(Q_BZ, xtal.uvw, Phi, masses)
    w_BZ_raw, _ = calc_freq_eig(D_BZ, sort_modes=False)
    w_BZ_cm = convert_frequencies(w_BZ_raw, 'cm-1')
    
    # Create histogram for DOS
    DOS, bins = np.histogram(np.real(w_BZ_cm.flatten()), bins=121)
    
    # Create plot
    print("\nCreating plot...")
    
    fig = plt.figure(figsize=(12, 8))
    
    # Layout: [left, bottom, width, height]
    ax_disp = plt.axes([0.08, 0.55, 0.68, 0.40])   # Dispersion
    ax_dos = plt.axes([0.78, 0.55, 0.18, 0.40])    # DOS
    ax_path = plt.axes([0.08, 0.35, 0.68, 0.15])   # q-components
    
    # Plot dispersion
    ax_disp.plot(range(len(q_path)), np.abs(w_cm), linewidth=2)
    
    # Add vertical lines at high-symmetry points
    for idx in segment_indices[1:-1]:
        ax_disp.axvline(x=idx, color='gray', linestyle='--', linewidth=0.5)
    
    ax_disp.set_ylabel('Frequency (cm⁻¹)', fontsize=14)
    ax_disp.set_xlim(0, len(q_path))
    ax_disp.set_xticks(segment_indices)
    ax_disp.set_xticklabels(labels)
    ax_disp.grid(True, alpha=0.3)
    ax_disp.set_title('AuTe₂ Phonon Dispersion', fontsize=16)
    
    # Plot DOS (share y-axis with dispersion)
    ax_dos.plot(DOS, bins[1:], '-', linewidth=2)
    ax_dos.set_ylim(ax_disp.get_ylim())
    ax_dos.set_xlabel('DOS', fontsize=10)
    ax_dos.set_yticklabels([])  # Remove y-tick labels
    ax_dos.grid(True, alpha=0.3)
    
    # Plot q-components
    ax_path.plot(range(len(q_path)), q_path[:, 0], 'r-', linewidth=1.5, label='q₁')
    ax_path.plot(range(len(q_path)), q_path[:, 1], 'g--', linewidth=1.5, label='q₂')
    ax_path.plot(range(len(q_path)), q_path[:, 2], 'b:', linewidth=2, label='q₃')
    
    for idx in segment_indices[1:-1]:
        ax_path.axvline(x=idx, color='gray', linestyle='--', linewidth=0.5)
    
    ax_path.set_ylabel('q (r.l.u.)', fontsize=11)
    ax_path.set_xlim(0, len(q_path))
    ax_path.set_xticks(segment_indices)
    ax_path.set_xticklabels([])
    ax_path.set_ylim(q_path.min() - 0.05, q_path.max() + 0.05)
    ax_path.legend(loc='upper right', fontsize=9, framealpha=0.9)
    ax_path.grid(True, alpha=0.3)
    
    plt.savefig('aute2_dispersion.png', dpi=150, bbox_inches='tight')
    print(f"\n✓ Saved plot to: aute2_dispersion.png")
    
    if show:
        plt.show(block=block)
    
    return fig, (ax_disp, ax_dos, ax_path), q_path, w_cm


if __name__ == "__main__":
    import sys
    import os
    
    sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
    
    # Plot dispersion (non-blocking by default)
    fig, axes, q_path, w_cm = plot_aute2_dispersion(block=False)
    
    print("\n" + "=" * 70)
    print("Dispersion calculation complete!")
    print("=" * 70)
    print("\nPlot window is open (non-blocking).")
    print("Close the plot window or press Ctrl+C to exit.")
    
    # Keep script running so plot stays visible
    try:
        plt.pause(0.1)  # Small pause to ensure plot displays
        input("\nPress Enter to exit...")
    except KeyboardInterrupt:
        print("\n\nExiting...")
    
    plt.close('all')
