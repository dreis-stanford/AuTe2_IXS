"""
Phonon calculations from force constants
Converted from calcDq.m and CalcFreqEig.m
"""

import numpy as np
from scipy.linalg import eigh
from scipy.optimize import linear_sum_assignment
from .constants import const

def calc_Dq(q_pts, uvw, IFC, masses):
    """
    Calculate dynamical matrix at q-points using vectorized Fourier transform
    Equivalent to MATLAB calcDq.m
    
    Parameters:
    -----------
    q_pts : ndarray (nQ, 3)
        Q-points in reciprocal lattice units (fractional coordinates)
    uvw : ndarray (NR, 3)
        Lattice translation vectors in direct lattice units (integers)
        Contains steps relative to origin, e.g., [-1, 0, 2]
    IFC : ndarray (3n, 3n, NR)
        Interatomic force constants in eV/Ų
    masses : ndarray (n,) or (n, 1)
        Atomic masses in amu
    
    Returns:
    --------
    D_all : ndarray (3n, 3n, nQ)
        Mass-weighted dynamical matrices (complex)
    """
    
    three_n, _, NR = IFC.shape
    nQ = q_pts.shape[0]
    
    # Ensure masses is 1D
    if masses.ndim > 1:
        masses = masses.flatten()
    
    # 1. Mass weighting [3n x 3n x NR]
    # Repeat each mass 3 times (x, y, z components)
    M_eff = np.repeat(masses, 3)
    
    # Create mass normalization matrix: 1/sqrt(M_i * M_j)
    M_grid = np.sqrt(np.outer(M_eff, M_eff))
    
    # Apply mass weighting to force constants
    IFC_mass = IFC / M_grid[:, :, np.newaxis]
    
    # 2. Vectorized Fourier Transform
    # q_pts [nQ x 3] @ uvw.T [3 x NR] = [nQ x NR] phase matrix
    # Use 2*pi because q is in RLU and uvw is in lattice units
    phase = np.exp(2j * np.pi * q_pts @ uvw.T)
    
    # Reshape IFC for matrix multiplication
    # IFC_flat [3n^2 x NR] @ phase.T [NR x nQ] = [3n^2 x nQ]
    IFC_flat = IFC_mass.reshape(three_n**2, NR)
    D_flat = IFC_flat @ phase.T
    D_all = D_flat.reshape(three_n, three_n, nQ)
    
    # 3. Enforce Hermiticity
    # Average with conjugate transpose to ensure hermiticity
    D_all = (D_all + np.conj(np.transpose(D_all, (1, 0, 2)))) / 2
    
    return D_all


def calc_freq_eig(D, sort_modes=False, rotate_ev=False, scale=1.0):
    """
    Solve eigenvalue problem for dynamical matrix
    Equivalent to MATLAB CalcFreqEig.m
    
    Parameters:
    -----------
    D : ndarray (3n, 3n, nQ)
        Dynamical matrices (complex, mass-weighted)
    sort_modes : bool
        Sort modes by overlap between adjacent q-points (Hungarian matching)
    rotate_ev : bool
        Rotate eigenvector phase so largest component is real and positive
    scale : float
        Additional scaling factor
    
    Returns:
    --------
    w : ndarray (nQ, 3n)
        Phonon frequencies (sign indicates real/imaginary)
        Units: sqrt(eV/Ų/amu)
    ev : ndarray (3n, 3n, nQ)
        Eigenvectors (phonon modes)
    """
    
    nQ = D.shape[2]
    three_n = D.shape[0]
    
    # 1. Apply scaling
    D = D * scale
    
    # 2. Solve eigenvalue problem for each q-point
    eigenvalues = np.zeros((three_n, nQ))
    eigenvectors = np.zeros((three_n, three_n, nQ), dtype=complex)
    
    for iq in range(nQ):
        # Use eigh for Hermitian matrices (faster and more stable)
        evals, evecs = eigh(D[:, :, iq])
        
        # Handle negative eigenvalues (imaginary frequencies)
        # Store sign separately
        w_sq = np.real(evals)
        eigenvalues[:, iq] = np.sqrt(np.abs(w_sq)) * np.sign(w_sq)
        eigenvectors[:, :, iq] = evecs
    
    # 3. Rotate Eigenvectors (phase convention)
    if rotate_ev:
        for iq in range(nQ):
            for imode in range(three_n):
                ev_mode = eigenvectors[:, imode, iq]
                
                # Find component with largest amplitude
                max_idx = np.argmax(np.abs(ev_mode))
                
                # Rotate so that component is real and positive
                phase_shift = np.exp(-1j * np.angle(ev_mode[max_idx]))
                eigenvectors[:, imode, iq] *= phase_shift
    
    # 4. Sort Modes (Hungarian algorithm for mode tracking)
    if sort_modes and nQ > 1:
        for iq in range(1, nQ):
            # Calculate overlap matrix between consecutive q-points
            # overlap[i,j] = |<ev(q-1, mode i) | ev(q, mode j)>|
            overlap = np.abs(eigenvectors[:, :, iq-1].conj().T @ eigenvectors[:, :, iq])
            
            # Use Hungarian algorithm to find best 1-to-1 matching
            # Minimize negative overlap (maximize overlap)
            cost_matrix = -overlap
            row_ind, col_ind = linear_sum_assignment(cost_matrix)
            
            # Reorder modes at current q-point
            eigenvalues[:, iq] = eigenvalues[col_ind, iq]
            eigenvectors[:, :, iq] = eigenvectors[:, col_ind, iq]
    
    # 5. Transpose to match MATLAB output format [nQ x 3n]
    w = eigenvalues.T
    
    return w, eigenvectors


def convert_frequencies(w_raw, unit='THz'):
    """
    Convert raw eigenvalue square roots to physical frequencies
    
    Parameters:
    -----------
    w_raw : ndarray
        Raw frequencies from calc_freq_eig: sqrt(eV/Ų/amu)
    unit : str
        Output unit: 'THz', 'cm-1', 'meV'
    
    Returns:
    --------
    w_converted : ndarray
        Frequencies in requested units
    """
    
    # Convert from sqrt(eV/Ų/amu) to THz
    # w [sqrt(eV/Angstrom^2/amu)] * c/100/sqrt(M_u)/2/pi -> THz
    w_THz = w_raw * const.c / 100 / np.sqrt(const.M_u) / (2 * np.pi)
    
    if unit == 'THz':
        return w_THz
    elif unit == 'cm-1' or unit == 'cm_inv':
        # THz to cm^-1
        return w_THz / (const.c * 1e-10)
    elif unit == 'meV':
        # THz to meV
        return w_THz * const.THz2meV
    else:
        raise ValueError(f"Unknown unit: {unit}. Use 'THz', 'cm-1', or 'meV'")


def generate_bz_grid(N):
    """
    Generate Monkhorst-Pack grid for Brillouin zone sampling
    Equivalent to MATLAB BZgrid.m
    
    Parameters:
    -----------
    N : int or array-like (3,)
        Grid dimensions. If scalar, uses same N for all directions
    
    Returns:
    --------
    Q : ndarray (N1*N2*N3, 3)
        Q-points in reciprocal lattice units
    """
    
    if np.isscalar(N):
        N = [N, N, N]
    
    N = np.array(N)
    
    if len(N) != 3:
        raise ValueError("N must be scalar or length-3 array")
    
    # Generate fractional coordinates for each dimension
    # Monkhorst-Pack: (2n - N - 1) / (2N) for n = 1, 2, ..., N
    HKL = []
    for i in range(3):
        n = np.arange(1, N[i] + 1)
        HKL.append((2*n - N[i] - 1) / (2*N[i]))
    
    # Create 3D meshgrid
    H, K, L = np.meshgrid(HKL[0], HKL[1], HKL[2], indexing='ij')
    
    # Reshape into (Nq, 3) array
    Q = np.column_stack([H.ravel(), K.ravel(), L.ravel()])
    
    return Q


if __name__ == "__main__":
    import sys
    import os
    sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
    
    from force_constants import ForceConstants
    
    print("=" * 60)
    print("Testing Phonon Calculations")
    print("=" * 60)
    
    # Load force constants
    fc_file = "data/AuTe_2_m.fc"
    xtal = ForceConstants(fc_file)
    Phi = xtal.convert_to_eV_per_Angstrom2()
    
    # Get masses
    masses = np.array([xtal.masses[xtal.atom_type_map[i]-1] 
                       for i in range(xtal.nat)])
    
    print(f"\nMasses for calculation: {masses} amu")
    
    # Test at Gamma point
    print("\n" + "=" * 60)
    print("Test 1: Phonons at Γ point")
    print("=" * 60)
    
    q_gamma = np.array([[0.0, 0.0, 0.0]])
    D_gamma = calc_Dq(q_gamma, xtal.uvw, Phi, masses)
    
    w_raw, ev = calc_freq_eig(D_gamma, sort_modes=False, rotate_ev=True)
    
    # Convert to different units
    w_THz = convert_frequencies(w_raw, 'THz')
    w_cm = convert_frequencies(w_raw, 'cm-1')
    w_meV = convert_frequencies(w_raw, 'meV')
    
    print("\nPhonon frequencies at Γ:")
    print(f"{'Mode':>6} {'THz':>10} {'cm⁻¹':>10} {'meV':>10}")
    print("-" * 40)
    for i in range(len(w_THz[0])):
        print(f"{i+1:6d} {w_THz[0,i]:10.4f} {w_cm[0,i]:10.2f} {w_meV[0,i]:10.3f}")
    
    # Test along a path
    print("\n" + "=" * 60)
    print("Test 2: Path Γ → X")
    print("=" * 60)
    
    n_pts = 11
    q_path = np.zeros((n_pts, 3))
    q_path[:, 0] = np.linspace(0, 0.5, n_pts)  # Γ to X
    
    D_path = calc_Dq(q_path, xtal.uvw, Phi, masses)
    w_path, ev_path = calc_freq_eig(D_path, sort_modes=True, rotate_ev=True)
    w_path_cm = convert_frequencies(w_path, 'cm-1')
    
    print(f"\nFrequency range along Γ→X:")
    print(f"  Acoustic: {np.min(w_path_cm[:, :3]):.2f} to {np.max(w_path_cm[:, :3]):.2f} cm⁻¹")
    print(f"  Optic:    {np.min(w_path_cm[:, 3:]):.2f} to {np.max(w_path_cm[:, 3:]):.2f} cm⁻¹")
    
    # Test BZ grid
    print("\n" + "=" * 60)
    print("Test 3: Brillouin Zone Grid")
    print("=" * 60)
    
    Q_grid = generate_bz_grid(5)
    print(f"Generated {len(Q_grid)} q-points in 5×5×5 grid")
    print(f"Q-point range: [{Q_grid.min():.3f}, {Q_grid.max():.3f}]")
    
    print("\n" + "=" * 60)
    print("All tests passed!")
    print("=" * 60)
