"""
Debye-Waller (atomic displacement) factors.

Mode 'phonon': anisotropic mean-square-displacement tensors U_k (Angstrom^2)
from a Brillouin-zone sum over the same force constants used for the
dispersion/IXS calculations, reusing the calc_Dq/calc_freq_eig/
generate_bz_grid pipeline (the same pipeline plot_dispersion.py uses for the
phonon DOS, just keeping the eigenvectors that the DOS calculation discards).
"""

import numpy as np
from .constants import const
from .phonons import calc_Dq, calc_freq_eig, convert_frequencies, generate_bz_grid

# <u^2> [Angstrom^2] = (hbar/(2*M*omega)) * (2n+1), with M in kg and omega in
# rad/s. Combining hbar (J.s), amu->kg, THz->rad/s (the 2*pi), and m^2->A^2
# (1e20) into a single prefactor (~0.505 A^2.THz) so that, with e_k already
# mass-normalized (e_k = ev_k / sqrt(M_k)):
#   <u_k,ab> = _DW_PREFACTOR/N_q * sum_{q,j} (2n+1)/f_THz * Re[e_k,a e_k,b*]
_DW_PREFACTOR = const.hbar * 1e8 / (4 * np.pi * const.amu)


def compute_U_tensors_phonon(xtal, Phi, masses, kT_THz, n_mesh=15):
    """
    Per-atom mean-square-displacement tensors U_k (Angstrom^2), via a
    Monkhorst-Pack BZ sum over the dynamical matrix at temperature kT_THz.

    Parameters:
    -----------
    xtal : ForceConstants
    Phi : ndarray
        Force constants in eV/Angstrom^2
    masses : ndarray (nat,)
        Atomic masses in amu
    kT_THz : float
        Temperature in THz
    n_mesh : int
        Monkhorst-Pack grid size (n_mesh x n_mesh x n_mesh), same
        generate_bz_grid() used for the DOS in plot_dispersion.py.

    Returns:
    --------
    U : ndarray (nat, 3, 3)
        Symmetric mean-square-displacement tensors, in the same Cartesian
        frame as xtal.b_l (and Q_cart = Q_rlu @ b_l).
    """
    masses = np.asarray(masses).flatten()
    nat = len(masses)

    q_grid = generate_bz_grid(n_mesh)
    nQ = q_grid.shape[0]
    D = calc_Dq(q_grid, xtal.uvw, Phi, masses)
    w_raw, ev = calc_freq_eig(D, sort_modes=False, rotate_ev=False)
    w_THz = convert_frequencies(w_raw, 'THz')  # (nQ, 3*nat)

    # Same acoustic-divergence floor as calc_ixs.
    min_w = 0.1  # THz
    w_safe = np.abs(w_THz)
    w_safe[w_safe < min_w] = min_w

    n_bose = 1.0 / (np.exp(w_safe / kT_THz) - 1.0)
    weight = (2 * n_bose + 1) / w_safe  # (nQ, 3*nat)

    U = np.zeros((nat, 3, 3))
    for k in range(nat):
        e_k = ev[3*k:3*k+3, :, :] / np.sqrt(masses[k])  # (3, 3*nat, nQ)
        U[k] = _DW_PREFACTOR / nQ * np.real(
            np.einsum('anq,bnq,nq->ab', e_k, np.conj(e_k), weight.T))

    return U


def debye_waller_factor(Q_cart, U):
    """
    Anisotropic Debye-Waller factor T(Q) = exp(-1/2 Q . U . Q), with Q_cart
    (Angstrom^-1, 2*pi convention) and U (Angstrom^2) in the same Cartesian
    frame. Reduces to exp(-8*pi^2 * U * sin^2(theta)/lambda^2) for isotropic
    U = U * I.
    """
    Q_cart = np.asarray(Q_cart)
    return np.exp(-0.5 * Q_cart @ U @ Q_cart)


# Literature ADPs for AuTe2 (calaverite), conventional C2/m cell, from
# Reithmayer et al., Acta Crystallographica B49, 6 (1993) -- same source as
# config.LATTICE_PARAMS. U_ij in the CIF/IUCr convention (Angstrom^2).
_CIF_CELL = dict(a=7.189, b=4.407, c=5.069, beta=89.96)  # Angstrom, degrees

_U_CIF_AU = np.array([
    [0.01520, 0.0,     0.00200],
    [0.0,     0.02050, 0.0],
    [0.00200, 0.0,     0.01390]])  # Au at (0, 0, 0)

_U_CIF_TE = np.array([
    [0.01670, 0.0,     0.00170],
    [0.0,     0.01030, 0.0],
    [0.00170, 0.0,     0.01040]])  # Te at (0.6884, 0, 0.2878)


def _conv_cell_vectors(a, b, c, beta):
    """
    Conventional C2/m cell as real-space row vectors (Angstrom), in the
    standard b-unique monoclinic convention (b || y, a and c in the xz-plane,
    tilted symmetrically about x/z by (90-beta)/2 so the angle between them
    is beta). This is the same convention T_PRIM_TO_CONV @ xtal.a_l produces
    for the DFT cell (alpha = gamma = 90 deg exactly), so the result lands in
    xtal's Cartesian frame to within <0.2 deg (negligible).
    """
    delta = np.radians((90.0 - beta) / 2)
    a_vec = a * np.array([np.cos(delta), 0.0, np.sin(delta)])
    b_vec = b * np.array([0.0, 1.0, 0.0])
    c_vec = c * np.array([np.sin(delta), 0.0, np.cos(delta)])
    return np.array([a_vec, b_vec, c_vec])


def _reciprocal_magnitudes(A):
    """|a*_i| (1/Angstrom, no 2*pi) for real-space cell rows A."""
    a_vec, b_vec, c_vec = A
    V = np.dot(a_vec, np.cross(b_vec, c_vec))
    a_star = np.array([
        np.linalg.norm(np.cross(b_vec, c_vec)),
        np.linalg.norm(np.cross(c_vec, a_vec)),
        np.linalg.norm(np.cross(a_vec, b_vec))]) / V
    return np.abs(a_star)


def _u_cif_to_cart(U_cif, A, D):
    """U_cif (CIF convention) -> Cartesian U tensor (Angstrom^2)."""
    return A.T @ D @ U_cif @ D @ A


def compute_U_tensors_cif(xtal):
    """
    Per-atom mean-square-displacement tensors U_k (Angstrom^2), from the
    fixed literature ADPs of Reithmayer et al., Acta Cryst. B49, 6 (1993).

    Material-specific (AuTe2 only): assumes the primitive-cell atom ordering
    (Au, Te1, Te2) of xtal.symbols/xtal.atom_type_map. Te1 and Te2 are the two
    C2/m-symmetry-equivalent Te sites; since the CIF's Te tensor has
    U12 = U23 = 0 (required by the m_y mirror, as both Au and Te sit at
    y = 0), the symmetry operation relating Te1 and Te2 (R = diag(-1,1,-1))
    leaves U unchanged, so both get the identical Cartesian U tensor.

    Returns:
    --------
    U : ndarray (3, 3, 3)
        Symmetric mean-square-displacement tensors, in the same Cartesian
        frame as xtal.b_l (and Q_cart = Q_rlu @ b_l).
    """
    if xtal.nat != 3:
        raise ValueError(
            f"compute_U_tensors_cif assumes the AuTe2 (Au, Te, Te) primitive "
            f"cell (nat=3), got nat={xtal.nat}")

    A = _conv_cell_vectors(**_CIF_CELL)
    D = np.diag(_reciprocal_magnitudes(A))

    U_Au = _u_cif_to_cart(_U_CIF_AU, A, D)
    U_Te = _u_cif_to_cart(_U_CIF_TE, A, D)

    return np.array([U_Au, U_Te, U_Te])
