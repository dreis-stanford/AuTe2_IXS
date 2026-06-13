"""
Interactive IXS analysis for single Q-points.

Material-generic analyzer (originally written for AuTe2, matching MATLAB
calcIXS_single_Q_AuTe2.m). AuTe2 is the default material; Silicon is
configured in single_q_analysis_si.py via the same class.
"""

import numpy as np
import sys
from . import config
from .constants import const
from .force_constants import ForceConstants
from .phonons import calc_Dq, calc_freq_eig, convert_frequencies
from .form_factors import calc_form_factor
from .ixs import calc_ixs, polarization_factor
from .sixcircle_interface import SixCircleInterface
from .aute2_structure import aute2_conv2prim_k, aute2_prim2conv_k
from .modulated_structure import ModulatedStructure


# ============================================================================
# Mode polarization classification (unified scheme)
# ============================================================================
# Strict L/T thresholds: meaningful for selection rules. On symmetry lines the
# (fixed) L-char numerics reach exact 0/1 to ~1e-6, so 0.98/0.02 is
# comfortably above numerical noise; tighten here if desired.
STRICT_L = 0.98   # L-char above this  -> 'L'    (longitudinal, selection-rule grade)
STRICT_T = 0.02   # L-char below this  -> 'T'    (transverse, selection-rule grade)
PRIMARY_L = 0.7   # L-char above this  -> 'M(L)' (mixed, primarily longitudinal)
PRIMARY_T = 0.3   # L-char below this  -> 'M(T)' (mixed, primarily transverse)
#                   anything else      -> 'M'    (mixed)


def _safe_pause(plt, fig, dt):
    """Pump the GUI event loop for `fig`, tolerating the user having closed
    the window early (figure already gone, or backend raises on a stale
    canvas) -- avoids crashing the REPL on plt.pause()/plt.close()."""
    try:
        if plt.fignum_exists(fig.number):
            plt.pause(dt)
    except Exception:
        pass


def classify_polarization(L):
    """
    Classify mode polarization from longitudinal character L in [0, 1].

    Returns 'L', 'T', 'M(L)', 'M(T)', 'M', or '--' (undefined, e.g. at
    reduced q = Gamma where L is NaN).
    """
    if np.isnan(L):
        return '--'
    if L > STRICT_L:
        return 'L'
    if L < STRICT_T:
        return 'T'
    if L > PRIMARY_L:
        return 'M(L)'
    if L < PRIMARY_T:
        return 'M(T)'
    return 'M'


class SingleQAnalyzer:
    """
    Analyze phonons and IXS at a single Q-point (material-generic).

    Defaults configure AuTe2 (monoclinic C2/m conv<->prim transforms, CDW
    satellites). Pass conv2prim/prim2conv and enable_satellites=False for
    other materials (see single_q_analysis_si.py for Silicon).
    """

    def __init__(self, xtal, Phi, masses, kT_THz,
                 conv2prim=None, prim2conv=None,
                 enable_satellites=True, material='AuTe2', sixc=None):
        """
        Initialize analyzer

        Parameters:
        -----------
        xtal : ForceConstants
            Crystal structure from force constants file
        Phi : ndarray
            Force constants in eV/Ų
        masses : ndarray
            Atomic masses in amu
        kT_THz : float
            Temperature in THz
        conv2prim, prim2conv : callable, optional
            Reciprocal-space coordinate transforms (default: AuTe2 C2/m)
        enable_satellites : bool
            Enable CDW satellite (H K L m) handling (default True; AuTe2)
        material : str
            Material name for printouts
        sixc : SixCircleInterface, optional
            If provided, used to compute the actual (tth, gam) scattering
            angles for each Q (via calculate_angles) so the IXS cross
            sections include the beam-polarization factor. If None, or if
            the reflection is inaccessible under the current frozen
            angles/limits, falls back to an in-plane (gam=0) elastic-Bragg
            estimate of tth from |Q|; for |Q| >= 4*pi/lambda this estimate
            saturates at tth=180 deg (polarization_factor=1, i.e. no
            correction).
        """
        self.xtal = xtal
        self.Phi = Phi
        self.masses = masses
        self.kT_THz = kT_THz
        self.material = material
        self.sixc = sixc
        self.conv2prim = conv2prim if conv2prim is not None else aute2_conv2prim_k
        self.prim2conv = prim2conv if prim2conv is not None else aute2_prim2conv_k
        self.mod_struct = ModulatedStructure() if enable_satellites else None

        # Per-atom element symbols and unique table labels (e.g. Au, Te1, Te2)
        self.symbols = [xtal.symbols[xtal.atom_type_map[i] - 1]
                        for i in range(xtal.nat)]
        counts = {s: self.symbols.count(s) for s in set(self.symbols)}
        seen = {}
        self.atom_labels = []
        for s in self.symbols:
            if counts[s] == 1:
                self.atom_labels.append(s)
            else:
                seen[s] = seen.get(s, 0) + 1
                self.atom_labels.append(f'{s}{seen[s]}')

        # Debye-Waller factor (see code/debye_waller.py)
        self.dw_U = None
        if config.DEBYE_WALLER_MODE == 'phonon':
            from .debye_waller import compute_U_tensors_phonon
            self.dw_U = compute_U_tensors_phonon(
                self.xtal, self.Phi, self.masses, self.kT_THz,
                n_mesh=config.DEBYE_WALLER_QMESH)
        elif config.DEBYE_WALLER_MODE == 'cif':
            if self.material == 'AuTe2':
                from .debye_waller import compute_U_tensors_cif
                self.dw_U = compute_U_tensors_cif(self.xtal)
            # else: the literature ADPs in compute_U_tensors_cif are AuTe2
            # (Reithmayer et al.) specific and don't apply to other
            # materials (e.g. Silicon) -- leave dw_U = None (no DW).
        elif config.DEBYE_WALLER_MODE != 'none':
            raise ValueError(
                f"Unknown DEBYE_WALLER_MODE: {config.DEBYE_WALLER_MODE!r}")

    def analyze(self, Q_input, coords='primitive', freq_unit='meV', print_results=True,
                print_detailed=False, tth_gam=None):
        """
        Analyze phonons and IXS at single Q-point

        Matches MATLAB calcIXS_single_Q_AuTe2.m

        Parameters:
        -----------
        Q_input : array-like (3,)
            Q-vector in specified coordinate system
        coords : str
            'primitive' or 'conventional'
        print_results : bool
            Print summary table
        print_detailed : bool
            Print detailed mode analysis
        threshold_L : float
            Longitudinal character threshold
        threshold_T : float
            Transverse character threshold
        tth_gam : tuple (tth, gam) in degrees, optional
            If given, use these directly for the beam-polarization factor
            instead of calling self.sixc.calculate_angles(Q_conv) -- avoids
            sixcircle's expensive nonlinear angle solve when the caller
            already knows (tth, gam) for this Q (e.g. analyze_array(), which
            gets them from forward kinematics via analyzer_array_offsets()).

        Returns:
        --------
        result : dict
            Dictionary with all analysis results
        """
        
        Q_input = np.asarray(Q_input)
        
        # Convert to primitive coordinates if needed
        if coords.lower() == 'conventional':
            Q_conv = Q_input
            Q_prim = aute2_conv2prim_k(Q_conv)
        else:
            Q_prim = Q_input
            Q_conv = aute2_prim2conv_k(Q_prim)
        
        result = {
            'Q_prim': Q_prim,
            'Q_conv': Q_conv
        }
        
        # Convert to Cartesian (2pi/A)
        Q_cart = Q_prim @ self.xtal.b_l
        Q_mag = np.linalg.norm(Q_cart)

        # For longitudinal character, use reduced q (phonon wavevector)
        result['Q_cart'] = Q_cart
        result['Q_mag'] = Q_mag
        
        # Find nearest reciprocal lattice vector
        # Search in Cartesian space to handle non-orthogonal reciprocal lattice
        Q_prim_rounded = np.round(Q_prim)
        Q_cart = Q_prim @ self.xtal.b_l
        
        # Search nearby integer combinations
        min_dist = np.inf
        best_G = Q_prim_rounded
        
        for dh in [-1, 0, 1]:
            for dk in [-1, 0, 1]:
                for dl in [-1, 0, 1]:
                    G_test = Q_prim_rounded + np.array([dh, dk, dl])
                    G_cart = G_test @ self.xtal.b_l
                    dist = np.linalg.norm(Q_cart - G_cart)
                    if dist < min_dist:
                        min_dist = dist
                        best_G = G_test
        
        G_prim = best_G
        Q_prim_reduced = Q_prim - G_prim
        
        Q_conv_G = aute2_prim2conv_k(G_prim)
        Q_conv_reduced = aute2_prim2conv_k(Q_prim_reduced)
        
        Q_cart_G = G_prim @ self.xtal.b_l
        Q_cart_reduced = Q_prim_reduced @ self.xtal.b_l
        Q_mag_reduced = np.linalg.norm(Q_cart_reduced)
        
        result['Q_reduced_prim'] = Q_prim_reduced
        result['Q_reduced_conv'] = Q_conv_reduced
        result['Q_reduced_cart'] = Q_cart_reduced
        result['G_prim'] = G_prim
        result['G_conv'] = Q_conv_G
        result['G_cart'] = Q_cart_G
        
        # Calculate dynamical matrix
        Q_for_calc = Q_prim_reduced.reshape(1, 3)
        Dq = calc_Dq(Q_for_calc, self.xtal.uvw, self.Phi, self.masses)
        
        # Solve for frequencies and eigenvectors
        w_raw, ev = calc_freq_eig(Dq, sort_modes=True, rotate_ev=True)
        
        # Convert frequencies
        w_THz = convert_frequencies(w_raw, 'THz').flatten()
        w_cm = convert_frequencies(w_raw, 'cm-1').flatten()
        w_meV = convert_frequencies(w_raw, 'meV').flatten()
        
        result['frequencies_THz'] = w_THz
        result['frequencies_cm'] = w_cm
        result['frequencies_meV'] = w_meV
        result['eigenvectors'] = ev
        
        nmodes = len(w_cm)
        
        # Calculate longitudinal/transverse character
        long_char = np.zeros(nmodes)
        Q_hat = None
        at_gamma = Q_mag_reduced < 1e-6
        
        if not at_gamma:
            Q_hat = Q_cart_reduced / Q_mag_reduced
            
            for imode in range(nmodes):
                ev_mode = ev[:, imode, 0]
                # DOF layout is atom-major [x1,y1,z1, x2,y2,z2, ...]
                # -> reshape to (nat, 3) so each row is one atom's (x,y,z)
                ev_reshaped = ev_mode.reshape(self.xtal.nat, 3)

                total_amp_sq = 0.0
                long_amp_sq = 0.0

                for iat in range(self.xtal.nat):
                    ev_atom = ev_reshaped[iat, :]

                    # Convert to physical displacement (undo mass weighting)
                    ev_atom_physical = ev_atom / np.sqrt(self.masses[iat])

                    # Sum |e_i|^2 (vdot conjugates; plain dot would give |sum e_i^2|)
                    atom_amp_sq = np.real(np.vdot(ev_atom_physical, ev_atom_physical))
                    total_amp_sq += atom_amp_sq

                    long_proj = np.dot(Q_hat, ev_atom_physical)
                    long_amp_sq += np.abs(long_proj)**2
                
                if total_amp_sq > 1e-12:
                    long_char[imode] = long_amp_sq / total_amp_sq
        else:
            long_char[:] = np.nan
        
        result['long_char'] = long_char
        result['Q_hat'] = Q_hat
        
        # Classify modes (unified scheme, see classify_polarization)
        pol_type = [classify_polarization(L) for L in long_char]
        result['pol_type'] = pol_type
        
        # Atomic participation
        atom_participation = np.zeros((nmodes, self.xtal.nat))
        
        for imode in range(nmodes):
            ev_mode = ev[:, imode, 0]
            # Atom-major layout -> (nat, 3), one row per atom
            ev_reshaped = ev_mode.reshape(self.xtal.nat, 3)

            physical_amp = np.zeros(self.xtal.nat)
            for iat in range(self.xtal.nat):
                ev_atom = ev_reshaped[iat, :]
                ev_atom_physical = ev_atom / np.sqrt(self.masses[iat])
                physical_amp[iat] = np.sum(np.abs(ev_atom_physical)**2)
            
            total_amp = np.sum(physical_amp)
            
            for iat in range(self.xtal.nat):
                atom_participation[imode, iat] = physical_amp[iat] / (total_amp + 1e-12)
        
        result['atom_participation'] = atom_participation
        
        # Calculate signed longitudinal components (Q*e for each atom)
        Q_cart = result['Q_cart']
        Q_mag = result['Q_mag']
        Q_hat = Q_cart_reduced / Q_mag_reduced if Q_mag_reduced > 1e-10 else np.zeros(3)
        
        longitudinal_signed = np.zeros((nmodes, self.xtal.nat))
        
        for imode in range(nmodes):
            ev_mode = ev[:, imode]
            total_ev_norm = np.linalg.norm(ev_mode)
            
            for iat in range(self.xtal.nat):
                e_atom = ev_mode[3*iat:3*iat+3]
                Q_dot_e = np.real(np.vdot(Q_hat, e_atom))
                longitudinal_signed[imode, iat] = 100 * Q_dot_e / (total_ev_norm + 1e-12)
        
        result['longitudinal_signed'] = longitudinal_signed
        
        # Calculate form factors (one per element, applied per atom)
        Q_sinThOverLambda = Q_mag / (4 * np.pi)

        # f(Q,E) = f0(Q) + f'(E) + i*f''(E) (xraylib, if available),
        # otherwise the real f0(Q) from the built-in Cromer-Mann coefficients.
        f_by_element = {s: calc_form_factor(Q_mag, s, config.ENERGY_KEV, scale=4*np.pi)
                        for s in set(self.symbols)}

        result['Q_sinThOverLambda'] = Q_sinThOverLambda
        result['form_factors'] = f_by_element

        # |Q| <= 4*pi/lambda is required for elastic scattering at this
        # wavelength to reach this Q at any angle; beyond that the
        # reflection is kinematically unreachable.
        result['q_reachable'] = Q_sinThOverLambda * config.WAVELENGTH <= 1.0

        # Effective per-atom form factors, including the Debye-Waller factor
        # T_k(Q) = exp(-1/2 Q.U_k.Q) if config.DEBYE_WALLER_MODE != 'none'.
        if self.dw_U is not None:
            from .debye_waller import debye_waller_factor
            dw_factors = np.array([debye_waller_factor(Q_cart, self.dw_U[k])
                                    for k in range(self.xtal.nat)])
        else:
            dw_factors = np.ones(self.xtal.nat)
        f_eff = np.array([f_by_element[self.symbols[k]] * dw_factors[k]
                           for k in range(self.xtal.nat)])
        result['debye_waller_factors'] = dw_factors
        result['form_factors_eff'] = f_eff

        # Prepare for IXS calculation
        fQ_matrix = f_eff.reshape(1, -1)

        # Calculate IXS
        Is, Ias, n, F, xs_info = calc_ixs(
            w_THz.reshape(1, -1), ev, Q_prim.reshape(1, 3),  # Use FULL Q_prim for IXS
            self.xtal.b_l, self.xtal.xs, fQ_matrix, self.masses, self.kT_THz,
            units='barn', per_steradian=True
        )
        
        # Beam-polarization factor P = 1 - sin^2(tth)*cos^2(gam) (horizontal
        # incident polarization, unanalyzed scattered polarization). Use the
        # actual (tth, gam) for this reflection from sixcircle when
        # available; otherwise fall back to the in-plane (gam=0) elastic
        # estimate tth = 2*asin(|Q|*lambda/(4*pi)) -- exact for |Q| within
        # reach of the Ewald sphere, and -> P=1 (no effect) as |Q| -> 4*pi/lambda.
        if tth_gam is not None:
            tth_pol, gam_pol = tth_gam
        else:
            angles = None
            if self.sixc is not None:
                try:
                    angles = self.sixc.calculate_angles(
                        (float(Q_conv[0]), float(Q_conv[1]), float(Q_conv[2])))
                except RuntimeError:
                    angles = None
            if angles is not None:
                tth_pol, gam_pol = angles['tth'], angles['gam']
            else:
                sin_half_tth = np.clip(Q_sinThOverLambda * config.WAVELENGTH, -1.0, 1.0)
                tth_pol = 2 * np.degrees(np.arcsin(sin_half_tth))
                gam_pol = 0.0
        pol_factor = polarization_factor(tth_pol, gam_pol)

        result['polarization_factor'] = pol_factor
        result['tth_pol'] = tth_pol
        result['gam_pol'] = gam_pol

        result['IXS_stokes'] = Is.flatten() * pol_factor
        result['IXS_antistokes'] = Ias.flatten() * pol_factor
        result['bose_factor'] = n.flatten()
        result['cross_section_info'] = xs_info
        
        # Print results
        if print_results:
            self._print_results(result, detailed=print_detailed, freq_unit=freq_unit)
        
        return result
    

    def analyze_array(self, Q_center, coords='conventional', freq_unit='meV', 
                     gap_h=60, gap_v=60, print_results=True):
        """
        Analyze phonons across analyzer array at BL43LXU
        
        Parameters:
        -----------
        Q_center : array-like (3,)
            Center Q position for analyzer array
        coords : str
            'primitive' or 'conventional'  
        freq_unit : str
            'meV', 'THz', or 'cm-1'
        gap_h, gap_v : float
            Horizontal and vertical analyzer gaps (microns)
        print_results : bool
            Print array table
            
        Returns:
        --------
        array_results : list of dict
            Results for each analyzer position
        """
        
        Q_center = np.asarray(Q_center)

        # Convert to conventional if needed
        if coords.lower() == 'primitive':
            Q_center_conv = aute2_prim2conv_k(Q_center)
        else:
            Q_center_conv = Q_center

        # Per-analyzer (dH, dK, dL) offsets for the arm centered on
        # Q_center_conv, computed from the real UB matrix and the BL43LXU
        # 7x4 analyzer-array geometry (see SixCircleInterface.analyzer_array_offsets).
        if self.sixc is None or self.sixc.sixc is None:
            raise RuntimeError(
                "analyze_array requires the sixcircle library to compute "
                "real per-analyzer Q positions (check config.SIXCIRCLE_PATH).")
        analyzers = self.sixc.analyzer_array_offsets(
            (float(Q_center_conv[0]), float(Q_center_conv[1]), float(Q_center_conv[2])))
        if analyzers is None:
            raise RuntimeError(
                f"Q={tuple(Q_center_conv)} is inaccessible under the current "
                f"frozen angles/limits, so no analyzer array can be calculated.")

        array_results = []
        
        # Suppress form factor messages during batch calculation
        import sys
        import io
        
        for name in sorted(analyzers.keys()):
            info = analyzers[name]
            dQ = np.array(info['dQ'])
            Q_analyzer = Q_center_conv + dQ

            # Analyze this Q point (suppress printing and form factor messages)
            old_stdout = sys.stdout
            sys.stdout = io.StringIO()
            try:
                result = self.analyze(Q_analyzer, coords='conventional',
                                     freq_unit=freq_unit, print_results=False,
                                     tth_gam=(info['tth'], info['gam']))
            finally:
                sys.stdout = old_stdout

            result['analyzer'] = name
            result['Q_offset'] = dQ
            result['dtth'] = info['dtth']
            result['dgam'] = info['dgam']
            result['xi'] = info['xi']
            result['zi'] = info['zi']
            array_results.append(result)
        
        if print_results:
            self._print_array_results(array_results, freq_unit)

        return array_results

    def analyze_grid(self, Q_conv_grid, tth_grid, gam_grid, freq_unit='meV'):
        """
        Vectorized phonon frequencies + Stokes IXS intensity over a grid of
        Q-points (e.g. a wider (dtth, dgam) region than the analyzer array
        spans, from SixCircleInterface.compute_angle_grid()).

        Unlike analyze(), this computes only frequencies and IXS_stokes (no
        L/T character, atom participation, etc.) and processes all grid
        points in a single batched calc_Dq/calc_freq_eig/calc_ixs call --
        the vectorization analyze_array() lacked when looping analyze() over
        many Q-points.

        Parameters
        ----------
        Q_conv_grid : ndarray (..., 3)
            Q-points in conventional coordinates, any leading shape.
        tth_grid, gam_grid : ndarray (...)
            Absolute (tth, gam) (degrees), same leading shape as
            Q_conv_grid, for the beam-polarization factor.
        freq_unit : str
            'meV', 'THz', or 'cm-1'

        Returns:
        --------
        dict with 'frequencies_<unit>' and 'IXS_stokes', each of shape
        (..., nmodes)
        """
        Q_conv_grid = np.asarray(Q_conv_grid)
        shape = Q_conv_grid.shape[:-1]
        Q_conv = Q_conv_grid.reshape(-1, 3)
        Q_prim = aute2_conv2prim_k(Q_conv)
        nQ = Q_prim.shape[0]

        Dq = calc_Dq(Q_prim, self.xtal.uvw, self.Phi, self.masses)
        w_raw, ev = calc_freq_eig(Dq, sort_modes=False, rotate_ev=False)
        w_THz = convert_frequencies(w_raw, 'THz')
        freq_key, freq_label = {'meV': ('frequencies_meV', 'meV'),
                                 'THz': ('frequencies_THz', 'THz'),
                                 'cm-1': ('frequencies_cm', 'cm-1')}[freq_unit]
        freq_out = convert_frequencies(w_raw, freq_label)
        nmodes = freq_out.shape[1]

        Q_cart = Q_prim @ self.xtal.b_l
        Q_mag = np.linalg.norm(Q_cart, axis=1)

        f_by_element = {s: calc_form_factor(Q_mag, s, config.ENERGY_KEV, scale=4 * np.pi)
                        for s in set(self.symbols)}

        if self.dw_U is not None:
            dw_factors = np.stack([
                np.exp(-0.5 * np.einsum('qi,ij,qj->q', Q_cart, self.dw_U[k], Q_cart))
                for k in range(self.xtal.nat)
            ], axis=1)  # (nQ, nat)
        else:
            dw_factors = np.ones((nQ, self.xtal.nat))

        f_eff = np.stack([f_by_element[self.symbols[k]]
                           for k in range(self.xtal.nat)], axis=1)  # (nQ, nat)
        fQ_matrix = f_eff * dw_factors

        Is, _, _, _, _ = calc_ixs(
            w_THz, ev, Q_prim, self.xtal.b_l, self.xtal.xs, fQ_matrix,
            self.masses, self.kT_THz, units='barn', per_steradian=True)

        pol_factor = polarization_factor(np.asarray(tth_grid).reshape(-1),
                                          np.asarray(gam_grid).reshape(-1))
        IXS_stokes = Is * pol_factor[:, np.newaxis]

        return {
            freq_key: freq_out.reshape(shape + (nmodes,)),
            'IXS_stokes': IXS_stokes.reshape(shape + (nmodes,)),
        }

    def analyze_array_grid(self, Q_center, coords='conventional', freq_unit='meV',
                            dtth_range=(-3.0, 3.0), dgam_range=(-3.0, 1.0),
                            n_tth=25, n_gam=25):
        """
        Phonon frequency/IXS maps over a (dtth, dgam) region, for the
        background of the array-command contour plots (see
        array_viz.plot_analyzer_array). Wider than the BL43LXU analyzer
        array itself spans (dtth in [-2.1, 2.1] deg, dgam in [-2.1, 0] deg).

        Parameters
        ----------
        Q_center : array-like (3,)
            Center Q position (same convention as analyze_array())
        coords : str
            'primitive' or 'conventional'
        freq_unit : str
            'meV', 'THz', or 'cm-1'
        dtth_range, dgam_range : (float, float)
            (min, max) angular offsets (degrees) from the array center to
            cover.
        n_tth, n_gam : int
            Grid resolution along dtth and dgam.

        Returns:
        --------
        grid : dict
            'dtth', 'dgam' (n_gam, n_tth): angular offsets from the array
            center (same convention as analyzer_array_offsets' dtth/dgam).
            'frequencies_<unit>', 'IXS_stokes' (n_gam, n_tth, nmodes).
        """
        Q_center = np.asarray(Q_center)
        if coords.lower() == 'primitive':
            Q_center_conv = aute2_prim2conv_k(Q_center)
        else:
            Q_center_conv = Q_center

        if self.sixc is None or self.sixc.sixc is None:
            raise RuntimeError(
                "analyze_array_grid requires the sixcircle library to compute "
                "real per-point Q positions (check config.SIXCIRCLE_PATH).")

        hkl = (float(Q_center_conv[0]), float(Q_center_conv[1]), float(Q_center_conv[2]))

        dtth = np.linspace(dtth_range[0], dtth_range[1], n_tth)
        dgam = np.linspace(dgam_range[0], dgam_range[1], n_gam)

        grid = self.sixc.compute_angle_grid(hkl, dtth, dgam)
        if grid is None:
            raise RuntimeError(
                f"Q={tuple(Q_center_conv)} is inaccessible under the current "
                f"frozen angles/limits, so no analyzer-array grid can be calculated.")

        result = self.analyze_grid(grid['Q_conv'], grid['tth'], grid['gam'], freq_unit=freq_unit)
        grid.update(result)
        return grid

    def _print_array_results(self, array_results, freq_unit='meV'):
        """Print compact analyzer array table: dQ offsets + freq/Pol pairs"""

        print('\n' + '='* 80)
        print(f'Analyzer Array Results - {len(array_results)} analyzers')
        print('Note: Analyzer (dH,dK,dL) offsets from Q_center computed from the real UB')
        print('      matrix and the BL43LXU 7x4 analyzer-array geometry.')
        print('='* 80)

        # Header
        freq_label = {'meV': 'meV', 'cm-1': 'cm^-1', 'THz': 'THz'}[freq_unit]

        print(f'{"#":>2s} {"dH":>5s} {"dK":>5s} {"dL":>5s}  '
              f'Freq ({freq_label}), modes 1-9 (low to high)')
        print('-'* 80)

        freq_key = {'meV': 'frequencies_meV',
                   'cm-1': 'frequencies_cm',
                   'THz': 'frequencies_THz'}[freq_unit]

        for result in array_results:
            num = result['analyzer'][1:]
            dQ = result['Q_offset']
            freqs = result[freq_key]

            freq_str = ' '.join(f'{f:5.1f}' for f in freqs)

            print(f'{num:>2s} {dQ[0]:5.2f} {dQ[1]:5.2f} {dQ[2]:5.2f}  '
                  f'{freq_str}')

        print('='*80 + '\n')

    @staticmethod
    def qe_projections(ev_mode, Q_hat, nat, noise_threshold=0.05):
        """
        Per-atom Q.epsilon projections for the mode table.

        Parameters:
        -----------
        ev_mode : ndarray (3*nat,) complex
            Mass-weighted eigenvector (polarization vector epsilon),
            atom-major layout [x1,y1,z1, x2,y2,z2, ...]
        Q_hat : ndarray (3,)
            Unit vector along full Q (Cartesian)
        nat : int
            Number of atoms in the cell
        noise_threshold : float
            Phases are reported as NaN where |Q.e| < noise_threshold
            (absolute scale; eigenvectors are normalized to 1, so the phase
            of a near-zero projection is numerical noise)

        Returns:
        --------
        magnitudes : ndarray (nat,)
            |Q_hat . epsilon_atom|
        phases_deg : ndarray (nat,)
            Phase relative to the atom with largest |Q.e|, in (-180, 180].
            NaN where the magnitude is below the noise threshold.

        Note: phases follow the dynamical-matrix convention (phase factor
        e^{2*pi*i q.R} over lattice vectors only), the same convention used
        by calc_ixs. They therefore include the physical phase advance of
        the wave across the atomic basis.
        """
        ev_reshaped = np.asarray(ev_mode).reshape(nat, 3)
        q_dot_e = np.array([np.vdot(Q_hat, ev_reshaped[iat, :])
                            for iat in range(nat)])

        magnitudes = np.abs(q_dot_e)
        ref_atom = int(np.argmax(magnitudes))
        phases_deg = np.degrees(np.angle(q_dot_e) - np.angle(q_dot_e[ref_atom]))
        phases_deg = np.mod(phases_deg + 180, 360) - 180
        phases_deg[magnitudes < noise_threshold] = np.nan
        return magnitudes, phases_deg

    def _print_results(self, result, detailed=False, freq_unit='meV'):
        """Print analysis results - compact format"""
        
        Q_conv = result['Q_conv']
        Q_prim = result['Q_prim']
        Q_cart = result['Q_cart']
        Q_mag = result['Q_mag']

        if not result['q_reachable']:
            Q_max = 4 * np.pi / config.WAVELENGTH
            print(f'\n|Q|={Q_mag:.2f} (2pi/A)  '
                  f'Q(conv): [{Q_conv[0]:.2f}, {Q_conv[1]:.2f}, {Q_conv[2]:.2f}] (r.l.u.)  '
                  f'Q(prim): [{Q_prim[0]:.2f}, {Q_prim[1]:.2f}, {Q_prim[2]:.2f}] (r.l.u.)')
            print(f'\n⚠ |Q|={Q_mag:.2f} A^-1 exceeds 4*pi/lambda = {Q_max:.2f} A^-1 '
                  f'(lambda={config.WAVELENGTH:.5f} A) -- unreachable by elastic '
                  f'scattering at this wavelength. Skipping IXS table.\n')
            return

        # Temperature and form factors
        kT_K = self.kT_THz * const.THz2meV / 1000 / 8.617333262e-5
        ff_str = '  '.join(f'f({s})={f:.1f}'
                           for s, f in sorted(result["form_factors"].items()))
        print(f'\nT={kT_K:.0f}K  '
              f'sin(θ)/λ={result["Q_sinThOverLambda"]:.4f} A^-1  {ff_str}')

        if self.dw_U is not None:
            dwf_str = '  '.join(f'{lbl}={d:.3f}' for lbl, d in
                                 zip(self.atom_labels, result['debye_waller_factors']))
            print(f'DWF(T={kT_K:.0f}K): {dwf_str}')

        # Calculate structure factor for elastic scattering
        # F(Q) = sum_atoms f_eff_atom(Q) * exp(2pii Q_prim * r_frac)
        F = 0.0 + 0.0j
        Q_prim_vec = Q_prim
        f_eff = result['form_factors_eff']

        for iat in range(self.xtal.nat):
            r_frac = self.xtal.xs[iat]  # Fractional coordinates
            # Sign convention matches calc_ixs (exp(-2j*pi*Q.r))
            phase = -2 * np.pi * np.dot(Q_prim_vec, r_frac)
            F += f_eff[iat] * np.exp(1j * phase)
        
        F_squared = np.abs(F)**2
        
        # Q info - 2 significant figures, two lines, just before table
        # Check if at Gamma
        Q_cart_reduced = result['Q_reduced_cart']
        Q_mag_reduced = np.linalg.norm(Q_cart_reduced)
        at_gamma = Q_mag_reduced < 1e-6
        
        gamma_note = "  (q at Gamma)" if at_gamma else ""
        
        print(f'\n|Q|={Q_mag:.2f} (2pi/A)  '
              f'Q(cart): [{Q_cart[0]:.2f}, {Q_cart[1]:.2f}, {Q_cart[2]:.2f}] (2pi/A){gamma_note}')
        print(f'Q(conv): [{Q_conv[0]:.2f}, {Q_conv[1]:.2f}, {Q_conv[2]:.2f}] (r.l.u.)  '
              f'Q(prim): [{Q_prim[0]:.2f}, {Q_prim[1]:.2f}, {Q_prim[2]:.2f}] (r.l.u.)')

        # Mode table with separators
        freq_label = {'meV': 'meV', 'cm-1': 'cm^-1', 'THz': 'THz'}[freq_unit]
        freq_data = {'meV': result['frequencies_meV'], 
                     'cm-1': result['frequencies_cm'],
                     'THz': result['frequencies_THz']}[freq_unit]
        
        atom_cols = '   '.join(f'{lbl}:Q*e,ph' for lbl in self.atom_labels)
        print('\n' + '='* 80)
        print(f'Mode Freq({freq_label:>3s})    {"IXS(S)":>8s}   {"IXS(AS)":>8s}  Pol  {atom_cols}')
        print('-'* 80)
        
        # Phonon modes
        w_meV = result['frequencies_meV']
        long_char = result['long_char']
        Is = result['IXS_stokes']
        Ias = result['IXS_antistokes']
        
        for i in range(len(w_meV)):
            # Determine polarization (unified scheme)
            pol = classify_polarization(long_char[i])
            
            # Format IXS values (0.00 to 99999.99)
            # Show --- for acoustic modes at or near Gamma point
            if at_gamma and freq_data[i] < 1.0:  # Acoustic mode at Gamma
                ixs_s_str = "---"
                ixs_as_str = "---"
            else:
                ixs_s = min(99999.99, max(0.0, Is[i]))
                ixs_as = min(99999.99, max(0.0, Ias[i]))
                ixs_s_str = f'{ixs_s:8.2f}'
                ixs_as_str = f'{ixs_as:8.2f}'
            
            # Calculate Q*e projections for each atom
            Q_hat_full = Q_cart / Q_mag if Q_mag > 1e-10 else np.zeros(3)

            ev_mode = result['eigenvectors'][:, i, 0]
            magnitudes, phases_deg = self.qe_projections(
                ev_mode, Q_hat_full, self.xtal.nat)

            # Format eigenvector info.
            # Blank everything at Q=0; blank phases where |Q.e| is numerical
            # noise, and for the degenerate acoustic modes at reduced q = Gamma
            # (eigenvectors are an arbitrary mixture there).
            acoustic_at_gamma = at_gamma and freq_data[i] < 1.0
            if Q_mag < 1e-10:
                ev_str = '    '.join(['     ---  '] * self.xtal.nat)
            else:
                ev_parts = []
                for iat in range(self.xtal.nat):
                    mag = magnitudes[iat]
                    if np.isnan(phases_deg[iat]) or acoustic_at_gamma:
                        ev_parts.append(f'{mag:4.2f},  --')
                    else:
                        ev_parts.append(f'{mag:4.2f},{phases_deg[iat]:+4.0f}')
                ev_str = '    '.join(ev_parts)
            
            print(f'{i+1:3d}  {freq_data[i]:7.2f}    {ixs_s_str:>8s}   {ixs_as_str:>8s}  {pol:>4s}    {ev_str}')
        
        # Elastic line with structure factor breakdown
        print('-'* 80)
        
        # Calculate Q*r for each atom (for elastic scattering)
        # Show e^(-iQ*r) even at Gamma (all = 1.0, phase = 0)
        if True:  # Always calculate
            elastic_parts = []
            # Show e^(-2pi*i Q*r) for each atom (amplitude and phase)
            # Sign convention matches calc_ixs (exp(-2j*pi*Q.r)) so phonon-row
            # and elastic-row phases are directly comparable.
            exponentials = np.zeros(self.xtal.nat, dtype=complex)

            for iat in range(self.xtal.nat):
                r_frac = self.xtal.xs[iat]
                phase = -2 * np.pi * np.dot(Q_prim, r_frac)
                exponentials[iat] = np.exp(1j * phase)
            
            # Find reference (first atom, or could use largest)
            ref_phase = np.angle(exponentials[0])
            
            for iat in range(self.xtal.nat):
                mag = np.abs(exponentials[iat])  # Always 1.0
                phase_rel = np.angle(exponentials[iat]) - ref_phase
                phase_deg = np.degrees(phase_rel)
                phase_deg = np.mod(phase_deg + 180, 360) - 180
                elastic_parts.append(f'{mag:3.1f},{phase_deg:+4.0f}')
            
            elastic_str = '    '.join(elastic_parts)
        
        print(f'  0      0.00    |F|^2 = {F_squared:8.0f}       El.  {elastic_str}')
        print('='* 80)
        
        #notes
        print('\nNotes: IXS units: barn/(unit cell*sr);  DFT uses unmodulated structure (no CDW).')
        print('\n' + '='* 80 + '\n')


# Material configurations for interactive_mode
MATERIALS = {
    'AuTe2': dict(
        fc_file="data/AuTe_2_m.fc",
        conv2prim=None, prim2conv=None,      # defaults = AuTe2 C2/m
        satellites=True,                      # CDW satellites (H K L m)
        sixcircle=True,                       # diffractometer configured for AuTe2
        name='AuTe2',
    ),
    'Si': dict(
        fc_file="data/Test__Silicon_dispersion/Qgrid_888/Cg.fc",
        conv2prim=None, prim2conv=None,       # set below (FCC), import here
        satellites=False,
        sixcircle=False,                      # sixcircle lattice is AuTe2-specific
        name='Silicon',
    ),
}


def interactive_mode(material='AuTe2'):
    """
    Interactive Q-point analysis mode (material-generic)
    Matches MATLAB interactive script
    """
    mat = MATERIALS[material]

    if material == 'Si':
        from .fcc_structure import fcc_conv2prim_k, fcc_prim2conv_k
        mat = dict(mat, conv2prim=fcc_conv2prim_k, prim2conv=fcc_prim2conv_k)

    print("=" * 80)
    print(f"  Interactive IXS Analysis for {mat['name']}")
    print("=" * 80)

    # Load force constants
    fc_file = mat['fc_file']
    print(f"\nLoading force constants from {fc_file}...")

    xtal = ForceConstants(fc_file)
    Phi = xtal.convert_to_eV_per_Angstrom2()

    # Get masses
    masses = np.array([xtal.masses[xtal.atom_type_map[i]-1]
                       for i in range(xtal.nat)])

    # Temperature (set config.TEMPERATURE to change)
    kT_THz = config.TEMPERATURE * const.kb * const.eV2THz

    print(f"Temperature: {config.TEMPERATURE:.1f} K ({kT_THz:.2f} THz)\n")

    # Initialize sixcircle interface once (reuse it for all angle calculations,
    # and so the analyzer can compute the beam-polarization factor per Q)
    sixc = None
    if mat['sixcircle']:
        try:
            try:
                from .sixcircle_interface import SixCircleInterface
            except ImportError:
                from sixcircle_interface import SixCircleInterface
            sixc = SixCircleInterface()
            status = "calculation only" if sixc.sixc is not None else "unavailable"
            print(f"✓ Sixcircle interface loaded ({status})")
        except Exception as e:
            print(f"⚠ Sixcircle not available: {e}")

    # Create analyzer
    analyzer = SingleQAnalyzer(xtal, Phi, masses, kT_THz,
                               conv2prim=mat['conv2prim'],
                               prim2conv=mat['prim2conv'],
                               enable_satellites=mat['satellites'],
                               material=mat['name'],
                               sixc=sixc)

    if analyzer.dw_U is not None:
        u_iso_str = '  '.join(
            f'{lbl}={np.trace(U)/3:.4f}'
            for lbl, U in zip(analyzer.atom_labels, analyzer.dw_U))
        print(f"Debye-Waller (phonon, U_iso, A^2): {u_iso_str}\n")

    # Default coordinate system
    coord_system = 'conventional'
    freq_unit = 'meV'

    print("Instructions:")
    print("  - Enter Q vector as three numbers (e.g., 0.5 0 0)")
    if mat['satellites']:
        print("  - For satellites, add order m: H K L m (e.g., 1 0 0 1 for first satellite)")
    print("  - Press Enter on empty line to quit")
    print("  - Type 'conv', 'prim', or 'cart' to change coordinate system")
    print("  - Type 'meV', 'THz', or 'invcm' to change frequency units")
    print("  - Type 'array' to calculate full analyzer array at last Q")
    if mat['sixcircle']:
        print("  - Type 'angles' to calculate diffractometer angles for last Q")
        print("  - Type 'viz' to show a 3D view of the scattering geometry for last Q")
        print("  - Type 'freeze' / 'limits' to view or change the angle mode and limits")
    print(f"  - Current system: {coord_system}\n")
    
    current_q = None
    
    while True:
        print('-' * 80)
        user_input = input(f'Enter Q ({coord_system}): ').strip()
        
        # Check if user wants to quit
        if not user_input:
            print("\nExiting...")
            break
        
        # Check for command keywords first
        if user_input.lower() == 'prim':
            coord_system = 'primitive'
            print("  → Switched to primitive coordinates\n")
            continue
        elif user_input.lower() == 'conv':
            coord_system = 'conventional'
            print("  → Switched to conventional coordinates\n")
            continue
        elif user_input.lower() == 'cart':
            coord_system = 'cartesian'
            print("  → Switched to Cartesian coordinates\n")
            continue
        elif user_input.lower() == 'mev':
            freq_unit = 'meV'
            print("  → Switched to meV\n")
            continue
        elif user_input.lower() in ['cm', 'cm-1', 'invcm']:
            freq_unit = 'cm-1'
            print("  → Switched to cm^-1\n")
            continue
        elif user_input.lower() == 'thz':
            freq_unit = 'THz'
            print(f"  → Switched to THz\n")
            continue
        elif user_input.lower() == 'array':
            if current_q is None:
                print("\n⚠ Enter a Q point first\n")
                continue
            fig = None
            try:
                import matplotlib.pyplot as plt
                Q_conv = analyzer.prim2conv(current_q) if coord_system == 'primitive' else current_q
                array_results = analyzer.analyze_array(Q_conv, coords='conventional', freq_unit=freq_unit)
                analyzer_qs = {str(int(r['analyzer'][1:])): r['Q_conv'] for r in array_results}
                grid = analyzer.analyze_array_grid(Q_conv, coords='conventional', freq_unit=freq_unit)

                from .array_viz import plot_analyzer_array
                Q_label = f"({Q_conv[0]:.3f}, {Q_conv[1]:.3f}, {Q_conv[2]:.3f})"
                fig = plot_analyzer_array(array_results, freq_unit=freq_unit, Q_label=Q_label, grid=grid)
                plt.show(block=False)
                _safe_pause(plt, fig, 0.1)

                print("Enter analyzer # for its IXS table (blank: redisplay table, 'q': back)\n")
                while True:
                    sel = input("Analyzer #: ").strip().lower()
                    _safe_pause(plt, fig, 0.05)
                    if sel in ('q', 'quit', 'b', 'back'):
                        break
                    if not sel:
                        analyzer._print_array_results(array_results, freq_unit)
                        continue
                    try:
                        num = str(int(sel))
                    except ValueError:
                        num = sel
                    if num not in analyzer_qs:
                        print(f"\n⚠ Unknown analyzer '{sel}'. "
                              f"Available: {', '.join(sorted(analyzer_qs, key=int))}\n")
                        continue
                    analyzer.analyze(analyzer_qs[num], coords='conventional', freq_unit=freq_unit)
            except Exception as e:
                print(f"\n✗ {e}\n")
                import traceback
                traceback.print_exc()
            finally:
                if fig is not None:
                    try:
                        import matplotlib.pyplot as plt
                        if plt.fignum_exists(fig.number):
                            plt.close(fig)
                    except Exception:
                        pass
            continue

        elif user_input.lower() == 'angles':
            if current_q is None:
                print("\n⚠ Enter a Q point first\n")
                continue
            if sixc is None:
                print("\n⚠ Sixcircle not available\n")
                continue
            try:
                Q_conv = analyzer.prim2conv(current_q) if coord_system == 'primitive' else current_q
                # sixcircle doesn't like numpy.float64
                hkl = (float(Q_conv[0]), float(Q_conv[1]), float(Q_conv[2]))
                print("\n" + "=" * 80)
                print("Diffractometer Angles (calculation only — no motion)")
                print("=" * 80)
                print(f"Q (conv): [{Q_conv[0]:.4f}, {Q_conv[1]:.4f}, {Q_conv[2]:.4f}]")
                print(sixc.describe_constraints(include_limits=False))
                a = sixc.calculate_angles(hkl)
                out_of_limits = False
                if a is None:
                    a = sixc.calculate_angles(hkl, ignore_limits=True)
                    out_of_limits = a is not None
                if a is None:
                    print("\n  ⚠ Inaccessible: no solution exists even with "
                          "limits relaxed to ±180°.")
                else:
                    if out_of_limits:
                        violations = [
                            f"{name}={a[name]:.3f} (limit [{lo:g}, {hi:g}])"
                            for name, (lo, hi) in config.ANGLE_LIMITS.items()
                            if not (lo <= a[name] <= hi)]
                        print("\n  ⚠ Outside current angle limits"
                              + (": " + ", ".join(violations) if violations else "."))
                        print("    Try 'limits' to widen them or 'freeze' to change mode.")
                    print("\n  " + "  ".join(
                        f"{k}={a[k]:8.3f}" for k in
                        ('th', 'chi', 'phi', 'mu')))
                    print("  " + "  ".join(
                        f"{k}={a[k]:8.3f}" for k in ('tth', 'gam')))
                    print("  " + "  ".join(
                        f"{k}={a[k]:8.3f}" for k in ('alpha', 'beta', 'omega')))
                    print(f"\n  SPEC:  mv tth {a['tth']:.4f} th {a['th']:.4f} "
                          f"chi {a['chi']:.4f} phi {a['phi']:.4f} "
                          f"mu {a['mu']:.4f} gam {a['gam']:.4f}")
                print("=" * 80 + "\n")
            except Exception as e:
                print(f"\n✗ {e}\n")
                import traceback
                traceback.print_exc()
            continue

        elif user_input.lower() in ('viz', 'geo', 'geometry'):
            if current_q is None:
                print("\n⚠ Enter a Q point first\n")
                continue
            if sixc is None:
                print("\n⚠ Sixcircle not available\n")
                continue
            try:
                Q_conv = analyzer.prim2conv(current_q) if coord_system == 'primitive' else current_q
                hkl = (float(Q_conv[0]), float(Q_conv[1]), float(Q_conv[2]))
                a = sixc.calculate_angles(hkl)
                if a is None:
                    print("\n  ⚠ Inaccessible within current angle limits — "
                          "nothing to draw.\n")
                    continue
                save_path = f"scattering_geometry_{hkl[0]:g}_{hkl[1]:g}_{hkl[2]:g}.png"
                # Launch the 3D view in its own process: it gets its own GUI
                # event loop (fully rotatable) and the window can be closed any
                # time without affecting this REPL, which continues immediately.
                import subprocess
                angles_arg = ",".join(
                    f"{a[k]:.6f}" for k in ('tth', 'th', 'chi', 'phi', 'mu', 'gam'))
                cmd = [sys.executable, '-m', 'code.geometry_viz',
                       '--angles', angles_arg,
                       '--hkl', f"{hkl[0]:g},{hkl[1]:g},{hkl[2]:g}",
                       '--save', save_path]
                ub = sixc.get_UB()
                if ub is not None:
                    cmd += ['--ub', ",".join(f"{x:.10g}" for x in ub.flatten())]
                try:
                    subprocess.Popen(cmd, stdout=subprocess.DEVNULL,
                                     stderr=subprocess.DEVNULL)
                    print(f"\n  Opened 3D geometry window (rotatable). "
                          f"Saved to {save_path}.\n")
                except Exception as e:
                    print(f"\n✗ Could not open geometry window: {e}\n")
            except Exception as e:
                print(f"\n✗ {e}\n")
                import traceback
                traceback.print_exc()
            continue

        elif user_input.lower() == 'freeze':
            # Show current frozen mode / values
            if sixc is None:
                print("\n⚠ Sixcircle not available\n")
                continue
            print("\n" + sixc.describe_constraints())
            print("\n  Usage: 'freeze 456'        set which angles are frozen")
            print("         'freeze mu=-0.17'   set a frozen angle's value")
            print("         codes: tth0 th1 chi2 phi3 mu4 gam5 omega6 azimuth7 alpha8 beta9\n")
            continue

        elif user_input.lower().startswith('freeze '):
            if sixc is None:
                print("\n⚠ Sixcircle not available\n")
                continue
            try:
                tokens = user_input.split()[1:]
                values = {}
                for tok in tokens:
                    if '=' in tok:
                        name, val = tok.split('=', 1)
                        values[name.strip()] = float(val)
                    elif tok.isdigit():
                        sixc.set_frozen(tok)
                    else:
                        raise ValueError(f"Cannot parse '{tok}'")
                if values:
                    sixc.set_frozen_values(**values)
                print("\n" + sixc.describe_constraints() + "\n")
            except Exception as e:
                print(f"\n✗ {e}\n")
            continue

        elif user_input.lower() == 'limits':
            if sixc is None:
                print("\n⚠ Sixcircle not available\n")
                continue
            print("\n" + sixc.describe_constraints())
            print("\n  Usage: 'limits tth 0 60'   set an angle's [lower, upper]")
            print("         'limits reset'       restore all to ±180°\n")
            continue

        elif user_input.lower().startswith('limits '):
            if sixc is None:
                print("\n⚠ Sixcircle not available\n")
                continue
            try:
                tokens = user_input.split()[1:]
                if len(tokens) == 1 and tokens[0].lower() == 'reset':
                    sixc.reset_limits()
                elif len(tokens) == 3:
                    angle, lo, hi = tokens[0], float(tokens[1]), float(tokens[2])
                    sixc.set_limit(angle, lo, hi)
                else:
                    raise ValueError("Usage: 'limits <angle> <lower> <upper>' or 'limits reset'")
                print("\n" + sixc.describe_constraints() + "\n")
            except Exception as e:
                print(f"\n✗ {e}\n")
            continue
        
        elif user_input.lower() == 'move':
            # This tool is calculation-only; it never drives the diffractometer.
            print("\n  This is a planning tool — it does not move the "
                  "diffractometer.")
            print("  Use 'angles' to get the positions (and a ready-to-paste "
                  "SPEC command)")
            print("  to run at the beamline.\n")
            continue
        elif user_input.lower().startswith('sixc '):
            # Pass command to sixcircle
            cmd = user_input[5:].strip()  # Remove 'sixc ' prefix
            if not cmd:
                print("\n⚠ Usage: sixc <command>\n")
                print("  Examples: sixc wh, sixc pa, sixc or_check\n")
                continue
            if sixc is None:
                print(f"\n⚠ Sixcircle not configured for {analyzer.material}\n")
                continue
            try:
                if sixc.simulation_mode:
                    print("\n⚠ Simulation mode - sixcircle commands not available\n")
                    continue
                if not sixc.sixc:
                    print("\n⚠ Sixcircle not loaded\n")
                    continue
                print("")
                # Import whitelist from config
                from config import SIXCIRCLE_ALLOWED_COMMANDS, SIXCIRCLE_COMMAND_HELP
                
                # Parse command
                cmd_clean = cmd.strip()
                cmd_name = cmd_clean.split('(')[0].strip()
                
                # Check if command is allowed
                if cmd_name not in SIXCIRCLE_ALLOWED_COMMANDS:
                    print(f"⚠ Command '{cmd_name}' not allowed")
                    print(SIXCIRCLE_COMMAND_HELP)
                    print("")
                    continue
                
                # Execute whitelisted command
                method_name, needs_args = SIXCIRCLE_ALLOWED_COMMANDS[cmd_name]
                method = getattr(sixc.sixc, method_name, None)
                
                if method is None:
                    print(f"⚠ Method '{method_name}' not found in sixcircle")
                    print("")
                    continue
                
                # Call method
                if needs_args:
                    print(f"⚠ Commands with arguments not yet implemented")
                    print(f"  '{cmd_name}' requires arguments")
                else:
                    result = method() if callable(method) else method
                    if result is not None:
                        print(result)
            except Exception as e:
                print(f"\n✗ Error: {e}\n")
                print("  Try: sixc wh() or sixc pa()\n")
            continue
        
        try:
            parts = user_input.split()
            
            if len(parts) == 3:
                # Standard Q input: H K L
                Q = np.array([float(x) for x in parts])
                satellite_order = 0
            elif len(parts) == 4:
                # Satellite Q input: H K L m
                if analyzer.mod_struct is None:
                    print(f'⚠ Satellites not configured for {analyzer.material}\n')
                    continue
                if coord_system != 'conventional':
                    # q_mod is defined in conventional r.l.u.
                    print("⚠ Satellite input (H K L m) requires 'conv' coordinates\n")
                    continue
                Q_main = np.array([float(x) for x in parts[:3]])
                satellite_order = int(parts[3])
                # Calculate satellite position
                q_mod = np.array(analyzer.mod_struct.q_mod)
                Q = Q_main + satellite_order * q_mod
                print(f'  Main reflection: [{Q_main[0]:.4f}, {Q_main[1]:.4f}, {Q_main[2]:.4f}]')
                print(f'  Satellite order m={satellite_order}')
                print(f'  Modulation vector: [{q_mod[0]:.4f}, {q_mod[1]:.4f}, {q_mod[2]:.4f}]')
                print(f'  Satellite Q = Q_main + {satellite_order}*q_mod')
                print(f'             = [{Q[0]:.4f}, {Q[1]:.4f}, {Q[2]:.4f}]\n')
            else:
                print('⚠ Error: Enter 3 numbers (H K L) or 4 numbers (H K L m) for satellites\n')
                continue
            
            # Handle Cartesian input
            if coord_system == 'cartesian':
                Q_prim = Q @ np.linalg.inv(xtal.b_l)
                print(f'Converted Cartesian Q = [{Q[0]:.4f}, {Q[1]:.4f}, {Q[2]:.4f}] (2pi/A)')
                print(f'       to Primitive q = [{Q_prim[0]:.4f}, {Q_prim[1]:.4f}, {Q_prim[2]:.4f}] (r.l.u.)\n')
                Q_input = Q_prim
                input_coords = 'primitive'
            else:
                Q_input = Q
                input_coords = coord_system
            
            # Analyze
            current_q = Q_input
            result = analyzer.analyze(Q_input, coords=input_coords,
                                     freq_unit=freq_unit,
                                     print_results=True, print_detailed=False)
            
        except ValueError as e:
            print(f'⚠ Error parsing input: {e}\n')
        except Exception as e:
            print(f'⚠ Error: {e}\n')
            import traceback
            traceback.print_exc()
    
    print("\n" + "=" * 80)
    print("Session ended.")
    print("=" * 80)


if __name__ == "__main__":
    import os
    sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
    
    # Test with single Q-point
    print("=" * 80)
    print("Testing Single Q-Point Analyzer")
    print("=" * 80)
    
    # Load data
    fc_file = "data/AuTe_2_m.fc"
    xtal = ForceConstants(fc_file)
    Phi = xtal.convert_to_eV_per_Angstrom2()
    
    masses = np.array([xtal.masses[xtal.atom_type_map[i]-1] 
                       for i in range(xtal.nat)])
    
    # Temperature (set config.TEMPERATURE to change)
    kT_THz = config.TEMPERATURE * const.kb * const.eV2THz

    # Create analyzer
    analyzer = SingleQAnalyzer(xtal, Phi, masses, kT_THz)
    
    # Test analysis at Q = (0.5, 0, 0) primitive
    print("\nTest: Q = (0.5, 0, 0) primitive")
    result = analyzer.analyze([0.5, 0.0, 0.0], coords='primitive')
    
    print("\n" + "=" * 80)
    print("Test complete! To run interactive mode:")
    print("  from single_q_analysis import interactive_mode")
    print("  interactive_mode()")
    print("=" * 80)
