"""
Interactive IXS analysis for single Q-points in AuTe2
Matches MATLAB calcIXS_single_Q_AuTe2.m and interactive script
"""

import numpy as np
import sys
from .constants import const
from .force_constants import ForceConstants
from .phonons import calc_Dq, calc_freq_eig, convert_frequencies
from .form_factors import CalcAtomicfQ
from .ixs import calc_ixs
from .sixcircle_interface import SixCircleInterface
from .aute2_structure import aute2_conv2prim_k, aute2_prim2conv_k
from .modulated_structure import ModulatedStructure


class SingleQAnalyzer:
    """
    Analyze phonons and IXS at a single Q-point
    Matches MATLAB calcIXS_single_Q_AuTe2.m
    """
    
    def __init__(self, xtal, Phi, masses, kT_THz):
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
        """
        self.xtal = xtal
        self.Phi = Phi
        self.masses = masses
        self.kT_THz = kT_THz
        self.mod_struct = ModulatedStructure()
        
    def analyze(self, Q_input, coords='primitive', freq_unit='meV', print_results=True, 
                print_detailed=False, threshold_L=0.7, threshold_T=0.3):
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
                ev_reshaped = ev_mode.reshape(3, self.xtal.nat)
                
                total_amp_sq = 0.0
                long_amp_sq = 0.0
                
                for iat in range(self.xtal.nat):
                    ev_atom = ev_reshaped[:, iat]
                    
                    # Convert to physical displacement (undo mass weighting)
                    ev_atom_physical = ev_atom / np.sqrt(self.masses[iat])
                    
                    atom_amp_sq = np.abs(np.dot(ev_atom_physical, ev_atom_physical))
                    total_amp_sq += atom_amp_sq
                    
                    long_proj = np.dot(Q_hat, ev_atom_physical)
                    long_amp_sq += np.abs(long_proj)**2
                
                if total_amp_sq > 1e-12:
                    long_char[imode] = long_amp_sq / total_amp_sq
        else:
            long_char[:] = np.nan
        
        result['long_char'] = long_char
        result['Q_hat'] = Q_hat
        
        # Classify modes
        pol_type = []
        for L in long_char:
            if np.isnan(L):
                pol_type.append('Undefined')
            elif L < threshold_T:
                pol_type.append('Transverse')
            elif L > threshold_L:
                pol_type.append('Longitudinal')
            else:
                pol_type.append('Mixed')
        
        result['pol_type'] = pol_type
        
        # Atomic participation
        atom_participation = np.zeros((nmodes, self.xtal.nat))
        
        for imode in range(nmodes):
            ev_mode = ev[:, imode, 0]
            ev_reshaped = ev_mode.reshape(3, self.xtal.nat)
            
            physical_amp = np.zeros(self.xtal.nat)
            for iat in range(self.xtal.nat):
                ev_atom = ev_reshaped[:, iat]
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
                longitudinal_signed[imode, iat] = 80 * Q_dot_e / (total_ev_norm + 1e-12)
        
        result['longitudinal_signed'] = longitudinal_signed
        
        # Calculate form factors
        Q_sinThOverLambda = Q_mag / (4 * np.pi)
        
        # For AuTe2: first atom is Au, next two are Te
        fAu = CalcAtomicfQ(Q_mag, 'Au', scale=4*np.pi, use_xraylib=False)
        fTe = CalcAtomicfQ(Q_mag, 'Te', scale=4*np.pi, use_xraylib=False)
        
        result['Q_sinThOverLambda'] = Q_sinThOverLambda
        result['form_factors'] = {'Au': fAu, 'Te': fTe}
        
        # Prepare for IXS calculation
        fQ_matrix = np.array([[fAu, fTe, fTe]])
        
        # Calculate IXS
        Is, Ias, n, F, xs_info = calc_ixs(
            w_THz.reshape(1, -1), ev, Q_prim.reshape(1, 3),  # Use FULL Q_prim for IXS
            self.xtal.b_l, self.xtal.xs, fQ_matrix, self.masses, self.kT_THz,
            units='barn', per_steradian=True
        )
        
        result['IXS_stokes'] = Is.flatten()
        result['IXS_antistokes'] = Ias.flatten()
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
        
        # Analyzer positions relative to center (from sixcircle ca6 output)
        # Row and column offsets in (dH, dK, dL) - approximate for now
        analyzers = {
            # Row 0 (a03-a09)
            'a03': (-0.142, -0.146, -0.006),
            'a04': (-0.095, -0.098, -0.007),
            'a05': (-0.048, -0.049, -0.035),
            'a06': (0.000, 0.000, 0.000),  # Center
            'a07': (0.048, 0.049, 0.034),
            'a08': (0.096, 0.099, 0.067),
            'a09': (0.144, 0.149, 0.099),
            
            # Row 1 (a14-a20)  
            'a14': (-0.197, -0.093, -0.006),
            'a15': (-0.150, -0.045, -0.007),
            'a16': (-0.103, 0.004, -0.034),
            'a17': (-0.055, 0.053, 0.001),
            'a18': (-0.007, 0.102, 0.035),
            'a19': (0.041, 0.152, 0.068),
            'a20': (0.089, 0.202, 0.100),
            
            # Row 2 (a25-a31)
            'a25': (-0.252, -0.040, -0.006),
            'a26': (-0.205, 0.009, -0.007),
            'a27': (-0.158, 0.057, -0.034),
            'a28': (-0.110, 0.106, 0.001),
            'a29': (-0.062, 0.156, 0.035),
            'a30': (-0.014, 0.205, 0.068),
            'a31': (0.034, 0.255, 0.100),
            
            # Row 3 (a35-a41)
            'a35': (-0.306, 0.014, -0.006),
            'a36': (-0.260, 0.062, -0.007),
            'a37': (-0.212, 0.111, -0.035),
            'a38': (-0.165, 0.160, 0.000),
            'a39': (-0.117, 0.209, 0.034),
            'a40': (-0.069, 0.259, 0.067),
            'a41': (-0.020, 0.309, 0.099),
        }
        
        Q_center = np.asarray(Q_center)
        
        # Convert to conventional if needed
        if coords.lower() == 'primitive':
            Q_center_conv = aute2_prim2conv_k(Q_center)
        else:
            Q_center_conv = Q_center
        
        array_results = []
        
        # Suppress form factor messages during batch calculation
        import sys
        import io
        
        for name in sorted(analyzers.keys()):
            dQ = np.array(analyzers[name])
            Q_analyzer = Q_center_conv + dQ
            
            # Analyze this Q point (suppress printing and form factor messages)
            old_stdout = sys.stdout
            sys.stdout = io.StringIO()
            try:
                result = self.analyze(Q_analyzer, coords='conventional', 
                                     freq_unit=freq_unit, print_results=False)
            finally:
                sys.stdout = old_stdout
            
            result['analyzer'] = name
            result['Q_offset'] = dQ
            array_results.append(result)
        
        if print_results:
            self._print_array_results(array_results, freq_unit)
        
        return array_results
    
    def _print_array_results(self, array_results, freq_unit='meV'):
        """Print simplified analyzer array table"""
        
        print('\n' + '='*80)
        print(f'Analyzer Array Results - {len(array_results)} analyzers')
        print('Note: Analyzer positions approximate, based on BL43LXU geometry at ~30 nm^-1')
        print('='*80)
        
        # Header
        freq_label = {'meV': 'meV', 'cm-1': 'cm^-1', 'THz': 'THz'}[freq_unit]
        
        print(f'{"Ana":>4s}  {"H":>7s} {"K":>7s} {"L":>7s}  {"   Frequencies (" + freq_label + ")":60s}  {"IXS (barn/uc*sr)"}')
        print('-'*80)
        
        freq_key = {'meV': 'frequencies_meV', 
                   'cm-1': 'frequencies_cm',
                   'THz': 'frequencies_THz'}[freq_unit]
        
        for result in array_results:
            Q = result['Q_conv']
            freqs = result[freq_key]
            ixs_stokes = result['IXS_stokes']
            
            # Show first 6 modes
            freq_str = '  '.join([f'{f:5.1f}' for f in freqs[:6]])
            ixs_str = '  '.join([f'{xs:5.2f}' if xs > 0.01 else ' ~0  ' 
                                for xs in ixs_stokes[:6]])
            
            print(f'{result["analyzer"]:>4s}  {Q[0]:7.3f} {Q[1]:7.3f} {Q[2]:7.3f}  '
                  f'{freq_str}  {ixs_str}')
        
        print('='* 80 + '\n')


    def _print_array_results(self, array_results, freq_unit='meV'):
        """Print simplified analyzer array table"""
        
        print('\n' + '='* 80)
        print(f'Analyzer Array Results - {len(array_results)} analyzers')
        print('Note: Analyzer positions approximate, based on BL43LXU geometry at ~30 nm^-1')
        print('='* 80)
        
        # Header
        freq_label = {'meV': 'meV', 'cm-1': 'cm^-1', 'THz': 'THz'}[freq_unit]
        
        print(f'{"Ana":>4s}  {"H":>7s} {"K":>7s} {"L":>7s}  {"   Frequencies (" + freq_label + ")":60s}  {"IXS (barn/uc*sr)"}')
        print('-'* 80)
        
        freq_key = {'meV': 'frequencies_meV', 
                   'cm-1': 'frequencies_cm',
                   'THz': 'frequencies_THz'}[freq_unit]
        
        for result in array_results:
            Q = result['Q_conv']
            freqs = result[freq_key]
            ixs_stokes = result['IXS_stokes']
            
            # Show first 6 modes
            freq_str = '  '.join([f'{f:5.1f}' for f in freqs[:6]])
            ixs_str = '  '.join([f'{xs:5.2f}' if xs > 0.01 else ' ~0  ' 
                                for xs in ixs_stokes[:6]])
            
            print(f'{result["analyzer"]:>4s}  {Q[0]:7.3f} {Q[1]:7.3f} {Q[2]:7.3f}  '
                  f'{freq_str}  {ixs_str}')
        
        print('='*80 + '\n')

    def _print_results(self, result, detailed=False, freq_unit='meV'):
        """Print analysis results - compact format"""
        
        Q_conv = result['Q_conv']
        Q_prim = result['Q_prim']
        Q_cart = result['Q_cart']
        Q_mag = result['Q_mag']
        
        # Temperature and form factors
        kT_K = self.kT_THz * const.THz2meV / 1000 / 8.617333262e-5
        print(f'\nT={kT_K:.0f}K  '
              f'sin(θ)/λ={result["Q_sinThOverLambda"]:.4f} A^-1  '
              f'f(Au)={result["form_factors"]["Au"]:.1f}  '
              f'f(Te)={result["form_factors"]["Te"]:.1f}')
        
        # Calculate structure factor for elastic scattering
        # F(Q) = sum_atoms f_atom(Q) * exp(2pii Q_prim * r_frac)
        F = 0.0 + 0.0j
        Q_prim_vec = Q_prim
        
        # AuTe2: atom 0 is Au, atoms 1,2 are Te
        f_Au = result["form_factors"]["Au"]
        f_Te = result["form_factors"]["Te"]
        
        for iat in range(self.xtal.nat):
            r_frac = self.xtal.xs[iat]  # Fractional coordinates
            phase = 2 * np.pi * np.dot(Q_prim_vec, r_frac)
            
            if iat == 0:  # Au
                F += f_Au * np.exp(1j * phase)
            else:  # Te
                F += f_Te * np.exp(1j * phase)
        
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
        
        print('\n' + '='* 80)
        print(f'Mode  Freq({freq_label:>3s})    IXS(S)   IXS(AS)    Pol   Au:Q*e,ph   Te1:Q*e,ph  Te2:Q*e,ph')
        print('-'* 80)
        
        # Phonon modes
        w_meV = result['frequencies_meV']
        long_char = result['long_char']
        Is = result['IXS_stokes']
        Ias = result['IXS_antistokes']
        
        for i in range(len(w_meV)):
            # Determine polarization
            L = long_char[i]
            if np.isnan(L):
                pol = 'M'
            elif L > 0.7:
                pol = 'L'
            elif L < 0.3:
                pol = 'T'
            elif L > 0.5:
                pol = 'M(L)'
            else:
                pol = 'M(T)'
            
            # Format IXS values (0.00 to 999.99)
            # Show --- for acoustic modes at or near Gamma point
            if at_gamma and freq_data[i] < 1.0:  # Acoustic mode at Gamma
                ixs_s_str = "   ---"
                ixs_as_str = "   ---"
            else:
                ixs_s = min(999.99, max(0.0, Is[i]))
                ixs_as = min(999.99, max(0.0, Ias[i]))
                ixs_s_str = f'{ixs_s:6.2f}'
                ixs_as_str = f'{ixs_as:6.2f}'
            
            # Calculate Q*e projections for each atom
            Q_hat_full = Q_cart / Q_mag if Q_mag > 1e-10 else np.zeros(3)
            
            ev_mode = result['eigenvectors'][:, i, 0]
            ev_reshaped = ev_mode.reshape(3, self.xtal.nat)
            
            # Calculate Q*e for each atom (complex)
            q_dot_e = np.zeros(self.xtal.nat, dtype=complex)
            for iat in range(self.xtal.nat):
                e_atom = ev_reshaped[:, iat]
                q_dot_e[iat] = np.vdot(Q_hat_full, e_atom)
            
            # Find reference atom (largest |Q*e|)
            magnitudes = np.abs(q_dot_e)
            ref_atom = np.argmax(magnitudes)
            ref_phase = np.angle(q_dot_e[ref_atom])
            
            # Calculate relative phases (degrees)
            phases = np.angle(q_dot_e) - ref_phase
            phases_deg = np.degrees(phases)
            phases_deg = np.mod(phases_deg + 180, 360) - 180
            
            # Format eigenvector info
            if Q_mag < 1e-10:
                ev_str = "      ---           ---            ---"
            else:
                ev_parts = []
                for iat in range(self.xtal.nat):
                    mag = magnitudes[iat]
                    phase = phases_deg[iat]
                    ev_parts.append(f'{mag:4.2f},{phase:+4.0f}')
                ev_str = '    '.join(ev_parts)
            
            print(f'{i+1:3d}   {freq_data[i]:7.2f}    {ixs_s_str:>6s}   {ixs_as_str:>6s}    {pol:>4s}     {ev_str}')
        
        # Elastic line with structure factor breakdown
        print('-'* 80)
        
        # Calculate Q*r for each atom (for elastic scattering)
        # Show e^(iQ*r) even at Gamma (all = 1.0, phase = 0)
        if True:  # Always calculate
            elastic_parts = []
            # Show e^(iQ*r) for each atom (amplitude and phase)
            exponentials = np.zeros(self.xtal.nat, dtype=complex)
            
            for iat in range(self.xtal.nat):
                r_frac = self.xtal.xs[iat]
                phase = 2 * np.pi * np.dot(Q_prim, r_frac)
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
        
        print(f'  0      0.00    |F|^2 = {F_squared:8.0f}         El.   {elastic_str}')
        print('='* 80)
        
        #notes
        print('\nNotes: IXS units: barn/(unit cell*sr);  DFT uses unmodulated structure (no CDW).')
        print('\n' + '='* 80 + '\n')
    
    def _format_xs(val):
        """Format cross-section value for display"""
        if val < 1e-20:
            return '    ~0   '
        elif val < 0.001:
            return f'{val:.2e}'
        elif val < 10:
            return f'{val:8.3f}'
        elif val < 10000:
            return f'{val:8.1f}'
        else:
            return f'{val:.2e}'


def interactive_mode():
    """
    Interactive Q-point analysis mode
    Matches MATLAB interactive script
    """
    
    print("=" * 80)
    print("  Interactive IXS Analysis for AuTe2")
    print("=" * 80)
    
    # Load force constants
    fc_file = "data/AuTe_2_m.fc"
    print(f"\nLoading force constants from {fc_file}...")
    
    xtal = ForceConstants(fc_file)
    Phi = xtal.convert_to_eV_per_Angstrom2()
    
    # Get masses
    masses = np.array([xtal.masses[xtal.atom_type_map[i]-1] 
                       for i in range(xtal.nat)])
    
    # Temperature
    kT_cm = 207  # cm^-1
    kT_THz = kT_cm * const.c * 80 / 1e12
    
    print(f"Temperature: {kT_cm:.1f} cm^-1 ({kT_THz:.2f} THz)\n")
    
    # Create analyzer
    analyzer = SingleQAnalyzer(xtal, Phi, masses, kT_THz)
    
    # Default coordinate system
    coord_system = 'conventional'
    freq_unit = 'meV'
    
    print("Instructions:")
    print("  - Enter Q vector as three numbers (e.g., 0.5 0 0)")
    print("  - For satellites, add order m: H K L m (e.g., 1 0 0 1 for first satellite)")
    print("  - Press Enter on empty line to quit")
    print("  - Type 'conv', 'prim', or 'cart' to change coordinate system")
    print("  - Type 'meV', 'THz', or 'invcm' to change frequency units")
    print("  - Type 'array' to calculate full analyzer array at last Q")
    print("  - Type 'angles' to calculate diffractometer angles for last Q")
    print(f"  - Current system: {coord_system}\n")
    
    current_q = None
    
    while True:
        print('─' * 80)
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
            try:
                Q_conv = aute2_prim2conv_k(current_q) if coord_system == 'primitive' else current_q
                analyzer.analyze_array(Q_conv, coords='conventional', freq_unit=freq_unit)
            except Exception as e:
                print(f"\n✗ {e}\n")
                import traceback
                traceback.print_exc()
            continue
            
        elif user_input.lower() == 'array':
            if current_q is None:
                print("\n⚠ Enter a Q point first\n")
                continue
            try:
                Q_conv = aute2_prim2conv_k(current_q) if coord_system == 'primitive' else current_q
                analyzer.analyze_array(Q_conv, coords='conventional', freq_unit=freq_unit)
            except Exception as e:
                print(f"\n✗ {e}\n")
                import traceback
                traceback.print_exc()
            continue
            
        elif user_input.lower() == 'angles':
            if current_q is None:
                print("\n⚠ Enter a Q point first\n")
                continue
            try:
                import sys
                import io
                
                # Suppress sixcircle verbose output
                old_stdout = sys.stdout
                sys.stdout = io.StringIO()
                
                try:
                    sixc = SixCircleInterface()
                finally:
                    sys.stdout = old_stdout
                
                print("\n" + "="* 80)
                print("Diffractometer Angles" + (" (SIMULATION)" if sixc.simulation_mode else ""))
                print("="* 80)
                Q_conv = aute2_prim2conv_k(current_q) if coord_system == 'primitive' else current_q
                print(f"Q (conv): [{Q_conv[0]:.4f}, {Q_conv[1]:.4f}, {Q_conv[2]:.4f}]")
                angles = sixc.move_to_hkl(tuple(Q_conv), check_only=True)
                print("\nAngles:")
                for key, val in angles.items():
                    print(f"  {key:5s} = {val:7.3f}deg")
                print("="* 80 + "\n")
            except Exception as e:
                print(f"\n✗ {e}\n")
            continue
        
        elif user_input.lower() == 'move':
            # Move diffractometer - REAL MOTION!
            if current_q is None:
                print("\n⚠ Enter a Q point first\n")
                continue
            try:
                sixc = SixCircleInterface()
                if sixc.simulation_mode:
                    print("\n⚠ In simulation mode - cannot move\n")
                    continue
                Q_conv = aute2_prim2conv_k(current_q) if coord_system == 'primitive' else current_q
                print("\n" + "="* 80)
                print("⚠ WARNING: Will MOVE diffractometer!")
                print("="* 80)
                print(f"Q (conv): [{Q_conv[0]:.4f}, {Q_conv[1]:.4f}, {Q_conv[2]:.4f}]")
                confirm = input("\nType 'yes' to confirm: ").strip().lower()
                if confirm != 'yes':
                    print("Cancelled\n")
                    continue
                print("\nMoving...")
                sixc.move_to_hkl(tuple(Q_conv), check_only=False)
                print("✓ Complete\n")
            except Exception as e:
                print(f"\n✗ {e}\n")
            continue
        elif user_input.lower().startswith('sixc '):
            # Pass command to sixcircle
            cmd = user_input[5:].strip()  # Remove 'sixc ' prefix
            if not cmd:
                print("\n⚠ Usage: sixc <command>\n")
                print("  Examples: sixc wh, sixc pa, sixc or_check\n")
                continue
            try:
                sixc = SixCircleInterface()
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
    
    # Temperature
    kT_cm = 207
    kT_THz = kT_cm * const.c * 80 / 1e12
    
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
