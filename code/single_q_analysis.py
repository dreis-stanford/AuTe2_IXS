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
        
        # Convert to Cartesian (2π/Å)
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
        
        if Q_mag_reduced > 1e-6:
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
        
        # Calculate signed longitudinal components (Q·e for each atom)
        Q_cart = result['Q_cart']
        Q_mag = result['Q_mag']
        Q_hat = Q_cart_reduced / Q_mag_reduced if Q_mag > 1e-10 else np.zeros(3)
        
        longitudinal_signed = np.zeros((nmodes, self.xtal.nat))
        
        for imode in range(nmodes):
            ev_mode = ev[:, imode]
            total_ev_norm = np.linalg.norm(ev_mode)
            
            for iat in range(self.xtal.nat):
                e_atom = ev_mode[3*iat:3*iat+3]
                Q_dot_e = np.real(np.vdot(Q_hat, e_atom))
                longitudinal_signed[imode, iat] = 100 * Q_dot_e / (total_ev_norm + 1e-12)
        
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
    
    def _print_results(self, result, detailed=False, freq_unit='meV'):
        """Print analysis results"""
        
        print('\n' + '╔' + '='*64 + '╗')
        print('║  IXS Analysis at Single Q Point - AuTe₂' + ' '*24 + '║')
        print('╚' + '='*64 + '╝\n')
        
        Q_conv = result['Q_conv']
        Q_prim = result['Q_prim']
        Q_cart = result['Q_cart']
        Q_mag = result['Q_mag']
        
        # Check if this is a satellite reflection
        q_mod = self.mod_struct.q_mod
        is_satellite = False
        satellite_order = 0
        Q_main = Q_conv
        
        # Try to determine if Q is near a satellite position
        for m in range(-3, 4):
            for h in range(-3, 4):
                for k in range(-3, 4):
                    for l in range(-3, 4):
                        Q_test = np.array([h, k, l]) + m * q_mod
                        if np.linalg.norm(Q_conv - Q_test) < 0.01:
                            is_satellite = (m != 0)
                            satellite_order = m
                            Q_main = np.array([h, k, l])
                            break
        
        if is_satellite:
            print(f'Q (conventional): [{Q_conv[0]:.4f}, {Q_conv[1]:.4f}, {Q_conv[2]:.4f}] (r.l.u.) ★ SATELLITE m={satellite_order:+d}')
            print(f'  Main peak:      [{Q_main[0]:.0f}, {Q_main[1]:.0f}, {Q_main[2]:.0f}]')
            print(f'  q_mod:          [{q_mod[0]:.4f}, {q_mod[1]:.4f}, {q_mod[2]:.4f}]')
        else:
            print(f'Q (conventional): [{Q_conv[0]:.4f}, {Q_conv[1]:.4f}, {Q_conv[2]:.4f}] (r.l.u.)')
        print(f'Q (primitive):    [{Q_prim[0]:.4f}, {Q_prim[1]:.4f}, {Q_prim[2]:.4f}] (r.l.u.)')
        print(f'Q (Cartesian):    [{Q_cart[0]:.4f}, {Q_cart[1]:.4f}, {Q_cart[2]:.4f}] (2π/Å)')
        print(f'|Q| = {Q_mag:.4f} (2π/Å)\n')
        
        Q_conv_G = result['G_conv']
        Q_prim_G = result['G_prim']
        Q_cart_G = result['G_cart']
        
        print('Nearest reciprocal lattice vector G:')
        print(f'  G (conventional): [{Q_conv_G[0]:.0f}, {Q_conv_G[1]:.0f}, {Q_conv_G[2]:.0f}] (r.l.u.)')
        print(f'  G (primitive):    [{Q_prim_G[0]:.0f}, {Q_prim_G[1]:.0f}, {Q_prim_G[2]:.0f}] (r.l.u.)')
        print(f'  G (Cartesian):    [{Q_cart_G[0]:.4f}, {Q_cart_G[1]:.4f}, {Q_cart_G[2]:.4f}] (2π/Å)\n')
        
        Q_conv_red = result['Q_reduced_conv']
        Q_prim_red = result['Q_reduced_prim']
        Q_cart_red = result['Q_reduced_cart']
        
        print('Reduced q = Q - G (first Brillouin zone):')
        print(f'  q (conventional): [{Q_conv_red[0]:.4f}, {Q_conv_red[1]:.4f}, {Q_conv_red[2]:.4f}] (r.l.u.)')
        print(f'  q (primitive):    [{Q_prim_red[0]:.4f}, {Q_prim_red[1]:.4f}, {Q_prim_red[2]:.4f}] (r.l.u.)')
        print(f'  q (Cartesian):    [{Q_cart_red[0]:.4f}, {Q_cart_red[1]:.4f}, {Q_cart_red[2]:.4f}] (2π/Å)')
        print(f'  |q| = {np.linalg.norm(Q_cart_red):.4f} (2π/Å)\n')
        
        # Temperature
        kT_K = self.kT_THz * const.THz2meV / 1000 / 8.617333262e-5
        print(f'Temperature: {self.kT_THz:.1f} THz ({kT_K:.1f} K)\n')
        
        # Form factors
        print('Form factors:')
        print(f'  sin(θ)/λ = {result["Q_sinThOverLambda"]:.4f} Å⁻¹')
        print(f'  f(Q) Au:  {result["form_factors"]["Au"]:.2f} (Z=79)')
        print(f'  f(Q) Te:  {result["form_factors"]["Te"]:.2f} (Z=52)\n')
        
        print(f'Cross-section units: {result["cross_section_info"]["units"]}')
        # DFT modulation note
        print("\nNote: DFT uses unmodulated AuTe₂ structure (no CDW in calculation).\n")
        
        # Mode table
        print('═'*95)
        # Choose frequency unit
        freq_label = {'meV': 'meV', 'cm-1': 'cm⁻¹', 'THz': 'THz'}[freq_unit]
        freq_data = {'meV': result['frequencies_meV'], 
                     'cm-1': result['frequencies_cm'],
                     'THz': result['frequencies_THz']}[freq_unit]
        
        print(f'Mode  Freq({freq_label:>4s})   L-char   IXS(S)    IXS(AS)   Au%   Au‖    Te1%  Te1‖   Te2%  Te2‖   Pol')
        print('─'*95)
        
        w_cm = result['frequencies_cm']
        w_meV = result['frequencies_meV']
        long_char = result['long_char']
        Is = result['IXS_stokes']
        Ias = result['IXS_antistokes']
        atom_part = result['atom_participation']
        long_sign = result['longitudinal_signed']
        pol_type = result['pol_type']
        
        for i in range(len(w_cm)):
            # Format IXS values
            ixs_s_str = self._format_xs(Is[i])
            ixs_as_str = self._format_xs(Ias[i])
            
            print(f'{i+1:2d}    {freq_data[i]:8.2f}     {long_char[i]:5.3f}   '
                  f'{ixs_s_str:9s}  {ixs_as_str:9s}  '
                  f'{atom_part[i,0]*100:4.1f} {long_sign[i,0]:+5.1f}  '
                  f'{atom_part[i,1]*100:4.1f} {long_sign[i,1]:+5.1f}  '
                  f'{atom_part[i,2]*100:4.1f} {long_sign[i,2]:+5.1f}  '
                  f'{pol_type[i]}')
        
        print('═'*95 + '\n')
    
    @staticmethod
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
    
    print("=" * 70)
    print("  Interactive IXS Analysis for AuTe₂")
    print("=" * 70)
    
    # Load force constants
    fc_file = "data/AuTe_2_m.fc"
    print(f"\nLoading force constants from {fc_file}...")
    
    xtal = ForceConstants(fc_file)
    Phi = xtal.convert_to_eV_per_Angstrom2()
    
    # Get masses
    masses = np.array([xtal.masses[xtal.atom_type_map[i]-1] 
                       for i in range(xtal.nat)])
    
    # Temperature
    kT_cm = 207  # cm⁻¹
    kT_THz = kT_cm * const.c * 100 / 1e12
    
    print(f"Temperature: {kT_cm:.1f} cm⁻¹ ({kT_THz:.2f} THz)\n")
    
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
    print("  - Type 'angles' to calculate diffractometer angles for last Q")
    print(f"  - Current system: {coord_system}\n")
    
    current_q = None
    
    while True:
        print('─' * 70)
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
            print("  → Switched to cm⁻¹\n")
            continue
        elif user_input.lower() == 'thz':
            freq_unit = 'THz'
            print(f"  → Switched to THz\n")
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
                
                print("\n" + "="*70)
                print("Diffractometer Angles" + (" (SIMULATION)" if sixc.simulation_mode else ""))
                print("="*70)
                Q_conv = aute2_prim2conv_k(current_q) if coord_system == 'primitive' else current_q
                print(f"Q (conv): [{Q_conv[0]:.4f}, {Q_conv[1]:.4f}, {Q_conv[2]:.4f}]")
                angles = sixc.move_to_hkl(tuple(Q_conv), check_only=True)
                print("\nAngles:")
                for key, val in angles.items():
                    print(f"  {key:5s} = {val:7.3f}°")
                print("="*70 + "\n")
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
                print("\n" + "="*70)
                print("⚠ WARNING: Will MOVE diffractometer!")
                print("="*70)
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
                print(f'Converted Cartesian Q = [{Q[0]:.4f}, {Q[1]:.4f}, {Q[2]:.4f}] (2π/Å)')
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
    
    print("\n" + "=" * 70)
    print("Session ended.")
    print("=" * 70)


if __name__ == "__main__":
    import os
    sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
    
    # Test with single Q-point
    print("=" * 70)
    print("Testing Single Q-Point Analyzer")
    print("=" * 70)
    
    # Load data
    fc_file = "data/AuTe_2_m.fc"
    xtal = ForceConstants(fc_file)
    Phi = xtal.convert_to_eV_per_Angstrom2()
    
    masses = np.array([xtal.masses[xtal.atom_type_map[i]-1] 
                       for i in range(xtal.nat)])
    
    # Temperature
    kT_cm = 207
    kT_THz = kT_cm * const.c * 100 / 1e12
    
    # Create analyzer
    analyzer = SingleQAnalyzer(xtal, Phi, masses, kT_THz)
    
    # Test analysis at Q = (0.5, 0, 0) primitive
    print("\nTest: Q = (0.5, 0, 0) primitive")
    result = analyzer.analyze([0.5, 0.0, 0.0], coords='primitive')
    
    print("\n" + "=" * 70)
    print("Test complete! To run interactive mode:")
    print("  from single_q_analysis import interactive_mode")
    print("  interactive_mode()")
    print("=" * 70)
