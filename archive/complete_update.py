"""
Complete update: add longitudinal components and unit switching
"""
import re

with open('code/single_q_analysis.py', 'r') as f:
    content = f.read()

# 1. Add longitudinal calculation after atom_participation
pattern = r"(result\['atom_participation'\] = atom_participation)"
replacement = r"""\1
        
        # Calculate signed longitudinal components (Q·e for each atom)
        Q_cart = result['Q_cart']
        Q_mag = result['Q_mag']
        Q_hat = Q_cart / Q_mag if Q_mag > 1e-10 else np.zeros(3)
        
        longitudinal_signed = np.zeros((nmodes, self.xtal.nat))
        
        for imode in range(nmodes):
            ev_mode = ev[:, imode]
            total_ev_norm = np.linalg.norm(ev_mode)
            
            for iat in range(self.xtal.nat):
                e_atom = ev_mode[3*iat:3*iat+3]
                Q_dot_e = np.real(np.vdot(Q_hat, e_atom))
                longitudinal_signed[imode, iat] = 100 * Q_dot_e / (total_ev_norm + 1e-12)
        
        result['longitudinal_signed'] = longitudinal_signed"""

content = re.sub(pattern, replacement, content)

# 2. Update analyze signature to include freq_unit
content = content.replace(
    "def analyze(self, Q_input, coords='primitive', print_results=True,",
    "def analyze(self, Q_input, coords='primitive', freq_unit='meV', print_results=True,"
)

# 3. Pass freq_unit to _print_results
content = content.replace(
    "self._print_results(result, detailed=print_detailed)",
    "self._print_results(result, detailed=print_detailed, freq_unit=freq_unit)"
)

# 4. Update _print_results signature
content = content.replace(
    "def _print_results(self, result, detailed=False):",
    "def _print_results(self, result, detailed=False, freq_unit='meV'):"
)

# 5. Update table header
old_header = "        print('Mode  Freq(cm⁻¹)  Freq(meV)   L-char   IXS(S)    IXS(AS)   Au%   Te1%  Te2%  Polarization')"
new_header = """        # Choose frequency unit
        freq_label = {'meV': 'meV', 'cm-1': 'cm⁻¹', 'THz': 'THz'}[freq_unit]
        freq_data = {'meV': result['frequencies_meV'], 
                     'cm-1': result['frequencies_cm'],
                     'THz': result['frequencies_THz']}[freq_unit]
        
        print(f'Mode  Freq({freq_label:>4s})   L-char   IXS(S)    IXS(AS)   Au%   Au‖    Te1%  Te1‖   Te2%  Te2‖   Pol')"""

content = content.replace(old_header, new_header)

# 6. Add long_sign and update print loop
old_section = """        atom_part = result['atom_participation']
        pol_type = result['pol_type']
        
        for i in range(len(w_cm)):
            # Format IXS values
            ixs_s_str = self._format_xs(Is[i])
            ixs_as_str = self._format_xs(Ias[i])
            
            print(f'{i+1:2d}    {w_cm[i]:8.2f}    {w_meV[i]:7.2f}     {long_char[i]:5.3f}   '
                  f'{ixs_s_str:9s}  {ixs_as_str:9s}  '
                  f'{atom_part[i,0]*100:4.1f}  {atom_part[i,1]*100:4.1f}  {atom_part[i,2]*100:4.1f}  '
                  f'{pol_type[i]}')"""

new_section = """        atom_part = result['atom_participation']
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
                  f'{pol_type[i]}')"""

content = content.replace(old_section, new_section)

# 7. Add freq_unit to interactive mode
# Initialize freq_unit
content = content.replace(
    "coord_system = 'primitive'\n    print(f'  - Current system: {coord_system}')",
    "coord_system = 'primitive'\n    freq_unit = 'meV'\n    print(f'  - Current: {coord_system} coords, {freq_unit} units')"
)

# Add unit switching commands
old_try = """        # Check for command keywords first
        if user_input.lower() == 'prim':
            coord_system = 'primitive'
            print("  → Switched to primitive coordinates\\n")
            continue
        elif user_input.lower() == 'conv':
            coord_system = 'conventional'
            print("  → Switched to conventional coordinates\\n")
            continue
        elif user_input.lower() == 'cart':
            coord_system = 'cartesian'
            print("  → Switched to Cartesian coordinates\\n")
            continue
        elif user_input.lower() == 'mev':
            freq_unit = 'meV'
            print("  → Switched to meV\\n")
            continue
        elif user_input.lower() in ['cm', 'cm-1', 'invcm']:
            freq_unit = 'cm-1'
            print("  → Switched to cm⁻¹\\n")
            continue
        elif user_input.lower() == 'thz':
            freq_unit = 'THz'
            print("  → Switched to THz\\n")
            continue
        
        try:"""

if "Check for command keywords first" not in content:
    content = content.replace(
        "        try:\n            Q = np.array([float(x) for x in user_input.split()])",
        old_try + "\n            Q = np.array([float(x) for x in user_input.split()])"
    )

# Pass freq_unit to analyze in interactive mode
content = content.replace(
    "result = analyzer.analyze(Q_input, coords=input_coords, \n                                     print_results=True, print_detailed=False)",
    "result = analyzer.analyze(Q_input, coords=input_coords,\n                                     freq_unit=freq_unit,\n                                     print_results=True, print_detailed=False)"
)

# Update instructions
old_instr = """    print("  - Type 'conv', 'prim', or 'cart' to change coordinate system")"""
new_instr = """    print("  - Type 'conv', 'prim', or 'cart' to change coordinate system")
    print("  - Type 'meV', 'THz', or 'invcm' to change frequency units")"""

content = content.replace(old_instr, new_instr)

with open('code/single_q_analysis.py', 'w') as f:
    f.write(content)

print("✓ Complete update applied")
print("  - Longitudinal components added")
print("  - Frequency unit switching added")
print("  - Interactive mode updated")
