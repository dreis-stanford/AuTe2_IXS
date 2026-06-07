import re

with open('code/single_q_analysis.py', 'r') as f:
    content = f.read()

# 1. Add longitudinal calculation after atom_participation
insert_pos = content.find("result['atom_participation'] = atom_participation")
if insert_pos > 0:
    insert_pos = content.find('\n', insert_pos) + 1
    
    longitudinal_code = """
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
                Q_dot_e = np.real(np.dot(Q_hat, e_atom))
                longitudinal_signed[imode, iat] = 100 * Q_dot_e / (total_ev_norm + 1e-12)
        
        result['longitudinal_signed'] = longitudinal_signed
"""
    content = content[:insert_pos] + longitudinal_code + content[insert_pos:]
    print("✓ Added longitudinal calculation")

# 2. Update header line
content = content.replace(
    "        print('Mode  Freq(cm⁻¹)  Freq(meV)   L-char   IXS(S)    IXS(AS)   Au%   Te1%  Te2%  Polarization')",
    "        print('Mode  Freq(cm⁻¹)  Freq(meV)   L-char   IXS(S)    IXS(AS)   Au%   Au‖    Te1%  Te1‖   Te2%  Te2‖   Pol')"
)
print("✓ Updated header")

# 3. Add long_sign extraction and update print statement
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
            
            print(f'{i+1:2d}    {w_cm[i]:8.2f}    {w_meV[i]:7.2f}     {long_char[i]:5.3f}   '
                  f'{ixs_s_str:9s}  {ixs_as_str:9s}  '
                  f'{atom_part[i,0]*100:4.1f} {long_sign[i,0]:+5.1f}  '
                  f'{atom_part[i,1]*100:4.1f} {long_sign[i,1]:+5.1f}  '
                  f'{atom_part[i,2]*100:4.1f} {long_sign[i,2]:+5.1f}  '
                  f'{pol_type[i]}')"""

content = content.replace(old_section, new_section)
print("✓ Updated print statement")

with open('code/single_q_analysis.py', 'w') as f:
    f.write(content)

print("\n✓ All updates complete!")
