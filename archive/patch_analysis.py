"""
Add signed longitudinal components to single_q_analysis.py
"""

import re

# Read the current file
with open('code/single_q_analysis.py', 'r') as f:
    content = f.read()

# Find the atom_participation calculation section (around line 177-194)
# We need to add longitudinal component calculation after it

# Pattern to find: after the atom_participation calculation
pattern = r"(result\['atom_participation'\] = atom_participation)"

# Replacement: add longitudinal calculation
replacement = r"""\1
        
        # Calculate signed longitudinal components (Q·e for each atom)
        # Use FULL Q (not reduced q) since that's what matters for scattering
        Q_cart = result['Q_cart']
        Q_mag = result['Q_mag']
        Q_hat = Q_cart / Q_mag if Q_mag > 1e-10 else np.zeros(3)
        
        longitudinal_signed = np.zeros((nmodes, self.xtal.nat))
        
        for imode in range(nmodes):
            ev_mode = ev[:, imode]
            
            for iat in range(self.xtal.nat):
                # Get displacement vector for this atom (in Cartesian)
                e_atom = ev_mode[3*iat:3*iat+3]
                
                # Signed projection onto Q: Q̂·e (real part)
                # Positive = along Q, Negative = opposite to Q
                Q_dot_e = np.real(np.dot(Q_hat, e_atom))
                
                # Normalize by total eigenvector norm for this mode
                total_ev_norm = np.linalg.norm(ev_mode)
                
                # Store as signed percentage-like value
                longitudinal_signed[imode, iat] = 100 * Q_dot_e / (total_ev_norm + 1e-12)
        
        result['longitudinal_signed'] = longitudinal_signed"""

content = re.sub(pattern, replacement, content)

# Now update the print format (around line 279)
# Change the header
old_header = "        print('Mode  Freq(cm⁻¹)  Freq(meV)   L-char   IXS(S)    IXS(AS)   Au%   Te1%  Te2%  Polarization')"
new_header = """        print('='*125)
        print('Mode  Freq(cm⁻¹)  Freq(meV)   L-char   IXS(S)    IXS(AS)   Au%   Au‖    Te1%  Te1‖   Te2%  Te2‖   Pol')"""

content = content.replace(old_header, new_header)

# Update the data printing loop (around line 282-294)
# Find the loop that prints each mode
old_print_pattern = r"(\s+)(print\(f\"{imode\+1:4d}\s+{w_cm\[imode\]:10\.2f}.*?\))"

# This is complex, let's just add a note for manual editing
print("File patched for longitudinal calculation.")
print("\nYou need to manually update the print statement around line 282-294")
print("to include longitudinal_signed values.")
print("\nLook for the line that prints mode data and add:")
print("  long_sign = result['longitudinal_signed']")
print("  Then in the print statement add: {long_sign[imode, 0]:+6.1f} etc.")

# Write the modified content
with open('code/single_q_analysis.py', 'w') as f:
    f.write(content)

print("\n✓ Patched: Added longitudinal_signed calculation")
print("✓ Backup saved to: code/single_q_analysis.py.backup")
