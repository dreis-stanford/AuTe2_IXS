with open('code/single_q_analysis_si.py', 'r') as f:
    lines = f.readlines()

# Find and fix the form factor calculation (around line 217-225)
for i, line in enumerate(lines):
    if 'fAu = CalcAtomicfQ' in line:
        # Replace the Au/Te form factor section
        lines[i] = '        fSi = CalcAtomicfQ(Q_mag, "Si", scale=4*np.pi, use_xraylib=False)\n'
        lines[i+1] = '        \n'
        lines[i+2] = '        result["form_factors"] = {"Si": fSi}\n'
        lines[i+3] = '        \n'
        lines[i+4] = '        # IXS cross-sections (use FULL Q, not reduced q!)\n'
        lines[i+5] = '        fQ_matrix = np.array([[fSi, fSi]])\n'
        # Remove extra Te lines
        del lines[i+6:i+9]
        break

# Fix the print statements (around line 288-289)
for i, line in enumerate(lines):
    if 'f(Q) Au:' in line:
        lines[i] = '        print(f\'  f(Q) Si:  {result["form_factors"]["Si"]:.2f} (Z=14)\\n\')\n'
        # Remove the Te line
        if i+1 < len(lines) and 'f(Q) Te:' in lines[i+1]:
            del lines[i+1]
        break

# Fix the header
for i, line in enumerate(lines):
    if 'Au%   Au‖    Te1%  Te1‖   Te2%  Te2‖' in line:
        lines[i] = line.replace('Au%   Au‖    Te1%  Te1‖   Te2%  Te2‖', 'Si1%  Si1‖   Si2%  Si2‖')
        break

# Fix title
for i, line in enumerate(lines):
    if 'AuTe₂' in line:
        lines[i] = line.replace('AuTe₂', 'Silicon')
        break

# Fix docstrings at top
for i, line in enumerate(lines):
    if 'AuTe2' in line and i < 20:
        lines[i] = line.replace('AuTe2', 'Silicon')

# Fix imports
for i, line in enumerate(lines):
    if 'from aute2_structure' in line:
        lines[i] = 'from fcc_structure import fcc_conv2prim_k, fcc_prim2conv_k\n'
        break

# Fix function calls
for i, line in enumerate(lines):
    if 'aute2_conv2prim_k' in line:
        lines[i] = line.replace('aute2_conv2prim_k', 'fcc_conv2prim_k')
    if 'aute2_prim2conv_k' in line:
        lines[i] = line.replace('aute2_prim2conv_k', 'fcc_prim2conv_k')

# Fix force constants file
for i, line in enumerate(lines):
    if 'fc_file = "data/AuTe_2_m.fc"' in line:
        lines[i] = '    fc_file = "data/Test__Silicon_dispersion/Qgrid_888/Cg.fc"\n'
        break

with open('code/single_q_analysis_si.py', 'w') as f:
    f.writelines(lines)

print("✓ Fixed Silicon version manually")
