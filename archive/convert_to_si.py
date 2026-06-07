"""Convert single_q_analysis.py to work for Silicon"""

with open('code/single_q_analysis_si.py', 'r') as f:
    content = f.read()

# 1. Update header/docstrings
content = content.replace('Interactive IXS analysis for single Q-points in AuTe2', 
                          'Interactive IXS analysis for single Q-points in Silicon')
content = content.replace('Matches MATLAB calcIXS_single_Q_AuTe2.m and interactive script',
                          'Based on single_q_analysis.py, adapted for Silicon FCC structure')
content = content.replace('MATLAB calcIXS_single_Q_AuTe2.m',
                          'single_q_analysis.py for AuTe2')

# 2. Change imports
content = content.replace('from aute2_structure import aute2_conv2prim_k, aute2_prim2conv_k',
                          'from fcc_structure import fcc_conv2prim_k, fcc_prim2conv_k')

# 3. Change coordinate transformation function calls
content = content.replace('aute2_conv2prim_k', 'fcc_conv2prim_k')
content = content.replace('aute2_prim2conv_k', 'fcc_prim2conv_k')

# 4. Update form factors (Silicon has 2 atoms, both Si)
old_form = """        # For AuTe2: first atom is Au, next two are Te
        fAu = CalcAtomicfQ(Q_mag, 'Au', scale=4*np.pi, use_xraylib=False)
        fTe = CalcAtomicfQ(Q_mag, 'Te', scale=4*np.pi, use_xraylib=False)
        
        # Form factor matrix
        result['form_factors'] = {'Au': fAu, 'Te': fTe}
        
        # IXS cross-sections (use FULL Q, not reduced q!)
        fQ_matrix = np.array([[fAu, fTe, fTe]])"""

new_form = """        # For Silicon: both atoms are Si
        fSi = CalcAtomicfQ(Q_mag, 'Si', scale=4*np.pi, use_xraylib=False)
        
        # Form factor matrix (both atoms same)
        result['form_factors'] = {'Si': fSi}
        
        # IXS cross-sections (use FULL Q, not reduced q!)
        fQ_matrix = np.array([[fSi, fSi]])"""

content = content.replace(old_form, new_form)

# 5. Update title
content = content.replace("IXS Analysis at Single Q Point - AuTe₂",
                          "IXS Analysis at Single Q Point - Silicon")

# 6. Update form factor printing
old_print = """        print(f'  f(Q) Au:  {result["form_factors"]["Au"]:.2f} (Z=79)')
        print(f'  f(Q) Te:  {result["form_factors"]["Te"]:.2f} (Z=52)\\n')"""

new_print = """        print(f'  f(Q) Si:  {result["form_factors"]["Si"]:.2f} (Z=14)\\n')"""

content = content.replace(old_print, new_print)

# 7. Update table header (Si1, Si2 instead of Au, Te1, Te2)
old_header = "print(f'Mode  Freq({freq_label:>4s})   L-char   IXS(S)    IXS(AS)   Au%   Au‖    Te1%  Te1‖   Te2%  Te2‖   Pol')"
new_header = "print(f'Mode  Freq({freq_label:>4s})   L-char   IXS(S)    IXS(AS)   Si1%  Si1‖   Si2%  Si2‖')"

content = content.replace(old_header, new_header)

# 8. Update force constants file in interactive mode
content = content.replace('fc_file = "data/AuTe_2_m.fc"',
                          'fc_file = "data/Test__Silicon_dispersion/Qgrid_888/Cg.fc"')

# 9. Update interactive mode title
content = content.replace('Interactive IXS Analysis for AuTe₂',
                          'Interactive IXS Analysis for Silicon')

with open('code/single_q_analysis_si.py', 'w') as f:
    f.write(content)

print("✓ Converted to Silicon version")
print("Changes:")
print("  - FCC coordinate transformations")
print("  - Silicon form factors (2 atoms)")
print("  - Updated headers and labels")
print("  - Si force constants file")
