with open('code/single_q_analysis_si.py', 'r') as f:
    lines = f.readlines()

# Find and replace the form factors section (around lines 215-225)
for i in range(len(lines)):
    if 'For Silicon: first atom is Au' in lines[i]:
        # Replace lines 217-224 with corrected Silicon version
        lines[i] = '        # For Silicon: both atoms are Si\n'
        lines[i+1] = '        fSi = CalcAtomicfQ(Q_mag, "Si", scale=4*np.pi, use_xraylib=False)\n'
        lines[i+2] = '        \n'
        lines[i+3] = '        result["Q_sinThOverLambda"] = Q_sinThOverLambda\n'
        lines[i+4] = '        result["form_factors"] = {"Si": fSi}\n'
        lines[i+5] = '        \n'
        lines[i+6] = '        # Prepare for IXS calculation\n'
        lines[i+7] = '        fQ_matrix = np.array([[fSi, fSi]])\n'
        break

with open('code/single_q_analysis_si.py', 'w') as f:
    f.writelines(lines)

print("✓ Fixed form factors section")
