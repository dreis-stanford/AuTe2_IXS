with open('code/single_q_analysis_si.py', 'r') as f:
    lines = f.readlines()

# Find and fix the print statement (around line 318-321)
for i, line in enumerate(lines):
    if "print(f'{i+1:2d}    {freq_data[i]:8.2f}" in line:
        # Replace the whole print statement with corrected 2-atom version
        lines[i] = "            print(f'{i+1:2d}    {freq_data[i]:8.2f}     {long_char[i]:5.3f}   '\n"
        lines[i+1] = "                  f'{ixs_s_str:9s}  {ixs_as_str:9s}  '\n"
        lines[i+2] = "                  f'{atom_part[i,0]*100:4.1f} {long_sign[i,0]:+5.1f}  '\n"
        lines[i+3] = "                  f'{atom_part[i,1]*100:4.1f} {long_sign[i,1]:+5.1f}')\n"
        # Delete any extra lines
        if i+4 < len(lines) and 'atom_part' not in lines[i+4]:
            # Keep the rest
            pass
        print(f"Fixed print statement at line {i+1}")
        break

with open('code/single_q_analysis_si.py', 'w') as f:
    f.writelines(lines)
