# Update the print section to include longitudinal components

with open('code/single_q_analysis.py', 'r') as f:
    lines = f.readlines()

# Find and update line 323 (the print statement for each mode)
for i, line in enumerate(lines):
    if i == 322:  # Line 323 (0-indexed)
        # Old line prints: atom_part[i,0]*100, atom_part[i,1]*100, atom_part[i,2]*100
        # New line should add longitudinal_signed values
        old = """            print(f'{i+1:2d}    {w_cm[i]:8.2f}    {w_meV[i]:7.2f}     {long_char[i]:5.3f}   '
                  f'{ixs_s_str:9s}  {ixs_as_str:9s}  '
                  f'{atom_part[i,0]*100:4.1f}  {atom_part[i,1]*100:4.1f}  {atom_part[i,2]*100:4.1f}  '
                  f'{pol_type[i]}')
"""
        new = """            print(f'{i+1:2d}    {w_cm[i]:8.2f}    {w_meV[i]:7.2f}     {long_char[i]:5.3f}   '
                  f'{ixs_s_str:9s}  {ixs_as_str:9s}  '
                  f'{atom_part[i,0]*100:4.1f} {long_sign[i,0]:+5.1f}  '
                  f'{atom_part[i,1]*100:4.1f} {long_sign[i,1]:+5.1f}  '
                  f'{atom_part[i,2]*100:4.1f} {long_sign[i,2]:+5.1f}  '
                  f'{pol_type[i]}')
"""
        lines[i] = new
        print(f"Updated line {i+1}")

# Also need to add long_sign = result['longitudinal_signed'] around line 317
for i, line in enumerate(lines):
    if 'atom_part = result' in line:
        # Add the longitudinal_signed extraction right after
        lines.insert(i+1, "        long_sign = result['longitudinal_signed']\n")
        print(f"Added long_sign extraction at line {i+2}")
        break

# Write back
with open('code/single_q_analysis.py', 'w') as f:
    f.writelines(lines)

print("\n✓ Updated print section to show longitudinal components")
