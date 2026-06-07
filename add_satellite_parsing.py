#!/usr/bin/env python3
"""
Update interactive_mode() to parse satellite order (4th parameter)
"""

filepath = 'code/single_q_analysis.py'

with open(filepath, 'r') as f:
    content = f.read()

# Find the Q parsing section in interactive_mode
old_parsing = """        try:
            Q = np.array([float(x) for x in user_input.split()])
            
            if len(Q) != 3:
                print('⚠ Error: Please enter exactly 3 numbers\\n')
                continue"""

new_parsing = """        try:
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
                print(f'             = [{Q[0]:.4f}, {Q[1]:.4f}, {Q[2]:.4f}]\\n')
            else:
                print('⚠ Error: Enter 3 numbers (H K L) or 4 numbers (H K L m) for satellites\\n')
                continue"""

if old_parsing in content:
    content = content.replace(old_parsing, new_parsing)
    print("✓ Updated Q parsing to support satellites")
else:
    print("⚠ Could not find Q parsing section - may need manual update")

# Update instructions
old_instructions = """    print("Instructions:")
    print("  - Enter Q vector as three numbers (e.g., 0.5 0 0)")
    print("  - Press Enter on empty line to quit")"""

new_instructions = """    print("Instructions:")
    print("  - Enter Q vector as three numbers (e.g., 0.5 0 0)")
    print("  - For satellites, add order m: H K L m (e.g., 1 0 0 1 for first satellite)")
    print("  - Press Enter on empty line to quit")"""

content = content.replace(old_instructions, new_instructions)

# Write back
with open(filepath, 'w') as f:
    f.write(content)

print("✓ Updated interactive mode")
print("\nUsage examples:")
print("  1 0 0       → Main Bragg peak (1,0,0)")
print("  1 0 0 1     → First satellite: (1,0,0) + q_mod")
print("  1 0 0 -1    → First satellite: (1,0,0) - q_mod")
print("  2 0 0 2     → Second satellite: (2,0,0) + 2*q_mod")
