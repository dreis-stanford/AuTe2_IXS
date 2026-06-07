with open('code/single_q_analysis.py', 'r') as f:
    lines = f.readlines()

# Find line with "coord_system = 'primitive'" in interactive_mode
for i, line in enumerate(lines):
    if i > 344 and "coord_system = 'primitive'" in line:  # After interactive_mode def
        # Add freq_unit on next line
        lines.insert(i+1, "    freq_unit = 'meV'\n")
        print(f"Added freq_unit = 'meV' at line {i+2}")
        break

with open('code/single_q_analysis.py', 'w') as f:
    f.writelines(lines)

print("✓ Added freq_unit initialization")
