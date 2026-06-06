with open('code/single_q_analysis.py', 'r') as f:
    lines = f.readlines()

# Find where coord_system is initialized in interactive_mode
for i, line in enumerate(lines):
    if "coord_system = 'primitive'" in line and 'interactive_mode' in ''.join(lines[max(0,i-20):i]):
        # Check if freq_unit is already on next line
        if i+1 < len(lines) and 'freq_unit' not in lines[i+1]:
            lines.insert(i+1, "    freq_unit = 'meV'\n")
            print(f"Added freq_unit initialization at line {i+2}")
            break
        elif 'freq_unit' in lines[i+1]:
            print("freq_unit already initialized")
            break

with open('code/single_q_analysis.py', 'w') as f:
    f.writelines(lines)

print("✓ Fixed freq_unit initialization")
