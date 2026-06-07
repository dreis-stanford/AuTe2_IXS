with open('code/single_q_analysis.py', 'r') as f:
    lines = f.readlines()

# Find and replace the instructions section
for i, line in enumerate(lines):
    if 'print("Instructions:")' in line or "print('Instructions:')" in line:
        # Replace the next few lines
        lines[i] = '    print("Instructions:")\n'
        lines[i+1] = '    print("  - Enter Q vector as three numbers (e.g., 0.5 0 0)"))\n'
        lines[i+2] = '    print("  - Press Enter on empty line to quit")\n'
        lines[i+3] = '    print("")\n'
        lines[i+4] = '    print("Commands:")\n'
        lines.insert(i+5, '    print("  - Coordinates: \'prim\', \'conv\', \'cart\'")\n')
        lines.insert(i+6, '    print("  - Freq units:  \'meV\', \'THz\', \'invcm\'")\n')
        lines.insert(i+7, '    print("")\n')
        # Update the Current system line
        for j in range(i+8, min(i+12, len(lines))):
            if 'Current system:' in lines[j]:
                lines[j] = f'    print(f"Current settings: {{coord_system}} coordinates, {{freq_unit}} units")\n'
                break
        break

with open('code/single_q_analysis.py', 'w') as f:
    f.writelines(lines)

print("✓ Fixed instructions")
