with open('code/single_q_analysis.py', 'r') as f:
    lines = f.readlines()

# Find and remove the first command parsing block (lines ~397-410)
new_lines = []
skip = False
skip_count = 0

for i, line in enumerate(lines):
    # Start skipping at "Check for coordinate system change"
    if "# Check for coordinate system change" in line:
        skip = True
        skip_count = 0
        continue
    
    # Stop skipping after the cartesian block
    if skip:
        skip_count += 1
        if skip_count > 15 or "# Try to parse Q vector" in line:
            skip = False
            # Don't include the comment line either
            if "# Try to parse Q vector" in line:
                continue
    
    if not skip:
        new_lines.append(line)

with open('code/single_q_analysis.py', 'w') as f:
    f.writelines(new_lines)

print("✓ Removed duplicate command parsing")
