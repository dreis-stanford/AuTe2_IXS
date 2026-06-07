with open('code/single_q_analysis.py', 'r') as f:
    content = f.read()

# Remove debug prints
content = content.replace(
    '        print(f"DEBUG: freq_unit = {freq_unit}")  # DEBUG\n',
    ''
)
content = content.replace(
    '            print(f"DEBUG: About to call analyze with freq_unit={freq_unit}")\n',
    ''
)
content = content.replace(
    '  → Switched to THz (DEBUG: freq_unit now = {freq_unit})',
    '  → Switched to THz'
)

with open('code/single_q_analysis.py', 'w') as f:
    f.write(content)

print("✓ Removed debug statements")
