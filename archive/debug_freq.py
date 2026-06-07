with open('code/single_q_analysis.py', 'r') as f:
    content = f.read()

# Add debug print right before freq_label is set
old = """        # Choose frequency unit
        freq_label = {'meV': 'meV', 'cm-1': 'cm⁻¹', 'THz': 'THz'}[freq_unit]"""

new = """        # Choose frequency unit
        print(f"DEBUG: freq_unit = {freq_unit}")  # DEBUG
        freq_label = {'meV': 'meV', 'cm-1': 'cm⁻¹', 'THz': 'THz'}[freq_unit]"""

content = content.replace(old, new)

with open('code/single_q_analysis.py', 'w') as f:
    f.write(content)

print("✓ Added debug print")
