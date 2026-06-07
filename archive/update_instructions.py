with open('code/single_q_analysis.py', 'r') as f:
    content = f.read()

# Find and update the instructions section
old_instructions = """Instructions:
  - Enter Q vector as three numbers (e.g., 0.5 0 0)
  - Press Enter on empty line to quit
  - Type 'conv', 'prim', or 'cart' to change coordinate system
  - Type 'meV', 'THz', or 'invcm' to change frequency units
  - Current system: {coord_system}, units: {freq_unit}"""

new_instructions = """Instructions:
  - Enter Q vector as three numbers (e.g., 0.5 0 0)
  - Press Enter on empty line to quit
  
Commands:
  - Coordinates: 'prim', 'conv', 'cart'
  - Freq units:  'meV', 'THz', 'invcm'
  
Current settings: {coord_system} coordinates, {freq_unit} units"""

content = content.replace(old_instructions, new_instructions)

with open('code/single_q_analysis.py', 'w') as f:
    f.write(content)

print("✓ Updated instructions")
