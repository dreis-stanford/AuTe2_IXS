with open('code/single_q_analysis.py', 'r') as f:
    content = f.read()

# Find the interactive loop and add unit/coord switching before the Q parsing
old_section = """        try:
            Q = np.array([float(x) for x in user_input.split()])
            
            if len(Q) != 3:
                print('⚠ Error: Please enter exactly 3 numbers\\n')
                continue"""

new_section = """        # Check for command keywords first
        if user_input.lower() == 'prim':
            coord_system = 'primitive'
            print("  → Switched to primitive coordinates\\n")
            continue
        elif user_input.lower() == 'conv':
            coord_system = 'conventional'
            print("  → Switched to conventional coordinates\\n")
            continue
        elif user_input.lower() == 'cart':
            coord_system = 'cartesian'
            print("  → Switched to Cartesian coordinates\\n")
            continue
        elif user_input.lower() == 'mev':
            freq_unit = 'meV'
            print("  → Switched to meV\\n")
            continue
        elif user_input.lower() in ['cm', 'cm-1', 'invcm']:
            freq_unit = 'cm-1'
            print("  → Switched to cm⁻¹\\n")
            continue
        elif user_input.lower() == 'thz':
            freq_unit = 'THz'
            print("  → Switched to THz\\n")
            continue
        
        try:
            Q = np.array([float(x) for x in user_input.split()])
            
            if len(Q) != 3:
                print('⚠ Error: Please enter exactly 3 numbers\\n')
                continue"""

content = content.replace(old_section, new_section)

# Add freq_unit initialization in interactive_mode
old_init = """    # Set coordinate system
    coord_system = 'primitive'
    print(f'  - Current system: {coord_system}')"""

new_init = """    # Set coordinate system and frequency unit
    coord_system = 'primitive'
    freq_unit = 'meV'
    print(f'  - Current system: {coord_system}, units: {freq_unit}')"""

content = content.replace(old_init, new_init)

# Pass freq_unit to analyze
old_analyze = """            result = analyzer.analyze(Q_input, coords=input_coords, 
                                     print_results=True, print_detailed=False)"""

new_analyze = """            result = analyzer.analyze(Q_input, coords=input_coords,
                                     freq_unit=freq_unit,
                                     print_results=True, print_detailed=False)"""

content = content.replace(old_analyze, new_analyze)

# Update instructions
old_instr = """    print('  - Type \\'conv\\', \\'prim\\', or \\'cart\\' to change coordinate system')"""
new_instr = """    print('  - Type \\'conv\\', \\'prim\\', or \\'cart\\' to change coordinate system')
    print('  - Type \\'meV\\', \\'THz\\', or \\'invcm\\' to change frequency units')"""

content = content.replace(old_instr, new_instr)

with open('code/single_q_analysis.py', 'w') as f:
    f.write(content)

print("✓ Added interactive unit switching")
