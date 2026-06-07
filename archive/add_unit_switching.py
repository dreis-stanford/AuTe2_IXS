with open('code/single_q_analysis.py', 'r') as f:
    content = f.read()

# Find the interactive_mode function
# Add freq_unit variable and switching logic

# 1. Add freq_unit initialization after coord_system initialization
old_init = """    # Default coordinate system
    coord_system = 'primitive'"""

new_init = """    # Default coordinate system and frequency unit
    coord_system = 'primitive'
    freq_unit = 'meV'  # Options: 'meV', 'cm-1', 'THz'"""

content = content.replace(old_init, new_init)

# 2. Add unit switching to the command processing
old_switch = """            if user_input.lower() == 'prim':
                coord_system = 'primitive'
                print("  → Switched to primitive coordinates")
                continue
            elif user_input.lower() == 'conv':
                coord_system = 'conventional'
                print("  → Switched to conventional coordinates")
                continue
            elif user_input.lower() == 'cart':
                coord_system = 'cartesian'
                print("  → Switched to Cartesian coordinates")
                continue"""

new_switch = """            if user_input.lower() == 'prim':
                coord_system = 'primitive'
                print("  → Switched to primitive coordinates")
                continue
            elif user_input.lower() == 'conv':
                coord_system = 'conventional'
                print("  → Switched to conventional coordinates")
                continue
            elif user_input.lower() == 'cart':
                coord_system = 'cartesian'
                print("  → Switched to Cartesian coordinates")
                continue
            elif user_input.lower() == 'mev':
                freq_unit = 'meV'
                print("  → Switched to meV")
                continue
            elif user_input.lower() in ['cm', 'cm-1', 'invcm']:
                freq_unit = 'cm-1'
                print("  → Switched to cm⁻¹")
                continue
            elif user_input.lower() == 'thz':
                freq_unit = 'THz'
                print("  → Switched to THz")
                continue"""

content = content.replace(old_switch, new_switch)

# 3. Pass freq_unit to analyze function
old_call = """            result = analyzer.analyze(Q_input, coords=coord_system,
                                     print_results=True, print_detailed=False)"""

new_call = """            result = analyzer.analyze(Q_input, coords=coord_system,
                                     freq_unit=freq_unit,
                                     print_results=True, print_detailed=False)"""

content = content.replace(old_call, new_call)

# 4. Update the prompt
old_prompt = """    print("\\nCommands: 'prim', 'conv', 'cart' to switch coordinates")
    print("          'q' or Ctrl+C to quit")"""

new_prompt = """    print("\\nCommands: 'prim', 'conv', 'cart' - switch coordinates")
    print("          'meV', 'THz', 'invcm' - switch frequency units")
    print("          'q' or Ctrl+C to quit")"""

content = content.replace(old_prompt, new_prompt)

# 5. Update analyze function signature
old_sig = """    def analyze(self, Q_input, coords='primitive', print_results=True, 
                print_detailed=False, threshold_L=0.7, threshold_T=0.3):"""

new_sig = """    def analyze(self, Q_input, coords='primitive', freq_unit='meV',
                print_results=True, print_detailed=False, 
                threshold_L=0.7, threshold_T=0.3):"""

content = content.replace(old_sig, new_sig)

# 6. Pass freq_unit to print function
old_print_call = """        if print_results:
            self._print_results(result, print_detailed)"""

new_print_call = """        if print_results:
            self._print_results(result, print_detailed, freq_unit)"""

content = content.replace(old_print_call, new_print_call)

# 7. Update _print_results signature and implementation
old_print_sig = """    def _print_results(self, result, print_detailed=False):"""
new_print_sig = """    def _print_results(self, result, print_detailed=False, freq_unit='meV'):"""

content = content.replace(old_print_sig, new_print_sig)

# 8. Update the frequency display in the table
old_freq_header = """        print('Mode  Freq(cm⁻¹)  Freq(meV)   L-char   IXS(S)    IXS(AS)   Au%   Au‖    Te1%  Te1‖   Te2%  Te2‖   Pol')"""

new_freq_header = """        # Choose which frequency to display
        freq_label = {'meV': 'meV', 'cm-1': 'cm⁻¹', 'THz': 'THz'}[freq_unit]
        freq_data = {'meV': result['frequencies_meV'], 
                     'cm-1': result['frequencies_cm'],
                     'THz': result['frequencies_THz']}[freq_unit]
        
        print(f'Mode  Freq({freq_label:>4s})   L-char   IXS(S)    IXS(AS)   Au%   Au‖    Te1%  Te1‖   Te2%  Te2‖   Pol')"""

content = content.replace(old_freq_header, new_freq_header)

# 9. Update the print loop to use freq_data
old_loop_print = """            print(f'{i+1:2d}    {w_cm[i]:8.2f}    {w_meV[i]:7.2f}     {long_char[i]:5.3f}   '"""

new_loop_print = """            print(f'{i+1:2d}    {freq_data[i]:8.2f}     {long_char[i]:5.3f}   '"""

content = content.replace(old_loop_print, new_loop_print)

with open('code/single_q_analysis.py', 'w') as f:
    f.write(content)

print("✓ Added frequency unit switching")
print("  Commands: meV, THz, invcm")
