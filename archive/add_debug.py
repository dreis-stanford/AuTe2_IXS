with open('code/single_q_analysis_si.py', 'r') as f:
    lines = f.readlines()

# Add debug print right before _print_results is called
for i, line in enumerate(lines):
    if "self._print_results(result, detailed=print_detailed, freq_unit=freq_unit)" in line:
        lines.insert(i, '            print(f"DEBUG: result keys = {list(result.keys())}")  # DEBUG\n')
        print(f"Added debug at line {i+1}")
        break

with open('code/single_q_analysis_si.py', 'w') as f:
    f.writelines(lines)
