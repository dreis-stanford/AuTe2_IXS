with open('code/single_q_analysis.py', 'r') as f:
    content = f.read()

# Fix the function signature
content = content.replace(
    'def _print_results(self, result, detailed=False):',
    'def _print_results(self, result, detailed=False, freq_unit="meV"):'
)

# Fix the call
content = content.replace(
    'self._print_results(result, print_detailed, freq_unit)',
    'self._print_results(result, detailed=print_detailed, freq_unit=freq_unit)'
)

with open('code/single_q_analysis.py', 'w') as f:
    f.write(content)

print("✓ Fixed function call")
