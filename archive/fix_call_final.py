with open('code/single_q_analysis.py', 'r') as f:
    content = f.read()

# Fix the call to _print_results
content = content.replace(
    "self._print_results(result, print_detailed)",
    "self._print_results(result, detailed=print_detailed, freq_unit=freq_unit)"
)

with open('code/single_q_analysis.py', 'w') as f:
    f.write(content)

print("✓ Fixed _print_results call to include freq_unit")
