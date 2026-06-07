with open('code/single_q_analysis.py', 'r') as f:
    content = f.read()

# Add debug at the point where freq_unit is changed
content = content.replace(
    """        elif user_input.lower() == 'thz':
            freq_unit = 'THz'
            print("  → Switched to THz\\n")
            continue""",
    """        elif user_input.lower() == 'thz':
            freq_unit = 'THz'
            print(f"  → Switched to THz (DEBUG: freq_unit now = {freq_unit})\\n")
            continue"""
)

# Add debug right before the analyze call
content = content.replace(
    """            result = analyzer.analyze(Q_input, coords=input_coords,
                                     freq_unit=freq_unit,""",
    """            print(f"DEBUG: About to call analyze with freq_unit={freq_unit}")
            result = analyzer.analyze(Q_input, coords=input_coords,
                                     freq_unit=freq_unit,"""
)

with open('code/single_q_analysis.py', 'w') as f:
    f.write(content)

print("✓ Added more debug")
