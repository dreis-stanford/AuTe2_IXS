# Read file
with open('code/single_q_analysis_si.py', 'r') as f:
    lines = f.readlines()

# Find line with Q_sinThOverLambda result and add form_factors after it
for i, line in enumerate(lines):
    if 'result["Q_sinThOverLambda"] = Q_sinThOverLambda' in line or "result['Q_sinThOverLambda'] = Q_sinThOverLambda" in line:
        # Insert form_factors line after
        if 'form_factors' not in lines[i+1]:
            lines.insert(i+1, '        result["form_factors"] = {"Si": fSi}\n')
        print(f"Added form_factors at line {i+2}")
        break

# Write back
with open('code/single_q_analysis_si.py', 'w') as f:
    f.writelines(lines)

print("✓ Fixed")
