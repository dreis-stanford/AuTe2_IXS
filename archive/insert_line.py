with open('code/single_q_analysis_si.py', 'r') as f:
    lines = f.readlines()

# Find the line with Q_sinThOverLambda and add form_factors right after
for i, line in enumerate(lines):
    if "result['Q_sinThOverLambda'] = Q_sinThOverLambda" in line:
        # Check if next line is NOT form_factors
        if i+1 < len(lines) and 'form_factors' not in lines[i+1]:
            lines.insert(i+1, '        result["form_factors"] = {"Si": fSi}\n')
            print(f"Inserted form_factors at line {i+2}")
        else:
            print(f"form_factors already exists at line {i+2}")
        break
else:
    print("Could not find Q_sinThOverLambda line!")

with open('code/single_q_analysis_si.py', 'w') as f:
    f.writelines(lines)
