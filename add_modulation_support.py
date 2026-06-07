#!/usr/bin/env python3
"""
Add modulation vector support to single_q_analysis.py
"""

import shutil
from datetime import datetime

# Backup
backup = f'code/single_q_analysis.py.backup_{datetime.now().strftime("%Y%m%d_%H%M%S")}'
shutil.copy2('code/single_q_analysis.py', backup)
print(f"✓ Backup: {backup}")

# Read file
with open('code/single_q_analysis.py', 'r') as f:
    lines = f.readlines()

# Find the imports section and add modulated_structure
for i, line in enumerate(lines):
    if 'from .aute2_structure import' in line:
        lines.insert(i+1, 'from .modulated_structure import ModulatedStructure\n')
        print("✓ Added modulated_structure import")
        break

# Find __init__ and add modulation vector
for i, line in enumerate(lines):
    if 'def __init__(self, xtal, Phi, masses, kT_THz):' in line:
        for j in range(i, min(i+20, len(lines))):
            if 'self.kT_THz = kT_THz' in lines[j]:
                lines.insert(j+1, '        self.mod_struct = ModulatedStructure()\n')
                print("✓ Added ModulatedStructure to __init__")
                break
        break

# Write back
with open('code/single_q_analysis.py', 'w') as f:
    f.writelines(lines)

print("✓ Done")
