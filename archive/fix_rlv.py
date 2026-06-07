#!/usr/bin/env python3
"""
Fix reciprocal lattice vector finding in single_q_analysis.py
"""

import shutil
from datetime import datetime

# Read the file
with open('code/single_q_analysis.py', 'r') as f:
    content = f.read()

# Create backup
backup = f'code/single_q_analysis.py.backup_{datetime.now().strftime("%Y%m%d_%H%M%S")}'
shutil.copy2('code/single_q_analysis.py', backup)
print(f"✓ Backup created: {backup}")

# Find and replace the old method
old_code = """        # Find nearest reciprocal lattice vector
        Q_prim_rounded = np.round(Q_prim)
        Q_prim_reduced = Q_prim - Q_prim_rounded
        
        Q_conv_G = aute2_prim2conv_k(Q_prim_rounded)
        Q_conv_reduced = aute2_prim2conv_k(Q_prim_reduced)
        
        Q_cart_G = Q_prim_rounded @ self.xtal.b_l
        Q_cart_reduced = Q_prim_reduced @ self.xtal.b_l
        Q_mag_reduced = np.linalg.norm(Q_cart_reduced)
        
        result['Q_reduced_prim'] = Q_prim_reduced
        result['Q_reduced_conv'] = Q_conv_reduced
        result['Q_reduced_cart'] = Q_cart_reduced
        result['G_prim'] = Q_prim_rounded
        result['G_conv'] = Q_conv_G
        result['G_cart'] = Q_cart_G"""

new_code = """        # Find nearest reciprocal lattice vector
        # Search in Cartesian space to handle non-orthogonal reciprocal lattice
        Q_prim_rounded = np.round(Q_prim)
        Q_cart = Q_prim @ self.xtal.b_l
        
        # Search nearby integer combinations
        min_dist = np.inf
        best_G = Q_prim_rounded
        
        for dh in [-1, 0, 1]:
            for dk in [-1, 0, 1]:
                for dl in [-1, 0, 1]:
                    G_test = Q_prim_rounded + np.array([dh, dk, dl])
                    G_cart = G_test @ self.xtal.b_l
                    dist = np.linalg.norm(Q_cart - G_cart)
                    if dist < min_dist:
                        min_dist = dist
                        best_G = G_test
        
        G_prim = best_G
        Q_prim_reduced = Q_prim - G_prim
        
        Q_conv_G = aute2_prim2conv_k(G_prim)
        Q_conv_reduced = aute2_prim2conv_k(Q_prim_reduced)
        
        Q_cart_G = G_prim @ self.xtal.b_l
        Q_cart_reduced = Q_prim_reduced @ self.xtal.b_l
        Q_mag_reduced = np.linalg.norm(Q_cart_reduced)
        
        result['Q_reduced_prim'] = Q_prim_reduced
        result['Q_reduced_conv'] = Q_conv_reduced
        result['Q_reduced_cart'] = Q_cart_reduced
        result['G_prim'] = G_prim
        result['G_conv'] = Q_conv_G
        result['G_cart'] = Q_cart_G"""

# Replace
content = content.replace(old_code, new_code)

# Write back
with open('code/single_q_analysis.py', 'w') as f:
    f.write(content)

print("✓ Fixed RLV finding code")
print("\nDone! Test with: python analyze_q.py")
