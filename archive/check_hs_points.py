import sys
sys.path.insert(0, 'code')
from fcc_structure import fcc_conv2prim_k
import numpy as np

# Standard FCC high-symmetry points in CONVENTIONAL coordinates
points_conv = {
    'G': np.array([0.0, 0.0, 0.0]),
    'X': np.array([0.0, 0.5, 0.5]),     
    'L': np.array([0.5, 0.5, 0.5]),
    'K': np.array([0.375, 0.75, 0.375]),
    'W': np.array([0.25, 0.75, 0.5]),
}

print("High-symmetry points:")
print(f"{'Point':<5} {'Conventional':<20} {'Primitive':<20}")
print("-" * 50)
for name, pt in points_conv.items():
    pt_prim = fcc_conv2prim_k(pt)
    print(f"{name:<5} {str(pt):<20} {str(pt_prim):<20}")
