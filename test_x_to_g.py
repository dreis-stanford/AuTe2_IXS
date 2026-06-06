import sys
sys.path.insert(0, 'code')
from fcc_structure import fcc_conv2prim_k
import numpy as np

# X point in conventional: (1, 1, 0) or equivalently (0, 0, 1)
X_conv_1 = np.array([1.0, 1.0, 0.0])
X_conv_2 = np.array([0.0, 0.0, 1.0])

# Convert to primitive
X_prim_1 = fcc_conv2prim_k(X_conv_1)
X_prim_2 = fcc_conv2prim_k(X_conv_2)

print("X point equivalents:")
print(f"  (1,1,0) conv -> {X_prim_1} prim")
print(f"  (0,0,1) conv -> {X_prim_2} prim")

# Path from (0,0,1) to (0,0,0) in conventional
print("\nX->G path in conventional (0,0,xi) with xi: 1->0")
for xi in [1.0, 0.75, 0.5, 0.25, 0.0]:
    q_conv = np.array([0.0, 0.0, xi])
    q_prim = fcc_conv2prim_k(q_conv)
    print(f"  (0, 0, {xi:.2f}) conv -> {q_prim} prim")
