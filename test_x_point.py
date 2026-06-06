import sys
sys.path.insert(0, 'code')
from force_constants import ForceConstants
from phonons import calc_Dq, calc_freq_eig, convert_frequencies
from fcc_structure import fcc_conv2prim_k
import numpy as np

fc_file = "data/Test__Silicon_dispersion/Qgrid_888/Cg.fc"
xtal = ForceConstants(fc_file)
Phi = xtal.convert_to_eV_per_Angstrom2()
masses = np.array([xtal.masses[xtal.atom_type_map[i]-1] for i in range(xtal.nat)])

# Try different X point definitions
x_points_conv = {
    'X1: (0, 0.5, 0.5)': np.array([0.0, 0.5, 0.5]),
    'X2: (1, 0, 0)': np.array([1.0, 0.0, 0.0]),
    'X3: (0.5, 0, 0)': np.array([0.5, 0.0, 0.0]),
}

print("Testing different X point definitions:")
print("=" * 60)
for name, x_conv in x_points_conv.items():
    x_prim = fcc_conv2prim_k(x_conv)
    q_path = x_prim.reshape(1, 3)
    D = calc_Dq(q_path, xtal.uvw, Phi, masses)
    w_raw, ev = calc_freq_eig(D)
    w_meV = convert_frequencies(w_raw, 'meV')
    print(f"\n{name}")
    print(f"  Conv: {x_conv}, Prim: {x_prim}")
    print(f"  Frequencies (meV): {w_meV[0]}")
