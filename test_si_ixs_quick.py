"""Quick test: Can we run IXS on Si?"""
import numpy as np
import sys
sys.path.insert(0, 'code')

from force_constants import ForceConstants
from phonons import calc_Dq, calc_freq_eig, convert_frequencies
from form_factors import CalcAtomicfQ
from ixs import calc_ixs

xtal = ForceConstants("data/Test__Silicon_dispersion/Qgrid_444/Cg.fc")
Phi = xtal.convert_to_eV_per_Angstrom2()
masses = np.array([xtal.masses[0], xtal.masses[0]])  # Both Si

Q_prim = np.array([[0.5, 0, 0]])
Dq = calc_Dq(Q_prim, xtal.uvw, Phi, masses)
w_raw, ev = calc_freq_eig(Dq, sort_modes=True)
w_THz = convert_frequencies(w_raw, 'THz')
w_cm = convert_frequencies(w_raw, 'cm-1')

# Form factors (both atoms are Si)
Q_cart = Q_prim @ xtal.b_l
Q_mag = np.linalg.norm(Q_cart)
fSi = CalcAtomicfQ(Q_mag, 'Si', scale=4*np.pi, use_xraylib=False)
fQ = np.array([[fSi, fSi]])

kT_THz = 6.2
Is, Ias, n, F, info = calc_ixs(w_THz, ev, Q_prim, xtal.b_l, xtal.xs, 
                                fQ, masses, kT_THz, units='barn')

print("IXS at X=(0.5,0,0):")
print(f"  Q_mag = {Q_mag:.4f} (2π/Å)")
print(f"  f(Si) = {fSi:.2f}")
print(f"\n  Mode   ω(cm⁻¹)    IXS(S)    IXS(AS)")
for i in range(6):
    print(f"  {i+1:4d}  {w_cm[0,i]:8.2f}  {Is[0,i]:9.3f}  {Ias[0,i]:9.3f}")
