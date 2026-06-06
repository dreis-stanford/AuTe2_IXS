with open('code/single_q_analysis.py', 'r') as f:
    content = f.read()

# Fix the calculation - handle complex arrays properly
old_calc = "                Q_dot_e = float(np.real(np.dot(Q_hat, e_atom)))"
new_calc = "                Q_dot_e = np.real(np.vdot(Q_hat, e_atom))  # vdot handles complex"

content = content.replace(old_calc, new_calc)

old_assign = "                longitudinal_signed[imode, iat] = 100 * np.real(Q_dot_e) / (total_ev_norm + 1e-12)"
new_assign = "                longitudinal_signed[imode, iat] = 100 * Q_dot_e / (total_ev_norm + 1e-12)"

content = content.replace(old_assign, new_assign)

with open('code/single_q_analysis.py', 'w') as f:
    f.write(content)

print("✓ Fixed complex number handling")
