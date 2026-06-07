with open('code/single_q_analysis.py', 'r') as f:
    content = f.read()

# Fix the Q_dot_e calculation - need to ensure it's a scalar
old_line = "                longitudinal_signed[imode, iat] = 100 * Q_dot_e / (total_ev_norm + 1e-12)"
new_line = "                longitudinal_signed[imode, iat] = 100 * np.real(Q_dot_e) / (total_ev_norm + 1e-12)"

content = content.replace(old_line, new_line)

# Also make sure Q_dot_e calculation returns scalar
old_calc = "                Q_dot_e = np.real(np.dot(Q_hat, e_atom))"
new_calc = "                Q_dot_e = float(np.real(np.dot(Q_hat, e_atom)))"

content = content.replace(old_calc, new_calc)

with open('code/single_q_analysis.py', 'w') as f:
    f.write(content)

print("✓ Fixed Q_dot_e scalar conversion")
