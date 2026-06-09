# Code Review — AuTe2 IXS Project
*Reviewed 2026-06-09. Three physics bugs verified numerically; details below.*

## Verified physics bugs (high priority)

### 1. sed script corrupted physics constants (`update_separators.sh`)
The script's global replacements (`s/\* *100/* 80/g` etc.) were meant for separator
lines but also hit numeric constants in `code/single_q_analysis.py`:

- **Lines 632 & 901:** `kT_THz = kT_cm * const.c * 80 / 1e12` — was `* 100`
  (m/s → cm/s). Temperature is 20% low (≈238 K instead of 298 K), so all Bose
  factors and Stokes/anti-Stokes ratios in the AuTe2 interactive analysis are wrong.
- **Line 234:** `longitudinal_signed = 80 * Q_dot_e / ...` — was `100 *` (percent scale).

Backups (`*.backup_20260608_084103`) confirm both were `100`. The Si version and
`ixs.py` still have the correct `* 100`. **Also fix or delete `update_separators.sh`
— rerunning it will re-corrupt the file.**

### 2. Eigenvector reshape transposed → broken L/T character (the known Si [100] bug)
DOF layout is atom-major `[x1,y1,z1, x2,y2,z2, ...]` (see `calc_Dq`'s
`np.repeat(masses, 3)` and `ixs.py`'s `ev[3*s:3*s+3]`). But the L-char and
atomic-participation code extracts per-atom vectors with:

```python
ev_reshaped = ev_mode.reshape(3, nat)   # WRONG: scrambles components
ev_atom = ev_reshaped[:, iat]           # gives [x1, z1, y2] for nat=2 !
```

Correct: `ev_mode.reshape(nat, 3)[iat, :]`.

Affected: `single_q_analysis.py` (~lines 161, 205, 530) and
`single_q_analysis_si.py` (lines 137, 185).

**Numerical verification** (Si, q = [0.1,0,0] primitive, numpy-only re-derivation):

| mode | ω (THz) | L-char current | L-char fixed |
|------|---------|----------------|--------------|
| TA×2 | 1.41    | 0.11 / 0.33    | 0.000        |
| LA   | 2.87    | 0.54           | **1.000**    |
| TO×2 | 15.13   | 0.33 / 0.56    | 0.000        |
| LO   | 15.13   | 0.12           | **1.000**    |

The fix gives the physically required clean L/T split; the current code gives the
fractional values described in KNOWN_ISSUES.md. This resolves that issue.

Related, same code blocks: `np.abs(np.dot(e, e))` on complex vectors computes
|Σeᵢ²|, not Σ|eᵢ|². Use `np.real(np.vdot(e, e))`. Masked at q-points where
`rotate_ev` makes eigenvectors nearly real, wrong in general (complex ev at
general q in a monoclinic cell).

### 3. Mass unit conversion off by 2.2% → all frequencies ~1.09% high
`force_constants.py` line 105: `self.masses[i] = m_val / const.M_u * 1e6`
(divides by 931.494). Quantum-Espresso `.fc` masses are in **Rydberg atomic units**
(unit = 2mₑ), so the correct divisor is amu/(2mₑ) = **911.444**.

Verification against the files: Si comes out as 28.085 amu with 911.444 (matches
to 10 digits); current code yields Si = 28.703, Au = 192.73, Te = 124.85 amu —
all 2.15% low. Uniform mass scaling doesn't change eigenvectors, but every
frequency is multiplied by √(931.494/911.444) → **+1.09% systematic error**, which
matters when comparing to measured IXS positions at 1.5 meV resolution.

## Other defects

- **Duplicate code in `single_q_analysis.py`:** `_print_array_results` defined twice
  (second definition silently wins); the `elif user_input.lower() == 'array'` branch
  duplicated (second is dead code).
- **`_format_xs(val)`** in the AuTe2 class lacks `@staticmethod` (Si version has it).
  Latent crash if ever called via `self`.
- **Satellite input bug (interactive):** `Q = Q_main + m*q_mod` adds the modulation
  vector (defined in conventional r.l.u.) to the typed Q regardless of the current
  coordinate system — wrong result if user is in `prim` or `cart` mode.
- **`angles` command in simulation mode:** calls `sixc.sixc.ca(...)` directly;
  `sixc.sixc` is `None` in simulation mode → AttributeError instead of the intended
  simulated-angle path (`move_to_hkl`/`_simulate_angles`).
- **Inconsistent L-char convention between materials:** AuTe2 version projects onto
  reduced q̂ (correct for phonon polarization); Si version uses full Q̂ including G —
  differs whenever G ≠ 0.
- **Phase-convention mismatch:** `ixs.py` uses exp(−2πi Q·r) for the one-phonon
  structure factor; the elastic |F|² in `_print_results` uses exp(+2πi Q·r).
  Harmless for intensities, but should be unified before adding interference terms.
- **`config.py` side effects at import:** prints banners and auto-loads
  `data/AuTe_2_m.fc` via a cwd-relative path — imports fail/misbehave when run from
  any other directory. Same cwd assumption in `_load_sixcircle`'s config-file
  renaming (a crash mid-load would leave `ini.conf.bak` behind).
- **Broken/orphaned modules:** `sixcircle_wrapper.py` references
  `config.PROJECT_ROOT` (doesn't exist); it and `sixcircle_minimal.py` use absolute
  `import config`, which fails inside the package. Decide: fix or archive.
- **Package name `code` shadows the Python stdlib `code` module** (used by pdb and
  interactive tooling). Worth renaming (e.g. `aute2_ixs`) before the package grows.
- **Hardcoded approximations:** analyzer-array offsets in `analyze_array` (TODO
  acknowledges; replace with `ca6()` output), dummy or0/or1 orientation angles in
  `sixcircle_interface._load_sixcircle`.

## Repo hygiene

~30 one-off `fix_*/debug_*/check_*` scripts at top level (beyond `archive/`);
10 backup variants of `single_q_analysis.py` inside `code/`; `fcc_structure.py.broken`,
`config.py.backup`, `*.bak`, `*~` files; `venv/` and `AuTe2_IXS.tar.gz` in the tree;
`tests/`, `docs/`, `notebooks/` directories empty while test scripts sit at top level.
Suggest: archive the one-offs, gitignore backups, move tests into `tests/` as pytest.

## What's in good shape

The computational core is clean and well-documented: vectorized `calc_Dq` with
Hermiticity enforcement, acoustic sum rule in the FC parser, Hungarian mode
tracking, both Cromer-Mann and xraylib form factors with explicit s = sinθ/λ
conventions, correct use of reduced q for dynamics and full Q for the IXS structure
factor, the two-structure (DFT vs experimental) design in config, and a whitelisted
command pass-through for sixcircle. The monoclinic conv↔prim transforms and the FCC
versions have built-in self-tests that pass.

## Suggested order of work

1. Fix the three verified bugs (sed constants, reshape+vdot, mass units); delete or
   constrain `update_separators.sh`.
2. Deduplicate `single_q_analysis.py`; fix satellite coordinate handling and the
   simulation-mode `angles` path.
3. Add a pytest suite locking in: Γ-point acoustic ω≈0, Si L-char = 1/0 along [100],
   Si frequencies vs known dispersion, conv↔prim round-trips, IXS Stokes/anti-Stokes
   ratio = (n+1)/n at known T.
4. Repo cleanup (archive one-offs, remove backups from `code/`).
5. Then the TODO list: polarization factor cos²(2θ), real `ca6()` analyzer positions,
   Debye-Waller.
