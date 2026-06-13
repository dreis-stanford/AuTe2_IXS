# Known Issues

### `setup_experiment()` setlat() blocking `input()` (live setup path only)
**Open, 2026-06-12.** `SixCircleInterface.setup_experiment()`
(`code/sixcircle_interface.py`) calls `self.sixc.setlat(...)` with 6 args,
which triggers `input('\nSample description...')` and hangs/`EOFError`s on
non-interactive runs. This path is **not** used by the interactive tools — only
by `code/test_sixcircle_integration.py` — so it doesn't affect `analyze_q.py`.
Fix when needed: pass a 7th `g_sample` arg to `setlat()` (the same blocking
pattern, now avoided in `_load_sixcircle`, which sets globals directly and
never calls the 6-arg `setlat`).

## Resolved

### `angles` command didn't produce angles
**Resolved 2026-06-12.** The interactive `angles` command now produces real
diffractometer angles in both modes. Root cause was a broken parallel
`_setup_simulation_UB`/`_simulate_angles` path that re-imported the sixcircle
library badly (no cwd config-file guard → `KeyError`/`sv_radi` NameError),
never set `g_sample` or the angle limits, and read a nonexistent
`sixcircle.A[...]` array (`ca()` only sets `outcaTTH/TH/CHI/PHI`).

Fix: since the sixcircle package is pure calculation (`mv()`/`br()` only
update module globals — no hardware I/O), the interface now loads the library
whenever it's available (even in "simulation mode") and computes angles via
`ca_s()`, which returns all six angles directly. The broken
`_setup_simulation_UB`/`_simulate_angles` were deleted. Frozen mode/values and
angle limits are now driven by `config.py` (`FROZEN_ANGLES`, `FROZEN_VALUES`,
`ANGLE_LIMITS`) and editable at runtime via the new `freeze`/`limits` commands.
Regression tests: `tests/test_sixcircle_angles.py`.

### Longitudinal Character Calculation - Silicon [100] Direction
**Resolved 2026-06-09.** Root cause: per-atom eigenvector extraction used
`ev_mode.reshape(3, nat)[:, iat]`, but the DOF layout is atom-major
`[x1,y1,z1, x2,y2,z2, ...]`, so x/y/z components were scrambled across atoms.
Fixed to `ev_mode.reshape(nat, 3)[iat, :]` in `single_q_analysis.py` and
`single_q_analysis_si.py` (L-char, atomic participation, and Q·e table).
Also replaced `np.abs(np.dot(e, e))` with `np.real(np.vdot(e, e))` for complex
eigenvectors. Si along [100]-type directions now gives clean L=1.0 / T=0.0.
Regression test: `tests/test_physics_regressions.py::test_si_lchar_clean_split`

### sed-corrupted constants in single_q_analysis.py
**Resolved 2026-06-09.** `update_separators.sh` globally replaced `* 100` with
`* 80`, corrupting the kT cm⁻¹→THz conversion (T was 20% low) and the Q·e
percent scale. Constants restored; the script now refuses to run.
Regression test: `tests/test_physics_regressions.py::test_kt_conversion_in_source`

### Mass unit conversion in force_constants.py
**Resolved 2026-06-09.** QE `.fc` masses are in Rydberg atomic units (2·mₑ);
the parser divided by M_u (931.494) instead of amu/(2mₑ) (911.444), making all
masses 2.15% low and all frequencies ~1.09% high. Now converts via
`m_val * 2*Me / M_u`. Si → 28.0855 amu, Au → 196.97, Te → 127.60.
Regression test: `tests/test_physics_regressions.py::test_si_mass_is_amu`
