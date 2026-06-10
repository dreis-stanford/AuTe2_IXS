# Known Issues

### `angles` command in simulation mode still doesn't produce angles
**Open, 2026-06-10.** The `angles` interactive command no longer crashes with
`AttributeError: 'NoneType' object has no attribute 'ca'` (fixed: it now
dispatches through `SixCircleInterface.move_to_hkl(hkl, check_only=True)`,
which calls `_simulate_angles`/`_setup_simulation_UB` in simulation mode).

However, `_setup_simulation_UB` (`code/sixcircle_interface.py`) itself appears
never to have been exercised end-to-end and has at least two bugs:

1. **`sv_radi` NameError on `import sixcircle_rqd`** — the external module's
   `init_rqd()` looks for `BL43XU_CONST.mac` in the cwd; without it, it falls
   into a broken `setincident(9)` fallback. *Worked around* by restoring
   `BL43XU_CONST.mac` to the project root (gitignored — it was moved to
   `archive/beamline_config/` in `5b2773e`). Side effect: this import also
   writes `previous.bl.conf` to the cwd (also gitignored).
2. **`g_sample` NameError / blocking `input()` in `sixcircle.setlat()`**
   (`_setup_simulation_UB` line ~342) — `_setup_simulation_UB` never sets
   `sixcircle.g_sample` (unlike `_load_sixcircle`, which sets it to
   `'AuTe2_exp'`), and calling `setlat()` with 6 args triggers
   `input('\nSample description...')`, which references the unset global and
   would also block on stdin. Needs `sixcircle.g_sample = ...` set first, and
   a 7th arg passed to `setlat()` to skip the prompt.

The same `setlat()` 6-arg blocking-`input()` pattern also affects the *live*
path: `SixCircleInterface.setup_experiment()` (`code/sixcircle_interface.py`
line ~146) calls `self.sixc.setlat(...)` with 6 args and hangs/`EOFError`s on
non-interactive runs (confirmed via `code/test_sixcircle_integration.py`,
pre-existing — not introduced by recent changes). Same fix needed: pass a 7th
`g_sample` arg.

Likely more issues beyond this in the `or0`/`or1` orientation setup. Tracked
under TODO #2 ("Fix Sixcircle/YAML Integration").

## Resolved

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
