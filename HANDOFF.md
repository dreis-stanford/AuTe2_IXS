# Session Handoff — 2026-06-09

## State
5 local commits today (eb12234..471522b), not yet pushed. All 14 tests in
tests/ pass (`pytest tests/ -v` from project root). See CODE_REVIEW.md for
the full review and KNOWN_ISSUES.md for resolved bugs.

## Sixcircle geometry sanity check — RESOLVED
Wrote check_sixcircle_geometry.py (uses scbasic.py from
~/Documents/MyPython/Others/sixcircle_1p85/ directly, no sixcircle module
side effects). Result: OVERALL PASS, no rotation-sense/lab-frame convention
mismatch.

- Test 1 (Si cubic, trivial mounting a||x b||y c||z): U == diag(1,-1,-1) as
  expected from beam-along--y geometry; ca() round-trips and verify_full
  passes exactly for (4,0,0), (0,4,0), (0,0,4).
- Test 2 (AuTe2, or0=(0,0,1) chi=90 / or1=(0,2,0) phi=90 from
  _load_sixcircle): or0/or1 round-trip tth/chi correctly; UB matrix is
  self-consistent (U has ~1.5° off-diagonal tilt, expected for a real
  mount).
- Found and fixed one real bug in code/verify_scattering.py: verify_full
  used k_in = k*[0,-1,0] (ignoring mu) while Q_lab included the R_mu
  rotation, so for AuTe2's mu=-0.1719° every reflection showed a spurious
  ~0.033 Å⁻¹ "direction error" (WARNING). Fixed k_in to
  k*[0,-cos(mu),-sin(mu)] and k_out_expected to include cos(gam)/sin(gam);
  all AuTe2 direction errors now <0.005 Å⁻¹ (PASS). 14/14 tests still pass.
- Conclusion: the hardcoded or0/or1 angles in
  code/sixcircle_interface.py _load_sixcircle are correct and consistent
  with scbasic.py conventions — no changes needed there.

Remaining minor item (not blocking): or0 (0,0,1) BASIC test (|Q| from
lattice vs from angle) is a WARNING at 3.55e-4 rel. error, from tth=6.45°
being rounded to 2 decimals vs exact Bragg angle 6.448°. Cosmetic only.

Fixed and verified today:
- sed-corrupted constants (kT * 80 -> * 100); update_separators.sh retired
- Eigenvector reshape bug (atom-major layout) -> L-char, participation,
  Q.e table all correct now; verified vs QE matdyn.f90 conventions and
  C2 symmetry (Au exactly 0/1 for q || b*)
- Mass units (Ry a.u. -> amu, divisor 911.444): all frequencies were 1.09% high
- Si L-char now uses reduced q (G-invariant)
- Unified classification: L > 0.95, T < 0.05 (strict, selection rules);
  M(L) > 0.7, M(T) < 0.3; else M. Constants STRICT_L etc. in
  code/single_q_analysis.py
- Si analyzer = thin wrapper over generic SingleQAnalyzer (one code path)

## Next task (agreed priority)
**Sixcircle geometry sanity check**: hand-built orientation matrix in
code/sixcircle_interface.py (_load_sixcircle has hardcoded dummy or0/or1
angles); verify ca(h,k,l) angles match expected beamline values. David
suspects rotation-sense / lab-frame convention mismatch between his
verify_scattering.py (code/verify_scattering.py) and sixcircle. Resources:
~/Documents/MyPython/Others/sixcircle_1p85/ has axis-definition PDFs,
math writeup (sixcircle_documentation/), and scbasic.py (rotation
matrices). Strategy: derive conventions from PDFs + scbasic.py, test with
cubic crystal + trivial mounting first, then AuTe2. Note: 'angles' command
crashes in simulation mode (sixc.sixc is None) - known, unfixed.

## Deferred (from CODE_REVIEW.md)
- Repo cleanup: ~30 one-off fix_*/debug_* scripts at top level, 10 backup
  copies of single_q_analysis.py in code/, broken sixcircle_wrapper.py /
  sixcircle_minimal.py (use absent config.PROJECT_ROOT / absolute imports)
- config.py side effects at import (cwd-relative fc path)
- Package dir 'code' shadows Python stdlib module
- TODO.md high-priority items: polarization factor cos^2(2theta), real
  ca6() analyzer positions (replace hardcoded offsets in analyze_array),
  Debye-Waller
- Old plots (aute2_dispersion.png etc.) predate the 1.09% frequency fix

## Environment notes
- Sandbox has no scipy/pytest; tests/conftest.py stubs scipy with numpy
  (works for these tests). On David's machine, real scipy/pytest are used.
- Stale git lock files can appear after sandbox commits; if git complains:
  rm .git/index.lock .git/HEAD.lock
