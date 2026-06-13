# AuTe2 IXS Analysis - TODO List

## Next Session Priorities (set 2026-06-13, end of day)

With the analyzer-array visualization (performance, wider-region grid,
linear-saturating IXS scale, linear-radius marker sizes) and the
geometry-viz UB refactor + top-down default view all done, the suggested
order is:

1. **Measurement planning** (#4) and **wiring up `code/q_optimizer.py`**
   (currently a standalone, unwired module - see README "Project
   Structure") - use it to help choose/rank Q-points for the measurement
   plan.
2. Verify the lab-frame axis conventions
   (`verify_scattering.rotate_Q_to_lab_frame`, sense of th/chi/phi
   rotations) against the sixcircle documentation (see #6a, second bullet).

## High Priority

### 1. Add Polarization Factor
- [ ] Implement cos²(2θ) polarization correction for σ→σ scattering
- [ ] Account for horizontal synchrotron polarization
- [ ] Add to IXS cross-section calculations
- [ ] Document effect on intensity predictions

### 2. Fix Sixcircle/YAML Integration
- [x] Fix the `angles` command so it produces real diffractometer angles
      (2026-06-12: load the library in all modes, compute via `ca_s()`;
      deleted the broken `_setup_simulation_UB`/`_simulate_angles`). See
      KNOWN_ISSUES.md "Resolved".
- [x] SPEC-like frozen mode/values and angle limits, configurable in
      `config.py` and editable at runtime via `freeze`/`limits` commands.
- [ ] Install yaml module (`pip install pyyaml`)
- [x] Integrate `ca6()` for accurate analyzer positions
      (2026-06-12: `SixCircleInterface.analyzer_array_offsets()` computes
      per-analyzer (dH,dK,dL) via the real UB matrix + BL43LXU 7x4 array
      geometry, matching ca6()'s math/naming convention)
- [x] Replace hardcoded analyzer offsets with real ca6() output
- [ ] `setup_experiment()` still uses a 6-arg `setlat()` that blocks on
      `input()` (live setup path; not used by the interactive tools) —
      see KNOWN_ISSUES.md.

### 3. Analyzer Array Improvements
- [x] Use actual ca6() output for Q-dependent analyzer positions

## Medium Priority

### 4. Measurement Planning
- [ ] Create measurement plan script for AuTe2
- [ ] Optimize Q-points for b-axis (CDW direction)
- [ ] Include satellite positions
- [ ] Calculate diffractometer accessibility
- [ ] Estimate measurement times

### 5. Documentation
- [ ] Add measurement strategy guide to README
- [ ] Document polarization effects
- [ ] Add examples of satellite analysis
- [ ] Create tutorial for analyzer array usage

### 6. Validation
- [ ] Compare with MATLAB results (if available)
- [ ] Verify satellite phonon calculations
- [ ] Cross-check IXS intensities

### 6a. Geometry Visualization (`viz` / `code/geometry_viz.py`)
- [x] Real-UB crystal frame (2026-06-12, `a1a5e14`): `SixCircleInterface.get_UB()`
      plus a `ub=`/`--ub` arg to `plot_scattering_geometry`/the `geometry_viz`
      CLI; surface normal = c* = `UB @ hkl`, real axes a,b,c = columns of
      `inv(UB).T`, both rotated to lab via `rotate_Q_to_lab_frame`. Idealized
      axes remain as a fallback when `ub=None`.
- [x] Default 3D view: view of the (horizontal) scattering plane oriented
      to show the front (mounted) crystal face
      (2026-06-13: `ax.view_init(elev=-90, azim=-90, roll=-90)`, superseding
      the earlier `elev=90, azim=-90, roll=90` top-down view from the same
      day). Incident beam (k_in, mu=0) points right (travels left-to-right);
      k_out points right-and-down for tth>0 (screen_x=-y_lab,
      screen_y=-x_lab, -z_lab toward viewer). The exposed crystal face (+c*)
      faces the viewer whenever n_lab_z<0 -- true across chi roughly
      -88..+75 deg for the real UB (checked at 7 representative reflections),
      so the front face stays visible across this range. Bulk/side faces are
      now opaque (alpha 0.30 -> 0.95) so they clearly sit behind the front
      face. The legend was removed; k_in/k_out/Q are now labeled directly on
      the plot (color-matched), and the surface-normal label shows its HKL,
      e.g. `n (001)`. The lab-frame reference triad moved to the new empty
      corner. Also fixed a z-ordering bug (`ax.computed_zorder = False`)
      where mplot3d's automatic depth-sorting could place k_in/k_out/Q
      behind the now-opaque crystal faces; draw order now guarantees they
      render on top.
- [ ] Double-check the lab-frame axis conventions (origin/sense of the
      th and tth rotations) against the sixcircle documentation
      (`~/Documents/MyPython/Others/sixcircle_1p85/sixcircle_documentation/`)
      to confirm `verify_scattering.rotate_Q_to_lab_frame` matches the real
      instrument, not just an internally-consistent convention. (Partial
      check 2026-06-13: the top-down view's k_in/k_out/tth geometry is
      self-consistent with `_kin_hat`/`_kout_hat`'s lab-frame definitions --
      no inconsistency found there -- but th/chi/phi rotation sense vs. the
      real instrument is still unverified.)

## Low Priority

### 7. Features
- [ ] Export results to file (CSV, JSON)
- [ ] Batch analysis mode (multiple Q-points from file)
- [ ] Plotting functions for dispersion

### 8. Code Quality
- [ ] Add unit tests for coordinate transformations
- [ ] Test RLV finding with edge cases
- [ ] Verify all satellite orders work correctly
- [ ] Add error handling for invalid inputs

## Completed ✓

- [x] Fix reciprocal lattice vector finding for monoclinic structures
- [x] Add satellite reflection support (H K L m input)
- [x] Suppress verbose sixcircle output
- [x] Clean up repository (archive temporary scripts)
- [x] Add analyzer array analysis
- [x] Add DFT modulation limitation note
- [x] Update documentation
- [x] Add Debye-Waller factor (2026-06-13: `code/debye_waller.py`,
      `config.DEBYE_WALLER_MODE` = 'none' | 'phonon' (BZ-summed anisotropic
      U from these force constants at `config.TEMPERATURE`) | 'cif' (fixed
      298K ADPs from Reithmayer et al., Acta Cryst. B49, 6 1993; AuTe2 only).
      Default is now 'cif'.)
- [x] Full energy-dependent form factors (2026-06-13: `calc_form_factor()` in
      `code/form_factors.py` returns f(Q,E) = f0(Q) + f'(E) + i*f''(E) via
      xraylib when available, wired into `single_q_analysis.py`. f''
      sign-convention vs. this codebase's exp(-2*pi*i*Q.r) phase convention
      verified against the optical theorem; see `tests/test_form_factors.py`.)
- [x] Analyzer-array visualization improvements (2026-06-13):
      - Performance: avoid redundant angle solves in the `array` command
        (~25x speedup, `cfb2d76`).
      - Wider-region contour background via `analyze_array_grid`/
        `SixCircleInterface.compute_angle_grid`, vectorized (`b5a9f9a`).
      - Explicit plot range (`dtth_range=(-3,3)`, `dgam_range=(-3,1)`
        degrees), linear IXS scale saturating at the 10th/90th percentiles,
        decluttered per-subplot axes with a single plot-range annotation in
        the figure title (`b6737f0`).
      - IXS marker size now linear in radius (not area), min marker area
        ~3 pt^2 / max ~215 pt^2, for better small-vs-large contrast.
- [x] Geometry-viz UB refactor and top-down default view (see #6a).

## Notes

**Polarization Factor:**
- For horizontal synchrotron polarization and σ→σ scattering
- Factor ≈ cos²(2θ) for modes in scattering plane
- Modes perpendicular to scattering plane suppressed
- Important for b-polarized CDW modes

**YAML Issue:**
- Currently running in simulation mode
- Need `pyyaml` for full sixcircle functionality
- Affects ca6() availability and diffractometer control

**CDW Strategy:**
- Modulation along b-axis, q_mod in (a*, c*) plane
- b-polarized modes at satellites are key
- Need to orient sample for optimal polarization selection
