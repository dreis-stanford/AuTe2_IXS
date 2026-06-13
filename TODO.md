# AuTe2 IXS Analysis - TODO List

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
- [ ] Default 3D view: orient the camera so the lab z-axis points mostly
      "up" on screen (currently an arbitrary default elevation/azimuth).
- [ ] Double-check the lab-frame axis conventions (origin/sense of the
      th and tth rotations) against the sixcircle documentation
      (`~/Documents/MyPython/Others/sixcircle_1p85/sixcircle_documentation/`)
      to confirm `verify_scattering.rotate_Q_to_lab_frame` matches the real
      instrument, not just an internally-consistent convention.

## Low Priority

### 7. Features
- [ ] Add Debye-Waller factor (currently disabled in config)
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
