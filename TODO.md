# AuTe2 IXS Analysis - TODO List

## High Priority

### 1. Add Polarization Factor
- [ ] Implement cos²(2θ) polarization correction for σ→σ scattering
- [ ] Account for horizontal synchrotron polarization
- [ ] Add to IXS cross-section calculations
- [ ] Document effect on intensity predictions

### 2. Fix Sixcircle/YAML Integration
- [ ] Install yaml module (`pip install pyyaml`)
- [ ] Test sixcircle loading in non-simulation mode
- [ ] Integrate `ca6()` for accurate analyzer positions
- [ ] Replace hardcoded analyzer offsets with real ca6() output
- [ ] Fix `_setup_simulation_UB` (`code/sixcircle_interface.py`) so the
      `angles` command produces real simulated angles — see
      KNOWN_ISSUES.md for the `g_sample`/`setlat()` bug and others
      likely lurking in the `or0`/`or1` setup

### 3. Analyzer Array Improvements
- [ ] Use actual ca6() output for Q-dependent analyzer positions
- [ ] Add coordinate system handling (ensure conventional units)
- [ ] Option to control number of modes displayed
- [ ] Option to set analyzer gaps interactively

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
