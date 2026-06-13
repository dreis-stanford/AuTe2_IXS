# AuTe2 IXS Experiment Planning

Python tools for planning Inelastic X-ray Scattering (IXS) experiments on AuTe₂ using DFT phonon predictions and six-circle diffractometer control.

## Overview

This project integrates:
- DFT phonon calculations to predict phonon dispersions and eigenvectors
- Sixcircle diffractometer interface for experiment planning at SPring-8 
- Modulated structure handling for CDW satellite reflections
- Q-point optimization module to help maximize information from measurements
  (implemented as a standalone class, not yet wired into the interactive tools)

## Quick Start

Interactive Analysis:
    python analyze_q.py
    python analyze_si.py

## Installation

Requirements: numpy scipy matplotlib

For SPring-8 users, edit code/config.py to set SIXCIRCLE_PATH

## Project Structure

- code/config.py - Configuration
- code/sixcircle_interface.py - Diffractometer (just calc. no control)  
- code/modulated_structure.py - Satellite (CDW) reflection handling, used by
  the `H K L m` input below
- code/debye_waller.py - Anisotropic Debye-Waller U tensors, from a phonon
  BZ sum (`'phonon'` mode) or fixed literature ADPs (`'cif'` mode); see
  `config.DEBYE_WALLER_MODE`
- code/form_factors.py - Atomic form factors f(Q), and full f(Q,E) =
  f0(Q) + f'(E) + i*f''(E) via xraylib (`calc_form_factor`)
- code/q_optimizer.py - Q-point optimization (standalone module, not yet
  wired into the interactive tools)
- analyze_q.py - AuTe2 launcher
- test_integration.py - Sixcircle integration tests; requires the external
  sixcircle package configured for a real instrument (not a general
  installation smoke test - use test_installation.py / pytest tests/ for that)

## Materials

AuTe2: Monoclinic C2/m, 3 atoms/cell, 9 phonon branches
Silicon: FCC diamond, 2 atoms/cell, 6 phonon branches

## Satellite Reflection Analysis

AuTe₂ has a CDW modulation with q_mod = (-0.4076, 0, 0.4479) at 298K.

Access satellites using 4-parameter input:
H K L m
where m is the satellite order.
Examples:
- `1 0 0` → Main Bragg peak
- `1 0 0 1` → First satellite: (1,0,0) + q_mod
- `1 0 0 -1` → First satellite: (1,0,0) - q_mod  
- `2 0 0 2` → Second order: (2,0,0) + 2×q_mod
Reference: Schutte et al., Acta Cryst. B44, 486 (1988)

## Interactive Commands
- `H K L [m]` - Analyze Q-point (optional m for satellites)
- `angles` - Calculate diffractometer angles for the last Q (calculation only;
  also prints a ready-to-paste SPEC `mv` command)
- `freeze` - Show/set the frozen-angle mode and held values, e.g.
  `freeze 456` (freeze mu, gam, omega) or `freeze mu=-0.17`
- `limits` - Show/set angle limits, e.g. `limits tth 0 60` (selects the
  positive-tth solution branch) or `limits reset`
- `conv`, `prim`, `cart` - Switch coordinate systems
- `meV`, `THz`, `invcm` - Switch frequency units
- `sixc <cmd>` - Pass command to sixcircle

Note: this is a planning tool — it never drives the diffractometer. Frozen
mode/values and angle limits start from `code/config.py` (`FROZEN_ANGLES`,
`FROZEN_VALUES`, `ANGLE_LIMITS`) and can be changed at runtime with the
`freeze`/`limits` commands above.

## Recent Updates

### Debye-Waller factor and dispersion-corrected form factors (2026-06-13)
- `config.DEBYE_WALLER_MODE` ('none' | 'phonon' | 'cif', default **'cif'**)
  selects an anisotropic Debye-Waller correction applied to the IXS and
  elastic structure factors:
  - `'phonon'`: U tensors from a Monkhorst-Pack BZ sum over these force
    constants at `config.TEMPERATURE` (`config.DEBYE_WALLER_QMESH`).
  - `'cif'`: fixed 298K ADPs for AuTe2 from Reithmayer et al., Acta Cryst.
    B49, 6 (1993); only applied when `material == 'AuTe2'`.
- `code/form_factors.calc_form_factor()` now returns the full energy-dependent
  form factor f(Q,E) = f0(Q) + f'(E) + i*f''(E) via xraylib (falls back to
  the real f0(Q) from Cromer-Mann coefficients if xraylib is unavailable).
  The sign of xraylib's `Fii` (anomalous f'') was verified against the
  optical theorem and against this codebase's exp(-2*pi*i*Q.r) phase
  convention (`code/ixs.py`) — the two sign flips cancel, giving the
  physically correct (absorptive) |F(Q)|^2. See the docstring of
  `calc_form_factor` and `tests/test_form_factors.py`.
- Fixed a stale-state issue where `pa()`/startup banners reported
  "Si(11 11 11)" (from a cached `previous.bl.conf`) even though calculations
  correctly used the configured `config.SI_ORDER` (12); the loader now calls
  `setorder(config.SI_ORDER)` so the cached state and banners stay in sync.

### Sixcircle integration
- Diffractometer interface for SPring-8 BL43LXU (1,184+ lines)
- Modulated structure and satellite handling
- Q-point optimization (standalone framework, see Project Structure)
- Simulation mode for testing
- Added or0/or1 orientation reflections used by sixcircle for angle calculations

### Output formatting improvements
- Standardized all separator lines to exactly 80 characters for consistent terminal display
- Reorganized output information for better readability:
  - Temperature, form factors, and scattering parameters in compact header 
  - IXS units and calculation notes moved to footer
- Enhanced phonon mode table with:
  - Clear polarization indicators (L=Longitudinal, T=Transverse, M=Mixed)
  - Eigenvector projections showing Q·e magnitude and relative phase for each atom
  - IXS cross-sections for both Stokes and anti-Stokes processes (not including polarization or DW)
  - Structure factor information for elastic scattering (not including pol. or DW)

### Physics fixes (2026-06-09, see CODE_REVIEW.md / KNOWN_ISSUES.md)
- Fixed sed-corrupted constants (kT scale factor, longitudinal_signed) in single_q_analysis.py
- Fixed eigenvector reshape (atom-major layout) - L-char, atomic participation,
  and Q·e tables are now correct (previously gave fractional/incorrect values)
- Fixed mass-unit conversion (Ry a.u. -> amu): all phonon frequencies were ~1.09% high
- Verified sixcircle geometry conventions and orientation matrix against scbasic.py

