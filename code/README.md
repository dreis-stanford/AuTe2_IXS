# AuTe2 IXS Analysis Code

## Main Analysis Scripts

**Interactive Analysis:**
- `single_q_analysis.py` - Single Q-point phonon and IXS analysis (use via `../analyze_q.py`)
- `single_q_analysis_si.py` - Silicon version (use via `../analyze_si.py`)

## Core Modules

**Phonon Calculations:**
- `force_constants.py` - Load force constants from DFT
- `phonons.py` - Calculate dynamical matrix, frequencies, eigenvectors

**IXS Cross-sections:**
- `ixs.py` - IXS cross-section calculation (production code)
- `form_factors.py` - Atomic form factors f0(Q), and full f(Q,E) =
  f0(Q) + f'(E) + i*f''(E) via xraylib (`calc_form_factor`)
- `debye_waller.py` - Anisotropic Debye-Waller U tensors (phonon BZ sum or
  literature ADPs; see `config.DEBYE_WALLER_MODE`)

**Crystal Structure:**
- `aute2_structure.py` - AuTe2 coordinate transformations (C2/m monoclinic)
- `fcc_structure.py` - FCC (Si) coordinate transformations
- `fcc_path.py` - FCC high-symmetry paths
- `modulated_structure.py` - Modulated structures

**Experimental Control:**
- `sixcircle_interface.py` - Diffractometer control interface
- `sixcircle_wrapper.py` - Wrapper functions
- `sixcircle_minimal.py` - Minimal implementation

**Utilities:**
- `constants.py` - Physical constants
- `config.py` - Configuration and structure loading

## Test/Debug Scripts

- `test_*.py` - Unit tests
- `check_*.py` - Verification scripts
- `debug_*.py` - Debugging tools
- `compare_*.py` - Comparison with MATLAB
- `verify_*.py` - Verification tools

## Plotting

- `plot_dispersion.py` - General dispersion plotting
- `plot_dispersion_aute2_working.py` - AuTe2 specific
- `plot_dispersion_si.py` - Silicon specific
