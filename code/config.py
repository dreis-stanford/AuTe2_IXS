"""
Configuration file for AuTe2 IXS project
"""

import os
from typing import Dict, Tuple

# ============================================================================
# PATHS
# ============================================================================

# Path to sixcircle code (Alfred's six-circle diffractometer code)
# Update this to point to where sixcircle is installed on your system
SIXCIRCLE_PATH = os.path.expanduser("~/Documents/MyPython/Others/sixcircle_1p85")

# Alternative: set to None if sixcircle is not available
# SIXCIRCLE_PATH = None

# Check if sixcircle exists
SIXCIRCLE_AVAILABLE = False
if SIXCIRCLE_PATH and os.path.exists(SIXCIRCLE_PATH):
    SIXCIRCLE_AVAILABLE = True
else:
    print(f"Note: Running in simulation mode (sixcircle not found)")
    print(f"      Looked at: {SIXCIRCLE_PATH}")
    print(f"      Update SIXCIRCLE_PATH in code/config.py if needed")

# ============================================================================
# EXPERIMENTAL PARAMETERS
# ============================================================================

# Beamline configuration
BEAMLINE = 43  # BL43LXU at SPring-8
SI_ORDER = 11  # Si(11,11,11) analyzer
INCIDENT_TYPE = 1  # 1 = M3 standard, 3 = KB micro-focus
FROZEN_ANGLES = "456"  # Freeze mu, gam, omega

# Analyzer slit gaps (micrometers)
DEFAULT_AGAP_V = 60  # Vertical gap
DEFAULT_AGAP_H = 60  # Horizontal gap

# ============================================================================
# SAMPLE PARAMETERS - AuTe2
# ============================================================================

# Lattice parameters for AuTe2 (to be refined from experiment)
# These are initial values - will be updated after orientation matrix refinement
LATTICE_PARAMS = {
    'a': 8.15,      # Angstroms
    'b': 4.47,
    'c': 5.08,
    'alpha': 90.0,  # degrees
    'beta': 90.0,
    'gamma': 90.0
}

# Modulation vector (in r.l.u. of average structure)
# This is the CDW/superstructure modulation
MODULATION_VECTOR = (0.28, 0.0, 0.0)  # (qh, qk, ql)

# Sample name/description
SAMPLE_NAME = "AuTe2_CDW"

# ============================================================================
# PHONON CALCULATION PARAMETERS
# ============================================================================

# DFT phonon data location
DFT_DATA_DIR = "data/dft"

# Energy unit for phonon frequencies
FREQ_UNIT = "meV"  # Options: "meV", "THz", "cm-1"
# Q-point coordinate system
Q_COORDINATE_SYSTEM = "conventional"  # Options: "primitive", "conventional", "cartesian"

# Q-point mesh for phonon interpolation
Q_MESH = (50, 50, 50)

# ============================================================================
# IXS CROSS-SECTION PARAMETERS
# ============================================================================

# Temperature for Debye-Waller factor (Kelvin)
TEMPERATURE = 300.0

# Energy resolution (meV) - typical for BL43LXU
ENERGY_RESOLUTION = 1.5

# Include Debye-Waller factor
USE_DEBYE_WALLER = False # Not yet implemented - TODO
# ============================================================================
# RECIPROCAL SPACE EXPLORATION
# ============================================================================

# Maximum |Q| to consider (nm^-1)
MAX_Q = 60.0

# Maximum satellite order to include
MAX_SATELLITE_ORDER = 2

# High-symmetry points for AuTe2 (orthorhombic)
HIGH_SYMMETRY_POINTS = {
    'Γ': (0.0, 0.0, 0.0),
    'X': (0.5, 0.0, 0.0),
    'Y': (0.0, 0.5, 0.0),
    'Z': (0.0, 0.0, 0.5),
    'S': (0.5, 0.5, 0.0),
    'T': (0.0, 0.5, 0.5),
    'U': (0.5, 0.0, 0.5),
    'R': (0.5, 0.5, 0.5),
}

# ============================================================================
# OUTPUT/PLOTTING
# ============================================================================

# Output directories
OUTPUT_DIR = "results"
FIGURE_DIR = "results/figures"
DATA_DIR = "data"

# Plotting defaults
PLOT_DPI = 300
PLOT_FORMAT = 'png'  # or 'pdf'

# ============================================================================
# HELPER FUNCTIONS
# ============================================================================

def get_lattice_parameters() -> Dict[str, float]:
    """Get current lattice parameters"""
    return LATTICE_PARAMS.copy()

def update_lattice_parameters(**kwargs):
    """
    Update lattice parameters
    
    Example:
        update_lattice_parameters(a=8.16, b=4.48)
    """
    for key, value in kwargs.items():
        if key in LATTICE_PARAMS:
            LATTICE_PARAMS[key] = value
        else:
            raise ValueError(f"Unknown lattice parameter: {key}")

def get_modulation_vector() -> Tuple[float, float, float]:
    """Get modulation vector"""
    return MODULATION_VECTOR

def update_modulation_vector(qh: float, qk: float, ql: float):
    """Update modulation vector"""
    global MODULATION_VECTOR
    MODULATION_VECTOR = (qh, qk, ql)

def print_config():
    """Print current configuration"""
    print("\n" + "="*70)
    print("AuTe2 IXS Project Configuration")
    print("="*70)
    
    print("\n[Paths]")
    print(f"  Sixcircle: {SIXCIRCLE_PATH}")
    print(f"  Available: {SIXCIRCLE_AVAILABLE}")
    
    print("\n[Experiment - BL{0}LXU]".format(BEAMLINE))
    print(f"  Si order: ({SI_ORDER}, {SI_ORDER}, {SI_ORDER})")
    print(f"  Incident type: {INCIDENT_TYPE}")
    print(f"  Frozen angles: {FROZEN_ANGLES}")
    print(f"  Analyzer gaps: V={DEFAULT_AGAP_V}μm, H={DEFAULT_AGAP_H}μm")
    
    print("\n[Sample - {0}]".format(SAMPLE_NAME))
    print(f"  a = {LATTICE_PARAMS['a']:.4f} Å")
    print(f"  b = {LATTICE_PARAMS['b']:.4f} Å")
    print(f"  c = {LATTICE_PARAMS['c']:.4f} Å")
    print(f"  α = {LATTICE_PARAMS['alpha']:.2f}°")
    print(f"  β = {LATTICE_PARAMS['beta']:.2f}°")
    print(f"  γ = {LATTICE_PARAMS['gamma']:.2f}°")
    
    print("\n[Modulation]")
    print(f"  q_mod = ({MODULATION_VECTOR[0]:.3f}, {MODULATION_VECTOR[1]:.3f}, {MODULATION_VECTOR[2]:.3f})")
    
    print("\n[Phonon Calculations]")
    print(f"  Frequency unit: {FREQ_UNIT}")
    print(f"  Q coordinate system: {Q_COORDINATE_SYSTEM}")
    print(f"  Q-mesh: {Q_MESH}")
    print(f"  Temperature: {TEMPERATURE} K")
    print(f"  Debye-Waller: {USE_DEBYE_WALLER} (not yet implemented)")    
    print("="*70 + "\n")

# Print config on import for visibility
if __name__ != "__main__":
    # Only print brief status on import
    status = "AVAILABLE" if SIXCIRCLE_AVAILABLE else "SIMULATION MODE"
    print(f"AuTe2 IXS Config loaded - Sixcircle: {status}")
