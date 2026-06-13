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
SIXCIRCLE_AVAILABLE = bool(SIXCIRCLE_PATH and os.path.exists(SIXCIRCLE_PATH))

# ============================================================================
# EXPERIMENTAL PARAMETERS
# ============================================================================

# Beamline configuration
BEAMLINE = 43  # BL43LXU at SPring-8
SI_ORDER = 12  # Si(12,12,12) analyzer
INCIDENT_TYPE = 1  # 1 = M3 standard, 3 = KB micro-focus

# Nominal Si lattice constant at 30°C (sixcircle's setorder() convention)
SI_LATTICE_CONSTANT = 5.431109  # Angstroms

# Incident wavelength/energy, set by the Si(n,n,n) analyzer backscattering
# condition lambda = 2*a_Si / (n*sqrt(3)) -- same formula as sixcircle_rqd's
# setorder(SI_ORDER). Single source of truth: sixcircle_interface.py and
# modulated_structure.py both read this instead of hardcoding lambda.
WAVELENGTH = 2 * SI_LATTICE_CONSTANT / (SI_ORDER * 3 ** 0.5)  # Angstroms
ENERGY_KEV = 12.39854 / WAVELENGTH

# ----------------------------------------------------------------------------
# Diffractometer modes, constraints, and limits (SPEC-like)
# ----------------------------------------------------------------------------
# NB: this layer is calculation-only. Nothing here moves the diffractometer;
# it only steers which angle solution sixcircle's ca()/ca6() returns.

# Which three angles are frozen, as a 3-digit code (sixcircle convention):
#   tth=0  th=1  chi=2  phi=3  mu=4  gam=5  omega=6  azimuth=7  alpha=8  beta=9
# "456" freezes mu, gam, omega (the usual IXS planning mode).
FROZEN_ANGLES = "456"

# Values (degrees) at which frozen angles are held. Only the angles named in
# FROZEN_ANGLES are actually applied; the rest are ignored. mu is held at the
# BL43LXU incident tilt by default.
FROZEN_VALUES = {
    'tth': 0.0, 'th': 0.0, 'chi': 0.0, 'phi': 0.0,
    'mu': -0.1719, 'gam': 0.0, 'omega': 0.0,
    'azimuth': 0.0, 'alpha': 0.0, 'beta': 0.0,
}

# Angle limits (degrees) as (lower, upper). ca()/ca6() reject solutions outside
# these, so narrowing them selects among multiple solutions (e.g. force tth>0).
ANGLE_LIMITS = {
    # tth defaults to >= 0 so the conventional top-face (tth>0) branch is
    # selected; sixcircle otherwise returns the negative-tth mirror solution
    # (bottom-face geometry) for reflections like (0,0,L>0). Widen with e.g.
    # `limits tth -180 180` to inspect the bottom-face branch.
    'tth':   (0.0, 180.0),
    'th':    (-180.0, 180.0),
    'chi':   (-180.0, 180.0),
    'phi':   (-180.0, 180.0),
    'mu':    (-180.0, 180.0),
    'gam':   (-180.0, 180.0),
    'alpha': (-180.0, 180.0),
    'beta':  (-180.0, 180.0),
}

# Surface normal / azimuthal reference, as an HKL direction. This is the
# sixcircle azimuth reference (g_haz, g_kaz, g_laz); the incidence/exit angles
# alpha and beta are measured against it. (0,0,1) is the c* reciprocal
# direction (the mounted face for AuTe2).
SURFACE_NORMAL = (0, 0, 1)

# Analyzer slit gaps (micrometers)
DEFAULT_AGAP_V = 60  # Vertical gap
DEFAULT_AGAP_H = 60  # Horizontal gap

# ============================================================================
# SIXCIRCLE COMMAND WHITELIST (for interactive mode security)
# ============================================================================

# Allowed sixcircle commands for interactive mode
# Format: 'command_name': ('method_name', needs_args)
SIXCIRCLE_ALLOWED_COMMANDS = {
    'wh': ('wh', False),      # where: show current position
    'pa': ('pa', False),      # show all parameters
    'ca': ('ca', False),      # show calculated positions
    'or0': ('or0', False),    # show first orientation reflection
    'or1': ('or1', False),    # show second orientation reflection
    'setlat': ('setlat', True),  # set lattice (requires args)
    'setaz': ('setaz', True),    # set azimuthal reference
    'freeze': ('freeze', True),  # freeze angles
    'unfreeze': ('unfreeze', True),  # unfreeze angles
}

# Help text for sixcircle commands
SIXCIRCLE_COMMAND_HELP = """
Available sixcircle commands:
  wh         - Show current diffractometer position
  pa         - Show all parameters
  ca         - Show calculated angle positions
  or0        - Show first orientation reflection
  or1        - Show second orientation reflection
  
Type 'sixc <command>' to execute
"""


# ============================================================================
# SAMPLE PARAMETERS - AuTe2
# ============================================================================

# Lattice parameters for AuTe2 (to be refined from experiment)
# These are initial values - will be updated after orientation matrix refinement
# Lattice parameters for AuTe2
# From Reithmayer et al., Acta Crystallographica B49, 6 (1993)
LATTICE_PARAMS = {
    'a': 7.189,      # Angstroms
    'b': 4.407,
    'c': 5.069,
    'alpha': 90.0,   # degrees
    'beta': 89.96,
    'gamma': 90.0
}

# Modulation vector (in r.l.u. of average structure)
# This is the CDW/superstructure modulation
MODULATION_VECTOR = (-0.4076, 0.0, 0.4479)  # From Schutte et al., Acta Cryst. B44, 486 (1988) at 298K  # (qh, qk, ql)

# Sample name/description
SAMPLE_NAME = "AuTe2_CDW"

# ============================================================================
# IXS CROSS-SECTION PARAMETERS
# ============================================================================

# Sample temperature (Kelvin), used for the Bose population factor in the
# IXS cross-section (see single_q_analysis.py: kT_THz derived from this).
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

def set_frozen_angles(code: str):
    """Set which three angles are frozen (3-digit sixcircle code, e.g. '456')."""
    global FROZEN_ANGLES
    code = str(code).strip()
    if not (code.isdigit() and len(code) == 3 and len(set(code)) == 3):
        raise ValueError(
            f"Frozen code must be three distinct digits 0-9, got {code!r}")
    FROZEN_ANGLES = code

def set_frozen_values(**kwargs):
    """Set the held values (degrees) of frozen angles, e.g. mu=-0.17, gam=0."""
    for key, value in kwargs.items():
        if key not in FROZEN_VALUES:
            raise ValueError(f"Unknown angle: {key}")
        FROZEN_VALUES[key] = float(value)

def set_angle_limit(angle: str, lower: float, upper: float):
    """Set the (lower, upper) limit in degrees for one angle."""
    if angle not in ANGLE_LIMITS:
        raise ValueError(f"Unknown angle: {angle}")
    if lower >= upper:
        raise ValueError(f"lower ({lower}) must be < upper ({upper})")
    ANGLE_LIMITS[angle] = (float(lower), float(upper))

def reset_angle_limits():
    """Reset all angle limits to the full ±180° range."""
    for angle in ANGLE_LIMITS:
        ANGLE_LIMITS[angle] = (-180.0, 180.0)


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
    print(f"  Wavelength: {WAVELENGTH:.5f} A  ({ENERGY_KEV:.3f} keV)")
    print(f"  Incident type: {INCIDENT_TYPE}")
    _ang = {'0':'tth','1':'th','2':'chi','3':'phi','4':'mu',
            '5':'gam','6':'omega','7':'azimuth','8':'alpha','9':'beta'}
    frozen_desc = ", ".join(
        f"{_ang[d]}={FROZEN_VALUES.get(_ang[d], 0.0):.4g}" for d in FROZEN_ANGLES)
    print(f"  Frozen angles: {FROZEN_ANGLES} ({frozen_desc})")
    print(f"  Analyzer gaps: V={DEFAULT_AGAP_V}μm, H={DEFAULT_AGAP_H}μm")
    non_default = {a: lu for a, lu in ANGLE_LIMITS.items() if lu != (-180.0, 180.0)}
    if non_default:
        lims = ", ".join(f"{a}∈[{lo:g},{hi:g}]" for a, (lo, hi) in non_default.items())
        print(f"  Angle limits: {lims}")
    else:
        print(f"  Angle limits: all ±180°")
    
    print("\n[Sample - {0}]".format(SAMPLE_NAME))
    print(f"  a = {LATTICE_PARAMS['a']:.4f} Å")
    print(f"  b = {LATTICE_PARAMS['b']:.4f} Å")
    print(f"  c = {LATTICE_PARAMS['c']:.4f} Å")
    print(f"  α = {LATTICE_PARAMS['alpha']:.2f}°")
    print(f"  β = {LATTICE_PARAMS['beta']:.2f}°")
    print(f"  γ = {LATTICE_PARAMS['gamma']:.2f}°")
    
    print("\n[Modulation]")
    print(f"  q_mod = ({MODULATION_VECTOR[0]:.3f}, {MODULATION_VECTOR[1]:.3f}, {MODULATION_VECTOR[2]:.3f})")

    print("\n[Thermal / IXS]")
    print(f"  Temperature: {TEMPERATURE} K")
    print(f"  Debye-Waller: {USE_DEBYE_WALLER} (not yet implemented)")
    print("="*70 + "\n")

