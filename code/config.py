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

# IMPORTANT: We maintain TWO structures for different purposes:
# 1. DFT structure: For phonon eigenvectors and one-phonon structure factors
# 2. Experimental structure: For Bragg/satellite positions and diffractometer geometry

# Which structure to use for what:
PHONON_STRUCTURE = 'DFT'          # Use DFT positions for phonon calculations
DIFFRACTION_STRUCTURE = 'EXPERIMENTAL'  # Use experimental for Bragg geometry


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

# Experimental atomic positions (fractional coordinates)
# C2/m space group, average structure
# From Reithmayer et al., Acta Crystallographica B49, 6 (1993)
EXPERIMENTAL_STRUCTURE = {
    'source': 'Reithmayer et al., Acta Cryst. B49, 6 (1993)',
    'lattice': {
        'a': 7.189,     # Angstroms
        'b': 4.407,
        'c': 5.069,
        'alpha': 90.0,  # degrees
        'beta': 89.96,
        'gamma': 90.0
    },
    'positions': {
        'Au': {
            'wyckoff': '2a',
            'frac_coords': [
                [0.0, 0.0, 0.0],
                [0.5, 0.5, 0.0],  # C-centering
            ],
            'mass': 196.97
        },
        'Te': {
            'wyckoff': '4i',
            'frac_coords': [
                [ 0.6884, 0.0,  0.2878],
                [-0.6884, 0.0, -0.2878],
                [ 0.1884, 0.5,  0.2878],
                [ 0.8116, 0.5,  0.7122],
            ],
            'mass': 127.60
        }
    }
}

# Modulation vector (in r.l.u. of average structure)
# This is the CDW/superstructure modulation
MODULATION_VECTOR = (-0.4076, 0.0, 0.4479)  # From Schutte et al., Acta Cryst. B44, 486 (1988) at 298K  # (qh, qk, ql)

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
    
    print("\n[Default Crystal]")
    print(f"  Using: {DEFAULT_CRYSTAL} structure")
    
    print("\n[Sample Orientation]")
    print(f"  Surface normal: {SURFACE_NORMAL} (conv. r.l.u.)")
    print(f"  Azimuthal ref:  {AZIMUTHAL_REFERENCE} (conv. r.l.u.)")
    
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


# ============================================================================
# AUTO-LOAD DFT STRUCTURE FROM FORCE CONSTANTS
# ============================================================================

def load_dft_structure_from_fc(fc_file='data/AuTe_2_m.fc'):
    """
    Auto-load DFT structure from force constants file
    
    Returns:
    --------
    dict with 'lattice' and 'positions' matching EXPERIMENTAL_STRUCTURE format
    """
    try:
        from .force_constants import ForceConstants
        import numpy as np
        
        print(f"Auto-loading DFT structure from {fc_file}...")
        fc = ForceConstants(fc_file)
        
        # Extract lattice parameters (already calculated in ForceConstants)
        structure = {
            'source': f'DFT (auto-loaded from {fc_file})',
            'lattice': {
                'a': fc.a,
                'b': fc.b,
                'c': fc.c,
                'alpha': fc.alpha,
                'beta': fc.beta,
                'gamma': fc.gamma
            },
            'positions': {}
        }
        
        # Extract atomic positions
        # fc.xs contains fractional coordinates [nat x 3]
        # fc.atom_type_map tells us which atom is which type
        # fc.symbols gives the element names
        
        for symbol in set(fc.symbols):
            structure['positions'][symbol] = {
                'frac_coords': [],
                'mass': None
            }
        
        for i in range(fc.nat):
            atom_type_idx = fc.atom_type_map[i] - 1  # Convert to 0-indexed (1,2 -> 0,1)
            symbol = fc.symbols[atom_type_idx]
            pos = fc.xs[i].tolist()
            structure['positions'][symbol]['frac_coords'].append(pos)
            structure['positions'][symbol]['mass'] = fc.masses[atom_type_idx]
        
        print(f"  ✓ Loaded DFT lattice: a={fc.a:.4f}, b={fc.b:.4f}, c={fc.c:.4f} Å")
        print(f"  ✓ Lattice angles: α={fc.alpha:.2f}°, β={fc.beta:.2f}°, γ={fc.gamma:.2f}°")
        for symbol, data in structure['positions'].items():
            n_atoms = len(data['frac_coords'])
            print(f"  ✓ Found {n_atoms} {symbol} atom(s), mass={data['mass']:.2f} amu")
        
        return structure
        
    except Exception as e:
        import traceback
        print(f"  ✗ Could not auto-load DFT structure: {e}")
        print("    Using placeholder DFT structure")
        traceback.print_exc()
        return None

# Try to auto-load DFT structure
print("\n" + "="*70)
print("Loading DFT Structure")
print("="*70)
_dft_loaded = load_dft_structure_from_fc()

if _dft_loaded is not None:
    DFT_STRUCTURE = _dft_loaded
else:
    # Placeholder if auto-load fails
    DFT_STRUCTURE = {
        'source': 'Placeholder - DFT structure not loaded',
        'lattice': {
            'a': 8.15, 'b': 4.47, 'c': 5.08,
            'alpha': 90.0, 'beta': 90.0, 'gamma': 90.0
        },
        'positions': {
            'Au': {'frac_coords': [[0.0, 0.0, 0.0], [0.5, 0.5, 0.0]], 'mass': 196.97},
            'Te': {'frac_coords': [[0.69, 0.0, 0.29], [0.31, 0.0, 0.71]], 'mass': 127.60}
        }
    }

print("="*70 + "\n")

# ============================================================================
# HELPER FUNCTIONS FOR STRUCTURE ACCESS
# ============================================================================

def get_phonon_structure():
    """Get structure to use for phonon/eigenvector calculations"""
    if PHONON_STRUCTURE == 'DFT':
        return DFT_STRUCTURE
    elif PHONON_STRUCTURE == 'EXPERIMENTAL':
        print("WARNING: Using experimental structure for phonons - eigenvectors may be inconsistent!")
        return EXPERIMENTAL_STRUCTURE
    else:
        raise ValueError(f"Unknown phonon structure: {PHONON_STRUCTURE}")

def get_diffraction_structure():
    """Get structure to use for Bragg peak/satellite calculations"""
    if DIFFRACTION_STRUCTURE == 'DFT':
        return DFT_STRUCTURE
    elif DIFFRACTION_STRUCTURE == 'EXPERIMENTAL':
        return EXPERIMENTAL_STRUCTURE
    else:
        raise ValueError(f"Unknown diffraction structure: {DIFFRACTION_STRUCTURE}")

def get_structure_info():
    """Print which structures are being used"""
    print("\nStructure Usage:")
    print(f"  Phonons/Eigenvectors: {PHONON_STRUCTURE}")
    print(f"    Source: {get_phonon_structure()['source']}")
    print(f"  Bragg/Diffractometer: {DIFFRACTION_STRUCTURE}")
    print(f"    Source: {get_diffraction_structure()['source']}")

# ============================================================================
# DEFAULT CRYSTAL SETTING
# ============================================================================

# Which crystal structure to use for calculations
# This affects which lattice is used for sixcircle and reciprocal space calculations
DEFAULT_CRYSTAL = 'EXPERIMENTAL'  # Options: 'DFT', 'EXPERIMENTAL'

def get_default_lattice():
    """
    Get the default lattice parameters for diffractometer calculations
    
    Returns:
    --------
    dict with lattice parameters
    """
    if DEFAULT_CRYSTAL == 'DFT':
        return DFT_STRUCTURE['lattice']
    elif DEFAULT_CRYSTAL == 'EXPERIMENTAL':
        return EXPERIMENTAL_STRUCTURE['lattice']
    else:
        raise ValueError(f"Unknown default crystal: {DEFAULT_CRYSTAL}")

def set_default_crystal(crystal: str):
    """
    Set which crystal structure to use as default
    
    Parameters:
    -----------
    crystal : 'DFT' or 'EXPERIMENTAL'
    """
    global DEFAULT_CRYSTAL
    if crystal.upper() not in ['DFT', 'EXPERIMENTAL']:
        raise ValueError("crystal must be 'DFT' or 'EXPERIMENTAL'")
    DEFAULT_CRYSTAL = crystal.upper()
    print(f"✓ Default crystal set to {DEFAULT_CRYSTAL}")

# ============================================================================
# SAMPLE ORIENTATION / AZIMUTHAL REFERENCE
# ============================================================================

# Surface normal (sample normal for thin film or preferred scattering geometry)
# Defined in conventional reciprocal lattice units
SURFACE_NORMAL = (0, 0, 1)  # Along c-axis (conventional)

# Azimuthal reference direction (in-plane reference)
# Defined in conventional reciprocal lattice units
AZIMUTHAL_REFERENCE = (0, 1, 0)  # Along b-axis (conventional)

# Note: For AuTe2 conventional cell (monoclinic C2/m):
#   α = 90° (b⊥c), β = 89.96° (a⊥c nearly), γ = 90° (a⊥b)
#   So b and c are perpendicular (suitable for azimuthal reference)

def get_surface_normal():
    """Get surface normal in conventional r.l.u."""
    return SURFACE_NORMAL

def get_azimuthal_reference():
    """Get azimuthal reference in conventional r.l.u."""
    return AZIMUTHAL_REFERENCE

def set_sample_orientation(surface_normal=None, azimuthal_ref=None):
    """
    Set sample orientation for diffractometer
    
    Parameters:
    -----------
    surface_normal : tuple (h, k, l) in conventional r.l.u.
    azimuthal_ref : tuple (h, k, l) in conventional r.l.u.
    """
    global SURFACE_NORMAL, AZIMUTHAL_REFERENCE
    
    if surface_normal is not None:
        SURFACE_NORMAL = tuple(surface_normal)
        print(f"✓ Surface normal set to {SURFACE_NORMAL}")
    
    if azimuthal_ref is not None:
        AZIMUTHAL_REFERENCE = tuple(azimuthal_ref)
        print(f"✓ Azimuthal reference set to {AZIMUTHAL_REFERENCE}")
    
    # Check if perpendicular
    import numpy as np
    if surface_normal is not None and azimuthal_ref is not None:
        dot = np.dot(surface_normal, azimuthal_ref)
        if abs(dot) < 1e-6:
            print("  ✓ Vectors are perpendicular")
        else:
            print(f"  ⚠ Warning: Vectors not perpendicular (dot product = {dot})")

