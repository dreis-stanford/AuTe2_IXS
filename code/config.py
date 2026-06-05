"""
Configuration file for AuTe2 IXS project
"""

import os

# Path to sixcircle code (Alfred's six-circle diffractometer code)
# Update this to point to where sixcircle is installed on your system
SIXCIRCLE_PATH = os.path.expanduser("~/Documents/myPython/Others/sixcircle_1p85")

# Alternative: set to None if sixcircle is not available
# SIXCIRCLE_PATH = None

# Check if sixcircle exists
if SIXCIRCLE_PATH and not os.path.exists(SIXCIRCLE_PATH):
    print(f"Warning: sixcircle not found at {SIXCIRCLE_PATH}")
    print("Set SIXCIRCLE_PATH=None in code/config.py if not needed")
