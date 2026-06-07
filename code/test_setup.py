"""
Test that everything is set up correctly
"""

import sys
import os

# Add code directory to path
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

import config
from .aute2_structure import AuTe2
import numpy as np
import xraylib

def test_imports():
    """Test that all imports work"""
    print("Testing imports...")
    
    # Test that config found sixcircle path
    if os.path.exists(config.SIXCIRCLE_PATH):
        print(f"✓ sixcircle found at {config.SIXCIRCLE_PATH}")
    else:
        print(f"✗ sixcircle NOT found at {config.SIXCIRCLE_PATH}")
        return False
    
    # Test scientific packages
    print("✓ numpy imported")
    
    # Test xraylib
    cs = xraylib.CS_Total(26, 10.0)  # Iron at 10 keV
    print(f"✓ xraylib working (Fe cross-section: {cs:.3e})")
    
    # Test h5py
    try:
        import h5py
        print("✓ h5py installed")
    except:
        print("✗ h5py missing")
        return False
    
    return True

def test_aute2():
    """Test AuTe2 class"""
    print("\nTesting AuTe2 class...")
    
    crystal = AuTe2()
    crystal.info()
    
    print("✓ AuTe2 class working")
    
    return True

def test_sixcircle_import():
    """Test that sixcircle can be imported when needed"""
    print("\nTesting sixcircle import...")
    
    try:
        # Import just the basic module, not rqd which has the bug
        import sys
        sys.path.insert(0, config.SIXCIRCLE_PATH)
        import sixcircle
        print("✓ sixcircle module can be imported")
        return True
    except Exception as e:
        print(f"✗ sixcircle import failed: {e}")
        return False

if __name__ == "__main__":
    print("=" * 50)
    print("AuTe2 IXS Project Setup Test")
    print("=" * 50)
    
    success = True
    success = test_imports() and success
    success = test_aute2() and success
    success = test_sixcircle_import() and success
    
    if success:
        print("\n" + "=" * 50)
        print("SUCCESS! Basic setup is working")
        print("=" * 50)
        print("\nNote: sixcircle_rqd has initialization issues")
        print("We'll work around this in your analysis code")
    else:
        print("\n" + "=" * 50)
        print("FAILED - check error messages above")
        print("=" * 50)
