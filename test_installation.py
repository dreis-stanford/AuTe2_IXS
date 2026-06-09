#!/usr/bin/env python
"""Test script to verify installation"""
import sys
import os

print("Testing AuTe2 IXS Installation")
print("=" * 60)

# Test Python version
print("\n1. Python version:", sys.version_info.major, ".", sys.version_info.minor)
if sys.version_info < (3, 8):
    print("   ERROR: Need Python 3.8+")
    sys.exit(1)

# Test packages
print("\n2. Checking packages...")
for pkg in ['numpy', 'scipy', 'matplotlib']:
    try:
        __import__(pkg)
        print(f"   ✓ {pkg}")
    except:
        print(f"   ✗ {pkg} MISSING - run: pip install {pkg}")
        sys.exit(1)

# Test project files
print("\n3. Checking files...")
for f in ['code/phonons.py', 'data/AuTe_2_m.fc']:
    if os.path.exists(f):
        print(f"   ✓ {f}")
    else:
        print(f"   ✗ {f} NOT FOUND")
        sys.exit(1)

# Test import (also exercises code/__init__.py, which auto-loads the
# DFT structure from data/AuTe_2_m.fc - expect a structure summary below)
print("\n4. Testing imports...")
from code.force_constants import ForceConstants
from code.phonons import calc_Dq
print("   ✓ All imports work")

print("\n✓ Installation successful!")
