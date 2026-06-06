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

# Test import
print("\n4. Testing imports...")
sys.path.insert(0, 'code')
from force_constants import ForceConstants
from phonons import calc_Dq
print("   ✓ All imports work")

print("\n✓ Installation successful!")
