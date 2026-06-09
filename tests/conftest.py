"""Pytest config: make project root importable; stub scipy if missing.

The stub only provides what these tests need (eigh). On a normal install
with scipy present, the real scipy is used.
"""
import os
import sys

PROJECT_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, PROJECT_ROOT)

try:
    import scipy.linalg  # noqa: F401
except ImportError:
    import types
    import numpy as np

    scipy = types.ModuleType('scipy')
    scipy_linalg = types.ModuleType('scipy.linalg')
    scipy_linalg.eigh = np.linalg.eigh
    scipy_optimize = types.ModuleType('scipy.optimize')
    scipy_optimize.linear_sum_assignment = None  # only used for multi-q sorting
    scipy.linalg, scipy.optimize = scipy_linalg, scipy_optimize
    sys.modules.update({'scipy': scipy, 'scipy.linalg': scipy_linalg,
                        'scipy.optimize': scipy_optimize})


def pytest_configure(config):
    # Force-constant files are loaded with cwd-relative paths (e.g. "data/...")
    os.chdir(PROJECT_ROOT)
