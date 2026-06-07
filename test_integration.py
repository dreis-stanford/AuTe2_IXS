#!/usr/bin/env python3
"""
Run integration tests for sixcircle interface
"""

import sys
from pathlib import Path

# Add code directory to path
sys.path.insert(0, str(Path(__file__).parent))

from code.test_sixcircle_integration import run_all_tests

if __name__ == '__main__':
    run_all_tests()
