#!/usr/bin/env python
"""Quick launcher for interactive Silicon analysis"""
import sys
import os

# Add code directory to path
sys.path.insert(0, os.path.join(os.path.dirname(__file__), 'code'))

from single_q_analysis_si import interactive_mode

if __name__ == "__main__":
    interactive_mode()
