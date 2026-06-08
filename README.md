# AuTe2 IXS Experiment Planning

Python tools for planning Inelastic X-ray Scattering (IXS) experiments on AuTe₂ using DFT phonon predictions and six-circle diffractometer control.

## Overview

This project integrates:
- DFT phonon calculations to predict phonon dispersions and eigenvectors
- Sixcircle diffractometer interface for experiment planning at SPring-8
- Modulated structure handling for CDW satellite reflections
- Q-point optimization to maximize information from measurements

## Quick Start

Interactive Analysis:
    python analyze_q.py
    python analyze_si.py
    python test_integration.py

## Installation

Requirements: numpy scipy matplotlib phonopy

For SPring-8 users, edit code/config.py to set SIXCIRCLE_PATH

## Project Structure

- code/config.py - Configuration
- code/sixcircle_interface.py - Diffractometer control  
- code/modulated_structure.py - Satellite reflections
- code/q_optimizer.py - Q-point optimization
- analyze_q.py - AuTe2 launcher
- test_integration.py - Integration tests

## Recent Updates

Sixcircle Integration - 1,184+ lines of new code
- Diffractometer interface for SPring-8 BL43LXU
- Modulated structure and satellite handling
- Q-point optimization
- Simulation mode for testing

## Materials

AuTe2: Monoclinic C2/m, 3 atoms/cell, 9 phonon branches
Silicon: FCC diamond, 2 atoms/cell, 6 phonon branches

## Satellite Reflection Analysis

AuTe₂ has a CDW modulation with q_mod = (-0.4076, 0, 0.4479) at 298K.

Access satellites using 4-parameter input:
H K L m
where m is the satellite order.
Examples:
- `1 0 0` → Main Bragg peak
- `1 0 0 1` → First satellite: (1,0,0) + q_mod
- `1 0 0 -1` → First satellite: (1,0,0) - q_mod  
- `2 0 0 2` → Second order: (2,0,0) + 2×q_mod
Reference: Schutte et al., Acta Cryst. B44, 486 (1988)
## Interactive Commands
- `H K L [m]` - Analyze Q-point (optional m for satellites)
- `angles` - Calculate diffractometer angles
- `conv`, `prim`, `cart` - Switch coordinate systems
- `meV`, `THz`, `invcm` - Switch frequency units
- `move` - Move diffractometer (with confirmation)
- `sixc <cmd>` - Pass command to sixcircle
