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
