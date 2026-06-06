# AuTe2 IXS Analysis Project

Python tools for calculating phonon dispersions and Inelastic X-ray Scattering (IXS) cross-sections.

## Materials

- AuTe2: Base-centered monoclinic (C2/m), 3 atoms/cell, 9 phonon branches
- Silicon: FCC diamond structure, 2 atoms/cell, 6 phonon branches

## Quick Start

### AuTe2 Interactive Analysis
python analyze_q.py

Commands: Enter Q as "0.5 0 0", type "prim"/"conv"/"cart" for coordinates, "meV"/"THz"/"invcm" for units

### Silicon Interactive Analysis  
python analyze_si.py

## Key Files

- code/single_q_analysis.py - Interactive IXS for AuTe2
- code/phonons.py - Phonon calculations
- code/ixs.py - IXS cross-sections
- analyze_q.py - AuTe2 launcher
- analyze_si.py - Silicon launcher

## Recent Work

Just added:
1. Signed longitudinal components (Q·e per atom) - shows direction along Q
2. Frequency unit switching (meV/THz/cm-1)

## To-Do

1. Force constant symmetry validation
2. Debye-Waller factors
3. 6-circle diffractometer angles
4. Investigate Si dispersion discrepancies

## For New Sessions

Tell me: "Continuing AuTe2 IXS project" and show "git log --oneline -5"
