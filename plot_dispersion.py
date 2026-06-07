#!/usr/bin/env python3
"""
Launcher for AuTe2 phonon dispersion plotting
"""

from code.plot_dispersion import plot_aute2_dispersion

if __name__ == '__main__':
    import sys
    
    # Parse optional arguments
    fc_file = 'data/AuTe_2_m.fc'
    if len(sys.argv) > 1:
        fc_file = sys.argv[1]
    
    plot_aute2_dispersion(fc_file=fc_file, show=True, block=True)
