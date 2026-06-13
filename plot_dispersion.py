#!/usr/bin/env python3
"""
Launcher for AuTe2 phonon dispersion plotting.

Opens the dispersion plot in its own (detached) process so the window stays
open and interactive — zoom/pan — while this terminal returns immediately.
Close the plot window any time; it doesn't affect this shell. A PNG is also
saved (aute2_dispersion.png).

Usage:
    python plot_dispersion.py [force_constants.fc]
    python plot_dispersion.py --wait   # block here until the window is closed
"""

import subprocess
import sys

if __name__ == '__main__':
    args = [a for a in sys.argv[1:] if a != '--wait']
    wait = '--wait' in sys.argv

    cmd = [sys.executable, '-m', 'code.plot_dispersion'] + args

    if wait:
        # Run in the foreground (blocks until the window is closed).
        subprocess.run(cmd)
    else:
        # Detached: the plot window stays open in its own process and this
        # terminal is free again right away.
        subprocess.Popen(cmd, stdout=subprocess.DEVNULL,
                         stderr=subprocess.DEVNULL)
        print("Opening dispersion plot in a separate window "
              "(this terminal is free; close the window any time).")
