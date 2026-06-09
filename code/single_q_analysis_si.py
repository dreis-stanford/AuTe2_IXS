"""
Interactive IXS analysis for single Q-points in Silicon.

Thin wrapper around the material-generic analyzer in single_q_analysis.py,
configured with FCC conv<->prim transforms and no CDW satellites. The
interactive ./analyze_si.py therefore behaves identically to ./analyze_q.py
(same commands, table format, and mode classification), except that
satellite input (H K L m) and the sixcircle diffractometer commands are
disabled for Silicon.
"""

from .single_q_analysis import (
    SingleQAnalyzer as _GenericAnalyzer,
    interactive_mode as _generic_interactive_mode,
    classify_polarization,  # re-export for convenience
)
from .fcc_structure import fcc_conv2prim_k, fcc_prim2conv_k


class SingleQAnalyzer(_GenericAnalyzer):
    """Single-Q analyzer configured for Silicon (FCC diamond)."""

    def __init__(self, xtal, Phi, masses, kT_THz):
        super().__init__(
            xtal, Phi, masses, kT_THz,
            conv2prim=fcc_conv2prim_k,
            prim2conv=fcc_prim2conv_k,
            enable_satellites=False,
            material='Silicon',
        )


def interactive_mode():
    """Interactive Q-point analysis for Silicon."""
    _generic_interactive_mode(material='Si')


if __name__ == "__main__":
    import sys
    import os
    sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
    interactive_mode()
