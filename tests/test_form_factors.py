"""Tests for code/form_factors.py.

Covers the sign-convention check for xraylib's dispersion-correction
imaginary part f''(E): see code/form_factors.py for the full derivation of
why calc_form_factor() uses Im(f) = +Fii(Z,E) (not -Fii) given this
codebase's exp(-2*pi*i*Q.r) phase convention (code/ixs.py).
"""
import numpy as np
import pytest
import scipy.constants as sc

from code import config
from code.constants import const
from code.form_factors import calc_form_factor, HAS_XRAYLIB

pytestmark = pytest.mark.skipif(not HAS_XRAYLIB, reason="xraylib not available")

import xraylib


def _f2_optical_theorem(symbol, A, E_keV):
    """Physically-positive f''(E), via the optical theorem applied to
    xraylib's own photoelectric cross section (CS_Photo) -- an
    independent check of xraylib's Fi/Fii tables."""
    Z = xraylib.SymbolToAtomicNumber(symbol)
    cs_photo = xraylib.CS_Photo(Z, E_keV)        # cm^2/g
    sigma_pe_A2 = cs_photo * A / sc.N_A * 1e16   # Angstrom^2 / atom
    hc_keV_A = const.hc / 1000.0                 # keV.Angstrom
    return sigma_pe_A2 * E_keV / (2 * const.re * hc_keV_A)


@pytest.mark.parametrize("symbol,A", [('Au', 196.96657), ('Te', 127.6)])
def test_xraylib_fii_sign_is_flipped_from_optical_theorem(symbol, A):
    """-Fii(Z,E) (not +Fii) is the physically-positive (absorptive) f''(E):
    it must match the optical-theorem f'' derived independently from
    xraylib's own photoelectric cross section (CS_Photo)."""
    Z = xraylib.SymbolToAtomicNumber(symbol)
    f2_optical = _f2_optical_theorem(symbol, A, config.ENERGY_KEV)
    f2_xraylib = xraylib.Fii(Z, config.ENERGY_KEV)

    assert f2_optical > 0
    assert f2_xraylib < 0
    assert -f2_xraylib == pytest.approx(f2_optical, rel=0.02)


def test_calc_form_factor_uses_xraylib_fii_sign_directly():
    """calc_form_factor's Im(f) = Fii(Z,E) as-is (the sign-flipped value
    relative to the physical f''). This is correct -- not a bug -- because
    it exactly cancels against this codebase's exp(-2*pi*i*Q.r) phase
    convention in code/ixs.py, giving |F(Q)|^2 equal to the standard
    (f0+f'+if'', f''>0, exp(+2*pi*i*Q.r)) convention. See
    code/form_factors.py for the derivation."""
    f = calc_form_factor(0.1927, 'Au', config.ENERGY_KEV, scale=1.0)
    Z = xraylib.SymbolToAtomicNumber('Au')
    assert np.imag(f) == pytest.approx(xraylib.Fii(Z, config.ENERGY_KEV))
