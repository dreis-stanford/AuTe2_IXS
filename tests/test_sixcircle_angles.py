"""Tests for the sixcircle angle-calculation interface (added 2026-06-12).

Covers the rewrite that made the `angles` command work end-to-end:
 - calculate_angles() returns all six angles via ca_s() (no stdout parsing,
   no dependence on a nonexistent sixcircle.A array);
 - angle limits select among solution branches (e.g. force tth > 0);
 - frozen-angle mode/values are honoured.

These require the external sixcircle library; the whole module is skipped if
it isn't available (e.g. on CI). Run from project root: pytest tests/ -v
"""
import copy
import numpy as np
import pytest

from code import config

pytestmark = pytest.mark.skipif(
    not config.SIXCIRCLE_AVAILABLE,
    reason="external sixcircle library not available")


@pytest.fixture
def sixc():
    """A calculation-only interface, with config constraints restored after."""
    from code.sixcircle_interface import SixCircleInterface
    saved_limits = copy.deepcopy(config.ANGLE_LIMITS)
    saved_frozen = config.FROZEN_ANGLES
    saved_values = copy.deepcopy(config.FROZEN_VALUES)
    s = SixCircleInterface(simulation_mode=True)
    if s.sixc is None:
        pytest.skip("sixcircle present but failed to load")
    try:
        yield s
    finally:
        config.ANGLE_LIMITS.update(saved_limits)
        config.FROZEN_ANGLES = saved_frozen
        config.FROZEN_VALUES.update(saved_values)
        s.apply_limits()
        s.apply_frozen()


def test_calculate_angles_returns_nine(sixc):
    a = sixc.calculate_angles((0, 0, 2))
    assert set(a) == {'tth', 'th', 'chi', 'phi', 'mu', 'gam', 'omega', 'alpha', 'beta'}
    assert all(isinstance(v, float) for v in a.values())


def test_th_is_half_tth_for_symmetric_reflection(sixc):
    # In the standard frozen mode, th tracks tth/2 (omega frozen at 0).
    a = sixc.calculate_angles((0, 0, 2))
    assert a['th'] == pytest.approx(a['tth'] / 2, abs=1e-3)


def test_default_selects_top_face_branch(sixc):
    """The default tth limit (>=0) selects the top-face (tth>0) solution."""
    a = sixc.calculate_angles((0, 0, 2))
    assert a['tth'] > 0


def test_negative_tth_branch_via_limits(sixc):
    """Widening the tth limit exposes the bottom-face (tth<0) mirror solution."""
    pos = sixc.calculate_angles((0, 0, 2))
    sixc.set_limit('tth', -60.0, 0.0)
    neg = sixc.calculate_angles((0, 0, 2))
    assert neg['tth'] < 0
    # Same |Q|, so |tth| is unchanged; only the branch differs.
    assert abs(neg['tth']) == pytest.approx(pos['tth'], abs=1e-3)


def test_frozen_value_is_honoured(sixc):
    sixc.set_frozen_values(mu=-0.5)
    a = sixc.calculate_angles((0, 0, 2))
    assert a['mu'] == pytest.approx(-0.5, abs=1e-6)


def test_set_frozen_updates_mode(sixc):
    sixc.set_frozen('045')
    assert config.FROZEN_ANGLES == '045'
    assert sixc.sixc.g_frozen == '045'


def test_inaccessible_reflection_returns_none(sixc):
    # |Q| far beyond what the wavelength can reach -> no solution.
    assert sixc.calculate_angles((50, 50, 50)) is None


def test_invalid_frozen_code_rejected(sixc):
    with pytest.raises(ValueError):
        sixc.set_frozen('44')   # not three distinct digits


def test_bad_limit_order_rejected(sixc):
    with pytest.raises(ValueError):
        sixc.set_limit('tth', 60.0, 0.0)   # lower >= upper
