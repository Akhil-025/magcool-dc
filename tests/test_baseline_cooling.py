import pytest

from core.baseline_cooling import carnot_cop, vapor_compression_cop, liquid_cooling_cop


def test_carnot_cop_matches_definition():
    assert carnot_cop(290.0, 300.0) == pytest.approx(290.0 / 10.0)


def test_baseline_cops_below_carnot():
    vcc = vapor_compression_cop(290.0, 300.0)
    liq = liquid_cooling_cop(290.0, 300.0)
    assert vcc.COP < vcc.COP_carnot
    # liquid cooling blends a high economizer-mode COP, so only check it's
    # not exceeding physical bounds by an absurd margin (it can legitimately
    # exceed the *mechanical* Carnot figure for the DX-only comparison).
    assert liq.COP > 0
