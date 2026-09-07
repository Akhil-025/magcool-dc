import math

from core.electrocaloric_cycle import (
    run_cycle,
    run_span_sweep,
    sensitivity_to_press_vs_paper_figure,
    MEASURED_DEVICE,
)


def test_measured_span_reproduces_paper_cop_exactly():
    """At the single calibration span (20.9K), COP_electrical should equal
    fraction_of_carnot * Carnot(span) by construction."""
    r = run_cycle(MEASURED_DEVICE["span_K"])
    expected = MEASURED_DEVICE["fraction_of_carnot"] * (r.T_cold / r.span_K)
    assert math.isclose(r.COP_electrical, expected, rel_tol=1e-9)
    assert r.is_measured_span
    assert not r.is_extrapolated
    assert r.measured_cooling_power_W == MEASURED_DEVICE["cooling_power_W"]


def test_other_spans_are_flagged_as_extrapolated():
    r = run_cycle(15.0)
    assert r.is_extrapolated
    assert not r.is_measured_span
    assert not r.is_far_extrapolation  # 15K is within [13.0, 20.9] bracketed range
    assert r.measured_cooling_power_W is None


def test_spans_below_13K_flagged_as_far_extrapolation():
    """Below Torello et al. 2020's 13.0K demonstrated span, there is no
    electrocaloric device evidence at all -- see module docstring."""
    r = run_cycle(10.0)
    assert r.is_extrapolated
    assert r.is_far_extrapolation


def test_cop_decreases_with_wider_span_at_fixed_efficiency_fraction():
    """Carnot COP = T_cold/span falls as span widens, so COP_electrical
    should too, holding fraction_of_carnot fixed."""
    r_narrow = run_cycle(5.0)
    r_wide = run_cycle(20.0)
    assert r_narrow.COP_electrical > r_wide.COP_electrical


def test_press_release_figure_gives_higher_cop_than_paper_figure():
    """64% (press) > 54% (paper's own abstract) -- this module defaults to
    the paper's figure; this test just locks in the direction of the
    flagged discrepancy."""
    s = sensitivity_to_press_vs_paper_figure()
    assert s["press_release_64pct"] > s["paper_54pct"]


def test_span_sweep_beats_matched_vcc_cop_at_the_measured_span():
    """Headline honest finding: at the ONE measured span (20.9K, nearest
    grid point 20K), this device's COP should exceed this repo's own VCC
    COP of ~6.11 at 20K span (comparison_table.csv)."""
    rows = run_span_sweep(spans_K=(20,), verbose=False)
    assert rows[0].COP_electrical > 6.11
