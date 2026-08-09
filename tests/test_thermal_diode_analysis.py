import pytest

from core.thermal_diode_analysis import (check_frequency_ceiling_claim,
                                           sweep_frequency_with_and_without_diode,
                                           demo_cycle_time_reduction)


def test_frequency_ceiling_finding_reports_no_internal_cap():
    finding = check_frequency_ceiling_claim(verbose=False)
    assert finding["amr_system_has_internal_frequency_cap"] is False
    assert finding["upper_bound_documented_as_mechanical_switching_limit"] is False
    assert finding["optimize_py_frequency_upper_bound_Hz"] > 0


def test_sweep_diode_never_beats_no_diode_baseline():
    """Cost-only accounting (see module honesty flag): diode-assisted
    COP_electrical must never exceed the no-diode baseline."""
    rows = sweep_frequency_with_and_without_diode(frequencies=(1.0, 2.0), verbose=False)
    for f, cop_base, cop_diode, delta_pct in rows:
        assert cop_diode <= cop_base + 1e-12
        assert delta_pct <= 1e-9


def test_demo_cycle_time_reduction_returns_fraction_between_zero_and_one():
    reduction = demo_cycle_time_reduction(verbose=False)
    assert 0.0 <= reduction <= 1.0