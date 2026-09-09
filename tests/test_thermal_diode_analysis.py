import pytest

from core.thermal_diode_analysis import (check_frequency_ceiling_claim,
                                           sweep_frequency_with_and_without_diode,
                                           demo_cycle_time_reduction,
                                           check_against_andrade_2024_benchmark,
                                           ANDRADE_2024_BENCHMARK)
from core.thermal_diode import DEFAULT_FERROFLUID_THERMAL_SWITCH


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


# ---- Andrade et al. (2024) benchmark check (VALIDATION UPDATE) ----

def test_ferrofluid_switch_sweep_never_beats_no_diode_baseline():
    """Same cost-only-accounting structural property as
    MechanicalContactDiode's sweep (no offsetting heat-transfer benefit
    is modeled), now checked for FerrofluidThermalSwitch too."""
    rows = sweep_frequency_with_and_without_diode(
        frequencies=(1.0, 2.0), diode=DEFAULT_FERROFLUID_THERMAL_SWITCH,
        verbose=False)
    for f, cop_base, cop_diode, delta_pct in rows:
        assert cop_diode <= cop_base + 1e-12
        assert delta_pct <= 1e-9


def test_andrade_2024_benchmark_digitized_correctly():
    assert ANDRADE_2024_BENCHMARK["mcm_mass_g"] == pytest.approx(7.0)
    assert ANDRADE_2024_BENCHMARK["symmetric_cycling_advantage_over_bare_Gd"] is False
    assert ANDRADE_2024_BENCHMARK["asymmetric_result_is_hardware_measured"] is False


def test_check_against_andrade_2024_benchmark_finds_qualitative_agreement():
    """This repo's model, having no offsetting heat-transfer benefit for
    the diode mechanism, structurally shows no symmetric-cycling
    advantage -- which is exactly Andrade et al. (2024)'s own reported
    symmetric-cycling finding. This is agreement on a NEGATIVE result,
    not a validated positive quantitative prediction (see function
    docstring)."""
    finding = check_against_andrade_2024_benchmark(verbose=False)
    assert finding["model_shows_no_symmetric_advantage"] is True
    assert finding["qualitative_agreement_with_benchmark_symmetric_finding"] is True
    assert finding["asymmetric_60pct_span_claim_reproduced_by_this_repo"] is False
