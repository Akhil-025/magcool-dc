"""
Regression tests for core/regime_crossover_analysis.py.

Pins the module's own headline finding (no COP crossover exists, and
refrigerant elimination does not rescue the total-emissions comparison at
this repo's own modeled COP gap) so a future, unrelated change elsewhere
in the AMR/loss-model/emissions stack that flips either conclusion gets
caught here rather than silently going unnoticed -- same discipline as
tests/test_hysteresis_sensitivity.py locking down that module's own null
result.
"""
import pytest

from core.regime_crossover_analysis import (
    run_cop_crossover_search, run_emissions_crossover_check,
    COP_SEARCH_SPANS_K, COP_SEARCH_VCC_ETA_RANGE,
)


def test_cop_crossover_search_runs_and_returns_expected_shape():
    result = run_cop_crossover_search(verbose=False)
    assert "rows" in result and "any_crossover_found" in result
    assert len(result["rows"]) == len(COP_SEARCH_SPANS_K)
    for row in result["rows"]:
        for key in ("span_K", "T_cold_K", "best_AMR_COP_electrical",
                    "best_AMR_design", "VCC_COPs_by_eta",
                    "worst_case_VCC_COP", "AMR_beats_worst_case_VCC"):
            assert key in row
        assert set(row["VCC_COPs_by_eta"].keys()) == set(COP_SEARCH_VCC_ETA_RANGE)


def test_cop_crossover_search_finds_no_crossover_at_current_model_state():
    """Locks down this module's own headline finding (see its top-level
    docstring): across every span/design combination searched, this
    repo's own AMR model does not beat vapor-compression's OWN worst-case
    setting. If a future change to core/amr_cycle.py, core/loss_model.py,
    or core/cascade.py's staging logic ever flips this, that is
    significant enough to be a deliberate, reviewed change -- not
    something that should happen silently."""
    result = run_cop_crossover_search(verbose=False)
    assert result["any_crossover_found"] is False, (
        "A COP crossover was found where none existed before -- this "
        "contradicts every real-world source checked when this module was "
        "built (Magnotherm Eclipse, Polaris, Ames Lab all report 'matches', "
        "not 'beats', vapor-compression). Verify this is a genuine model "
        "improvement, not a bug, before treating it as good news.")
    # Every individual row should also independently agree.
    assert all(not row["AMR_beats_worst_case_VCC"] for row in result["rows"])


def test_best_amr_cop_decreases_with_span_within_documented_tolerance():
    """Sanity check independent of the crossover question: best achievable
    AMR_COP_electrical should decrease as span widens (more temperature
    lift per unit of magnetic work is intrinsically harder), matching the
    same qualitative trend results/comparison_table.csv already shows.
    A violation here would suggest a bug in the search itself (e.g. a
    design that should remain feasible at small span silently dropping
    out), not a physics finding worth reporting on its own."""
    result = run_cop_crossover_search(verbose=False)
    cops = [row["best_AMR_COP_electrical"] for row in result["rows"]]
    # Allow for the fact this is a best-of-grid search, not a smooth
    # analytic function, but the overall trend across 8 spans should still
    # be clearly decreasing (final value well under the first).
    assert cops[-1] < cops[0]
    assert cops[-1] < 0.5 * cops[0]


def test_emissions_crossover_check_runs_and_reproduces_documented_ratio():
    """Locks down the ~3.6x AMR-vs-VCC total-emissions ratio at the
    Eclipse-derived default operating point (amr_cop=1.76, vcc_cop=6.66,
    0.4kW) documented in this module's own top-level docstring."""
    result = run_emissions_crossover_check(verbose=False)
    assert "amr" in result and "vcc" in result
    assert result["amr_wins_total_emissions"] is False
    assert result["ratio_amr_to_vcc_total_emissions"] == pytest.approx(3.6, abs=1.0)


def test_emissions_crossover_check_at_a_closer_cop_gap_can_still_favor_vcc():
    """Even a much smaller (more favorable to AMR) COP gap than the
    default Eclipse-derived one should not automatically flip the total-
    emissions comparison, since core/emissions.py's own default
    refrigerant-leak assumptions are small relative to typical operational
    emissions at any realistic COP gap -- confirms the module's own
    'operational emissions dominate' framing holds beyond just the one
    default data point, not only at the specific default arguments."""
    result = run_emissions_crossover_check(capacity_kW=0.4, amr_cop=5.0,
                                            vcc_cop=6.0, verbose=False)
    assert result["amr_wins_total_emissions"] is False
