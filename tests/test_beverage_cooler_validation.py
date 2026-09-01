"""
Regression tests for core/beverage_cooler_validation.py.

Pins (a) the literature reference constants read from the two source
publications, so a future edit can't silently substitute a different
figure, and (b) the two checks' own observed model results within a
generous tolerance band, so an unrelated change elsewhere in the AMR/
loss-model stack that materially shifts these numbers gets noticed here
-- same discipline as tests/test_giguere_validation.py.
"""
import pytest

from core.beverage_cooler_validation import (
    run_eclipse_directional_check, run_polaris_second_law_validation,
    ECLIPSE_REPORTED_ENERGY_SAVING_PCT, ECLIPSE_T_COLD_C, ECLIPSE_T_AMBIENT_C,
    POLARIS_FIELD_T, POLARIS_SPAN_K, POLARIS_PLUGIN_COP,
    POLARIS_SECOND_LAW_EFF_PCT, POLARIS_SPECIFIC_COOLING_W_PER_KG,
)


def test_eclipse_reference_values_match_the_source():
    """Read directly from naturalrefrigerants.com's coverage of Magnotherm's
    ATMOsphere Europe Summit 2025 presentation, independently corroborated
    by EIT RawMaterials and refindustry.com/HAUSER -- pin them."""
    assert ECLIPSE_REPORTED_ENERGY_SAVING_PCT == 15.0
    assert ECLIPSE_T_COLD_C == pytest.approx(4.5)
    # T_ambient is this repo's own assumption (not independently reported by
    # any source), not a literature figure -- just confirm it stays in a
    # physically sane retail-floor range if someone edits it.
    assert 15.0 <= ECLIPSE_T_AMBIENT_C <= 30.0


def test_polaris_reference_values_match_the_paper():
    """Read directly from Liang, Pickett, Hermann et al., "Polaris: From
    Laboratory Prototypes to Market-Ready Sustainable Magnetic Beverage
    Coolers," Applied Thermal Engineering / ScienceDirect (2025) -- pin
    them so a future edit can't silently substitute a different figure."""
    assert POLARIS_FIELD_T == pytest.approx(0.8)
    assert POLARIS_SPAN_K == pytest.approx(15.0)
    assert POLARIS_PLUGIN_COP == pytest.approx(1.0)
    assert POLARIS_SECOND_LAW_EFF_PCT == pytest.approx(5.4)
    assert POLARIS_SPECIFIC_COOLING_W_PER_KG == pytest.approx(131.0)


def test_eclipse_check_runs_and_returns_expected_keys():
    """Does not assert model-vs-literature agreement (see this function's
    own honesty flag -- Eclipse's internal design is proprietary, so this
    repo's own default parameters are not expected to match it closely);
    only that the function runs cleanly end-to-end and returns a complete,
    well-formed result."""
    result = run_eclipse_directional_check()
    for key in ("T_cold_K", "T_hot_K", "span_K", "AMR_Qc_W", "AMR_COP_electrical",
                "AMR_n_stages", "VCC_COP", "model_predicted_saving_pct",
                "reported_saving_pct", "source"):
        assert key in result
    assert result["reported_saving_pct"] == 15.0
    assert result["T_cold_K"] == pytest.approx(ECLIPSE_T_COLD_C + 273.15)


def test_polaris_check_second_law_efficiency_is_within_documented_tolerance():
    """Locks down the ~+17.5% relative error between this repo's own model
    (COP_electrical / COP_carnot at Gd/0.8T/15K span, using
    staged_baseline_result()'s already-parasitic-loss-aware COP_electrical,
    NOT its exergy_eff field -- see run_polaris_second_law_validation()'s
    own honesty flag for why) and Polaris's own peer-reviewed, directly
    reported 5.4% second-law efficiency, so a future change to
    core/amr_cycle.py, core/loss_model.py, or core/cascade.py's staging
    logic that materially shifts this comparison gets caught here."""
    result = run_polaris_second_law_validation()
    assert result["AMR_Qc_W"] > 0, (
        "Model returned infeasible (Qc=0) at Polaris's own reported 0.8T/15K "
        "operating point -- see run_polaris_second_law_validation()'s own "
        "max_stages docstring note; this default should be feasible.")
    assert result["model_second_law_eff_pct"] is not None
    rel_err = (result["model_second_law_eff_pct"] - POLARIS_SECOND_LAW_EFF_PCT) \
        / POLARIS_SECOND_LAW_EFF_PCT
    assert rel_err == pytest.approx(0.175, abs=0.25), (
        f"Model second-law efficiency {result['model_second_law_eff_pct']:.1f}% vs. "
        f"Polaris's reported {POLARIS_SECOND_LAW_EFF_PCT}% drifted well outside the "
        f"documented ~+17.5% relative error band (got {rel_err:+.1%}) -- re-check "
        "whether this is a genuine model change or a regression.")


def test_polaris_check_returns_expected_keys():
    result = run_polaris_second_law_validation()
    for key in ("T_cold_K", "span_K", "field_T", "AMR_Qc_W", "AMR_n_stages",
                "AMR_COP_electrical", "COP_carnot", "model_second_law_eff_pct",
                "reported_second_law_eff_pct", "reported_plugin_COP", "source"):
        assert key in result
    assert result["field_T"] == pytest.approx(POLARIS_FIELD_T)
    assert result["span_K"] == pytest.approx(POLARIS_SPAN_K)


def test_polaris_second_law_eff_not_confused_with_magnetic_only_exergy_eff():
    """Regression guard for the specific mistake this module's own honesty
    flag warns against: comparing Polaris's plug-in figure against
    StagedBaselineResult.exergy_eff (magnetic-work-only, excludes parasitic
    losses) instead of COP_electrical/COP_carnot gives a ~53% figure --
    roughly 10x too optimistic and not a meaningful match to Polaris's
    5.4%. This test confirms the module's ACTUAL reported number stays far
    below that wrong-quantity magnitude, catching a future accidental
    swap back to the wrong field."""
    result = run_polaris_second_law_validation()
    assert result["model_second_law_eff_pct"] < 20.0, (
        "model_second_law_eff_pct is implausibly high -- check that "
        "run_polaris_second_law_validation() is still using COP_electrical, "
        "not COP (magnetic-only) or exergy_eff, in its second-law-efficiency "
        "calculation.")
