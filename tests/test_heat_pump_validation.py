"""Regression tests for core/heat_pump_validation.py."""
import pytest

from core.heat_pump_validation import (
    run_ames_lab_architecture_check,
    AMES_BASELINE_SPD_W_PER_KG, AMES_OPTIMIZED_SPD_W_PER_KG,
    AMES_PROJECTED_MAX_SPD_W_PER_KG, AMES_PERFORMANCE_RANGE_W,
)


def test_ames_lab_reference_values_match_the_paper():
    """Slaughter, Griffith, Czernuszewicz & Pecharsky, Applied Energy 377
    (2025) 124696 -- pin the reference figures so a future edit can't
    silently substitute a different one."""
    assert AMES_BASELINE_SPD_W_PER_KG == pytest.approx(5.9)
    assert AMES_OPTIMIZED_SPD_W_PER_KG == pytest.approx(81.3)
    assert AMES_PROJECTED_MAX_SPD_W_PER_KG == pytest.approx(114.0)
    assert AMES_PERFORMANCE_RANGE_W == (37.0, 43500.0)
    # Optimized SPD must be a genuine improvement over baseline, and the
    # projected ceiling must be >= the already-demonstrated optimized value
    # -- both are physically-required orderings, not just literature trivia.
    assert AMES_OPTIMIZED_SPD_W_PER_KG > AMES_BASELINE_SPD_W_PER_KG
    assert AMES_PROJECTED_MAX_SPD_W_PER_KG >= AMES_OPTIMIZED_SPD_W_PER_KG


def test_ames_lab_check_runs_and_returns_expected_keys():
    result = run_ames_lab_architecture_check()
    for key in ("T_cold_K", "span_K", "AMR_Qc_W", "AMR_n_stages",
                "AMR_COP_electrical", "model_specific_cooling_power_w_per_kg_MCM",
                "ames_baseline_SPD_w_per_kg_whole_device",
                "ames_optimized_SPD_w_per_kg_whole_device",
                "ames_projected_max_SPD_w_per_kg_whole_device",
                "ames_performance_range_W", "source"):
        assert key in result
    assert result["AMR_Qc_W"] > 0, (
        "Model returned infeasible (Qc=0) at the default residential-heat-"
        "pump-scale operating point -- this default should be feasible.")
    assert result["model_specific_cooling_power_w_per_kg_MCM"] is not None


def test_model_specific_cooling_power_is_not_silently_compared_to_ames_SPD():
    """Regression guard for the specific mistake this module's own honesty
    flag warns against: this repo's MCM-only specific cooling power is
    structurally larger than Ames Lab's own whole-device SPD (different
    mass denominators), so the two should never converge to a similar
    magnitude by coincidence -- if they ever do, something about the
    default design parameters or the SPD calculation itself likely
    changed in a way worth re-examining before trusting this module's own
    'these are not comparable' framing."""
    result = run_ames_lab_architecture_check()
    assert result["model_specific_cooling_power_w_per_kg_MCM"] > \
        2 * result["ames_optimized_SPD_w_per_kg_whole_device"], (
        "model_specific_cooling_power_w_per_kg_MCM is no longer clearly "
        "larger than Ames Lab's own whole-device SPD -- re-check whether "
        "this module's own 'not directly comparable, different "
        "denominators' framing is still accurate.")
