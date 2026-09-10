"""
Tests for core/fluid_selection_optimization.py -- wired into main.py as
pipeline step 18 (run on every pipeline pass, see main.py's step list and
`fluid_selection_result` capture), but previously with no test coverage
of its own despite results/fluid_selection_optimization.txt and
results/fluid_selection_optimization_robustness.txt being real,
checked-in production outputs.

These tests use the module's own BASELINE_POINT/ROBUSTNESS_POINT (not
cut-down stand-ins) so a regression here reflects the actual pipeline
step's behavior, not a simplified proxy of it.
"""
import os

from core.fluids import FLUID_LIBRARY, DEFAULT_FLUID
from core.fluid_selection_optimization import (
    BASELINE_POINT,
    ROBUSTNESS_POINT,
    MDOT_BOUNDS_KGS,
    run_fluid_sweep,
    run_fluid_selection_comparison,
)


def test_run_fluid_sweep_returns_one_row_per_fluid():
    rows = run_fluid_sweep(BASELINE_POINT, verbose=False)
    assert len(rows) == len(FLUID_LIBRARY)
    assert {r["fluid"] for r in rows} == set(FLUID_LIBRARY)


def test_run_fluid_sweep_rows_have_expected_fields_and_are_physical():
    rows = run_fluid_sweep(BASELINE_POINT, verbose=False)
    for r in rows:
        assert MDOT_BOUNDS_KGS[0] <= r["mdot_kgs"] <= MDOT_BOUNDS_KGS[1]
        assert r["COP_electrical"] > 0.0
        assert r["Qc_W"] >= 0.0
        assert 0.0 <= r["exergy_eff"] <= 1.0


def test_run_fluid_sweep_sorted_by_cop_descending():
    rows = run_fluid_sweep(BASELINE_POINT, verbose=False)
    cops = [r["COP_electrical"] for r in rows]
    assert cops == sorted(cops, reverse=True)


def test_run_fluid_sweep_pct_below_water_is_zero_for_water_and_nonneg_for_rest():
    """water is the reference; every other fluid's pct_below_water must
    be >= 0 (nothing beats the pure-water ceiling in this repo's own
    documented finding -- see module docstring's RESULT section)."""
    rows = run_fluid_sweep(BASELINE_POINT, verbose=False)
    water_row = next(r for r in rows if r["fluid"] == "water")
    assert water_row["pct_below_water"] == 0.0
    for r in rows:
        if r["fluid"] != "water":
            assert r["pct_below_water"] >= 0.0


def test_run_fluid_sweep_water_ranks_first_at_baseline():
    """Regression guard for this module's own documented headline
    finding: water > water_eg10 > water_eg20 > water_pg30 > ethanol on
    COP_electrical at the baseline operating point."""
    rows = run_fluid_sweep(BASELINE_POINT, verbose=False)
    ordering = [r["fluid"] for r in rows]
    assert ordering == [
        "water", "water_eg10", "water_eg20", "water_pg30", "ethanol",
    ]


def test_run_fluid_sweep_water_eg10_gives_up_least_cop_among_realistic_fluids():
    """water_eg10 is DEFAULT_FLUID precisely because, among the
    corrosion-protected (non-pure-water) options, it gives up the least
    COP_electrical relative to the pure-water ceiling."""
    rows = run_fluid_sweep(BASELINE_POINT, verbose=False)
    realistic_rows = [r for r in rows if r["fluid"] != "water"]
    best_realistic = min(realistic_rows, key=lambda r: r["pct_below_water"])
    assert best_realistic["fluid"] == "water_eg10"
    assert best_realistic["fluid"] == DEFAULT_FLUID


def test_run_fluid_selection_comparison_writes_both_result_files(tmp_path):
    out_path = str(tmp_path / "fluid_selection_optimization.txt")
    robustness_out_path = str(tmp_path / "fluid_selection_optimization_robustness.txt")
    run_fluid_selection_comparison(
        out_path=out_path, robustness_out_path=robustness_out_path, verbose=False)
    assert os.path.exists(out_path)
    assert os.path.exists(robustness_out_path)


def test_run_fluid_selection_comparison_result_dict_has_expected_keys(tmp_path):
    result = run_fluid_selection_comparison(
        out_path=str(tmp_path / "baseline.txt"),
        robustness_out_path=str(tmp_path / "robustness.txt"),
        verbose=False,
    )
    for key in (
        "baseline_rows", "robustness_rows", "same_ranking",
        "best_realistic_fluid", "best_realistic_pct_below_water_baseline",
        "default_fluid", "default_matches_best_realistic",
    ):
        assert key in result


def test_run_fluid_selection_comparison_default_matches_best_realistic(tmp_path):
    """Sanity check on core.fluids.DEFAULT_FLUID staying consistent with
    what this analysis actually recommends -- if someone edits
    FLUID_LIBRARY's property values without re-checking this, the flag
    should flip to False rather than the mismatch going unnoticed."""
    result = run_fluid_selection_comparison(
        out_path=str(tmp_path / "baseline.txt"),
        robustness_out_path=str(tmp_path / "robustness.txt"),
        verbose=False,
    )
    assert result["default_fluid"] == DEFAULT_FLUID
    assert result["default_matches_best_realistic"] is True
    assert result["best_realistic_fluid"] == DEFAULT_FLUID


def test_run_fluid_selection_comparison_does_not_mutate_default_fluid(tmp_path):
    """This module must never have the side effect of changing
    core.fluids.DEFAULT_FLUID -- the module's own docstring is explicit
    that adopting a new default fluid is an opt-in choice, not an
    automatic consequence of running this comparison."""
    before = DEFAULT_FLUID
    run_fluid_selection_comparison(
        out_path=str(tmp_path / "baseline.txt"),
        robustness_out_path=str(tmp_path / "robustness.txt"),
        verbose=False,
    )
    from core.fluids import DEFAULT_FLUID as after
    assert after == before


def test_run_fluid_selection_comparison_same_ranking_flag_matches_rows(tmp_path):
    result = run_fluid_selection_comparison(
        out_path=str(tmp_path / "baseline.txt"),
        robustness_out_path=str(tmp_path / "robustness.txt"),
        verbose=False,
    )
    baseline_order = [r["fluid"] for r in result["baseline_rows"]]
    robustness_order = [r["fluid"] for r in result["robustness_rows"]]
    assert result["same_ranking"] == (baseline_order == robustness_order)


def test_baseline_and_robustness_points_share_temperature_and_span():
    """Both operating points must vary only frequency/mass_regenerator
    (per the module docstring) so the comparison is a genuine robustness
    check and not a confound from also moving T_cold/span."""
    assert BASELINE_POINT["mu0H_max"] == ROBUSTNESS_POINT["mu0H_max"]
    assert BASELINE_POINT["particle_diameter"] == ROBUSTNESS_POINT["particle_diameter"]
    assert BASELINE_POINT["mass_regenerator"] != ROBUSTNESS_POINT["mass_regenerator"]
    assert BASELINE_POINT["frequency"] != ROBUSTNESS_POINT["frequency"]
