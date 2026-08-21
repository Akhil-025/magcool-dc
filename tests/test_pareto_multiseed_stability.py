"""
Tests for core/pareto_multiseed_stability.py. Uses tiny pop_size/n_gen and
only 2 seeds to keep this a fast smoke test: checks the harness runs
end-to-end and produces a self-consistent summary dict -- NOT that any
particular COP/material-share number comes out a specific way (that's the
actual research question, answered by running the module at production
settings, not by this test suite).
"""
import pytest

from core.pareto_multiseed_stability import run_pareto_multiseed_stability_check

_POP, _GEN = 8, 4


def test_multiseed_check_runs_end_to_end(tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)
    import os
    os.makedirs("results", exist_ok=True)
    result = run_pareto_multiseed_stability_check(
        seeds=(1, 2), pop_size=_POP, n_gen=_GEN, verbose=False)
    assert "per_seed" in result and "summary" in result
    assert len(result["per_seed"]) == 2


def test_summary_has_expected_keys_and_consistent_types(tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)
    import os
    os.makedirs("results", exist_ok=True)
    result = run_pareto_multiseed_stability_check(
        seeds=(1, 2), pop_size=_POP, n_gen=_GEN, verbose=False)
    summary = result["summary"]
    if summary is None:
        pytest.skip("both seeds returned an empty front at this tiny pop_size/n_gen")
    for key in ("best_COP_electrical_mean", "best_COP_electrical_std",
                "knee_COP_electrical_mean", "knee_material_consistent_across_seeds",
                "material_share_stats"):
        assert key in summary
    assert isinstance(summary["knee_material_consistent_across_seeds"], bool)
    # material shares for each family should sum to ~100% per seed on average
    total_share = sum(s["mean_pct"] for s in summary["material_share_stats"].values())
    assert 99.0 <= total_share <= 101.0
