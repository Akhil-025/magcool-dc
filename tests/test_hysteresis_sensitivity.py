"""
Tests for core/hysteresis_sensitivity.py's ON/OFF Pareto-front A/B
comparison. Uses much smaller pop_size/n_gen than the module's own
default (which is itself already reduced from optimize.py's production
default -- see that module's docstring honesty flag #1) to keep this a
fast smoke test: it checks that the harness runs end-to-end, restores the
mutated module-level constants correctly regardless of outcome, and
produces a self-consistent output file/dict -- NOT that the ON/OFF
material-composition difference is large or in any particular direction
(that's the actual research question, answered by running the module at
its own default or higher, not by this test suite).
"""
import os
import pytest

from core.first_order_mce import (
    GD5SI2GE2_FIRST_ORDER, LAFESIH_FIRST_ORDER, MNFEPSI_FIRST_ORDER,
)
from core.hysteresis_sensitivity import (
    run_hysteresis_sensitivity, _set_all_hysteresis, _material_counts,
    run_hysteresis_multiseed_stability_check,
)

_POP, _GEN = 14, 5  # small but >= n_var=7-ish reference-direction count concerns
                     # aren't relevant here (pymoo just warns, doesn't fail);
                     # chosen purely to keep this test fast (~1-2s/material).


def _hysteresis_snapshot():
    return {
        "Gd5Si2Ge2": GD5SI2GE2_FIRST_ORDER.hysteresis_loss_J_per_kg,
        "La(Fe,Si)13Hy": LAFESIH_FIRST_ORDER.hysteresis_loss_J_per_kg,
        "Mn-Fe-P-Si": MNFEPSI_FIRST_ORDER.hysteresis_loss_J_per_kg,
    }


def test_run_hysteresis_sensitivity_restores_original_values(tmp_path):
    """Regardless of outcome, the module-level *_FIRST_ORDER constants'
    hysteresis_loss_J_per_kg must be restored to their pre-call values
    after run_hysteresis_sensitivity() returns -- this is what makes it
    safe to call from within a larger test/analysis session without
    leaking mutated global state into unrelated later code."""
    before = _hysteresis_snapshot()
    out_path = str(tmp_path / "hysteresis_sensitivity.txt")
    run_hysteresis_sensitivity(pop_size=_POP, n_gen=_GEN, seed=1, out_path=out_path)
    after = _hysteresis_snapshot()
    assert before == after
    for v in before.values():
        assert v > 0.0, "fixture assumes the Phase 16 placeholders are nonzero by default"


def test_run_hysteresis_sensitivity_restores_values_even_on_exception(monkeypatch, tmp_path):
    """The try/finally in run_hysteresis_sensitivity() must restore
    original values even if the second (hysteresis-OFF) optimization run
    raises."""
    before = _hysteresis_snapshot()

    import core.hysteresis_sensitivity as hs
    calls = {"n": 0}
    real_run_optimization = hs.optimize.run_optimization

    def _flaky_run_optimization(*args, **kwargs):
        calls["n"] += 1
        if calls["n"] == 2:
            raise RuntimeError("simulated failure on the OFF run")
        return real_run_optimization(*args, **kwargs)

    monkeypatch.setattr(hs.optimize, "run_optimization", _flaky_run_optimization)

    out_path = str(tmp_path / "hysteresis_sensitivity.txt")
    with pytest.raises(RuntimeError, match="simulated failure"):
        run_hysteresis_sensitivity(pop_size=_POP, n_gen=_GEN, seed=1, out_path=out_path)

    after = _hysteresis_snapshot()
    assert before == after


def test_run_hysteresis_sensitivity_writes_output_file(tmp_path):
    out_path = str(tmp_path / "sub" / "hysteresis_sensitivity.txt")
    result = run_hysteresis_sensitivity(pop_size=_POP, n_gen=_GEN, seed=1, out_path=out_path)
    assert os.path.isfile(out_path)
    with open(out_path) as f:
        content = f.read()
    assert "OFF (pre-Ph16)" in content
    assert "ON (Ph16)" in content
    assert "rows_on" in result and "rows_off" in result
    assert "counts_on" in result and "counts_off" in result


def test_set_all_hysteresis_helper_zeroes_all_three_constants():
    before = _hysteresis_snapshot()
    try:
        _set_all_hysteresis(0.0)
        assert GD5SI2GE2_FIRST_ORDER.hysteresis_loss_J_per_kg == 0.0
        assert LAFESIH_FIRST_ORDER.hysteresis_loss_J_per_kg == 0.0
        assert MNFEPSI_FIRST_ORDER.hysteresis_loss_J_per_kg == 0.0
    finally:
        GD5SI2GE2_FIRST_ORDER.hysteresis_loss_J_per_kg = before["Gd5Si2Ge2"]
        LAFESIH_FIRST_ORDER.hysteresis_loss_J_per_kg = before["La(Fe,Si)13Hy"]
        MNFEPSI_FIRST_ORDER.hysteresis_loss_J_per_kg = before["Mn-Fe-P-Si"]


def test_material_counts_helper():
    rows = [{"material": "A"}, {"material": "B"}, {"material": "A"}]
    counts = _material_counts(rows)
    assert counts == {"A": 2, "B": 1}

# ---------------------------------------------------------------------------
# run_hysteresis_multiseed_stability_check() -- closes the "Open item" flagged
# in this module's own honesty flag #1 and in ROADMAP.md's Phase 16 entry.
# ---------------------------------------------------------------------------

def test_multiseed_stability_check_restores_original_values(tmp_path):
    before = _hysteresis_snapshot()
    out_path = str(tmp_path / "stability.txt")
    run_hysteresis_multiseed_stability_check(seeds=(1, 2), pop_size=_POP,
                                              n_gen=_GEN, out_path=out_path)
    after = _hysteresis_snapshot()
    assert before == after


def test_multiseed_stability_check_returns_one_row_per_seed(tmp_path):
    out_path = str(tmp_path / "stability.txt")
    result = run_hysteresis_multiseed_stability_check(
        seeds=(1, 2, 3), pop_size=_POP, n_gen=_GEN, out_path=out_path)
    assert len(result["per_seed"]) == 3
    assert [s["seed"] for s in result["per_seed"]] == [1, 2, 3]
    for s in result["per_seed"]:
        assert 0.0 <= s["lafesih_frac_off"] <= 1.0
        assert 0.0 <= s["lafesih_frac_on"] <= 1.0
    assert isinstance(result["stable"], bool)


def test_multiseed_stability_check_writes_output_file_and_cleans_scratch(tmp_path):
    out_path = str(tmp_path / "sub" / "stability.txt")
    run_hysteresis_multiseed_stability_check(
        seeds=(1,), pop_size=_POP, n_gen=_GEN, out_path=out_path)
    assert os.path.isfile(out_path)
    with open(out_path) as f:
        content = f.read()
    assert "RESULT:" in content
    assert "seed" in content
    # scratch CSVs the function writes internally must not be left behind
    assert not os.path.exists("results/_scratch_hysteresis_multiseed_on.csv")
    assert not os.path.exists("results/_scratch_hysteresis_multiseed_off.csv")


def test_multiseed_stability_check_stable_flag_matches_per_seed_data(tmp_path):
    out_path = str(tmp_path / "stability.txt")
    result = run_hysteresis_multiseed_stability_check(
        seeds=(1, 2), pop_size=_POP, n_gen=_GEN, out_path=out_path)
    expected_stable = all(
        s["lafesih_frac_on"] >= s["lafesih_frac_off"] - 1e-9
        for s in result["per_seed"])
    assert result["stable"] == expected_stable