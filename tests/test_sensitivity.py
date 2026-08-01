"""
Tests for core/sensitivity.py's Sobol global sensitivity analysis. Uses a
small n_base throughout to keep the suite fast; assertions target the
qualitative, seed-robust findings the module's own docstring/Sobol-stage
write-up rests on, not tight numeric values (Sobol indices from a small
sample are inherently noisy).
"""
import pytest

from core.sensitivity import run_sobol, model_cop, PROBLEM

N_BASE = 16  # small but enough for the qualitative checks below to be stable


def test_problem_definition_has_five_bounded_parameters():
    assert PROBLEM["num_vars"] == 5
    assert len(PROBLEM["names"]) == 5
    assert len(PROBLEM["bounds"]) == 5
    for lo, hi in PROBLEM["bounds"]:
        assert lo < hi


def test_model_cop_constant_loss_mode_is_positive():
    # midpoint of each parameter's bounds
    params = [0.5 * (lo + hi) for lo, hi in PROBLEM["bounds"]]
    cop = model_cop(params, use_state_dependent_losses=False)
    assert cop > 0


def test_model_cop_state_dependent_mode_is_positive():
    params = [0.5 * (lo + hi) for lo, hi in PROBLEM["bounds"]]
    cop = model_cop(params, use_state_dependent_losses=True)
    assert cop > 0


def test_constant_loss_model_shows_near_zero_sensitivity_to_field_frequency_flow(tmp_path):
    """Central Sobol-stage finding: with a constant parasitic_fraction,
    COP_electrical is algebraically independent of mu0H_max/frequency/
    fluid_mdot (they only change Qc, which cancels out of the COP
    formula), so their total-order Sobol indices should be ~0."""
    Si = run_sobol(n_base=N_BASE, use_state_dependent_losses=False, seed=1,
                    out_path=str(tmp_path / "sobol_constant_scratch.txt"))
    names = PROBLEM["names"]
    for pname in ["mu0H_max_T", "frequency_Hz", "fluid_mdot_kgs"]:
        st = Si["ST"][names.index(pname)]
        assert st < 0.05


def test_state_dependent_model_restores_sensitivity_to_field_frequency_flow(tmp_path):
    """Companion finding: once eddy/pumping/base-overhead losses are made
    state-dependent, mu0H_max/frequency/fluid_mdot should carry real
    (nonzero, seed-robust) sensitivity, unlike in the constant-loss case."""
    Si = run_sobol(n_base=N_BASE, use_state_dependent_losses=True, seed=1,
                    out_path=str(tmp_path / "sobol_state_scratch.txt"))
    names = PROBLEM["names"]
    total_st = sum(Si["ST"][names.index(p)]
                   for p in ["mu0H_max_T", "frequency_Hz", "fluid_mdot_kgs"])
    assert total_st > 0.05


def test_state_dependent_model_ignores_parasitic_fraction_parameter(tmp_path):
    """parasitic_fraction is documented as used only by the constant-loss
    formulation; with the state-dependent loss model it should be a
    structural no-op, giving an exactly-zero (not just small) ST."""
    Si = run_sobol(n_base=N_BASE, use_state_dependent_losses=True, seed=1,
                    out_path=str(tmp_path / "sobol_state_scratch2.txt"))
    names = PROBLEM["names"]
    assert Si["ST"][names.index("parasitic_fraction")] == pytest.approx(0.0, abs=1e-9)


def test_run_sobol_writes_output_file(tmp_path):
    out_path = tmp_path / "sobol_results.txt"
    run_sobol(n_base=N_BASE, out_path=str(out_path), use_state_dependent_losses=False)
    assert out_path.exists()
    text = out_path.read_text()
    assert "Sobol sensitivity analysis" in text
    for name in PROBLEM["names"]:
        assert name in text
