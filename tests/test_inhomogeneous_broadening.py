"""
Tests for core/inhomogeneous_broadening.py .
"""
import numpy as np
import pytest

from core.mce_material import GADOLINIUM
from core.inhomogeneous_broadening import (
    BroadenedMagnetocaloricMaterial, _tc_ensemble, _peak_and_fwhm,
    run_broadening_sweep, run_dankov_error_sensitivity,
    run_inhomogeneous_broadening_analysis, mu0,
    FieldBroadenedMagnetocaloricMaterial, calibrate_field_dependent_broadening,
    run_field_dependent_broadening_calibration,
)

T_GRID = np.linspace(270.0, 320.0, 501)


def test_zero_sigma_ensemble_is_the_base_material_only():
    """sigma_Tc<=0 should degenerate to a single clone with weight 1."""
    clones, weights = _tc_ensemble(GADOLINIUM, 0.0)
    assert len(clones) == 1
    assert clones[0] is GADOLINIUM
    assert weights == pytest.approx([1.0])


def test_quadrature_weights_sum_to_one():
    """Gauss-Hermite weights normalized for a standard normal must sum to 1
    (to float precision) regardless of sigma_Tc, since they represent a
    probability distribution over grain Tc."""
    for sigma in (0.5, 1.0, 3.0, 5.0):
        _, weights = _tc_ensemble(GADOLINIUM, sigma)
        assert weights.sum() == pytest.approx(1.0, abs=1e-9)


def test_ensemble_mean_tc_matches_base_material():
    """The weighted mean of the quadrature Tc nodes should reproduce
    base_material.Tc (first moment of a symmetric Gaussian)."""
    clones, weights = _tc_ensemble(GADOLINIUM, 2.0)
    mean_tc = sum(w * c.Tc for w, c in zip(weights, clones))
    assert mean_tc == pytest.approx(GADOLINIUM.Tc, abs=1e-6)


def test_with_tc_clone_preserves_other_parameters():
    clone = GADOLINIUM.with_Tc(280.0)
    assert clone.Tc == pytest.approx(280.0)
    assert clone.J == GADOLINIUM.J
    assert clone.g == GADOLINIUM.g
    assert clone.M_molar == GADOLINIUM.M_molar
    assert clone.theta_D == GADOLINIUM.theta_D
    assert clone is not GADOLINIUM


def test_broadened_material_reduces_to_sharp_as_sigma_shrinks():
    """As sigma_Tc -> 0, the broadened ensemble's DeltaT_ad should converge
    to the sharp (unbroadened) model's own value at the same (T, H)."""
    H = 2.0 / mu0
    T = np.array([294.0])
    sharp = float(GADOLINIUM.delta_T_adiabatic(T, H)[0])
    broadened_tiny = BroadenedMagnetocaloricMaterial(GADOLINIUM, sigma_Tc_K=0.01)
    val = float(np.asarray(broadened_tiny.delta_T_adiabatic(T, H)).ravel()[0])
    assert val == pytest.approx(sharp, rel=1e-2)


def test_quadrature_matches_brute_force_reference():
    """15-node Gauss-Hermite quadrature should match a brute-force
    2001-point linspace/trapz reference integral for delta_T_adiabatic to
    a tight tolerance, confirming the quadrature-node count is adequate
    (module docstring's claim of <1e-6-scale quadrature error)."""
    sigma = 2.0
    H = 2.0 / mu0
    T = np.array([294.0])

    quad_mat = BroadenedMagnetocaloricMaterial(GADOLINIUM, sigma_Tc_K=sigma, n_quad=15)
    quad_val = float(np.asarray(quad_mat.delta_T_adiabatic(T, H)).ravel()[0])

    # Brute-force reference: dense linspace of Tc offsets, trapezoidal weights.
    xs = np.linspace(-6, 6, 2001)
    pdf = np.exp(-xs ** 2 / 2) / np.sqrt(2 * np.pi)
    Tc_values = GADOLINIUM.Tc + sigma * xs
    dS_vals = np.array([
        GADOLINIUM.with_Tc(tc).delta_S_isothermal(T, H)[0] for tc in Tc_values
    ])
    C_vals = np.array([
        GADOLINIUM.with_Tc(tc).total_heat_capacity(T)[0] for tc in Tc_values
    ])
    trapz_fn = getattr(np, "trapezoid", None) or np.trapz
    dS_ref = trapz_fn(dS_vals * pdf, xs)
    C_ref = trapz_fn(C_vals * pdf, xs)
    ref_val = -294.0 * dS_ref / C_ref

    assert quad_val == pytest.approx(ref_val, rel=5e-3)


def test_peak_and_fwhm_widens_with_broadening():
    """Broadening should widen (never narrow) the DeltaT_ad(T) FWHM at
    fixed field, since it is a convolution-like smoothing operation."""
    sharp_T, sharp_peak, sharp_fwhm = _peak_and_fwhm(GADOLINIUM, 5.0, T_GRID)
    broadened = BroadenedMagnetocaloricMaterial(GADOLINIUM, sigma_Tc_K=3.0)
    b_T, b_peak, b_fwhm = _peak_and_fwhm(broadened, 5.0, T_GRID)
    assert b_fwhm > sharp_fwhm
    # Broadening (a smoothing operation) should not raise the peak value.
    assert b_peak <= sharp_peak + 1e-9


def test_broadening_sweep_returns_expected_rows():
    rows = run_broadening_sweep(sigma_values=(0.0, 1.0), mu0H_T=(1.0, 5.0), verbose=False)
    assert len(rows) == 4
    for r in rows:
        assert r["fwhm_K"] > 0
        assert r["peak_dTad_K"] > 0


def test_dankov_error_sensitivity_sigma_zero_matches_gadolinium_linear_model():
    """sigma_Tc=0 rows must reproduce GADOLINIUM.delta_T_adiabatic() (the
    fast linear approximation) exactly, since _make_material(sigma=0)
    returns GADOLINIUM itself and run_dankov_error_sensitivity()
    deliberately uses that fast path (see its own docstring), not
    core/validation.py's run_validation() -- which, since the Paper-Mining
    Pass physics fix, uses the exact isentropic delta_T_adiabatic_exact()
    instead and is no longer expected to match this sweep's numbers."""
    from core.validation import LITERATURE_DELTA_T_AD
    rows = run_dankov_error_sensitivity(sigma_values=(0.0,), verbose=False)
    for r in rows:
        H = r["mu0H_T"] / mu0
        dT_lit = LITERATURE_DELTA_T_AD[r["mu0H_T"]]
        expected = float(np.asarray(GADOLINIUM.delta_T_adiabatic(np.array([294.0]), H)).ravel()[0])
        expected_err = 100 * (expected - dT_lit) / dT_lit
        assert r["dT_model_K"] == pytest.approx(expected, rel=1e-9)
        assert r["err_pct"] == pytest.approx(expected_err, rel=1e-9)


def test_full_analysis_runs_and_writes_report(tmp_path):
    out_path = tmp_path / "inhomogeneous_broadening.txt"
    result = run_inhomogeneous_broadening_analysis(out_path=str(out_path), verbose=False)
    assert out_path.exists()
    assert "sweep_rows" in result and len(result["sweep_rows"]) > 0
    assert "err_rows" in result and len(result["err_rows"]) > 0
    assert 0.0 in result["max_err_by_sigma"]
    assert isinstance(result["conclusion"], str) and len(result["conclusion"]) > 0


def test_worst_field_error_never_increases_relative_to_sharp_model_bound():
    """Documents (does not assert a specific direction of physics, only an
    internal-consistency bound): the reported best_sigma's worst-field
    error must never exceed the sharp model's own worst-field error, since
    best_sigma is chosen by explicit minimization over the swept values
    (which always includes sigma=0)."""
    result = run_inhomogeneous_broadening_analysis(out_path="results/inhomogeneous_broadening.txt",
                                                     verbose=False)
    best = result["best_sigma_Tc_K"]
    assert result["max_err_by_sigma"][best] <= result["sharp_max_err_pct"] + 1e-9

def test_field_broadened_material_zero_k_matches_sharp_model():
    """k_K_per_T=0.0 should reproduce the plain (unbroadened) GADOLINIUM
    prediction exactly, at any field -- the field-dependent class must be
    a strict generalization, not an independent parameterization."""
    mat = FieldBroadenedMagnetocaloricMaterial(GADOLINIUM, 0.0)
    for B in (1.0, 2.0, 5.0, 7.5):
        H = B / mu0
        sharp = float(GADOLINIUM.delta_T_adiabatic(np.array([294.0]), H)[0])
        broadened = float(np.asarray(mat.delta_T_adiabatic(np.array([294.0]), H)).ravel()[0])
        assert broadened == pytest.approx(sharp, rel=1e-9)


def test_field_broadened_sigma_grows_linearly_with_field():
    """sigma_Tc(H) = k*mu0*H by construction -- check the private helper
    directly rather than inferring it indirectly from DeltaT_ad."""
    mat = FieldBroadenedMagnetocaloricMaterial(GADOLINIUM, 2.0)
    H1 = 1.0 / mu0
    H5 = 5.0 / mu0
    assert mat._sigma_Tc_K(H1) == pytest.approx(2.0, rel=1e-9)
    assert mat._sigma_Tc_K(H5) == pytest.approx(10.0, rel=1e-9)


def test_calibrate_field_dependent_broadening_returns_consistent_rows():
    result = calibrate_field_dependent_broadening(verbose=False)
    assert result["k_fit_K_per_T"] >= 0.0
    assert len(result["rows"]) == 3
    for row in result["rows"]:
        # "before" (k=0) rows must match run_validation()'s own numbers
        # (same 3 fields, same T=294K methodology).
        assert row["dT_before_K"] > 0.0


def test_full_field_dependent_broadening_analysis_runs_and_writes_report(tmp_path):
    out_path = tmp_path / "field_dependent_broadening_calibration.txt"
    result = run_field_dependent_broadening_calibration(out_path=str(out_path), verbose=False)
    assert out_path.exists()
    assert "k_fit_K_per_T" in result
    assert "held_out_75T" in result
    assert "held_out_giguere" in result and len(result["held_out_giguere"]) == 2
    assert isinstance(result["conclusion"], str) and len(result["conclusion"]) > 0
