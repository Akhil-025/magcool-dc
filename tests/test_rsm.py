"""
Tests for core/rsm.py's quadratic response-surface surrogate for AMR
cooling capacity (Qc). Uses small train/test sample counts throughout so
the suite stays fast -- these are structural/smoke checks, not a
recalibration of the surrogate's held-out accuracy.
"""
import numpy as np

from core.rsm import (
    sample_design, evaluate, build_quadratic_features, fit_rsm,
    VAR_NAMES, BOUNDS,
)


def test_sample_design_respects_bounds():
    X = sample_design(n=200, seed=1)
    assert X.shape == (200, len(VAR_NAMES))
    for i, (lo, hi) in enumerate(BOUNDS):
        assert X[:, i].min() >= lo
        assert X[:, i].max() <= hi


def test_sample_design_is_reproducible_given_seed():
    X1 = sample_design(n=50, seed=7)
    X2 = sample_design(n=50, seed=7)
    assert np.array_equal(X1, X2)


def test_evaluate_returns_nonnegative_cooling_capacity():
    X = sample_design(n=20, seed=3)
    Y = evaluate(X)
    assert Y.shape == (20,)
    assert np.all(Y >= 0)
    assert np.all(np.isfinite(Y))


def test_build_quadratic_features_shape_and_names():
    """5 variables -> 1 intercept + 5 linear + 5 squared + C(5,2)=10
    pairwise interaction terms = 21 columns/names total."""
    X = sample_design(n=10, seed=1)
    Phi, feat_names = build_quadratic_features(X, VAR_NAMES)
    expected_cols = 1 + len(VAR_NAMES) + len(VAR_NAMES) + 10
    assert Phi.shape == (10, expected_cols)
    assert len(feat_names) == expected_cols
    assert feat_names[0] == "intercept"
    assert np.array_equal(Phi[:, 0], np.ones(10))


def test_build_quadratic_features_squared_terms_match_input():
    X = sample_design(n=10, seed=1)
    Phi, feat_names = build_quadratic_features(X, VAR_NAMES)
    sq_idx = feat_names.index(f"{VAR_NAMES[0]}^2")
    assert np.allclose(Phi[:, sq_idx], X[:, 0] ** 2)


def test_fit_rsm_achieves_reasonable_held_out_fit(tmp_path):
    """A smoke-level accuracy bar, not a tight regression: the quadratic
    surrogate should explain a substantial majority of held-out variance
    in Qc, well above a random/constant-predictor baseline (R^2=0)."""
    out_path = tmp_path / "rsm_coefficients.txt"
    coeffs, feat_names, r2, rmse = fit_rsm(
        out_path=str(out_path), n_train=150, n_test=50)
    assert out_path.exists()
    assert np.isfinite(r2)
    assert r2 > 0.5
    assert rmse >= 0
    assert len(coeffs) == len(feat_names)


def test_fit_rsm_output_file_lists_all_coefficients(tmp_path):
    out_path = tmp_path / "rsm_coefficients.txt"
    coeffs, feat_names, _, _ = fit_rsm(
        out_path=str(out_path), n_train=120, n_test=40)
    text = out_path.read_text()
    for name in feat_names:
        assert name in text
