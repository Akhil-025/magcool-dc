"""
Tests for core/nanocomposite_material.py .
"""
import numpy as np
import pytest

from core.first_order_mce import (
    lafesih_composition_tuned_material, LAFESIH_TC_MIN_K, LAFESIH_TC_MAX_K,
)
from core.cascade import _target_composition_for_peak, _peak_temperature, run_cascade
from core.nanocomposite_material import (
    WeightedMaterialEnsemble, nanocomposite_tuned_material, NANOCOMPOSITE_FAMILY,
    NANOCOMPOSITE_SPREAD_K, NANOCOMPOSITE_TC_MIN_K, NANOCOMPOSITE_TC_MAX_K,
    run_robustness_check,
)

mu0 = 4 * np.pi * 1e-7


def test_weighted_ensemble_rejects_bad_weights():
    m1 = lafesih_composition_tuned_material(285.0)
    m2 = lafesih_composition_tuned_material(290.0)
    with pytest.raises(ValueError):
        WeightedMaterialEnsemble([m1, m2], [0.5, 0.6])
    with pytest.raises(ValueError):
        WeightedMaterialEnsemble([m1, m2], [1.0])


def test_weighted_ensemble_tc_is_weighted_mean():
    m1 = lafesih_composition_tuned_material(280.0)
    m2 = lafesih_composition_tuned_material(290.0)
    ens = WeightedMaterialEnsemble([m1, m2], [0.25, 0.75])
    assert ens.Tc == pytest.approx(0.25 * 280.0 + 0.75 * 290.0)


def test_weighted_ensemble_single_material_reduces_to_that_material():
    """A degenerate 1-material, weight-1.0 ensemble must reproduce that
    material's own delta_T_adiabatic exactly."""
    m = lafesih_composition_tuned_material(285.0)
    ens = WeightedMaterialEnsemble([m], [1.0])
    T = np.array([290.0, 295.0, 300.0])
    H = 2.0 / mu0
    np.testing.assert_allclose(
        np.asarray(ens.delta_T_adiabatic(T, H)).ravel(),
        np.asarray(m.delta_T_adiabatic(T, H)).ravel(),
    )


def test_weighted_ensemble_delta_T_adiabatic_is_linear_combination():
    """delta_T_adiabatic mixing is explicitly a weighted SUM of each
    phase's own delta_T_adiabatic (see class docstring) -- verify this
    directly against a hand-computed combination, not just via the class."""
    m1 = lafesih_composition_tuned_material(283.0)
    m2 = lafesih_composition_tuned_material(288.0)
    m3 = lafesih_composition_tuned_material(293.0)
    weights = (0.2, 0.5, 0.3)
    ens = WeightedMaterialEnsemble([m1, m2, m3], weights)
    T = np.array([294.0])
    H = 2.0 / mu0
    expected = sum(w * float(np.asarray(m.delta_T_adiabatic(T, H)).ravel()[0])
                    for w, m in zip(weights, [m1, m2, m3]))
    actual = float(np.asarray(ens.delta_T_adiabatic(T, H)).ravel()[0])
    assert actual == pytest.approx(expected, rel=1e-9)


def test_weighted_ensemble_hysteresis_loss_is_weighted_average():
    m1 = lafesih_composition_tuned_material(283.0)
    m2 = lafesih_composition_tuned_material(293.0)
    ens = WeightedMaterialEnsemble([m1, m2], [0.4, 0.6])
    expected = 0.4 * m1.hysteresis_loss_J_per_kg + 0.6 * m2.hysteresis_loss_J_per_kg
    assert ens.hysteresis_loss_J_per_kg == pytest.approx(expected)


def test_nanocomposite_tuned_material_builds_three_phases_at_expected_tc():
    mat = nanocomposite_tuned_material(290.0, spread_K=4.0)
    assert len(mat.materials) == 3
    tcs = sorted(m.Tc for m in mat.materials)
    assert tcs == pytest.approx([286.0, 290.0, 294.0])


def test_nanocomposite_tc_window_is_tightened_by_spread():
    assert NANOCOMPOSITE_TC_MIN_K == pytest.approx(LAFESIH_TC_MIN_K + NANOCOMPOSITE_SPREAD_K)
    assert NANOCOMPOSITE_TC_MAX_K == pytest.approx(LAFESIH_TC_MAX_K - NANOCOMPOSITE_SPREAD_K)


def test_nanocomposite_family_stays_within_lafesih_window_at_boundary():
    """Building a nanocomposite at the family's own documented tc_min/tc_max
    boundary must not raise (its outermost phase should land exactly at
    LAFESIH's own boundary, not beyond it)."""
    nanocomposite_tuned_material(NANOCOMPOSITE_TC_MIN_K)
    nanocomposite_tuned_material(NANOCOMPOSITE_TC_MAX_K)


def test_nanocomposite_peak_temperature_is_monotonic_in_center_tc():
    """Mirrors the monotonicity check core/cascade.py's own
    _target_composition_for_peak() docstring says was verified numerically
    for GD_FAMILY/LAFESIH_FAMILY -- required for the bracketed root-finder
    to converge correctly for NANOCOMPOSITE_FAMILY too."""
    centers = np.linspace(NANOCOMPOSITE_TC_MIN_K + 5, NANOCOMPOSITE_TC_MAX_K - 5, 6)
    peaks = [_peak_temperature(nanocomposite_tuned_material(c), 2.0) for c in centers]
    assert np.all(np.diff(peaks) > 0)


def test_target_composition_for_peak_converges_for_nanocomposite_family():
    T_target = 296.15
    tc = _target_composition_for_peak(T_target, 2.0, NANOCOMPOSITE_FAMILY)
    mat = NANOCOMPOSITE_FAMILY.tuned_fn(tc)
    peak_T = _peak_temperature(mat, 2.0)
    assert peak_T == pytest.approx(T_target, abs=0.5)


def test_nanocomposite_underperforms_perfectly_tuned_single_phase_at_design_point():
    """Documents the expected trade-off (blending costs peak height) at
    the design span itself, using the SAME Tc for both -- a regression
    guard on the qualitative direction of the effect, not just its
    existence."""
    Tc = 285.4
    nano = nanocomposite_tuned_material(Tc)
    single = lafesih_composition_tuned_material(Tc)
    T = np.array([296.15])
    H = 2.0 / mu0
    nano_dT = float(np.asarray(nano.delta_T_adiabatic(T, H)).ravel()[0])
    single_dT = float(np.asarray(single.delta_T_adiabatic(T, H)).ravel()[0])
    assert nano_dT <= single_dT


def test_run_robustness_check_returns_well_formed_result(tmp_path):
    out_path = tmp_path / "nanocomposite_robustness.txt"
    result = run_robustness_check(out_path=str(out_path), verbose=False)
    assert out_path.exists()
    assert len(result["rows"]) == 4
    assert isinstance(result["conclusion"], str) and len(result["conclusion"]) > 0
    design_rows = [r for r in result["rows"] if r["is_design_span"]]
    assert len(design_rows) == 1


def test_run_cascade_accepts_nanocomposite_material_directly():
    """core.cascade.run_cascade (and, through it, AMRSystem) must accept a
    WeightedMaterialEnsemble the same way it accepts any other material --
    no special-casing required, confirming the minimal interface
    (delta_T_adiabatic + optional hysteresis_loss_J_per_kg) is sufficient."""
    mat = nanocomposite_tuned_material(285.4)
    result = run_cascade(291.15, 10.0, 1, material=mat, mu0H_max=2.0, mass_per_stage=5.0)
    assert result["feasible"]
    assert result["Qc_W"] > 0