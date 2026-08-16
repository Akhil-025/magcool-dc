"""
Tests for core/antiperovskite_material.py (Phase 24: Ga1-xCMn3+x
composition-tunable antiperovskite family).
"""
import numpy as np
import pytest

from core.antiperovskite_material import (
    GA1XCMN3X_REF, GA1XCMN3X_TC_MIN_K, GA1XCMN3X_TC_MAX_K,
    ga1xcmn3x_composition_tuned_material,
)


def test_reference_material_matches_wang_et_al_x007():
    """GA1XCMN3X_REF is the x=0.07 composition -- Tc should match Wang et
    al. (2009)'s own reported 296.5K exactly (it is the literature value
    itself, not a fit)."""
    assert GA1XCMN3X_REF.Tc == pytest.approx(296.5)


def test_tc_window_matches_measured_endpoints():
    assert GA1XCMN3X_TC_MIN_K == pytest.approx(250.0)
    assert GA1XCMN3X_TC_MAX_K == pytest.approx(323.5)


@pytest.mark.parametrize("tc", [250.0, 281.5, 296.5, 323.5])
def test_tuned_material_reproduces_measured_points(tc):
    """The four Tc values Wang et al. actually measured (x=0/0.06/0.07/
    0.08) should be reproduced exactly by construction (Tc is the only
    tuned parameter)."""
    mat = ga1xcmn3x_composition_tuned_material(tc)
    assert mat.Tc == pytest.approx(tc)


def test_tuned_material_preserves_other_parameters():
    """Only Tc (and its derived Weiss constant) should change -- J, g,
    M_molar, theta_D, n_atoms_per_fu stay fixed at the REF calibration
    (same convention as MagnetocaloricMaterial.with_Tc() itself, and the
    other composition_tuned_material() helpers in this repo)."""
    mat = ga1xcmn3x_composition_tuned_material(310.0)
    assert mat.J == GA1XCMN3X_REF.J
    assert mat.g == GA1XCMN3X_REF.g
    assert mat.M_molar == pytest.approx(GA1XCMN3X_REF.M_molar)
    assert mat.theta_D == GA1XCMN3X_REF.theta_D
    assert mat.n_atoms_per_fu == GA1XCMN3X_REF.n_atoms_per_fu


def test_out_of_range_tc_raises():
    with pytest.raises(ValueError):
        ga1xcmn3x_composition_tuned_material(GA1XCMN3X_TC_MIN_K - 1.0)
    with pytest.raises(ValueError):
        ga1xcmn3x_composition_tuned_material(GA1XCMN3X_TC_MAX_K + 1.0)


def test_boundary_tc_values_do_not_raise():
    ga1xcmn3x_composition_tuned_material(GA1XCMN3X_TC_MIN_K)
    ga1xcmn3x_composition_tuned_material(GA1XCMN3X_TC_MAX_K)


def test_zero_hysteresis_second_order_material_has_no_hysteresis_field():
    """This family reuses core.mce_material.MagnetocaloricMaterial (the
    second-order/mean-field class), which -- unlike
    core.first_order_mce.FirstOrderMCEMaterial -- has no
    hysteresis_loss_J_per_kg field at all, consistent with Wang et al.'s
    own directly-reported finding of no observable hysteresis for this
    composition series."""
    mat = GA1XCMN3X_REF
    assert not hasattr(mat, "hysteresis_loss_J_per_kg")


def test_delta_T_adiabatic_finite_and_peaks_near_tc():
    """Sanity check: the tuned material's DeltaT_ad(T) should be finite
    across a T range spanning its Tc and roughly peak somewhere in that
    neighborhood (mean-field/second-order materials peak close to Tc,
    unlike the first-order Landau families' well-documented +10-11K
    offset -- see module docstring)."""
    mat = ga1xcmn3x_composition_tuned_material(296.5)
    mu0 = 4 * np.pi * 1e-7
    H = 2.0 / mu0
    T = np.linspace(270.0, 320.0, 101)
    dTad = mat.delta_T_adiabatic(T, H)
    assert np.all(np.isfinite(dTad))
    T_peak = T[int(np.argmax(dTad))]
    assert abs(T_peak - mat.Tc) < 15.0


def test_named_override():
    mat = ga1xcmn3x_composition_tuned_material(300.0, name="custom label")
    assert mat.name == "custom label"
    assert mat.Tc == pytest.approx(300.0)
