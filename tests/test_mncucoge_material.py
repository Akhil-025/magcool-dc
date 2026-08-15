"""
Tests for core/first_order_mce.py's Phase 25 addition: MNCUCOGE_FIRST_ORDER /
mncucoge_composition_tuned_material() (Mn1-xCuxCoGe magnetostructural family).
"""
import numpy as np
import pytest

from core.first_order_mce import (
    MNCUCOGE_FIRST_ORDER, MNCUCOGE_TC_MIN_K, MNCUCOGE_TC_MAX_K,
    mncucoge_composition_tuned_material,
)


def test_reference_material_tc_matches_samanta_x080():
    """MNCUCOGE_FIRST_ORDER is the x=0.080 composition -- Tc should match
    Samanta et al. (2012)'s own reported 302K exactly (it is the
    literature value itself)."""
    assert MNCUCOGE_FIRST_ORDER.Tc == pytest.approx(302.0)


def test_calibration_reproduces_literature_peak_entropy_change():
    """The (A,B,C) Landau coefficients were grid-searched to reproduce
    Samanta et al.'s own |DeltaS_M|=52.5 J/(kg K) at 5T target. Confirm
    the peak (scanned across T, not evaluated only at nominal Tc, per
    this module's own documented convention -- see __main__ block) lands
    within 1% of that target."""
    mu0 = 4 * np.pi * 1e-7
    H = 5.0 / mu0
    Ts = np.linspace(280.0, 330.0, 251)
    dS = MNCUCOGE_FIRST_ORDER.delta_S_isothermal(Ts, H, 0.0)
    peak = np.max(np.abs(dS))
    assert peak == pytest.approx(52.5, rel=0.01)


def test_peak_lands_above_nominal_tc():
    """Same qualitative model behavior as the module's honesty flag #3
    for the other three families: the Landau-model peak should sit above
    (not at) the nominal Tc."""
    mu0 = 4 * np.pi * 1e-7
    H = 5.0 / mu0
    Ts = np.linspace(280.0, 330.0, 251)
    dS = MNCUCOGE_FIRST_ORDER.delta_S_isothermal(Ts, H, 0.0)
    T_peak = Ts[int(np.argmax(np.abs(dS)))]
    assert T_peak > MNCUCOGE_FIRST_ORDER.Tc


def test_tc_window_bounds():
    assert MNCUCOGE_TC_MIN_K == pytest.approx(291.0)
    assert MNCUCOGE_TC_MAX_K == pytest.approx(316.0)


@pytest.mark.parametrize("tc", [291.0, 302.0, 316.0])
def test_tuned_material_hits_requested_tc(tc):
    mat = mncucoge_composition_tuned_material(tc)
    assert mat.Tc == pytest.approx(tc)


def test_tuned_material_preserves_other_parameters():
    mat = mncucoge_composition_tuned_material(305.0)
    assert mat.J == MNCUCOGE_FIRST_ORDER.J
    assert mat.g == MNCUCOGE_FIRST_ORDER.g
    assert mat.M_molar == pytest.approx(MNCUCOGE_FIRST_ORDER.M_molar)
    assert mat.theta_D == MNCUCOGE_FIRST_ORDER.theta_D
    assert mat.n_atoms_per_fu == MNCUCOGE_FIRST_ORDER.n_atoms_per_fu
    assert mat.A == MNCUCOGE_FIRST_ORDER.A
    assert mat.B == MNCUCOGE_FIRST_ORDER.B
    assert mat.C == MNCUCOGE_FIRST_ORDER.C
    assert mat.hysteresis_loss_J_per_kg == MNCUCOGE_FIRST_ORDER.hysteresis_loss_J_per_kg


def test_out_of_range_tc_raises():
    with pytest.raises(ValueError):
        mncucoge_composition_tuned_material(MNCUCOGE_TC_MIN_K - 1.0)
    with pytest.raises(ValueError):
        mncucoge_composition_tuned_material(MNCUCOGE_TC_MAX_K + 1.0)


def test_boundary_tc_values_do_not_raise():
    mncucoge_composition_tuned_material(MNCUCOGE_TC_MIN_K)
    mncucoge_composition_tuned_material(MNCUCOGE_TC_MAX_K)


def test_delta_T_adiabatic_finite():
    mat = mncucoge_composition_tuned_material(302.0)
    mu0 = 4 * np.pi * 1e-7
    H = 2.0 / mu0
    T = np.linspace(280.0, 330.0, 101)
    dTad = mat.delta_T_adiabatic(T, H)
    assert np.all(np.isfinite(dTad))


def test_named_override():
    mat = mncucoge_composition_tuned_material(300.0, name="custom label")
    assert mat.name == "custom label"
    assert mat.Tc == pytest.approx(300.0)
