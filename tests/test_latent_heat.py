"""
Tests for core/first_order_mce.py's Phase 26 addition: latent_heat_J_per_kg /
latent_heat_width_K / latent_heat_capacity() / _field_dependent_transition_T() /
GD5SI2GE2_FIRST_ORDER_LATENT_HEAT, and the corresponding
core/giguere_validation.py::run_latent_heat_validation().
"""
import numpy as np
import pytest

from core.first_order_mce import (
    GD5SI2GE2_FIRST_ORDER, GD5SI2GE2_FIRST_ORDER_LATENT_HEAT,
)

mu0 = 4 * np.pi * 1e-7


def test_default_materials_unaffected_by_latent_heat_addition():
    """Every pre-existing FirstOrderMCEMaterial instance defaults to
    latent_heat_J_per_kg=0.0 -- delta_T_adiabatic() must be bit-for-bit
    unchanged (this is the module's own stated design goal for this
    field)."""
    Ts = np.linspace(260.0, 300.0, 801)
    dT = GD5SI2GE2_FIRST_ORDER.delta_T_adiabatic(Ts, 7.0 / mu0)
    peak = float(np.max(dT))
    # Documented/verified value prior to this phase's addition.
    assert peak == pytest.approx(24.169816296242782, rel=1e-9)


def test_latent_heat_capacity_zero_when_default():
    T = np.linspace(200.0, 350.0, 50)
    C = GD5SI2GE2_FIRST_ORDER.latent_heat_capacity(T, T_center=276.0)
    assert np.all(C == 0.0)


def test_latent_heat_capacity_integrates_to_L():
    """Gaussian normalization sanity check: integrating latent_heat_capacity
    over a wide T range should recover latent_heat_J_per_kg."""
    mat = GD5SI2GE2_FIRST_ORDER_LATENT_HEAT
    T = np.linspace(mat.Tc - 50.0, mat.Tc + 50.0, 20001)
    C = mat.latent_heat_capacity(T, T_center=mat.Tc)
    trapz_fn = getattr(np, "trapezoid", None) or np.trapz
    integral = trapz_fn(C, T)
    assert integral == pytest.approx(mat.latent_heat_J_per_kg, rel=1e-3)


def test_field_dependent_transition_T_shifts_with_field():
    """The transition location at 7T should sit well above the zero-field
    Tc for GD5SI2GE2_FIRST_ORDER -- this is the whole reason the latent
    heat spike must NOT be fixed at self.Tc (see module docstring)."""
    mat = GD5SI2GE2_FIRST_ORDER
    T_transition = mat._field_dependent_transition_T(7.0 / mu0)
    assert T_transition > mat.Tc + 5.0


def test_latent_heat_reduces_peak_dTad_at_7T():
    """The literature-grounded latent-heat instance should show a real,
    substantial reduction in peak DeltaT_ad at 7T relative to the raw
    model -- confirmed value from this phase's own validation run."""
    Ts = np.linspace(260.0, 300.0, 1601)
    dT_raw = GD5SI2GE2_FIRST_ORDER.delta_T_adiabatic(Ts, 7.0 / mu0)
    dT_lh = GD5SI2GE2_FIRST_ORDER_LATENT_HEAT.delta_T_adiabatic(Ts, 7.0 / mu0)
    peak_raw = float(np.max(dT_raw))
    peak_lh = float(np.max(dT_lh))
    assert peak_lh < peak_raw
    assert peak_lh == pytest.approx(19.011424164678456, rel=1e-6)
    # Real improvement, but NOT a full fix to Giguere's 10.0K direct target.
    assert peak_lh > 10.0


def test_latent_heat_instance_parameters():
    mat = GD5SI2GE2_FIRST_ORDER_LATENT_HEAT
    assert mat.latent_heat_J_per_kg == pytest.approx(276.0 * 18.0)
    assert mat.latent_heat_width_K == pytest.approx(5.0 / 2.3548, rel=1e-4)
    # Everything else should be inherited unchanged from GD5SI2GE2_FIRST_ORDER.
    assert mat.Tc == GD5SI2GE2_FIRST_ORDER.Tc
    assert mat.A == GD5SI2GE2_FIRST_ORDER.A
    assert mat.B == GD5SI2GE2_FIRST_ORDER.B
    assert mat.C == GD5SI2GE2_FIRST_ORDER.C
    assert mat.hysteresis_loss_J_per_kg == GD5SI2GE2_FIRST_ORDER.hysteresis_loss_J_per_kg


def test_total_heat_capacity_matches_lattice_only_when_no_latent_heat():
    T = np.linspace(260.0, 300.0, 21)
    lattice_only = GD5SI2GE2_FIRST_ORDER.lattice_heat_capacity(T)
    total = GD5SI2GE2_FIRST_ORDER.total_heat_capacity(T, H_final=7.0 / mu0)
    assert np.allclose(lattice_only, total)


def test_total_heat_capacity_exceeds_lattice_near_transition_with_latent_heat():
    mat = GD5SI2GE2_FIRST_ORDER_LATENT_HEAT
    H = 7.0 / mu0
    T_transition = mat._field_dependent_transition_T(H)
    T = np.array([T_transition])
    lattice_only = mat.lattice_heat_capacity(T)
    total = mat.total_heat_capacity(T, H_final=H)
    assert total[0] > lattice_only[0]


def test_run_latent_heat_validation_runs_and_returns_expected_keys():
    from core.giguere_validation import run_latent_heat_validation
    result = run_latent_heat_validation(verbose=False)
    assert "peak_dTad_7T_raw_K" in result
    assert "peak_dTad_7T_latent_heat_K" in result
    assert "gap_closed_pct" in result
    assert result["peak_dTad_7T_latent_heat_K"] < result["peak_dTad_7T_raw_K"]
    assert 0.0 < result["gap_closed_pct"] < 100.0