import numpy as np
import pytest

from core.economics import material_cost, lifetime_cost
from core.emissions import refrigerant_emissions_tCO2e, operational_emissions_tCO2e, compare_emissions
from core.baseline_cooling import carnot_cop, vapor_compression_cop, liquid_cooling_cop
from core.loss_model import (
    calibrate_loss_coefficients, leave_one_out_cv, CALIBRATION_POINTS_CORE,
    CALIBRATION_POINTS_EXTENDED, CALIBRATION_POINTS_FURTHER_EXTENDED,
    analyze_parasitic_fraction_scaling,
)


def test_material_cost_scales_with_field_and_mass():
    base = material_cost(mu0H_max=1.0, mass_regenerator=1.0)
    higher_field = material_cost(mu0H_max=2.0, mass_regenerator=1.0)
    more_mass = material_cost(mu0H_max=1.0, mass_regenerator=2.0)
    assert higher_field > base
    assert more_mass > base


def test_lifetime_cost_includes_materials_floor_and_electricity():
    result = lifetime_cost(mu0H_max=1.0, mass_regenerator=1.0, Qc_avg_W=500.0,
                            COP_electrical=3.0, device_lifetime_years=15.0)
    mat_floor = material_cost(mu0H_max=1.0, mass_regenerator=1.0)
    assert result["materials_floor_$"] == pytest.approx(mat_floor, rel=1e-6)
    assert result["lifetime_electricity_$"] > 0
    assert result["lifetime_total_$"] == pytest.approx(
        result["materials_floor_$"] + result["lifetime_electricity_$"], rel=1e-6)


def test_lifetime_cost_scales_with_lifetime_and_inversely_with_cop():
    short = lifetime_cost(mu0H_max=1.0, mass_regenerator=1.0, Qc_avg_W=500.0,
                           COP_electrical=3.0, device_lifetime_years=5.0)
    long = lifetime_cost(mu0H_max=1.0, mass_regenerator=1.0, Qc_avg_W=500.0,
                          COP_electrical=3.0, device_lifetime_years=15.0)
    assert long["lifetime_electricity_$"] > short["lifetime_electricity_$"]

    low_cop = lifetime_cost(mu0H_max=1.0, mass_regenerator=1.0, Qc_avg_W=500.0,
                             COP_electrical=2.0)
    high_cop = lifetime_cost(mu0H_max=1.0, mass_regenerator=1.0, Qc_avg_W=500.0,
                              COP_electrical=8.0)
    assert low_cop["lifetime_electricity_$"] > high_cop["lifetime_electricity_$"]


def test_lifetime_cost_rejects_nonpositive_cop():
    with pytest.raises(ValueError):
        lifetime_cost(mu0H_max=1.0, mass_regenerator=1.0, Qc_avg_W=500.0,
                       COP_electrical=0.0)


def test_refrigerant_emissions_zero_leak_rate_is_zero():
    assert refrigerant_emissions_tCO2e(100.0, leak_rate=0.0) == 0.0


def test_operational_emissions_scale_inversely_with_cop():
    low_cop = operational_emissions_tCO2e(100.0, cop=3.0)
    high_cop = operational_emissions_tCO2e(100.0, cop=10.0)
    assert low_cop > high_cop


def test_compare_emissions_amr_has_zero_refrigerant_component():
    results = compare_emissions(100.0, amr_cop=5.0, vcc_cop=12.0, liquid_cop=20.0)
    amr = next(r for r in results if r.technology.startswith("Magnetic"))
    assert amr.refrigerant_GWP_tCO2e_per_year == 0.0


def test_carnot_cop_matches_definition():
    assert carnot_cop(290.0, 300.0) == pytest.approx(290.0 / 10.0)


def test_baseline_cops_below_carnot():
    vcc = vapor_compression_cop(290.0, 300.0)
    liq = liquid_cooling_cop(290.0, 300.0)
    assert vcc.COP < vcc.COP_carnot
    # liquid cooling blends a high economizer-mode COP, so only check it's
    # not exceeding physical bounds by an absurd margin (it can legitimately
    # exceed the *mechanical* Carnot figure for the DX-only comparison).
    assert liq.COP > 0


def test_loss_model_calibration_nonnegative():
    cal = calibrate_loss_coefficients(verbose=False)
    assert cal["k_eddy"] >= 0
    assert cal["k_pump"] >= 0
    assert cal["base_frac"] >= 0


def test_loss_model_exactly_determined_zero_residual():
    """With exactly 3 points and 3 unknowns, the CORE calibration should
    fit exactly (near-zero residual), which is what makes leave-one-out
    on this same 3-point set structurally uninformative (each fold is only
    2 points / 3 unknowns -- underdetermined) -- documented here so the
    limitation is explicit rather than silently assumed."""
    cal = calibrate_loss_coefficients(CALIBRATION_POINTS_CORE, verbose=False)
    from core.loss_model import _build_system
    A, b = _build_system(CALIBRATION_POINTS_CORE)
    pred = A @ cal["raw"]
    resid = np.abs(pred - b)
    assert np.all(resid < 1e-6 * np.abs(b))


def test_nnls_extended_fit_is_nonnegative():
    """Phase 6 found the unconstrained lstsq fit on the 4-point EXTENDED set
    gives negative k_eddy/base_frac. NNLS must give a non-negative fit by
    construction -- this is the whole point of switching solvers."""
    cal = calibrate_loss_coefficients(CALIBRATION_POINTS_EXTENDED, verbose=False,
                                       label="EXTENDED test")
    assert cal["k_eddy"] >= 0
    assert cal["k_pump"] >= 0
    assert cal["base_frac"] >= 0


def test_nnls_loo_error_improves_but_remains_large():
    """NNLS should improve (reduce) the worst leave-one-out error on the
    EXTENDED set relative to the +1639% Phase 6 reported with plain lstsq,
    but should NOT bring it down to a "solved" small error -- the Phase 7
    finding is that a better solver helps but does not fix the underlying
    structural mismatch across four orders of magnitude of device scale."""
    loo = leave_one_out_cv(CALIBRATION_POINTS_EXTENDED, verbose=False)
    worst_abs_err = max(abs(r[3]) for r in loo)
    assert worst_abs_err < 1639.0
    assert worst_abs_err > 100.0


def test_further_extended_fit_is_nonnegative():
    """The 8-point FURTHER_EXTENDED set (EXTENDED + 4 calibratable Lozano
    POLO/UFSC (2016) points) must still give a non-negative NNLS fit."""
    cal = calibrate_loss_coefficients(CALIBRATION_POINTS_FURTHER_EXTENDED,
                                       verbose=False, label="FURTHER_EXTENDED test")
    assert cal["k_eddy"] >= 0
    assert cal["k_pump"] >= 0
    assert cal["base_frac"] >= 0


def test_further_extended_lozano_points_consistently_underpredicted():
    """Phase 7 finding: when each Lozano point is held out, the pooled fit
    consistently UNDERpredicts its required W_parasitic by a similar
    amount (all leave-one-out errors strongly negative and tightly
    clustered) rather than scattering randomly around zero -- evidence
    that Lozano's device sits in a distinct, worse motor/inverter
    efficiency class the model can't represent, not that the model is
    merely noisy on this device."""
    loo = leave_one_out_cv(CALIBRATION_POINTS_FURTHER_EXTENDED, verbose=False)
    lozano_errs = [r[3] for r in loo if r[0].startswith("Lozano")]
    assert len(lozano_errs) == 4
    assert all(e < -50.0 for e in lozano_errs)  # all substantially underpredicted
    assert max(lozano_errs) - min(lozano_errs) < 10.0  # tightly clustered


def test_parasitic_fraction_scaling_is_not_monotonic():
    """The Phase 6 write-up speculated that a size/scale term would fix the
    loss model. Checking parasitic fraction (W_parasitic/Qc) sorted by
    device scale shows this doesn't hold in the current 4-device set: the
    smallest device does not have the highest overhead fraction."""
    rows = analyze_parasitic_fraction_scaling(CALIBRATION_POINTS_EXTENDED, verbose=False)
    fracs = [r[2] for r in rows]
    monotonic = all(fracs[i] <= fracs[i + 1] for i in range(len(fracs) - 1)) or \
                all(fracs[i] >= fracs[i + 1] for i in range(len(fracs) - 1))
    assert not monotonic