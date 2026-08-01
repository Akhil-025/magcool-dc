import numpy as np

from core.loss_model import (
    calibrate_loss_coefficients, leave_one_out_cv, CALIBRATION_POINTS_CORE,
    CALIBRATION_POINTS_EXTENDED, CALIBRATION_POINTS_FURTHER_EXTENDED,
    analyze_parasitic_fraction_scaling,
)


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
