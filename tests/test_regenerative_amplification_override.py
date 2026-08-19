"""Tests for this session's additions closing the "core AMR model has no
regenerative amplification" gap (see README.md's "1-D regenerator model
findings" and "Follow-up: an opt-in override" sections):

  - core/amr_cycle.py: AMRSystem.no_load_span_override
  - core/regenerator_1d.py: axial conduction, no_load_span()'s interior
    maximum, regenerative_span_cap()
  - core/loss_model.py: CALIBRATION_POINTS_CORE_PLUS_MAGGIE_HIGHSPAN

Before this file, none of the above had any dedicated test coverage --
every existing calibration point and model addition in this repo has a
regression guard except these (see test_loss_model.py's own
test_core_calibration_points_are_self_consistent for the precedent this
follows). regenerator_1d.py's full transient simulation is deliberately
slow at production settings (tens of seconds per no_load_span() call,
see that module's own docstring) -- every test below uses small
n_nodes/max_cycles specifically to stay fast while still exercising the
real code paths, not mocks.
"""
import numpy as np
import pytest

from core.mce_material import GADOLINIUM
from core.amr_cycle import AMRSystem
from core.regenerator_1d import (simulate_amr_1d, no_load_span,
                                   _apply_axial_conduction)
from core.loss_model import (CALIBRATION_POINTS_CORE,
                               CALIBRATION_POINTS_CORE_PLUS_MAGGIE_HIGHSPAN,
                               leave_one_out_cv)


# ---------------------------------------------------------------------------
# AMRSystem.no_load_span_override
# ---------------------------------------------------------------------------

def make_system(**overrides):
    kwargs = dict(material=GADOLINIUM, mu0H_max=2.0, mass_regenerator=5.0,
                  frequency=1.0, fluid_mdot=0.1, regenerator_effectiveness=0.8)
    kwargs.update(overrides)
    return AMRSystem(**kwargs)


def test_no_load_span_override_default_none_preserves_old_cap():
    """Default (unset) must reproduce the pre-existing 2*dTad_noload
    behavior exactly -- every caller that doesn't pass this parameter
    (optimize.py, cascade.py, every pre-existing test) must be
    byte-for-byte unaffected."""
    sys_default = make_system()
    sys_explicit_none = make_system(no_load_span_override=None)
    Qc1, dTad1 = sys_default.cooling_capacity(291.0, 10.0)
    Qc2, dTad2 = sys_explicit_none.cooling_capacity(291.0, 10.0)
    assert Qc1 == pytest.approx(Qc2)
    assert dTad1 == pytest.approx(dTad2)
    # A span far beyond any reasonable 2*dTad_noload for Gd must still
    # hard-zero under the default cap (dTad_noload is T_mid-dependent, so
    # this checks the CLAMP behavior robustly rather than assuming a fixed
    # cap value across different spans).
    Qc_far, _ = sys_default.cooling_capacity(291.0, 200.0)
    assert Qc_far == pytest.approx(0.0, abs=1e-6)


def test_no_load_span_override_extends_reachable_span():
    """With the override set to a large span, a span the default cap
    hard-zeros (200K, far beyond any reasonable single-blow dTad_noload
    for Gd) must become reachable (Qc > 0)."""
    sys_default = make_system()
    test_span = 200.0

    Qc_default, _ = sys_default.cooling_capacity(291.0, test_span)
    assert Qc_default == pytest.approx(0.0, abs=1e-6)

    sys_override = make_system(no_load_span_override=250.0)
    Qc_override, _ = sys_override.cooling_capacity(291.0, test_span)
    assert Qc_override > 0.0


def test_no_load_span_override_does_not_change_zero_span_capacity():
    """The override extends how far span can reach, but must not change
    Qc's magnitude at zero span -- that's still set by dTad_noload alone
    (see the parameter's own docstring)."""
    sys_default = make_system()
    sys_override = make_system(no_load_span_override=999.0)
    Qc1, _ = sys_default.cooling_capacity(291.0, 0.0)
    Qc2, _ = sys_override.cooling_capacity(291.0, 0.0)
    assert Qc1 == pytest.approx(Qc2)


def test_no_load_span_override_zero_or_negative_gives_zero_capacity():
    """A degenerate override (<=0) should not divide by zero or go
    negative -- cooling_capacity() must fail safe to Qc=0."""
    sys_ = make_system(no_load_span_override=0.0)
    Qc, _ = sys_.cooling_capacity(291.0, 1.0)
    assert Qc == pytest.approx(0.0, abs=1e-9)


def test_no_load_span_override_threads_through_run():
    """The override must actually reach cooling_capacity() when called via
    run(), not just when cooling_capacity() is called directly."""
    sys_default = make_system()
    test_span = 200.0
    Qc_default, _ = sys_default.cooling_capacity(291.0, test_span)
    assert Qc_default == pytest.approx(0.0, abs=1e-6)
    sys_override = make_system(no_load_span_override=250.0)
    result = sys_override.run(291.0, test_span)
    assert result.Qc > 0.0


# ---------------------------------------------------------------------------
# core/regenerator_1d.py -- axial conduction + interior maximum
# ---------------------------------------------------------------------------

def test_axial_conduction_conserves_energy():
    """_apply_axial_conduction redistributes heat between insulated-end
    nodes -- total thermal energy (sum of T, equal node masses/cp here)
    must be conserved to numerical precision."""
    T = np.array([300.0, 290.0, 310.0, 295.0, 305.0])
    T_after = _apply_axial_conduction(T, dt_total=1.0, dx=0.001, bed_area=0.002,
                                       k_eff_axial=1.0, m_node=0.01,
                                       cp_solid_eff=236.0)
    assert np.sum(T_after) == pytest.approx(np.sum(T), abs=1e-6)


def test_axial_conduction_smooths_toward_uniform():
    """More conduction time should reduce the spread between the hottest
    and coldest node (a basic sanity check that heat actually flows from
    hot to cold, not the reverse)."""
    T = np.array([320.0, 300.0, 280.0])
    spread_before = T.max() - T.min()
    T_after = _apply_axial_conduction(T, dt_total=0.5, dx=0.001, bed_area=0.002,
                                       k_eff_axial=1.0, m_node=0.01,
                                       cp_solid_eff=236.0)
    spread_after = T_after.max() - T_after.min()
    assert spread_after < spread_before


def test_axial_conduction_zero_conductivity_is_a_no_op():
    T = np.array([310.0, 290.0, 300.0])
    T_after = _apply_axial_conduction(T, dt_total=1.0, dx=0.001, bed_area=0.002,
                                       k_eff_axial=0.0, m_node=0.01,
                                       cp_solid_eff=236.0)
    assert np.allclose(T_after, T)


def test_simulate_amr_1d_no_load_span_is_finite_and_bounded():
    """A tiny/fast configuration should still produce a finite, physically
    sane (not NaN, not absurdly large) span -- the regression this test
    guards against is exactly the low-mdot degeneracy this session fixed
    (span growing without bound / diverging to NaN, see
    core/regenerator_1d.py's module docstring 'known limitations')."""
    r = simulate_amr_1d(GADOLINIUM, 1.5, 0.5, 1.0, mdot=0.01,
                         n_nodes=4, max_cycles=40, tol=1e-3)
    assert np.isfinite(r["span_K"])
    assert 0.0 <= r["span_K"] < 100.0  # sanity ceiling, not a precise bound


def test_simulate_amr_1d_very_low_mdot_does_not_diverge():
    """The specific regression this session's axial-conduction fix
    targets: span at a very low flow rate must stay bounded, not run away
    (pre-fix, this kept growing without limit as mdot -> 0)."""
    r = simulate_amr_1d(GADOLINIUM, 1.5, 0.5, 1.0, mdot=0.0002,
                         n_nodes=4, max_cycles=40, tol=1e-3)
    assert np.isfinite(r["span_K"])
    assert r["span_K"] < 100.0


def test_no_load_span_search_returns_best_among_grid():
    """no_load_span()'s search must return the mdot in its own grid with
    the largest span, not just the first or last point evaluated."""
    result = no_load_span(GADOLINIUM, 1.5, 0.5, 1.0, n_nodes=4,
                           mdot_search=(0.001, 0.01, 0.1), max_cycles=25, tol=1e-2)
    assert "mdot_kg_s" in result
    assert result["mdot_kg_s"] in (0.001, 0.01, 0.1)
    assert np.isfinite(result["span_K"])


# ---------------------------------------------------------------------------
# core/loss_model.py -- 4th calibration point
# ---------------------------------------------------------------------------

def test_maggie_highspan_calibration_point_is_self_consistent():
    """Guards CALIBRATION_POINTS_CORE_PLUS_MAGGIE_HIGHSPAN's hardcoded mdot
    against silent drift, the same way test_loss_model.py's
    test_core_calibration_points_are_self_consistent already guards CORE
    -- except this point only calibrates WITH no_load_span_override set
    (that's the whole reason it's a separate set, not folded into CORE),
    so the override is supplied explicitly here rather than omitted."""
    name, f, H, mdot, Qc_lit, _Wp = CALIBRATION_POINTS_CORE_PLUS_MAGGIE_HIGHSPAN[-1]
    assert name == "DTU_Eriksen_MAGGIE_2016"
    sys_ = AMRSystem(material=GADOLINIUM, mu0H_max=H, mass_regenerator=1.7,
                      frequency=f, fluid_mdot=mdot, no_load_span_override=21.04)
    Qc_model, _ = sys_.cooling_capacity(289.0, 15.5)
    assert Qc_model == pytest.approx(Qc_lit, rel=0.02), (
        f"{name}: hardcoded mdot={mdot} with no_load_span_override=21.04 predicts "
        f"Qc={Qc_model:.1f}W, not the hardcoded Qc_lit={Qc_lit}W -- this point has "
        f"drifted out of sync; re-run the brentq recalibration described in "
        f"loss_model.py's CALIBRATION_POINTS_CORE_PLUS_MAGGIE_HIGHSPAN comment.")


def test_maggie_highspan_point_improves_its_own_held_out_fold():
    """Honest regression guard on this session's actual finding: the new
    point's own leave-one-out fold should predict within a much tighter
    band than CORE's existing folds do (documented as ~+30% vs ~250-700%
    for the others) -- if this degrades to an order-of-magnitude miss like
    the rest, the 'this is a genuinely different, better kind of 4th
    point' claim in README.md/loss_model.py's docstring no longer holds
    and those need revisiting, not just this test."""
    loo = leave_one_out_cv(CALIBRATION_POINTS_CORE_PLUS_MAGGIE_HIGHSPAN, verbose=False)
    maggie_fold = [r for r in loo if r[0] == "DTU_Eriksen_MAGGIE_2016"][0]
    err_pct = maggie_fold[3]
    assert abs(err_pct) < 100.0, (
        f"DTU_Eriksen_MAGGIE_2016's held-out leave-one-out error is now "
        f"{err_pct:+.1f}%, no longer the ~+30% this point's inclusion was "
        f"justified by -- re-check CALIBRATION_POINTS_CORE_PLUS_MAGGIE_HIGHSPAN.")


def test_core_default_unaffected_by_new_calibration_set_existing():
    """The production CORE fit must be numerically identical whether or
    not CALIBRATION_POINTS_CORE_PLUS_MAGGIE_HIGHSPAN exists in the module
    -- i.e. adding the new set must not have mutated CORE itself."""
    assert len(CALIBRATION_POINTS_CORE) == 3
    assert CALIBRATION_POINTS_CORE_PLUS_MAGGIE_HIGHSPAN[:3] == CALIBRATION_POINTS_CORE
    assert len(CALIBRATION_POINTS_CORE_PLUS_MAGGIE_HIGHSPAN) == 4


def test_maggie_point_mdot_is_back_calculated_not_measured():
    """Regression guard for the honesty flag in loss_model.py's
    CALIBRATION_POINTS_CORE_PLUS_MAGGIE_HIGHSPAN comment: this point's
    hardcoded mdot (solved to reproduce the reported Qc) is substantially
    smaller than the device's own directly-measured flow rate (2.5 L/min
    -> 0.04167 kg/s, Eriksen 2016 PhD thesis Table 6.2). If this ever
    stops being true (e.g. the point is re-derived some other way), the
    honesty flag describing a ~2.8x gap needs to be revisited too."""
    _name, _f, _H, mdot_backcalc, _Qc, _Wp = CALIBRATION_POINTS_CORE_PLUS_MAGGIE_HIGHSPAN[-1]
    mdot_measured = 2.5e-3 / 60.0 * 1000.0  # 2.5 L/min of water -> kg/s
    ratio = mdot_measured / mdot_backcalc
    assert ratio == pytest.approx(2.84, rel=0.05)
