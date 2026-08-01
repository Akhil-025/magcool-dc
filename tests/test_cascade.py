import numpy as np
import pytest

from core.cascade import (
    run_cascade, run_graded_cascade, GD_FAMILY, LAFESIH_FAMILY,
    _target_composition_for_peak, _peak_temperature,
    validate_astronautics_graded_bed,
)
from core.mce_material import GADOLINIUM
from core.first_order_mce import (
    composition_tuned_material, lafesih_composition_tuned_material,
    GIANT_MCE_TC_MIN_K, GIANT_MCE_TC_MAX_K, LAFESIH_TC_MIN_K, LAFESIH_TC_MAX_K,
)


def test_run_cascade_more_stages_more_capacity():
    """More stages sharing the same total span should generally cost less
    work per stage (smaller span/stage -> easier lift), so 3-stage COP
    should not be dramatically worse than 1-stage for a fixed total span."""
    r1 = run_cascade(291.15, 12.0, 1, material=GADOLINIUM, mass_per_stage=5.0)
    r3 = run_cascade(291.15, 12.0, 3, material=GADOLINIUM, mass_per_stage=5.0)
    assert r1["feasible"] and r3["feasible"]
    assert r1["Qc_W"] > 0 and r3["Qc_W"] > 0


def test_peak_temperature_is_near_material_tc():
    """Sanity check for the two-pass coarse/fine _peak_temperature search
    added in Phase 9 for speed: the located peak should be within a few K
    of the material's own Tc (exact offset varies with field/model, but it
    should not be wildly off, e.g. not landing outside the search window
    or at a spurious low-T artifact)."""
    mat = composition_tuned_material(250.0, apply_giguere_correction=False)
    peak_T = _peak_temperature(mat, 2.0, T_range=(100.0, 330.0))
    assert abs(peak_T - mat.Tc) < 20.0


def test_target_composition_for_peak_gd_family_converges():
    """The brentq-based root-finder (Phase 9 replacement for the original
    fixed-point iteration, which was found to fail for the narrower
    LAFESIH_FAMILY transition) must still correctly solve for GD_FAMILY --
    this is the regression this test guards against (an earlier version of
    the Phase 9 change searched down to T=-20K, hit this Landau model's
    low-temperature DeltaT_ad numerical artifact, and returned garbage
    Tc~568K)."""
    T_target = 293.0
    Tc = _target_composition_for_peak(T_target, 2.0, GD_FAMILY)
    assert GIANT_MCE_TC_MIN_K <= Tc <= GIANT_MCE_TC_MAX_K
    mat = composition_tuned_material(Tc, apply_giguere_correction=True)
    peak_T = _peak_temperature(mat, 2.0, T_range=(100.0, GIANT_MCE_TC_MAX_K + 40.0))
    assert abs(peak_T - T_target) < 0.5


def test_target_composition_for_peak_lafesih_family_converges():
    """Same check for LAFESIH_FAMILY (Phase 9 addition) -- this is the
    specific case the original fixed-point iteration could NOT reliably
    solve (it left dTad at ~0.6K instead of ~21K at some stages due to the
    family's much narrower transition; see _target_composition_for_peak's
    docstring)."""
    T_target = 313.25
    Tc = _target_composition_for_peak(T_target, 1.44, LAFESIH_FAMILY)
    assert LAFESIH_TC_MIN_K <= Tc <= LAFESIH_TC_MAX_K
    mat = lafesih_composition_tuned_material(Tc)
    mu0 = 4 * np.pi * 1e-7
    dT_at_target = mat.delta_T_adiabatic(np.array([T_target]), 1.44 / mu0)[0]
    # this is the actual failure mode Phase 9 fixed: dT_at_target collapsing
    # to ~0.6K (span_fraction clamps Qc to 0) instead of landing near the
    # true peak (~21K for this material at 1.44T)
    assert dT_at_target > 15.0


def test_run_graded_cascade_gd_family_default_matches_explicit_family():
    """family=None must reproduce family=GD_FAMILY exactly (backward
    compatibility with the pre-Phase-9 API, which only supported the Gd
    family and took apply_giguere_correction directly)."""
    r_default = run_graded_cascade(291.15, 10.0, 3, mass_per_stage=5.0)
    r_explicit = run_graded_cascade(291.15, 10.0, 3, mass_per_stage=5.0, family=GD_FAMILY)
    assert r_default["feasible"] == r_explicit["feasible"]
    assert r_default["Qc_W"] == pytest.approx(r_explicit["Qc_W"], rel=1e-9)
    assert r_default["COP_cascade"] == pytest.approx(r_explicit["COP_cascade"], rel=1e-9)


def test_run_graded_cascade_gd_family_small_span_feasible_in_range():
    """At a small span near the ASHRAE data-center range, every stage's
    needed composition should stay within the documented giant-MCE window
    (no fallback to plain Gd needed) -- this is the "38 cells" majority
    case reported by the module's own __main__ sweep."""
    r = run_graded_cascade(291.15, 6.0, 2, mass_per_stage=5.0)
    assert r["feasible"]
    assert r["n_stages_out_of_range"] == 0
    assert r["Qc_W"] > 0
    assert r["COP_cascade"] > 0


def test_run_graded_cascade_lafesih_family_astronautics_range():
    """The 6-layer La(Fe,Si)13Hy grading used for the Astronautics device
    (Phase 9) should be fully in-range and feasible at that device's own
    operating point."""
    r = run_graded_cascade(305.0, 11.0, 6, mu0H_max=1.44, mass_per_stage=1.52 / 6,
                            frequency=4.0, family=LAFESIH_FAMILY)
    assert r["feasible"]
    assert r["n_stages_out_of_range"] == 0
    assert r["Qc_W"] > 0
    assert r["COP_cascade"] > 0
    for s in r["stage_info"]:
        assert s["in_range"]


def test_composition_tuned_material_out_of_range_raises():
    with pytest.raises(ValueError):
        composition_tuned_material(GIANT_MCE_TC_MAX_K + 50.0)
    with pytest.raises(ValueError):
        lafesih_composition_tuned_material(LAFESIH_TC_MIN_K - 50.0)


def test_validate_astronautics_graded_bed_reproduces_reported_qc_and_close_cop():
    """The Phase 9 headline result: calibrating fluid_mdot to reproduce the
    literature Qc=2502W, the 6-layer graded La(Fe,Si)13Hy bed's predicted
    COP should be within the same order-of-magnitude error the rest of
    validation_system.py sees for other devices (well under 50%), NOT the
    "no calibration found" outcome the single-layer LAFESIH_FIRST_ORDER
    material gave this same device."""
    r = validate_astronautics_graded_bed()
    assert r["feasible"]
    assert r["Qc_W"] == pytest.approx(2502.0, rel=1e-3)
    assert abs(r["COP_error_pct"]) < 50.0