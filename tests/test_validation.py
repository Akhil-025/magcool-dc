"""
Tests for core/validation.py, including the Paper-Mining Pass Part 2
extensions: the Giguere et al. (1999) pure-Gd cross-check (§2) and the
Dan'kov et al. (1998) Curie-point field-shift check (§3).
"""
import numpy as np
import pytest

from core.validation import (
    run_validation, run_giguere_gd_extension, run_curie_shift_check,
    LITERATURE_DELTA_T_AD, GIGUERE_GD_CROSSCHECK,
    DANKOV_CURIE_SHIFT_RATE_K_PER_T, DANKOV_CURIE_SHIFT_FIELD_RANGE_T,
    mu0,
)


def test_run_validation_structure():
    """Basic structural regression guard: right number of rows, positive
    predictions, error percentages consistent with the returned values."""
    rows = run_validation(verbose=False)
    assert len(rows) == len(LITERATURE_DELTA_T_AD)
    for B, dT_lit, dT_model, err_pct in rows:
        assert dT_model > 0
        assert err_pct == pytest.approx(100 * (dT_model - dT_lit) / dT_lit)


def test_run_validation_physics_fix_reduces_error_at_every_field():
    """The whole point of the Paper-Mining Pass physics fix (exact
    isentropic method + electronic term + fitted grain-Tc-broadening,
    core/inhomogeneous_broadening.py's GADOLINIUM_CALIBRATED): compared
    against the TRUE pre-fix baseline -- the original linear
    approximation with no electronic term and no broadening (sigma_Tc=0)
    -- every field's |error| against Dan'kov et al. (1998) should be
    smaller, and the worst remaining error should be well under half of
    what it was. (Note: GADOLINIUM itself now also carries the electronic
    term, since that fix applies unconditionally -- so the pre-fix
    baseline below is reconstructed explicitly with gamma=0 rather than
    read off plain GADOLINIUM, to isolate the FULL before/after
    comparison rather than an intermediate step.)"""
    import dataclasses
    from core.mce_material import GADOLINIUM
    GADOLINIUM_PRE_FIX = dataclasses.replace(GADOLINIUM, sommerfeld_gamma_J_per_molK2=0.0)

    rows = run_validation(verbose=False)
    old_errs = {}
    for B, dT_lit in LITERATURE_DELTA_T_AD.items():
        H = B / mu0
        dT_old = float(GADOLINIUM_PRE_FIX.delta_T_adiabatic(np.array([294.0]), H)[0])
        old_errs[B] = abs(100 * (dT_old - dT_lit) / dT_lit)
    worst_old = max(old_errs.values())
    worst_new = 0.0
    for B, dT_lit, dT_model, err_pct in rows:
        assert abs(err_pct) < old_errs[B]
        worst_new = max(worst_new, abs(err_pct))
    assert worst_new < 0.5 * worst_old


def test_giguere_reference_values_match_the_paper():
    """Read directly from Giguere et al. (1999)'s own pure-Gd methods-section
    cross-check paragraph -- pin them so a future edit can't silently
    substitute a different number."""
    assert GIGUERE_GD_CROSSCHECK[5.0]["range_K"] == (10.5, 11.5)
    assert GIGUERE_GD_CROSSCHECK[7.0]["range_K"] == (12.0, 13.0)


def test_giguere_gd_extension_uses_same_fixed_temperature_methodology():
    """run_giguere_gd_extension() must evaluate dTad at the SAME fixed
    T=294K AND the same calibrated model run_validation() uses
    (core/inhomogeneous_broadening.py's GADOLINIUM_CALIBRATED, Paper-
    Mining Pass physics fix), for direct comparability -- not a
    peak-scanned value or an uncalibrated model, either of which would
    silently change the comparison."""
    from core.inhomogeneous_broadening import GADOLINIUM_CALIBRATED
    mu0_local = 4 * 3.141592653589793 * 1e-7
    rows = run_giguere_gd_extension(verbose=False)
    for B, lo, hi, dT_model, err_pct, in_range in rows:
        expected = GADOLINIUM_CALIBRATED.delta_T_adiabatic_exact(294.0, B / mu0_local)
        assert dT_model == pytest.approx(expected, rel=1e-9)


def test_giguere_gd_extension_5T_falls_within_range_after_physics_fix():
    """With GADOLINIUM_CALIBRATED (exact isentropic method + electronic
    term + fitted grain-Tc-broadening -- see run_validation()'s
    docstring), the 5T prediction now falls INSIDE Giguere et al.'s
    independently-reported 10.5-11.5K range -- a genuine held-out
    success, since none of that fit used this dataset. 7T remains above
    Giguere et al.'s range, reported as such."""
    rows = {row[0]: row for row in run_giguere_gd_extension(verbose=False)}
    B5, lo5, hi5, dT5, err5, in_range5 = rows[5.0]
    assert in_range5
    assert lo5 <= dT5 <= hi5

    B7, lo7, hi7, dT7, err7, in_range7 = rows[7.0]
    assert not in_range7
    assert dT7 > hi7


def test_curie_shift_check_field_range_matches_dankov():
    """Dan'kov et al. (1998) report the ~6 K/T rate specifically over
    2-7.5T -- pin the checked range to that."""
    assert DANKOV_CURIE_SHIFT_FIELD_RANGE_T == (2.0, 7.5)
    assert DANKOV_CURIE_SHIFT_RATE_K_PER_T == 6.0


def test_curie_shift_check_returns_one_peak_per_field():
    result = run_curie_shift_check(verbose=False)
    assert len(result["fields_T"]) == len(result["peak_Ts_K"])
    assert len(result["fields_T"]) >= 10
    assert min(result["fields_T"]) == pytest.approx(2.0)
    assert max(result["fields_T"]) == pytest.approx(7.5)


def test_curie_shift_check_documents_the_actual_null_finding():
    """This is the genuine, documented finding (see run_curie_shift_check's
    own docstring): the model's peak temperature does NOT shift with field
    at anywhere near the literature's ~6 K/T rate -- it's pinned within a
    fraction of a Kelvin across the whole 2-7.5T range. Locks this down as
    a regression test so a future, unrelated change to mce_material.py
    can't silently start passing (or silently drift further) without
    someone noticing and updating the docstring/ROADMAP accordingly."""
    result = run_curie_shift_check(verbose=False)
    assert abs(result["fitted_slope_K_per_T"]) < 0.5
    assert abs(result["fitted_slope_K_per_T"]) < 0.1 * result["literature_slope_K_per_T"]
    peak_Ts = result["peak_Ts_K"]
    assert max(peak_Ts) - min(peak_Ts) < 0.5


def test_curie_shift_check_peaks_stay_near_nominal_tc():
    """Whatever the (lack of) field dependence, the peaks located should
    still sit near the material's own Tc=294K, not off in some spurious
    region of the search bounds."""
    from core.mce_material import GADOLINIUM
    result = run_curie_shift_check(verbose=False)
    for Tp in result["peak_Ts_K"]:
        assert abs(Tp - GADOLINIUM.Tc) < 5.0