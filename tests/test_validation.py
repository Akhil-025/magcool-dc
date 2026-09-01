"""
Tests for core/validation.py, including the Paper-Mining Pass Part 2
extensions: the Giguere et al. (1999) pure-Gd cross-check (§2) and the
Dan'kov et al. (1998) Curie-point field-shift check (§3).
"""
import pytest

from core.validation import (
    run_validation, run_giguere_gd_extension, run_curie_shift_check,
    LITERATURE_DELTA_T_AD, GIGUERE_GD_CROSSCHECK,
    DANKOV_CURIE_SHIFT_RATE_K_PER_T, DANKOV_CURIE_SHIFT_FIELD_RANGE_T,
)


def test_run_validation_reproduces_existing_dankov_errors():
    """Regression guard for the pre-existing Dan'kov et al. checks --
    unrelated to this pass, but should not have been disturbed by it."""
    rows = run_validation(verbose=False)
    assert len(rows) == len(LITERATURE_DELTA_T_AD)
    for B, dT_lit, dT_model, err_pct in rows:
        assert dT_model > 0
        assert err_pct == pytest.approx(100 * (dT_model - dT_lit) / dT_lit)


def test_giguere_reference_values_match_the_paper():
    """Read directly from Giguere et al. (1999)'s own pure-Gd methods-section
    cross-check paragraph -- pin them so a future edit can't silently
    substitute a different number."""
    assert GIGUERE_GD_CROSSCHECK[5.0]["range_K"] == (10.5, 11.5)
    assert GIGUERE_GD_CROSSCHECK[7.0]["range_K"] == (12.0, 13.0)


def test_giguere_gd_extension_uses_same_fixed_temperature_methodology():
    """run_giguere_gd_extension() must evaluate dTad at the SAME fixed
    T=294K run_validation() uses, for direct comparability -- not a
    peak-scanned value, which would silently change the comparison."""
    from core.mce_material import GADOLINIUM
    import numpy as np
    mu0 = 4 * 3.141592653589793 * 1e-7
    rows = run_giguere_gd_extension(verbose=False)
    for B, lo, hi, dT_model, err_pct, in_range in rows:
        expected = float(GADOLINIUM.delta_T_adiabatic(np.array([294.0]), B / mu0)[0])
        assert dT_model == pytest.approx(expected, rel=1e-9)


def test_giguere_gd_extension_model_overestimates_relative_to_giguere_range():
    """Documents the actual (honest, not-hidden) finding: a model
    calibrated to Dan'kov et al.'s 5T value (12.3K, Phase 32-corrected --
    see LITERATURE_DELTA_T_AD's own note) overestimates relative to
    Giguere et al.'s LOWER Gd range (10.5-11.5K) -- a modest residual
    cross-paper discrepancy, plausibly genuine sample/technique variation,
    not a bug, and should stay reported rather than silently reconciled.
    Note: this assertion's direction (model overestimates) is UNCHANGED by
    the Phase 32 correction -- what changed is the model's reported
    behavior at 5T specifically against LITERATURE_DELTA_T_AD (flipped
    from apparently underestimating by ~7.5% under the old, incorrect
    14.6K reference to overestimating by ~9.8% under the corrected 12.3K
    one -- see run_validation()'s own test for that comparison)."""
    rows = run_giguere_gd_extension(verbose=False)
    for B, lo, hi, dT_model, err_pct, in_range in rows:
        assert dT_model > hi
        assert not in_range
        assert err_pct > 0


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