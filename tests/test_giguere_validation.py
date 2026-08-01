"""
Regression tests for core/giguere_validation.py.

Locks down the ~2.4x discrepancy between the first-order Landau model's
predicted DeltaT_ad and Giguere et al.'s (1999) directly measured value,
so a future recalibration of first_order_mce.py's (A, B, C) or theta_D
can't silently drift this correction factor without a test noticing.
"""
import pytest

from core.giguere_validation import (
    run_validation, DTAD_CORRECTION_FACTOR,
    GIGUERE_DIRECT_DTAD_7T, GIGUERE_CLAUSIUS_CLAPEYRON_DTAD_7T,
    GIGUERE_INDIRECT_MAXWELL_DTAD_7T, _model_peak_dTad,
)


def test_giguere_reference_values_match_the_paper():
    """These are read directly from Giguere et al. (1999)'s Fig. 3 and
    Clausius-Clapeyron cross-check -- pin them so a future edit can't
    silently substitute a different number."""
    assert GIGUERE_DIRECT_DTAD_7T == 10.0
    assert GIGUERE_CLAUSIUS_CLAPEYRON_DTAD_7T == pytest.approx(9.9)
    assert GIGUERE_INDIRECT_MAXWELL_DTAD_7T == pytest.approx(14.9)


def test_model_overestimates_direct_measurement_by_documented_factor():
    """The module's docstring reports a ~2.4x overestimate of the model's
    peak DeltaT_ad at 7T vs. Giguere's direct measurement. Pin this to a
    generous-but-bounded range so silent drift (e.g. from an unrelated
    change to first_order_mce.py's calibration) gets caught."""
    _, dT_model = _model_peak_dTad(7.0)
    overestimate_factor = dT_model / GIGUERE_DIRECT_DTAD_7T
    assert overestimate_factor == pytest.approx(2.4, rel=0.15)


def test_model_overestimate_exceeds_papers_own_indirect_vs_direct_gap():
    """Central honesty-flag claim: the model's overestimate vs. Giguere's
    DIRECT measurement should be worse than the ~1.49x gap the paper
    itself found between its own indirect (Maxwell-relation) and direct
    values -- the two effects are additive, not duplicates."""
    _, dT_model = _model_peak_dTad(7.0)
    overestimate_vs_direct = dT_model / GIGUERE_DIRECT_DTAD_7T
    papers_own_overestimate = GIGUERE_INDIRECT_MAXWELL_DTAD_7T / GIGUERE_DIRECT_DTAD_7T
    assert overestimate_vs_direct > papers_own_overestimate


def test_dtad_correction_factor_brings_model_to_direct_measurement():
    """DTAD_CORRECTION_FACTOR is derived so that model_dTad * factor ==
    the directly measured value at 7T -- this is the whole point of the
    correction, so it should hold exactly (up to floating point), not
    approximately."""
    _, dT_model = _model_peak_dTad(7.0)
    assert dT_model * DTAD_CORRECTION_FACTOR == pytest.approx(
        GIGUERE_DIRECT_DTAD_7T, rel=1e-9)


def test_dtad_correction_factor_is_between_zero_and_one():
    """The model overestimates DeltaT_ad, so the correction that brings it
    down to the direct measurement must shrink it (0 < factor < 1)."""
    assert 0.0 < DTAD_CORRECTION_FACTOR < 1.0


def test_run_validation_writes_output_and_returns_matching_summary(tmp_path):
    out_path = tmp_path / "giguere_validation.txt"
    summary = run_validation(out_path=str(out_path), verbose=False)
    assert out_path.exists()
    assert "Giguere" in out_path.read_text()

    assert summary["direct_dTad_7T_K"] == GIGUERE_DIRECT_DTAD_7T
    assert summary["indirect_dTad_7T_K"] == GIGUERE_INDIRECT_MAXWELL_DTAD_7T
    assert summary["overestimate_factor_vs_direct"] == pytest.approx(
        summary["model_dTad_7T_K"] / GIGUERE_DIRECT_DTAD_7T)
