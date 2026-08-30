"""Phase 31 addition: core/uncertainty_propagation.py had no dedicated
test file even though it is wired into main.py's pipeline (see README's
Tier-1 test-coverage gap). n_draws is kept small (20-40) throughout so
this file runs in well under a second, rather than the module's own
default of 2000 draws intended for a real report."""

from core.uncertainty_propagation import (
    monte_carlo_cop_uncertainty,
    uncertainty_band_across_spans,
)


def test_monte_carlo_returns_a_sane_confidence_band():
    result = monte_carlo_cop_uncertainty(n_draws=30, seed=1, verbose=False)
    assert result.COP_electrical_mean > 0
    assert result.COP_electrical_std >= 0
    assert result.COP_electrical_p05 <= result.COP_electrical_mean <= result.COP_electrical_p95
    # NOTE (genuine finding, not a test-tolerance workaround): Qc's p05/p95
    # band is essentially zero-width here (differs only at ~1e-13, i.e.
    # floating-point noise) because Qc (cooling capacity) is set by the
    # AMR's thermal/NTU/mass-flow calculation, which this codebase's
    # architecture does NOT route the calibrated loss coefficients
    # through -- only COP_electrical (Qc / input electrical power) is
    # sensitive to the perturbed calibration. So a small floating-point
    # epsilon, not a strict "<=", is the correct check here; a strict
    # "Qc_p05 <= Qc_mean <= Qc_p95" occasionally fails on rounding alone.
    # Worth flagging to the project authors: the "Qc confidence band" in
    # uncertainty_band_across_spans() is not actually informative -- only
    # the COP_electrical band carries real calibration uncertainty.
    assert result.Qc_p05 <= result.Qc_mean + 1e-6
    assert result.Qc_mean <= result.Qc_p95 + 1e-6
    assert 0.0 <= result.fraction_failed_draws <= 1.0


def test_monte_carlo_is_reproducible_given_the_same_seed():
    a = monte_carlo_cop_uncertainty(n_draws=25, seed=42, verbose=False)
    b = monte_carlo_cop_uncertainty(n_draws=25, seed=42, verbose=False)
    assert a.COP_electrical_mean == b.COP_electrical_mean
    assert a.COP_electrical_std == b.COP_electrical_std
    assert a.Qc_mean == b.Qc_mean


def test_monte_carlo_different_seeds_give_different_draws():
    a = monte_carlo_cop_uncertainty(n_draws=25, seed=1, verbose=False)
    b = monte_carlo_cop_uncertainty(n_draws=25, seed=2, verbose=False)
    # Not a strict requirement of correctness, but if two different RNG
    # seeds produced bit-identical means it would suggest the seed isn't
    # actually being threaded through to the perturbation draws.
    assert a.COP_electrical_mean != b.COP_electrical_mean


def test_monte_carlo_wider_assumed_uncertainty_gives_wider_or_equal_confidence_band():
    tight = monte_carlo_cop_uncertainty(n_draws=200, measurement_uncertainty_frac=0.05,
                                         seed=7, verbose=False)
    wide = monte_carlo_cop_uncertainty(n_draws=200, measurement_uncertainty_frac=0.30,
                                        seed=7, verbose=False)
    tight_width = tight.COP_electrical_p95 - tight.COP_electrical_p05
    wide_width = wide.COP_electrical_p95 - wide.COP_electrical_p05
    assert wide_width >= tight_width


def test_monte_carlo_result_reports_the_requested_n_draws_and_frac():
    result = monte_carlo_cop_uncertainty(n_draws=17, measurement_uncertainty_frac=0.22, verbose=False)
    assert result.n_draws == 17
    assert result.measurement_uncertainty_frac == 0.22


def test_qc_is_essentially_invariant_to_calibration_uncertainty_unlike_cop():
    """Documents a genuine architectural finding surfaced while writing
    these tests: Qc (cooling capacity) does not depend on the calibrated
    loss-model coefficients in this codebase, only COP_electrical does
    (Qc / input electrical power). So perturbing the calibration inputs
    changes the COP distribution meaningfully but leaves Qc's distribution
    essentially at its point estimate (std ~1e-13, floating-point noise
    only). This means uncertainty_band_across_spans()'s Qc_p05/Qc_p95
    columns are not actually informative uncertainty bands -- only its
    COP_electrical_p05/p95 columns are. Worth a one-line caveat in the
    module's own report output if this is surfaced in a paper."""
    result = monte_carlo_cop_uncertainty(n_draws=100, measurement_uncertainty_frac=0.30,
                                          seed=3, verbose=False)
    cop_relative_spread = (result.COP_electrical_p95 - result.COP_electrical_p05) / result.COP_electrical_mean
    qc_relative_spread = (result.Qc_p95 - result.Qc_p05) / result.Qc_mean
    assert cop_relative_spread > 0.01     # COP band is meaningfully wide
    assert qc_relative_spread < 1e-6      # Qc band is not


def test_uncertainty_band_across_spans_returns_one_row_per_span():
    spans = [5, 10, 15]
    rows = uncertainty_band_across_spans(spans_K=spans, n_draws=20, verbose=False)
    assert len(rows) == len(spans)
    assert [r["span_K"] for r in rows] == spans


def test_uncertainty_band_across_spans_rows_have_expected_keys():
    rows = uncertainty_band_across_spans(spans_K=[10], n_draws=20, verbose=False)
    for key in ("span_K", "COP_electrical_mean", "COP_electrical_p05",
                "COP_electrical_p95", "Qc_mean_W", "Qc_p05_W", "Qc_p95_W"):
        assert key in rows[0]
