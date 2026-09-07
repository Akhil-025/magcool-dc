"""
Tests for core/mnfepsi_doped_hysteresis_speculative.py.

This module is explicitly SPECULATIVE (see its own docstring) -- these
tests do not assert that the fitted k or the resulting J/kg estimates are
"correct" in any physical sense (there is no ground truth to check them
against; that is the entire reason the module exists as a flagged
estimate rather than a new MNFEPSI_FIRST_ORDER default). What they DO
check: the fit against Zhang et al.'s own 5-point table is reproducible
and internally consistent, the two V-doped estimates come out in the
qualitative range the module's own bottom-line claim commits to (well
under MNFEPSI_FIRST_ORDER's 25.0 J/kg placeholder), and the module does
not mutate anything in core/first_order_mce.py as a side effect.
"""
import pytest

from core.first_order_mce import MNFEPSI_FIRST_ORDER
from core.mnfepsi_doped_hysteresis_speculative import (
    _ZHANG_TABLE, _fit_k, K_MEAN, K_STD, K_LEAST_SQUARES,
    VONFE_ESTIMATE, VONMN_ESTIMATE, run_estimate,
)


def test_fit_k_reproducible_and_matches_module_level_constants():
    """_fit_k() is a pure function of the hardcoded _ZHANG_TABLE -- calling
    it again must reproduce exactly the K_MEAN/K_STD/K_LEAST_SQUARES
    values computed at import time."""
    k_mean, k_std, k_ls = _fit_k()
    assert k_mean == pytest.approx(K_MEAN)
    assert k_std == pytest.approx(K_STD)
    assert k_ls == pytest.approx(K_LEAST_SQUARES)


def test_fit_k_matches_documented_values():
    """The module docstring quotes mean k=0.151 (std 0.028) and
    least-squares k=0.162 -- check the live fit still lands there,
    so a silent edit to _ZHANG_TABLE doesn't drift the printed claim."""
    assert K_MEAN == pytest.approx(0.151, abs=0.001)
    assert K_STD == pytest.approx(0.028, abs=0.001)
    assert K_LEAST_SQUARES == pytest.approx(0.162, abs=0.001)


def test_zhang_table_has_five_points():
    assert len(_ZHANG_TABLE) == 5


@pytest.mark.parametrize("est", [VONFE_ESTIMATE, VONMN_ESTIMATE])
def test_estimates_are_far_below_placeholder_and_flagged_as_extrapolation(est):
    """Both module docstring bottom-line claims to check for real: (a)
    each estimate is well under MNFEPSI_FIRST_ORDER's 25.0 J/kg
    placeholder, and (b) each is correctly flagged as a far
    extrapolation (T_hys well below Zhang et al.'s own measured range)."""
    assert 0.0 < est.estimate_J_per_kg_mean_k < 5.0
    assert 0.0 < est.estimate_J_per_kg_ls_k < 5.0
    assert est.estimate_J_per_kg_mean_k < MNFEPSI_FIRST_ORDER.hysteresis_loss_J_per_kg
    assert est.is_far_extrapolation is True


def test_mean_k_and_least_squares_k_give_consistent_order_of_magnitude():
    """The two fitting methods (mean-of-ratios vs least-squares-through-
    origin) shouldn't diverge wildly from each other for the same input --
    checks the module's own claimed ~5-10% spread between the two."""
    for est in (VONFE_ESTIMATE, VONMN_ESTIMATE):
        ratio = est.estimate_J_per_kg_ls_k / est.estimate_J_per_kg_mean_k
        assert 0.9 < ratio < 1.2


def test_run_estimate_does_not_mutate_mnfepsi_first_order():
    """This is a speculative, unwired module -- it must never change the
    validated MNFEPSI_FIRST_ORDER.hysteresis_loss_J_per_kg default as a
    side effect of running."""
    before = MNFEPSI_FIRST_ORDER.hysteresis_loss_J_per_kg
    run_estimate(verbose=False)
    after = MNFEPSI_FIRST_ORDER.hysteresis_loss_J_per_kg
    assert before == after


def test_run_estimate_returns_expected_shape():
    result = run_estimate(verbose=False)
    assert set(result.keys()) == {"k_mean", "k_std", "k_least_squares", "estimates"}
    assert len(result["estimates"]) == 2
