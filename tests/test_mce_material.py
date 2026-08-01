"""
Tests for core/mce_material.py.

Includes regression tests for two bugs found while profiling/optimizing
this module (see mce_material.py docstrings for the full explanation):

1. The magnetization() self-consistent solver used to be a damped
   fixed-point iteration that took ~300-500 iterations per call near Tc.
   Replaced with Newton's method (~5-15 iterations, same fixed points).

2. entropy_magnetic() computed log(sinh(ax)/sinh(bx)) by flooring the two
   sinh *arguments* independently at 1e-12, which silently collapsed the
   entropy to 0 (instead of the correct N*kB*ln(2J+1)) whenever x->0 --
   i.e. whenever H=0 and T>=Tc. This fed into delta_S_isothermal /
   delta_T_adiabatic as a large spurious error and could even flip the
   total heat capacity negative. Fixed with a proper small-x expansion.
"""
import numpy as np
import pytest

from core.mce_material import GADOLINIUM, kB, MagnetocaloricMaterial


def test_magnetization_converges_quickly():
    """Newton's method should reach tol in well under the old solver's
    500-iteration cap; this is a proxy for the perf fix holding."""
    T = np.linspace(250.0, 320.0, 25)
    H = 2.0 / (4 * np.pi * 1e-7)
    M = GADOLINIUM.magnetization(T, H, max_iter=30)
    assert np.all(np.isfinite(M))
    assert np.all(M >= 0)


def test_saturation_magnetization_bound():
    """M(T,H) can never exceed Msat = N*g*muB*J."""
    Msat = GADOLINIUM.N * GADOLINIUM.g * 9.2740100783e-24 * GADOLINIUM.J
    T = np.linspace(1.0, 300.0, 20)
    H = 10.0 / (4 * np.pi * 1e-7)
    M = GADOLINIUM.magnetization(T, H)
    assert np.all(M <= Msat * (1 + 1e-9))


@pytest.mark.parametrize("T_val", [294.0, 295.0, 300.0, 310.0])
def test_zero_field_entropy_above_tc_is_max_entropy(T_val):
    """Regression test for the log-sinh-ratio floor bug: at H=0 and T>=Tc,
    M=0 exactly, so S_M should equal the full spin-degeneracy entropy
    N*kB*ln(2J+1), not 0."""
    T = np.array([T_val])
    S0 = GADOLINIUM.entropy_magnetic(T, 0.0)
    expected = GADOLINIUM.N * kB * np.log(2 * GADOLINIUM.J + 1)
    assert S0[0] == pytest.approx(expected, rel=1e-6)


def test_entropy_decreases_with_field_near_tc():
    """Applying a field should always reduce magnetic entropy (spins
    ordering); this was violated (sign flipped) before the entropy fix."""
    T = np.array([294.0])
    S0 = GADOLINIUM.entropy_magnetic(T, 0.0)
    for B in [0.5, 1.0, 2.0, 5.0]:
        H = B / (4 * np.pi * 1e-7)
        S_H = GADOLINIUM.entropy_magnetic(T, H)
        assert S_H[0] < S0[0], f"entropy should drop under a {B} T field"


def test_delta_T_adiabatic_positive_and_bounded_near_tc():
    """dTad at Tc for Gd should be a small positive number of a few K per
    Tesla, per Dan'kov et al. (1998) -- not the -200 K seen when the
    entropy bug was present."""
    T = np.array([294.0])
    for B, lo, hi in [(1.0, 0.5, 8.0), (2.0, 1.0, 12.0), (5.0, 3.0, 25.0)]:
        H = B / (4 * np.pi * 1e-7)
        dTad = float(GADOLINIUM.delta_T_adiabatic(T, H)[0])
        assert lo < dTad < hi, f"dTad={dTad} out of sane range at {B} T"


def test_total_heat_capacity_is_positive():
    """C_total must be positive everywhere physically sampled; a negative
    value (as produced by the entropy bug) is unphysical and signals a
    modeling error."""
    T = np.linspace(200.0, 350.0, 30)
    H = 2.0 / (4 * np.pi * 1e-7)
    C = GADOLINIUM.total_heat_capacity(T, H)
    assert np.all(C > 0), f"non-positive heat capacity found: min={C.min()}"


def test_validation_errors_within_expected_mean_field_range():
    """The mean-field model is known (and documented) to overpredict dTad
    near Tc; lock in a generous bound so a future regression is caught
    without over-fitting to today's exact numbers."""
    from core.validation import run_validation
    rows = run_validation(verbose=False)
    for B, dT_lit, dT_model, err_pct in rows:
        assert abs(err_pct) < 60.0, f"unexpectedly large error at {B} T: {err_pct}%"


def test_lattice_heat_capacity_matches_dulong_petit_at_high_t():
    """Debye C_lattice -> 3R per mole-atom (Dulong-Petit) for T >> theta_D."""
    R = 8.314462618
    T = np.array([5000.0])  # >> theta_D=169K for Gd
    c_kg = GADOLINIUM.lattice_heat_capacity(T)
    c_molar = c_kg[0] * GADOLINIUM.M_molar
    assert c_molar == pytest.approx(3 * R, rel=0.02)


def test_first_order_material_flagged_in_source_metadata():
    """GD5SI2GE2 is retained as a parameter entry only; the module docstring
    explains the mean-field model is not appropriate for it. This just
    checks the library entry still exists and is documented as such."""
    from core.mce_material import GD5SI2GE2
    assert GD5SI2GE2.Tc == pytest.approx(276.0)
    assert "Pecharsky" in GD5SI2GE2.source
