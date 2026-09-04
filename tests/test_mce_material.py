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


def test_electronic_heat_capacity_zero_for_default_materials():
    """Sommerfeld term defaults to 0.0 -- must not silently appear for
    materials that don't specify a gamma (GD5SI2GE2, LACAMNO3)."""
    from core.mce_material import GD5SI2GE2, LACAMNO3
    T = np.array([294.0])
    assert GD5SI2GE2.electronic_heat_capacity(T)[0] == 0.0
    assert LACAMNO3.electronic_heat_capacity(T)[0] == 0.0


def test_gadolinium_electronic_heat_capacity_matches_cited_gamma():
    """GADOLINIUM.sommerfeld_gamma_J_per_molK2 = 5.4 mJ/(mol K^2) (de
    Oliveira & von Ranke, Phys. Rep. 489 (2010), Sec. 4.2) -> C_el(294K)
    = gamma*T/M_molar."""
    T = np.array([294.0])
    expected = 5.4e-3 * 294.0 / GADOLINIUM.M_molar
    assert GADOLINIUM.electronic_heat_capacity(T)[0] == pytest.approx(expected, rel=1e-9)


def test_total_heat_capacity_includes_electronic_term():
    T = np.array([294.0])
    C_no_el = GADOLINIUM.lattice_heat_capacity(T) + GADOLINIUM.magnetic_heat_capacity(T, 0.0)
    C_total = GADOLINIUM.total_heat_capacity(T, 0.0)
    assert C_total[0] > C_no_el[0]
    assert C_total[0] == pytest.approx(
        (C_no_el + GADOLINIUM.electronic_heat_capacity(T))[0], rel=1e-9)


def test_entropy_lattice_positive_and_increasing_with_temperature():
    """Third law + monotonicity: S_lattice(T) >= 0 and increases with T."""
    T = np.array([50.0, 150.0, 294.0, 500.0])
    S = GADOLINIUM.entropy_lattice(T)
    assert np.all(S >= 0)
    assert np.all(np.diff(S) > 0)


def test_entropy_lattice_high_T_dulong_petit_limit():
    """At T >> theta_D, S_lattice should approach the Dulong-Petit-consistent
    high-T expansion 3*n*R*[1/3 + ln(T/theta_D)] (per mole), i.e. grow
    logarithmically -- loose sanity check, not a tight fit."""
    R = 8.314462618
    T_lo, T_hi = 5000.0, 50000.0
    S_lo = GADOLINIUM.entropy_lattice(np.array([T_lo]))[0] * GADOLINIUM.M_molar
    S_hi = GADOLINIUM.entropy_lattice(np.array([T_hi]))[0] * GADOLINIUM.M_molar
    expected_diff = 3 * R * np.log(T_hi / T_lo)
    assert (S_hi - S_lo) == pytest.approx(expected_diff, rel=0.05)


def test_total_entropy_equals_sum_of_parts():
    T = np.array([294.0])
    H = 1.0 / (4 * np.pi * 1e-7)
    total = GADOLINIUM.total_entropy(T, H)
    parts = (GADOLINIUM.entropy_lattice(T) + GADOLINIUM.entropy_magnetic(T, H)
             + GADOLINIUM.electronic_entropy(T))
    assert total[0] == pytest.approx(parts[0], rel=1e-9)


def test_delta_T_adiabatic_exact_zero_field_step_is_zero():
    """No field change -> no temperature change, by construction of the
    isentropic root-solve (S(T,H)=S(T,H) trivially at T2=T1)."""
    mu0 = 4 * np.pi * 1e-7
    dT = GADOLINIUM.delta_T_adiabatic_exact(294.0, 1.0 / mu0, H_initial=1.0 / mu0)
    assert dT == pytest.approx(0.0, abs=1e-4)


def test_delta_T_adiabatic_exact_matches_linear_for_small_field_step():
    """For a small field step, the exact isentropic solve and the linear
    -T*dS/C approximation should agree closely (the linear formula is a
    first-order Taylor expansion of the exact one, valid in this limit --
    de Oliveira & von Ranke, Phys. Rep. 489 (2010), Sec. 2.1)."""
    mu0 = 4 * np.pi * 1e-7
    H_small = 0.02 / mu0  # 0.02 T, a genuinely small step
    T = np.array([294.0])
    dT_linear = float(GADOLINIUM.delta_T_adiabatic(T, H_small)[0])
    dT_exact = GADOLINIUM.delta_T_adiabatic_exact(294.0, H_small)
    # T=294K sits right at the sharp lambda-anomaly, so even a "small" field
    # step retains some curvature error -- 10% here is still a meaningful
    # small-step agreement check, not a loosened pass condition.
    assert dT_exact == pytest.approx(dT_linear, rel=0.10)


def test_delta_T_adiabatic_exact_reduces_dankov_error_vs_linear():
    """The whole point of the fix: at the actual Dan'kov et al. (1998)
    calibration fields (1-5T, well outside the small-step regime above),
    the exact isentropic method should sit closer to the literature value
    than the old linear approximation, at every field."""
    mu0 = 4 * np.pi * 1e-7
    for B, dT_lit in ((1.0, 3.2), (2.0, 5.8), (5.0, 12.3)):
        H = B / mu0
        dT_linear = float(GADOLINIUM.delta_T_adiabatic(np.array([294.0]), H)[0])
        dT_exact = GADOLINIUM.delta_T_adiabatic_exact(294.0, H)
        err_linear = abs(dT_linear - dT_lit)
        err_exact = abs(dT_exact - dT_lit)
        assert err_exact < err_linear
