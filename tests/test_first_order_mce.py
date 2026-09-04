"""
Unit tests for core/first_order_mce.py.

Previously this module was only exercised indirectly via main.py's
integration run (through giant_mce_analysis.py / cascade.py). These tests
target the Landau model itself: the calibrated GD5SI2GE2_FIRST_ORDER /
LAFESIH_FIRST_ORDER materials, the composition-tuning helpers and their
documented Tc bounds, and the Giguere dTad_correction wiring.
"""
import numpy as np
import pytest

from core.first_order_mce import (
    GD5SI2GE2_FIRST_ORDER, LAFESIH_FIRST_ORDER, MNFEPSI_FIRST_ORDER,
    MNCUCOGE_FIRST_ORDER,
    composition_tuned_material, lafesih_composition_tuned_material,
    mnfepsi_composition_tuned_material,
    GIANT_MCE_TC_MIN_K, GIANT_MCE_TC_MAX_K, LAFESIH_TC_MIN_K, LAFESIH_TC_MAX_K,
    MNFEPSI_TC_MIN_K, MNFEPSI_TC_MAX_K,
)


def _peak_dS_and_dT(material, T_range, mu0H_tesla):
    mu0 = 4 * np.pi * 1e-7
    Ts = np.linspace(*T_range, 401)
    H = mu0H_tesla / mu0
    dS = material.delta_S_isothermal(Ts, H)
    i_peak = int(np.argmin(dS))
    dT = material.delta_T_adiabatic(Ts, H)
    return float(dS[i_peak]), float(dT[i_peak])


def test_gd5si2ge2_peak_entropy_change_matches_calibration_target():
    """The module docstring's calibration target is peak |DeltaS_M| ~ 18
    J/(kg K) at 5T, found by scanning T near Tc rather than evaluating at
    the fixed nominal Tc (the field-shifted peak sits a few K above Tc)."""
    dS_peak, _ = _peak_dS_and_dT(GD5SI2GE2_FIRST_ORDER, (260.0, 300.0), 5.0)
    assert dS_peak == pytest.approx(-18.0, rel=0.15)


def test_gd5si2ge2_peak_entropy_change_is_negative_and_grows_with_field():
    """Magnetization increases entropy order, so DeltaS_M on magnetizing
    must be negative, and a larger applied field should produce a larger
    (more negative) peak entropy change."""
    dS_1T, _ = _peak_dS_and_dT(GD5SI2GE2_FIRST_ORDER, (260.0, 300.0), 1.0)
    dS_5T, _ = _peak_dS_and_dT(GD5SI2GE2_FIRST_ORDER, (260.0, 300.0), 5.0)
    assert dS_1T < 0 and dS_5T < 0
    assert dS_5T < dS_1T


def test_lafesih_peak_entropy_change_matches_calibration_target():
    """Calibration target for La(Fe0.90Si0.10)13H1.1: peak |DeltaS_M| ~ 31
    J/(kg K) at 5T near Tc=287K (Fujieda, Fujita & Fukamichi 2002)."""
    dS_peak, _ = _peak_dS_and_dT(LAFESIH_FIRST_ORDER, (272.0, 312.0), 5.0)
    assert dS_peak == pytest.approx(-31.0, rel=0.1)


def test_mnfepsi_peak_entropy_change_matches_calibration_target():
    """Calibration target for Mn0.68Fe1.22P0.62Si0.38: peak |DeltaS_M| ~ 17.6
    J/(kg K) at 2T near Tc=331.2K (Hanggai et al., Acta Materialia 302
    (2026) 121677) -- NOTE this target is at 2T, not the 5T used for
    GD5SI2GE2_FIRST_ORDER/LAFESIH_FIRST_ORDER."""
    dS_peak, _ = _peak_dS_and_dT(MNFEPSI_FIRST_ORDER, (316.0, 356.0), 2.0)
    assert dS_peak == pytest.approx(-17.6, rel=0.1)


def test_dTad_correction_defaults_to_uncorrected_for_calibrated_materials():
    """GD5SI2GE2_FIRST_ORDER, LAFESIH_FIRST_ORDER and MNFEPSI_FIRST_ORDER
    are the base calibrations documented in the module -- none should
    silently carry an applied Giguere correction; only
    composition_tuned_material() with apply_giguere_correction=True
    should."""
    assert GD5SI2GE2_FIRST_ORDER.dTad_correction == 1.0
    assert LAFESIH_FIRST_ORDER.dTad_correction == 1.0
    assert MNFEPSI_FIRST_ORDER.dTad_correction == 1.0


def test_composition_tuned_material_applies_giguere_correction_by_default():
    """apply_giguere_correction defaults to True and should scale
    delta_T_adiabatic relative to the uncorrected (correction=1.0) case,
    without touching delta_S_isothermal (per the field's own docstring:
    the correction is NOT applied to entropy, only to DeltaT_ad)."""
    from core.giguere_validation import DTAD_CORRECTION_FACTOR

    corrected = composition_tuned_material(276.0, apply_giguere_correction=True)
    uncorrected = composition_tuned_material(276.0, apply_giguere_correction=False)
    assert corrected.dTad_correction == pytest.approx(DTAD_CORRECTION_FACTOR)
    assert uncorrected.dTad_correction == 1.0

    T = np.array([278.0])
    H = 5.0 / (4 * np.pi * 1e-7)
    dTad_corrected = corrected.delta_T_adiabatic(T, H)[0]
    dTad_uncorrected = uncorrected.delta_T_adiabatic(T, H)[0]
    assert dTad_corrected == pytest.approx(
        dTad_uncorrected * DTAD_CORRECTION_FACTOR, rel=1e-9)

    dS_corrected = corrected.delta_S_isothermal(T, H)[0]
    dS_uncorrected = uncorrected.delta_S_isothermal(T, H)[0]
    assert dS_corrected == pytest.approx(dS_uncorrected, rel=1e-9)


def test_composition_tuned_material_rejects_out_of_range_tc():
    """Tc targets outside the documented giant-MCE window for the
    Gd5(SixGe1-x)4(-Ga) family must raise, not silently extrapolate --
    this is the explicit guard the module docstring describes, including
    the specific case of Gd5Si4 (Tc=335K), which is a normal second-order
    ferromagnet, not part of this giant-MCE family."""
    with pytest.raises(ValueError):
        composition_tuned_material(GIANT_MCE_TC_MIN_K - 1.0)
    with pytest.raises(ValueError):
        composition_tuned_material(GIANT_MCE_TC_MAX_K + 1.0)
    with pytest.raises(ValueError):
        composition_tuned_material(335.0)  # Gd5Si4, explicitly out of family


def test_composition_tuned_material_accepts_boundary_tc_values():
    lo = composition_tuned_material(GIANT_MCE_TC_MIN_K)
    hi = composition_tuned_material(GIANT_MCE_TC_MAX_K)
    assert lo.Tc == GIANT_MCE_TC_MIN_K
    assert hi.Tc == GIANT_MCE_TC_MAX_K


def test_lafesih_composition_tuned_material_rejects_out_of_range_tc():
    with pytest.raises(ValueError):
        lafesih_composition_tuned_material(LAFESIH_TC_MIN_K - 1.0)
    with pytest.raises(ValueError):
        lafesih_composition_tuned_material(LAFESIH_TC_MAX_K + 1.0)


def test_lafesih_composition_tuned_material_never_applies_giguere_correction():
    """No Giguere-style cross-check exists for this material family (per
    the docstring), so dTad_correction should always pass through
    LAFESIH_FIRST_ORDER's own (uncorrected, 1.0) default."""
    mat = lafesih_composition_tuned_material(300.0)
    assert mat.dTad_correction == LAFESIH_FIRST_ORDER.dTad_correction == 1.0


def test_mnfepsi_composition_tuned_material_rejects_out_of_range_tc():
    with pytest.raises(ValueError):
        mnfepsi_composition_tuned_material(MNFEPSI_TC_MIN_K - 1.0)
    with pytest.raises(ValueError):
        mnfepsi_composition_tuned_material(MNFEPSI_TC_MAX_K + 1.0)


def test_mnfepsi_composition_tuned_material_never_applies_giguere_correction():
    """No Giguere-style cross-check exists for this material family (per
    the docstring), so dTad_correction should always pass through
    MNFEPSI_FIRST_ORDER's own (uncorrected, 1.0) default."""
    mat = mnfepsi_composition_tuned_material(310.0)
    assert mat.dTad_correction == MNFEPSI_FIRST_ORDER.dTad_correction == 1.0


def test_mnfepsi_composition_tuned_material_shifts_tc_only():
    mat = mnfepsi_composition_tuned_material(300.0)
    assert mat.Tc == 300.0
    assert mat.A == MNFEPSI_FIRST_ORDER.A
    assert mat.B == MNFEPSI_FIRST_ORDER.B
    assert mat.C == MNFEPSI_FIRST_ORDER.C
    assert mat.theta_D == MNFEPSI_FIRST_ORDER.theta_D
    assert mat.M_molar == MNFEPSI_FIRST_ORDER.M_molar


def test_composition_tuned_material_shifts_tc_only():
    """Per the documented simplifying assumption, only Tc should change
    between a composition-tuned material and the base calibration -- A, B,
    C, theta_D, M_molar and n_atoms_per_fu are all held fixed."""
    mat = composition_tuned_material(250.0, apply_giguere_correction=False)
    assert mat.Tc == 250.0
    assert mat.A == GD5SI2GE2_FIRST_ORDER.A
    assert mat.B == GD5SI2GE2_FIRST_ORDER.B
    assert mat.C == GD5SI2GE2_FIRST_ORDER.C
    assert mat.theta_D == GD5SI2GE2_FIRST_ORDER.theta_D
    assert mat.M_molar == GD5SI2GE2_FIRST_ORDER.M_molar


def test_lattice_heat_capacity_positive_and_increases_with_temperature():
    """Debye lattice heat capacity should be positive everywhere and
    monotonically increasing well below the Debye temperature (standard
    T^3-ish low-T behavior)."""
    Ts = np.array([100.0, 200.0, 276.0])
    C = GD5SI2GE2_FIRST_ORDER.lattice_heat_capacity(Ts)
    assert np.all(C > 0)
    assert C[0] < C[1] < C[2]


def test_equilibrium_m_is_odd_in_field_at_fixed_tau():
    """At fixed reduced temperature, reversing the sign of the reduced
    field should flip the sign of the equilibrium magnetization (the
    Landau free energy here has no explicit field-odd bias term beyond
    -h*m)."""
    tau = 0.5  # deep in the ordered phase
    h = 0.05
    m_pos = GD5SI2GE2_FIRST_ORDER._equilibrium_m(tau, h)
    m_neg = GD5SI2GE2_FIRST_ORDER._equilibrium_m(tau, -h)
    assert m_pos == pytest.approx(-m_neg, abs=1e-6)


def test_delta_S_isothermal_accepts_array_input():
    T = np.linspace(260.0, 300.0, 11)
    dS = GD5SI2GE2_FIRST_ORDER.delta_S_isothermal(T, 5.0 / (4 * np.pi * 1e-7))
    assert dS.shape == T.shape
    assert np.all(np.isfinite(dS))


# --- hysteresis_loss_J_per_kg -------------------------------

def test_hysteresis_loss_default_is_zero_for_bare_dataclass():
    """A FirstOrderMCEMaterial instance created without specifying
    hysteresis_loss_J_per_kg must default to 0.0 (dataclass default) --
    this is what keeps every FirstOrderMCEMaterial instance predating
     (e.g. any built directly in a test or notebook without the
    new field) behaviorally identical to before ."""
    from core.first_order_mce import FirstOrderMCEMaterial
    mat = FirstOrderMCEMaterial(name="bare", Tc=280.0, A=10.0, B=-4.0, C=8.0,
                                 J=3.5, g=2.0, M_molar=0.157, theta_D=200.0,
                                 n_atoms_per_fu=9)
    assert mat.hysteresis_loss_J_per_kg == 0.0


def test_three_first_order_families_have_nonzero_hysteresis_loss():
    """The three calibrated first-order constants should each carry a
    positive (nonzero)  hysteresis-loss placeholder -- see each
    constant's own block comment in core/first_order_mce.py for the
    literature source and honesty flags on the exact value."""
    assert GD5SI2GE2_FIRST_ORDER.hysteresis_loss_J_per_kg > 0.0
    assert LAFESIH_FIRST_ORDER.hysteresis_loss_J_per_kg > 0.0
    assert MNFEPSI_FIRST_ORDER.hysteresis_loss_J_per_kg > 0.0


def test_composition_tuned_material_inherits_hysteresis_loss():
    """composition_tuned_material() (Gd5(SixGe1-x)4 family) must carry
    hysteresis_loss_J_per_kg through unchanged from the base
    GD5SI2GE2_FIRST_ORDER constant -- this was a real bug fixed during
     implementation (the field was silently dropping to its 0.0
    dataclass default for every tuned instance before the fix, which
    would have made every material core.optimize.py/core.cascade.py
    actually use hysteresis-free regardless of the base constant's own
    value)."""
    tc = 285.0
    mat = composition_tuned_material(tc)
    assert mat.hysteresis_loss_J_per_kg == GD5SI2GE2_FIRST_ORDER.hysteresis_loss_J_per_kg


def test_lafesih_composition_tuned_material_inherits_hysteresis_loss():
    tc = 285.0
    mat = lafesih_composition_tuned_material(tc)
    assert mat.hysteresis_loss_J_per_kg == LAFESIH_FIRST_ORDER.hysteresis_loss_J_per_kg


def test_mnfepsi_composition_tuned_material_inherits_hysteresis_loss():
    tc = 310.0
    mat = mnfepsi_composition_tuned_material(tc)
    assert mat.hysteresis_loss_J_per_kg == MNFEPSI_FIRST_ORDER.hysteresis_loss_J_per_kg


def test_hysteresis_loss_is_mutable_field_not_frozen():
    """core.hysteresis_sensitivity.py's ON/OFF A/B comparison depends on
    being able to mutate the module-level *_FIRST_ORDER constants'
    hysteresis_loss_J_per_kg in place and restore it afterward -- this
    requires FirstOrderMCEMaterial to NOT be a frozen dataclass. Guards
    against a future refactor silently freezing the class and breaking
    that diagnostic."""
    original = LAFESIH_FIRST_ORDER.hysteresis_loss_J_per_kg
    try:
        LAFESIH_FIRST_ORDER.hysteresis_loss_J_per_kg = 999.0
        assert LAFESIH_FIRST_ORDER.hysteresis_loss_J_per_kg == 999.0
    finally:
        LAFESIH_FIRST_ORDER.hysteresis_loss_J_per_kg = original
    assert LAFESIH_FIRST_ORDER.hysteresis_loss_J_per_kg == original


# ---------------------------------------------------------------------------
# _equilibrium_m()'s h==0 closed-form fast path (this session's speedup for
# main.py's step 11f / core.cascade.run_graded_cascade()): guards against a
# silent regression of the ONE thing that matters -- m**2 (the only quantity
# any caller ever uses, via delta_S_isothermal()'s s = -0.5*A*m**2) must
# match a from-scratch, independent np.roots()-based reference computation
# to numerical precision, for every h==0 call the fast path intercepts. The
# raw sign of m is NOT checked (see _equilibrium_m()'s own docstring: it is
# a genuine, physically meaningless degeneracy at h=0 that the pre-fast-path
# implementation also had no defined convention for).
# ---------------------------------------------------------------------------

def _reference_equilibrium_m(mat, tau, h):
    """Independent re-implementation of the ORIGINAL (pre-fast-path)
    np.roots()-based algorithm, kept local to this test rather than
    imported, so this test can't be satisfied by a bug that breaks both
    the fast path and the "reference" identically."""
    coeffs = [mat.C, 0, mat.B, 0, mat.A * (tau - 1), -h]
    roots = np.roots(coeffs)
    real_roots = roots[np.abs(roots.imag) < 1e-6].real
    real_roots = real_roots[np.abs(real_roots) <= 1.5]
    if len(real_roots) == 0:
        return 0.0

    def f(m):
        return (0.5 * mat.A * (tau - 1) * m ** 2 + 0.25 * mat.B * m ** 4
                + (mat.C / 6) * m ** 6 - h * m)
    vals = [f(m) for m in real_roots]
    return real_roots[int(np.argmin(vals))]


@pytest.mark.parametrize("material", [
    GD5SI2GE2_FIRST_ORDER, LAFESIH_FIRST_ORDER, MNFEPSI_FIRST_ORDER, MNCUCOGE_FIRST_ORDER,
])
def test_equilibrium_m_h_zero_fast_path_matches_reference_msquared(material):
    for tau in np.linspace(0.80, 1.20, 41):
        fast = material._equilibrium_m(tau, 0.0)
        ref = _reference_equilibrium_m(material, tau, 0.0)
        assert fast ** 2 == pytest.approx(ref ** 2, abs=1e-9), (
            f"{material.name} at tau={tau}: fast path m**2={fast**2!r} != "
            f"reference m**2={ref**2!r}")


def test_equilibrium_m_h_nonzero_path_unchanged():
    """h != 0 must still go through the original general quintic solve,
    completely untouched by the fast path -- checked by comparing directly
    against the same independent reference implementation."""
    mat = GD5SI2GE2_FIRST_ORDER
    for tau in np.linspace(0.85, 1.15, 15):
        for h in (0.001, 0.05, 0.2, -0.1):
            actual = mat._equilibrium_m(tau, h)
            ref = _reference_equilibrium_m(mat, tau, h)
            assert actual == pytest.approx(ref, abs=1e-9)


def test_delta_S_isothermal_end_to_end_finite_after_fast_path():
    """Full delta_S_isothermal() output (what every real caller actually
    uses) must stay finite and sensibly-signed (entropy change from
    field-on is <=0 for a normal MCE material) across each family's own
    transition region."""
    for mat in (GD5SI2GE2_FIRST_ORDER, LAFESIH_FIRST_ORDER, MNFEPSI_FIRST_ORDER):
        T = np.linspace(mat.Tc - 20, mat.Tc + 20, 50)
        dS = mat.delta_S_isothermal(T, H_final=4.0 / (4e-7 * np.pi), H_initial=0.0)
        assert np.all(np.isfinite(dS))
        assert np.all(dS <= 1e-9)