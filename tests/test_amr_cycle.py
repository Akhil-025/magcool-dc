import pytest

from core.mce_material import GADOLINIUM
from core.amr_cycle import (AMRSystem, BLOW_FRACTION_MASCHE,
                             _blow_fraction_multiplier, CYCLE_TYPE_FACTORS)
from core.loss_model import StateDependentLossModel


def make_system(**overrides):
    kwargs = dict(material=GADOLINIUM, mu0H_max=2.0, mass_regenerator=5.0,
                  frequency=1.0, fluid_mdot=0.1, regenerator_effectiveness=0.8)
    kwargs.update(overrides)
    return AMRSystem(**kwargs)


def test_cooling_capacity_nonnegative():
    sys_ = make_system()
    Qc, dTad = sys_.cooling_capacity(T_cold=291.0, T_span=10.0)
    assert Qc >= 0.0


def test_cooling_capacity_zero_for_very_large_span():
    """Qc must fall to 0 once the imposed span swamps the no-load dTad."""
    sys_ = make_system()
    Qc_at_limit, _ = sys_.cooling_capacity(T_cold=291.0, T_span=200.0)
    assert Qc_at_limit == pytest.approx(0.0, abs=1e-6)


def test_electrical_cop_never_exceeds_ideal_cop():
    """Adding parasitic losses can only reduce COP, never increase it."""
    sys_ = make_system()
    result = sys_.run(T_cold=291.0, T_span=10.0)
    assert result.COP_electrical <= result.COP + 1e-9


def test_exergy_efficiency_between_zero_and_one():
    sys_ = make_system()
    result = sys_.run(T_cold=291.0, T_span=10.0)
    assert 0.0 <= result.exergy_eff <= 1.0 + 1e-9


def test_state_dependent_loss_model_increases_with_frequency_and_field():
    lm = StateDependentLossModel()
    low = lm.parasitic_power(frequency=0.5, mu0H=1.0, mdot=0.1, Qc=100.0)
    high = lm.parasitic_power(frequency=4.0, mu0H=3.0, mdot=0.1, Qc=100.0)
    assert high > low


def test_ntu_thermal_model_lets_mass_affect_cooling_capacity():
    """With use_ntu_thermal_model=True, more regenerator mass should improve
    (or at least not worsen) effectiveness / cooling capacity, unlike the
    constant-effectiveness model where mass has zero effect."""
    small = make_system(mass_regenerator=1.0, use_ntu_thermal_model=True)
    large = make_system(mass_regenerator=10.0, use_ntu_thermal_model=True)
    Qc_small, _ = small.cooling_capacity(291.0, 10.0)
    Qc_large, _ = large.cooling_capacity(291.0, 10.0)
    assert Qc_large >= Qc_small


def test_default_blow_fraction_is_symmetric_and_backward_compatible():
    """blow_fraction defaults to 0.5 (symmetric blow), which must return a
    multiplier of exactly 1.0 -- i.e. every pre-existing (pre-blow-fraction)
    result is reproduced exactly unless the caller opts in."""
    mult = _blow_fraction_multiplier(0.5, value_at_low=70.0, value_at_peak=330.0)
    assert mult == pytest.approx(1.0, rel=1e-9)


def test_blow_fraction_reproduces_masche_relative_qc_swing():
    """At the two Masche et al. (2022) reported blow fractions, Qc should
    scale by the same ~4.7x relative factor the paper reports between
    25.0% and 41.6% blow fraction (70W -> 330W at fixed T_span/U/f)."""
    sys_low = make_system(blow_fraction=BLOW_FRACTION_MASCHE["low"]["blow_fraction"])
    sys_peak = make_system(blow_fraction=BLOW_FRACTION_MASCHE["best_found"]["blow_fraction"])
    Qc_low, _ = sys_low.cooling_capacity(291.0, 10.0)
    Qc_peak, _ = sys_peak.cooling_capacity(291.0, 10.0)
    expected_ratio = (BLOW_FRACTION_MASCHE["best_found"]["Qc_W"]
                       / BLOW_FRACTION_MASCHE["low"]["Qc_W"])
    assert Qc_peak / Qc_low == pytest.approx(expected_ratio, rel=1e-6)


def test_blow_fraction_at_best_found_beats_symmetric_default():
    """The paper's best-found blow fraction (41.6%) should outperform the
    model's pre-existing implicit symmetric (50%) assumption, matching the
    paper's own qualitative finding that flow-profile asymmetry is a real,
    exploitable lever."""
    sys_peak = make_system(blow_fraction=BLOW_FRACTION_MASCHE["best_found"]["blow_fraction"])
    sys_symmetric = make_system(blow_fraction=0.5)
    r_peak = sys_peak.run(291.0, 10.0)
    r_symmetric = sys_symmetric.run(291.0, 10.0)
    assert r_peak.Qc > r_symmetric.Qc
    assert r_peak.COP_electrical > r_symmetric.COP_electrical


def test_amr_run_is_reasonably_fast():
    """Regression guard for the Newton-solver performance fix: 100
    sequential AMR evaluations (used heavily by optimize.py / sensitivity.py
    / rsm.py) should take well under a second, not the ~6s the old damped
    fixed-point solver needed."""
    import time
    lm = StateDependentLossModel()
    t0 = time.time()
    for _ in range(100):
        sys_ = make_system(loss_model=lm, use_ntu_thermal_model=True)
        sys_.run(291.0, 10.0)
    elapsed = time.time() - t0
    assert elapsed < 2.0, f"AMR evaluations took {elapsed:.2f}s for 100 calls (expected <2s)"


# --- Phase 16: hysteresis loss wiring ----------------------------------

def test_hysteresis_power_is_zero_for_gadolinium():
    """GADOLINIUM (core.mce_material.MagnetocaloricMaterial) has no
    hysteresis_loss_J_per_kg attribute at all -- getattr's default must
    make _hysteresis_power_W() return exactly 0.0, so every pre-Phase-16
    Gd-based result in this test file (make_system()'s own default
    material) is completely unaffected by the Phase 16 addition."""
    sys_ = make_system(mass_regenerator=8.0, frequency=2.0)
    assert sys_._hysteresis_power_W() == 0.0


def test_hysteresis_power_matches_formula_for_first_order_material():
    """W_hys = hysteresis_loss_J_per_kg * mass_regenerator * frequency,
    exactly -- direct formula check, independent of the rest of the AMR
    cycle solve."""
    from core.first_order_mce import lafesih_composition_tuned_material
    mat = lafesih_composition_tuned_material(285.0)
    sys_ = make_system(material=mat, mass_regenerator=6.0, frequency=3.0)
    expected = mat.hysteresis_loss_J_per_kg * 6.0 * 3.0
    assert sys_._hysteresis_power_W() == pytest.approx(expected)


def test_hysteresis_power_scales_linearly_with_mass_and_frequency():
    from core.first_order_mce import lafesih_composition_tuned_material
    mat = lafesih_composition_tuned_material(285.0)
    base = make_system(material=mat, mass_regenerator=5.0, frequency=1.0)
    double_mass = make_system(material=mat, mass_regenerator=10.0, frequency=1.0)
    double_freq = make_system(material=mat, mass_regenerator=5.0, frequency=2.0)
    assert double_mass._hysteresis_power_W() == pytest.approx(2 * base._hysteresis_power_W())
    assert double_freq._hysteresis_power_W() == pytest.approx(2 * base._hysteresis_power_W())


def test_hysteresis_increases_w_parasitic_and_reduces_cop_electrical():
    """End-to-end wiring check via run(): switching a first-order
    material's hysteresis_loss_J_per_kg from its literature-placeholder
    value to 0.0, with every other parameter held fixed, must (a) leave
    Qc and W (magnetic work) unchanged, since hysteresis is accounted for
    purely as an ADDITIONAL parasitic electrical load, not folded into
    the ideal-cycle thermodynamics, and (b) strictly increase
    W_parasitic and strictly decrease COP_electrical when hysteresis is
    turned on."""
    from core.first_order_mce import lafesih_composition_tuned_material
    mat = lafesih_composition_tuned_material(285.3)
    lm = StateDependentLossModel()
    sys_ = AMRSystem(material=mat, mu0H_max=1.105, mass_regenerator=14.82,
                      frequency=1.197, fluid_mdot=0.4131,
                      regenerator_effectiveness=0.9, loss_model=lm,
                      blow_fraction=0.41, particle_diameter=1.6991e-3,
                      bed_cross_section_area=0.002)

    original_hyst = mat.hysteresis_loss_J_per_kg
    assert original_hyst > 0.0, "test fixture assumes a nonzero placeholder"
    try:
        r_on = sys_.run(291.0, 10.0)
        mat.hysteresis_loss_J_per_kg = 0.0
        r_off = sys_.run(291.0, 10.0)
    finally:
        mat.hysteresis_loss_J_per_kg = original_hyst

    assert r_on.Qc == pytest.approx(r_off.Qc)
    assert r_on.W_mag == pytest.approx(r_off.W_mag)
    assert r_on.W_parasitic > r_off.W_parasitic
    assert r_on.COP_electrical < r_off.COP_electrical
    expected_delta = original_hyst * 14.82 * 1.197
    assert (r_on.W_parasitic - r_off.W_parasitic) == pytest.approx(expected_delta, rel=1e-6)


def test_hysteresis_power_included_in_no_loss_model_path_too():
    """core.cascade.py's _single_stage() baseline helper builds an
    AMRSystem WITHOUT a loss_model (using the constant parasitic_fraction
    fallback). Phase 16 deliberately adds hysteresis power in run()
    UNCONDITIONALLY (outside the `if self.loss_model is not None` branch)
    so this path is not silently missed -- regression guard for that
    specific wiring choice."""
    from core.first_order_mce import lafesih_composition_tuned_material
    mat = lafesih_composition_tuned_material(285.0)
    sys_ = AMRSystem(material=mat, mu0H_max=1.5, mass_regenerator=6.0,
                      frequency=2.0, fluid_mdot=0.3,
                      regenerator_effectiveness=0.85, loss_model=None,
                      parasitic_fraction=0.1)
    result = sys_.run(291.0, 10.0)
    expected_hyst = mat.hysteresis_loss_J_per_kg * 6.0 * 2.0
    expected_base = sys_.parasitic_fraction * result.Qc
    assert result.W_parasitic == pytest.approx(expected_base + expected_hyst)

# --- Phase 17: AMR cycle topology (cycle_type) ---

def test_cycle_type_defaults_to_brayton_and_is_backward_compatible():
    """Omitting cycle_type must give byte-for-byte the same result as
    explicitly passing cycle_type='brayton' -- the required backward-
    compatibility guarantee for this Phase 17 addition (same pattern as
    particle_diameter=None in Phase 15 and blow_fraction=0.5 earlier)."""
    default_sys = make_system()
    assert default_sys.cycle_type == "brayton"
    explicit_sys = make_system(cycle_type="brayton")
    r_default = default_sys.run(T_cold=291.0, T_span=10.0)
    r_explicit = explicit_sys.run(T_cold=291.0, T_span=10.0)
    assert r_default == r_explicit


def test_invalid_cycle_type_raises_value_error():
    with pytest.raises(ValueError):
        make_system(cycle_type="stirling")


def test_cycle_type_factors_brayton_is_identity():
    assert CYCLE_TYPE_FACTORS["brayton"] == {"qc_multiplier": 1.0, "eta_uplift": 1.0}


def test_cycle_type_ordering_carnot_ge_ericsson_ge_brayton():
    """Per the qualitative ranking this Phase 17 addition targets (see
    CYCLE_TYPE_FACTORS's docstring): at fixed span/eps/field/frequency,
    Carnot-like should show the highest cooling capacity AND second-law
    (exergy) efficiency, Brayton-like (this model's pre-Phase-17 default)
    the lowest, Ericsson-like in between."""
    results = {}
    for ct in ("brayton", "ericsson", "carnot"):
        sys_ = make_system(cycle_type=ct)
        results[ct] = sys_.run(T_cold=291.0, T_span=10.0)

    assert results["carnot"].Qc >= results["ericsson"].Qc >= results["brayton"].Qc
    assert (results["carnot"].exergy_eff >= results["ericsson"].exergy_eff
            >= results["brayton"].exergy_eff)
    # A strict inequality somewhere confirms the multipliers are actually
    # doing something, not just a no-op ordering by coincidence.
    assert results["carnot"].Qc > results["brayton"].Qc


def test_cycle_type_carnot_exergy_eff_still_bounded():
    """The uplifted eta_2nd_law must still respect the model's existing
    np.clip(..., 0.02, 0.95) ceiling -- Carnot-like is an idealized
    reference, not a claim of exceeding the model's own efficiency bound."""
    sys_ = make_system(cycle_type="carnot", regenerator_effectiveness=0.99)
    result = sys_.run(T_cold=291.0, T_span=5.0)
    assert 0.0 <= result.exergy_eff <= 1.0 + 1e-9


def test_cycle_type_qc_multiplier_applied_directly():
    """cooling_capacity() should scale linearly with qc_multiplier at a
    fixed operating point, independent of run()'s downstream loss
    accounting -- a direct regression guard on the wiring, not just the
    ordering above."""
    brayton_sys = make_system(cycle_type="brayton")
    ericsson_sys = make_system(cycle_type="ericsson")
    Qc_brayton, _ = brayton_sys.cooling_capacity(291.0, 10.0)
    Qc_ericsson, _ = ericsson_sys.cooling_capacity(291.0, 10.0)
    expected_ratio = (CYCLE_TYPE_FACTORS["ericsson"]["qc_multiplier"]
                       / CYCLE_TYPE_FACTORS["brayton"]["qc_multiplier"])
    assert Qc_ericsson == pytest.approx(Qc_brayton * expected_ratio, rel=1e-9)
