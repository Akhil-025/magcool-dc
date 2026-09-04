import pytest

from core.mce_material import GADOLINIUM
from core.amr_cycle import (AMRSystem, BLOW_FRACTION_MASCHE,
                             _blow_fraction_multiplier, CYCLE_TYPE_FACTORS)
from core.loss_model import StateDependentLossModel
from core.thermal_diode import MechanicalContactDiode


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


# --- hysteresis loss wiring ----------------------------------

def test_hysteresis_power_is_zero_for_gadolinium():
    """GADOLINIUM (core.mce_material.MagnetocaloricMaterial) has no
    hysteresis_loss_J_per_kg attribute at all -- getattr's default must
    make _hysteresis_power_W() return exactly 0.0, so every previous
    Gd-based result in this test file (make_system()'s own default
    material) is completely unaffected by the addition."""
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
    fallback). deliberately adds hysteresis power in run()
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

# --- AMR cycle topology (cycle_type) ---

def test_cycle_type_defaults_to_brayton_and_is_backward_compatible():
    """Omitting cycle_type must give byte-for-byte the same result as
    explicitly passing cycle_type='brayton' -- the required backward-
    compatibility guarantee for this addition (same pattern as
    particle_diameter=None in and blow_fraction=0.5 earlier)."""
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
    """Per the qualitative ranking this addition targets (see
    CYCLE_TYPE_FACTORS's docstring): at fixed span/eps/field/frequency,
    Carnot-like should show the highest cooling capacity AND second-law
    (exergy) efficiency, Brayton-like (this model's previous default)
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


# ---------------------------------------------------------------------------
# thermal-diode-assisted AMRSystem (core/thermal_diode.py)
# ---------------------------------------------------------------------------

def test_thermal_diode_default_none_reproduces_pre_phase18_numbers():
    """thermal_diode=None (the default) must give IDENTICAL results to a
    system built without the parameter at all -- the same backward-
    compatibility guarantee gave particle_diameter/cycle_type."""
    sys_with_default = make_system()
    sys_explicit_none = make_system(thermal_diode=None)
    r1 = sys_with_default.run(291.0, 10.0)
    r2 = sys_explicit_none.run(291.0, 10.0)
    assert r1.W_parasitic == pytest.approx(r2.W_parasitic)
    assert r1.COP_electrical == pytest.approx(r2.COP_electrical)


def test_thermal_diode_adds_switching_power_to_parasitic():
    diode = MechanicalContactDiode(forward_conductance_W_K=5.0,
                                    reverse_conductance_W_K=0.5,
                                    actuation_energy_J_per_cycle=1.0)
    sys_no_diode = make_system()
    sys_diode = make_system(thermal_diode=diode)
    r_no_diode = sys_no_diode.run(291.0, 10.0)
    r_diode = sys_diode.run(291.0, 10.0)
    expected_switch_power = diode.switching_power_W(sys_diode.f)
    assert r_diode.W_parasitic == pytest.approx(
        r_no_diode.W_parasitic + expected_switch_power)


def test_thermal_diode_never_changes_Qc_or_W_mag():
    """ deliberately adds a cost-only term -- Qc and W_mag must be
    bit-for-bit unaffected by thermal_diode (see core/thermal_diode.py's
    honesty flag: no heat-transfer benefit is modeled)."""
    diode = MechanicalContactDiode(forward_conductance_W_K=5.0,
                                    reverse_conductance_W_K=0.5,
                                    actuation_energy_J_per_cycle=1.0)
    sys_no_diode = make_system()
    sys_diode = make_system(thermal_diode=diode)
    r_no_diode = sys_no_diode.run(291.0, 10.0)
    r_diode = sys_diode.run(291.0, 10.0)
    assert r_diode.Qc == pytest.approx(r_no_diode.Qc)
    assert r_diode.W_mag == pytest.approx(r_no_diode.W_mag)


def test_thermal_diode_zero_actuation_energy_is_a_noop():
    diode = MechanicalContactDiode(forward_conductance_W_K=5.0,
                                    reverse_conductance_W_K=0.5,
                                    actuation_energy_J_per_cycle=0.0)
    sys_no_diode = make_system()
    sys_diode = make_system(thermal_diode=diode)
    r_no_diode = sys_no_diode.run(291.0, 10.0)
    r_diode = sys_diode.run(291.0, 10.0)
    assert r_diode.W_parasitic == pytest.approx(r_no_diode.W_parasitic)


def test_thermal_diode_stacks_additively_with_hysteresis():
    """ (hysteresis) and (diode switching) must combine
    additively in W_parasitic, each independent of the other -- both are
    added unconditionally in run() via the same accounting pattern."""
    from core.first_order_mce import LAFESIH_FIRST_ORDER
    diode = MechanicalContactDiode(forward_conductance_W_K=5.0,
                                    reverse_conductance_W_K=0.5,
                                    actuation_energy_J_per_cycle=1.0)
    sys_neither = AMRSystem(material=GADOLINIUM, mu0H_max=2.0,
                             mass_regenerator=5.0, frequency=1.0,
                             fluid_mdot=0.1, regenerator_effectiveness=0.8)
    sys_hyst_only = AMRSystem(material=LAFESIH_FIRST_ORDER, mu0H_max=2.0,
                               mass_regenerator=5.0, frequency=1.0,
                               fluid_mdot=0.1, regenerator_effectiveness=0.8)
    sys_diode_only = AMRSystem(material=GADOLINIUM, mu0H_max=2.0,
                                mass_regenerator=5.0, frequency=1.0,
                                fluid_mdot=0.1, regenerator_effectiveness=0.8,
                                thermal_diode=diode)
    sys_both = AMRSystem(material=LAFESIH_FIRST_ORDER, mu0H_max=2.0,
                          mass_regenerator=5.0, frequency=1.0,
                          fluid_mdot=0.1, regenerator_effectiveness=0.8,
                          thermal_diode=diode)
    base = sys_neither.run(291.0, 10.0).W_parasitic
    hyst_extra = sys_hyst_only.run(291.0, 10.0).W_parasitic - base
    diode_extra = sys_diode_only.run(291.0, 10.0).W_parasitic - base
    both = sys_both.run(291.0, 10.0).W_parasitic
    assert both == pytest.approx(base + hyst_extra + diode_extra)

def test_cooling_capacity_span_sweep_is_monotonically_nonincreasing():
    """The whole point of the clamp: Qc_W must never increase as span
    increases, across the requested spans, even where raw Qc (unclamped)
    is not monotonic."""
    sys_ = make_system()
    spans = [2.0, 5.0, 8.0, 11.0, 14.0, 17.0]
    rows = sys_.cooling_capacity_span_sweep(291.0, spans)
    Qc_values = [r["Qc_W"] for r in rows]
    assert all(Qc_values[i] >= Qc_values[i + 1] - 1e-9 for i in range(len(Qc_values) - 1))


def test_cooling_capacity_span_sweep_clamp_never_exceeds_raw():
    sys_ = make_system()
    rows = sys_.cooling_capacity_span_sweep(291.0, [3.0, 6.0, 9.0, 12.0])
    for r in rows:
        assert r["Qc_W"] <= r["Qc_raw_W"] + 1e-9


def test_cooling_capacity_span_sweep_matches_raw_call_when_already_monotonic():
    """At small spans, well below any near-Tc discontinuity, raw
    cooling_capacity() should already be monotonic, so the clamp should
    be a no-op (Qc_W == Qc_raw_W) there."""
    sys_ = make_system()
    rows = sys_.cooling_capacity_span_sweep(291.0, [1.0, 2.0, 3.0])
    for r in rows:
        assert r["Qc_W"] == pytest.approx(r["Qc_raw_W"], rel=1e-6)


def test_cooling_capacity_span_sweep_dense_grid_catches_gap_between_sparse_spans():
    """If the reopening excursion sits entirely between two widely-spaced
    requested spans, the clamp must still catch it via the internal dense
    grid -- a running-min over only the sparse points would miss it."""
    sys_ = make_system()
    # Only two requested spans, straddling a wide gap -- if the internal
    # dense grid weren't used, a naive running-min over just these two
    # points could under-clamp relative to the same case scanned densely.
    sparse_rows = sys_.cooling_capacity_span_sweep(291.0, [2.0, 18.0])
    dense_rows = sys_.cooling_capacity_span_sweep(291.0, list(range(2, 19)))
    sparse_at_18 = next(r for r in sparse_rows if r["span_K"] == 18.0)
    dense_at_18 = next(r for r in dense_rows if r["span_K"] == 18.0)
    # The densely-scanned clamp can only be <= the sparse one (more points
    # seen along the way can only lower the running minimum further).
    assert dense_at_18["Qc_W"] <= sparse_at_18["Qc_W"] + 1e-9


def test_cooling_capacity_span_sweep_reproduces_tusek_reopening_clamp():
    """Direct reproduction of the exact Tusek AMR(A) V*=0.95 case
    documented in LIMITATIONS.md and diagnosed by
    core/validation_system.py's diagnose_qc_feasibility_reopening(): raw
    Qc at span=12.23K reopens to a large positive value (physically
    backwards, since Qc must not increase relative to smaller spans
    already evaluated); cooling_capacity_span_sweep() must clamp it to
    0W, matching what the model already reports at the neighboring
    spans in this same reopened window."""
    from core.validation_system import _calibrate_mdot, _t_cold_for_row
    row = {"device_group": "Tusek_fig10_AMRA_Vstar0.95",
           "span_K": "7.26", "Qc_W": "5.27",
           "mu0H_T": "1.15", "mass_MCM_kg": "0.1763", "frequency_Hz": "0.3",
           "material": "Gd (packed bed - AMR A, 0.1mm parallel plates)"}
    mdot_cal, sys_ = _calibrate_mdot(row)
    t_cold = _t_cold_for_row(row)

    Qc_raw, _ = sys_.cooling_capacity(t_cold, 12.23)
    assert Qc_raw > 5.0  # confirms the raw reopening artifact is still present

    rows = sys_.cooling_capacity_span_sweep(t_cold, [7.26, 12.23, 14.75])
    row_1223 = next(r for r in rows if r["span_K"] == 12.23)
    assert row_1223["Qc_raw_W"] == pytest.approx(Qc_raw, rel=1e-9)
    assert row_1223["Qc_W"] == pytest.approx(0.0, abs=1e-6)


def test_cooling_capacity_span_sweep_does_not_change_plain_cooling_capacity():
    """Additive-only guarantee: cooling_capacity() itself must be
    byte-for-byte unaffected by cooling_capacity_span_sweep() existing or
    having been called."""
    sys_ = make_system()
    before = sys_.cooling_capacity(291.0, 10.0)
    sys_.cooling_capacity_span_sweep(291.0, [2.0, 6.0, 10.0, 14.0])
    after = sys_.cooling_capacity(291.0, 10.0)
    assert before == after
