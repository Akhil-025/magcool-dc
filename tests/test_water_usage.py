"""Phase 31 addition: core/water_usage.py had no dedicated test file even
though it is wired into main.py's pipeline (see README's Tier-1 test-
coverage gap). These tests mirror the style and level of test_emissions.py
and test_economics.py -- direct function calls, plain asserts, no fixtures."""

from core.water_usage import (
    annual_water_liters,
    compare_water_usage,
    WUE_L_PER_KWH_BY_REJECTION_CLASS,
)


def test_annual_water_liters_scales_linearly_with_capacity():
    small = annual_water_liters(50.0, cop_electrical=5.0, rejection_class="evaporative_tower")
    large = annual_water_liters(100.0, cop_electrical=5.0, rejection_class="evaporative_tower")
    assert large == small * 2


def test_annual_water_liters_dry_rejection_is_far_lower_than_evaporative():
    dry = annual_water_liters(100.0, cop_electrical=5.0, rejection_class="dry_air_cooled")
    evap = annual_water_liters(100.0, cop_electrical=5.0, rejection_class="evaporative_tower")
    assert dry < evap
    # WUE_L_PER_KWH_BY_REJECTION_CLASS says dry is ~36x lower (0.05 vs 1.8 L/kWh) --
    # check the ratio directly so this test breaks loudly if the benchmark
    # constants are ever edited without updating the module docstring.
    ratio = WUE_L_PER_KWH_BY_REJECTION_CLASS["evaporative_tower"] / WUE_L_PER_KWH_BY_REJECTION_CLASS["dry_air_cooled"]
    assert evap / dry == ratio


def test_annual_water_liters_is_independent_of_cop():
    # Documented explicitly in the module: WUE is a property of the
    # rejection technology (evaporative vs. dry), not of how much
    # electricity the cooling system itself consumes -- so COP must NOT
    # change the water-liters figure for a fixed rejection_class.
    low_cop = annual_water_liters(100.0, cop_electrical=3.0, rejection_class="evaporative_tower")
    high_cop = annual_water_liters(100.0, cop_electrical=15.0, rejection_class="evaporative_tower")
    assert low_cop == high_cop


def test_annual_water_liters_unknown_rejection_class_raises():
    try:
        annual_water_liters(100.0, cop_electrical=5.0, rejection_class="not_a_real_class")
        assert False, "expected a KeyError for an unrecognized rejection_class"
    except KeyError:
        pass


def test_compare_water_usage_returns_three_technologies_in_order():
    results = compare_water_usage(capacity_kW_IT=100.0, amr_cop=4.6, vcc_cop=3.2, liquid_cop=4.0)
    labels = [r.technology for r in results]
    assert labels == ["Magnetic (AMR)", "Vapor-compression", "Liquid cooling"]


def test_compare_water_usage_amr_default_rejection_class_is_dry_air_cooled():
    results = compare_water_usage(capacity_kW_IT=100.0, amr_cop=4.6, vcc_cop=3.2, liquid_cop=4.0)
    amr = results[0]
    assert amr.rejection_class == "dry_air_cooled"
    # AMR's own reported figure must be internally consistent with the
    # per-kWh benchmark table for whichever class it was assigned.
    assert amr.WUE_L_per_kWh_IT == WUE_L_PER_KWH_BY_REJECTION_CLASS["dry_air_cooled"]


def test_compare_water_usage_amr_uses_less_water_than_vcc_baseline_under_defaults():
    # This is the module's own headline claim (see write_water_usage_report's
    # "AMR ... uses X% less annual water" line) -- assert it holds under the
    # documented default rejection-class assignments, not just print it.
    results = compare_water_usage(capacity_kW_IT=100.0, amr_cop=4.6, vcc_cop=3.2, liquid_cop=4.0)
    amr = next(r for r in results if r.technology.startswith("Magnetic"))
    vcc = next(r for r in results if r.technology == "Vapor-compression")
    assert amr.annual_water_liters < vcc.annual_water_liters


def test_compare_water_usage_respects_caller_supplied_rejection_class_override():
    # The module docstring explicitly says a caller can override the
    # default assignment (e.g. an AMR paired with an evaporative tower by
    # design choice) -- check that override is actually honored, not
    # silently ignored.
    default_amr = compare_water_usage(100.0, amr_cop=4.6, vcc_cop=3.2, liquid_cop=4.0)[0]
    overridden_amr = compare_water_usage(
        100.0, amr_cop=4.6, vcc_cop=3.2, liquid_cop=4.0,
        amr_rejection_class="evaporative_tower")[0]
    assert overridden_amr.rejection_class == "evaporative_tower"
    assert overridden_amr.annual_water_liters > default_amr.annual_water_liters


def test_annual_water_liters_zero_capacity_is_zero():
    assert annual_water_liters(0.0, cop_electrical=5.0, rejection_class="dry_air_cooled") == 0.0
