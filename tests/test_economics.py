import pytest

from core.economics import material_cost, lifetime_cost


def test_material_cost_scales_with_field_and_mass():
    base = material_cost(mu0H_max=1.0, mass_regenerator=1.0)
    higher_field = material_cost(mu0H_max=2.0, mass_regenerator=1.0)
    more_mass = material_cost(mu0H_max=1.0, mass_regenerator=2.0)
    assert higher_field > base
    assert more_mass > base


def test_lifetime_cost_includes_materials_floor_and_electricity():
    result = lifetime_cost(mu0H_max=1.0, mass_regenerator=1.0, Qc_avg_W=500.0,
                            COP_electrical=3.0, device_lifetime_years=15.0)
    mat_floor = material_cost(mu0H_max=1.0, mass_regenerator=1.0)
    assert result["materials_floor_$"] == pytest.approx(mat_floor, rel=1e-6)
    assert result["lifetime_electricity_$"] > 0
    assert result["lifetime_total_$"] == pytest.approx(
        result["materials_floor_$"] + result["lifetime_electricity_$"], rel=1e-6)


def test_lifetime_cost_scales_with_lifetime_and_inversely_with_cop():
    short = lifetime_cost(mu0H_max=1.0, mass_regenerator=1.0, Qc_avg_W=500.0,
                           COP_electrical=3.0, device_lifetime_years=5.0)
    long = lifetime_cost(mu0H_max=1.0, mass_regenerator=1.0, Qc_avg_W=500.0,
                          COP_electrical=3.0, device_lifetime_years=15.0)
    assert long["lifetime_electricity_$"] > short["lifetime_electricity_$"]

    low_cop = lifetime_cost(mu0H_max=1.0, mass_regenerator=1.0, Qc_avg_W=500.0,
                             COP_electrical=2.0)
    high_cop = lifetime_cost(mu0H_max=1.0, mass_regenerator=1.0, Qc_avg_W=500.0,
                              COP_electrical=8.0)
    assert low_cop["lifetime_electricity_$"] > high_cop["lifetime_electricity_$"]


def test_lifetime_cost_rejects_nonpositive_cop():
    with pytest.raises(ValueError):
        lifetime_cost(mu0H_max=1.0, mass_regenerator=1.0, Qc_avg_W=500.0,
                       COP_electrical=0.0)


# =============================================================================
# Phase 19: geometric (Halbach-cylinder) magnet-mass cost term
# =============================================================================
from core.economics import (
    bom_cost,
    bom_cost_geometric,
    geometric_magnet_mass_kg,
    full_system_cost_estimate_geometric,
)


def test_geometric_magnet_mass_increases_with_field():
    low = geometric_magnet_mass_kg(mu0H_max=1.0, mass_regenerator=5.0)
    high = geometric_magnet_mass_kg(mu0H_max=2.0, mass_regenerator=5.0)
    assert high > low


def test_geometric_magnet_mass_increases_with_regenerator_mass():
    small = geometric_magnet_mass_kg(mu0H_max=1.5, mass_regenerator=2.0)
    large = geometric_magnet_mass_kg(mu0H_max=1.5, mass_regenerator=8.0)
    assert large > small


def test_bom_cost_geometric_return_shape_matches_bom_cost():
    """bom_cost_geometric() must return the same key set as bom_cost() so
    downstream callers (e.g. optimize.py's cost_index()) can treat them
    interchangeably via the same ["materials_bom_total_$"] lookup."""
    flat = bom_cost(1.5, 5.0, "Gd")
    geom = bom_cost_geometric(1.5, 5.0, "Gd")
    assert set(flat.keys()) == set(geom.keys())
    assert geom["materials_bom_total_$"] > 0


def test_bom_cost_geometric_diverges_from_flat_ratio_at_high_field():
    """At a high enough field the super-linear geometric relation must
    produce a materially different (larger, at fields well above 2T)
    total BOM cost than the flat per-Tesla ratio -- otherwise Phase 19
    would not have changed anything."""
    flat_low = bom_cost(1.0, 5.0, "Gd")["materials_bom_total_$"]
    geom_low = bom_cost_geometric(1.0, 5.0, "Gd")["materials_bom_total_$"]
    flat_high = bom_cost(3.0, 5.0, "Gd")["materials_bom_total_$"]
    geom_high = bom_cost_geometric(3.0, 5.0, "Gd")["materials_bom_total_$"]
    # ratio of geometric to flat cost should grow with field (super-linear
    # vs. linear), not stay constant
    ratio_low = geom_low / flat_low
    ratio_high = geom_high / flat_high
    assert ratio_high > ratio_low


def test_full_system_cost_estimate_geometric_applies_multiplier():
    from core.economics import NON_MATERIALS_COST_MULTIPLIER
    result = full_system_cost_estimate_geometric(1.5, 5.0, "Gd")
    expected = result["materials_bom_total_$"] * NON_MATERIALS_COST_MULTIPLIER
    assert result["full_system_cost_estimate_$"] == pytest.approx(expected, rel=1e-6)
    assert result["non_materials_multiplier"] == NON_MATERIALS_COST_MULTIPLIER


def test_amorphous_material_cost_performance_note_is_qualitative_only():
    """Phase 22 item 3: the amorphous-materials note must be a non-empty
    qualitative string, and must NOT silently be wired into
    MCM_COST_PER_KG_BY_FAMILY (no numeric $/kg is sourced for this repo's
    corpus -- see the note's own text and core/economics.py's section
    docstring for why)."""
    from core.economics import (
        amorphous_material_cost_performance_note, MCM_COST_PER_KG_BY_FAMILY,
    )
    note = amorphous_material_cost_performance_note()
    assert isinstance(note, str) and len(note) > 100
    assert "amorphous" in note.lower()
    assert "amorphous" not in " ".join(MCM_COST_PER_KG_BY_FAMILY.keys()).lower()


def test_amorphous_material_cost_performance_note_is_stable_across_calls():
    from core.economics import amorphous_material_cost_performance_note
    assert amorphous_material_cost_performance_note() == amorphous_material_cost_performance_note()

# =============================================================================
# Phase 31: non-materials sensitivity band, Rowe (2011) VCC cross-check,
# and the GE&R/CEC (2024) commercial MCM price reality check
# =============================================================================
from core.economics import (
    full_system_cost_estimate_range,
    NON_MATERIALS_COST_MULTIPLIER_LOW,
    NON_MATERIALS_COST_MULTIPLIER_MID,
    NON_MATERIALS_COST_MULTIPLIER_HIGH,
    rowe2011_vcc_compressor_cost_cross_check,
    commercial_mcm_price_reality_check,
    VAPOR_COMPRESSION,
)


def test_non_materials_multiplier_band_is_ordered():
    assert NON_MATERIALS_COST_MULTIPLIER_LOW < NON_MATERIALS_COST_MULTIPLIER_MID
    assert NON_MATERIALS_COST_MULTIPLIER_MID < NON_MATERIALS_COST_MULTIPLIER_HIGH


def test_full_system_cost_estimate_range_ordered_and_consistent_with_bom():
    result = full_system_cost_estimate_range(1.5, 5.0, "Gd")
    assert (result["full_system_cost_low_$"]
            < result["full_system_cost_mid_$"]
            < result["full_system_cost_high_$"])
    assert result["full_system_cost_mid_$"] == pytest.approx(
        result["materials_bom_total_$"] * NON_MATERIALS_COST_MULTIPLIER_MID, rel=1e-6)


def test_full_system_cost_estimate_range_mid_matches_point_estimate():
    """The MID band value must equal the existing full_system_cost_estimate()'s
    default point estimate exactly (same multiplier, same bom_cost() call) --
    this is an additive function, not a change in existing behavior."""
    from core.economics import full_system_cost_estimate
    point = full_system_cost_estimate(1.5, 5.0, "Gd")
    ranged = full_system_cost_estimate_range(1.5, 5.0, "Gd")
    assert ranged["full_system_cost_mid_$"] == pytest.approx(
        point["full_system_cost_estimate_$"], rel=1e-6)


def test_rowe2011_vcc_cross_check_reports_both_sources_unchanged():
    result = rowe2011_vcc_compressor_cost_cross_check()
    assert result["ashrae_derived_capex_usd_per_kw"] == VAPOR_COMPRESSION.capex_per_kw_cooling
    low, high = result["rowe2011_compressor_capex_usd_per_kw_range"]
    assert low == pytest.approx(500.0)
    assert high == pytest.approx(1900.0)
    assert low < high


def test_commercial_mcm_price_reality_check_does_not_mutate_working_price():
    from core.economics import MCM_COST_PER_KG_BY_FAMILY, COST_MCM_PER_KG
    before = dict(MCM_COST_PER_KG_BY_FAMILY)
    result = commercial_mcm_price_reality_check("Gd")
    assert MCM_COST_PER_KG_BY_FAMILY == before  # unchanged, see section docstring
    assert result["working_price_usd_per_kg"] == COST_MCM_PER_KG
    assert result["target_scaled_ratio_vs_working_price"] > 1.0


def test_commercial_mcm_price_reality_check_unrecognized_family_falls_back_to_gd():
    from core.economics import COST_MCM_PER_KG
    result = commercial_mcm_price_reality_check("NotARealFamily")
    assert result["working_price_usd_per_kg"] == COST_MCM_PER_KG


# =============================================================================
# Phase 32: bottom-up, market-catalog-sourced non-materials BOM
# =============================================================================
from core.economics import (
    bottom_up_non_materials_bom,
    full_system_cost_estimate_bottom_up,
    cross_check_full_system_cost_methods,
    exact_hx_duty_multiplier,
    HX_COST_PER_KW_RANGE,
    PUMP_COST_PER_KW_RANGE,
    MOTOR_COST_PER_KW_RANGE,
    DRIVE_COST_PER_KW_RANGE,
)


def test_exact_hx_duty_multiplier_matches_energy_balance():
    # multiplier = 2 + 1/COP by construction; check directly and at a
    # couple of representative COP values
    assert exact_hx_duty_multiplier(5.26) == pytest.approx(2.0 + 1.0 / 5.26)
    assert exact_hx_duty_multiplier(2.0) == pytest.approx(2.5)
    assert exact_hx_duty_multiplier(1.0) == pytest.approx(3.0)


def test_exact_hx_duty_multiplier_rejects_nonpositive_cop():
    with pytest.raises(ValueError):
        exact_hx_duty_multiplier(0.0)


def test_bottom_up_non_materials_bom_defaults_to_exact_multiplier():
    from core.economics import HX_COST_PER_KW_RANGE
    COP = 4.0
    result = bottom_up_non_materials_bom(Qc_avg_W=1000.0, COP_electrical=COP)
    expected_hx_kW = 1.0 * exact_hx_duty_multiplier(COP)
    expected_mid = HX_COST_PER_KW_RANGE[1] * expected_hx_kW
    assert result["heat_exchangers_$"][1] == pytest.approx(expected_mid, rel=1e-6)


def test_bottom_up_non_materials_bom_flat_override_still_works():
    result = bottom_up_non_materials_bom(Qc_avg_W=1000.0, COP_electrical=4.0,
                                           hx_duty_multiplier=2.0)
    assert result["heat_exchangers_$"][1] == pytest.approx(
        HX_COST_PER_KW_RANGE[1] * 1.0 * 2.0, rel=1e-6)


def test_component_cost_ranges_are_ordered_low_mid_high():
    for r in (HX_COST_PER_KW_RANGE, PUMP_COST_PER_KW_RANGE,
              MOTOR_COST_PER_KW_RANGE, DRIVE_COST_PER_KW_RANGE):
        assert r[0] < r[1] < r[2]


def test_bottom_up_non_materials_bom_scales_with_power():
    small = bottom_up_non_materials_bom(Qc_avg_W=500.0, COP_electrical=5.0)
    large = bottom_up_non_materials_bom(Qc_avg_W=5000.0, COP_electrical=5.0)
    assert (large["non_materials_bom_total_mid_$"]
            > small["non_materials_bom_total_mid_$"])
    # controls/enclosure allowance is fixed, not power-scaled
    assert small["controls_and_enclosure_$"] == large["controls_and_enclosure_$"]


def test_bottom_up_non_materials_bom_bands_are_ordered():
    result = bottom_up_non_materials_bom(Qc_avg_W=2000.0, COP_electrical=4.0)
    assert (result["non_materials_bom_total_low_$"]
            < result["non_materials_bom_total_mid_$"]
            < result["non_materials_bom_total_high_$"])


def test_bottom_up_non_materials_bom_rejects_nonpositive_cop():
    with pytest.raises(ValueError):
        bottom_up_non_materials_bom(Qc_avg_W=500.0, COP_electrical=0.0)


def test_full_system_cost_estimate_bottom_up_includes_materials_and_non_materials():
    result = full_system_cost_estimate_bottom_up(1.5, 5.0, Qc_avg_W=500.0,
                                                    COP_electrical=5.0)
    expected_mid = (result["materials_bom_total_$"]
                     + result["non_materials_breakdown"]["non_materials_bom_total_mid_$"])
    assert result["full_system_cost_bottom_up_mid_$"] == pytest.approx(expected_mid, rel=1e-6)
    assert (result["full_system_cost_bottom_up_low_$"]
            < result["full_system_cost_bottom_up_mid_$"]
            < result["full_system_cost_bottom_up_high_$"])


def test_cross_check_reports_both_methods_and_a_ratio():
    result = cross_check_full_system_cost_methods(1.5, 5.0, Qc_avg_W=500.0,
                                                     COP_electrical=5.0)
    assert result["borrowed_multiplier_method"]["mid_$"] > 0
    assert result["bottom_up_component_method"]["mid_$"] > 0
    assert result["borrowed_vs_bottom_up_mid_ratio"] > 0
    # documented finding: the borrowed-multiplier method reports a HIGHER
    # full-system cost than the bare-component bottom-up method at this
    # repo's own representative operating points (see function docstring)
    assert (result["borrowed_multiplier_method"]["mid_$"]
            > result["bottom_up_component_method"]["mid_$"])


# =============================================================================
# Phase 33: MAGNET_TO_MCM_MASS_RATIO_PER_TESLA cross-checked against 11
# real reported devices (Rowe 2011, Table 1)
# =============================================================================
from core.economics import (
    rowe2011_magnet_mass_ratio_cross_check,
    ROWE2011_DEVICE_MAGNET_MCM_DATA,
    MAGNET_TO_MCM_MASS_RATIO_PER_TESLA,
)


def test_rowe2011_cross_check_covers_all_eleven_devices():
    result = rowe2011_magnet_mass_ratio_cross_check()
    assert result["n_devices"] == len(ROWE2011_DEVICE_MAGNET_MCM_DATA) == 11
    assert len(result["per_device"]) == 11


def test_rowe2011_cross_check_ratios_are_positive_and_ordered():
    result = rowe2011_magnet_mass_ratio_cross_check()
    assert result["min_mass_ratio_per_tesla"] > 0
    assert result["min_mass_ratio_per_tesla"] <= result["median_mass_ratio_per_tesla"]
    assert result["median_mass_ratio_per_tesla"] <= result["max_mass_ratio_per_tesla"]


def test_rowe2011_cross_check_does_not_mutate_existing_constant():
    """See section docstring: this cross-check must never silently change
    the load-bearing MAGNET_TO_MCM_MASS_RATIO_PER_TESLA constant that
    material_cost()/bom_cost()/lifetime_cost() etc. all depend on."""
    from core import economics
    before = economics.MAGNET_TO_MCM_MASS_RATIO_PER_TESLA
    rowe2011_magnet_mass_ratio_cross_check()
    assert economics.MAGNET_TO_MCM_MASS_RATIO_PER_TESLA == before == MAGNET_TO_MCM_MASS_RATIO_PER_TESLA


def test_rowe2011_cross_check_reports_current_module_value():
    result = rowe2011_magnet_mass_ratio_cross_check()
    assert result["current_module_value"] == MAGNET_TO_MCM_MASS_RATIO_PER_TESLA


def test_rowe2011_per_device_mass_ratio_matches_manual_calculation():
    """Spot-check one device's arithmetic directly against the formula
    documented in the section docstring."""
    result = rowe2011_magnet_mass_ratio_cross_check()
    okamura = next(d for d in result["per_device"] if d["device"] == "Okamura")
    expected = (3.38 / 0.8) * (7.45 / 7.9) / 1.0
    assert okamura["mass_ratio_per_tesla"] == pytest.approx(expected, rel=1e-2)


# =============================================================================
# Phase 34: act on the Phase 33 finding -- update the working default AND
# keep the legacy value directly usable
# =============================================================================
from core.economics import (
    material_cost as _material_cost_p34,
    bom_cost as _bom_cost_p34,
    MAGNET_TO_MCM_MASS_RATIO_PER_TESLA_BJORK2011_LEGACY,
    MAGNET_TO_MCM_MASS_RATIO_PER_TESLA_ROWE2011_MEDIAN,
    compare_legacy_and_updated_magnet_ratio,
)


def test_working_default_equals_rowe2011_median_constant():
    assert MAGNET_TO_MCM_MASS_RATIO_PER_TESLA == MAGNET_TO_MCM_MASS_RATIO_PER_TESLA_ROWE2011_MEDIAN


def test_rowe2011_median_constant_matches_cross_check_function():
    """The hardcoded MAGNET_TO_MCM_MASS_RATIO_PER_TESLA_ROWE2011_MEDIAN
    constant must stay in sync with what rowe2011_magnet_mass_ratio_cross_check()
    actually computes from the raw device data -- if the device data ever
    changes, this test catches the drift."""
    result = rowe2011_magnet_mass_ratio_cross_check()
    assert MAGNET_TO_MCM_MASS_RATIO_PER_TESLA_ROWE2011_MEDIAN == pytest.approx(
        result["median_mass_ratio_per_tesla"], rel=1e-6)


def test_legacy_constant_unchanged_from_original_value():
    assert MAGNET_TO_MCM_MASS_RATIO_PER_TESLA_BJORK2011_LEGACY == 3.0


def test_material_cost_default_uses_updated_ratio():
    default = _material_cost_p34(1.5, 5.0)
    legacy = _material_cost_p34(1.5, 5.0, mass_ratio_per_tesla=MAGNET_TO_MCM_MASS_RATIO_PER_TESLA_BJORK2011_LEGACY)
    assert default > legacy  # updated ratio (13.47) > legacy ratio (3.0)


def test_material_cost_legacy_override_reproduces_pre_phase34_value():
    """Pinned pre-Phase-34 value: material_cost(1.5, 5.0) with the old flat
    3.0 ratio. If this ever changes, something broke reproducibility of
    this module's historical numbers."""
    legacy = _material_cost_p34(1.5, 5.0, mass_ratio_per_tesla=MAGNET_TO_MCM_MASS_RATIO_PER_TESLA_BJORK2011_LEGACY)
    assert legacy == pytest.approx(1000.0, rel=1e-6)  # 3.0*1.5*5=22.5kg*$40 + 5kg*$20


def test_bom_cost_accepts_mass_ratio_override():
    default = _bom_cost_p34(1.5, 5.0, "Gd")
    legacy = _bom_cost_p34(1.5, 5.0, "Gd", mass_ratio_per_tesla=MAGNET_TO_MCM_MASS_RATIO_PER_TESLA_BJORK2011_LEGACY)
    assert default["magnet_mass_kg"] > legacy["magnet_mass_kg"]


def test_compare_legacy_and_updated_magnet_ratio_reports_both():
    result = compare_legacy_and_updated_magnet_ratio(1.5, 5.0, "Gd")
    assert result["legacy_bjork2011"]["mass_ratio_per_tesla"] == MAGNET_TO_MCM_MASS_RATIO_PER_TESLA_BJORK2011_LEGACY
    assert result["updated_rowe2011_median"]["mass_ratio_per_tesla"] == MAGNET_TO_MCM_MASS_RATIO_PER_TESLA
    assert (result["updated_rowe2011_median"]["materials_bom_total_$"]
            > result["legacy_bjork2011"]["materials_bom_total_$"])
    assert result["materials_bom_total_ratio"] > 1.0


def test_compare_legacy_and_updated_magnet_ratio_matches_direct_bom_cost_calls():
    result = compare_legacy_and_updated_magnet_ratio(2.0, 8.0, "Gd")
    direct_legacy = _bom_cost_p34(2.0, 8.0, "Gd", mass_ratio_per_tesla=MAGNET_TO_MCM_MASS_RATIO_PER_TESLA_BJORK2011_LEGACY)
    direct_updated = _bom_cost_p34(2.0, 8.0, "Gd")
    assert result["legacy_bjork2011"]["materials_bom_total_$"] == direct_legacy["materials_bom_total_$"]
    assert result["updated_rowe2011_median"]["materials_bom_total_$"] == direct_updated["materials_bom_total_$"]


def test_rowe2011_cross_check_note_reflects_phase34_update():
    """The cross-check's own note should describe the update having
    happened, not merely flag an unresolved discrepancy (Phase 33's
    original framing, superseded by Phase 34)."""
    result = rowe2011_magnet_mass_ratio_cross_check()
    assert "legacy_value" in result
    assert result["legacy_value"] == MAGNET_TO_MCM_MASS_RATIO_PER_TESLA_BJORK2011_LEGACY
    assert "Phase 34" in result["note"]
