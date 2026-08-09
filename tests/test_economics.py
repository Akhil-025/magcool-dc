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