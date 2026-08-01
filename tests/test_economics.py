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
