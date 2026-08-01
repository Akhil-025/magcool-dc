from core.emissions import refrigerant_emissions_tCO2e, operational_emissions_tCO2e, compare_emissions


def test_refrigerant_emissions_zero_leak_rate_is_zero():
    assert refrigerant_emissions_tCO2e(100.0, leak_rate=0.0) == 0.0


def test_operational_emissions_scale_inversely_with_cop():
    low_cop = operational_emissions_tCO2e(100.0, cop=3.0)
    high_cop = operational_emissions_tCO2e(100.0, cop=10.0)
    assert low_cop > high_cop


def test_compare_emissions_amr_has_zero_refrigerant_component():
    results = compare_emissions(100.0, amr_cop=5.0, vcc_cop=12.0, liquid_cop=20.0)
    amr = next(r for r in results if r.technology.startswith("Magnetic"))
    assert amr.refrigerant_GWP_tCO2e_per_year == 0.0
