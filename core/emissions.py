"""
emissions.py
=============
Refrigerant-free GWP/emissions comparison.

Magnetic (AMR) cooling uses a single-phase heat transfer fluid (typically
water/glycol) and no HFC/HFO refrigerant. This distinguishes it from
vapor-compression and, where mechanical chilling is required, liquid
cooling systems, both of which use refrigerants with non-zero global
warming potential (GWP) and associated leakage risk.

Refrigerant GWP values (IPCC AR5, 100-year GWP):
    R-410A   : 2088   (common CRAC/DX refrigerant)
    R-134a   : 1430   (common chiller refrigerant)
    R-32     : 675    (lower-GWP DX alternative)
    R-1234ze : 7      (near-zero-GWP HFO)

Leak rate assumption: 2–10% per year is commonly cited for commercial
HVAC/refrigeration systems (US EPA GreenChill; ASHRAE Standard 15). This
module uses a clearly identified mid-range assumption of 4% per year.

Representative refrigerant charge: approximately 0.3–0.5 kg per kW of
cooling capacity for packaged DX/chiller equipment. A default value of
0.4 kg/kW is used here.
"""

from dataclasses import dataclass

REFRIGERANT_GWP = {
    "R-410A": 2088,
    "R-134a": 1430,
    "R-32": 675,
    "R-1234ze": 7,
}

CHARGE_KG_PER_KW = 0.4
DEFAULT_LEAK_RATE_PER_YEAR = 0.04
CO2_PER_KWH_GRID = 0.71   # kg CO2/kWh, representative Indian grid average
                            # (CEA CO2 Baseline Database order-of-magnitude;
                            # cite the latest CEA figure directly in the paper
                            # rather than relying on this placeholder)


@dataclass
class EmissionsResult:
    technology: str
    refrigerant_GWP_tCO2e_per_year: float
    operational_CO2_tCO2e_per_year: float
    total_tCO2e_per_year: float


def refrigerant_emissions_tCO2e(capacity_kW, refrigerant="R-410A",
                                  leak_rate=DEFAULT_LEAK_RATE_PER_YEAR):
    charge_kg = CHARGE_KG_PER_KW * capacity_kW
    leaked_kg_per_year = charge_kg * leak_rate
    gwp = REFRIGERANT_GWP[refrigerant]
    return leaked_kg_per_year * gwp / 1000.0  # tCO2e


def operational_emissions_tCO2e(capacity_kW, cop, annual_hours=8760,
                                  avg_load_fraction=0.7,
                                  co2_per_kwh=CO2_PER_KWH_GRID):
    cooling_kWh = capacity_kW * annual_hours * avg_load_fraction
    electricity_kWh = cooling_kWh / cop if cop > 0 else float("inf")
    return electricity_kWh * co2_per_kwh / 1000.0  # tCO2e


def compare_emissions(capacity_kW, amr_cop, vcc_cop, liquid_cop,
                        vcc_refrigerant="R-410A", liquid_refrigerant="R-134a"):
    results = []
    results.append(EmissionsResult(
        "Magnetic (AMR) - no refrigerant", 0.0,
        operational_emissions_tCO2e(capacity_kW, amr_cop),
        operational_emissions_tCO2e(capacity_kW, amr_cop)))
    vcc_refrig = refrigerant_emissions_tCO2e(capacity_kW, vcc_refrigerant)
    vcc_op = operational_emissions_tCO2e(capacity_kW, vcc_cop)
    results.append(EmissionsResult("Vapor-compression", vcc_refrig, vcc_op, vcc_refrig + vcc_op))
    liq_refrig = refrigerant_emissions_tCO2e(capacity_kW, liquid_refrigerant)
    liq_op = operational_emissions_tCO2e(capacity_kW, liquid_cop)
    results.append(EmissionsResult("Liquid cooling", liq_refrig, liq_op, liq_refrig + liq_op))
    return results


if __name__ == "__main__":
    print("Annual emissions comparison, 100 kW cooling capacity, illustrative COPs")
    print("(operational emissions dominate at any realistic COP -- refrigerant leak")
    print(" emissions are a secondary but nonzero, categorically-avoided term for AMR)")
    print("=" * 80)
    for r in compare_emissions(100.0, amr_cop=5.0, vcc_cop=12.0, liquid_cop=20.0):
        print(f"{r.technology:<32} refrigerant={r.refrigerant_GWP_tCO2e_per_year:7.2f} "
              f"tCO2e/yr  operational={r.operational_CO2_tCO2e_per_year:8.2f} tCO2e/yr  "
              f"total={r.total_tCO2e_per_year:8.2f} tCO2e/yr")
    print("\nNote: at these representative COPs, operational emissions dominate. "
        "If AMR has a lower COP than the baseline technologies, its operational "
        "emissions will also be higher. The refrigerant-free design remains a "
        "genuine environmental benefit (zero direct refrigerant leakage and no "
        "refrigerant phase-out liability), but it does not by itself guarantee "
        "lower total emissions. Both operational and refrigerant-related "
        "emissions should therefore be reported separately.")