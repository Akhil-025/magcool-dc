"""
emissions.py
=============
Phase 5: refrigerant-free GWP/emissions comparison. Magnetic (AMR) cooling
uses a single-phase heat transfer fluid (typically water/glycol) and no
HFC/HFO refrigerant at all -- a categorical difference from vapor-compression
and (indirectly, via its chiller-plant mechanical-cooling hours) liquid
cooling, both of which carry refrigerant charge with nonzero GWP and leak
risk. This module quantifies that avoided-emissions case, since the Phase
1-4 COP results do not, on their own, favor magnetic cooling.

Refrigerant GWP values (AR5 100-year GWP, IPCC):
    R-410A  : 2088   (common CRAC/DX refrigerant)
    R-134a  : 1430    (common chiller refrigerant)
    R-32     : 675      (lower-GWP DX alternative, increasingly used)
    R-1234ze : 7          (HFO, near-zero GWP, emerging in chiller plants)

Leak rate assumption: 2-10%/year of refrigerant charge is a commonly cited
range for commercial HVAC/refrigeration systems (US EPA GreenChill program
guidance; ASHRAE Standard 15 discussions of annual leakage). This module
uses a mid-range 4%/year default, clearly labeled as an assumption.

Refrigerant charge per kW cooling: ~0.3-0.5 kg/kW is representative for
packaged DX/chiller equipment (manufacturer datasheets, order-of-magnitude;
not a specific product spec). Default 0.4 kg/kW.
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
    print("\nNote: at these representative COPs, operational emissions dominate and "
          "AMR's lower COP (from Phase 1-4) makes its OPERATIONAL emissions the "
          "highest of the three -- the refrigerant-free case is a real, categorical "
          "benefit (zero leak risk, zero phase-out/replacement liability under "
          "F-gas-type regulation) but it does NOT overturn the emissions comparison "
          "on its own unless AMR's COP is closed relative to the baselines. Report "
          "both numbers plainly; don't let the qualitative refrigerant-free story "
          "imply a quantitative emissions win it hasn't earned yet.")
