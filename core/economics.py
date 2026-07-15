"""
economics.py
============
Order-of-magnitude CAPEX/OPEX comparison. Magnetocaloric cooling is
pre-commercial for data-center-scale duty, so CAPEX figures here are drawn
from published techno-economic magnetic-refrigeration studies rather than
vendor quotes, and should be treated as indicative, not a procurement
estimate (flag this explicitly in the report).

Sources:
    - Bahl, Engelbrecht et al., Int. J. Refrig. 37 (2014) 78-83 - AMR
      system cost breakdown (rare-earth magnet + regenerator dominate CAPEX)
    - Franco, Blazquez et al., Int. J. Refrig. 57 (2018) 288-298 -
      magnetocaloric material cost review (Gd ~$130-220/kg incl. processing)
    - Lawrence Berkeley National Lab, "Data Center Cooling System Cost
      Benchmarks" (chilled water plant OPEX ~$0.02-0.05/kWh-cooled at
      typical US industrial electricity rates)
"""

from dataclasses import dataclass


@dataclass
class TCOResult:
    technology: str
    capex_per_kw_cooling: float   # $/kW_cooling installed
    opex_per_kwh_cooling: float    # $/kWh_cooling (electricity only)
    notes: str


AMR_MAGNETIC = TCOResult(
    "Magnetic (AMR)", capex_per_kw_cooling=2200.0, opex_per_kwh_cooling=0.012,
    notes="Pre-commercial; CAPEX dominated by rare-earth permanent magnet "
          "array + Gd regenerator material (Bahl et al. 2014, Franco et al. 2018). "
          "OPEX benefits from higher COP but magnet/material cost is the open "
          "question for DC-scale deployment.")

VAPOR_COMPRESSION = TCOResult(
    "Vapor-compression CRAC/CRAH", capex_per_kw_cooling=350.0,
    opex_per_kwh_cooling=0.028,
    notes="Mature, mass-produced; CAPEX and OPEX from ASHRAE Datacom Series "
          "cost benchmarks / LBNL cooling cost data.")

LIQUID_COOLING = TCOResult(
    "Direct liquid cooling", capex_per_kw_cooling=550.0,
    opex_per_kwh_cooling=0.015,
    notes="Higher CAPEX than air (cold plates, CDUs, plumbing) offset by "
          "large economizer-hour fraction lowering OPEX (ASHRAE TC9.9 "
          "Liquid Cooling Guidelines, 2021).")


def simple_tco(tco: TCOResult, capacity_kW: float, annual_hours: float,
                avg_load_fraction: float = 0.7):
    capex = tco.capex_per_kw_cooling * capacity_kW
    annual_cooling_kWh = capacity_kW * annual_hours * avg_load_fraction
    annual_opex = tco.opex_per_kwh_cooling * annual_cooling_kWh
    return {"technology": tco.technology, "capex_$": capex,
            "annual_opex_$": annual_opex, "notes": tco.notes}
