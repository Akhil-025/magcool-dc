"""
economics.py
============
Order-of-magnitude CAPEX/OPEX comparison.

Phase 6 update: replaced the earlier placeholder $175/kg (loosely
attributed to "Franco et al. 2018 order-of-magnitude") with the actual
per-kg costs used in a dedicated magnetic-refrigerator cost-optimization
study, found via a targeted search this session: Bjørk, Bahl & Smith,
"Determining the minimum mass and cost of a magnetic refrigerator", Int. J.
Refrigeration 34 (2011) 1805-1816 -- $40/kg for NdFeB (N42, 1.2-1.3T
remanence) permanent magnet material, $20/kg for the magnetocaloric
material (Gd). Their own worked examples (100W/20K device: 0.8kg magnet +
0.3kg Gd; 50W/30K device: 0.15kg magnet + 0.04kg Gd) show magnet mass
running roughly 2.7-3.75x the MCM mass and increasing with field --
approximated here as magnet_mass ~= 3.0 * mu0H_max[T] * mass_regenerator,
a rough fit to their two reported points, NOT a validated scaling law.

Sources:
    - Bjork, Bahl & Smith, Int. J. Refrig. 34 (2011) 1805-1816 -- magnet/MCM
      unit costs and worked mass examples (used directly above)
    - Bahl, Engelbrecht et al., Int. J. Refrig. 37 (2014) 78-83 -- AMR
      system cost breakdown context
    - Lawrence Berkeley National Lab, "Data Center Cooling System Cost
      Benchmarks" (chilled water plant OPEX ~$0.02-0.05/kWh-cooled at
      typical US industrial electricity rates)
"""

from dataclasses import dataclass

COST_MCM_PER_KG = 20.0          # $/kg, Bjork et al. 2011
COST_MAGNET_PER_KG = 40.0        # $/kg, Bjork et al. 2011 (NdFeB N42)
MAGNET_TO_MCM_MASS_RATIO_PER_TESLA = 3.0  # rough fit to Bjork et al.'s two
                                             # worked examples, see docstring


def material_cost(mu0H_max, mass_regenerator):
    """Bottom-up magnet + MCM material cost, $ (Bjork et al. 2011 unit costs
    and mass-ratio approximation -- see module docstring). This is a
    materials-only FLOOR, not full system cost (excludes heat exchangers,
    pumps, motor/drive, controls, enclosure -- Bahl et al. 2014 note these
    dominate total AMR system cost, materials are a minority share, but no
    specific multiplier is cited here pending a bottom-up BOM in Phase 7)."""
    magnet_mass = MAGNET_TO_MCM_MASS_RATIO_PER_TESLA * mu0H_max * mass_regenerator
    return COST_MAGNET_PER_KG * magnet_mass + COST_MCM_PER_KG * mass_regenerator


@dataclass
class TCOResult:
    technology: str
    capex_per_kw_cooling: float   # $/kW_cooling installed
    opex_per_kwh_cooling: float    # $/kWh_cooling (electricity only)
    notes: str


AMR_MAGNETIC = TCOResult(
    "Magnetic (AMR)", capex_per_kw_cooling=2200.0, opex_per_kwh_cooling=0.012,
    notes="Pre-commercial; this $/kW figure is a rough placeholder, NOT yet "
          "derived from material_cost() above -- use material_cost(mu0H, "
          "mass_regenerator) directly with a specific optimize.py Pareto "
          "design for a grounded, materials-only cost floor (excludes heat "
          "exchangers/pumps/motor/controls, which Bahl et al. 2014 note "
          "dominate total system cost). Reconciling this placeholder with "
          "the bottom-up figure is a Phase 7 item.")

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
