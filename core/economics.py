"""
economics.py
============
Order-of-magnitude CAPEX/OPEX comparison.

Material costs are based on the magnetic refrigerator cost-optimization
study by Bjørk, Bahl & Smith, "Determining the minimum mass and cost of a
magnetic refrigerator", Int. J. Refrigeration 34 (2011) 1805-1816.

The model uses:
    - $40/kg for NdFeB (N42, 1.2-1.3 T remanence) permanent magnet material
    - $20/kg for gadolinium (Gd) magnetocaloric material

Their worked examples (100 W / 20 K device: 0.8 kg magnet + 0.3 kg Gd;
50 W / 30 K device: 0.15 kg magnet + 0.04 kg Gd) indicate that magnet mass
is typically several times the magnetocaloric material mass and increases
with magnetic field strength. This behaviour is approximated here as

    magnet_mass ≈ 3.0 × mu0H_max[T] × mass_regenerator

which provides a rough fit to the published examples but should not be
interpreted as a validated scaling law.

Sources:
    - Bjørk, Bahl & Smith, Int. J. Refrig. 34 (2011) 1805-1816 — magnet and
      magnetocaloric material costs and worked mass examples
    - Bjørk, Bahl & Nielsen, "The lifetime cost of a magnetic refrigerator",
      Int. J. Refrig. 63 (2016) 48-62 -- Phase 7 addition: this follow-up
      study by the same group adds device OPERATING cost (electricity over
      the device lifetime) to the same materials-only building cost. It
      explicitly assumes $0.10/kWh electricity (US/China/India-representative;
      many European countries are higher) and states that, like the 2011
      paper, "actual manufacturing, transportation, maintenance and
      auxiliary systems are ignored" -- so it does NOT resolve the
      HX/pump/motor/controls capital-cost gap either; that remains open
      (see `lifetime_cost()` docstring).
    - Bahl, Engelbrecht et al., Int. J. Refrig. 37 (2014) 78-83 — AMR
      system cost breakdown context
    - Lawrence Berkeley National Laboratory, "Data Center Cooling System
      Cost Benchmarks" — representative chilled-water OPEX
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
    specific multiplier is used here pending development of a detailed
    bottom-up bill-of-materials (BOM) model."""
    magnet_mass = MAGNET_TO_MCM_MASS_RATIO_PER_TESLA * mu0H_max * mass_regenerator
    return COST_MAGNET_PER_KG * magnet_mass + COST_MCM_PER_KG * mass_regenerator


ELECTRICITY_PRICE_PER_KWH = 0.10  # $/kWh, Bjørk, Bahl & Nielsen (2016), US/
                                    # China/India-representative; many
                                    # European countries are higher -- see
                                    # module docstring


def lifetime_cost(mu0H_max, mass_regenerator, Qc_avg_W, COP_electrical,
                   device_lifetime_years=15.0, capacity_factor=1.0,
                   electricity_price_per_kwh=ELECTRICITY_PRICE_PER_KWH):
    """Phase 7 addition, complementing (not replacing) `material_cost()`.

    Combines the materials-only building-cost floor with the device's
    lifetime OPERATING cost (electricity), following the methodology of
    Bjørk, Bahl & Nielsen (2016) -- see module docstring. `capacity_factor`
    lets the caller apply an average-load derating (analogous to
    `simple_tco()`'s `avg_load_fraction`) if the device doesn't run
    continuously at Qc_avg_W.

    THIS STILL DOES NOT CLOSE THE ROADMAP'S "full-system cost" GAP: the
    2016 source paper explicitly states that, like the 2011 paper, actual
    manufacturing, transportation, maintenance and auxiliary systems
    (i.e. the heat exchangers/pumps/motor-drive/controls/enclosure
    hardware itself, as opposed to the energy to run it) are ignored. No
    published $ breakdown for that hardware was found for this codebase
    (searched; nothing beyond materials-cost and now electricity-cost
    studies was located) -- so a full bottom-up BOM model for that
    remaining hardware stays a genuinely open Phase 7 item, not resolved
    by this function.

    Returns a dict with the materials floor, lifetime electricity cost,
    and their sum, so callers can see each piece separately rather than
    a single opaque total."""
    materials_floor = material_cost(mu0H_max, mass_regenerator)
    if COP_electrical <= 0:
        raise ValueError("COP_electrical must be positive")
    electrical_power_kW = (Qc_avg_W / COP_electrical) * capacity_factor / 1000.0
    hours_per_year = 24 * 365
    annual_electricity_cost = (electrical_power_kW * hours_per_year
                                * electricity_price_per_kwh)
    lifetime_electricity_cost = annual_electricity_cost * device_lifetime_years
    return {
        "materials_floor_$": round(materials_floor, 2),
        "lifetime_electricity_$": round(lifetime_electricity_cost, 2),
        "lifetime_total_$": round(materials_floor + lifetime_electricity_cost, 2),
        "device_lifetime_years": device_lifetime_years,
        "note": "Excludes HX/pump/motor/controls/enclosure hardware CAPEX -- "
                "see function docstring. Not a full-system cost.",
    }


@dataclass
class TCOResult:
    technology: str
    capex_per_kw_cooling: float   # $/kW_cooling installed
    opex_per_kwh_cooling: float    # $/kWh_cooling (electricity only)
    notes: str


AMR_MAGNETIC = TCOResult(
    "Magnetic (AMR)", capex_per_kw_cooling=2200.0, opex_per_kwh_cooling=0.012,
    notes="Pre-commercial; this $/kW figure is a rough placeholder and is not "
        "derived directly from material_cost(). Use "
        "material_cost(mu0H, mass_regenerator) with a specific design to "
        "estimate a materials-only cost floor. This excludes heat "
        "exchangers, pumps, motor/drive, controls and enclosure, which "
        "Bahl et al. (2014) identify as major contributors to total AMR "
        "system cost. A detailed bottom-up cost model is left for future work.")

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