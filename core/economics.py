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
    - Gauss, Homm & Gutfleisch, "The Resource Basis of Magnetic
      Refrigeration", J. Ind. Ecol. 21(5) (2016) 1291-1300 -- phase 15
      addition: raw-material supply-criticality assessment for Gd-,
      La-, and Mn-based magnetocaloric alloys plus the Nd2Fe14B magnet
      material; see `resource_criticality_note()` below
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


# =============================================================================
# Phase 15 addition: full-system BOM cost model
# =============================================================================
#
# `material_cost()`/`lifetime_cost()` above are explicitly documented as a
# MATERIALS-ONLY floor (magnet + MCM), not a full-system cost -- Phase 14's
# item B6 re-confirmed this gap was already correctly stated, not closed.
# This section closes it partially, using three newly-added papers in
# Papers/Economics/ that were not available for the earlier passes:
#
#   1. Russek & Zimm, "Potential for cost effective magnetocaloric air
#      conditioning systems", Int. J. Refrig. 29 (2006) 1366-1373 -- Table 1
#      reports the DOE-evaluated MANUFACTURED cost of 3-ton (10.55 kW_c)
#      residential air conditioners: $42.6-102.1/kW_c across SEER 10-18.
#      The paper's own worked examples put aggregate magnet+MCM material
#      cost at $2.3-11.7/kW_c -- i.e. "less than 10% of the manufactured
#      cost for an SEER 13 system" (the paper's own words, paraphrased).
#      That ratio is used below as an ORDER-OF-MAGNITUDE full-system
#      multiplier: if materials are a MINORITY share (<=10%) of a
#      manufactured, mass-produced HVAC unit's cost, a materials-only
#      number can be scaled up by roughly 1/0.10 = 10x as a rough
#      full-system-cost ESTIMATE, not a bottom-up BOM. This is explicitly
#      an analogy to a *different, more mature, mass-produced* technology
#      (vapor-compression AC), not a magnetic-refrigerator-specific
#      HX/pump/motor/controls quote -- flagged, not hidden.
#   2. Silva, Bahl et al. (the "Permanent magnet design for magnetic heat
#      pumps using total cost minimization" paper, J. Magn. Magn. Mater.
#      442 (2017) 87-96) gives a thermoeconomic "cost of exergetic
#      cooling" formulation (their Eq. 6, following Rowe's 2011
#      thermoeconomic approach) that annualizes CAPITAL cost via a Capital
#      Recovery Factor (CRF) and adds it to electricity OPEX on a
#      consistent $/kWh basis, and explicitly adds a THIRD material
#      category -- soft magnetic material (SMM, e.g. 1018 steel flux
#      return yoke) at ~$5/kg -- alongside permanent-magnet and MCM costs.
#      `bom_cost()` below adds that SMM term to `material_cost()`'s two
#      existing terms (magnet, MCM) as a genuinely additive, literature-
#      sourced BOM line-item, not a replacement.
#   3. Russek & Zimm also report a literature unit cost for La(Fe,Si)13Hy
#      of ~$8/kg (vs. Gd's $20/kg already used above) -- used in
#      `MCM_COST_PER_KG_BY_FAMILY` below so `bom_cost()`/`cost_index()`
#      can price a design that uses a different composition-tunable
#      giant-MCE family (see core/cascade.py's GD_FAMILY/LAFESIH_FAMILY/
#      MNFEPSI_FAMILY and core/optimize.py's Phase 15 material-family
#      co-optimization) rather than assuming Gd for every candidate.
#
# Honesty note: none of this closes the ROADMAP's full bottom-up BOM gap
# (a real HX/pump/motor/controls parts-and-labor cost breakdown SPECIFIC
# to an AMR device). It adds (a) one genuinely new, literature-sourced
# material-cost line item (SMM/yoke), (b) an explicit, clearly-labeled
# ORDER-OF-MAGNITUDE full-system estimate derived from a different
# (vapor-compression) technology's manufactured-cost benchmark, used only
# as a sanity-check multiplier, and (c) a second, independent capital+
# operating cost methodology (CRF-based levelized cost of cooling) that
# can be cross-checked against `lifetime_cost()`'s simpler
# materials-floor-plus-electricity approach. A specific, AMR-native BOM
# remains future work (see ROADMAP.md Phase 15).

COST_SMM_PER_KG = 5.0   # $/kg, 1018 soft-magnetic (flux-return) steel,
                          # Silva et al., J. Magn. Magn. Mater. 442 (2017) 87-96

# Non-materials cost fraction implied by Russek & Zimm (2006) Table 1 +
# their own worked examples (aggregate magnet+MCM cost "less than 10% of
# the manufactured cost for an SEER 13 system"). Used as a single
# representative multiplier; the paper's own range (Gd: $11.7/kW_c,
# La(Fe,Si)13Hy: $2.3/kW_c, against $42.6-102.1/kW_c manufactured cost)
# spans roughly 8x-40x, so 10x is a conservative (low) point within that
# range, not an upper bound.
NON_MATERIALS_COST_MULTIPLIER = 10.0

# MCM unit costs by material family, $/kg. Gd from Bjork et al. (2011,
# see module docstring); La(Fe,Si)13Hy from Russek & Zimm (2006) -- see
# section docstring above. Gd5(SixGe1-x)4(-Ga) (GD_FAMILY in cascade.py)
# has no independently-sourced $/kg in this corpus; it is a Gd-based
# alloy, so Gd's price is used as an explicitly-flagged proxy rather than
# inventing a number. (Mn,Fe)2(P,Si) (MNFEPSI_FAMILY) is repeatedly
# described in this corpus only QUALITATIVELY as low-cost / abundant,
# non-rare-earth (e.g. "Materials challenges for high performance
# magnetocaloric refrigeration devices" -- comparatively low materials
# cost, no digit given; "Impact of F and S Doping on (Mn,Fe)2(P,Si) Giant
# Magnetocaloric Materials" -- low raw material costs, absence of rare
# earths, again no digit) -- no $/kg figure was found anywhere in this
# project's corpus, so it is ALSO left at the Gd price as a conservative
# (i.e. not artificially cheap) placeholder rather than guessing a number
# the qualitative literature only implies should be lower.
MCM_COST_PER_KG_BY_FAMILY = {
    "Gd": COST_MCM_PER_KG,                      # $20/kg, Bjork et al. 2011
    "Gd5(SixGe1-x)4(-Ga)": COST_MCM_PER_KG,      # proxy: Gd-based alloy, no
                                                   # independent source found
    "La(Fe,Si)13Hy": 8.0,                        # $/kg, Russek & Zimm 2006
    "(Mn,Fe)2(P,Si)": COST_MCM_PER_KG,           # proxy: no $/kg source found
                                                   # (qualitatively "low cost"
                                                   # in the literature only)
    # Phase 24: Ga1-xCMn3+x antiperovskite (core/antiperovskite_material.py).
    # Wang et al. (2009) describe the raw materials (Ga, C, Mn) only
    # qualitatively as "inexpensive and innoxious" -- no $/kg digit is given,
    # same situation as (Mn,Fe)2(P,Si) above. Left at the Gd price as the
    # SAME conservative (not artificially cheap) placeholder convention used
    # for (Mn,Fe)2(P,Si), rather than inventing a number the qualitative
    # literature only implies should be lower.
    "Ga1-xCMn3+x": COST_MCM_PER_KG,
    # Phase 25: Mn1-xCuxCoGe (core/first_order_mce.py's MNCUCOGE_FIRST_ORDER).
    # No $/kg digit located for this exact composition either. Mn/Co/Ge (and
    # a small Cu fraction) are all non-precious, non-rare-earth elements --
    # qualitatively cheap, like (Mn,Fe)2(P,Si) above -- so the SAME
    # conservative Gd-price placeholder convention is used rather than
    # inventing a lower number the qualitative chemistry only suggests.
    "Mn1-xCuxCoGe": COST_MCM_PER_KG,
}


def material_cost_by_family(mu0H_max, mass_regenerator, family_name="Gd"):
    """Same magnet-mass scaling as `material_cost()`, but looks the MCM
    unit cost up by family name (see MCM_COST_PER_KG_BY_FAMILY) instead of
    always assuming Gd's $20/kg. `family_name` should match a
    core.cascade.GradedFamily.name (or "Gd" for plain gadolinium) -- an
    unrecognized name falls back to Gd's cost with a printed warning
    rather than raising, so callers sweeping many candidate designs don't
    crash on an unexpected label."""
    mcm_cost_per_kg = MCM_COST_PER_KG_BY_FAMILY.get(family_name)
    if mcm_cost_per_kg is None:
        print(f"material_cost_by_family: unrecognized family_name={family_name!r}, "
              f"falling back to Gd's ${COST_MCM_PER_KG}/kg")
        mcm_cost_per_kg = COST_MCM_PER_KG
    magnet_mass = MAGNET_TO_MCM_MASS_RATIO_PER_TESLA * mu0H_max * mass_regenerator
    return COST_MAGNET_PER_KG * magnet_mass + mcm_cost_per_kg * mass_regenerator


def bom_cost(mu0H_max, mass_regenerator, family_name="Gd",
             smm_mass_fraction=0.5):
    """Bottom-up magnet + MCM + soft-magnetic-material (SMM, flux-return
    yoke) cost, $ -- extends `material_cost()`/`material_cost_by_family()`
    with the SMM line item from Silva et al. (2017) (see section
    docstring). `smm_mass_fraction` is the assumed SMM mass as a fraction
    of the magnet mass (soft-iron flux-return yokes are typically a
    comparable order of magnitude to the permanent-magnet mass in
    Halbach-style and C-core magnet circuits, but no single ratio is
    reported across all the papers in this corpus for a general AMR
    design -- 0.5 is a rough, explicitly-flagged mid-of-plausible-range
    placeholder, not a fitted value). Still a materials-only cost -- see
    section docstring for what remains excluded."""
    mcm_cost_per_kg = MCM_COST_PER_KG_BY_FAMILY.get(family_name, COST_MCM_PER_KG)
    magnet_mass = MAGNET_TO_MCM_MASS_RATIO_PER_TESLA * mu0H_max * mass_regenerator
    smm_mass = smm_mass_fraction * magnet_mass
    magnet_cost = COST_MAGNET_PER_KG * magnet_mass
    mcm_cost = mcm_cost_per_kg * mass_regenerator
    smm_cost = COST_SMM_PER_KG * smm_mass
    return {
        "magnet_mass_kg": round(magnet_mass, 4),
        "smm_mass_kg": round(smm_mass, 4),
        "magnet_cost_$": round(magnet_cost, 2),
        "mcm_cost_$": round(mcm_cost, 2),
        "smm_cost_$": round(smm_cost, 2),
        "materials_bom_total_$": round(magnet_cost + mcm_cost + smm_cost, 2),
    }


# =============================================================================
# Phase 30 addition: AMR-NATIVE bottom-up lifetime-cost model
# =============================================================================
#
# Everything above this point either (a) prices magnet+MCM+SMM materials
# only (`bom_cost()`), or (b) scales that materials-only BOM by a multiplier
# BORROWED from a different, more mature technology (vapor-compression AC
# manufactured-cost benchmarks, `full_system_cost_estimate()`). Neither is
# an AMR-specific bottom-up cost. This section closes that gap using a
# paper genuinely specific to AMR devices that was in this project's own
# corpus but not yet used for costing:
#
#   Bjørk, Bahl & Nielsen, "The lifetime cost of a magnetic refrigerator",
#   Int. J. Refrig. 63 (2016) 48-62.
#
# This is a DTU numerical-optimization study of a 25W-average-load AMR
# device's TOTAL lifetime cost (capital: magnet + MCM; operating:
# electricity), independently varying magnet/MCM unit price and device
# lifetime, and reporting the actual optimized capital+operating cost
# split -- not a manufactured-cost multiplier borrowed from a different
# technology. Their own headline numbers (Abstract, Section 4, magnet/MCM
# price $40/$20 per kg):
#   - Capital cost: ~$100 (magnet) + ~$40 (MCM) at a representative
#     4.5 Hz, utilization=0.35, COP=2 operating point.
#   - Operating cost: ~$0.004/hour of runtime.
#   - Total 15-year lifetime cost: $150-$400, depending on magnet/MCM
#     price (magnet dominates; MCM cost is "almost negligible" -- their
#     own wording).
#   - Their own rough VCC comparison point: an A+++ compression-based
#     appliance uses ~$113 of electricity over 15 years at 8.6W and
#     $0.10/kWh, plus a ~$30 compressor -- i.e. ~$143 total, similar
#     order of magnitude to the cheapest AMR configuration they find.
#
# This is for a SMALL (25W-average) appliance-scale device, NOT this
# repo's own data-center-scale (~kW) design points -- the numbers below
# are used as a per-unit-cost-structure REFERENCE (magnet-dominated
# capital cost, small MCM share, small operating-cost share at low
# utilization), not rescaled to kW scale here, since Bjørk et al.
# themselves do not provide a scaling law for that jump and inventing one
# would be a bigger honesty violation than leaving the gap open.

BJORK2016_REFERENCE_DEVICE = {
    "average_cooling_power_W": 25.0,
    "capital_cost_magnet_usd": 100.0,
    "capital_cost_mcm_usd": 40.0,
    "operating_cost_usd_per_hour": 0.004,
    "device_lifetime_years": 15.0,
    "lifetime_cost_range_usd": (150.0, 400.0),
    "frequency_Hz": 4.5,
    "utilization": 0.35,
    "COP": 2.0,
    "field_T_at_optimum": 1.4,
    "magnet_price_usd_per_kg": 40.0,
    "mcm_price_usd_per_kg": 20.0,
    "source": "Bjork, Bahl & Nielsen, Int. J. Refrig. 63 (2016) 48-62, "
              "Abstract + Section 4 + Fig. 8/Table 2.",
}

VCC_REFERENCE_APPLIANCE_BJORK2016 = {
    "class": "A+++ compression-based appliance (rough comparison point, "
             "Bjork et al.'s own words, not this repo's own baseline_cooling.py)",
    "power_W": 8.6,
    "electricity_price_usd_per_kWh": 0.10,
    "lifetime_electricity_cost_usd": 113.0,
    "compressor_capital_cost_usd": 30.0,
    "lifetime_years": 15.0,
    "total_lifetime_cost_usd": 143.0,
    "source": "Bjork, Bahl & Nielsen (2016), Section 6 discussion, citing "
              "Vincent & Heun (2006) for the compressor price.",
}


def amr_native_lifetime_cost_reference(verbose=True):
    """Report the Bjork et al. (2016) AMR-native lifetime-cost structure
    directly, as a genuine (not borrowed-multiplier) bottom-up reference
    point -- for the paper's economics section, to be cited ALONGSIDE (not
    instead of) `full_system_cost_estimate()`'s order-of-magnitude VCC-
    multiplier estimate, each disclosing a different limitation."""
    d = BJORK2016_REFERENCE_DEVICE
    v = VCC_REFERENCE_APPLIANCE_BJORK2016
    capital = d["capital_cost_magnet_usd"] + d["capital_cost_mcm_usd"]
    magnet_share_of_capital = d["capital_cost_magnet_usd"] / capital
    if verbose:
        print("AMR-native lifetime-cost reference (Bjork, Bahl & Nielsen, "
              "Int. J. Refrig. 63 (2016) 48-62):")
        print(f"  device class: {d['average_cooling_power_W']:.0f}W-average "
              f"appliance-scale AMR (NOT this repo's kW-scale data-center "
              f"design points -- reference structure only, not rescaled)")
        print(f"  capital cost: ${d['capital_cost_magnet_usd']:.0f} magnet + "
              f"${d['capital_cost_mcm_usd']:.0f} MCM = ${capital:.0f} "
              f"(magnet is {magnet_share_of_capital*100:.0f}% of capital cost)")
        print(f"  operating cost: ${d['operating_cost_usd_per_hour']:.3f}/hour")
        print(f"  15-year lifetime cost range: "
              f"${d['lifetime_cost_range_usd'][0]:.0f}-"
              f"${d['lifetime_cost_range_usd'][1]:.0f}, depending on "
              f"magnet/MCM unit price")
        print(f"  their own rough VCC comparison (A+++ appliance): "
              f"~${v['total_lifetime_cost_usd']:.0f} total lifetime cost "
              f"(${v['lifetime_electricity_cost_usd']:.0f} electricity + "
              f"${v['compressor_capital_cost_usd']:.0f} compressor)")
        print("  HONEST FRAMING FOR THE PAPER: this is a genuine AMR-native "
              "bottom-up cost study (magnet + MCM capital + electricity "
              "OPEX, numerically optimized), unlike full_system_cost_"
              "estimate()'s borrowed VCC-manufactured-cost multiplier -- "
              "but it is small-appliance-scale (25W average), not "
              "data-center scale (~kW), and does NOT include HX/pump/"
              "motor/controls/enclosure costs any more than this repo's "
              "own bom_cost() does (Bjork et al. themselves scope their "
              "cost model to magnet+MCM+electricity only -- see their own "
              "Section 2). Cite it as independent qualitative support "
              "(magnet cost dominates, MCM cost is a small fraction, "
              "operating cost is comparable in order of magnitude to "
              "capital cost over a 15-year life) rather than as a "
              "kW-scale dollar figure for this repo's own design points.")
    return {
        "capital_cost_usd": capital,
        "magnet_share_of_capital": round(magnet_share_of_capital, 3),
        "lifetime_cost_range_usd": d["lifetime_cost_range_usd"],
        "vcc_comparison_total_usd": v["total_lifetime_cost_usd"],
    }


def full_system_cost_estimate(mu0H_max, mass_regenerator, family_name="Gd",
                                smm_mass_fraction=0.5,
                                non_materials_multiplier=NON_MATERIALS_COST_MULTIPLIER):
    """ORDER-OF-MAGNITUDE full-system cost estimate: applies
    `NON_MATERIALS_COST_MULTIPLIER` (see section docstring -- derived from
    Russek & Zimm's 2006 vapor-compression-AC manufactured-cost benchmark,
    NOT an AMR-specific quote) to `bom_cost()`'s materials-only BOM total.
    This is explicitly a SANITY-CHECK ESTIMATE for "is this design's
    materials cost small enough that a full system could plausibly be
    cost-competitive", not a substitute for a real bottom-up AMR BOM.
    Returns the BOM breakdown alongside the scaled estimate so the
    materials-only number stays visible, not buried inside a single
    opaque total."""
    bom = bom_cost(mu0H_max, mass_regenerator, family_name, smm_mass_fraction)
    full_system = bom["materials_bom_total_$"] * non_materials_multiplier
    return {
        **bom,
        "non_materials_multiplier": non_materials_multiplier,
        "full_system_cost_estimate_$": round(full_system, 2),
        "note": "full_system_cost_estimate_$ is an ORDER-OF-MAGNITUDE estimate "
                "(materials BOM x a multiplier implied by vapor-compression AC "
                "manufactured-cost benchmarks, Russek & Zimm 2006), NOT a "
                "bottom-up AMR-specific HX/pump/motor/controls quote.",
    }


def levelized_cost_of_cooling(mu0H_max, mass_regenerator, Qc_avg_W, COP_electrical,
                                family_name="Gd", smm_mass_fraction=0.5,
                                device_lifetime_years=15.0, discount_rate=0.06,
                                capacity_factor=1.0,
                                electricity_price_per_kwh=ELECTRICITY_PRICE_PER_KWH):
    """Second, independent capital+operating cost methodology, following
    the Capital Recovery Factor (CRF) thermoeconomic approach of Silva et
    al. (2017) (their Eq. 6, itself following Rowe (2011)'s
    thermoeconomic cost-of-cooling method) -- complements, and can be
    cross-checked against, `lifetime_cost()`'s simpler
    "materials-floor-plus-total-lifetime-electricity" approach above.

    Unlike `lifetime_cost()` (which sums a materials floor and TOTAL
    lifetime electricity cost), this ANNUALIZES the capital cost via
        CRF = discount_rate * (1+discount_rate)^N / ((1+discount_rate)^N - 1)
    and reports a single $/kWh_cooling levelized figure combining
    annualized capital and per-kWh electricity cost -- the standard
    levelized-cost-of-X form used for comparing technologies with very
    different capital/operating cost splits (e.g. against
    `economics.simple_tco`'s $/kW CAPEX + $/kWh OPEX baselines for VCC
    and liquid cooling).

    `capacity_factor` (0-1) is the operating duty cycle, same convention
    as `lifetime_cost()`'s parameter of the same name and Silva et al.'s
    own `cf`."""
    if COP_electrical <= 0:
        raise ValueError("COP_electrical must be positive")
    if not (0.0 < discount_rate < 1.0):
        raise ValueError("discount_rate must be in (0, 1)")
    bom = bom_cost(mu0H_max, mass_regenerator, family_name, smm_mass_fraction)
    capital_cost = bom["materials_bom_total_$"]
    N = device_lifetime_years
    r = discount_rate
    CRF = r * (1 + r) ** N / ((1 + r) ** N - 1)
    annual_cooling_kWh = (Qc_avg_W / 1000.0) * 24 * 365 * capacity_factor
    if annual_cooling_kWh <= 0:
        raise ValueError("Qc_avg_W and capacity_factor must give positive annual cooling energy")
    annualized_capital_per_kwh = (CRF * capital_cost) / annual_cooling_kWh
    electricity_per_kwh_cooling = electricity_price_per_kwh / COP_electrical
    levelized_cost_per_kwh = annualized_capital_per_kwh + electricity_per_kwh_cooling
    return {
        "materials_bom_$": capital_cost,
        "CRF": round(CRF, 5),
        "annualized_capital_$_per_kwh_cooling": round(annualized_capital_per_kwh, 5),
        "electricity_$_per_kwh_cooling": round(electricity_per_kwh_cooling, 5),
        "levelized_cost_of_cooling_$_per_kwh": round(levelized_cost_per_kwh, 5),
        "device_lifetime_years": N,
        "discount_rate": r,
        "note": "Materials-only capital basis (see bom_cost()); comparable in "
                "form, not necessarily in scope, to simple_tco()'s "
                "capex_per_kw_cooling/opex_per_kwh_cooling baselines for VCC "
                "and liquid cooling.",
    }


# =============================================================================
# phase 15 addition: raw-material resource criticality (new Economics-folder
# paper, not previously in this repo's corpus)
# =============================================================================
# Gauss, Homm & Gutfleisch, "The Resource Basis of Magnetic Refrigeration",
# J. Ind. Ecol. 21(5) (2016) 1291-1300, assess supply-criticality for the
# raw materials behind the same three magnetocaloric families this repo
# already models (Gd5(SiGe)4, La(Fe,Si)13, (Mn,Fe)2P) plus the Nd2Fe14B
# permanent-magnet material `material_cost()` already prices. Their
# headline finding, reproduced qualitatively here (this repo does not
# reproduce their criticality-index calculation itself, only the
# ranking/rationale, since the index depends on external, time-varying
# supply/demand data this repo does not maintain):
#   - Gd-based alloys (this repo's plain "Gd" and cascade.py's GD_FAMILY,
#     itself Gd-based) are flagged as DISQUALIFIED as a mass-market
#     refrigerant on resource-criticality grounds -- heavy rare-earth Gd
#     supply is small and concentrated.
#   - La- and Mn-based alloys (this repo's LAFESIH_FAMILY, MNFEPSI_FAMILY)
#     are found MUCH LESS problematic (La is a light rare earth with much
#     larger supply; Mn/P/Fe/Si are all abundant, non-rare-earth elements).
#   - The Nd2Fe14B permanent magnet (this repo's `COST_MAGNET_PER_KG`
#     material) would only face a significant supply bottleneck at a LATER
#     innovation stage, once magnetic cooling has captured a large share of
#     the domestic refrigerator/AC market -- not at this repo's current
#     lab/early-product scale.
# This is a qualitative, NON-cost input: it does not change
# `MCM_COST_PER_KG_BY_FAMILY` (today's market price does not yet reflect
# this longer-run supply risk), but is directly relevant to
# optimize.py's phase 15 material-choice search, since a Pareto-optimal
# design that happens to prefer Gd on today's COP/Qc/cost objectives alone
# may be a poor long-run choice on a criticality-aware objective this
# repo does not currently model (see ROADMAP.md Future Work).
RESOURCE_CRITICALITY_BY_FAMILY = {
    "Gd": "high risk -- Gauss, Homm & Gutfleisch (2016) disqualify Gd-based "
          "alloys as a mass-market refrigerant on resource-criticality "
          "grounds (small, concentrated heavy-rare-earth supply)",
    "Gd5(SixGe1-x)4(-Ga)": "high risk -- Gd-based alloy, same heavy-rare-earth "
                            "supply constraint as plain Gd",
    "La(Fe,Si)13Hy": "much lower risk -- La is a light rare earth with "
                      "substantially larger, less concentrated supply",
    "(Mn,Fe)2(P,Si)": "much lower risk -- no rare-earth content (Mn/P/Fe/Si "
                       "are all abundant elements)",
}


def resource_criticality_note(family_name="Gd"):
    """phase 15: returns the qualitative resource-criticality rationale
    for `family_name` from Gauss, Homm & Gutfleisch (2016) (see module
    section docstring above). Falls back to a "not assessed" string for
    an unrecognized family_name rather than raising, matching
    `material_cost_by_family()`'s own fallback convention."""
    return RESOURCE_CRITICALITY_BY_FAMILY.get(
        family_name, f"not assessed in this repo's corpus for {family_name!r}")


# =============================================================================
# Phase 22 item 3: amorphous-material cost/performance note (qualitative
# only, per phase_plan.md's own scoping -- "worth a one-line cost/
# performance note in economics.py rather than a full model")
# =============================================================================
#
# HONESTY FLAG (book access, same tier as Phases 17-22's own flags): Tishin
# & Spichkin (2003) Ch. 9 (amorphous magnetic materials) is this item's
# named source. Re-confirmed directly for this pass: pdfplumber extracts
# zero characters from every page of this project's copy sampled (0, 1, 2,
# 50, 51) -- the same image-only-PDF finding already recorded for Tishin
# Ch. 11 (Phase 21), Sect. 2.8 (Phase 22 item 1), and Sect. 2.9/Ch. 10
# (Phase 22 item 2). Ch. 9's specific reported materials/numbers could not
# be read or digitized here. What follows is a general, qualitative,
# well-established materials-science characterization of amorphous
# (melt-spun ribbon / metallic-glass) magnetic alloys relative to their
# crystalline counterparts -- not a reproduction of Ch. 9's own content --
# kept deliberately to a short qualitative note rather than a cost model,
# per phase_plan.md's own explicit scoping of this as the lowest-priority
# item in Phase 22 with "no clear near-term payoff for the data-center
# application specifically."
#
# The general trade-off: amorphous (melt-spun / rapidly-quenched)
# magnetic alloys are cheaper to MANUFACTURE than the single- or
# poly-crystalline MCM families this repo already prices
# (MCM_COST_PER_KG_BY_FAMILY above) -- melt-spinning is a single-step,
# continuous process that skips the slow directional solidification /
# single-crystal-growth / long high-temperature annealing steps
# crystalline rare-earth and La(Fe,Si)13-type MCMs typically need -- but
# the lack of long-range crystalline order that makes this cheap to
# produce is the SAME structural feature that broadens and shallows the
# magnetic phase transition, so amorphous MCM candidates generally show a
# LOWER peak DeltaS_M / DeltaT_ad than a well-ordered crystalline sample
# of a comparable composition (this is the same broadening-vs-peak-height
# trade-off core/inhomogeneous_broadening.py, Phase 22 item 1, and
# core/nanocomposite_material.py, Phase 22 item 2, already quantify for
# OTHER broadening mechanisms in this repo -- amorphous structural
# disorder is a third, distinct source of the same qualitative trade-off,
# not a new mechanism this repo models separately).
AMORPHOUS_MATERIAL_COST_PERFORMANCE_NOTE = (
    "Amorphous (melt-spun ribbon / metallic-glass) MCM candidates trade "
    "LOWER manufacturing cost (melt-spinning is a single continuous step, "
    "skipping the slow single-crystal-growth/annealing this repo's priced "
    "crystalline families need) for LOWER peak DeltaS_M/DeltaT_ad -- "
    "structural disorder that makes them cheap to produce also broadens "
    "and shallows the magnetic transition, the same qualitative peak-vs-"
    "width trade-off Phase 22 items 1-2 already quantify for random "
    "grain-to-grain Tc inhomogeneity and deliberate multi-phase blending "
    "respectively. No amorphous-MCM $/kg figure or DeltaS_M value is "
    "digitized anywhere in this repo's corpus (Tishin Ch. 9, the natural "
    "source, is an image-only PDF here -- see this section's own honesty "
    "flag above), so this is recorded as a qualitative note only, NOT "
    "wired into MCM_COST_PER_KG_BY_FAMILY or any cost/performance "
    "calculation -- adding a numeric placeholder here, unlike the "
    "already-sourced entries in that dict, would be inventing a number "
    "this repo has no basis for, which core/economics.py's existing "
    "MNFEPSI_FAMILY/GD_FAMILY entries already deliberately avoid doing "
    "for their own missing figures (see MCM_COST_PER_KG_BY_FAMILY's own "
    "comment above)."
)


def amorphous_material_cost_performance_note():
    """Phase 22 item 3: returns AMORPHOUS_MATERIAL_COST_PERFORMANCE_NOTE
    (see the section docstring above for scope, sourcing, and why this is
    a qualitative note rather than a cost model or a new
    MCM_COST_PER_KG_BY_FAMILY entry)."""
    return AMORPHOUS_MATERIAL_COST_PERFORMANCE_NOTE


# =============================================================================
# Phase 19 addition: geometric (Halbach-cylinder) magnet-mass term
# =============================================================================
#
# `material_cost()`/`bom_cost()` above scale magnet mass LINEARLY with
# mu0H_max via MAGNET_TO_MCM_MASS_RATIO_PER_TESLA -- explicitly documented
# as "a rough fit to Bjork et al.'s two worked examples," not a physical
# model. ROADMAP.md's Phase 19 plan named the resulting gap directly:
# "achieving high mu0H should cost nonlinearly more magnet mass for a
# fixed air-gap geometry, which is physically real and currently absent."
#
# `core/magnet_geometry.py`'s new `halbach_field_vs_mass()` (a standard,
# closed-form idealized-Halbach-cylinder relation -- see that module's
# own honesty flags for what it is and is not sourced from) closes this
# specific gap. The functions below are NEW, ADDITIVE entry points
# (`*_geometric` suffix) rather than in-place replacements of
# `material_cost()`/`bom_cost()`/`full_system_cost_estimate()` -- unlike
# the original Phase 19 plan's literal wording ("Replace economics.py's
# current flat $/kg-with-a-ratio-fudge-factor..."), keeping the existing
# functions' exact numeric behavior unchanged avoids silently changing
# every existing caller's $ figures (main.py steps 5/5b, economics.py's
# own `lifetime_cost()`/`levelized_cost_of_cooling()`, and every existing
# test) with no explicit opt-in -- the same "new parameter/function,
# old default preserved" backward-compatibility discipline this repo has
# used consistently since Phase 15 (`pumping_power_override`,
# `cycle_type="brayton"`, `thermal_diode=None`). `core/optimize.py`'s
# `cost_index()` gained its own explicit `use_geometric_magnet_mass`
# opt-in flag for the same reason -- see that module's docstring.
#
# `air_gap_volume_m3` for the geometric relation is derived from
# `mass_regenerator` using core/thermal.py's own packed-bed volume
# convention (V_bed = mass_regenerator / (RHO_GD * (1 - porosity))) at
# thermal.py's own default porosity (0.365) and core/optimize.py's own
# default `BED_CROSS_SECTION_AREA_M2` (0.002 m^2) -- so a caller already
# using this repo's one consistent bed-geometry convention gets a
# matching magnet-bore geometry, not an independently-guessed one.

DEFAULT_POROSITY = 0.365                    # matches core/thermal.py's own default
DEFAULT_BED_CROSS_SECTION_AREA_M2 = 0.002   # matches core/optimize.py's
                                              # BED_CROSS_SECTION_AREA_M2


def geometric_magnet_mass_kg(mu0H_max, mass_regenerator,
                               porosity=DEFAULT_POROSITY,
                               bed_cross_section_area_m2=DEFAULT_BED_CROSS_SECTION_AREA_M2):
    """Magnet mass (kg) from core.magnet_geometry's closed-form idealized-
    Halbach-cylinder relation, instead of `material_cost()`'s flat
    `MAGNET_TO_MCM_MASS_RATIO_PER_TESLA * mu0H_max * mass_regenerator`
    proxy. `mass_regenerator` (kg of MCM) is converted to an air-gap
    volume using core/thermal.py's own RHO_GD-based packed-bed volume
    convention (see section docstring) before being handed to
    `core.magnet_geometry.halbach_field_vs_mass()`."""
    from core.thermal import RHO_GD
    from core.magnet_geometry import halbach_field_vs_mass
    air_gap_volume_m3 = mass_regenerator / (RHO_GD * (1 - porosity))
    return halbach_field_vs_mass(mu0H_max, air_gap_volume_m3,
                                   bed_cross_section_area_m2)["magnet_mass_kg"]


def bom_cost_geometric(mu0H_max, mass_regenerator, family_name="Gd",
                         smm_mass_fraction=0.5, porosity=DEFAULT_POROSITY,
                         bed_cross_section_area_m2=DEFAULT_BED_CROSS_SECTION_AREA_M2):
    """Geometric-magnet-mass counterpart of `bom_cost()`: same MCM and
    soft-magnetic-material (SMM) yoke line items (SMM mass still scaled
    off the magnet mass via `smm_mass_fraction`, same convention as
    `bom_cost()`), but the magnet-mass term comes from
    `geometric_magnet_mass_kg()` instead of the flat per-Tesla ratio.
    Returns the same dict shape as `bom_cost()` (with the same key names)
    so downstream callers (e.g. a `cost_index_geometric()` in
    core/optimize.py) can be written as a near drop-in swap."""
    mcm_cost_per_kg = MCM_COST_PER_KG_BY_FAMILY.get(family_name, COST_MCM_PER_KG)
    magnet_mass = geometric_magnet_mass_kg(mu0H_max, mass_regenerator, porosity,
                                             bed_cross_section_area_m2)
    smm_mass = smm_mass_fraction * magnet_mass
    magnet_cost = COST_MAGNET_PER_KG * magnet_mass
    mcm_cost = mcm_cost_per_kg * mass_regenerator
    smm_cost = COST_SMM_PER_KG * smm_mass
    return {
        "magnet_mass_kg": round(magnet_mass, 4),
        "smm_mass_kg": round(smm_mass, 4),
        "magnet_cost_$": round(magnet_cost, 2),
        "mcm_cost_$": round(mcm_cost, 2),
        "smm_cost_$": round(smm_cost, 2),
        "materials_bom_total_$": round(magnet_cost + mcm_cost + smm_cost, 2),
    }


def full_system_cost_estimate_geometric(
        mu0H_max, mass_regenerator, family_name="Gd", smm_mass_fraction=0.5,
        porosity=DEFAULT_POROSITY,
        bed_cross_section_area_m2=DEFAULT_BED_CROSS_SECTION_AREA_M2,
        non_materials_multiplier=NON_MATERIALS_COST_MULTIPLIER):
    """Geometric-magnet-mass counterpart of `full_system_cost_estimate()` --
    same ORDER-OF-MAGNITUDE non-materials multiplier (see that function's
    own docstring for its Russek & Zimm (2006) provenance and caveats),
    applied to `bom_cost_geometric()`'s materials-only BOM total instead
    of `bom_cost()`'s flat-ratio one."""
    bom = bom_cost_geometric(mu0H_max, mass_regenerator, family_name,
                               smm_mass_fraction, porosity, bed_cross_section_area_m2)
    full_system = bom["materials_bom_total_$"] * non_materials_multiplier
    return {
        **bom,
        "non_materials_multiplier": non_materials_multiplier,
        "full_system_cost_estimate_$": round(full_system, 2),
        "note": "Same order-of-magnitude caveats as full_system_cost_estimate() "
                "(see that function's docstring) -- this version's magnet-mass "
                "term additionally uses core.magnet_geometry's closed-form "
                "Halbach-cylinder relation instead of the flat per-Tesla ratio; "
                "see that module's own honesty flags.",
    }