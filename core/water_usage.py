"""
water_usage.py
================
 addition.

emissions.py already splits GWP into refrigerant + operational terms and
explicitly flags that AMR's refrigerant-free design "is a genuine
environmental benefit... but does not by itself guarantee lower total
emissions." Water is the analogous second axis reviewers now expect
alongside PUE (pue_annualized.py) and emissions (emissions.py): 2025-2026
trade coverage of both Magnotherm ("no reliance on local water supply")
and the broader data-center water-scarcity debate make this a citable,
quantifiable comparison this repo did not previously make.

This module does NOT introduce new thermodynamics. It converts each
technology's own already-computed COP (from comparison_table.csv /
pue_annualized.py) into an annual water-consumption figure via the
industry-standard Water Usage Effectiveness (WUE) metric, using
technology-appropriate WUE benchmarks from current trade/industry
reporting (Water Usage Effectiveness = liters of water consumed per kWh
of IT energy, Green Grid, 2011).

Key modeling choice, stated explicitly: AMR is water-COOLED (uses a
water/glycol heat-transfer fluid, per this repo's own emissions.py
docstring) but is NOT water-CONSUMING in the WUE sense -- it rejects heat
to a dry radiator/air-cooled heat exchanger loop (per Magnotherm's own
"Data centers consuming less energy to cool, with no reliance on local
water supply" claim, see commercial_landscape.py), not an open
evaporative cooling tower. This distinction -- water as a closed-loop
HEAT-TRANSFER MEDIUM (near-zero net consumption) vs. water CONSUMED BY
EVAPORATION to reject heat to atmosphere (the actual quantity WUE
measures) -- is the source of AMR's water advantage and is stated
explicitly here rather than left implicit, since "AMR uses water too"
is a common point of confusion this module exists partly to head off.

Sources:
  - Green Grid, "Water Usage Effectiveness (WUE): A Green Grid Data
    Center Sustainability Metric" (2011) -- WUE definition,
    WUE = Annual Site Water Usage (L) / IT Equipment Energy (kWh).
  - Industry-average WUE ~1.8 L/kWh, dominated by evaporative cooling-
    tower makeup water (80-90% of direct consumption per multiple 2025-
    2026 industry sources); best-in-class facilities with reclaimed-
    water/closed-loop or air-side-economized designs report 0.02-0.7
    L/kWh (e.g. Meta 0.26, Microsoft 0.49, AWS 0.19 L/kWh, all FY2024-
    2025 disclosures). Direct-liquid-cooled/dry-cooled closed-loop
    facilities report near-zero (Microsoft's new zero-water-evaporation
    chip-level liquid-cooling design, announced August 2024).
  - Magnotherm's own public claim (naturalrefrigerants.com, 2026; see
    core/commercial_landscape.py's COMMERCIAL_CLAIMS) that Stellar
    involves "no reliance on local water supply" -- a vendor claim, not
    an independently audited WUE figure, flagged as such below exactly
    the same way commercial_landscape.py flags its other vendor claims.

HONEST FRAMING: WUE benchmarks above are DATA-CENTER-FACILITY-LEVEL
figures (mixing IT load, humidification, and whichever cooling
technology that whole facility uses), not a like-for-like measurement of
"WUE contributed by AMR specifically" vs. "WUE contributed by a VCC/
evaporative-tower system specifically" at the SAME facility. This module
treats them as representative TECHNOLOGY-CLASS references (dry/closed-
loop rejection = near-zero, evaporative-tower rejection = ~1.8 L/kWh
industry average), the same simplification this repo's own
pue_annualized.py already applies to its ASHRAE climate-bin framing --
not a claim that swapping only the compressor for an AMR at an existing
evaporative-tower facility would, by itself, move that facility's WUE to
zero (a real facility may still use towers for the AMR's own heat
rejection loop depending on design choice; a dry/air-cooled heat
rejection loop, not the AMR unit itself, is what eliminates the water
term, and either VCC or AMR could in principle be paired with either
heat-rejection strategy).
"""

from dataclasses import dataclass
from typing import List

# Liters of water consumed per kWh of ELECTRICAL input to the cooling
# system (not IT load -- see annual_water_liters() below for how this is
# combined with COP to get a per-technology figure), by rejection
# technology class. These are current (2025-2026) industry-reported WUE
# figures for FACILITIES predominantly using each rejection strategy, used
# here as representative per-technology references (see module docstring
# "HONEST FRAMING" for the exact scope of this simplification).
WUE_L_PER_KWH_BY_REJECTION_CLASS = {
    "evaporative_tower": 1.8,
    # Industry-average WUE, ~80-90% of which is evaporative cooling-tower
    # makeup water (naturalrefrigerants/multiple 2025-2026 trade sources;
    # EESI; Green Grid methodology). This repo's VAPOR_COMPRESSION and
    # LIQUID_COOLING baselines (economics.py) are both modeled here as
    # using this rejection class by default -- see per-technology notes
    # below for why liquid cooling's OWN number is set separately.
    "closed_loop_liquid": 0.05,
    # Representative of a well-designed closed-loop direct-liquid-cooling
    # system with a dry/air-cooled or reclaimed-water heat-rejection loop
    # (e.g. best-in-class hyperscaler figures 0.02-0.26 L/kWh; Microsoft's
    # announced zero-water-evaporation chip-level liquid-cooling design,
    # Aug 2024, approaches the low end of this). NOT zero -- a small,
    # non-zero placeholder for humidification/incidental losses shared
    # across the whole facility regardless of the primary cooling
    # technology, since WUE as defined includes those terms too.
    "dry_air_cooled": 0.05,
    # Same near-zero placeholder as closed_loop_liquid, for a heat-
    # rejection loop with no evaporative component at all (dry radiator/
    # air-cooled heat exchanger) -- the rejection class this module
    # assigns to AMR by default, per Magnotherm's own "no reliance on
    # local water supply" claim (see module docstring) and the physical
    # reasoning that an AMR's own water/glycol loop is closed and
    # non-evaporative unless the SYSTEM DESIGNER chooses to reject its
    # heat via an evaporative tower (a design choice, not a property of
    # the magnetocaloric cycle itself -- stated explicitly, not implied).
}


@dataclass
class WaterUsageResult:
    technology: str
    rejection_class: str
    COP_electrical: float
    WUE_L_per_kWh_IT: float
    annual_water_liters: float
    annual_water_liters_per_kW_IT: float


def annual_water_liters(capacity_kW_IT, cop_electrical, rejection_class,
                         annual_hours=8760, avg_load_fraction=0.7):
    """Annual water CONSUMPTION (liters) implied by a technology's own
    already-computed electrical COP and its assigned rejection_class WUE
    benchmark, at a given IT capacity and duty cycle. NOTE: unlike PUE
    (which scales with 1/COP, i.e. WORSE COP means MORE electricity),
    WUE as defined here is a per-kWh-of-IT-energy water intensity
    independent of COP -- it is a property of the REJECTION TECHNOLOGY
    (does heat leave via evaporation or not), not of how much electricity
    the cooling system itself consumes. COP is therefore NOT a direct
    input to the liters/kWh_IT figure -- it is retained in the returned
    dataclass purely for cross-reference against the PUE/emissions
    modules' own COP-indexed reporting, so a reader can look up all three
    (COP, PUE, water) for the same technology/operating point in one
    place without recomputing anything."""
    annual_it_kwh = capacity_kW_IT * annual_hours * avg_load_fraction
    wue = WUE_L_PER_KWH_BY_REJECTION_CLASS[rejection_class]
    return annual_it_kwh * wue


def compare_water_usage(capacity_kW_IT, amr_cop, vcc_cop, liquid_cop,
                          amr_rejection_class="dry_air_cooled",
                          vcc_rejection_class="evaporative_tower",
                          liquid_rejection_class="closed_loop_liquid",
                          annual_hours=8760, avg_load_fraction=0.7) -> List[WaterUsageResult]:
    """Same three-technology comparison shape as emissions.compare_emissions()
    and pue_annualized.pue_comparison_table() -- takes the SAME COPs those
    modules already computed (from comparison_table.csv / step 4's
    baseline sweep) rather than recomputing anything, so all three
    (emissions, PUE, water) reports in a paper's results section are
    talking about the exact same operating point.

    Default rejection_class assignments reflect the MOST COMMON current
    industry deployment pattern for each technology (see
    WUE_L_PER_KWH_BY_REJECTION_CLASS and module docstring), not a
    thermodynamic necessity -- a caller modeling a specific real facility
    should override these (e.g. a VCC system with a dry-cooler instead of
    a tower, or an AMR paired with an evaporative tower by design choice,
    both change which number applies)."""
    results = []
    for label, cop, rc in [
        ("Magnetic (AMR)", amr_cop, amr_rejection_class),
        ("Vapor-compression", vcc_cop, vcc_rejection_class),
        ("Liquid cooling", liquid_cop, liquid_rejection_class),
    ]:
        liters = annual_water_liters(capacity_kW_IT, cop, rc, annual_hours,
                                      avg_load_fraction)
        wue = WUE_L_PER_KWH_BY_REJECTION_CLASS[rc]
        results.append(WaterUsageResult(
            technology=label, rejection_class=rc, COP_electrical=cop,
            WUE_L_per_kWh_IT=wue, annual_water_liters=round(liters, 1),
            annual_water_liters_per_kW_IT=round(
                liters / capacity_kW_IT if capacity_kW_IT > 0 else float("nan"), 1),
        ))
    return results


def write_water_usage_report(path="results/water_usage_comparison.txt",
                              capacity_kW_IT=100.0, amr_cop=4.63, vcc_cop=3.2,
                              liquid_cop=4.0):
    """capacity_kW_IT/amr_cop/vcc_cop/liquid_cop default to illustrative
    placeholders (100 kW-IT, 4.63/3.2/4.0 -- matching pue_annualized.py's
    own historical __main__ example point) ONLY when the caller doesn't
    supply real ones. main.py's stage 16 passes the actual step-4
    baseline-sweep representative_row's COPs and a facility-scale
    capacity ( fix -- these were previously hardcoded
    unconditionally, disconnecting this report from the rest of the
    pipeline's own computed operating point)."""
    is_default = (capacity_kW_IT == 100.0 and amr_cop == 4.63
                  and vcc_cop == 3.2 and liquid_cop == 4.0)
    lines = []
    lines.append(f"Water-usage (WUE) comparison, {capacity_kW_IT:.0f} kW-IT facility"
                  + ("  [ILLUSTRATIVE -- default placeholder COPs, not a computed "
                     "operating point]" if is_default else ""))
    lines.append("=" * 78)
    if is_default:
        lines.append("(replace amr_cop/vcc_cop/liquid_cop with the actual step-4 "
                      "baseline-sweep COPs at the representative span before citing "
                      "in the paper -- these are illustrative defaults matching "
                      "pue_annualized.py's own __main__ example point)")
    else:
        lines.append("(amr_cop/vcc_cop/liquid_cop below are the actual step-4 "
                      "baseline-sweep COPs at the representative span, same basis "
                      "run_economics()/run_emissions() already use -- see main.py's "
                      "_run_paper_strengthening_additions())")
    lines.append("")
    results = compare_water_usage(capacity_kW_IT, amr_cop, vcc_cop, liquid_cop)
    for r in results:
        lines.append(f"{r.technology:<20} rejection={r.rejection_class:<18} "
                      f"WUE={r.WUE_L_per_kWh_IT:5.2f} L/kWh_IT "
                      f"annual water={r.annual_water_liters:>12,.0f} L "
                      f"({r.annual_water_liters_per_kW_IT:>8,.0f} L/kW_IT/yr)")
    amr_liters = results[0].annual_water_liters
    baseline_liters = [r.annual_water_liters for r in results[1:]]
    if baseline_liters and min(baseline_liters) > 0:
        reduction_pct = 100 * (1 - amr_liters / min(baseline_liters))
        lines.append("")
        lines.append(f"AMR (dry/air-cooled rejection, default) uses "
                      f"{reduction_pct:.0f}% less annual water than the better "
                      f"of the two evaporative/closed-loop baselines above, "
                      f"under the default rejection-class assignments.")
    lines.append("")
    lines.append("HONEST FRAMING FOR THE PAPER: this compares TECHNOLOGY-CLASS "
                  "WUE references (dry rejection vs. evaporative-tower rejection "
                  "vs. closed-loop liquid), not a like-for-like measurement at "
                  "one real facility -- see this module's own docstring "
                  "'HONEST FRAMING' section. AMR's water advantage in this "
                  "framing comes from being assignable to a DRY heat-rejection "
                  "loop (a design choice enabled by, not uniquely guaranteed by, "
                  "its refrigerant-free/water-glycol-loop architecture -- "
                  "Magnotherm's own 'no reliance on local water supply' claim "
                  "is a vendor statement, not an independently audited WUE "
                  "figure, see core/commercial_landscape.py), not from the "
                  "magnetocaloric cycle itself consuming less water than a "
                  "compressor cycle would under the SAME rejection strategy.")
    text = "\n".join(lines)
    with open(path, "w") as f:
        f.write(text)
    print(text)
    print(f"\nWrote {path}")


if __name__ == "__main__":
    write_water_usage_report()
