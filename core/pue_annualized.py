"""
pue_annualized.py
==================
Phase 24 addition.

Two additions data-center-engineering reviewers consistently ask for and
`comparison_table.csv` (a single design-point-at-a-time COP comparison)
does not provide:

1. PUE framing. COP is a thermodynamicist's metric; Power Usage
   Effectiveness (PUE = Total_facility_power / IT_power) is the metric a
   data-center engineering audience actually uses. This module converts a
   cooling-system COP into its *cooling-only contribution to PUE*
   (PUE_cooling_only = 1 + P_cooling/P_IT), holding all non-cooling
   overhead (power distribution, lighting, etc.) fixed and separately
   reported, per the standard PUE decomposition (Green Grid, 2016;
   ASHRAE Datacom Series).

2. Annualized / part-load comparison. A single design-point COP at one
   fixed span can understate or overstate a technology's real annual
   advantage, because outdoor-air-coupled cooling loads vary with climate
   and time of year. This module builds a simple bin-weighted annual
   energy comparison across a representative set of ASHRAE-climate-zone
   outdoor temperature bins, re-using this repo's own per-span COP
   functions (AMRSystem, vapor_compression_cop, liquid_cooling_cop) rather
   than introducing a new thermodynamic model.

Both additions are deliberately simple (bin-weighted averages, not a full
8760-hour TMY simulation) -- adequate to show whether the headline
5-20K-span conclusion is robust across a representative annual load
profile, not a claim of hourly-resolution accuracy.

Sources for climate bin structure and IT/cooling load-fraction defaults:
  - ASHRAE TC9.9, "Thermal Guidelines for Data Processing Environments",
    5th ed. (2021) -- W-class recommended/allowable envelopes.
  - Green Grid, "PUE: A Comprehensive Examination of the Metric" (2012)
    and 2016 addendum -- standard PUE decomposition and reporting practice.
  - Lawrence Berkeley National Laboratory (LBNL), "Data Center Energy
    Usage" reports -- representative non-cooling overhead PUE component
    (~0.10-0.15 of PUE from power distribution/lighting, used here as a
    fixed placeholder, flagged as such).
"""

from dataclasses import dataclass
from typing import List

import numpy as np

from core.baseline_cooling import vapor_compression_cop, liquid_cooling_cop
from core.amr_cycle import AMRSystem
from core.mce_material import GADOLINIUM
from core.loss_model import StateDependentLossModel

# --- 1. PUE framing -------------------------------------------------------

NON_COOLING_PUE_OVERHEAD = 0.12
"""Fixed placeholder for non-cooling PUE overhead (power distribution,
lighting, UPS losses), per LBNL data-center energy-usage survey figures.
NOT re-derived here -- flagged explicitly as a placeholder rather than a
site-specific measurement, exactly as this repo's own economics.py flags
its full-system BOM cost as order-of-magnitude."""


@dataclass
class PUEResult:
    technology: str
    COP: float
    PUE_cooling_only: float     # 1 + 1/COP
    PUE_total_estimate: float   # + fixed non-cooling overhead


def cop_to_pue(technology: str, COP: float,
               non_cooling_overhead: float = NON_COOLING_PUE_OVERHEAD) -> PUEResult:
    """PUE = (IT_power + cooling_power + other_overhead_power) / IT_power.
    Cooling power per unit IT power = 1/COP (COP defined as Qc/W_electrical,
    Qc taken here as equal to IT heat load, the standard data-center
    assumption for steady state). Non-cooling overhead is a separately
    reported additive term, NOT folded into COP."""
    if COP <= 0:
        pue_cool = float("inf")
    else:
        pue_cool = 1.0 + 1.0 / COP
    return PUEResult(technology, COP, round(pue_cool, 3),
                      round(pue_cool + non_cooling_overhead, 3))


def pue_comparison_table(comparison_rows: List[dict]) -> List[dict]:
    """Given comparison_table.csv-style rows (AMR_COP_electrical,
    VaporCompression_COP, LiquidCooling_COP per span), return the same rows
    augmented with PUE_cooling_only / PUE_total_estimate for each technology.
    """
    out = []
    for r in comparison_rows:
        amr_pue = cop_to_pue("AMR", r["AMR_COP_electrical"])
        vcc_pue = cop_to_pue("VaporCompression", r["VaporCompression_COP"])
        liq_pue = cop_to_pue("LiquidCooling", r["LiquidCooling_COP"])
        out.append({
            "span_K": r["span_K"],
            "AMR_PUE_cooling_only": amr_pue.PUE_cooling_only,
            "AMR_PUE_total_estimate": amr_pue.PUE_total_estimate,
            "VCC_PUE_cooling_only": vcc_pue.PUE_cooling_only,
            "VCC_PUE_total_estimate": vcc_pue.PUE_total_estimate,
            "Liquid_PUE_cooling_only": liq_pue.PUE_cooling_only,
            "Liquid_PUE_total_estimate": liq_pue.PUE_total_estimate,
        })
    return out


# --- 2. Annualized / part-load climate-weighted comparison ---------------

@dataclass
class ClimateBin:
    label: str
    T_outdoor_C: float
    annual_hours_fraction: float   # fraction of the year at ~this outdoor T


# Representative, deliberately coarse 6-bin annual outdoor-temperature
# profile for a mixed/moderate US climate zone (ASHRAE zone 4A-like,
# e.g. representative of a large fraction of US data-center siting),
# built from typical-year percentile temperature distributions reported in
# ASHRAE/DOE climate-zone summaries. NOT a full 8760-hour TMY dataset --
# a coarse approximation adequate for a directional annualized check.
REPRESENTATIVE_CLIMATE_PROFILE_4A = [
    ClimateBin("very cold (<0C)", -5.0, 0.08),
    ClimateBin("cold (0-10C)", 5.0, 0.20),
    ClimateBin("mild (10-18C)", 14.0, 0.27),
    ClimateBin("warm (18-24C)", 21.0, 0.25),
    ClimateBin("hot (24-30C)", 27.0, 0.15),
    ClimateBin("very hot (>30C)", 33.0, 0.05),
]


def _amr_cop_at(T_cold_K, span_K, mu0H_max=2.0, mass_regenerator=5.0,
                 frequency=2.0, fluid_mdot=0.08, regenerator_effectiveness=0.85):
    sys_ = AMRSystem(GADOLINIUM, mu0H_max=mu0H_max, mass_regenerator=mass_regenerator,
                      frequency=frequency, fluid_cp=4186.0, fluid_mdot=fluid_mdot,
                      regenerator_effectiveness=regenerator_effectiveness,
                      loss_model=StateDependentLossModel(), use_ntu_thermal_model=True)
    return sys_.run(T_cold_K, span_K).COP_electrical


def annualized_energy_comparison(T_it_setpoint_C=27.0,
                                  climate_profile: List[ClimateBin] = None,
                                  economizer_below_C=18.0,
                                  verbose=True):
    """Bin-weighted annual comparison of AMR vs. vapor-compression vs.
    liquid cooling, using this repo's own per-span COP models evaluated at
    the span implied by each climate bin's outdoor temperature (i.e. the
    heat-rejection span the cooling system must cover that bin).

    Below `economizer_below_C`, liquid cooling is assumed to run in
    free/economizer mode (very high effective COP, per ASHRAE TC9.9
    W-class envelopes and this repo's own baseline_cooling.py docstring);
    vapor-compression and AMR do NOT get an economizer credit here (neither
    is modeled with a free-cooling bypass mode in this repo), which is
    itself a conservative-for-AMR, generous-for-liquid-cooling framing
    worth stating explicitly in the paper's limitations.
    """
    profile = climate_profile or REPRESENTATIVE_CLIMATE_PROFILE_4A
    T_it_K = T_it_setpoint_C + 273.15
    rows = []
    for b in profile:
        raw_span_K = max(1.0, T_it_setpoint_C - b.T_outdoor_C + 8.0)
        # +8K floor term represents the minimum approach/compressor-lift the
        # system still needs even when outdoor air is cool, since none of
        # this repo's baseline correlations model a literal zero-span limit.
        # AMR span is additionally clamped to this repo's own validated
        # ASHRAE 5-20K envelope (comparison_table.csv/main.py's own sweep
        # range) -- above 20K the 0-D model is already known (README/
        # regenerator_1d.py) to hit a structural span cap and return
        # COP_electrical=0.0, which is a documented model limitation, not a
        # new finding, so those bins are excluded from the AMR annual
        # average and flagged rather than silently zeroed or crashed on.
        span_K = raw_span_K
        T_cold_K = T_it_K - span_K
        amr_span_feasible = span_K <= 20.0
        amr_cop = _amr_cop_at(T_cold_K, min(span_K, 20.0)) if amr_span_feasible else 0.0
        vcc = vapor_compression_cop(T_cold_K, T_it_K)
        liq = liquid_cooling_cop(T_cold_K, T_it_K)
        liq_cop = 25.0 if b.T_outdoor_C < economizer_below_C else liq.COP
        rows.append({
            "bin": b.label, "T_outdoor_C": b.T_outdoor_C,
            "hours_fraction": b.annual_hours_fraction,
            "span_K": round(span_K, 1),
            "AMR_COP": round(amr_cop, 2),
            "AMR_span_feasible": amr_span_feasible,
            "VCC_COP": round(vcc.COP, 2),
            "Liquid_COP": round(liq_cop, 2),
        })

    def weighted_mean_input_power_per_unit_cooling(cop_key, rows_subset, norm_hours):
        # weight by 1/COP (electrical input per unit cooling load), the
        # physically correct quantity to time-average for annual energy,
        # NOT COP itself (averaging COP directly overweights low-load
        # high-COP hours).
        return sum((r["hours_fraction"] / norm_hours) * (1.0 / r[cop_key])
                   for r in rows_subset if r[cop_key] > 0)

    amr_rows_feasible = [r for r in rows if r["AMR_span_feasible"]]
    amr_hours_covered = sum(r["hours_fraction"] for r in amr_rows_feasible)
    amr_hours_excluded = 1.0 - amr_hours_covered
    amr_annual = weighted_mean_input_power_per_unit_cooling(
        "AMR_COP", amr_rows_feasible, amr_hours_covered) if amr_hours_covered > 0 else float("nan")
    vcc_annual = weighted_mean_input_power_per_unit_cooling("VCC_COP", rows, 1.0)
    liq_annual = weighted_mean_input_power_per_unit_cooling("Liquid_COP", rows, 1.0)
    amr_effective_annual_cop = (1.0 / amr_annual) if amr_annual == amr_annual and amr_annual > 0 else float("nan")
    vcc_effective_annual_cop = 1.0 / vcc_annual
    liq_effective_annual_cop = 1.0 / liq_annual

    if verbose:
        print(f"{'bin':<20}{'T_out(C)':>10}{'hrs frac':>10}{'span(K)':>9}"
              f"{'AMR COP':>10}{'VCC COP':>10}{'Liq COP':>10}")
        for r in rows:
            print(f"{r['bin']:<20}{r['T_outdoor_C']:>10.1f}{r['hours_fraction']:>10.2f}"
                  f"{r['span_K']:>9.1f}{r['AMR_COP']:>10.2f}{r['VCC_COP']:>10.2f}"
                  f"{r['Liquid_COP']:>10.2f}")
        print()
        print(f"Annual (bin-weighted, power-averaged) effective COP:")
        if amr_hours_excluded > 1e-9:
            print(f"  AMR:              {amr_effective_annual_cop:.2f}  "
                  f"(computed over {amr_hours_covered*100:.0f}% of annual hours; "
                  f"{amr_hours_excluded*100:.0f}% excluded -- span exceeds this "
                  f"repo's validated 20K AMR envelope / hits the documented "
                  f"0-D structural span cap, see regenerator_1d.py)")
        else:
            print(f"  AMR:              {amr_effective_annual_cop:.2f}")
        print(f"  Vapor-compression:{vcc_effective_annual_cop:.2f}")
        print(f"  Liquid cooling:   {liq_effective_annual_cop:.2f} "
              f"(includes economizer-mode credit below {economizer_below_C}C)")
        print()
        print("HONEST FRAMING FOR THE PAPER: this is a coarse 6-bin "
              "approximation of one representative mixed climate (ASHRAE "
              "zone 4A-like), not an 8760-hour TMY simulation, and gives "
              "liquid cooling an economizer credit that AMR/VCC do not "
              "receive here (neither is modeled with a free-cooling bypass "
              "in this repo). AMR's annual figure additionally only covers "
              "the fraction of hours whose implied span falls inside this "
              "repo's own validated 5-20K envelope -- colder-climate bins "
              "with a larger implied span are excluded rather than forced "
              "through a model already known not to extrapolate there. "
              "This is adequate to check whether the single-design-point "
              "conclusion in comparison_table.csv survives across a "
              "representative annual load profile, but should NOT be "
              "reported as a precise annual-energy or PUE prediction "
              "without a full TMY-driven simulation and a wider-span-"
              "capable AMR model.")

    return {
        "rows": rows,
        "AMR_effective_annual_COP": (round(amr_effective_annual_cop, 3)
                                      if amr_effective_annual_cop == amr_effective_annual_cop else None),
        "AMR_annual_hours_fraction_covered": round(amr_hours_covered, 3),
        "VCC_effective_annual_COP": round(vcc_effective_annual_cop, 3),
        "Liquid_effective_annual_COP": round(liq_effective_annual_cop, 3),
    }


def write_pue_annualized_report(path="results/pue_annualized_analysis.txt",
                                  amr_cop=4.63, vcc_cop=3.2, liquid_cop=4.0):
    """PUE-framing COPs default to illustrative placeholders (4.63/3.2/4.0,
    matching this module's own historical __main__ example point) ONLY
    when the caller doesn't supply real ones. Callers with access to the
    actual step-4 baseline-sweep representative_row (main.py does) should
    pass its AMR_COP_electrical/VaporCompression_COP/LiquidCooling_COP
    values instead -- using the placeholders unconditionally was flagged
    (Phase 31 follow-up) as reporting a PUE comparison disconnected from
    this repo's own computed operating point. The annualized/part-load
    section below is unaffected -- it has always computed its own COPs
    directly via AMRSystem/vapor_compression_cop/liquid_cooling_cop, not
    from these placeholders."""
    import io, contextlib
    buf = io.StringIO()
    with contextlib.redirect_stdout(buf):
        print("PUE framing (non_cooling_overhead placeholder = "
              f"{NON_COOLING_PUE_OVERHEAD})")
        for name, cop in [("AMR", amr_cop), ("VaporCompression", vcc_cop),
                           ("LiquidCooling", liquid_cop)]:
            r = cop_to_pue(name, cop)
            print(f"  {name:<18} COP={cop:.2f}  PUE_cooling_only={r.PUE_cooling_only:.3f}"
                  f"  PUE_total_estimate={r.PUE_total_estimate:.3f}")
        print()
        print("Annualized / part-load climate-weighted comparison")
        print("-" * 55)
        annualized_energy_comparison(verbose=True)
    with open(path, "w") as f:
        f.write(buf.getvalue())
    print(buf.getvalue())
    print(f"Wrote {path}")


if __name__ == "__main__":
    write_pue_annualized_report()
