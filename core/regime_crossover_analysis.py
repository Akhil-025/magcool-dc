"""
regime_crossover_analysis.py
==============================
Phase 34: the direct answer to "where, if anywhere, does this repo's OWN
model show magnetic cooling beating conventional cooling" -- not by
picking one favorable-looking point (as core.beverage_cooler_validation
and core.heat_pump_validation both deliberately do NOT do; both of those
modules check REALISM against real deployments, not superiority), but by
a systematic search across span and baseline-technology quality, using
ONLY this repo's own already-tested functions (core.cascade.
staged_baseline_result, core.baseline_cooling.vapor_compression_cop,
core.emissions.compare_emissions) -- no new physics, no new loss model,
just a wider, more honest sweep of the same model already used for
results/comparison_table.csv.

HONEST HEADLINE RESULT (do not skip this before reading the functions
below): there is NO region found. Two independent checks, both null:

1. COP CROSSOVER SEARCH (run_cop_crossover_search()): swept span from 3K
   to 30K, baseline vapor-compression second-law efficiency from 0.25
   (representative of small/simple/residential-grade compressors, below
   this repo's own core.baseline_cooling documented 0.35-0.45 data-center
   range -- see that module's own docstring) up to 0.55 (well-optimized
   chilled-water), AND a broad grid over this repo's own AMR design
   freedoms (mass_regenerator, frequency, fluid_mdot, mu0H_max) at EACH
   span. Confirmed: AMR_COP_electrical does not exceed VCC_COP at ANY
   combination tried, including VCC's own WORST-case (eta=0.25) setting
   at the smallest span tested (3K) -- at span=3K the best achievable
   AMR_COP_electrical found by this repo's own model, across dozens of
   design combinations, saturates at ~5.9 (a real, structural ceiling in
   this repo's own StateDependentLossModel/AMRSystem stack, not a search
   failure -- confirmed by many different (frequency, mdot, mu0H)
   combinations independently converging to the identical value), against
   VCC's worst-case COP of ~24.3 at the same span. The gap does not close
   as span narrows toward the beverage-cooler/heat-pump regime the way an
   optimistic reading of results/comparison_table.csv's own "gets
   relatively better as span narrows" trend might suggest extrapolating
   to -- it narrows in RATIO terms but a factor-of-4-plus gap remains
   even at the most favorable point found.

2. EMISSIONS CROSS-CHECK (run_emissions_crossover_check()): does
   eliminating refrigerant entirely (this repo's own core.emissions.py,
   already built, already honest about this exact risk in its own
   docstring: "if AMR has a lower COP than the baseline technologies, its
   operational emissions will also be higher... [refrigerant emissions do
   not by themselves] guarantee lower total emissions") rescue the
   comparison at beverage-cooler scale, where refrigerant charge is a
   larger fraction of a small system? Checked directly using the SAME
   COPs run_eclipse_directional_check() itself computed (AMR=1.76,
   VCC=6.66 at 0.4kW): NO -- AMR's total annual emissions come out ~3.6x
   HIGHER than vapor-compression's at this scale in this repo's own
   model, because operational (energy-driven) emissions dominate at this
   COP gap by roughly two orders of magnitude over the refrigerant-leak
   term, exactly as core.emissions.py's own docstring already warned.

WHY REPORT A NULL RESULT AS ITS OWN MODULE (same discipline as
core.hysteresis_sensitivity.run_hysteresis_paired_significance_test()'s
own null finding, or core.loss_model.run_core_plus_tusek_multipoint_diagnostic()'s
correction of an initially-too-rosy read): the real-world sources this
repo's own core.beverage_cooler_validation and core.heat_pump_validation
modules check against (Magnotherm Eclipse, Polaris, Ames Lab) ALL,
independently and consistently, frame their own results the same way --
"matches" vapor-compression on cost/weight/power-density/COP, "eliminates
refrigerant emissions" -- NONE of the three peer-reviewed/press sources
checked this session claims to BEAT vapor-compression on COP or total
emissions. This module's own null result is therefore not a modeling
failure to be explained away; it is the SAME conclusion the real
literature has independently reached, arrived at here through this
repo's own model rather than by reading it off someone else's press
release. That is a genuinely useful, publication-honest thing for this
repo to be able to say for itself, in its own words, with its own numbers
-- "not yet competitive on performance, real value is
parity-plus-decarbonization-optionality, corroborated independently
by the deployed literature" is a defensible, citable claim; "we found a
market where it wins" would not be, on this evidence.
"""

from core.cascade import staged_baseline_result
from core.baseline_cooling import vapor_compression_cop
from core.emissions import compare_emissions
from core.mce_material import GADOLINIUM


COP_SEARCH_SPANS_K = (3.0, 5.0, 8.0, 10.0, 15.0, 20.0, 25.0, 30.0)
COP_SEARCH_VCC_ETA_RANGE = (0.25, 0.35, 0.42, 0.55)
# 0.25 extends BELOW this repo's own core.baseline_cooling documented
# 0.35-0.45 data-center range, representing smaller/simpler/residential-
# grade compressors (no economizer, smaller compressor efficiency) --
# this repo has no dedicated citation for this specific lower figure, so
# it is used here as a deliberately generous (favorable-to-AMR) lower
# bound, not asserted as a specific literature value the way the
# 0.35-0.45 range is.
COP_SEARCH_MASS_KG = (0.5, 2.0, 5.0, 10.0)
COP_SEARCH_FREQ_HZ = (1.0, 2.0, 4.0)
COP_SEARCH_MDOT_KG_S = (0.02, 0.05, 0.1, 0.2)
COP_SEARCH_MU0H_T = (1.5, 2.0, 3.0)
T_AMBIENT_K = 295.0  # 21.85C, representative indoor/ambient reference


def run_cop_crossover_search(verbose=True):
    """For each span in COP_SEARCH_SPANS_K, finds the BEST achievable
    AMR_COP_electrical across every (mass, frequency, mdot, mu0H)
    combination in the grids above (a genuine, generous search of this
    repo's own AMR design freedom, not a single fixed design point), and
    compares it against vapor-compression at every eta_2nd_law in
    COP_SEARCH_VCC_ETA_RANGE, including VCC's own least-favorable
    (lowest-COP) setting. Reports every span where AMR's best achievable
    COP_electrical beats VCC's WORST-CASE COP -- if this list comes back
    empty (as it does at the time of writing -- see this module's own
    top-level docstring), that is a genuine, actively-searched-for null
    result, not merely "no test was run for narrow spans"."""
    rows = []
    for span in COP_SEARCH_SPANS_K:
        T_cold = T_AMBIENT_K - span
        best_amr_cop = 0.0
        best_design = None
        for mass in COP_SEARCH_MASS_KG:
            for freq in COP_SEARCH_FREQ_HZ:
                for mdot in COP_SEARCH_MDOT_KG_S:
                    for mu0H in COP_SEARCH_MU0H_T:
                        amr = staged_baseline_result(
                            T_cold, span, material=GADOLINIUM, mu0H_max=mu0H,
                            mass_regenerator=mass, frequency=freq, fluid_mdot=mdot,
                            regenerator_effectiveness=0.85, max_stages=4)
                        if amr.Qc > 0 and amr.COP_electrical > best_amr_cop:
                            best_amr_cop = amr.COP_electrical
                            best_design = (mass, freq, mdot, mu0H)
        vcc_cops = {eta: vapor_compression_cop(T_cold, T_AMBIENT_K, eta_2nd_law=eta).COP
                    for eta in COP_SEARCH_VCC_ETA_RANGE}
        worst_case_vcc_cop = min(vcc_cops.values())
        amr_beats_worst_case_vcc = best_amr_cop > worst_case_vcc_cop
        rows.append({
            "span_K": span, "T_cold_K": T_cold,
            "best_AMR_COP_electrical": best_amr_cop, "best_AMR_design": best_design,
            "VCC_COPs_by_eta": vcc_cops, "worst_case_VCC_COP": worst_case_vcc_cop,
            "AMR_beats_worst_case_VCC": amr_beats_worst_case_vcc,
        })
        if verbose:
            vcc_str = "  ".join(f"eta={e:.2f}:COP={c:.2f}" for e, c in vcc_cops.items())
            print(f"span={span:>5.1f}K  best_AMR_COP_electrical={best_amr_cop:6.2f} "
                  f"(mass={best_design[0]}kg freq={best_design[1]}Hz mdot={best_design[2]}kg/s "
                  f"mu0H={best_design[3]}T)  |  VCC: {vcc_str}  |  "
                  f"AMR beats VCC's OWN worst case: {amr_beats_worst_case_vcc}")

    any_crossover = any(r["AMR_beats_worst_case_VCC"] for r in rows)
    if verbose:
        print()
        if any_crossover:
            winning_spans = [r["span_K"] for r in rows if r["AMR_beats_worst_case_VCC"]]
            print(f"CROSSOVER FOUND at span(s): {winning_spans} -- re-verify this is not a "
                  "search or unit error before reporting it anywhere, since it contradicts "
                  "every real-world source checked this session (see this module's own "
                  "top-level docstring).")
        else:
            print("NO CROSSOVER FOUND at any span/design combination searched, including "
                  "against vapor-compression's own least-favorable setting. This matches "
                  "(not contradicts) every real-world source checked this session -- see "
                  "this module's own top-level docstring for why a null result here is "
                  "itself the useful, reportable finding.")
    return {"rows": rows, "any_crossover_found": any_crossover}


def run_emissions_crossover_check(capacity_kW=0.4, amr_cop=1.76, vcc_cop=6.66,
                                    liquid_cop=15.0, verbose=True):
    """Does eliminating refrigerant entirely rescue the comparison on TOTAL
    emissions, even though (per run_cop_crossover_search() above) it does
    not on COP? Uses this repo's own already-built, already-tested
    core.emissions.compare_emissions() directly -- no new emissions model.

    Defaults (capacity_kW=0.4, amr_cop=1.76, vcc_cop=6.66) reuse the exact
    numbers core.beverage_cooler_validation.run_eclipse_directional_check()
    itself computes at its own default parameters, so this check is
    testing the SAME operating point that module already validated the
    model's realism at, not a new invented one."""
    results = compare_emissions(capacity_kW, amr_cop, vcc_cop, liquid_cop)
    by_tech = {r.technology: r for r in results}
    amr_r = by_tech["Magnetic (AMR) - no refrigerant"]
    vcc_r = by_tech["Vapor-compression"]
    amr_wins_total_emissions = amr_r.total_tCO2e_per_year < vcc_r.total_tCO2e_per_year

    if verbose:
        print(f"Emissions cross-check at {capacity_kW}kW (AMR_COP={amr_cop}, VCC_COP={vcc_cop}):")
        for r in results:
            print(f"  {r.technology:<32} refrigerant={r.refrigerant_GWP_tCO2e_per_year:.4f} "
                  f"operational={r.operational_CO2_tCO2e_per_year:.4f} "
                  f"total={r.total_tCO2e_per_year:.4f} tCO2e/yr")
        ratio = amr_r.total_tCO2e_per_year / vcc_r.total_tCO2e_per_year
        print(f"AMR total emissions are {ratio:.1f}x vapor-compression's at this COP gap -- "
              f"{'AMR wins on total emissions' if amr_wins_total_emissions else 'refrigerant elimination does NOT rescue the total-emissions comparison'} "
              "at this operating point, because operational (energy-driven) emissions "
              "dominate the refrigerant-leak term by roughly two orders of magnitude here "
              "-- exactly as core/emissions.py's own docstring already warns.")

    return {"amr": amr_r, "vcc": vcc_r,
            "amr_wins_total_emissions": amr_wins_total_emissions,
            "ratio_amr_to_vcc_total_emissions": amr_r.total_tCO2e_per_year / vcc_r.total_tCO2e_per_year}


if __name__ == "__main__":
    run_cop_crossover_search()
    print()
    run_emissions_crossover_check()
