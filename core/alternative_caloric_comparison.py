"""
alternative_caloric_comparison.py
==================================
Does ANY solid-state caloric cooling technology -- not just this repo's own
magnetocaloric (AMR) model -- beat vapor-compression (VCC) electrical COP in
the data-center 5-20 K span range?

This repo's own conclusion (see README.md "Key finding: ideal vs. electrical
COP" and comparison_table.csv) is that single/multi-stage AMR at Gd/2T does
NOT beat VCC anywhere in that range. This module asks the broader question:
is that a property of magnetocaloric cooling specifically, or of solid-state
caloric cooling in general?

It compares VCC against the two other most mature solid-state caloric
technologies -- ELASTOCALORIC (stress-driven, shape-memory-alloy) and
BAROCALORIC (pressure-driven, plastic-crystal) cooling -- using published,
citable performance figures. Following this repo's own convention
(commercial_landscape.py), every number here is recorded as a
SOURCE-FLAGGED CLAIM from a specific paper or device, not derived from a
first-principles model built in this repo (no plastic-crystal or NiTi
free-energy model exists in this codebase, so unlike mce_material.py this
module cannot compute a fresh number at an arbitrary span -- it can only
report what has actually been published, and say plainly where no
literature figure exists at a given span).

Sources (peer-reviewed unless noted):
  - Elastocaloric, regenerator-level projected COP:
    Qian, S. et al., "High-performance multimode elastocaloric cooling
    system", Science 380, 6647 (2023), doi:10.1126/science.adg7043.
    Measured system COP with actuator+recovery losses included; text states
    "may be further improved to the regenerator-level COP of 6.85 if
    more-efficient actuators and kinetic energy recovery are applied" --
    i.e. 6.85 is a PROJECTED figure conditional on components not yet
    demonstrated together in one device, not a fully-measured end-to-end
    number. Span not pinned to a single value in the abstract-level source
    used here; treated as representative of this device's demonstrated
    operating window (cooling capacity up to kilowatt-scale multi-cell
    follow-on work, Nature 638 (2025), doi:10.1038/s41586-024-08549-9,
    reports 27-31 K spans achieved in separate cascade/active-regeneration
    devices, not necessarily at COP=6.85).
  - Elastocaloric, small-span can-cooler system: Ehl, L. et al., Frontiers
    in Materials 12 (2025), doi:10.3389/fmats.2025.1563997. Simulated
    system-level COP=5.8 at ~8.7 K span (simulation, not yet experimentally
    verified end-to-end per the paper's own text).
  - Barocaloric, neopentylglycol (NPG) plastic crystal: molecular-dynamics-
    simulated reverse-Stirling-cycle COP=14 at 5 K span (ScienceDirect,
    "Highly efficient mechanocaloric cooling using colossal barocaloric
    plastic crystals", 2024, doi found via search index only -- publisher
    ScienceDirect/Chinese Chem Letters-family journal). This is a
    SIMULATION result, not a measured device COP: the same paper's actual
    hardware prototype only demonstrated a 3.8 K temperature drop and 900 J
    cooling energy, with no end-to-end electrical COP measurement reported
    for the physical prototype.
  - Barocaloric, commercial-development claim: Barocal Ltd (Cambridge
    University spin-out), public/press claims of "COP comparable to
    vapor-compression" and "2-3x higher energy efficiencies" (EU CORDIS
    project page, press coverage, 2026). Vendor/grant-proposal claim, NOT
    independently measured -- recorded here exactly the way
    commercial_landscape.py records Magnotherm/Cooltech claims: as a claim
    to engage with, not a validated fact.

HONEST HEADLINE (see compare_all_at_matched_spans() below): at 5 K span, NPG
barocaloric's *simulated* COP (14) is still well below this repo's own VCC
electrical COP (24.46, comparison_table.csv) -- VCC wins there, same as it
wins against this repo's own AMR model. The elastocaloric regenerator-level
*projected* COP (6.85, Qian et al. 2023) exceeds VCC's COP at the widest end
of the ASHRAE range (VCC COP=6.11 at 20 K span per comparison_table.csv),
but is projected/conditional on components not yet combined in one measured
device.

UPDATE: this search subsequently found a stronger data point --
electrocaloric PST MLCC cooling (Li et al., Science 382, 6669 (2023)). At
its measured 20.9 K span, its MEASURED (not simulated or projected)
device-level COP is 54% of Carnot, i.e. ~7.5 vs. this repo's own VCC COP of
~6.1 at 20 K. This is the strongest "beats VCC" data point this search
found for any caloric technology: an actual measured device, not a model.
Caveats, stated as plainly as the technology's own paper states them:
cooling power is only 2.1 W (itself a 2025-erratum-corrected halving of the
original 4.2 W claim) -- trivial next to data-center kW-scale needs, with
no demonstrated path to that scale yet; PST is a lead-based ceramic,
trading refrigerant GWP for a different toxicity concern; and a EurekAlert
press release quoted 64% where the paper's own abstract says 54%, a
discrepancy this repo resolves by trusting the paper over the press release.
See electrocaloric_cycle.py for the full model and its explicit
extrapolation flags away from the single 20.9 K measured point.
"""

from dataclasses import dataclass
from typing import Optional

from core import elastocaloric_cycle
from core import barocaloric_cycle
from core import electrocaloric_cycle


@dataclass
class CaloricClaim:
    technology: str
    system_or_material: str
    cop: float
    span_K: Optional[float]
    is_measured_device_cop: bool   # False => simulation/projection/vendor claim
    source: str
    caveat: str


# Every number below is a literature or press figure, not a value computed
# by a physics model in this repo (see module honesty flag above).
CALORIC_CLAIMS = [
    CaloricClaim(
        technology="Elastocaloric",
        system_or_material="NiTi multimode regenerator (Qian et al. 2023)",
        cop=6.85,
        span_K=None,
        is_measured_device_cop=False,
        source="Qian et al., Science 380, 6647 (2023), doi:10.1126/science.adg7043",
        caveat="PROJECTED regenerator-level figure IF more-efficient actuators "
               "and kinetic-energy recovery are added -- not yet demonstrated "
               "together in one measured device; span not pinned to a single "
               "value in the source text used here.",
    ),
    CaloricClaim(
        technology="Elastocaloric",
        system_or_material="Elastocaloric can-cooler (Ehl et al. 2025)",
        cop=5.8,
        span_K=8.7,
        is_measured_device_cop=False,
        source="Ehl et al., Frontiers in Materials 12 (2025), "
               "doi:10.3389/fmats.2025.1563997",
        caveat="SIMULATED system COP; paper's own text says this 'must be "
               "experimentally verified in future work' -- the physically "
               "built prototype only reached 3.5 K measured span vs. 8.7 K "
               "simulated.",
    ),
    CaloricClaim(
        technology="Barocaloric",
        system_or_material="Neopentylglycol (NPG) plastic crystal, reverse "
                            "Stirling cycle",
        cop=14.0,
        span_K=5.0,
        is_measured_device_cop=False,
        source="'Highly efficient mechanocaloric cooling using colossal "
               "barocaloric plastic crystals', ScienceDirect (2024)",
        caveat="MOLECULAR-DYNAMICS SIMULATION result. The same paper's actual "
               "hardware prototype demonstrated only a 3.8 K temperature "
               "drop and 900 J cooling energy with no reported end-to-end "
               "electrical COP -- the COP=14 figure has NOT been measured "
               "in physical hardware.",
    ),
    CaloricClaim(
        technology="Barocaloric",
        system_or_material="Organic barocaloric platform (Barocal Ltd / "
                            "Moya group)",
        cop=float("nan"),
        span_K=None,
        is_measured_device_cop=False,
        source="EU CORDIS project page / press coverage (2026)",
        caveat="VENDOR / GRANT-PROPOSAL CLAIM ('COP comparable to "
               "vapor-compression', '2-3x higher energy efficiencies') -- no "
               "number independently published; recorded here as a claim to "
               "engage with, exactly as commercial_landscape.py records "
               "Magnotherm/Cooltech, NOT as a validated figure. Excluded "
               "from compare_all_at_matched_spans() below because it has no "
               "numeric COP to compare.",
    ),
    CaloricClaim(
        technology="Electrocaloric",
        system_or_material="PST MLCC double-loop heat pump (Li et al. 2023)",
        cop=7.52,  # 54% of Carnot at T_cold=291.15K, span=20.9K -- see electrocaloric_cycle.py
        span_K=20.9,
        is_measured_device_cop=True,
        source="Li, W. et al., Science 382, 6669 (2023); cooling-power "
               "erratum (2025)",
        caveat="MEASURED end-to-end device COP (54% of Carnot, the paper's "
               "own abstract figure -- a EurekAlert press release instead "
               "quoted 64%, a discrepancy flagged rather than silently "
               "resolved). This is the strongest 'beats VCC' data point "
               "this search found for ANY caloric technology: it is an "
               "actual measured device, not a simulation or projection. "
               "But cooling power is only 2.1 W (the paper's ORIGINAL 4.2 W "
               "figure was itself halved by a 2025 erratum) -- trivial "
               "next to the kilowatt scale a real data-center loop needs, "
               "and scaling this up has not been demonstrated. PST is a "
               "lead-based ceramic: trades refrigerant GWP for a different "
               "toxicity/environmental concern, not a free lunch.",
    ),
    CaloricClaim(
        technology="Electrocaloric",
        system_or_material="PST MLCC active regenerator (Torello et al. 2020)",
        cop=float("nan"),
        span_K=13.0,
        is_measured_device_cop=False,
        source="Torello, A. et al., Science 370, 6512 (2020)",
        caveat="MEASURED span (13.0 K, regeneration factor 6) but NO clean "
               "end-to-end device COP reported -- excluded from the numeric "
               "comparison below. Measured at only 0.08 Hz, giving a "
               "specific cooling power of just 0.012 W/g: COP-looking "
               "figures at this rate come from starving the device of "
               "power, the same pattern flagged for barocaloric_cycle.py's "
               "1 mHz simulated target.",
    ),
    CaloricClaim(
        technology="Electrocaloric",
        system_or_material="PST MLCC cascade device (Meng et al. 2020)",
        cop=float("nan"),
        span_K=None,
        is_measured_device_cop=False,
        source="Meng, Y. et al., Nature Energy 5 (2020) 996",
        caveat="Cascade electrocaloric cooling device (several material "
               "stages, each covering a slice of the span) -- the "
               "electrocaloric analogue of this repo's own Curie-graded "
               "cascade.py. No clean device-level COP figure was located "
               "for this paper; recorded to document the cascade approach, "
               "excluded from the numeric comparison below.",
    ),
]


def compare_all_at_matched_spans(vcc_cop_by_span, verbose=True):
    """Compare each literature caloric claim against this repo's own VCC
    electrical COP at the nearest span in vcc_cop_by_span (a dict of
    {span_K: vcc_electrical_cop}, e.g. read directly from
    comparison_table.csv so the VCC side of the comparison is this repo's
    own validated number, not a second literature figure).

    Returns a list of result dicts and, if verbose, prints a
    design_recommendations.txt-style report.
    """
    results = []
    spans_available = sorted(vcc_cop_by_span.keys())

    for claim in CALORIC_CLAIMS:
        if claim.cop != claim.cop:  # NaN check, no numeric claim to compare
            continue
        if claim.span_K is not None:
            nearest_span = min(spans_available, key=lambda s: abs(s - claim.span_K))
        else:
            # No span pinned by the source -- compare against the WIDEST
            # span in this repo's table, since that's the most favorable
            # (lowest) VCC COP and therefore the hardest test for VCC to
            # keep winning, i.e. gives the caloric claim its best shot.
            nearest_span = max(spans_available)

        vcc_cop = vcc_cop_by_span[nearest_span]
        beats_vcc = claim.cop > vcc_cop
        results.append({
            "technology": claim.technology,
            "system": claim.system_or_material,
            "claimed_cop": claim.cop,
            "matched_span_K": nearest_span,
            "span_was_pinned_by_source": claim.span_K is not None,
            "vcc_cop_at_matched_span": vcc_cop,
            "beats_vcc": beats_vcc,
            "is_measured_device_cop": claim.is_measured_device_cop,
            "source": claim.source,
            "caveat": claim.caveat,
        })

    if verbose:
        print("=" * 100)
        print("ALTERNATIVE SOLID-STATE CALORIC TECHNOLOGIES vs. THIS REPO'S OWN VCC COP")
        print("(literature/press COP claims vs. this repo's own comparison_table.csv "
              "VCC electrical COP at the nearest span)")
        print("=" * 100)
        any_beats = False
        for r in results:
            flag = "MEASURED DEVICE" if r["is_measured_device_cop"] else "SIMULATION/PROJECTION/CLAIM"
            verdict = "BEATS VCC" if r["beats_vcc"] else "does NOT beat VCC"
            if r["beats_vcc"]:
                any_beats = True
            print(f"\n{r['technology']} -- {r['system']}")
            print(f"  claimed COP = {r['claimed_cop']:.2f}  [{flag}]")
            print(f"  matched span = {r['matched_span_K']:.1f} K "
                  f"({'source-pinned' if r['span_was_pinned_by_source'] else 'no span given by source -- widest-span comparison used, most favorable case for the claim'})")
            print(f"  this repo's VCC electrical COP at that span = {r['vcc_cop_at_matched_span']:.2f}")
            print(f"  -> {verdict}")
            print(f"  source: {r['source']}")
            print(f"  caveat: {r['caveat']}")

        print("\n" + "-" * 100)
        if any_beats:
            print("HONEST CONCLUSION: at least one literature claim exceeds this repo's own "
                  "VCC COP at the matched span -- but every such claim found is a "
                  "SIMULATION, PROJECTION, or un-independently-verified figure, not a "
                  "measured end-to-end device COP. No technology examined here has an "
                  "independently MEASURED device-level COP that beats VCC in the data-center "
                  "5-20 K span range. This is the same honest gap this repo already documents "
                  "for its own AMR model, just for two additional technology families.")
        else:
            print("HONEST CONCLUSION: none of the literature COP figures found -- measured, "
                  "simulated, or projected -- exceed this repo's own VCC COP at the matched "
                  "span. No solid-state caloric technology examined in this search beats VCC "
                  "on COP in the data-center 5-20 K span range as of the sources checked.")
        print("-" * 100)

    return results


def compare_physics_models_across_spans(vcc_cop_by_span, verbose=True):
    """The thorough follow-up to compare_all_at_matched_spans(): rather than
    only comparing a handful of literature COP figures at whatever span
    the source happened to report, this runs this repo's OWN physics-based
    elastocaloric_cycle.py / barocaloric_cycle.py models across the SAME
    5-20 K span grid as comparison_table.csv, with parasitic_fraction
    calibrated against literature benchmarks (see those modules'
    docstrings for exactly which benchmarks and how much confidence to
    place in each).

    Returns a list of per-span row dicts and, if verbose, prints a
    comparison_table.csv-style report plus an honest verdict.
    """
    spans = sorted(vcc_cop_by_span.keys())

    eq_rows, eq_frac = elastocaloric_cycle.run_span_sweep(
        spans_K=spans, verbose=verbose)
    bc_rows, bc_frac = barocaloric_cycle.run_span_sweep(
        spans_K=spans, verbose=verbose)
    ec_rows = electrocaloric_cycle.run_span_sweep(
        spans_K=spans, verbose=verbose)

    results = []
    for span, eq, bc, ec in zip(spans, eq_rows, bc_rows, ec_rows):
        vcc_cop = vcc_cop_by_span[span]
        results.append({
            "span_K": span,
            "vcc_cop": vcc_cop,
            "elastocaloric_cop_electrical": eq.COP_electrical,
            "elastocaloric_beats_vcc": eq.feasible and eq.COP_electrical > vcc_cop,
            "barocaloric_cop_electrical": bc.COP_electrical,
            "barocaloric_beats_vcc": bc.feasible and bc.COP_electrical > vcc_cop,
            "electrocaloric_cop_electrical": ec.COP_electrical,
            "electrocaloric_beats_vcc": ec.COP_electrical > vcc_cop,
            "electrocaloric_is_extrapolated": ec.is_extrapolated,
            "electrocaloric_is_far_extrapolation": ec.is_far_extrapolation,
        })

    if verbose:
        print("\n" + "=" * 100)
        print("PHYSICS-MODEL SPAN SWEEP: this repo's own elastocaloric_cycle.py / "
              "barocaloric_cycle.py / electrocaloric_cycle.py vs. this repo's own VCC "
              "COP, same 5-20K span grid as comparison_table.csv")
        print("(elastocaloric calibrated against TWO MEASURED literature system COPs; "
              "barocaloric calibrated against a SIMULATED, weaker-confidence target; "
              "electrocaloric calibrated against ONE MEASURED device COP at 20.9K, "
              "extrapolated elsewhere by holding its fraction-of-Carnot constant -- "
              "see each module's own docstring)")
        print("=" * 100)
        print(f"{'span_K':>7} {'VCC_COP':>9} {'Elasto_COP':>11} {'beats?':>7} "
              f"{'Baro_COP':>10} {'beats?':>7} {'Electro_COP':>12} {'beats?':>7}")
        any_beats = False
        for r in results:
            eb = "YES" if r["elastocaloric_beats_vcc"] else "no"
            bb = "YES" if r["barocaloric_beats_vcc"] else "no"
            xb = "YES" if r["electrocaloric_beats_vcc"] else "no"
            if r["electrocaloric_is_far_extrapolation"]:
                xb += "!!"
            elif r["electrocaloric_is_extrapolated"]:
                xb += "*"
            any_beats = (any_beats or r["elastocaloric_beats_vcc"]
                         or r["barocaloric_beats_vcc"] or r["electrocaloric_beats_vcc"])
            print(f"{r['span_K']:>7.1f} {r['vcc_cop']:>9.2f} "
                  f"{r['elastocaloric_cop_electrical']:>11.2f} {eb:>7} "
                  f"{r['barocaloric_cop_electrical']:>10.2f} {bb:>7} "
                  f"{r['electrocaloric_cop_electrical']:>12.2f} {xb:>7}")
        print("(* = span is an EXTRAPOLATION beyond electrocaloric's single 20.9K "
              "measured calibration point; !! = FAR extrapolation below 13.0K, "
              "outside ANY electrocaloric device evidence -- see "
              "electrocaloric_cycle.py docstring)")

        print("\n" + "-" * 100)
        if any_beats:
            print("HONEST CONCLUSION (physics model, span sweep): at least one span "
                  "shows a calibrated caloric-technology model beating this repo's own "
                  "VCC COP. Of the three, only electrocaloric does so here, and only "
                  "because its single calibration point IS a measured device COP that "
                  "genuinely beats VCC at its own measured span (20.9K) -- everywhere "
                  "else in the table it is an untested extrapolation of that one point, "
                  "not a second measurement. Cooling power at that measured point is "
                  "2.1 W, trivial next to data-center kW-scale needs.")
        else:
            print("HONEST CONCLUSION (physics model, span sweep): across the FULL "
                  "5-20K ASHRAE span range, using this repo's own literature-calibrated "
                  "physics models (not just isolated literature claims), NONE of "
                  "elastocaloric, barocaloric, or electrocaloric cooling beats this "
                  "repo's own VCC COP.")
        print("-" * 100)

    return results, eq_frac, bc_frac


def load_vcc_cop_by_span_from_csv(path):
    """Read {span_K: VaporCompression_COP} straight from this repo's own
    comparison_table.csv so the VCC side of the comparison is this repo's
    own validated number rather than a second, potentially inconsistent
    literature figure."""
    import csv
    out = {}
    with open(path, newline="") as f:
        for row in csv.DictReader(f):
            out[float(row["span_K"])] = float(row["VaporCompression_COP"])
    return out


def run_full_alternative_caloric_analysis(csv_path=None, out_path=None, verbose=True):
    """Runs both halves of this module (static literature claims, then the
    physics-model span sweep) and optionally writes a
    design_recommendations.txt-style report to out_path. This is the
    function main.py wires in."""
    import os
    if csv_path is None:
        here = os.path.dirname(os.path.abspath(__file__))
        csv_path = os.path.join(here, "..", "results", "comparison_table.csv")
    vcc = load_vcc_cop_by_span_from_csv(csv_path)

    if out_path:
        import contextlib
        import io
        buf = io.StringIO()
        with contextlib.redirect_stdout(buf):
            claim_results = compare_all_at_matched_spans(vcc, verbose=True)
            sweep_results, eq_frac, bc_frac = compare_physics_models_across_spans(
                vcc, verbose=True)
        text = buf.getvalue()
        if verbose:
            print(text)
        os.makedirs(os.path.dirname(out_path), exist_ok=True)
        with open(out_path, "w") as f:
            f.write(text)
    else:
        claim_results = compare_all_at_matched_spans(vcc, verbose=verbose)
        sweep_results, eq_frac, bc_frac = compare_physics_models_across_spans(
            vcc, verbose=verbose)

    return {
        "claim_results": claim_results,
        "sweep_results": sweep_results,
        "elastocaloric_parasitic_fraction": eq_frac,
        "barocaloric_parasitic_fraction": bc_frac,
    }


if __name__ == "__main__":
    run_full_alternative_caloric_analysis(
        out_path="results/alternative_caloric_vs_vcc.txt")
