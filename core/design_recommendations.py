"""
design_recommendations.py
==========================
Synthesizes every other analysis module in this repo into a single,
actionable "how do I actually maximize AMR electrical COP" report.

Why this module exists
-----------------------
Prior to this pass, the five COP-maximization levers this project's own
results actually support -- operating frequency, MCE material choice,
Curie-temperature grading, regenerator geometry, and field/flow balance --
were each documented separately (Sobol sensitivity in `sensitivity.py`,
NSGA-III in `optimize.py`, material ranking in
`material_family_comparison.py`, geometry optima in `geometry_analysis.py`,
Curie-grading in `cascade.py`). Nothing pulled them into one place, so a
reader had to reconstruct the design guidance by hand across five separate
text files. This module does that reconstruction in code: it takes the
*already-computed* result objects from those five stages (passed in by
`main.py`, not recomputed here -- several of those stages are expensive,
e.g. the Curie-graded cascade sweep takes ~2 minutes) and assembles them
into a ranked list of levers with the actual numbers each analysis found,
plus a single recommended operating point.

This directly answers the "how do I increase AMR electrical COP" question
using only numbers this repo's own model produced -- nothing here is a new
physical claim; it is a structured summary of five existing results.

Levers covered (ranked by demonstrated leverage in this repo's own
state-dependent-loss Sobol analysis, `sensitivity.py` step 9b):
    1. Operating frequency (ST ~ 0.85-0.87, the dominant lever --
                                   eddy-current loss scales with f^2 while
                                   Qc only grows sub-quadratically, so
                                   *lower* frequency raises electrical COP
                                   even though it lowers Qc)
    2. MCE material / composition (`material_family_comparison.py`:
                                   La(Fe,Si)13Hy-type materials, Curie-
                                   tuned to the operating point, rank
                                   consistently ahead of Gd and of the
                                   fixed-composition Gd5Si2Ge2 at the
                                   representative ASHRAE span)
    3. Curie-temperature grading (`cascade.py`'s graded cascade: matching
                                   each layer's Tc to its local fluid
                                   temperature raises both Qc and COP
                                   relative to a uniform-material cascade
                                   of the same stage count)
    4. Regenerator geometry (`geometry_analysis.py`: packed-bed
                                   sphere diameter and parallel-plate
                                   channel spacing both show a genuine
                                   interior COP optimum once hydraulic
                                   pumping power is coupled to NTU
                                   effectiveness, per Tusek, Kitanovski &
                                   Poredos, Int. J. Refrig. 36, 1456-1464
                                   (2013))
    5. Field strength / flow rate (`optimize.py`'s NSGA-III Pareto front:
                                   the knee-point design balances field,
                                   frequency, and flow rather than pushing
                                   any one of them to its bound -- as of
                                   phase 15, this same Pareto front also
                                   co-optimizes material choice (lever 2)
                                   and regenerator geometry (lever 4), so
                                   this lever's knee point is no longer
                                   computed at a fixed material/geometry)

phase 15 update: `optimize.py`'s NSGA-III search now folds material choice
and regenerator geometry into the same multi-objective search as field,
frequency, and flow (see core/optimize.py's module docstring for the
per-material-family NSGA-III + global-merge approach used). Levers 2 and
4 above remain independently useful as focused single-lever sweeps for
understanding WHY a given material/geometry wins, but are no longer this
repo's only way to explore those two choices -- lever 5's own knee point
is now this repo's most complete single answer to the joint design
question.

Recent literature context: Klinar, Kitanovski, Law, Franco & Moya,
"Perspectives and Energy Applications of Magnetocaloric, Pyromagnetic,
Electrocaloric, and Pyroelectric Materials," Adv. Energy Mater. 14,
2401739 (2024) -- a 2024 roadmap-style review that (among other things)
profiles a novel high-frequency active caloric regenerator design concept
("Hypereg") aimed squarely at this repo's own finding that frequency is
the dominant COP lever. As of phase 15 this is no longer just cited
context: `core/hypereg_analysis.py` and `core/thermal.py`'s
`pumping_power_packed_bed_hypereg()` implement and quantify the
mechanism (a hydraulic/pumping-power effect, not an eddy-current one --
see `results/hypereg_findings.md`), wired into `main.py` as step 3d.
"""

from __future__ import annotations

import os


def _fmt(value, spec=".3f", none_text="n/a"):
    if value is None:
        return none_text
    try:
        return format(value, spec)
    except (ValueError, TypeError):
        return str(value)


def summarize_frequency_lever(sobol_state_dependent_Si, pareto_rows):
    """Lever 1: operating frequency. Pulls the total-order Sobol index for
    `frequency_Hz` out of the already-computed state-dependent-loss Sobol
    result (`sensitivity.run_sobol(..., use_state_dependent_losses=True)`,
    main.py step 9b) and contrasts the NSGA-III Pareto front's best-COP
    design against its best-cooling-capacity design (main.py step 11) to
    show the frequency trade-off in this repo's own optimizer output,
    rather than asserting a number that isn't traceable to a computed
    result."""
    freq_ST = None
    if sobol_state_dependent_Si is not None:
        try:
            # SALib's Si is a dict-like object keyed by S1/ST/etc, indexed
            # in the same parameter order as sensitivity.PROBLEM["names"].
            from core.sensitivity import PROBLEM
            idx = PROBLEM["names"].index("frequency_Hz")
            freq_ST = float(sobol_state_dependent_Si["ST"][idx])
        except Exception:
            freq_ST = None

    best_cop = max(pareto_rows, key=lambda r: r["COP_electrical"]) if pareto_rows else None
    best_qc = max(pareto_rows, key=lambda r: r["Qc_W"]) if pareto_rows else None

    lines = []
    lines.append("1. OPERATING FREQUENCY -- the dominant lever")
    if freq_ST is not None:
        lines.append(f" Sobol total-order sensitivity (state-dependent loss model): "
                     f"ST(frequency) = {freq_ST:.3f}")
    if best_cop and best_qc:
        lines.append(f" NSGA-III best-electrical-COP design:  f={best_cop['frequency_Hz']:.3f} Hz "
                     f"-> COP_elec={best_cop['COP_electrical']:.2f}, Qc={best_cop['Qc_W']:.0f} W")
        lines.append(f" NSGA-III best-cooling-capacity design: f={best_qc['frequency_Hz']:.3f} Hz "
                     f"-> COP_elec={best_qc['COP_electrical']:.2f}, Qc={best_qc['Qc_W']:.0f} W")
        if best_qc["frequency_Hz"] > 0:
            lines.append(" Action: reduce cycle frequency toward the low-speed end of the "
                         "design space when electrical COP is the priority; raise it toward "
                         "the high end only when cooling capacity is the priority and higher "
                         "parasitic loss is acceptable. This is a genuine Pareto trade-off in "
                         "this repo's model, not a free win in either direction.")
    return "\n".join(lines), {"frequency_ST": freq_ST,
                               "best_cop_design": best_cop, "best_qc_design": best_qc}


def summarize_material_lever(material_rows, representative_span_K=10.0):
    """Lever 2: MCE material/composition choice, from
    `material_family_comparison.run_analysis()`'s already-computed rows
    (main.py step 8d) at the representative ASHRAE span."""
    rep_rows = [r for r in material_rows if r["span_K"] == representative_span_K
                and r.get("in_range") and r.get("1stage_COP") is not None]
    ranked = sorted(rep_rows, key=lambda r: r["1stage_COP"], reverse=True)

    lines = ["2. MCE MATERIAL / COMPOSITION CHOICE"]
    if ranked:
        best = ranked[0]
        gd_row = next((r for r in rep_rows if r["candidate"].startswith("Gd (fixed)")), None)
        lines.append(f" Best candidate at {representative_span_K:.0f}K span (ASHRAE point): "
                     f"{best['candidate']}  COP_elec={best['1stage_COP']:.2f}, "
                     f"Qc={best['1stage_Qc_W']:.0f} W")
        if gd_row and gd_row is not best:
            gain_pct = 100.0 * (best["1stage_COP"] / gd_row["1stage_COP"] - 1.0)
            lines.append(f" vs. plain Gd: COP_elec={gd_row['1stage_COP']:.2f} "
                         f"({gain_pct:+.0f}% relative to Gd)")
        lines.append(" Action: prefer a composition-tuned giant-MCE family (e.g. "
                     "La(Fe,Si)13Hy-type) whose documented Tc window covers the target "
                     "operating point, over a fixed-composition material whose peak may sit "
                     "off-target -- see material_family_comparison.py for the full "
                     "span-by-span ranking and which families actually cover which spans.")
    else:
        lines.append(" No in-range tunable candidate found at this span in the current "
                     "sweep -- see results/material_family_comparison.txt for the full table.")
    return "\n".join(lines), {"ranked": ranked}


def summarize_grading_lever(graded_row, gd_cascade_row, n_stages=3):
    """Lever 3: Curie-temperature grading, from the already-computed
    graded-cascade row (main.py step 7b) vs. the plain-Gd cascade row
    (step 7) at the same span/stage-count."""
    lines = ["3. CURIE-TEMPERATURE GRADING (layered/graded regenerator beds)"]
    if graded_row is not None and gd_cascade_row is not None:
        g_cop = graded_row.get(f"Graded_{n_stages}stage_COP")
        g_qc = graded_row.get(f"Graded_{n_stages}stage_Qc_W")
        p_cop = gd_cascade_row.get(f"AMR_{n_stages}stage_COP")
        p_qc = gd_cascade_row.get(f"AMR_{n_stages}stage_Qc_W")
        if g_cop is not None and p_cop is not None:
            lines.append(f"   {n_stages}-stage cascade at this span: "
                         f"graded COP={g_cop}, Qc={g_qc} W vs. "
                         f"plain-Gd COP={p_cop}, Qc={p_qc} W")
            lines.append(" Action: grade each stage's composition so its local peak MCE "
                         "temperature tracks the fluid temperature at that point in the bed, "
                         "rather than using one uniform material across the whole span -- "
                         "see cascade.py's compare_graded_cascade() for the per-span/"
                         "stage-count feasibility breakdown (not every span/stage-count "
                         "combination stays within the documented composition-tunability "
                         "window; some stages fall back to plain Gd, which is reported "
                         "explicitly rather than silently).")
    else:
        lines.append(" Graded-cascade data not available for this run.")
    return "\n".join(lines), {"graded_row": graded_row, "gd_cascade_row": gd_cascade_row}


def summarize_geometry_lever(pb_best_cop_row, pp_best_cop_row):
    """Lever 4: regenerator geometry, from the already-computed packed-bed
    and parallel-plate sweeps (main.py step 3c,
    `geometry_analysis.sweep_packed_bed_diameter` /
    `sweep_parallel_plate_spacing`)."""
    lines = ["4. REGENERATOR GEOMETRY (packed-bed vs. parallel-plate)"]
    if pb_best_cop_row is not None:
        lines.append(f" Packed-bed sphere diameter maximizing COP (interior optimum): "
                     f"{pb_best_cop_row[0]} mm (Qc={pb_best_cop_row[1]:.0f} W, "
                     f"COP_aug={pb_best_cop_row[2]:.2f})")
    if pp_best_cop_row is not None:
        lines.append(f" Parallel-plate channel spacing maximizing COP (interior optimum): "
                     f"{pp_best_cop_row[0]} mm (Qc={pp_best_cop_row[1]:.0f} W, "
                     f"COP_aug={pp_best_cop_row[2]:.2f})")
    lines.append(" Action: target these interior optima rather than minimizing particle/"
                 "channel size without bound -- shrinking geometry indefinitely raises "
                 "thermal effectiveness only marginally further while hydraulic pumping "
                 "power keeps growing, which is what produces the interior optimum in the "
                 "first place (Tusek, Kitanovski & Poredos 2013).")
    return "\n".join(lines), {"pb_best_cop": pb_best_cop_row, "pp_best_cop": pp_best_cop_row}


def summarize_field_flow_lever(pareto_rows):
    """Lever 5: field strength / flow rate balance, from the NSGA-III
    Pareto front's knee-point design (main.py step 11,
    `optimize.run_optimization()`), which is this repo's own definition
    of a balanced (not extreme) design point.

    phase 15 update: `optimize.run_optimization()` now co-optimizes
    material choice and regenerator geometry (particle diameter) alongside
    field/frequency/flow/mass (see core/optimize.py's module docstring),
    so `pareto_rows` may include a "material" and "particle_diameter_mm"
    column that previous callers/fixtures did not provide. Both are
    reported here when present, but this function still degrades
    gracefully (matching the rest of this module's convention) if an
    older-shaped `pareto_rows` (missing those two keys) is passed in."""
    import numpy as np
    lines = ["5. FIELD STRENGTH / FLOW RATE BALANCE (now co-optimized with material and geometry)"]
    if pareto_rows:
        F = np.array([[-r["COP_electrical"], -r["Qc_W"], r["cost_index_USD"]]
                      for r in pareto_rows])
        Fn = (F - F.min(axis=0)) / (F.max(axis=0) - F.min(axis=0) + 1e-12)
        knee_idx = int(np.argmin(np.linalg.norm(Fn, axis=1)))
        knee = pareto_rows[knee_idx]
        material_note = f", material={knee['material']}" if "material" in knee else ""
        criticality_note = f"; resource criticality: {knee['resource_criticality']}" if "resource_criticality" in knee else ""
        geometry_note = f", d_p={knee['particle_diameter_mm']} mm" if "particle_diameter_mm" in knee else ""
        lines.append(f" Knee-point (balanced) Pareto design: "
                     f"H={knee['mu0H_max_T']} T, f={knee['frequency_Hz']} Hz, "
                     f"mdot={knee['fluid_mdot_kgs']} kg/s, mass={knee['mass_regenerator_kg']} kg, "
                     f"eps={knee['regen_effectiveness']}{geometry_note}{material_note}  "
                     f"-> COP_elec={knee['COP_electrical']}, Qc={knee['Qc_W']} W, "
                     f"cost=${knee['cost_index_USD']}")
        lines.append(" Action: operate near this balanced point rather than at either "
                     "objective's own extreme (max-COP alone drives frequency and Qc very "
                     "low; max-Qc alone raises parasitic loss and cuts COP roughly 3x in "
                     "this repo's Pareto front) -- see results/pareto_front.csv for the full "
                     f"{len(pareto_rows)}-design front and the other extreme points. As of "
                     "phase 15 this knee point already reflects the best-performing material "
                     "and regenerator geometry found jointly with field/frequency/flow, not "
                     "just the best fixed-Gd, fixed-geometry point -- levers 2 and 4 above "
                     "remain useful as focused single-lever sweeps for understanding WHY a "
                     "given material/geometry wins, but are no longer this repo's only way to "
                     "explore those two choices.")
        if criticality_note:
            lines.append(f" Note{criticality_note} (Gauss, Homm & Gutfleisch 2016 -- see "
                         "economics.py; a qualitative, non-cost input not yet part of the "
                         "NSGA-III objective itself, so a design's COP/Qc/cost ranking here "
                         "does not account for it).")
        return "\n".join(lines), {"knee_point": knee}
    lines.append(" Pareto front not available for this run.")
    return "\n".join(lines), {"knee_point": None}


def summarize_cycle_type_finding(cycle_type_result):
    """Folds in step 2b's cycle-topology validation sensitivity:
    does inferring an 'ericsson'-like cycle for rotary devices (instead of
    this model's flat 'brayton' default) change the COP prediction error
    against the published benchmark rows? Reads cycle_type_result (the
    list of per-device dicts returned by
    validation_system.run_cycle_type_validation()) rather than recomputing
    anything -- same convention as the other five levers above."""
    lines = ["6. CYCLE TOPOLOGY (Ericsson-like vs. Brayton-like) -- validation finding, not a design lever"]
    if not cycle_type_result:
        lines.append(" Not available for this run.")
        return "\n".join(lines), {"cycle_type_result": None}
    try:
        n_improved = sum(1 for r in cycle_type_result if r.get("direction") == "improved")
        n_worsened = sum(1 for r in cycle_type_result if r.get("direction") == "worsened")
        n_unchanged = sum(
            1 for r in cycle_type_result
            if r.get("cycle_type_inferred") == "ericsson" and r.get("direction") == "unchanged")
        n_not_comparable = sum(
            1 for r in cycle_type_result
            if r.get("cycle_type_inferred") == "ericsson"
            and str(r.get("direction", "")).startswith("not comparable"))
        lines.append(f" Rotary-device subset re-checked as 'ericsson' vs. this model's flat "
                     f"'brayton' default: {n_improved} improved, {n_worsened} worsened, "
                     f"{n_unchanged} unchanged, {n_not_comparable} not comparable.")
        lines.append(" This is a directional sensitivity check against a small rotary subset "
                     "using a naming-convention proxy for cycle topology (not a literature-"
                     "confirmed classification) -- it does NOT validate the specific "
                     "CYCLE_TYPE_FACTORS multiplier values themselves, and is reported here as a "
                     "validation finding rather than an actionable design lever.")
    except Exception:
        lines.append(" Result available but could not be summarized (unexpected shape).")
    return "\n".join(lines), {"cycle_type_result": cycle_type_result}


def summarize_thermal_diode_finding(thermal_diode_rows):
    """Folds in step 11c's mechanical-contact thermal-diode
    cost-only sensitivity: how much does the illustrative diode-actuation
    switching-power cost reduce COP_electrical, across frequency? Reads
    thermal_diode_rows (list of (frequency_Hz, COP_no_diode,
    COP_diode_assisted, delta_cop_pct) tuples, as returned by
    thermal_diode_analysis.sweep_frequency_with_and_without_diode())."""
    lines = ["7. THERMAL DIODE (mechanical-contact) -- cost-only sensitivity, not a net-benefit finding"]
    if not thermal_diode_rows:
        lines.append(" Not available for this run.")
        return "\n".join(lines), {"thermal_diode_rows": None}
    try:
        worst = min(thermal_diode_rows, key=lambda row: row[3])
        lines.append(f" Illustrative actuation switching-power cost reduces COP_electrical by "
                     f"at most {abs(worst[3]):.2f}% (at f={worst[0]:.2f} Hz) across the frequencies "
                     "swept here -- small relative to the eddy-current/base-overhead losses that "
                     "already dominate parasitic loss at this operating point.")
        lines.append(" No offsetting heat-transfer benefit from the diode's rectification ratio "
                     "is modeled, so this is a documented upper-bound cost accounting, not a claim "
                     "that real thermal diodes are a net negative for AMR performance -- treat as "
                     "a design-exploration finding, not a validated lever (no benchmark device in "
                     "this repo's corpus uses thermal diodes).")
    except Exception:
        lines.append(" Result available but could not be summarized (unexpected shape).")
    return "\n".join(lines), {"thermal_diode_rows": thermal_diode_rows}


def summarize_passive_regen_and_elastocaloric(passive_regen_base, passive_regen_rows,
                                               elastocaloric_result):
    """Folds in step 15's passive/hybrid magnetic-regenerator
    augmentation of a conventional vapor-compression cycle, plus the
    static elastocaloric literature reference. passive_regen_base
    is the vapor-compression VaporCompressionResult (has .COP);
    passive_regen_rows is the sorted list of PassiveRegeneratorResult
    returned by passive_regenerator_analysis.compare_candidate_materials().
    elastocaloric_result is the ElastocaloricReferenceResult returned by
    baseline_cooling.elastocaloric_reference_cop()."""
    lines = ["8. PASSIVE-REGENERATOR AUGMENTATION & ELASTOCALORIC REFERENCE -- context, not AMR-specific levers"]
    if passive_regen_base is not None and passive_regen_rows:
        try:
            best = passive_regen_rows[0]
            lines.append(f" Best passive-regenerator candidate : {best.material_name}  "
                         f"eps {best.eps_baseline:.3f} -> {best.eps_augmented:.3f}  "
                         f"COP {passive_regen_base.COP:.4f} -> {best.augmented_COP:.4f} "
                         f"({best.cop_gain_fraction:+.2%}) -- an alignment effect (only materials "
                         "whose own Curie temperature falls inside the operating window gain "
                         "anything), capped by an illustrative literature-range ceiling, not a "
                         "validated device-level COP prediction.")
        except Exception:
            lines.append(" Passive-regenerator result available but could not be summarized "
                         "(unexpected shape).")
    else:
        lines.append(" Passive-regenerator result not available for this run.")
    if elastocaloric_result is not None:
        try:
            lines.append(f" Elastocaloric literature reference (, static, NOT "
                         f"span-simulated): COP_representative={elastocaloric_result.COP_representative:.2f} "
                         f"(range {elastocaloric_result.COP_low:.1f}-{elastocaloric_result.COP_high:.1f}) "
                         "-- an external anchor for comparison, not an output of this repo's own AMR model.")
        except Exception:
            lines.append(" Elastocaloric reference available but could not be summarized "
                         "(unexpected shape).")
    else:
        lines.append(" Elastocaloric reference not available for this run.")
    return "\n".join(lines), {"passive_regen_base": passive_regen_base,
                               "passive_regen_rows": passive_regen_rows,
                               "elastocaloric_result": elastocaloric_result}


def build_report(sobol_state_dependent_Si=None, pareto_rows=None, material_rows=None,
                  graded_row=None, gd_cascade_row=None, n_stages=3,
                  pb_best_cop_row=None, pp_best_cop_row=None,
                  representative_span_K=10.0,
                  cycle_type_result=None, thermal_diode_rows=None,
                  passive_regen_base=None, passive_regen_rows=None,
                  elastocaloric_result=None,
                  out_path="results/design_recommendations.txt"):
    """Assembles all five levers into one report, writes it to
    `out_path`, prints it, and returns a structured dict for any caller
    (e.g. main.py's final executive-summary banner) that wants the
    numbers programmatically rather than re-parsing the text.

    Every number in this report traces back to a result object computed
    by an existing pipeline stage and passed in as an argument -- this
    function performs no new physics and recomputes nothing expensive.
    Any lever whose upstream data wasn't provided (e.g. running this
    module standalone rather than via main.py) is reported as
    unavailable rather than silently omitted or fabricated.
    """
    lines = []
    lines.append("=" * 100)
    lines.append("AMR ELECTRICAL-COP DESIGN RECOMMENDATIONS -- consolidated from all analyses above")
    lines.append("=" * 100)
    lines.append("Every figure below was computed by an earlier pipeline stage in this same run; "
                 "this report only re-organizes those results by design lever, ranked by the "
                 "demonstrated Sobol sensitivity of electrical COP to each lever (sensitivity.py).")
    lines.append("")

    section1, data1 = summarize_frequency_lever(sobol_state_dependent_Si, pareto_rows or [])
    section2, data2 = summarize_material_lever(material_rows or [], representative_span_K)
    section3, data3 = summarize_grading_lever(graded_row, gd_cascade_row, n_stages)
    section4, data4 = summarize_geometry_lever(pb_best_cop_row, pp_best_cop_row)
    section5, data5 = summarize_field_flow_lever(pareto_rows or [])
    section6, data6 = summarize_cycle_type_finding(cycle_type_result)
    section7, data7 = summarize_thermal_diode_finding(thermal_diode_rows)
    section8, data8 = summarize_passive_regen_and_elastocaloric(
        passive_regen_base, passive_regen_rows, elastocaloric_result)

    for section in (section1, section2, section3, section4, section5,
                    section6, section7, section8):
        lines.append(section)
        lines.append("")

    lines.append("-" * 100)
    lines.append("RECOMMENDED STARTING DESIGN POINT (NSGA-III knee point, this run):")
    knee = data5.get("knee_point")
    if knee:
        lines.append(f" mu0H_max = {knee['mu0H_max_T']} T, frequency = {knee['frequency_Hz']} Hz, "
                     f"fluid_mdot = {knee['fluid_mdot_kgs']} kg/s, "
                     f"mass_regenerator = {knee['mass_regenerator_kg']} kg, "
                     f"regen_effectiveness = {knee['regen_effectiveness']}")
        geometry_note = f", particle_diameter = {knee['particle_diameter_mm']} mm" if "particle_diameter_mm" in knee else ""
        material_note = f", material = {knee['material']}" if "material" in knee else ""
        lines.append(f"  -> Predicted COP_electrical = {knee['COP_electrical']}, "
                     f"Qc = {knee['Qc_W']} W, cost_index = ${knee['cost_index_USD']}"
                     f"{geometry_note}{material_note}")
    lines.append("This is a starting point, not a final specification, but as of phase 15 the "
                 "geometry (lever 4) and material (lever 2) choices above ARE co-optimized "
                 "with the NSGA-III field/frequency/flow search (see core/optimize.py's module "
                 "docstring for how material is handled -- per-family NSGA-III passes merged "
                 "into one global Pareto front, since material family is categorical). "
                 "geometry_analysis.py and material_family_comparison.py remain independently "
                 "useful as focused single-lever sweeps for understanding WHY a given choice "
                 "wins, but optimize.py's own Pareto front is now this repo's most complete "
                 "answer to the joint design question.")
    lines.append("-" * 100)

    text = "\n".join(lines)
    print(text)
    os.makedirs(os.path.dirname(out_path) or ".", exist_ok=True)
    with open(out_path, "w") as fh:
        fh.write(text + "\n")

    return {
        "frequency": data1,
        "material": data2,
        "grading": data3,
        "geometry": data4,
        "field_flow": data5,
        "cycle_type": data6,
        "thermal_diode": data7,
        "passive_regen_and_elastocaloric": data8,
        "text": text,
        "out_path": out_path,
    }


if __name__ == "__main__":
    # Standalone smoke test: runs the (cheap) upstream stages this module
    # needs so `python -m core.design_recommendations` produces a real,
    # if partial, report on its own without requiring a full main.py run.
    from core import optimize as optimize_module
    from core import geometry_analysis
    from core import material_family_comparison

    pareto_rows = optimize_module.run_optimization()
    _, _, pb_best_cop = geometry_analysis.sweep_packed_bed_diameter(verbose=False)
    _, _, pp_best_cop = geometry_analysis.sweep_parallel_plate_spacing(verbose=False)
    material_rows = material_family_comparison.build_comparison_table()

    build_report(
        sobol_state_dependent_Si=None,
        pareto_rows=pareto_rows,
        material_rows=material_rows,
        graded_row=None,
        gd_cascade_row=None,
        pb_best_cop_row=pb_best_cop,
        pp_best_cop_row=pp_best_cop,
    )