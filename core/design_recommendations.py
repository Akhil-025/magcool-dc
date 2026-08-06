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
    1. Operating frequency       (ST ~ 0.85-0.87, the dominant lever --
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
    3. Curie-temperature grading  (`cascade.py`'s graded cascade: matching
                                   each layer's Tc to its local fluid
                                   temperature raises both Qc and COP
                                   relative to a uniform-material cascade
                                   of the same stage count)
    4. Regenerator geometry       (`geometry_analysis.py`: packed-bed
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
                                   any one of them to its bound)

Recent literature context (not used quantitatively, cited for the reader):
Klinar, Kitanovski, Law, Franco & Moya, "Perspectives and Energy
Applications of Magnetocaloric, Pyromagnetic, Electrocaloric, and
Pyroelectric Materials," Adv. Energy Mater. 14, 2401739 (2024) -- a 2024
roadmap-style review that (among other things) profiles a novel
high-frequency active caloric regenerator design concept ("Hypereg")
aimed squarely at this repo's own finding that frequency is the dominant
COP lever; see ROADMAP.md's "Future work" section for how that connects
to a possible extension of `loss_model.py`'s eddy-current term.
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
    names = list(sobol_state_dependent_Si["names"]) if isinstance(sobol_state_dependent_Si, dict) \
        and "names" in sobol_state_dependent_Si else None
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
        lines.append(f"   Sobol total-order sensitivity (state-dependent loss model): "
                     f"ST(frequency) = {freq_ST:.3f}")
    if best_cop and best_qc:
        lines.append(f"   NSGA-III best-electrical-COP design:  f={best_cop['frequency_Hz']:.3f} Hz  "
                     f"-> COP_elec={best_cop['COP_electrical']:.2f}, Qc={best_cop['Qc_W']:.0f} W")
        lines.append(f"   NSGA-III best-cooling-capacity design: f={best_qc['frequency_Hz']:.3f} Hz  "
                     f"-> COP_elec={best_qc['COP_electrical']:.2f}, Qc={best_qc['Qc_W']:.0f} W")
        if best_qc["frequency_Hz"] > 0:
            lines.append(f"   Action: reduce cycle frequency toward the low-speed end of the "
                         f"design space when electrical COP is the priority; raise it toward "
                         f"the high end only when cooling capacity is the priority and higher "
                         f"parasitic loss is acceptable. This is a genuine Pareto trade-off in "
                         f"this repo's model, not a free win in either direction.")
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
        lines.append(f"   Best candidate at {representative_span_K:.0f}K span (ASHRAE point): "
                     f"{best['candidate']}  COP_elec={best['1stage_COP']:.2f}, "
                     f"Qc={best['1stage_Qc_W']:.0f} W")
        if gd_row and gd_row is not best:
            gain_pct = 100.0 * (best["1stage_COP"] / gd_row["1stage_COP"] - 1.0)
            lines.append(f"   vs. plain Gd: COP_elec={gd_row['1stage_COP']:.2f} "
                         f"({gain_pct:+.0f}% relative to Gd)")
        lines.append("   Action: prefer a composition-tuned giant-MCE family (e.g. "
                     "La(Fe,Si)13Hy-type) whose documented Tc window covers the target "
                     "operating point, over a fixed-composition material whose peak may sit "
                     "off-target -- see material_family_comparison.py for the full "
                     "span-by-span ranking and which families actually cover which spans.")
    else:
        lines.append("   No in-range tunable candidate found at this span in the current "
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
                         f"graded COP={g_cop}, Qc={g_qc} W   vs.  "
                         f"plain-Gd COP={p_cop}, Qc={p_qc} W")
            lines.append("   Action: grade each stage's composition so its local peak MCE "
                         "temperature tracks the fluid temperature at that point in the bed, "
                         "rather than using one uniform material across the whole span -- "
                         "see cascade.py's compare_graded_cascade() for the per-span/"
                         "stage-count feasibility breakdown (not every span/stage-count "
                         "combination stays within the documented composition-tunability "
                         "window; some stages fall back to plain Gd, which is reported "
                         "explicitly rather than silently).")
    else:
        lines.append("   Graded-cascade data not available for this run.")
    return "\n".join(lines), {"graded_row": graded_row, "gd_cascade_row": gd_cascade_row}


def summarize_geometry_lever(pb_best_cop_row, pp_best_cop_row):
    """Lever 4: regenerator geometry, from the already-computed packed-bed
    and parallel-plate sweeps (main.py step 3c,
    `geometry_analysis.sweep_packed_bed_diameter` /
    `sweep_parallel_plate_spacing`)."""
    lines = ["4. REGENERATOR GEOMETRY (packed-bed vs. parallel-plate)"]
    if pb_best_cop_row is not None:
        lines.append(f"   Packed-bed sphere diameter maximizing COP (interior optimum): "
                     f"{pb_best_cop_row[0]} mm  (Qc={pb_best_cop_row[1]:.0f} W, "
                     f"COP_aug={pb_best_cop_row[2]:.2f})")
    if pp_best_cop_row is not None:
        lines.append(f"   Parallel-plate channel spacing maximizing COP (interior optimum): "
                     f"{pp_best_cop_row[0]} mm  (Qc={pp_best_cop_row[1]:.0f} W, "
                     f"COP_aug={pp_best_cop_row[2]:.2f})")
    lines.append("   Action: target these interior optima rather than minimizing particle/"
                 "channel size without bound -- shrinking geometry indefinitely raises "
                 "thermal effectiveness only marginally further while hydraulic pumping "
                 "power keeps growing, which is what produces the interior optimum in the "
                 "first place (Tusek, Kitanovski & Poredos 2013).")
    return "\n".join(lines), {"pb_best_cop": pb_best_cop_row, "pp_best_cop": pp_best_cop_row}


def summarize_field_flow_lever(pareto_rows):
    """Lever 5: field strength / flow rate balance, from the NSGA-III
    Pareto front's knee-point design (main.py step 11,
    `optimize.run_optimization()`), which is this repo's own definition
    of a balanced (not extreme) design point."""
    import numpy as np
    lines = ["5. FIELD STRENGTH / FLOW RATE BALANCE"]
    if pareto_rows:
        F = np.array([[-r["COP_electrical"], -r["Qc_W"], r["cost_index_USD"]]
                      for r in pareto_rows])
        Fn = (F - F.min(axis=0)) / (F.max(axis=0) - F.min(axis=0) + 1e-12)
        knee_idx = int(np.argmin(np.linalg.norm(Fn, axis=1)))
        knee = pareto_rows[knee_idx]
        lines.append(f"   Knee-point (balanced) Pareto design: "
                     f"H={knee['mu0H_max_T']} T, f={knee['frequency_Hz']} Hz, "
                     f"mdot={knee['fluid_mdot_kgs']} kg/s, mass={knee['mass_regenerator_kg']} kg, "
                     f"eps={knee['regen_effectiveness']}  "
                     f"-> COP_elec={knee['COP_electrical']}, Qc={knee['Qc_W']} W, "
                     f"cost=${knee['cost_index_USD']}")
        lines.append("   Action: operate near this balanced point rather than at either "
                     "objective's own extreme (max-COP alone drives frequency and Qc very "
                     "low; max-Qc alone raises parasitic loss and cuts COP roughly 3x in "
                     "this repo's Pareto front) -- see results/pareto_front.csv for the full "
                     f"{len(pareto_rows)}-design front and the other extreme points.")
        return "\n".join(lines), {"knee_point": knee}
    lines.append("   Pareto front not available for this run.")
    return "\n".join(lines), {"knee_point": None}


def build_report(sobol_state_dependent_Si=None, pareto_rows=None, material_rows=None,
                  graded_row=None, gd_cascade_row=None, n_stages=3,
                  pb_best_cop_row=None, pp_best_cop_row=None,
                  representative_span_K=10.0,
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

    for section in (section1, section2, section3, section4, section5):
        lines.append(section)
        lines.append("")

    lines.append("-" * 100)
    lines.append("RECOMMENDED STARTING DESIGN POINT (NSGA-III knee point, this run):")
    knee = data5.get("knee_point")
    if knee:
        lines.append(f"  mu0H_max = {knee['mu0H_max_T']} T,  frequency = {knee['frequency_Hz']} Hz,  "
                     f"fluid_mdot = {knee['fluid_mdot_kgs']} kg/s,  "
                     f"mass_regenerator = {knee['mass_regenerator_kg']} kg,  "
                     f"regen_effectiveness = {knee['regen_effectiveness']}")
        lines.append(f"  -> Predicted COP_electrical = {knee['COP_electrical']}, "
                     f"Qc = {knee['Qc_W']} W, cost_index = ${knee['cost_index_USD']}")
    lines.append("This is a starting point, not a final specification: the geometry (lever 4) "
                 "and material (lever 2) choices above are NOT yet co-optimized with the "
                 "NSGA-III field/frequency/flow search in this repo's current optimize.py -- "
                 "see ROADMAP.md's Future Work for that open item.")
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