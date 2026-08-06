"""
material_family_comparison.py
==============================
Answers "which material is actually best for this application?" directly,
instead of leaving it implicit across several separate analyses.

giant_mce_analysis.py already compares Gd against ONE fixed-composition
giant-MCE material (Gd5Si2Ge2, Tc=276K) and shows it collapses to ~0 in the
ASHRAE range because its Tc is fixed and off-target. But three families in
this codebase are composition-TUNABLE (core/cascade.py's GD_FAMILY,
LAFESIH_FAMILY, MNFEPSI_FAMILY, each built on a first_order_mce.py Landau
material + a composition_tuned_material()-style factory) -- so the fair
question is not "does the one fixed composition work here" but "can EACH
family be tuned to a composition that works here, and if so, how does it
actually perform against Gd and against each other".

This module runs all five candidates --
    1. Gd                      (fixed, core.mce_material.GADOLINIUM)
    2. Gd5Si2Ge2 (fixed comp.) (core.first_order_mce.GD5SI2GE2_FIRST_ORDER)
    3. Gd5(SixGe1-x)4(-Ga)     (GD_FAMILY, composition-tuned per span)
    4. La(Fe,Si)13Hy           (LAFESIH_FAMILY, composition-tuned per span)
    5. (Mn,Fe)2(P,Si)          (MNFEPSI_FAMILY, composition-tuned per span)
through the SAME fixed ASHRAE operating point(s) and the SAME cascade logic
(core.cascade.run_cascade, 1-4 stage), and outputs a ranked table of Qc,
COP, and whether each tunable family's own documented Tc window
(GIANT_MCE_TC_MIN_K/_MAX_K etc. in first_order_mce.py) actually covers the
composition needed to hit that point.

Each tunable family is composition-tuned so its OWN peak DeltaT_ad lands at
the operating point's midpoint temperature T_mid = T_cold + span/2 -- NOT
at its raw Tc parameter, since giant_mce_analysis.py's
landau_peak_offset_K() documents a real, independently-confirmed +10-11K
peak-vs-Tc offset for this Landau formulation (peak sits ABOVE Tc). The
Tc needed to hit a given T_mid is solved with the same bracketed
root-finder cascade.py's Curie-graded cascade already relies on
(cascade._target_composition_for_peak), so this module adds no new
numerics -- it only reuses that machinery outside the graded-cascade
context, as a single-composition-per-span comparison instead of a
6-stage-graded one.

If the solved Tc falls outside a family's documented window, that family
is reported as NOT COVERING the ASHRAE range at that span and falls back
to its fallback_material (plain Gd) for the Qc/COP columns -- the same
fallback behavior compare_graded_cascade() already uses per-stage.
"""

import csv

from core.mce_material import GADOLINIUM
from core.first_order_mce import GD5SI2GE2_FIRST_ORDER
from core.cascade import (
    run_cascade, GD_FAMILY, LAFESIH_FAMILY, MNFEPSI_FAMILY,
    _target_composition_for_peak,
)
from core.baseline_cooling import vapor_compression_cop, liquid_cooling_cop

MU0H_MAX = 2.0
MASS_PER_STAGE = 5.0
STAGE_COUNTS = (1, 2, 3, 4)
T_COLD_C = 18.0
T_COLD_K = T_COLD_C + 273.15
SPANS_K = (5.0, 10.0, 15.0, 20.0)   # ASHRAE-range representative sweep
REPRESENTATIVE_SPAN_K = 10.0        # matches main.py's steps 5/6/9/11 anchor

TUNABLE_FAMILIES = (GD_FAMILY, LAFESIH_FAMILY, MNFEPSI_FAMILY)


def _tuned_candidate(family, T_mid_K, mu0H_max=MU0H_MAX):
    """Returns (material, tc_used_K, in_range: bool) for `family` tuned so
    its own peak DeltaT_ad lands at T_mid_K. Falls back to
    family.fallback_material (plain Gd) if the required Tc sits outside
    the family's documented tunability window."""
    tc = _target_composition_for_peak(T_mid_K, mu0H_max, family)
    in_range = family.tc_min <= tc <= family.tc_max
    if not in_range:
        return family.fallback_material, tc, False
    return family.tuned_fn(tc), tc, True


def _eval_cascade(material, T_cold_K, span_K, mu0H_max=MU0H_MAX,
                   mass_per_stage=MASS_PER_STAGE, stage_counts=STAGE_COUNTS):
    """Runs core.cascade.run_cascade for each stage count -- the same
    cascade logic run_cascade_comparison() (main.py step 7) and
    giant_mce_analysis.py's single-stage AMRSystem calls both build on,
    just applied here across all five candidates uniformly."""
    out = {}
    for n in stage_counts:
        res = run_cascade(T_cold_K, span_K, n, material=material,
                           mu0H_max=mu0H_max, mass_per_stage=mass_per_stage)
        out[f"{n}stage_COP"] = round(res["COP_cascade"], 3) if res["feasible"] else None
        out[f"{n}stage_Qc_W"] = round(res["Qc_W"], 1) if res["feasible"] else None
    return out


def build_comparison_table(spans_K=SPANS_K, T_cold_K=T_COLD_K):
    """Builds one row per (candidate, span). Fixed-composition candidates
    (Gd, Gd5Si2Ge2) use the same material at every span; tunable families
    are re-tuned per span to their own best composition for that span's
    T_mid, so each family is evaluated at its own fairest shot rather than
    penalized for a single fixed composition (that comparison, for
    Gd5Si2Ge2 specifically, already exists and is intentionally kept
    separate -- see cascade_comparison_giant_mce.csv / fig20)."""
    rows = []
    for span in spans_K:
        T_mid = T_cold_K + span / 2.0

        rows.append({"candidate": "Gd (fixed)", "span_K": span, "T_mid_K": round(T_mid, 1),
                     "Tc_used_K": GADOLINIUM.Tc, "tc_window": "n/a (not tunable)",
                     "in_range": True,
                     **_eval_cascade(GADOLINIUM, T_cold_K, span)})

        rows.append({"candidate": "Gd5Si2Ge2 (fixed comp.)", "span_K": span, "T_mid_K": round(T_mid, 1),
                     "Tc_used_K": GD5SI2GE2_FIRST_ORDER.Tc, "tc_window": "n/a (not tunable)",
                     "in_range": True,
                     **_eval_cascade(GD5SI2GE2_FIRST_ORDER, T_cold_K, span)})

        for family in TUNABLE_FAMILIES:
            material, tc_used, in_range = _tuned_candidate(family, T_mid)
            rows.append({
                "candidate": f"{family.name} (tuned)", "span_K": span, "T_mid_K": round(T_mid, 1),
                "Tc_used_K": round(tc_used, 1),
                "tc_window": f"{family.tc_min:.1f}-{family.tc_max:.1f}K",
                "in_range": in_range,
                **_eval_cascade(material, T_cold_K, span),
            })
    return rows


def run_analysis(out_csv="results/material_family_comparison.csv",
                  out_txt="results/material_family_comparison.txt"):
    rows = build_comparison_table()

    fieldnames = ["candidate", "span_K", "T_mid_K", "Tc_used_K", "tc_window", "in_range"]
    fieldnames += [f"{n}stage_COP" for n in STAGE_COUNTS] + [f"{n}stage_Qc_W" for n in STAGE_COUNTS]
    with open(out_csv, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)

    lines = []
    lines.append("Four-way material family comparison at the ASHRAE operating point")
    lines.append(f"(T_cold={T_COLD_C:.0f}C={T_COLD_K:.2f}K, spans={list(SPANS_K)}K, "
                 f"mu0H={MU0H_MAX:.1f}T, {MASS_PER_STAGE:.0f}kg/stage, 1-4 stage cascade)")
    lines.append("=" * 100)
    lines.append("Tunable families (GD_FAMILY, LAFESIH_FAMILY, MNFEPSI_FAMILY) are re-tuned per span")
    lines.append("so their OWN peak lands at that span's T_mid = T_cold + span/2 (see module docstring")
    lines.append("for the +10-11K peak-vs-Tc offset this accounts for). Falls back to plain Gd where")
    lines.append("the required Tc sits outside the family's documented tunability window.")
    lines.append("")

    header = (f"{'candidate':<30}{'span_K':>7}{'T_mid_K':>9}{'Tc_used_K':>11}  "
              f"{'tc_window':<19}{'in_range':>9}{'1stg_COP':>10}{'1stg_Qc_W':>11}")
    lines.append(header)
    for r in rows:
        lines.append(f"{r['candidate']:<30}{r['span_K']:>7.0f}{r['T_mid_K']:>9.1f}"
                     f"{r['Tc_used_K']:>11.1f}  {r['tc_window']:<19}{str(r['in_range']):>9}"
                     f"{'--' if r['1stage_COP'] is None else r['1stage_COP']:>10}"
                     f"{'--' if r['1stage_Qc_W'] is None else r['1stage_Qc_W']:>11}")
    lines.append("")

    # Ranked summary at the representative span (matches steps 5/6/9/11's anchor).
    # Fallback candidates (in_range=False) are literally running plain Gd under a
    # different label -- rank only the genuinely distinct candidates, and list
    # fallbacks separately so they don't read as an independent result next to
    # the real Gd row they're identical to.
    rep_rows = [r for r in rows if r["span_K"] == REPRESENTATIVE_SPAN_K]
    rankable = [r for r in rep_rows if r["in_range"]]
    fallen_back = [r for r in rep_rows if not r["in_range"]]
    ranked = sorted(
        [r for r in rankable if r["1stage_COP"] is not None],
        key=lambda r: r["1stage_COP"], reverse=True,
    )
    infeasible = [r for r in rankable if r["1stage_COP"] is None]

    lines.append(f"RANKED (1-stage COP_electrical, at representative span={REPRESENTATIVE_SPAN_K:.0f}K):")
    for i, r in enumerate(ranked, 1):
        lines.append(f"  {i}. {r['candidate']:<26} COP={r['1stage_COP']:.2f}  Qc={r['1stage_Qc_W']:.1f}W")
    for r in infeasible:
        lines.append(f"  --  {r['candidate']:<26} INFEASIBLE at this span")
    if fallen_back:
        lines.append("")
        lines.append("  Not independently ranked (Tc window doesn't cover this point -> fell back")
        lines.append("  to plain Gd, so these are identical to the Gd row above, not a distinct result):")
        for r in fallen_back:
            lines.append(f"      {r['candidate']:<26} == Gd (fixed) at this span")
    lines.append("")

    T_hot_rep = T_COLD_K + REPRESENTATIVE_SPAN_K
    vcc = vapor_compression_cop(T_COLD_K, T_hot_rep)
    liq = liquid_cooling_cop(T_COLD_K, T_hot_rep)
    lines.append(f"For reference, baselines at this operating point: "
                 f"VCC COP={vcc.COP:.2f}, Liquid COP={liq.COP:.2f}")
    lines.append("")
    lines.append("CONCLUSION:")
    if ranked:
        best = ranked[0]
        lines.append(f"Best-performing AMR candidate at the representative point: {best['candidate']} "
                     f"(1-stage COP_electrical={best['1stage_COP']:.2f}, Qc={best['1stage_Qc_W']:.1f}W).")
    lines.append("This does not by itself establish that any AMR candidate beats vapor-compression or")
    lines.append("liquid cooling on COP -- see the baselines above and comparison_table.csv/step 4 for")
    lines.append("that separate, already-established conclusion. What this table adds is a fair,")
    lines.append("apples-to-apples ranking AMONG the giant-MCE material options themselves, including")
    lines.append("whether each family's documented tunability window can even reach the ASHRAE point.")

    with open(out_txt, "w") as f:
        f.write("\n".join(lines) + "\n")
    print("\n".join(lines))
    print(f"\nWrote {out_csv}")
    print(f"Wrote {out_txt}")
    return rows


if __name__ == "__main__":
    run_analysis()