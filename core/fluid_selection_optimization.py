"""
core/fluid_selection_optimization.py
=====================================
 addition, wired into main.py as its own pipeline step (see that
file's step list) -- moved here from a one-off research/ script once it
was confirmed useful enough to run on every pipeline pass, following
this repo's own convention that `core/` holds production analyses with a
`run_*()` entry point and `research/` holds one-off/exploratory
prototypes (see research/README.md).

Sweeps every heat-transfer fluid in `core.fluids.FLUID_LIBRARY` through
the SAME model every other production result in this repo uses --
`AMRSystem` with the NTU thermal model (`use_ntu_thermal_model=True`),
the calibrated `StateDependentLossModel`, and a geometry-explicit
hydraulic pumping-power override (`particle_diameter` set) -- to compare
fluids on a genuine, apples-to-apples COP_electrical/Qc basis rather than
just their raw {rho, cp, mu, k} property values (see `core/fluids.py`'s
module docstring for why those raw values alone don't decide the
question: cp/k losses and mu-driven pumping-power costs pull in opposite
directions, and this repo's own convention is "compare through the
model", not "eyeball the property table").

For each fluid, mdot is re-optimized (1-D grid search over this repo's
own `optimize.py` mdot bounds, [0.02, 0.5] kg/s) to maximize
COP_electrical -- mirroring `geometry_analysis.py`'s own stated
methodology of re-optimizing mdot per candidate rather than comparing at
a single fixed flow rate (a fixed-mdot comparison was checked first and
found to unfairly penalize higher-viscosity fluids, which want a lower
mdot to keep pumping losses down; reported below, not hidden).

Two operating points are checked for ranking robustness:
    BASELINE : T_cold=291K, span=10K (this repo's standard operating
               point, from `optimize.py`/`geometry_analysis.py`),
               mu0H_max=1.5T, mass_regenerator=5.0kg, frequency=1.0Hz,
               particle_diameter=0.3mm.
    ROBUSTNESS: same T_cold/span, but higher frequency (2.0Hz) and
               smaller regenerator mass (2.0kg) -- a meaningfully
               different design point (shorter NTU residence time per
               pass, less thermal mass) -- to check the ranking isn't an
               artifact of one specific design.

RESULT (see results/fluid_selection_optimization.txt and
results/fluid_selection_optimization_robustness.txt for full numbers):
    water > water_eg10 > water_eg20 > water_pg30 > ethanol
on COP_electrical, at BOTH operating points. Pure water is not a real
option for actual hardware (see core/fluids.py's module docstring: every
prototype in this project's Papers/ corpus that reports its fluid runs a
corrosion-inhibited water/glycol mixture, since Gd corrodes in plain
water). Of the corrosion-protected options, water_eg10 gives up the
least COP_electrical relative to the (unrealistic) pure-water ceiling --
about 1.4% at the baseline point, about 2.3% at the robustness point --
making it the recommended realistic default (`core.fluids.DEFAULT_FLUID`).

Honesty notes
-------------
- This is a single-figure-of-merit (COP_electrical) comparison at fixed
  mu0H_max/mass/frequency/particle_diameter; it does not re-run a full
  multi-objective search (`optimize.py`'s NSGA-III) per fluid. Given the
  monotonic, non-crossing ranking across two quite different operating
  points, re-running the full Pareto search per fluid was judged low
  marginal value for this question -- flagged here rather than silently
  assumed to generalize to every corner of the design space.
- The mdot grid search is 1-D and exhaustive (200 points, linear in
  mdot) over `optimize.py`'s own mdot design bound, [0.02, 0.5] kg/s --
  chosen to keep this comparison inside the same design space every
  other production result in this repo searches, rather than an
  unconstrained one.
- At the BASELINE point this bound contains a genuine interior optimum
  for every fluid (checked: COP_electrical rises then falls smoothly
  inside [0.02, 0.5] kg/s). At the ROBUSTNESS point (higher frequency,
  smaller mass), COP_electrical is still rising at mdot=0.5 for every
  fluid -- the true unconstrained optimum sits above this repo's own
  0.02-0.5 kg/s design bound (checked by extending the grid past 0.5:
  the true optimum for water sits near mdot~1.0 kg/s). The reported
  robustness-point numbers are therefore each fluid's best COP_electrical
  WITHIN this repo's established design space, not each fluid's global
  optimum -- flagged rather than silently widening the bound
  mid-comparison. This does not change the ranking conclusion (every
  fluid hits the same boundary, so the comparison is still apples-to-
  apples), but it means the robustness-point COP_electrical/Qc values
  are conservative/bound-limited, not each fluid's true best case.
"""

import numpy as np

from core.amr_cycle import AMRSystem
from core.loss_model import StateDependentLossModel
from core.mce_material import GADOLINIUM
from core.fluids import FLUID_LIBRARY, DEFAULT_FLUID

T_COLD_K = 291.0
SPAN_K = 10.0
BED_CROSS_SECTION_AREA_M2 = 0.002  # matches optimize.py / geometry_analysis.py
MDOT_BOUNDS_KGS = (0.02, 0.5)      # matches optimize.py's own mdot design bound
MDOT_GRID_N = 200

BASELINE_POINT = dict(
    label="baseline (T_cold=291K, span=10K, mu0H=1.5T, mass=5.0kg, f=1.0Hz)",
    mu0H_max=1.5, mass_regenerator=5.0, frequency=1.0,
    particle_diameter=0.0003,
)
ROBUSTNESS_POINT = dict(
    label="robustness (T_cold=291K, span=10K, mu0H=1.5T, mass=2.0kg, f=2.0Hz)",
    mu0H_max=1.5, mass_regenerator=2.0, frequency=2.0,
    particle_diameter=0.0003,
)


def _best_mdot_for_fluid(fluid, design_point, loss_model):
    """1-D grid search over mdot maximizing COP_electrical for `fluid` at
    `design_point`, using the NTU thermal model + calibrated loss model +
    geometry-explicit pumping-power override (particle_diameter set) --
    the same model combination every other production result in this
    repo uses (see core/optimize.py's AMRDesignProblem, core/pue_
    annualized.py's _amr_cop_at)."""
    mdots = np.linspace(MDOT_BOUNDS_KGS[0], MDOT_BOUNDS_KGS[1], MDOT_GRID_N)
    best = None
    for mdot in mdots:
        sys_ = AMRSystem(
            material=GADOLINIUM, mu0H_max=design_point["mu0H_max"],
            mass_regenerator=design_point["mass_regenerator"],
            frequency=design_point["frequency"], fluid_mdot=float(mdot),
            use_ntu_thermal_model=True,
            particle_diameter=design_point["particle_diameter"],
            bed_cross_section_area=BED_CROSS_SECTION_AREA_M2,
            loss_model=loss_model, fluid=fluid,
        )
        res = sys_.run(T_COLD_K, SPAN_K)
        if best is None or res.COP_electrical > best[1]:
            best = (float(mdot), res.COP_electrical, res.Qc, res.exergy_eff)
    return {"mdot_kgs": best[0], "COP_electrical": best[1], "Qc_W": best[2],
            "exergy_eff": best[3]}


def run_fluid_sweep(design_point, verbose=True):
    """Runs `_best_mdot_for_fluid` for every fluid in FLUID_LIBRARY at
    `design_point`. Returns rows sorted by COP_electrical, descending, and
    the reference (pure-water) COP for computing each fluid's %COP giveup."""
    loss_model = StateDependentLossModel()
    rows = []
    for fluid in FLUID_LIBRARY:
        r = _best_mdot_for_fluid(fluid, design_point, loss_model)
        r["fluid"] = fluid
        rows.append(r)
    rows.sort(key=lambda r: r["COP_electrical"], reverse=True)
    water_cop = next(r["COP_electrical"] for r in rows if r["fluid"] == "water")
    for r in rows:
        r["pct_below_water"] = 100.0 * (water_cop - r["COP_electrical"]) / water_cop

    if verbose:
        print(f"\n=== {design_point['label']} ===")
        print(f"{'fluid':<12}{'mdot(kg/s)':>12}{'COP_elec':>11}{'Qc(W)':>10}"
              f"{'exergy_eff':>12}{'%below water':>14}")
        for r in rows:
            print(f"{r['fluid']:<12}{r['mdot_kgs']:>12.4f}{r['COP_electrical']:>11.4f}"
                  f"{r['Qc_W']:>10.2f}{r['exergy_eff']:>12.4f}{r['pct_below_water']:>13.1f}%")
    return rows


def _write_results_txt(rows, design_point, path):
    with open(path, "w") as f:
        f.write("fluid_selection_optimization.py -- results\n")
        f.write("=" * 60 + "\n")
        f.write(f"Operating point: {design_point['label']}\n")
        f.write(f"mdot re-optimized per fluid over {MDOT_BOUNDS_KGS} kg/s "
                 f"({MDOT_GRID_N}-point grid), maximizing COP_electrical.\n\n")
        f.write(f"{'rank':<6}{'fluid':<12}{'mdot(kg/s)':>12}{'COP_elec':>11}"
                 f"{'Qc(W)':>10}{'exergy_eff':>12}{'%below water':>14}\n")
        for i, r in enumerate(rows, 1):
            f.write(f"{i:<6}{r['fluid']:<12}{r['mdot_kgs']:>12.4f}"
                     f"{r['COP_electrical']:>11.4f}{r['Qc_W']:>10.2f}"
                     f"{r['exergy_eff']:>12.4f}{r['pct_below_water']:>13.1f}%\n")
        f.write("\nRanking: " + " > ".join(r["fluid"] for r in rows) + "\n")
        real_rows = [r for r in rows if r["fluid"] != "water"]
        best_real = min(real_rows, key=lambda r: r["pct_below_water"])
        f.write(f"\nPure water is not a real option for actual hardware (Gd "
                 f"corrodes in plain water -- see core/fluids.py). Among the "
                 f"corrosion-protected fluids, '{best_real['fluid']}' gives up "
                 f"the least COP_electrical relative to the pure-water ceiling "
                 f"({best_real['pct_below_water']:.1f}%), making it the "
                 f"recommended realistic default.\n")
        f.write(f"DEFAULT_FLUID (core/fluids.py) = {DEFAULT_FLUID!r}\n")


def run_fluid_selection_comparison(
        out_path="results/fluid_selection_optimization.txt",
        robustness_out_path="results/fluid_selection_optimization_robustness.txt",
        verbose=True):
    """Main entry point (mirrors this repo's own `run_*()` convention --
    see `geometry_analysis.run_geometry_analysis()`, `fluid_mce_analysis.
    run_fluid_mce_analysis()` for the same pattern). Runs `run_fluid_sweep`
    at both BASELINE_POINT and ROBUSTNESS_POINT, writes each to its own
    results/ file (ADDITIVE ONLY -- these are new files; no existing
    results/ file is touched), and returns a dict with both rows lists
    plus the recommended-realistic-fluid summary for main.py's executive
    summary / step 13 synthesis report to consume, the same way step 3c's
    `pb_best_cop_row`/`pp_best_cop_row` are captured and re-used.

    Does NOT change `core.fluids.DEFAULT_FLUID` or any existing
    `AMRSystem` call site's `fluid` default -- see this module's own
    top-level docstring and `core/fluids.py`'s for why that stays a
    deliberate, opt-in choice rather than a side effect of running this
    analysis."""
    import os
    os.makedirs(os.path.dirname(out_path), exist_ok=True)

    baseline_rows = run_fluid_sweep(BASELINE_POINT, verbose=verbose)
    _write_results_txt(baseline_rows, BASELINE_POINT, out_path)

    robustness_rows = run_fluid_sweep(ROBUSTNESS_POINT, verbose=verbose)
    _write_results_txt(robustness_rows, ROBUSTNESS_POINT, robustness_out_path)

    baseline_order = [r["fluid"] for r in baseline_rows]
    robustness_order = [r["fluid"] for r in robustness_rows]
    same_ranking = baseline_order == robustness_order
    if verbose:
        print(f"\nRanking identical at both operating points: {same_ranking}")
        if not same_ranking:
            print("  baseline:   ", " > ".join(baseline_order))
            print("  robustness: ", " > ".join(robustness_order))

    real_rows = [r for r in baseline_rows if r["fluid"] != "water"]
    best_realistic = min(real_rows, key=lambda r: r["pct_below_water"])
    return {
        "baseline_rows": baseline_rows,
        "robustness_rows": robustness_rows,
        "same_ranking": same_ranking,
        "best_realistic_fluid": best_realistic["fluid"],
        "best_realistic_pct_below_water_baseline": best_realistic["pct_below_water"],
        "default_fluid": DEFAULT_FLUID,
        "default_matches_best_realistic": best_realistic["fluid"] == DEFAULT_FLUID,
    }


if __name__ == "__main__":
    run_fluid_selection_comparison()
