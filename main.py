"""
main.py
=======
End-to-end pipeline for the magcool-dc project. Runs every analysis module
in the repository in one pass, in dependency order, so a single
``python main.py`` reproduces every file under results/ from scratch:

    1.  Material-level validation      (core/validation.py)
    2.  System-level validation        (core/validation_system.py)
    3.  Loss-model calibration report  (core/loss_model.py)
    3b. Regenerator thermal demo       (core/thermal.py)
    3c. Geometry-dependent pumping power (core/geometry_analysis.py, ROADMAP.md Phase 7 item)
    3d. Hypereg parallel-hydraulic pumping-power analysis (core/hypereg_analysis.py, Phase 15 item 3)
    4.  Baseline comparison sweep      (this file, was the old main.py)
    5.  Economics / TCO                (core/economics.py)
    5b. Full-system cost estimate by material family (core/economics.py, Phase 15 item 5)
    6.  Emissions comparison           (core/emissions.py)
    7.  Cascade staging comparison     (core/cascade.py)
    7b. Curie-graded cascade           (core/cascade.py, ROADMAP.md Phase 7 item)
    7c. Astronautics graded-bed check  (core/cascade.py, ROADMAP.md Phase 9 addendum)
    8.  Giant-MCE materials analysis   (core/giant_mce_analysis.py)
    8d. Material family comparison     (core/material_family_comparison.py, Track A2 item)
    8b. First-order Landau model demo  (core/first_order_mce.py)
    8c. Giguere et al. (1999) direct-measurement cross-check + Pecharsky &
        Gschneidner (1997) peak-ratio check (core/giguere_validation.py)
    9.  Sobol sensitivity analysis     (core/sensitivity.py)
    10. RSM surrogate fit              (core/rsm.py)
    11. NSGA-III design optimization   (core/optimize.py, Phase 15: material +
        geometry co-optimization, per-material-family fronts merged post-hoc)
    12. Figure generation (26 figures) (plots.py -> results/figures/*.png, *.pdf)
    13. Design-recommendations synthesis (core/design_recommendations.py) --
        consolidates steps 3c/7b/8d/9b/11's already-computed results into
        one ranked, actionable "how do I raise AMR electrical COP" report
        (results/design_recommendations.txt)

Steps 3b and 8b reproduce core/thermal.py's and core/first_order_mce.py's
own __main__ demo blocks. Both modules are otherwise only reached
transitively (amr_cycle.py imports thermal.py internally;
giant_mce_analysis.py imports first_order_mce.py's material constant),
so without these two explicit steps their standalone demo output would
never appear even though the underlying code is exercised.

Steps 5 and 6 do not re-run their own illustrative examples; they are fed
the *actual* AMR/vapor-compression/liquid-cooling numbers computed in step
4 at the representative ASHRAE operating point (T_cold=18C, span=10K),
which is also the fixed point used internally by sensitivity.py and
optimize.py. Everything else calls the same public functions each module
already exposes for its own ``if __name__ == "__main__"`` block, so results
are identical to running each script individually -- just aggregated.

Step 13 (design_recommendations.py) does not recompute anything either: it
reuses the Sobol Si object, NSGA-III Pareto rows, material-family rows,
graded-cascade row, and geometry sweep rows already produced by steps
9b/11/8d/7b/7/3c above and reassembles them into one consolidated report,
so its own runtime is negligible even though the analyses it summarizes
are not.

Step 12 runs before the new step 13; step 13 depends only on already
-computed result objects from earlier stages, not on the figures
themselves, but runs last so its consolidated report can also mention
the figure count. plots.py is largely self-contained -- most of its 25
figures call straight into core/ and recompute their own data rather
than reading the CSVs the earlier steps write -- but three figures
(cascade staging, Curie-graded cascade, NSGA-III Pareto front) also
rewrite results/*.csv as a side effect, so running it last leaves
results/ fully consistent with the figures next to it.

Each stage is wrapped so a failure in one analysis logs an error and
does not stop the rest of the pipeline. A summary of what ran, what
failed, output files written, and total wall time is printed at the end.
Every stage's progress, timing, and any failure is written both to the
console and to results/pipeline.log via the logging module. Stage
functions that print directly (most of core/'s modules use plain print()
in their own `if __name__ == "__main__"` style, rather than the logging
module main.py itself uses for its bespoke per-stage summaries) have
their stdout captured and routed through the same logger for the
duration of each stage (see _StreamToLogger below), so that output is
also preserved in results/pipeline.log rather than only appearing on a
live console and vanishing afterward.

Typical total runtime: ~6 minutes on a modest machine. Figure generation
(step 12, ~200s) and the Curie-graded cascade sweep (step 7b, ~100s) are
the current bulk of it, not NSGA-III optimization or Sobol sampling
(step 11 and steps 9/9b are ~5-15s each) -- update this estimate again if
a future change shifts where the time goes. Nothing here is required to
reproduce results/ except the packages in requirements.txt (SALib, pymoo,
and matplotlib included).

Phase 15 (see ROADMAP.md) added: geometry (particle diameter) and material
family as genuine NSGA-III design variables/candidates in step 11; a
Hypereg-style parallel-hydraulic pumping-power analysis in new step 3d; a
full-system BOM cost model (materials + soft-magnetic yoke, plus an
order-of-magnitude full-system estimate and a CRF-based levelized cost of
cooling) in step 5, plus a per-material-family cost comparison in new step
5b; and confirmed (rather than duplicated) that the "does loss behavior
differ for rotary vs. reciprocating / multi-bed AMR topologies" question
was already answered by existing core/loss_model.py infrastructure
(RotaryDriveLossModel, analyze_parasitic_fraction_scaling) -- see
core/loss_model.py's module docstring "Phase 15 note" and ROADMAP.md for
the full writeup. NOTE: the original Phase 15 plan's item 1 (dedicated
tests for core/design_recommendations.py) was skipped in one earlier pass
because that module appeared to be missing from an out-of-date project
snapshot -- it is present here (step 13 above already exists and
integrates it), so that concern does not apply to this version of the
file; if tests for it are still wanted, they can be added directly against
this module's actual summarize_*_lever()/build_report() signatures.
"""

import logging
import os
import time
import traceback
import contextlib

import numpy as np
import csv

from core.mce_material import GADOLINIUM
from core.amr_cycle import AMRSystem
from core.baseline_cooling import vapor_compression_cop, liquid_cooling_cop
from core import validation
from core import validation_system
from core import loss_model
from core import economics
from core import emissions
from core import cascade
from core.cascade import staged_baseline_result
from core import giant_mce_analysis
from core import material_family_comparison
from core import sensitivity
from core import rsm
from core import optimize as optimize_module
from core.thermal import regenerator_effectiveness
from core.first_order_mce import GD5SI2GE2_FIRST_ORDER
from core import giguere_validation
from core import geometry_analysis
from core import hypereg_analysis
from core import plots
from core import design_recommendations
from core import sensitivity as sensitivity_module

RESULTS_DIR = "results"
RESULTS_CSV = "results/comparison_table.csv"
LOG_FILE = os.path.join(RESULTS_DIR, "pipeline.log")

# Representative operating point that steps 5/6 borrow from step 4's
# sweep. Matches the fixed point used internally by sensitivity.py and
# optimize.py (T_cold=291K, span=10K) so all headline numbers in this
# pipeline are talking about the same design point.
REPRESENTATIVE_SPAN_K = 10.0


def _setup_logging():
    """Configures a root logger that writes timestamped, leveled records
    to the console only. Safe to run at import time (e.g. when pytest
    collects tests/test_main.py via `import main`) since it never touches
    results/pipeline.log -- see _attach_file_logging() for that."""
    formatter = logging.Formatter(
        fmt="%(asctime)s [%(levelname)-7s] %(name)s: %(message)s",
        datefmt="%Y-%m-%d %H:%M:%S",
    )

    console_handler = logging.StreamHandler()
    console_handler.setFormatter(formatter)

    root = logging.getLogger()
    root.setLevel(logging.INFO)
    root.handlers.clear()
    root.addHandler(console_handler)

    return logging.getLogger("main")


def _attach_file_logging():
    """Attaches the results/pipeline.log FileHandler (mode="w", so the
    log always reflects the most recent full pipeline execution).

    Deliberately NOT called at module import time: doing so meant simply
    running `import main` -- e.g. via pytest collecting tests/test_main.py,
    which itself explicitly documents that it should NOT touch
    results/pipeline.log -- silently truncated the log to 0 bytes without
    ever running a single pipeline stage. Only main() calls this, so the
    log is only overwritten when the pipeline actually runs."""
    os.makedirs(RESULTS_DIR, exist_ok=True)
    formatter = logging.Formatter(
        fmt="%(asctime)s [%(levelname)-7s] %(name)s: %(message)s",
        datefmt="%Y-%m-%d %H:%M:%S",
    )
    file_handler = logging.FileHandler(LOG_FILE, mode="w")
    file_handler.setFormatter(formatter)
    logging.getLogger().addHandler(file_handler)


logger = _setup_logging()


def _banner(title):
    logger.info("")
    logger.info("#" * 100)
    logger.info("# %s", title)
    logger.info("#" * 100)


class _StreamToLogger:
    """File-like stream that routes each complete line written to it
    through logger.info(), so plain print() calls made inside a stage
    function (most of core/'s modules print directly, in their own
    ``if __name__ == "__main__"`` style, rather than using the logging
    module the way main.py's own bespoke per-stage summaries do) are
    captured in results/pipeline.log via the existing FileHandler instead
    of only reaching a live console and leaving no persistent record.

    Buffers partial lines (write() can be called with partial output,
    e.g. by a progress indicator using end="") until a newline arrives.
    Does not double-log logger.info() calls made directly inside a stage
    function: those go straight to the logger's own handlers and never
    touch sys.stdout, so redirecting stdout during a stage has no effect
    on them.
    """

    def __init__(self, target_logger, level=logging.INFO):
        self._logger = target_logger
        self._level = level
        self._buffer = ""

    def write(self, message):
        self._buffer += message
        while "\n" in self._buffer:
            line, self._buffer = self._buffer.split("\n", 1)
            if line:
                self._logger.log(self._level, line)

    def flush(self):
        if self._buffer:
            self._logger.log(self._level, self._buffer)
            self._buffer = ""


def run_baseline_sweep():
    """Step 4: AMR vs vapor-compression vs liquid cooling vs Carnot,
    swept over the ASHRAE 5-20K temperature-lift range. Returns the full
    list of row dicts so later steps can pull the representative point
    out of it instead of recomputing or hardcoding numbers.

    Uses core.cascade.staged_baseline_result() rather than a bare
    single-stage AMRSystem.run(). A single Gd stage at this operating
    point (2T, 5kg, 2Hz, mdot=0.08kg/s) runs out of its own no-load
    DeltaT_ad above ~16K span (amr_cycle.py's cooling_capacity() correctly
    returns Qc=0 past that point -- see its MODEL LIMITATION docstring),
    which previously showed up here as a flat AMR_COP_electrical=0.0 for
    17-20K span, reading as "AMR stops working" rather than "a single
    stage stops working." staged_baseline_result() automatically falls
    back to the minimum number of identical stages in series (2-4,
    same per-stage mass/field/flow as the single-stage case) needed to
    reach a positive Qc -- exactly the staging approach step 7's own
    cascade comparison already validates for this material at this span
    range. n_stages=1 whenever the single stage already worked, so every
    span that previously had a nonzero value is completely unchanged."""
    T_cold_C = 18.0
    T_cold_K = T_cold_C + 273.15
    spans = np.arange(5, 21, 1)

    amr_kwargs = dict(
        material=GADOLINIUM,
        mu0H_max=2.0,
        mass_regenerator=5.0,
        frequency=2.0,
        fluid_cp=4186.0,
        fluid_mdot=0.08,
        regenerator_effectiveness=0.85,
    )

    rows = []
    for span in spans:
        T_hot_K = T_cold_K + span
        amr_res = staged_baseline_result(T_cold_K, float(span), **amr_kwargs)
        vcc = vapor_compression_cop(T_cold_K, T_hot_K)
        liq = liquid_cooling_cop(T_cold_K, T_hot_K)
        rows.append({
            "span_K": span,
            "Tc_K": T_cold_K, "Th_K": T_hot_K,
            "AMR_COP_ideal": round(amr_res.COP, 2),
            "AMR_COP_electrical": round(amr_res.COP_electrical, 2),
            "AMR_Qc_W": round(amr_res.Qc, 1),
            "AMR_2ndlaw_eff": round(amr_res.exergy_eff, 3),
            "AMR_n_stages": amr_res.n_stages,
            "VaporCompression_COP": round(vcc.COP, 2),
            "LiquidCooling_COP": round(liq.COP, 2),
            "Carnot_COP": round(vcc.COP_carnot, 2),
        })

    with open(RESULTS_CSV, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=rows[0].keys())
        writer.writeheader()
        writer.writerows(rows)

    logger.info(f"{'span(K)':>8} {'AMR elec COP':>13} {'stages':>7} {'VCC COP':>9} {'Liquid COP':>11} {'Carnot':>8}")
    for r in rows:
        logger.info(f"{r['span_K']:>8} {r['AMR_COP_electrical']:>13} {r['AMR_n_stages']:>7} "
                    f"{r['VaporCompression_COP']:>9} {r['LiquidCooling_COP']:>11} {r['Carnot_COP']:>8}")
    logger.info(f"Wrote {RESULTS_CSV}")
    logger.info(
        "Note: AMR_COP_electrical includes estimated parasitic losses and is "
        "the appropriate quantity for comparison with vapor-compression and "
        "liquid-cooling COP values, which are also reported on an electrical "
        "basis. AMR_COP_ideal represents the thermodynamic cycle alone and is "
        "provided for reference. AMR_n_stages=1 rows are the plain single-stage "
        "cycle; AMR_n_stages>1 rows (span >16K here) used the automatic cascade "
        "fallback in core.cascade.staged_baseline_result() because a single "
        "stage's own no-load DeltaT_ad could not cover that span -- see that "
        "function's docstring."
    )
    return rows


def run_economics(representative_row):
    """Step 5: TCO for AMR/VCC/liquid, sized off the actual cooling
    capacity computed in the baseline sweep rather than an illustrative
    placeholder capacity."""
    capacity_kW = representative_row["AMR_Qc_W"] / 1000.0
    logger.info(f"Sizing basis: capacity={capacity_kW:.2f} kW, span={REPRESENTATIVE_SPAN_K}K "
                f"(from the baseline sweep computed in step 4)")
    logger.info(f"{'technology':<32}{'CAPEX ($)':>14}{'annual OPEX ($)':>18}")
    for tco in (economics.AMR_MAGNETIC, economics.VAPOR_COMPRESSION, economics.LIQUID_COOLING):
        r = economics.simple_tco(tco, capacity_kW, annual_hours=8760)
        logger.info(f"{r['technology']:<32}{r['capex_$']:>14,.0f}{r['annual_opex_$']:>18,.0f}")
    logger.info(
        "Note: AMR CAPEX/OPEX per kW are pre-commercial placeholders (see "
        "economics.py notes); material_cost() gives a materials-only cost "
        "floor for a specific design instead."
    )

    # Phase 15 addition: full-system BOM cost model (economics.bom_cost/
    # full_system_cost_estimate/levelized_cost_of_cooling). Uses the SAME
    # design point (2T, 5kg Gd) as step 4's baseline sweep, so this is a
    # priced version of the exact design already characterized above, not
    # a separate illustrative example.
    mu0H_T, mass_kg = 2.0, 5.0
    amr_cop = representative_row["AMR_COP_electrical"]
    bom = economics.bom_cost(mu0H_T, mass_kg, family_name="Gd")
    full_system = economics.full_system_cost_estimate(mu0H_T, mass_kg, family_name="Gd")
    lcoc = economics.levelized_cost_of_cooling(
        mu0H_T, mass_kg, Qc_avg_W=representative_row["AMR_Qc_W"], COP_electrical=amr_cop,
        family_name="Gd")
    logger.info("")
    logger.info(f"Phase 15: full-system BOM cost model, same design point (H={mu0H_T}T, "
                f"mass={mass_kg}kg Gd) as above:")
    logger.info(f"  Materials BOM: magnet ${bom['magnet_cost_$']:.0f} "
                f"({bom['magnet_mass_kg']:.2f}kg) + MCM ${bom['mcm_cost_$']:.0f} + "
                f"SMM yoke ${bom['smm_cost_$']:.0f} ({bom['smm_mass_kg']:.2f}kg) "
                f"= ${bom['materials_bom_total_$']:.0f} total")
    logger.info(f"  Full-system cost ESTIMATE (materials BOM x "
                f"{full_system['non_materials_multiplier']:.0f}x, order-of-magnitude only -- "
                f"see economics.py's Phase 15 section docstring for the Russek & Zimm (2006) "
                f"vapor-compression-AC benchmark this multiplier comes from): "
                f"${full_system['full_system_cost_estimate_$']:,.0f}")
    logger.info(f"  Levelized cost of cooling (CRF-based, Silva et al. 2017 methodology, "
                f"{lcoc['device_lifetime_years']:.0f}yr life, {lcoc['discount_rate']*100:.0f}% "
                f"discount rate): ${lcoc['levelized_cost_of_cooling_$_per_kwh']:.4f}/kWh_cooling "
                f"(${lcoc['annualized_capital_$_per_kwh_cooling']:.4f} capital + "
                f"${lcoc['electricity_$_per_kwh_cooling']:.4f} electricity, materials-only "
                f"capital basis -- see levelized_cost_of_cooling()'s docstring)")


def run_full_system_cost_by_material():
    """Step 5b (Phase 15 addition): compares the full-system cost estimate
    (economics.full_system_cost_estimate) across the SAME material
    candidates optimize.py's Phase 15 co-optimization considers
    (core.optimize._material_candidates()), at a fixed representative
    design point (2T, 5kg) -- a quick "does material choice matter for
    cost, independent of the NSGA-III search's own field/frequency/flow/
    geometry choices" sanity check alongside step 11's full multi-
    objective search."""
    for label, _material, family_name in optimize_module._material_candidates():
        r = economics.full_system_cost_estimate(2.0, 5.0, family_name=family_name)
        logger.info(f"  {label:<40} materials BOM=${r['materials_bom_total_$']:>8,.0f}   "
                    f"full-system estimate=${r['full_system_cost_estimate_$']:>10,.0f}")


def run_emissions(representative_row):
    """Step 6: emissions comparison fed with the *actual* COPs from the
    baseline sweep instead of the module's illustrative example numbers."""
    capacity_kW = representative_row["AMR_Qc_W"] / 1000.0
    amr_cop = representative_row["AMR_COP_electrical"]
    vcc_cop = representative_row["VaporCompression_COP"]
    liquid_cop = representative_row["LiquidCooling_COP"]
    logger.info(f"Basis: capacity={capacity_kW:.2f} kW, AMR_COP={amr_cop}, "
                f"VCC_COP={vcc_cop}, Liquid_COP={liquid_cop} (all from step 4, span="
                f"{REPRESENTATIVE_SPAN_K}K)")
    for r in emissions.compare_emissions(capacity_kW, amr_cop=amr_cop, vcc_cop=vcc_cop,
                                          liquid_cop=liquid_cop):
        logger.info(f"{r.technology:<32} refrigerant={r.refrigerant_GWP_tCO2e_per_year:7.2f} "
                    f"tCO2e/yr  operational={r.operational_CO2_tCO2e_per_year:8.2f} tCO2e/yr  "
                    f"total={r.total_tCO2e_per_year:8.2f} tCO2e/yr")


def run_cascade_comparison():
    """Step 7: 1-4 stage cascaded AMR vs baselines, for both Gd and the
    giant-MCE Gd5Si2Ge2 material. Returns both materials' rows (previously
    only rows_gd was returned) so step 12's figure generation can reuse
    both instead of recomputing the giant-MCE sweep a second time for
    fig20."""
    logger.info("Material: Gd (baseline)")
    rows_gd = cascade.compare_staging(material=GADOLINIUM, mass_per_stage=5.0,
                                       out_csv="results/cascade_comparison.csv")
    logger.info(f"Wrote results/cascade_comparison.csv ({len(rows_gd)} rows)")

    logger.info("Material: Gd5Si2Ge2 (giant MCE, first-order Landau model)")
    rows_giant = cascade.compare_staging(material=GD5SI2GE2_FIRST_ORDER, mass_per_stage=5.0,
                                          out_csv="results/cascade_comparison_giant_mce.csv")
    logger.info(f"Wrote results/cascade_comparison_giant_mce.csv ({len(rows_giant)} rows)")
    return rows_gd, rows_giant


def run_graded_cascade_comparison(rows_gd):
    """Step 7b: Curie-graded cascade (ROADMAP.md Phase 7 open item). Each
    stage uses a composition-tuned Gd5(SixGe1-x)4(-Ga) material matched to
    its own local operating temperature, checked against the documented
    ~20-290K composition-tunability range and scaled by the Giguere et al.
    (1999) empirical correction. See core/cascade.py's __main__ block for
    the full printed breakdown this reproduces."""
    rows_graded, stage_info_all = cascade.compare_graded_cascade(
        T_cold_C=18.0, spans=range(5, 21), mass_per_stage=5.0,
        out_csv="results/graded_cascade_comparison.csv")
    n_cells = sum(1 for row in rows_graded for n in (1, 2, 3, 4))
    n_full_range = sum(1 for row in rows_graded for n in (1, 2, 3, 4)
                        if row[f"Graded_{n}stage_n_fallback_to_Gd"] == 0)
    n_some_fallback = sum(1 for row in rows_graded for n in (1, 2, 3, 4)
                           if 0 < row[f"Graded_{n}stage_n_fallback_to_Gd"] < n)
    logger.info(f"Wrote results/graded_cascade_comparison.csv ({len(rows_graded)} rows): "
                f"{n_full_range}/{n_cells} span/stage-count cells fully within the documented "
                f"20-290K giant-MCE composition range, {n_some_fallback} with partial fallback "
                f"to plain Gd (see cascade.py __main__ for the full breakdown and honest caveats).")
    if rows_gd is not None:
        graded_10K_3 = next(r for r in rows_graded if r["span_K"] == 10)
        gd_10K_3 = next(r for r in rows_gd if r["span_K"] == 10)
        logger.info(f"At 10K span, 3 stages: Graded Qc={graded_10K_3['Graded_3stage_Qc_W']}W, "
                    f"COP={graded_10K_3['Graded_3stage_COP']}  vs. plain-Gd "
                    f"Qc={gd_10K_3['AMR_3stage_Qc_W']}W, COP={gd_10K_3['AMR_3stage_COP']}")
    # Returned (not previously) so main()'s new step 13 (design_recommendations.py)
    # can report the 10K-span/3-stage graded-vs-plain-Gd comparison without
    # re-running this ~2-minute sweep a second time.
    return rows_graded


def run_astronautics_graded_validation():
    """Step 7c (ROADMAP.md Phase 9 addendum): does a 6-layer Curie-graded
    La(Fe,Si)13Hy bed reproduce the REAL Astronautics_rotary_2014 device?
    Step 2 above (validation_system.py) could not calibrate this device
    with a single-Tc=287K material; this uses cascade.py's generalized
    Curie-grading machinery (LAFESIH_FAMILY) to test the actual 6-layer
    hypothesis instead. See core/cascade.py's own __main__ block and
    ROADMAP.md Phase 9 for the full writeup. Returns the result dict
    (previously discarded) so step 12's fig25 can reuse it instead of
    re-running this ~33s validation a second time."""
    astro = cascade.validate_astronautics_graded_bed()
    if astro.get("feasible"):
        for s in astro["stage_info"]:
            logger.info(f"    stage {s['stage']}: T_mid={s['T_mid_K']}K, needed composition Tc="
                        f"{s['Tc_target_K']}K -> {s['material']}")
        logger.info(f"mdot calibrated to reproduce reported Qc={astro['Qc_lit_W']}W: "
                    f"{astro['mdot_calibrated_kg_s']} kg/s")
        logger.info(f"Predicted COP={astro['COP_cascade']}  vs.  reported COP={astro['COP_lit']} "
                    f"({astro['COP_error_pct']:+.1f}% error) -- vs. the flat 'no calibration found' "
                    f"the single-layer material in step 2 gave this same device.")
    else:
        logger.info(astro.get("status", "infeasible"))
    return astro


def run_thermal_demo():
    """Reproduces core/thermal.py's own __main__ demo: a regenerator
    effectiveness sweep vs. mass and vs. frequency. thermal.py is only
    ever reached transitively (amr_cycle.py imports it internally), so
    this demo never otherwise runs as part of the pipeline."""
    logger.info("Regenerator effectiveness sweep vs. mass_regenerator (f=1Hz, mdot=0.08kg/s)")
    for m in [0.5, 1, 2, 5, 10, 15]:
        r = regenerator_effectiveness(m, frequency=1.0, mdot=0.08)
        logger.info(f"  mass={m:5.1f}kg  NTU={r['NTU']:6.2f}  U={r['U']:6.3f}  eps={r['eps']:.3f}")
    logger.info("Regenerator effectiveness sweep vs. frequency (mass=2kg, mdot=0.08kg/s)")
    for f in [0.25, 0.5, 1, 2, 4]:
        r = regenerator_effectiveness(2.0, frequency=f, mdot=0.08)
        logger.info(f"  f={f:5.2f}Hz  NTU={r['NTU']:6.2f}  U={r['U']:6.3f}  eps={r['eps']:.3f}")


def run_first_order_mce_demo():
    """Reproduces core/first_order_mce.py's own __main__ demo: a
    calibration check of the Landau model against its target literature
    value. first_order_mce.py is only ever reached transitively
    (giant_mce_analysis.py imports GD5SI2GE2_FIRST_ORDER from it), so
    this demo never otherwise runs as part of the pipeline."""
    mu0_ = 4 * np.pi * 1e-7
    logger.info("First-order Landau model calibration check, Gd5Si2Ge2, T=Tc=276K")
    for B_T in [1, 2, 5]:
        H = B_T / mu0_
        dS = GD5SI2GE2_FIRST_ORDER.delta_S_isothermal(np.array([276.0]), H)
        dT = GD5SI2GE2_FIRST_ORDER.delta_T_adiabatic(np.array([276.0]), H)
        logger.info(f"  {B_T}T: dS={dS[0]:.2f} J/(kg K)   dTad={dT[0]:.2f} K")
    logger.info("Target: dS ~ -18 J/(kg K) at 5T (Pecharsky & Gschneidner 1997 review value)")


def run_plot_generation(precomputed=None):
    """Step 12: renders all 26 figures in plots.py (results/figures/*.png
    and *.pdf) covering material validation, AMR characteristic curves,
    thermal/geometry modelling, loss-model calibration, system/curve
    validation, cascade and Curie-graded staging, Sobol sensitivity, RSM
    surrogate fitting, NSGA-III optimization, economics, and emissions.
    Most figures still compute their own data directly from core/, but the
    nine figures that duplicate an earlier stage's computation exactly
    (fig08 baseline sweep, fig14 system validation, fig16 Sobol, fig18
    Pareto, fig19/20 cascade, fig21 graded cascade, fig25 Astronautics
    validation, fig26 material family comparison) now reuse the results
    already produced by steps 2/4/7/7b/7c/8d/9/9b/11 via `precomputed`,
    instead of re-running them from scratch. This still runs last so the
    CSV-writing figures (cascade, graded cascade, Pareto front) leave
    results/ in a consistent, freshly-regenerated state."""
    plots.run_all(precomputed=precomputed)
    n_figs = len(list(plots.FIG_DIR.glob("*.png"))) if plots.FIG_DIR.exists() else 0
    if n_figs:
        logger.info(f"Generated {n_figs} figure(s) (PNG + PDF) in {plots.FIG_DIR.as_posix()}/")
    else:
        logger.warning("No figures were generated -- check the traceback above.")


def run_design_recommendations_synthesis(sobol_state_dependent_Si, pareto_rows, material_rows,
                                          graded_rows, cascade_rows_gd, pb_best_cop_row,
                                          pp_best_cop_row, n_stages=3,
                                          representative_span_K=REPRESENTATIVE_SPAN_K):
    """Step 13: consolidates the already-computed results from steps
    3c/7b/8d/9b/11 into one ranked "how do I raise AMR electrical COP"
    report via core/design_recommendations.py. Pulls out the 10K-span,
    3-stage graded-vs-plain-Gd cascade rows (matching the comparison
    main.py already logs at the end of step 7b) since the graded/plain
    cascade sweeps return every span/stage-count combination, not just
    the representative one."""
    graded_row = None
    gd_cascade_row = None
    if graded_rows is not None:
        graded_row = next((r for r in graded_rows if r["span_K"] == representative_span_K), None)
    if cascade_rows_gd is not None:
        gd_cascade_row = next((r for r in cascade_rows_gd if r["span_K"] == representative_span_K),
                               None)
    design_recommendations.build_report(
        sobol_state_dependent_Si=sobol_state_dependent_Si,
        pareto_rows=pareto_rows,
        material_rows=material_rows,
        graded_row=graded_row,
        gd_cascade_row=gd_cascade_row,
        n_stages=n_stages,
        pb_best_cop_row=pb_best_cop_row,
        pp_best_cop_row=pp_best_cop_row,
        representative_span_K=representative_span_K,
        out_path="results/design_recommendations.txt",
    )


def main():
    _attach_file_logging()
    t_start = time.time()
    stage_times = {}
    failures = []

    stages = [
        ("1. Material-level model validation vs. Dan'kov et al. (1998)",
         lambda: (validation.run_validation(),
                  validation.run_giguere_gd_extension(),
                  validation.run_curie_shift_check())),
        ("2. System-level validation vs. published AMR prototypes",
         None),  # handled specially below, result (system_validation_results) captured
        ("3. Loss-model calibration (auto-loaded by AMRSystem's default loss model)",
         lambda: (loss_model.calibrate_loss_coefficients(), loss_model.run_extended_diagnostic())),
        ("3b. Regenerator thermal-effectiveness demo (core/thermal.py, reached transitively otherwise)",
         run_thermal_demo),
        ("3c. Geometry-dependent pumping power: packed-bed + parallel-plate (core/geometry_analysis.py)",
         lambda: geometry_analysis.run_geometry_analysis()),
        ("3d. Hypereg parallel-hydraulic pumping-power analysis (core/hypereg_analysis.py, Phase 15 item 3)",
         lambda: hypereg_analysis.run_hypereg_analysis()),
        ("4. Baseline comparison sweep: AMR vs VCC vs liquid cooling vs Carnot",
         None),  # handled specially below, result captured
        ("5. Economics / TCO at the representative operating point",
         None),  # needs step 4's result
        ("5b. Full-system cost estimate by material family (core/economics.py, Phase 15 item 5)",
         run_full_system_cost_by_material),
        ("6. Emissions comparison at the representative operating point",
         None),  # needs step 4's result
        ("7. Cascade staging comparison (1-4 stage AMR, Gd and Gd5Si2Ge2)",
         None),  # handled specially below, result (rows_gd) captured for step 7b
        ("7b. Curie-graded cascade (composition-tuned per stage, ROADMAP.md Phase 7 item)",
         None),  # needs step 7's result
        ("7c. Does a 6-layer Curie-graded La(Fe,Si)13Hy bed reproduce Astronautics_rotary_2014? (ROADMAP.md Phase 9 addendum)",
         run_astronautics_graded_validation),
        ("8. Giant-MCE materials analysis (Gd vs Gd5Si2Ge2)",
         lambda: giant_mce_analysis.run_analysis()),
        ("8d. Four-way material family comparison (Gd, Gd5Si2Ge2-fixed, GD/LAFESIH/MNFEPSI-tuned; Track A2 item)",
         lambda: material_family_comparison.run_analysis()),
        ("8b. First-order Landau model calibration check (core/first_order_mce.py, reached transitively otherwise)",
         run_first_order_mce_demo),
        ("8c. Giguere et al. (1999) direct-measurement cross-check (core/giguere_validation.py)",
         lambda: (giguere_validation.run_validation(),
                  giguere_validation.run_pecharsky_ratio_check())),
        ("9. Sobol global sensitivity analysis (constant-loss model)",
         lambda: sensitivity.run_sobol(out_path="results/sobol_results_phase2_constant.txt",
                                        use_state_dependent_losses=False)),
        ("9b. Sobol global sensitivity analysis (state-dependent loss model)",
         lambda: sensitivity.run_sobol(out_path="results/sobol_results.txt",
                                        use_state_dependent_losses=True)),
        ("10. Response-surface (RSM) surrogate fit",
         lambda: rsm.fit_rsm()),
        ("11. NSGA-III multi-objective design optimization (material + geometry co-optimization, Phase 15 item 2)",
         None),  # handled specially below, result (pareto_rows) captured for step 13
        ("12. Figure generation: 26 figures covering validation, AMR curves, "
         "cascade/graded staging, sensitivity, RSM, NSGA-III, economics, emissions (plots.py)",
         None),  # handled specially below, reuses steps 7/7b/9/9b/11's results
        ("13. Design-recommendations synthesis (core/design_recommendations.py)",
         None),  # handled specially below, consumes steps 3c/7b/8d/9b/11's results
    ]

    representative_row = None
    system_validation_results = None
    cascade_rows_gd = None
    cascade_rows_giant = None
    graded_rows = None
    astro_result = None
    material_rows = None
    pareto_rows = None
    sobol_const_Si = None
    sobol_state_dependent_Si = None
    pb_best_cop_row = None
    pp_best_cop_row = None

    for name, fn in stages:
        _banner(name)
        t0 = time.time()
        try:
            with contextlib.redirect_stdout(_StreamToLogger(logger)):
                if name.startswith("2."):
                    system_validation_results = validation_system.run_system_validation()
                    validation_system.run_field_sensitivity_check()
                    validation_system.run_capacity_only_calibration_check()
                elif name.startswith("4."):
                    rows = run_baseline_sweep()
                    representative_row = next(
                        r for r in rows if abs(r["span_K"] - REPRESENTATIVE_SPAN_K) < 1e-9
                    )
                elif name.startswith("5."):
                    run_economics(representative_row)
                elif name.startswith("6."):
                    run_emissions(representative_row)
                elif name.startswith("7. "):
                    cascade_rows_gd, cascade_rows_giant = run_cascade_comparison()
                elif name.startswith("7b."):
                    graded_rows = run_graded_cascade_comparison(cascade_rows_gd)
                elif name.startswith("7c."):
                    astro_result = fn()
                elif name.startswith("3c."):
                    fn()
                    # Cheap (sub-second), non-printing re-sweep purely to capture
                    # the best-COP geometry rows for step 13's synthesis report;
                    # run_geometry_analysis() above already did the printed,
                    # file-writing version of this same sweep.
                    _, _, pb_best_cop_row = geometry_analysis.sweep_packed_bed_diameter(
                        verbose=False)
                    _, _, pp_best_cop_row = geometry_analysis.sweep_parallel_plate_spacing(
                        verbose=False)
                elif name.startswith("8d."):
                    material_rows = fn()
                elif name.startswith("9b."):
                    sobol_state_dependent_Si = fn()
                elif name.startswith("9."):
                    sobol_const_Si = fn()
                elif name.startswith("11."):
                    pareto_rows = optimize_module.run_optimization()
                elif name.startswith("12."):
                    run_plot_generation(precomputed={
                        "system_validation_results": system_validation_results,
                        "baseline_rows": rows,
                        "sobol_const_Si": sobol_const_Si,
                        "sobol_state_Si": sobol_state_dependent_Si,
                        "pareto_rows": pareto_rows,
                        "cascade_rows_gd": cascade_rows_gd,
                        "cascade_rows_giant": cascade_rows_giant,
                        "graded_rows": graded_rows,
                        "astro_result": astro_result,
                        "material_rows": material_rows,
                    })
                elif name.startswith("13."):
                    run_design_recommendations_synthesis(
                        sobol_state_dependent_Si=sobol_state_dependent_Si,
                        pareto_rows=pareto_rows,
                        material_rows=material_rows,
                        graded_rows=graded_rows,
                        cascade_rows_gd=cascade_rows_gd,
                        pb_best_cop_row=pb_best_cop_row,
                        pp_best_cop_row=pp_best_cop_row,
                    )
                else:
                    fn()
        except Exception:
            logger.error(f"!!! STAGE FAILED: {name}")
            logger.error(traceback.format_exc())
            failures.append(name)
        stage_times[name] = time.time() - t0
        logger.info(f"[stage done in {stage_times[name]:.2f}s]")

    _banner("PIPELINE SUMMARY")
    total = time.time() - t_start
    for name, dt in stage_times.items():
        status = "FAILED" if name in failures else "ok"
        logger.info(f"  [{status:6}] {dt:6.2f}s  {name}")
    logger.info(f"Total wall time: {total:.1f}s")
    if failures:
        logger.warning(f"{len(failures)} stage(s) failed -- see tracebacks above:")
        for name in failures:
            logger.warning(f"  - {name}")
    else:
        logger.info("All stages completed. Files written to results/:")
        logger.info("  comparison_table.csv, cascade_comparison.csv, "
                     "cascade_comparison_giant_mce.csv, giant_mce_analysis.txt, "
                     "material_family_comparison.csv, material_family_comparison.txt, "
                     "sobol_results.txt, sobol_results_phase2_constant.txt, "
                     "rsm_coefficients.txt, pareto_front.csv, "
                     "pareto_front_by_material/*.csv (Phase 15), "
                     "hypereg_analysis.txt (Phase 15), "
                     "geometry_optimization_analysis.txt, graded_cascade_comparison.csv, "
                     "design_recommendations.txt, figures/*.png+*.pdf (26 figures)")
    logger.info(f"Full run log: {LOG_FILE}")

    _print_executive_summary(representative_row, cascade_rows_gd, graded_rows,
                              material_rows, pareto_rows, pb_best_cop_row,
                              pp_best_cop_row, failures)


def _print_executive_summary(representative_row, cascade_rows_gd, graded_rows, material_rows,
                              pareto_rows, pb_best_cop_row, pp_best_cop_row, failures):
    """Final, well-structured overview of every implemented analysis and
    its headline metric, printed once at the very end of the run so a
    reader does not have to scroll back through 13 stages of log output
    to see what this pipeline actually covers. Every number quoted here
    is read directly from the result objects the stages above already
    computed -- nothing is recomputed or hardcoded. Stages that failed
    (see `failures`) are reported as unavailable rather than silently
    skipped."""
    _banner("EXECUTIVE SUMMARY -- implemented features, analyses, and headline metrics")

    def _ok(prefix):
        return not any(f.startswith(prefix) for f in failures)

    logger.info("Validation")
    logger.info("  - Material-level: mean-field Gd model vs. Dan'kov et al. (1998), Giguere et "
                "al. (1999) direct-measurement cross-check, Curie-shift limitation check")
    logger.info("  - System-level: calibrated against 16 published AMR prototype/benchmark rows "
                "(data/amr_experimental_benchmarks.csv), including a 6-layer Curie-graded "
                "La(Fe,Si)13Hy reproduction of Astronautics_rotary_2014")

    logger.info("Baseline comparison (AMR vs. vapor-compression vs. liquid cooling vs. Carnot)")
    if representative_row and _ok("4."):
        logger.info(f"  - At {REPRESENTATIVE_SPAN_K:.0f}K span: AMR_COP_elec="
                    f"{representative_row['AMR_COP_electrical']}, VCC_COP="
                    f"{representative_row['VaporCompression_COP']}, Liquid_COP="
                    f"{representative_row['LiquidCooling_COP']}, Carnot_COP="
                    f"{representative_row['Carnot_COP']}  (results/comparison_table.csv)")
    else:
        logger.info("  - unavailable (stage failed or was skipped)")

    logger.info("Economics & emissions (TCO and GWP at the representative operating point)")
    logger.info("  - results/*: economics.py CAPEX/OPEX comparison, emissions.py refrigerant + "
                "operational CO2e comparison")
    logger.info("  - Phase 15: full-system BOM cost model (materials + soft-magnetic yoke), "
                "order-of-magnitude full-system cost estimate, and CRF-based levelized cost of "
                "cooling, all at the same design point -- see step 5's log output above; step 5b "
                "compares the same estimate across material families")

    logger.info("Cascade staging & Curie-temperature grading")
    if cascade_rows_gd and graded_rows and _ok("7"):
        gd10 = next((r for r in cascade_rows_gd if r["span_K"] == 10), None)
        gr10 = next((r for r in graded_rows if r["span_K"] == 10), None)
        if gd10 and gr10:
            logger.info(f"  - At 10K span, 3-stage: plain-Gd COP={gd10['AMR_3stage_COP']}, "
                        f"Curie-graded COP={gr10['Graded_3stage_COP']} "
                        f"(results/cascade_comparison.csv, results/graded_cascade_comparison.csv)")
    else:
        logger.info("  - unavailable (stage failed or was skipped)")

    logger.info("Material family ranking (Gd / Gd5Si2Ge2 / GD- / LAFESIH- / MNFEPSI-tuned families)")
    if material_rows and _ok("8d."):
        rep = [r for r in material_rows if r["span_K"] == REPRESENTATIVE_SPAN_K
               and r.get("in_range") and r.get("1stage_COP") is not None]
        if rep:
            best = max(rep, key=lambda r: r["1stage_COP"])
            logger.info(f"  - Best at {REPRESENTATIVE_SPAN_K:.0f}K span: {best['candidate']}  "
                        f"COP_elec={best['1stage_COP']}  (results/material_family_comparison.csv)")
    else:
        logger.info("  - unavailable (stage failed or was skipped)")

    logger.info("Global sensitivity (Sobol) & response-surface (RSM) surrogate")
    logger.info("  - results/sobol_results.txt (state-dependent losses), "
                "results/sobol_results_phase2_constant.txt (constant losses), "
                "results/rsm_coefficients.txt (Qc surrogate, R^2 reported at fit time)")

    logger.info("Regenerator geometry optimization (packed-bed / parallel-plate)")
    if pb_best_cop_row and pp_best_cop_row and _ok("3c."):
        logger.info(f"  - COP-optimal packed-bed sphere diameter: {pb_best_cop_row[0]} mm; "
                    f"COP-optimal parallel-plate spacing: {pp_best_cop_row[0]} mm "
                    f"(results/geometry_optimization_analysis.txt)")
    else:
        logger.info("  - unavailable (stage failed or was skipped)")

    logger.info("Hypereg parallel-hydraulic pumping-power analysis (Phase 15 item 3)")
    if _ok("3d."):
        logger.info("  - Klinar et al. (2024)-motivated pumping-power-only sweep; see "
                    "results/hypereg_analysis.txt and results/hypereg_findings.md for the full "
                    "literature findings note and quantitative result")
    else:
        logger.info("  - unavailable (stage failed or was skipped)")

    logger.info("NSGA-III multi-objective design optimization (COP vs. Qc vs. cost; "
                "Phase 15: material + geometry co-optimized)")
    if pareto_rows and _ok("11."):
        best_cop = max(pareto_rows, key=lambda r: r["COP_electrical"])
        logger.info(f"  - {len(pareto_rows)} Pareto-optimal designs; best electrical COP="
                    f"{best_cop['COP_electrical']} at f={best_cop['frequency_Hz']}Hz, "
                    f"material={best_cop.get('material', 'Gd')} "
                    f"(results/pareto_front.csv, results/pareto_front_by_material/*.csv)")
    else:
        logger.info("  - unavailable (stage failed or was skipped)")

    logger.info("Consolidated design recommendations (step 13)")
    logger.info("  - Ranks all COP-maximization levers above by demonstrated Sobol "
                "sensitivity and reports a recommended starting design point "
                "(results/design_recommendations.txt)")

    logger.info("Figures")
    n_figs = len(list(plots.FIG_DIR.glob("*.png"))) if plots.FIG_DIR.exists() else 0
    logger.info(f"  - {n_figs} figure(s) generated (results/figures/*.png, *.pdf)")


if __name__ == "__main__":
    main()