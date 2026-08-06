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
    4.  Baseline comparison sweep      (this file, was the old main.py)
    5.  Economics / TCO                (core/economics.py)
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
    11. NSGA-III design optimization   (core/optimize.py)
    12. Figure generation (26 figures) (plots.py -> results/figures/*.png, *.pdf)

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

Step 12 runs last. plots.py is largely self-contained -- most of its 25
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
from core import plots

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
    giant-MCE Gd5Si2Ge2 material."""
    logger.info("Material: Gd (baseline)")
    rows_gd = cascade.compare_staging(material=GADOLINIUM, mass_per_stage=5.0,
                                       out_csv="results/cascade_comparison.csv")
    logger.info(f"Wrote results/cascade_comparison.csv ({len(rows_gd)} rows)")

    logger.info("Material: Gd5Si2Ge2 (giant MCE, first-order Landau model)")
    rows_giant = cascade.compare_staging(material=GD5SI2GE2_FIRST_ORDER, mass_per_stage=5.0,
                                          out_csv="results/cascade_comparison_giant_mce.csv")
    logger.info(f"Wrote results/cascade_comparison_giant_mce.csv ({len(rows_giant)} rows)")
    return rows_gd


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


def run_astronautics_graded_validation():
    """Step 7c (ROADMAP.md Phase 9 addendum): does a 6-layer Curie-graded
    La(Fe,Si)13Hy bed reproduce the REAL Astronautics_rotary_2014 device?
    Step 2 above (validation_system.py) could not calibrate this device
    with a single-Tc=287K material; this uses cascade.py's generalized
    Curie-grading machinery (LAFESIH_FAMILY) to test the actual 6-layer
    hypothesis instead. See core/cascade.py's own __main__ block and
    ROADMAP.md Phase 9 for the full writeup."""
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


def run_plot_generation():
    """Step 12: renders all 26 figures in plots.py (results/figures/*.png
    and *.pdf) covering material validation, AMR characteristic curves,
    thermal/geometry modelling, loss-model calibration, system/curve
    validation, cascade and Curie-graded staging, Sobol sensitivity, RSM
    surrogate fitting, NSGA-III optimization, economics, and emissions.
    Most figures recompute their own data directly from core/ rather than
    re-reading the CSVs written above, but this still runs last so the
    CSV-writing figures (cascade, graded cascade, Pareto front) leave
    results/ in a consistent, freshly-regenerated state."""
    plots.run_all()
    n_figs = len(list(plots.FIG_DIR.glob("*.png"))) if plots.FIG_DIR.exists() else 0
    if n_figs:
        logger.info(f"Generated {n_figs} figure(s) (PNG + PDF) in {plots.FIG_DIR.as_posix()}/")
    else:
        logger.warning("No figures were generated -- check the traceback above.")


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
         lambda: (validation_system.run_system_validation(),
                  validation_system.run_field_sensitivity_check(),
                  validation_system.run_capacity_only_calibration_check())),
        ("3. Loss-model calibration (auto-loaded by AMRSystem's default loss model)",
         lambda: (loss_model.calibrate_loss_coefficients(), loss_model.run_extended_diagnostic())),
        ("3b. Regenerator thermal-effectiveness demo (core/thermal.py, reached transitively otherwise)",
         run_thermal_demo),
        ("3c. Geometry-dependent pumping power: packed-bed + parallel-plate (core/geometry_analysis.py)",
         lambda: geometry_analysis.run_geometry_analysis()),
        ("4. Baseline comparison sweep: AMR vs VCC vs liquid cooling vs Carnot",
         None),  # handled specially below, result captured
        ("5. Economics / TCO at the representative operating point",
         None),  # needs step 4's result
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
        ("11. NSGA-III multi-objective design optimization",
         lambda: optimize_module.run_optimization()),
        ("12. Figure generation: 26 figures covering validation, AMR curves, "
         "cascade/graded staging, sensitivity, RSM, NSGA-III, economics, emissions (plots.py)",
         run_plot_generation),
    ]

    representative_row = None
    cascade_rows_gd = None

    for name, fn in stages:
        _banner(name)
        t0 = time.time()
        try:
            with contextlib.redirect_stdout(_StreamToLogger(logger)):
                if name.startswith("4."):
                    rows = run_baseline_sweep()
                    representative_row = next(
                        r for r in rows if abs(r["span_K"] - REPRESENTATIVE_SPAN_K) < 1e-9
                    )
                elif name.startswith("5."):
                    run_economics(representative_row)
                elif name.startswith("6."):
                    run_emissions(representative_row)
                elif name.startswith("7. "):
                    cascade_rows_gd = run_cascade_comparison()
                elif name.startswith("7b."):
                    run_graded_cascade_comparison(cascade_rows_gd)
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
                     "geometry_optimization_analysis.txt, figures/*.png+*.pdf (26 figures)")
    logger.info(f"Full run log: {LOG_FILE}")


if __name__ == "__main__":
    main()