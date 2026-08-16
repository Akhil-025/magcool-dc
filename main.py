"""
main.py
=======
End-to-end pipeline for the magcool-dc project. Runs every analysis module
in the repository in one pass, in dependency order, so a single
``python main.py`` reproduces every file under results/ from scratch:

    1.  Material-level validation      (core/validation.py)
    1b. Inhomogeneous/polycrystalline Tc-broadening sensitivity (core/
        inhomogeneous_broadening.py, Phase 22 item 1) -- runs right after
        step 1 since it directly extends that same Dan'kov et al. (1998)
        comparison
    2.  System-level validation        (core/validation_system.py)
    2b. Cycle-type validation sensitivity (core/validation_system.py, Phase 17)
    3.  Loss-model calibration report  (core/loss_model.py)
    3b. Regenerator thermal demo       (core/thermal.py)
    3c. Geometry-dependent pumping power (core/geometry_analysis.py, ROADMAP.md Phase 7 item)
    3d. Hypereg parallel-hydraulic pumping-power analysis (core/hypereg_analysis.py, Phase 15 item 3)
    3e. Geometry-explicit intragranular eddy-current loss + pump/motor
        efficiency demo (core/thermal.py, core/amr_cycle.py, Phase 27/28)
    4.  Baseline comparison sweep      (this file, was the old main.py)
    5.  Economics / TCO                (core/economics.py)
    5b. Full-system cost estimate by material family (core/economics.py, Phase 15 item 5)
    6.  Emissions comparison           (core/emissions.py)
    7.  Cascade staging comparison     (core/cascade.py)
    7b. Curie-graded cascade           (core/cascade.py, ROADMAP.md Phase 7 item)
    7c. Astronautics graded-bed check  (core/cascade.py, ROADMAP.md Phase 9 addendum)
    7d. Extending the graded-bed structural fix to the remaining STRUCTURAL
        devices (core/cascade.py, Paper-Mining Pass review item 1 follow-up)
    7e. Does the Giguere DeltaT_ad correction narrow the Astronautics graded-bed
        COP error? (core/cascade.py, Paper-Mining Pass review item 2)
    8.  Giant-MCE materials analysis   (core/giant_mce_analysis.py)
    8d. Material family comparison     (core/material_family_comparison.py, Track A2 item)
    8b. First-order Landau model demo  (core/first_order_mce.py)
    8c. Giguere et al. (1999) direct-measurement cross-check + Pecharsky &
        Gschneidner (1997) peak-ratio check + latent-heat Cp spike
        (core/giguere_validation.py, Phase 26)
    9.  Sobol sensitivity analysis     (core/sensitivity.py)
    10. RSM surrogate fit              (core/rsm.py)
    11. NSGA-III design optimization   (core/optimize.py, Phase 15: material +
        geometry co-optimization, per-material-family fronts merged post-hoc)
    11b. Hysteresis sensitivity (ON/OFF Pareto-front A/B check, Phase 16)
        (core/hysteresis_sensitivity.py)
    11c. Thermal-diode sensitivity study, mechanical-contact active thermal
        diode (core/thermal_diode.py, core/thermal_diode_analysis.py, Phase 18)
    11d. Magnet-geometry (Halbach-cylinder) field-vs-mass cost model
        (core/magnet_geometry.py, Phase 19)
    11e. Magnet-geometry Pareto sensitivity, production-settings multi-seed
        stability check (core/magnet_geometry.py, Paper-Mining Pass review item 4)
    11f. Layered/graded-bed NSGA-III co-optimization (core/optimize.py,
        core/cascade.py, Phase 29) -- reduced-effort pipeline pass
        (n_layers=1-3, pop_size=20, n_gen=10); the full 1-6 layer,
        production-settings version is directly callable but not run
        here, purely for pipeline-runtime reasons
    12. Figure generation (34 figures) (plots.py -> results/figures/*.png, *.pdf)
    13. Design-recommendations synthesis (core/design_recommendations.py) --
        consolidates steps 3c/7b/8d/9b/11's already-computed results into
        one ranked, actionable "how do I raise AMR electrical COP" report
        (results/design_recommendations.txt)
    14. Magnetocaloric fluids (ferrofluid/MR suspension) as an alternative
        working-body class (core/fluid_mce_cycle.py, core/fluid_mce_analysis.py,
        Phase 20) -- design-exploration/comparison tool, not a validated
        feature (see that module's own honesty flags)
    15. Passive/hybrid magnetic regenerator: does loading a conventional
        (vapor-compression) gas cycle's internal regenerator with a
        magnetocaloric material's own Curie-point heat-capacity anomaly
        raise its COP? (core/baseline_cooling.py's augmented_regenerator_cop()
        + passive_regenerator_augmentation(), core/passive_regenerator_analysis.py,
        Phase 21) -- design-exploration/comparison tool with an illustrative,
        literature-range-anchored effectiveness-to-COP ceiling, not a
        validated feature (see that module's own honesty flags)

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
the figure count. plots.py is largely self-contained -- most of its 34
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

Phase 16 (see ROADMAP.md) added a quantified thermal-hysteresis loss term
(core/first_order_mce.py's new hysteresis_loss_J_per_kg field,
core/amr_cycle.py's AMRSystem._hysteresis_power_W()) for first-order MCE
materials, which Phase 15's NSGA-III material selection had no visibility
into at all. New step 11b (core/hysteresis_sensitivity.py) reruns
optimize.run_optimization() twice at identical pop_size/n_gen/seed --
once with each first-order candidate's hysteresis loss at its Phase-16
literature-placeholder value, once forced to 0.0 (exactly reproducing
pre-Phase-16 behavior) -- and reports how the merged Pareto front's
material composition shifts as a result (results/hysteresis_sensitivity.txt).
This is deliberately run at a smaller pop_size/n_gen than step 11's own
production settings to keep pipeline runtime reasonable; see that
module's docstring honesty flag #1 before treating its output as a
settled, publication-quality answer rather than a directional check.

Phase 17 (see ROADMAP.md) added an AMR cycle-topology switch to
core/amr_cycle.py's AMRSystem (`cycle_type`: "brayton" [default,
pre-Phase-17 behavior], "ericsson", "carnot" -- see that module's
CYCLE_TYPE_FACTORS for the honesty flag on these being illustrative,
qualitatively-ordered multipliers, not a digitization of Kitanovski et
al.'s own closed-form Sect. 4.1.1-4.1.4 relations, which this project's
copy of that book does not include). New step 2b
(validation_system.run_cycle_type_validation()) reruns the existing
system-level COP validation with each rotary-drive benchmark device's
cycle_type inferred as "ericsson" instead of the flat "brayton" default,
and reports whether the per-device COP prediction error shrinks. A full
NSGA-III categorical cycle_type search (mirroring how Phase 15 handled
material family) was deliberately NOT added in this pass -- see
ROADMAP.md's Phase 17 entry for why.

A follow-up pass (see ROADMAP.md's Phase 17 entry, "closed after the
fact") threaded cycle_type through core/cascade.py's `run_graded_cascade()`/
`validate_astronautics_graded_bed()` (Phase 17's own "did NOT do" item),
and step 7c now also reports the same brayton-vs-ericsson comparison for
the Astronautics_rotary_2014 graded-bed reproduction that step 2b already
gave DTU_Eriksen_rotary_Gd_2015 -- a genuine, single-device null result
(ericsson does not narrow this device's much larger -81.1% error), stated
plainly rather than omitted because it disagrees with the other rotary
device's result.

Phase 18 (see ROADMAP.md) added a narrowly-scoped mechanical-contact
active thermal diode model (core/thermal_diode.py's MechanicalContactDiode,
core/amr_cycle.py's AMRSystem `thermal_diode` parameter, default None =
pre-Phase-18 behavior unchanged). New step 11c
(core/thermal_diode_analysis.py) is a cost-only sensitivity study, NOT a
validated feature: it (1) directly checks and reports that this repo's
model has no internal mechanical-switching frequency ceiling for a
diode-assisted design to relax (the Phase 18 plan's own premise for that
part of the item), and (2) sweeps the (illustrative, unbenchmarked --
see that module's honesty flag) diode actuation switching-power cost
against frequency at the representative operating point. No offsetting
heat-transfer benefit from the diode's rectification_ratio is modeled
(no closed-form relation for Sect. 6.2.4 was available to digitize --
this project's copy of Kitanovski et al. does not include Ch. 6), and no
AMR benchmark device in this repo's corpus uses thermal diodes, so this
step is explicitly a design-exploration tool rather than a validated
result -- see ROADMAP.md's Phase 18 entry for the full scoping
discussion and what was deliberately not built.

Phase 19 (see ROADMAP.md) added core/magnet_geometry.py: a standard,
closed-form idealized-Halbach-cylinder relation for magnet mass vs.
field, usable (as an opt-in, not a default -- see that module's own
docstring for why) as a replacement for economics.py's pre-Phase-19 flat
per-Tesla magnet-mass ratio with a genuinely nonlinear (super-linear-in-
field) one. New step 11d runs `run_magnet_geometry_analysis()` (a cheap,
deterministic cost-per-Kelvin sweep, no NSGA-III) and a reduced-
resolution `run_geometric_cost_pareto_sensitivity()` A/B Pareto-front
comparison (flat vs. geometric magnet-mass cost term), the same
controlled-A/B pattern step 11b already established for Phase 16. This
pass also found and flagged (see magnet_geometry.py's HONESTY FLAG #2)
that the Phase 19 plan's own citation for the qualitative "~2 T is a
cost/performance sweet spot" claim pointed at the wrong Bjørk et al.
paper; the simple proxy this step checks that claim against does NOT
independently reproduce it -- reported honestly rather than massaged to
agree.

Phase 20 (see ROADMAP.md) added core/fluid_mce_cycle.py: FerrofluidMCESystem,
a NEW SIBLING to AMRSystem (not a parameter on it, per the plan's own
scoping) modeling magnetocaloric fluids (ferrofluid/magnetorheological
suspension) as an alternative working-body class -- a continuous
flow-through loop with no packed regenerator bed, using standard
Krieger-Dougherty suspension rheology and Darcy-Weisbach pipe pumping
power instead of core/thermal.py's packed-bed correlations, and a
mixture-heat-capacity dilution model for the suspension's own effective
adiabatic temperature change (see that module's own docstring for why
neither project book was needed, or available, for any of this). New
step 14 (core/fluid_mce_analysis.py) runs a particle-volume-fraction
sweep (the viscosity-vs-phi/MCE-intensity-vs-phi tradeoff the plan named
directly) and a comparison, at this architecture's own favorable (small)
span, against a representative solid AMRSystem and against
core/baseline_cooling.py's liquid-cooling/vapor-compression references.
No usable magnetocaloric-fluid-as-working-body benchmark was found in
this pass's literature search (see that module's HONESTY FLAG #2, which
also documents a distinct, adjacent ferrofluid-AS-THERMAL-SWITCH
technology this module deliberately does not conflate with its own
topic) -- this step is a design-exploration/comparison tool, not a
validated result, the same disposition Phase 18 gave
core/thermal_diode_analysis.py. Its headline finding is a genuine,
unforced one: fluid dilution combined with this architecture's lack of
regeneration collapses the usable span to a fraction of a Kelvin up to a
couple of Kelvin at realistic particle loadings in this repo's own
model, dramatically less than solid AMR achieves at the same field/flow
-- reported directly rather than hidden, along with a secondary,
narrower finding (the fluid system's own COP_electrical can exceed solid
AMR's at that shared tiny span, though both trail the liquid-cooling and
vapor-compression baselines there) that this step's own writeup is
careful not to overstate.

Phase 21 (see ROADMAP.md) added core/baseline_cooling.py's
`passive_regenerator_augmentation()` / `augmented_regenerator_cop()`: a
PASSIVE (not actively magnetized/demagnetized) magnetic regenerator
loaded into a conventional gas cycle's own internal regenerator, whose
Curie-point heat-capacity anomaly (core/mce_material.py's own
`total_heat_capacity()`, already computed for every other analysis in
this repo) can raise regenerator effectiveness (reusing
core/thermal.py's existing `regenerator_effectiveness()`, given a new
backward-compatible `cp_solid` override) relative to the same material's
own lattice-only heat capacity -- cheap, per the plan's own framing,
because it recombines two things this repo already computes rather than
adding new physics or new benchmark data. New step 15
(core/passive_regenerator_analysis.py) compares Gd, Gd5Si2Ge2, and
La0.7Ca0.3MnO3 as candidate passive-regenerator materials at the
representative ASHRAE operating point: Gd (Tc=294K, inside the
[291.15K, 301.15K] window) gives the largest boost, while the other two
materials (Curie temperatures well outside the window) give none --
directly confirming the plan's own "alignment" framing in this repo's
own model. See that module's docstring for a real implementation
correction made along the way (an early version compared each
candidate's full heat capacity against one shared flat reference value
and thereby mostly ranked materials by bulk lattice heat capacity, not
magnetic alignment; the fix compares each material against its OWN
lattice-only baseline instead). The effectiveness-to-COP mapping itself
is an illustrative, literature-range-anchored ceiling (generic
internal-heat-exchanger COP-improvement figures from the refrigeration
literature), not a fitted or digitized coefficient -- Tishin & Spichkin
(2003)'s own passive-regenerator chapter (Ch. 11) could not be digitized
for this pass, same book-access limitation already documented for
Phase 20.
"""

import logging
import os
import time
import traceback
import contextlib

import numpy as np
import csv

from core.mce_material import GADOLINIUM
from core.baseline_cooling import (vapor_compression_cop, liquid_cooling_cop,
                                    elastocaloric_reference_cop)
from core import validation
from core import validation_system
from core import loss_model
from core import economics
from core import emissions
from core import cascade
from core.cascade import staged_baseline_result
from core import giant_mce_analysis
from core import material_family_comparison
from core import nanocomposite_material
from core import sensitivity
from core import rsm
from core import optimize as optimize_module
from core.thermal import regenerator_effectiveness
from core.first_order_mce import GD5SI2GE2_FIRST_ORDER
from core import giguere_validation
from core import geometry_analysis
from core import hypereg_analysis
from core import hysteresis_sensitivity
from core import thermal_diode_analysis
from core import magnet_geometry
from core import fluid_mce_analysis
from core import passive_regenerator_analysis
from core import inhomogeneous_broadening
from core import plots
from core import design_recommendations

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

    # Phase 23: a single STATIC literature reference point (not a function
    # of span -- see core/baseline_cooling.py's own Phase 23 honesty flag
    # for why), computed once and repeated on every row so it plots as a
    # flat comparison line/column alongside the span-dependent AMR/VCC/
    # liquid/Carnot figures, exactly the way Carnot is already a per-row
    # reference figure but elastocaloric here is deliberately NOT.
    elasto = elastocaloric_reference_cop()

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
            "Elastocaloric_COP_ref": round(elasto.COP_representative, 2),
        })

    with open(RESULTS_CSV, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=rows[0].keys())
        writer.writeheader()
        writer.writerows(rows)

    logger.info(f"{'span(K)':>8} {'AMR elec COP':>13} {'stages':>7} {'VCC COP':>9} {'Liquid COP':>11} {'Carnot':>8} {'Elasto ref':>10}")
    for r in rows:
        logger.info(f"{r['span_K']:>8} {r['AMR_COP_electrical']:>13} {r['AMR_n_stages']:>7} "
                    f"{r['VaporCompression_COP']:>9} {r['LiquidCooling_COP']:>11} {r['Carnot_COP']:>8} "
                    f"{r['Elastocaloric_COP_ref']:>10}")
    logger.info(f"Elastocaloric reference COP={elasto.COP_representative:.2f} "
                f"(range {elasto.COP_low:.1f}-{elasto.COP_high:.1f}), a flat literature "
                f"anchor repeated across every span -- NOT a span-dependent simulation. "
                f"Source: {elasto.source_note}")
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
    logger.info("")
    logger.info("Phase 22 item 3 note (qualitative only, no candidate priced above uses "
                "amorphous MCM data -- see core/economics.py's own section docstring):")
    logger.info(f"  {economics.amorphous_material_cost_performance_note()}")


def run_emissions(representative_row):
    """Step 6: emissions comparison fed with the *actual* COPs from the
    baseline sweep instead of the module's illustrative example numbers.

    Also logs the same comparison at emissions.FACILITY_SCALE_KW (1 MW):
    compare_emissions() is linear in capacity, so the technology ranking
    and ratio are identical at either scale, but the representative
    operating point's own capacity (~1.3 kW here) makes the absolute
    tCO2e/yr numbers small enough that the ratio is easy to lose next to
    them -- the facility-scale row makes the same conclusion easier to
    read (see emissions.py's own "Reporting scale" docstring note)."""
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

    logger.info(f"Same comparison at facility scale ({emissions.FACILITY_SCALE_KW:.0f} kW = "
                f"{emissions.FACILITY_SCALE_KW / 1000.0:.1f} MW, same COPs -- linear model, "
                "so this is the same ratio read at a scale where it's easier to see):")
    facility_rows = emissions.compare_emissions(
        emissions.FACILITY_SCALE_KW, amr_cop=amr_cop, vcc_cop=vcc_cop, liquid_cop=liquid_cop)
    for r in facility_rows:
        logger.info(f"{r.technology:<32} refrigerant={r.refrigerant_GWP_tCO2e_per_year:9.2f} "
                    f"tCO2e/yr  operational={r.operational_CO2_tCO2e_per_year:10.2f} tCO2e/yr  "
                    f"total={r.total_tCO2e_per_year:10.2f} tCO2e/yr")
    amr_total = facility_rows[0].total_tCO2e_per_year
    baseline_totals = [r.total_tCO2e_per_year for r in facility_rows[1:] if r.total_tCO2e_per_year > 0]
    if baseline_totals:
        best_baseline = min(baseline_totals)
        logger.info(f"Ratio (AMR total / best baseline total): {amr_total / best_baseline:.2f}x "
                    "-- identical at either scale by construction of the linear model.")


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
        # ROADMAP.md Phase 17 follow-up: closes that phase's own "did NOT
        # do" item (cycle_type was never threaded through cascade.py).
        # Astronautics_rotary_2014 is the one device with the largest COP
        # error on record, and is itself the naming-convention "rotary"
        # case infer_cycle_type_for_device() would flag -- so check
        # directly whether cycle_type="ericsson" (which narrowed
        # DTU_Eriksen_rotary_Gd_2015's error in step 2b) does the same
        # here, reusing the already-computed brayton `astro` result above
        # rather than recomputing it.
        ericsson = cascade.validate_astronautics_graded_bed(cycle_type="ericsson")
        if ericsson.get("feasible"):
            logger.info(f"Phase 17 follow-up: same graded bed under cycle_type='ericsson' "
                        f"(rotary-device naming-convention proxy, see validation_system."
                        f"infer_cycle_type_for_device()): COP={ericsson['COP_cascade']} "
                        f"({ericsson['COP_error_pct']:+.1f}% error) vs. brayton's "
                        f"{astro['COP_error_pct']:+.1f}% -- "
                        + ("ericsson narrows the error, a second data point in the same "
                           "direction as DTU_Eriksen_rotary_Gd_2015's step 2b result."
                           if abs(ericsson["COP_error_pct"]) < abs(astro["COP_error_pct"]) else
                           "ericsson does NOT narrow the error here, unlike "
                           "DTU_Eriksen_rotary_Gd_2015's step 2b result -- the remaining "
                           "gap is dominated by other documented issues (single-Tc-vs-6-"
                           "real-layers approximation, the ~2.4x DeltaT_ad overestimate in "
                           "giguere_validation.py), not by cycle topology."))
        astro["ericsson_followup"] = ericsson
    else:
        logger.info(astro.get("status", "infeasible"))
    return astro


def run_remaining_structural_devices_graded_validation():
    """Step 7d (calibration_failure_diagnostics.txt / Paper-Mining Pass
    review item 1, follow-up to step 7c): step 2c's diagnose_calibration_failure()
    classified 11/11 NO-CALIBRATION-FOUND benchmark rows as STRUCTURAL (span
    exceeds 2*dTad_noload for every mdot), not a search-space artifact -- see
    results/calibration_failure_diagnostics.txt. Step 7c already showed a
    6-layer Curie-graded La(Fe,Si)13Hy bed reproduces Astronautics_rotary_2014,
    the largest of these. This step extends the graded-bed structural fix to
    the other most data-center-relevant STRUCTURAL rows:

      - DTU_MagQueen_2018: genuinely a 10-layer graded bed per its own
        source note (like Astronautics) -- validate_magqueen_graded_bed().
        mass_MCM_kg is unreported, so run_magqueen_mass_sensitivity()
        sweeps it instead of assuming one value.
      - DTU_Eriksen_MAGGIE_2016: also genuinely graded, and unlike every
        other row in this set its real per-layer alloy compositions AND
        their measured Curie temperatures are directly reported (Eriksen
        et al. 2015 IJR paper + Eriksen 2016 PhD thesis, both now in this
        repo's Papers/) -- validate_maggie_real_graded_bed() uses those
        REAL Tc's directly (via GADOLINIUM.with_Tc()), not a composition
        search against a hypothetical tunable family. This closes the
        item explicitly deferred as future work in this function's
        previous pass.
      - Risoe_DTU_Gd_2011: real hardware is a single plain-Gd bed (NOT
        reported as graded) -- validate_risoe_dtu_graded_bed() tests the
        separate, explicitly-labeled HYPOTHETICAL question of whether a
        graded redesign could close the gap, not a reproduction claim.
      - Cooltech_2013_rotary: same hypothetical-redesign caveat as Risoe,
        plus it is a capacity-only row (no COP_lit) -- validate_cooltech_graded_bed()
        checks Qc feasibility only, and mass is again swept
        (run_cooltech_mass_sensitivity()) since it too is unreported."""
    out = {}
    logger.info("--- DTU_MagQueen_2018 (real 10-layer graded bed; mass unreported -> swept) ---")
    out["magqueen_mass_sensitivity"] = cascade.run_magqueen_mass_sensitivity()

    logger.info("--- DTU_Eriksen_MAGGIE_2016 (REAL 4-composition Gd/Gd-Y graded bed -- "
                "actual reported Curie temperatures, not a composition search) ---")
    maggie_result = cascade.run_maggie_span_sensitivity()
    out["maggie"] = maggie_result
    maggie = maggie_result["maggie"]
    companion = maggie_result["companion_2015"]
    if maggie.get("feasible"):
        logger.info(f"MAGGIE (15.5K span, 0.61Hz): Qc={maggie['Qc_W']}W (target "
                    f"{maggie['Qc_lit_W']}W), COP_cascade={maggie['COP_cascade']}  vs. "
                    f"reported COP={maggie['COP_lit']} ({maggie['COP_error_pct']:+.1f}% error) "
                    "-- vs. the flat 'no calibration found' the single-Tc Gd model gave this "
                    "row in step 2.")
    if companion.get("feasible"):
        logger.info(f"Companion DTU_Eriksen_rotary_Gd_2015 (10.2K span, 0.75Hz, SAME physical "
                    f"prototype) under the same real 4-layer model: COP_cascade="
                    f"{companion['COP_cascade']} vs. reported COP={companion['COP_lit']} "
                    f"({companion['COP_error_pct']:+.1f}% error) -- notably WORSE than step 2's "
                    "own single-Tc Gd approximation for this same row (-2.1% error there). The "
                    "real 4-layer model turns MAGGIE's own row from uncalibratable to "
                    "calibrated, but at this companion point it trades away accuracy the "
                    "simpler single-Tc approximation already had -- a genuine, stated trade-off, "
                    "not a strict improvement in every respect.")

    logger.info("--- Risoe_DTU_Gd_2011 (HYPOTHETICAL graded-redesign test -- real device is a "
                "single plain-Gd bed, not reported as graded) ---")
    risoe = cascade.validate_risoe_dtu_graded_bed()
    out["risoe_dtu"] = risoe
    if risoe.get("feasible"):
        logger.info(f"6-stage graded redesign: Qc={risoe['Qc_W']}W (target {risoe['Qc_lit_W']}W), "
                    f"COP_cascade={risoe['COP_cascade']}  vs. reported COP={risoe['COP_lit']} "
                    f"({risoe['COP_error_pct']:+.1f}% error)  -- "
                    f"{risoe['n_stages_out_of_range']}/6 stages fell back to plain Gd "
                    "(their needed Tc exceeded GD_FAMILY's documented 20-290K range). "
                    + ("Qc target reached but COP collapsed to 0 because at least one stage's "
                       "own Qc fell to 0 (span_fraction clamp) -- the structural Qc gap the "
                       "single-Tc model could not reach IS closed by splitting into 6 stages, "
                       "but a real efficient redesign at this span/field would need those "
                       "hot-side stages covered by an actually-tunable material this family's "
                       "documented range does not reach, not by plain Gd."
                       if risoe["COP_cascade"] == 0.0 else
                       "Unlike the pure-Qc-feasibility outcome elsewhere in this set, this "
                       "config also lands a non-zero cascade COP."))
    else:
        logger.info(risoe.get("status", "infeasible"))

    logger.info("--- Cooltech_2013_rotary (HYPOTHETICAL graded-redesign test, capacity-only "
                "row -- no COP_lit; mass also unreported -> swept) ---")
    out["cooltech_mass_sensitivity"] = cascade.run_cooltech_mass_sensitivity()

    logger.info("CONCLUSION: the graded-bed STRUCTURE (splitting one large span across several "
                "stages, each handling a small local span its own material can reach) closes "
                "the Qc-feasibility gap for every device checked here, at every swept mass -- "
                "confirming step 7c's Astronautics finding generalizes as a MECHANISM. Whether "
                "that translates into a genuinely EFFICIENT (positive-COP) design depends on "
                "whether the per-stage material is REAL or hypothetical: DTU_MagQueen_2018 and "
                "DTU_Eriksen_MAGGIE_2016 both use real, literature-reported graded compositions "
                "and both land non-zero, genuinely-informative COP errors (-92% and -69% "
                "respectively -- large, but real numbers now exist where none did before); "
                "Risoe_DTU's 6-stage HYPOTHETICAL redesign, by contrast, only gets 2/6 stages "
                "within GD_FAMILY's documented range at this device's span, so its COP result "
                "should be read as a Qc-feasibility finding, not an efficiency prediction. The "
                "MAGGIE companion check also surfaces a genuine trade-off: the real 4-layer "
                "model that newly calibrates MAGGIE's own 15.5K row makes the ALREADY-working "
                "10.2K companion row's COP error measurably worse (-2.1% -> -46.9%) -- a single "
                "graded-bed model does not uniformly improve every operating point of the same "
                "physical device, and that is reported directly rather than cherry-picked away.")
    return out


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


def run_eddy_and_pump_efficiency_demo():
    """Phase 27/28 demo. Neither core.thermal.intragranular_eddy_power()
    (Phase 27) nor AMRSystem.pump_motor_efficiency (Phase 28) has its own
    dedicated validation/analysis module (both are opt-in refinements to
    existing machinery, not new candidate materials or standalone
    analyses) -- this reproduces the same honest findings ROADMAP.md's
    own Phase 27/28 entries document, so they appear in a fresh
    `python main.py` run rather than only existing in this repo's own
    test suite and prior interactive session."""
    from core.thermal import intragranular_eddy_power
    from core.loss_model import StateDependentLossModel
    from core.amr_cycle import AMRSystem
    from core.mce_material import GADOLINIUM

    logger.info("Phase 27: geometry-explicit intragranular eddy-current loss vs. the "
                "CORE-calibrated support-structure k_eddy term, at 2T/2Hz")
    lm = StateDependentLossModel()
    k_eddy_term_W = lm.k_eddy * 2.0 ** 2 * 2.0 ** 2
    logger.info(f"  support-structure k_eddy term at 2T/2Hz: {k_eddy_term_W:.3f} W")
    for d_mm in (0.07, 0.17, 0.5, 100.0):
        W = intragranular_eddy_power(2.0, 2.0, particle_diameter=d_mm / 1000.0,
                                      mass_regenerator=2.0)
        ratio = W / k_eddy_term_W if k_eddy_term_W > 0 else float("nan")
        logger.info(f"  d_p={d_mm:7.2f}mm  intragranular eddy={W:10.6f} W  "
                    f"(ratio to support-structure term: {ratio:.6f})")
    logger.info("  -> negligible at realistic packed-bed particle sizes (0.07-0.17mm) -- "
                "the mechanism is real and wired in, but does not meaningfully change any "
                "existing result at sub-mm scale (see ROADMAP.md Phase 27).")

    logger.info("Phase 28: pump/motor efficiency opt-in -- idealized (default, "
                "efficiency=1.0) vs. literature-grounded (efficiency=0.6) at a "
                "representative geometry-explicit AMRSystem")
    sys_ideal = AMRSystem(material=GADOLINIUM, mu0H_max=2.0, mass_regenerator=2.0,
                           frequency=2.0, fluid_mdot=0.1, particle_diameter=0.0005,
                           use_ntu_thermal_model=True, loss_model=lm,
                           pump_motor_efficiency=1.0)
    sys_real = AMRSystem(material=GADOLINIUM, mu0H_max=2.0, mass_regenerator=2.0,
                          frequency=2.0, fluid_mdot=0.1, particle_diameter=0.0005,
                          use_ntu_thermal_model=True, loss_model=lm,
                          pump_motor_efficiency=AMRSystem.PUMP_MOTOR_EFFICIENCY_LITERATURE)
    r_ideal = sys_ideal.run(290.0, 10.0)
    r_real = sys_real.run(290.0, 10.0)
    logger.info(f"  idealized (eff=1.0):   COP_electrical={r_ideal.COP_electrical:.4f}")
    logger.info(f"  literature (eff=0.6):  COP_electrical={r_real.COP_electrical:.4f}")
    logger.info("  -> default AMRSystem/optimize.py behavior is UNCHANGED (eff=1.0); "
                "the literature value is available as an explicit opt-in, not silently "
                "applied to any production Pareto front (see ROADMAP.md Phase 28).")


def run_plot_generation(precomputed=None):
    """Step 12: renders all 34 figures in plots.py (results/figures/*.png
    and *.pdf) covering material validation, AMR characteristic curves,
    thermal/geometry modelling, loss-model calibration, system/curve
    validation, cascade and Curie-graded staging, Sobol sensitivity, RSM
    surrogate fitting, NSGA-III optimization, economics, emissions, and
    (figs 27-34) the Phase 16-22 sensitivity studies: Tc-broadening,
    nanocomposite off-design robustness, thermal-diode actuation cost,
    magnetocaloric-fluid volume fraction, passive-regenerator alignment,
    rotary-device cycle-type validation, hysteresis-loss and
    magnet-geometry Pareto-front sensitivity.
    Most figures still compute their own data directly from core/, but the
    eleven figures that duplicate an earlier stage's computation exactly
    (fig08 baseline sweep, fig14 system validation, fig16 Sobol, fig18
    Pareto, fig19/20 cascade, fig21 graded cascade, fig25 Astronautics
    validation, fig26 material family comparison, fig33 hysteresis A/B,
    fig34 magnet-geometry A/B) now reuse the results already produced by
    steps 2/4/7/7b/7c/8d/9/9b/11/11b/11d via `precomputed`, instead of
    re-running them from scratch. fig33/fig34 were previously the biggest
    offenders here — each re-ran a full ~9s NSGA-III optimization
    (fig34 ran it twice, once per FLAT/GEOMETRIC variant) that step
    11b/11d had already just computed one stage earlier; reusing
    `hysteresis_result`/`magnet_geometry_result` removes roughly 25-30s
    from a full `python main.py` run. This still runs last so the
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
                                          representative_span_K=REPRESENTATIVE_SPAN_K,
                                          cycle_type_result=None, thermal_diode_rows=None):
    """Step 13: consolidates the already-computed results from steps
    2b/3c/7b/8d/9b/11/11c into one ranked "how do I raise AMR electrical
    COP" report via core/design_recommendations.py. Pulls out the 10K-span,
    3-stage graded-vs-plain-Gd cascade rows (matching the comparison
    main.py already logs at the end of step 7b) since the graded/plain
    cascade sweeps return every span/stage-count combination, not just
    the representative one.

    Also folds in three findings that previously lived only in their own
    result files and never made it into this consolidated report: cycle-
    type validation (step 2b), thermal-diode cost sensitivity (step 11c),
    and passive/hybrid-regenerator augmentation (step 21/Phase 21) plus
    the static elastocaloric literature reference (Phase 23). cycle_type_
    result and thermal_diode_rows are passed in (steps 2b/11c both run
    before this one). The passive-regenerator and elastocaloric numbers
    are cheap (sub-second, no file I/O) to recompute directly here rather
    than plumbing them through from step 15, which runs AFTER this
    report is built -- same convention as step 3c's own cheap geometry
    re-sweep above."""
    graded_row = None
    gd_cascade_row = None
    if graded_rows is not None:
        graded_row = next((r for r in graded_rows if r["span_K"] == representative_span_K), None)
    if cascade_rows_gd is not None:
        gd_cascade_row = next((r for r in cascade_rows_gd if r["span_K"] == representative_span_K),
                               None)

    # Cheap (sub-second), non-printing, non-file-writing recomputation --
    # see the docstring note above for why these two are computed here
    # rather than passed in from their own (later-running, or text-only)
    # pipeline stages.
    passive_regen_base, passive_regen_rows = passive_regenerator_analysis.compare_candidate_materials(
        verbose=False)
    elastocaloric_result = elastocaloric_reference_cop()

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
        cycle_type_result=cycle_type_result,
        thermal_diode_rows=thermal_diode_rows,
        passive_regen_base=passive_regen_base,
        passive_regen_rows=passive_regen_rows,
        elastocaloric_result=elastocaloric_result,
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
                  validation.run_curie_shift_check(),
                  validation.run_curie_shift_check_v2())),
        ("1b. Inhomogeneous/polycrystalline Tc-broadening sensitivity "
         "(core/inhomogeneous_broadening.py, Phase 22 item 1)",
         lambda: inhomogeneous_broadening.run_inhomogeneous_broadening_analysis()),
        ("2. System-level validation vs. published AMR prototypes",
         None),  # handled specially below, result (system_validation_results) captured
        ("2b. Cycle-type (Ericsson-like vs. Brayton-like) validation sensitivity "
         "(core/validation_system.py, Phase 17)",
         lambda: validation_system.run_cycle_type_validation()),
        ("2c. Calibration-failure root-cause diagnostic: search-space artifact vs. "
         "structural limitation (core/validation_system.py, Paper-Mining Pass review item 1)",
         lambda: validation_system.run_calibration_failure_diagnostics()),
        ("3. Loss-model calibration (auto-loaded by AMRSystem's default loss model)",
         lambda: (loss_model.calibrate_loss_coefficients(), loss_model.run_extended_diagnostic())),
        ("3b. Regenerator thermal-effectiveness demo (core/thermal.py, reached transitively otherwise)",
         run_thermal_demo),
        ("3c. Geometry-dependent pumping power: packed-bed + parallel-plate (core/geometry_analysis.py)",
         lambda: geometry_analysis.run_geometry_analysis()),
        ("3d. Hypereg parallel-hydraulic pumping-power analysis (core/hypereg_analysis.py, Phase 15 item 3)",
         lambda: hypereg_analysis.run_hypereg_analysis()),
        ("3e. Geometry-explicit intragranular eddy-current loss + pump/motor "
         "efficiency demo (core/thermal.py, core/amr_cycle.py, Phase 27/28)",
         run_eddy_and_pump_efficiency_demo),
        ("4. Baseline comparison sweep: AMR vs VCC vs liquid cooling vs Carnot",
         None),  # handled specially below, result captured
        ("5. Economics / TCO at the representative operating point",
         None),  # needs step 4's result
        ("5b. Full-system cost estimate by material family (core/economics.py, Phase 15 item 5) + amorphous-material cost/performance note (Phase 22 item 3)",
         run_full_system_cost_by_material),
        ("6. Emissions comparison at the representative operating point",
         None),  # needs step 4's result
        ("7. Cascade staging comparison (1-4 stage AMR, Gd and Gd5Si2Ge2)",
         None),  # handled specially below, result (rows_gd) captured for step 7b
        ("7b. Curie-graded cascade (composition-tuned per stage, ROADMAP.md Phase 7 item)",
         None),  # needs step 7's result
        ("7c. Does a 6-layer Curie-graded La(Fe,Si)13Hy bed reproduce Astronautics_rotary_2014? (ROADMAP.md Phase 9 addendum)",
         run_astronautics_graded_validation),
        ("7d. Extending the graded-bed structural fix to the remaining STRUCTURAL "
         "devices (DTU_MagQueen_2018, Risoe_DTU_Gd_2011, Cooltech_2013_rotary) "
         "(calibration_failure_diagnostics.txt / Paper-Mining Pass review item 1 follow-up)",
         run_remaining_structural_devices_graded_validation),
        ("7e. Does the Giguere DeltaT_ad correction (giguere_validation.py, ~2.42x "
         "overestimate factor) narrow the Astronautics_rotary_2014 graded-bed "
         "-81.1% COP error? (core/cascade.py, Paper-Mining Pass review item 2)",
         lambda: cascade.run_astronautics_giguere_correction_sensitivity()),
        ("8. Giant-MCE materials analysis (Gd vs Gd5Si2Ge2)",
         lambda: giant_mce_analysis.run_analysis()),
        ("8d. Eight-way material family comparison (Gd, Gd5Si2Ge2-fixed, GD/LAFESIH/MNFEPSI/GA1XCMN3X/MNCUCOGE-tuned, LAFESIH-nanocomposite; Track A2 item + Phase 22 item 2 + Phase 24/25)",
         lambda: material_family_comparison.run_analysis()),
        ("8e. Nanocomposite off-design robustness check (core/nanocomposite_material.py, Phase 22 item 2 follow-up)",
         lambda: nanocomposite_material.run_robustness_check()),
        ("8b. First-order Landau model calibration check (core/first_order_mce.py, reached transitively otherwise)",
         run_first_order_mce_demo),
        ("8c. Giguere et al. (1999) direct-measurement cross-check + Pecharsky & "
         "Gschneidner (1997) peak-ratio check + latent-heat Cp spike (core/giguere_validation.py, "
         "Phase 26)",
         lambda: (giguere_validation.run_validation(),
                  giguere_validation.run_pecharsky_ratio_check(),
                  giguere_validation.run_latent_heat_validation())),
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
        ("11b. Hysteresis sensitivity: does Phase 16's thermal-hysteresis loss change the "
         "Phase 15 material-selection result? (core/hysteresis_sensitivity.py, Phase 16)",
         None),  # handled specially below, result (hysteresis_result) captured for the executive summary
        ("11c. Thermal-diode cost-only sensitivity (upper bound on switching-power "
         "overhead, NOT a net-benefit finding): mechanical-contact active thermal diode "
         "(core/thermal_diode.py, core/thermal_diode_analysis.py, Phase 18)",
         lambda: thermal_diode_analysis.run_thermal_diode_analysis()),
        ("11d. Magnet-geometry (Halbach-cylinder) field-vs-mass cost model "
         "(core/magnet_geometry.py, Phase 19)",
         None),  # handled specially below, result (magnet_geometry_result) captured
                 # for the executive summary
        ("11e. Magnet-geometry Pareto sensitivity: production-settings, multi-seed "
         "stability check (core/magnet_geometry.py, Paper-Mining Pass review item 4)",
         lambda: magnet_geometry.run_magnet_geometry_multiseed_stability_check()),
        ("11f. Layered/graded-bed NSGA-III co-optimization (core/optimize.py, "
         "core/cascade.py, Phase 29) -- reduced pop_size/n_gen/n_layers_range vs. "
         "the function's own full-quality defaults, purely to keep this pipeline "
         "stage's own runtime bounded (see run_layered_optimization()'s own "
         "docstring for the full 1-6 layer, pop_size=40/n_gen=25 version, callable "
         "directly for a dedicated deep run)",
         None),  # handled specially below, result (layered_pareto_rows) captured
                 # for the executive summary
        ("12. Figure generation: 34 figures covering validation, AMR curves, "
         "cascade/graded staging, sensitivity, RSM, NSGA-III, economics, emissions, "
         "Tc-broadening, nanocomposite robustness, thermal-diode, fluid-MCE, passive "
         "regenerator, cycle-type, hysteresis and magnet-geometry Pareto sensitivity (plots.py)",
         None),  # handled specially below, reuses steps 7/7b/9/9b/11's results
        ("13. Design-recommendations synthesis (core/design_recommendations.py)",
         None),  # handled specially below, consumes steps 3c/7b/8d/9b/11's results
        ("14. Magnetocaloric fluids (ferrofluid/MR suspension) as an alternative "
         "working-body class (core/fluid_mce_cycle.py, core/fluid_mce_analysis.py, "
         "Phase 20)",
         None),  # handled specially below, result (fluid_mce_result) captured
                 # for the executive summary
        ("15. Passive/hybrid magnetic regenerator augmentation of a conventional "
         "gas cycle (core/baseline_cooling.py, core/passive_regenerator_analysis.py, "
         "Phase 21)",
         None),  # handled specially below, result (passive_regen_result) captured
                 # for the executive summary
    ]

    representative_row = None
    system_validation_results = None
    cascade_rows_gd = None
    cascade_rows_giant = None
    graded_rows = None
    astro_result = None
    material_rows = None
    pareto_rows = None
    layered_pareto_rows = None
    sobol_const_Si = None
    sobol_state_dependent_Si = None
    pb_best_cop_row = None
    pp_best_cop_row = None
    hysteresis_result = None
    magnet_geometry_result = None
    magnet_geometry_multiseed_result = None
    fluid_mce_result = None
    passive_regen_result = None
    cycle_type_result = None
    thermal_diode_rows = None
    curie_shift_v2_result = None
    astronautics_giguere_result = None

    for name, fn in stages:
        _banner(name)
        t0 = time.time()
        try:
            with contextlib.redirect_stdout(_StreamToLogger(logger)):
                if name.startswith("1."):
                    curie_shift_v2_result = fn()[-1]  # (run_validation(), run_giguere_gd_extension(),
                                                        # run_curie_shift_check(), run_curie_shift_check_v2())
                elif name.startswith("2."):
                    system_validation_results = validation_system.run_system_validation()
                    validation_system.run_field_sensitivity_check()
                    validation_system.run_capacity_only_calibration_check()
                    # ROADMAP.md Group A completion pass: full 3-point Qc(span)
                    # curve-shape check for Tusek AMR(A) V*=0.95 (the one digitized
                    # curve whose anchor point calibrates -- see that function's
                    # docstring for why only this AMR/V* combo is checked here).
                    # Reads data/tusek_ate2013_figs/fig10_data.csv directly rather
                    # than going through the CSV benchmark rows, since
                    # run_curve_validation() (called separately by step 12's fig15)
                    # only ever compares one companion point, not the full
                    # digitized curve.
                    validation_system.run_tusek_multipoint_curve_validation()
                elif name.startswith("2b."):
                    cycle_type_result = validation_system.run_cycle_type_validation()
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
                elif name.startswith("7e."):
                    astronautics_giguere_result = fn()
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
                elif name.startswith("11b."):
                    hysteresis_result = hysteresis_sensitivity.run_hysteresis_sensitivity()
                elif name.startswith("11c."):
                    fn()
                    # Cheap (sub-second), non-printing re-sweep purely to capture
                    # the diode-vs-no-diode COP rows for step 13's synthesis
                    # report -- same convention as step 3c's geometry re-sweep
                    # above; run_thermal_diode_analysis() already did the
                    # printed, file-writing version of this same sweep.
                    thermal_diode_rows = thermal_diode_analysis.sweep_frequency_with_and_without_diode(
                        verbose=False)
                elif name.startswith("11e."):
                    magnet_geometry_multiseed_result = fn()
                elif name.startswith("11f."):
                    layered_pareto_rows = optimize_module.run_layered_optimization(
                        n_layers_range=(1, 2, 3), pop_size=20, n_gen=10, seed=1,
                        out_csv="results/layered_pareto_front.csv",
                        per_n_layers_out_dir="results/layered_pareto_front_by_n")
                elif name.startswith("11d."):
                    magnet_geometry.run_magnet_geometry_analysis()
                    magnet_geometry_result = magnet_geometry.run_geometric_cost_pareto_sensitivity()
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
                        "hysteresis_result": hysteresis_result,
                        "magnet_geometry_result": magnet_geometry_result,
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
                        cycle_type_result=cycle_type_result,
                        thermal_diode_rows=thermal_diode_rows,
                    )
                elif name.startswith("14."):
                    fluid_mce_result = fluid_mce_analysis.run_fluid_mce_analysis()
                elif name.startswith("15."):
                    passive_regen_result = passive_regenerator_analysis.run_passive_regenerator_analysis()
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
                     "hysteresis_sensitivity.txt (Phase 16), "
                     "cycle_type_validation.txt (Phase 17), "
                     "calibration_failure_diagnostics.txt (Paper-Mining Pass review item 1), "
                     "thermal_diode_analysis.txt (Phase 18), "
                     "magnet_geometry_analysis.txt, "
                     "magnet_geometry_pareto_sensitivity.txt, "
                     "pareto_front_magnet_flat.csv, pareto_front_magnet_geometric.csv "
                     "(Phase 19), "
                     "magnet_geometry_multiseed_stability.txt "
                     "(Paper-Mining Pass review item 4), "
                     "fluid_mce_analysis.txt (Phase 20), "
                     "passive_regenerator_analysis.txt (Phase 21), "
                     "geometry_optimization_analysis.txt, graded_cascade_comparison.csv, "
                     "design_recommendations.txt, figures/*.png+*.pdf (34 figures)")
    logger.info(f"Full run log: {LOG_FILE}")

    _print_executive_summary(representative_row, cascade_rows_gd, graded_rows,
                              material_rows, pareto_rows, pb_best_cop_row,
                              pp_best_cop_row, hysteresis_result,
                              magnet_geometry_result, fluid_mce_result,
                              passive_regen_result, failures, curie_shift_v2_result, 
                              astronautics_giguere_result, layered_pareto_rows, magnet_geometry_multiseed_result)


def _print_executive_summary(representative_row, cascade_rows_gd, graded_rows, material_rows,
                              pareto_rows, pb_best_cop_row, pp_best_cop_row,
                              hysteresis_result, magnet_geometry_result, fluid_mce_result,
                              passive_regen_result, failures, curie_shift_v2_result, 
                              astronautics_giguere_result, layered_pareto_rows, magnet_geometry_multiseed_result):
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
    if curie_shift_v2_result and _ok("1."):
        status = curie_shift_v2_result["status"]
        if status == "fit_succeeded":
            logger.info(f"  - Phenomenological Curie-shift patch (step 1, Paper-Mining Pass "
                        f"review item 3): fit SUCCEEDED, curie_shift_K_per_T="
                        f"{curie_shift_v2_result['curie_shift_K_per_T']:.2f} reproduces a "
                        f"{curie_shift_v2_result['fitted_slope_K_per_T']:.2f} K/T peak-shift "
                        f"rate on GADOLINIUM_FIELD_SHIFTED (plain GADOLINIUM unchanged)")
        else:
            logger.info("  - Phenomenological Curie-shift patch (step 1, Paper-Mining Pass "
                        "review item 3): implemented and tested, confirmed NOT to reproduce "
                        "Dan'kov et al.'s ~6 K/T shift for a specific, code-confirmed structural "
                        "reason (the H=0 reference entropy DeltaS_M(T,H) is measured against is "
                        "structurally blind to a field-only Tc shift) -- see this function's own "
                        "printed diagnostic above and core/validation.py's "
                        "calibrate_curie_shift()/diagnose_curie_shift_block() docstrings")
    if _ok("2c."):
        logger.info("  - Calibration-failure root cause (step 2c, Paper-Mining Pass review item "
                     "1): checks, per NO-CALIBRATION-FOUND row, whether span already exceeds "
                     "2*dTad_noload (structural -- no mdot bound can fix it) or whether a much "
                     "wider mdot search closes the gap (search-space artifact) -- see "
                     "results/calibration_failure_diagnostics.txt for the per-device breakdown "
                     "and which classification the largest devices (Astronautics, DTU MagQueen, "
                     "Risoe/DTU, Cooltech) fall into")

    logger.info("Baseline comparison (AMR vs. vapor-compression vs. liquid cooling vs. Carnot)")
    if representative_row and _ok("4."):
        logger.info(f"  - At {REPRESENTATIVE_SPAN_K:.0f}K span: AMR_COP_elec="
                    f"{representative_row['AMR_COP_electrical']}, VCC_COP="
                    f"{representative_row['VaporCompression_COP']}, Liquid_COP="
                    f"{representative_row['LiquidCooling_COP']}, Carnot_COP="
                    f"{representative_row['Carnot_COP']}  (results/comparison_table.csv)")
        logger.info(f"  - Elastocaloric (Phase 23, static literature reference, NOT "
                    f"span-simulated -- see core/baseline_cooling.py's own honesty flag): "
                    f"COP_ref={representative_row['Elastocaloric_COP_ref']}")
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
    if astronautics_giguere_result and _ok("7e."):
        base_err = astronautics_giguere_result["baseline"].get("COP_error_pct")
        corr_err = astronautics_giguere_result["corrected"].get("COP_error_pct")
        improves = astronautics_giguere_result["correction_improves"]
        logger.info(f"  - Astronautics_rotary_2014 Giguere-correction sensitivity (step 7e, "
                    f"Paper-Mining Pass review item 2): uncorrected COP_error={base_err}% -> "
                    f"Giguere-corrected COP_error={corr_err}% "
                    f"({'narrows the error' if improves else 'does NOT narrow the error'}) "
                    f"-- {'consistent with' if not improves else 'contrary to'} the existing "
                    "honesty flag that DTAD_CORRECTION_FACTOR was never shown to transfer from "
                    "Gd5Si2Ge2 to La(Fe,Si)13Hy; see core/cascade.py's "
                    "run_astronautics_giguere_correction_sensitivity() docstring")

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

    logger.info("Layered/graded-bed NSGA-III co-optimization (n_layers=1-3, reduced "
                "pop_size/n_gen for pipeline runtime; Phase 29)")
    if layered_pareto_rows and _ok("11f."):
        best_cascade_cop = max(layered_pareto_rows, key=lambda r: r["COP_cascade"])
        n_layers_counts = {}
        for r in layered_pareto_rows:
            n_layers_counts[r["n_layers"]] = n_layers_counts.get(r["n_layers"], 0) + 1
        logger.info(f"  - {len(layered_pareto_rows)} globally non-dominated designs across "
                    f"n_layers in {sorted(n_layers_counts)}; best cascade COP="
                    f"{best_cascade_cop['COP_cascade']} at n_layers={best_cascade_cop['n_layers']} "
                    f"(results/layered_pareto_front.csv, "
                    f"results/layered_pareto_front_by_n/n_layers_*.csv). n_layers "
                    f"representation: {dict(sorted(n_layers_counts.items()))} -- see "
                    f"run_layered_optimization()'s own docstring for the full 1-6 layer, "
                    f"production-settings version (not run here for pipeline-runtime reasons).")
    else:
        logger.info("  - unavailable (stage failed or was skipped)")

    logger.info("Hysteresis sensitivity: does Phase 16's thermal-hysteresis loss change the "
                "Phase 15 material-selection result? (step 11b, Phase 16)")
    if hysteresis_result and _ok("11b."):
        counts_on = hysteresis_result["counts_on"]
        counts_off = hysteresis_result["counts_off"]
        rows_on = hysteresis_result["rows_on"]
        rows_off = hysteresis_result["rows_off"]
        lafesih_on = sum(n for label, n in counts_on.items() if "La(Fe,Si)13Hy" in label)
        lafesih_off = sum(n for label, n in counts_off.items() if "La(Fe,Si)13Hy" in label)
        frac_on = lafesih_on / len(rows_on) if rows_on else 0.0
        frac_off = lafesih_off / len(rows_off) if rows_off else 0.0
        logger.info(f"  - La(Fe,Si)13Hy share of merged front: {frac_off:.0%} (hysteresis OFF, "
                    f"pre-Phase-16) -> {frac_on:.0%} (hysteresis ON, Phase 16) -- see "
                    "results/hysteresis_sensitivity.txt and that module's docstring honesty "
                    "flags 1-2 before treating this as a settled, publication-quality answer "
                    "rather than a directional sensitivity check")
    else:
        logger.info("  - unavailable (stage failed or was skipped)")

    logger.info("Thermal-diode cost-only sensitivity (upper bound on switching-power "
                "overhead, not a net-benefit finding): mechanical-contact active thermal "
                "diode (step 11c, Phase 18)")
    if _ok("11c."):
        logger.info("  - Cost-only, unbenchmarked design-exploration study (see "
                    "results/thermal_diode_analysis.txt and core/thermal_diode.py's docstring "
                    "honesty flag): confirms this repo's model has no internal mechanical-"
                    "switching frequency ceiling for a diode to relax, and quantifies the "
                    "small COP_electrical cost of an illustrative diode actuation-switching-"
                    "power term -- no benchmark device in this repo's corpus uses thermal "
                    "diodes, so this is not a validated feature")
    else:
        logger.info("  - unavailable (stage failed or was skipped)")

    logger.info("Magnet-geometry (Halbach-cylinder) field-vs-mass cost model "
                "(step 11d, Phase 19)")
    if magnet_geometry_result and _ok("11d."):
        rows_flat = magnet_geometry_result["rows_flat"]
        rows_geom = magnet_geometry_result["rows_geometric"]
        mean_flat = (sum(r["mu0H_max_T"] for r in rows_flat) / len(rows_flat)
                     if rows_flat else float("nan"))
        mean_geom = (sum(r["mu0H_max_T"] for r in rows_geom) / len(rows_geom)
                     if rows_geom else float("nan"))
        logger.info(f"  - Merged Pareto front mean mu0H_max_T: {mean_flat:.2f} T (flat "
                    f"per-Tesla magnet-mass ratio) -> {mean_geom:.2f} T (Phase 19 "
                    "geometric Halbach-cylinder magnet-mass relation) -- see "
                    "results/magnet_geometry_pareto_sensitivity.txt and "
                    "core/magnet_geometry.py's docstring honesty flags (incl. a citation "
                    "correction found in this pass) before treating this as a settled, "
                    "publication-quality answer rather than a directional sensitivity check")
    else:
        logger.info("  - unavailable (stage failed or was skipped)")

    if magnet_geometry_multiseed_result and _ok("11e."):
        per_seed = magnet_geometry_multiseed_result["per_seed"]
        stable = magnet_geometry_multiseed_result["stable"]
        seed_summary = "; ".join(
            f"seed{s['seed']}: {s['mean_flat_T']:.2f}T->{s['mean_geometric_T']:.2f}T"
            for s in per_seed)
        logger.info(f"  - Magnet-geometry Pareto sensitivity, production-settings multi-seed "
                    f"stability check (step 11e, Paper-Mining Pass review item 4): "
                    f"{'STABLE' if stable else 'NOT STABLE'} ({seed_summary}) -- "
                    f"{'the geometric cost term consistently pulls mean field down at every seed checked' if stable else 'at least one seed reverses direction, so the mean-field effect is seed-dependent / not reliably signed'} "
                    "-- see results/magnet_geometry_multiseed_stability.txt")

    logger.info("Consolidated design recommendations (step 13)")
    logger.info("  - Ranks all COP-maximization levers above by demonstrated Sobol "
                "sensitivity and reports a recommended starting design point "
                "(results/design_recommendations.txt)")

    logger.info("Figures")
    n_figs = len(list(plots.FIG_DIR.glob("*.png"))) if plots.FIG_DIR.exists() else 0
    logger.info(f"  - {n_figs} figure(s) generated (results/figures/*.png, *.pdf)")

    logger.info("Magnetocaloric fluids (ferrofluid/MR suspension) as an alternative "
                "working-body class (step 14, Phase 20)")
    if fluid_mce_result and _ok("14."):
        sweep = fluid_mce_result["sweep"]
        comp = fluid_mce_result["comparison"]
        best = sweep["best_row"]
        interior_note = ("a genuine interior optimum" if sweep["interior_optimum_found"]
                          else "no interior optimum (monotonic over the swept range)")
        logger.info(f"  - Volume-fraction sweep: COP_electrical best at phi={best['phi']:.2f} "
                    f"({interior_note}) -- results/fluid_mce_analysis.txt")
        logger.info(f"  - At its own favorable span ({comp['span_K']:.2f} K, phi="
                    f"{comp['phi']:.2f}): ferrofluid Qc={comp['fluid_MCE']['Qc_W']:.2f}W "
                    f"vs. solid AMR Qc={comp['solid_AMR']['Qc_W']:.2f}W -- fluid dilution "
                    "plus this architecture's lack of regeneration (see "
                    "core/fluid_mce_cycle.py's own docstring) collapses usable span far "
                    "below solid AMR's, a design-exploration finding, not a validated "
                    "benchmark-backed result (no fluid-MCE-as-working-body device was "
                    "found in this pass's own literature search -- see that module's "
                    "HONESTY FLAG #2)")
    else:
        logger.info("  - unavailable (stage failed or was skipped)")

    logger.info("Passive/hybrid magnetic regenerator: augmenting a conventional gas "
                "cycle's own regenerator with a magnetocaloric material's Curie-point "
                "heat-capacity anomaly (step 15, Phase 21)")
    if passive_regen_result and _ok("15."):
        best = passive_regen_result["candidate_results"][0]
        base_cop = passive_regen_result["base"].COP
        logger.info(f"  - Best candidate at the representative operating point: "
                    f"{best.material_name}  eps {best.eps_baseline:.3f} -> "
                    f"{best.eps_augmented:.3f}  COP {base_cop:.3f} -> "
                    f"{best.augmented_COP:.3f} ({best.cop_gain_fraction:+.2%}) -- "
                    "results/passive_regenerator_analysis.txt")
        logger.info("  - Alignment effect confirmed directly (not assumed): candidate "
                    "materials whose own Curie temperature falls outside the operating "
                    "window show ~0% gain, since delta_eps is clipped at 0 by "
                    "construction (see core/baseline_cooling.py's own docstring). The "
                    "effectiveness-to-COP mapping is an illustrative, literature-range-"
                    "anchored ceiling, not a fitted or digitized coefficient -- see that "
                    "module's Phase 21 honesty flag before treating this as a validated "
                    "device-level COP prediction")
    else:
        logger.info("  - unavailable (stage failed or was skipped)")


if __name__ == "__main__":
    main()