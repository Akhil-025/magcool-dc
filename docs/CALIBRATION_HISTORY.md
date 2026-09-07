# Calibration issues: history, fixes, and open items

This consolidates what was previously spread across six separate,
date-stamped working notes (`CALIBRATION_ERROR_ROOT_CAUSE_ANALYSIS.md`,
`CALIBRATION_ISSUES_AND_FIXES_2026-09-02.md`,
`CALIBRATION_ISSUES_STATUS_2026-09-02_PART2.md`/`_PART3.md`,
`WHY_NOT_ALL_FIXED_2026-09-02.md`, and `PAPERS_REVIEW_AND_PLAN_2026-09-02.md`)
into one document. See `LIMITATIONS.md` for the full technical derivations
behind each item below; this is the consolidated summary.

## Genuinely fixed

**Loss-accounting bug in graded/layered-bed reproductions.**
`core/cascade.py`'s `run_explicit_material_cascade()` / `run_graded_cascade()`
(and the three functions that reproduce real multi-layer devices:
`validate_maggie_real_graded_bed()`, `validate_astronautics_graded_bed()`,
`validate_magqueen_graded_bed()`) modeled a real single-bed device — one
set of magnets, one pump, one motor, N Curie-graded material layers packed
in one housing — as N *independent* full AMR systems, each separately
paying its own eddy, pumping, and baseline-overhead loss and summing all N
totals. A real graded bed only has one of each of those shared overhead
sources. The error got monotonically worse with layer count, a clean
fingerprint of an N-times-overcounted shared term:

| Device | Layers | COP error (before) | COP error (after) |
|---|---|---|---|
| MAGGIE | 4 | −69.3% | **−16.0%** |
| Astronautics | 6 | −81.1% | **+0.9%** |
| MagQueen | 10 | −92.2% | −47.4% |

Fix: a `shared_hardware` parameter. Magnetic work is still summed per
stage (legitimate — different Curie-shifted layers really do require
different thermodynamic work), but parasitic overhead is computed once
from the device's own aggregate operating point instead of once per
layer. The three real-device validation functions use it; the default
stays `False` so nothing else in the codebase changes. MagQueen's
residual −47% isn't attributable to this bug — its CSV row already flags
unreported mass (assumed) and derived (not directly measured) Qc/COP, a
data-provenance issue rather than a code bug.

**Qc(span) feasibility-reopening artifact (+787% on the Tušek case).**
`AMRSystem.cooling_capacity_span_sweep()` now applies a running-minimum
monotonicity clamp. Confirmed against the specific case that surfaced it.

**Hysteresis-share seed instability.** Re-run at n=20 seeds; the earlier
apparent reversal is statistically indistinguishable from zero, so this
is closed as resolved rather than merely re-confirmed as noisy.

**Tušek benchmark used the wrong regenerator correlation.** The
Tusek_singlebed_Gd_2010_spanceiling device is a parallel-plate
regenerator (0.1mm plate spacing, 0.25mm plate thickness, porosity
0.2564, per Tušek, Kitanovski, Zupan, Prebil, Poredoš, *Appl. Therm.
Eng.* 53 (2013) 57-66, Table 1), but `core/regenerator_1d.py` was calling
the packed-bed correlation with packed-bed default geometry for every
device regardless of actual type. `core/thermal.py` already had a real,
literature-sourced parallel-plate correlation
(`regenerator_effectiveness_parallel_plate()`, Nickolay & Martin 2002 /
Tušek et al. 2013 Eq. 4); it just wasn't wired into the transient
simulator. Added a `geometry="parallel_plate"` option (default stays
`"packed_bed"`, so every existing caller is unaffected) and re-ran with
the device's real, source-verified geometry: peak span went from ~1.4K
(−92% to −97% error) to ~16.5K (−16% error), using an existing
correlation and real device dimensions pulled from the source paper, not
a tuned parameter. A separate, small (~11%) discrepancy between the
model's idealized porosity (0.2857) and the paper's directly-measured
value (0.2564) is noted in the code as a real but minor remaining
imprecision.

**Gd near-Tc overprediction magnitude.** Three combined, literature-sourced
fixes: the exact isentropic ΔT_ad definition (superseding the small-field-
step linear approximation), a Sommerfeld electronic entropy term, and a
fitted polycrystalline grain-Curie-temperature broadening. Dan'kov et al.
errors went from +48.9%/+29.2%/+9.8% to +4.2%/−1.1%/−8.1%. The residual
few-percent error is consistent with genuine short-range-correlation
physics that mean-field theory omits (see below).

## Investigated, not resolved — and why

**Curie-point field-shift rate.** Reported literature rate is ~6 K/T; the
model predicts ~0. Two literature-grounded fix attempts were tried and
ruled out with evidence rather than abandoned without cause: an Oguchi
pair-cluster correction (doesn't beat mean-field on real data), and a
critical-exponent "universal curve" scaling (Dan'kov et al.'s own 3 points
imply a field-scaling exponent n≈0.95, and neither candidate theory —
mean-field n=0.667, 3D Heisenberg n=0.637 — is close). Closing this needs
primary data this repo doesn't have (a full digitized Dan'kov ΔT_ad(T)
curve, or their exact sample geometry for a demagnetization check), not a
third fit to the same three points. Reference/exploratory code:
`research/oguchi_pair_cluster_prototype.py`.

**Field-dependent Tc-broadening.** Motivated by a direct re-read of
Dan'kov, Tishin, Pecharsky & Gschneidner (1998), whose own prose (not just
its three calibration numbers) states that field broadens *and* shifts the
λ-type heat-capacity peak. The shift half is already modeled and already
shown unable to reproduce the reported rate (see above). The broadening
half is different — and grows with field, unlike the earlier constant-
sigma_Tc sweep, which correctly found a constant broadening can't fix a
field-dependent error pattern. Fitting `sigma_Tc(mu0*H) = k * mu0*H`
against the same three Dan'kov calibration points converges to k=0.000: no
broadening at all is preferred, because broadening only ever lowers
DeltaT_ad near the peak, and the sharp model's error is already small at
5T (+9.8%) and large at 1T (+48.9%) — any k>0 overcorrects the small error
long before making a dent in the large one. The underlying physics (field
broadens the peak) is real and literature-sourced; a linear-in-H form
calibrated to 3 points simply doesn't have the right shape to exploit it.
Reported as a negative result rather than searched until something looked
better. Nothing in the default pipeline uses this class — it's an
explicit, opt-in diagnostic (`core/inhomogeneous_broadening.py`'s
`FieldBroadenedMagnetocaloricMaterial`), unimported elsewhere.

**No regenerative amplification in the 0-D model.** `AMRSystem.
cooling_capacity()` caps achievable span at `2 x dTad_noload` — the
material's own single-blow adiabatic ΔT — which has nothing to do with the
regenerator mechanism that lets real devices reach much larger spans by
building up a temperature gradient over many cycles. Cross-checked against
literature independent of the benchmark set's own failures: the DTU
"MAGGIE" thesis reports a directly-measured 29.2 K no-load span for the
same 1.13 T / 1.7 kg hardware whose other two operating points are already
in this repo's CSV, versus the model's own best-case structural ceiling of
14.1 K at that field — a 2.1x gap on a data point that needed no `mdot`
back-calculation at all. `analyze_regenerative_amplification_gap()`
(`core/validation_system.py`, wired into `main.py` as step 2d) makes the
gap's size measurable across the whole benchmark set: 10 of 17 rows with a
well-defined `dTad_noload` exceed the model's own span cap, with a ratio
(actual span / model ceiling) of 1.04x-13.97x, median 1.63x; the
independently-sourced MAGGIE point lands at 4.39x, in the middle of that
range — a genuine cross-check, since it didn't feed the ratio statistics.
A real fix needs either an NTU/utilization-based semi-analytical formula
or a transient 1-D blow-by-blow AMR solver, both standard in the
numerical-modeling literature already in this repo's `Papers/` folder —
not a patch to the single-blow formula. `core/regenerator_1d.py`'s 1-D
transient simulator is a real step in that direction (see its own module
docstring for the still-open discretization-convergence question); it
remains diagnostic-only by default rather than wired into
`cooling_capacity()`, since forcing it in before it's grid-converged would
trade an honest "not yet" for a confidently wrong number.

**Regenerator grid convergence.** Confirming the parallel-plate geometry
fix above wasn't a fluke required checking convergence: running n_nodes=40
for 12,000 cycles (instead of 1,500) shows span settling at 1.04K, flat
over the last 10 cycles — not converging toward the n_nodes=20 answer
(16.5K) but to a genuinely different, smaller steady state. That rules out
"needs more time to settle" and confirms a real discretization
inconsistency. Four hypotheses were tested and eliminated (operator-
splitting order, axial-conduction magnitude, NTU magnitude, grid
resolution/cycle count); root cause remains open. A real fix needs a
properly re-derived, grid-independent spatial discretization for the
regenerator's fluid-solid coupling, verified against a known analytical or
reference solution — not a different formula chosen because the benchmark
number looks better. Reference/exploratory code:
`research/regenerator_1d_v2_coupled_prototype.py`.

**Loss-model device-data gaps (leave-one-out error 333%-1639% depending on
device).** `StateDependentLossModel`'s three coefficients are fit to
exactly three CORE calibration points — an exactly-determined system, not
a real regression, so a held-out point can miss by hundreds of percent.
This explains Tušek's (one of the three fitting points itself, flagged as
possibly using an inflated field value) −76% COP error and Okamura's
(whose frequency is a guessed 1.0 Hz placeholder feeding the f²-scaling
eddy term) +50% error. Closing this needs new experimental devices with
independently-reported component-level efficiency breakdowns; a review of
this repo's paper corpus (below) found none exist there. No amount of
further modeling effort manufactures data that doesn't exist.

**Opposite-direction material-level error.** `first_order_mce.py`'s Landau
model overestimates first-order materials' single-blow ΔT_ad by ~2.4x
versus Giguère et al.'s direct measurement (already caught by
`giguere_validation.py`). This doesn't reliably cancel the near-Tc issue
above — different mechanism, different level of the model.

## Paper-mining review (Papers/ corpus)

A full re-inventory of the `Papers/` folder (42 PDFs, 40 distinct papers,
cross-checked against `Literature_Review.md`, `CITATION_AUDIT_PHASE30.md`,
and every citation string already in `core/*.py`) found that 37 of the 40
papers were already cited and incorporated by earlier work. Two
housekeeping notes: one exact duplicate PDF under a misleading filename
(the "nine-layer" file is actually a byte-for-byte copy of the Astronautics
paper already used as the `Astronautics_rotary_2014` calibration point),
and one unrelated 1997 semiconductor-physics paper that appears to be an
accidental include and wasn't used.

Only three papers were genuinely new territory. De Oliveira & von Ranke,
*Phys. Rep.* 489 (2010) 89-159 independently confirms, in its own Gd
calculation, that mean-field theory's omission of short-range spin
correlations is what drives the near-Tc overprediction and the
Curie-shift-rate mismatch above — motivating the field-dependent
broadening investigation. Tishin & Spichkin (2003) remains fully
image-based (no extractable text layer) wherever it's cited as a gap.
Law, Moreno-Ramírez, Díaz-García & Franco, *J. Appl. Phys.* 133, 040903
(2023) was reviewed and did not change any conclusion above.

## What this means practically

Every number `results/` and `main.py`'s pipeline actually report reflects
the fixed items above — nothing unresolved is silently smoothed over.
The remaining open items are real, honestly characterized, backed by
ruled-out alternatives rather than "not tried," and the exploratory code
written while investigating them is preserved under `research/` rather
than discarded. Closing what's left needs one of two substantial,
multi-day undertakings, each with its own validation cycle:

1. An Oguchi/Bethe-Peierls pair-cluster model for Gd — a first-principles
   attempt at the near-Tc and Curie-shift issues, validated against
   Dan'kov and Giguère from the start.
2. A properly re-derived, grid-independent regenerator discretization —
   checked against a known reference case before trusting it on Tušek.

Neither is a same-session patch; both would need to be built with their
own validation cycle and their own honesty flags, the same way every other
addition in this repository has been.
