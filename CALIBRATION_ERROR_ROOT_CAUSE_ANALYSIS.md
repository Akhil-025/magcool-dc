# Why the calibration/COP errors happen — root-cause analysis and fixes

Scope: `Mag/` codebase, `results/` outputs, and the cited literature in `Papers/`.
This covers two distinct, independently-confirmed problems, one now fixed in code,
one quantified but left open (with a clear reason why).

---

## 1. FIXED: a loss-accounting bug that was wrecking COP for every graded/layered-bed reproduction

**Where:** `core/cascade.py` — `run_explicit_material_cascade()`, `run_graded_cascade()`,
and the three functions that call them to reproduce real multi-layer devices:
`validate_maggie_real_graded_bed()`, `validate_astronautics_graded_bed()`,
`validate_magqueen_graded_bed()`.

**What was wrong:** these functions model a real single-bed device (one set of magnets,
one pump, one motor, N Curie-graded material layers packed inside) as if it were N
*independent* full AMR systems, each separately computing and paying its own
`W_eddy = k_eddy·f²·mu0H²`, `W_pump = k_pump·mdot²`, and `W_base = base_frac·Qc`, then
summing all N stages' totals. A real graded bed only has **one** of each of those
overhead sources — the eddy loss and baseline electronics overhead don't multiply by
the number of material layers packed into one housing, and the fluid flow rate
(`fluid_mdot`) was already being passed identically, unscaled, to every stage, meaning
it was always intended to represent one shared flow path.

**Evidence this was a real bug, not just "the model is imprecise":** COP error got
monotonically worse as the number of layers increased — a clean fingerprint of an
N-times-overcounted shared term:

| Device | Layers | COP error (before) |
|---|---|---|
| MAGGIE | 4 | −69.3% |
| Astronautics | 6 | −81.1% |
| MagQueen | 10 | −92.2% |

**The fix:** added a `shared_hardware` parameter. When `True`, magnetic work
(`W_mag`) is still summed per stage (that's legitimate — different Curie-shifted
layers really do require different thermodynamic work), but the parasitic overhead
is computed **once**, from the device's own aggregate `(frequency, mu0H, mdot,
Qc_target)`, instead of once per layer. Wired the three real-device validation
functions to use it (they're reproducing one physical device each, so this is
unambiguous). Left the *default* at `False`, so nothing else in the codebase
changes — `run_cascade()`/`compare_staging()` (genuinely-separate cascade
refrigeration modules) and the exploratory `compare_graded_cascade()` /
`design_recommendations.txt` lever-3 numbers are untouched.

**Result after the fix:**

| Device | Before | After |
|---|---|---|
| MAGGIE | −69.3% | **−16.0%** |
| Astronautics | −81.1% | **+0.9%** |
| MagQueen | −92.2% | −47.4% |

Astronautics and MAGGIE are now solidly in line with the rest of the benchmark
set's accuracy. MagQueen's residual −47% is not attributable to this bug (it barely
moves under `cycle_type="ericsson"` either, −44.6%) — the CSV row's own notes already
flag that MagQueen's mass is unreported (assumed), and its Qc/COP are *derived*, not
directly measured, from the paper's heat-pump-mode numbers. That's a data-provenance
problem, not a code bug, and I didn't try to paper over it with a further correction.

All 27 existing `tests/test_cascade.py` tests still pass, plus `tests/test_validation.py`.

A side effect worth knowing about: `run_astronautics_giguere_correction_sensitivity()`
(the function that tests whether the Giguère `ΔT_ad` correction narrows Astronautics'
error) now starts from the fixed +0.9% baseline instead of the old −81.1%, and
applying the correction moves it to −0.3%. I updated that function's and the
corresponding `main.py` step's stale docstring text to say so.

---

## 2. QUANTIFIED, NOT FIXED: the core AMR model has no regenerative amplification

**Where:** `core/amr_cycle.py` — `AMRSystem.cooling_capacity()`.

**What's wrong:** the model caps the achievable temperature span at
`2 × dTad_noload`, where `dTad_noload` is the *magnetocaloric material's own
single-blow adiabatic ΔT* at the operating field and mid-bed temperature. That's
physically just the temperature change from one magnetize/demagnetize pulse — it
has nothing to do with the *regenerator* part of "Active Magnetic Regenerator,"
which is what actually lets real devices reach large spans: the packed bed builds
up a temperature gradient over many cycles, similar to a counter-flow cascade. A
material with only a 2–4 K single-blow ΔT can support a 20–30 K device span this
way. The current model can't represent that mechanism at all, by construction — it's
a single-node (0-D), single-blow evaluation, not a spatial/temporal regenerator
simulation.

I confirmed this independently against your own cited literature, not just the
benchmark set's own failures: the DTU "MAGGIE" thesis (already in `Papers/`) reports
a directly-measured **29.2 K no-load span** for the same 1.13 T / 1.7 kg hardware
whose two other operating points are already in your CSV. Your model's own best-case
structural ceiling at that field is **14.1 K** — a 2.1× gap, on a data point that
required no `mdot` back-calculation at all (unlike almost every other row). I added
this as a new CSV row, `DTU_Eriksen_MAGGIE_2016_noloadspan`, since it was flagged as
available-but-unused in an existing row's note.

I added `analyze_regenerative_amplification_gap()` to `core/validation_system.py`
(writes `results/regenerative_amplification_diagnostic.txt`, and is now wired into
`main.py` as step 2d) to make the size of this gap visible and reproducible across
every benchmark row, not just the ones that already fail to calibrate:

- **10 of 17** rows with a well-defined `dTad_noload` exceed the model's own
  structural span cap.
- Ratio (`actual span / model's own 2×dTad_noload ceiling`) ranges **1.04×–13.97×**,
  median **1.63×**.
- The independently-sourced MAGGIE no-load point lands at **4.39×** — right in the
  middle of that range, which is a genuine cross-check, not a circular one, since it
  didn't come from any calibration search.
- Two rows (`Astronautics_rotary_2014`, `DTU_MagQueen_2018`) show `dTad_noload ≈ 0`
  — a *different*, compounding problem: their fixed `T_COLD_LAFESIH_K` assumption
  puts the mid-bed temperature off the single-composition La(Fe,Si)₁₃Hy material's
  real MCE peak. These two need to be read separately from the "clean" ratio table,
  which I did (flagged as "near-zero / ratio undefined" rather than folded into the
  statistics as spuriously huge numbers).

**Why I didn't try to fix the formula itself:** your own `ROADMAP.md` (item A3)
already considered this territory and explicitly declined to invent an unsourced
smoothing function for the same reason I'd run into — no literature source in this
corpus gives the exact functional form for how span scales with NTU/utilization/
material ΔT_ad in a real regenerator. (Their A3 decision was actually about a
narrower question — smoothing the linear formula's corner — not the ceiling's scale;
I want to be precise that this is a related but larger open item, not a re-litigation
of a decision they already made correctly on its own terms.) The technically correct
fix is a real regenerator model — either an NTU/utilization-based semi-analytical
formula or a transient 1-D blow-by-blow AMR solver, both standard in the numerical-
modeling papers already in your `Papers/AMR Theory and Modeling/` folder — not a
patch to the existing single-blow formula. That's a substantial modeling project,
not a same-pass fix, so I scoped this contribution to making the gap's size
rigorously measurable instead of guessing at a replacement formula.

**Downstream reach:** every part of the pipeline that calls `cooling_capacity()`
inherits this — `comparison_table.csv`, both `cascade_comparison*.csv` files,
`design_recommendations.txt`, the NSGA-III Pareto front, and the Sobol sensitivity
analysis. Fix #1 (shared hardware) repairs the *work/COP* side for the three graded
real-device reproductions specifically; this section is about the *capacity/span*
side of the single-stage kernel underneath literally everything else in the
codebase, graded-bed or not.

---

## Secondary, smaller contributors (found, not changed)

- **Loss model has zero degrees of freedom.** `StateDependentLossModel`'s three
  coefficients are fit to exactly three CORE calibration points — an exactly-
  determined system, not a real regression. Your own leave-one-out check already
  shows this can miss a held-out point by ~680%. This explains why Tušek (one of
  the three fitting points itself, and flagged as possibly using an inflated field
  value) shows −76% COP error, and why Okamura (whose frequency is a guessed 1.0 Hz
  placeholder feeding the f²-scaling eddy term) shows +50%. Both are pre-existing,
  already-documented issues; I didn't touch the loss-model fit itself.
- **Opposite-direction material-level error.** `first_order_mce.py`'s Landau model
  overestimates first-order materials' single-blow ΔT_ad by ~2.4× versus Giguère et
  al.'s direct measurement (already caught in your own `giguere_validation.py`).
  This doesn't reliably cancel problem #2 above — it's a different mechanism at a
  different level of the model.

---

## Files changed / added

- `core/cascade.py` — `shared_hardware` fix (problem #1).
- `core/validation_system.py` — new `analyze_regenerative_amplification_gap()`
  (problem #2, diagnostic only).
- `main.py` — wired the new diagnostic in as step 2d; updated step 7e's and
  `run_astronautics_giguere_correction_sensitivity()`'s docstrings to stop quoting
  the now-superseded −81.1% baseline.
- `data/amr_experimental_benchmarks.csv` — added the MAGGIE 29.2 K no-load-span row.
- `results/regenerative_amplification_diagnostic.txt` — new output from the
  diagnostic above.
