# Known Limitations & Open Items -- Consolidated Ledger

This file exists for one reason: this project's own convention is to flag
every limitation, calibration gap, and unresolved finding explicitly rather
than smooth it over -- but by Phase 31 those flags were scattered across
README.md, ROADMAP.md, several dedicated `*_ROOT_CAUSE_ANALYSIS.md` /
`CITATION_AUDIT_*.md` files, module docstrings, and individual
`results/*.txt` diagnostic outputs. A judge, reviewer, or new contributor
asking "what doesn't work yet?" had to hunt across a dozen files to get a
complete answer. This file is that answer, in one place, with a pointer to
the authoritative source for each item. It does not replace those sources
(which contain the full diagnosis and reasoning) -- it is an index into them.

Nothing below is new information; every item was already flagged somewhere
in the repository before this file existed. Several were verified,
extended, or fixed directly during Phase 31 (noted inline).

---

## 1. Model-accuracy limitations (physics)

### 1.1 Mean-field Gd model overpredicts near the Curie point
Mean-field (Brillouin/Weiss) theory systematically overpredicts Gd's
adiabatic temperature change close to Tc=294K, worst at low field and
shrinking as field increases (+48.9% at 1T, +29.2% at 2T, +9.8% at 5T vs.
Dan'kov et al. 1998). This is a textbook, well-known limitation of
mean-field theory near a critical point (it ignores short-range spin
correlations) -- not a bug, and not fixable without a materially
different (e.g. renormalization-group) theory.
- Phase 35 correction: the 5T literature reference value was previously
  wrong (14.6K, making the model APPEAR to underestimate by -7.5% at 5T
  only -- an odd sign flip inconsistent with the theory above). A
  pixel-calibrated re-digitization of Dan'kov et al.'s own Fig. 10 found
  the correct value is 12.3K, not 14.6K -- see Section 5 below for the
  full methodology and cross-validation. The corrected table (all three
  fields overestimate, shrinking monotonically with field) is now
  internally consistent with the stated mean-field-near-Tc theory.
- Source: `core/mce_material.py`, `core/validation.py`'s
  `LITERATURE_DELTA_T_AD` (with the Phase 35 digitization note inline).

### 1.2 First-order Landau model overestimates giant-MCE ΔTad
The 6th-order Landau free-energy expansion used for Gd5Si2Ge2,
La(Fe,Si)13Hy, and (Mn,Fe)2(P,Si) overestimates ΔTad by roughly 2.4x
against Giguere et al.'s (1999) direct measurement.
- Source: `core/giguere_validation.py`; `tests/test_giant_mce_analysis.py`.
- Status: open. Re-fitting the 6th-order coefficients against Giguere's
  data specifically (rather than Pecharsky & Gschneidner's) is the
  concretely-scoped next step, not attempted in Phase 31.

### 1.3 1-D regenerator: direction-inconsistent error, even after two fix passes
`core/regenerator_1d.py`'s multi-cycle transient model undershoots on two
of three benchmark devices and overshoots on the third:
- Tusek: -96.9%, Lozano: +111.1%, DTU/MAGGIE: -62.0%
  (see `results/regenerator_1d_validation.txt` for the full table).
- Phase 31 fix attempted: replaced the previous ad hoc axial-conductivity
  multiplier with a real, citable Maxwell-Eucken packed-bed
  composite-conductivity model (`core/regenerator_1d.py`'s
  `_packed_bed_effective_axial_conductivity()`). **Honest result: this did
  NOT improve accuracy** -- it made the undershoot slightly WORSE on two
  rows while leaving the overshooting row essentially unchanged, because
  the more rigorously-justified conductivity value is *higher* than the
  old placeholder at these bed porosities (more damping, not less). This
  demonstrates the remaining gap is not simply "the conductivity constant
  is a bit off" -- a legitimately different, textbook-derived constant
  moves the same two rows further the same wrong way. The
  regenerator-effectiveness/NTU coupling or the single-blow reference
  itself are the more likely remaining suspects.
- Status: open, and NOT wired into `core/amr_cycle.py`'s
  `cooling_capacity()` (the function every other part of this codebase
  actually uses) until resolved.
- Source: `core/regenerator_1d.py` module docstring ("Known limitations"),
  `validate_against_benchmarks()`'s own printed/written caveat text.

### 1.4 Loss-model calibration has (at best) zero degrees of freedom
`StateDependentLossModel`'s 3 coefficients (k_eddy, k_pump, base_frac) are
fit to exactly 3 CORE calibration points -- an exactly-determined system,
not a regression. Leave-one-out cross-validation on the 4-point EXTENDED
set shows ~333% error predicting the held-out smallest device (improved
from ~680% after the Phase 31 Tusek point correction below, but still an
order-of-magnitude miss).
- Source: `core/loss_model.py` module docstring; `leave_one_out_cv()`,
  `run_extended_diagnostic()` (both directly re-runnable).
- Status: open. Needs 1-2 more independently-sourced, well-documented AMR
  devices added as genuine (not guessed) calibration points to become an
  actual regression with a reportable residual -- deliberately NOT
  fabricated in Phase 31 to close this gap artificially.

### 1.5 Hysteresis-share reversal finding is not stable across NSGA-III seeds
The Phase 16 finding that thermal hysteresis loss shifts material-family
representation on the Pareto front (La(Fe,Si)13Hy share 88%->100%) does
NOT reproduce consistently across random seeds at production NSGA-III
settings -- one seed shows the share moving the opposite direction.
- Source: `hysteresis_multiseed_stability.txt`,
  `results/hysteresis_sensitivity.txt`.
- Status: open. Needs either many more seeds to report a stable
  distribution, or should be reported as "direction-dependent on seed",
  not as a point-estimate finding.

---

## 2. Calibration-data provenance (now resolved in Phase 31, documented for the record)

### 2.1 Tusek calibration point corrected to the paper-verified field/point
Previously, `core/loss_model.py`'s `CALIBRATION_POINTS_CORE` used a
stopgap Tusek point (1.69T, 0.196kg, 0.25Hz, guessed Qc=6.5W) because the
paper-verified operating point (1.15T, 0.1763kg, 0.3Hz) did not calibrate
against an unverified guessed (span, Qc) pair. A rigorous pixel-calibrated
digitization of the source paper's actual Figs. 10-11 had already been
completed (`data/tusek_ate2013_figs/`) but was never wired into the
calibration.
- **Fixed in Phase 31**: `CALIBRATION_POINTS_CORE` now uses the verified
  point (span=7.26K, Qc=5.27W, COP=5.38 at the correct 1.15T/0.1763kg/
  0.3Hz), recomputed via the same brentq procedure used for every other
  CORE point. This calibrates cleanly (mdot_cal=0.007351 kg/s, matching
  the CSV's own independent cross-check to 2 sig figs).
- Every downstream number this touched was re-verified and updated:
  `core/loss_model.py`'s docstring narrative (leave-one-out errors,
  parasitic-fraction ranking), `tests/test_loss_model.py`'s hardcoded
  self-consistency dictionaries, and `results/regenerator_1d_validation.txt`.
- Full regression run: 130+ tests across every module depending on
  `StateDependentLossModel`'s default coefficients, all passing after the
  coefficient shift (k_eddy 30.52->31.25, base_frac 0.0484->0.0387).
- **Not yet regenerated**: several `results/*.txt` diagnostic files
  reference "Tusek" in contexts that may still show pre-fix numbers if
  they were generated by a one-off script rather than a function called
  fresh each pipeline run (`results/cycle_type_validation.txt`,
  `results/design_recommendations.txt`,
  `results/geometry_optimization_analysis.txt`,
  `results/regenerative_amplification_diagnostic.txt`). Re-run
  `python main.py` (or the specific analysis module) to refresh these
  before citing their exact numbers in a report.

---

## 3. Software robustness (found and fixed in Phase 31)

### 3.1 ProcessPoolExecutor hangs indefinitely in restricted environments
`core/cascade.py`'s parallel mass-sensitivity sweeps
(`run_magqueen_mass_sensitivity`, `run_cooltech_mass_sensitivity`) used
bare `future.result()` with no timeout, and several other pool call sites
used `pool.shutdown(wait=True)` with no timeout either. In a sandboxed
environment with unreliable process forking, this hung the calling
process indefinitely with no exception ever raised -- reproduced directly
via this project's own test suite
(`test_magqueen_mass_sensitivity_parallel_matches_sequential` and its
cooltech analog).
- **Fixed in Phase 31**: every `ProcessPoolExecutor` call site in
  `core/cascade.py` now goes through one of three new helpers
  (`_pool_map_or_none`, `_pool_submit_all_or_none`, `_safe_pool_shutdown`)
  that bound every wait to a hard timeout and fall back to the
  pre-existing sequential code path on any failure. A "poisoned executor"
  flag prevents a broken pool from being retried on every iteration of a
  `brentq` root-finding loop (an early version of this fix multiplied
  the wait instead of bounding it -- caught and corrected before landing).
- Verified: both previously-hanging tests now complete (87.5s and 44s
  respectively, bounded instead of unbounded); full downstream regression
  suite (95+ tests across every module importing `core.cascade`) passing.

### 3.2 Missing test coverage for 5 pipeline-wired modules
`core/regenerator_1d.py`, `uncertainty_propagation.py`, `water_usage.py`,
`pue_annualized.py`, and `commercial_landscape.py` were all imported and
called from `main.py` but had no dedicated `tests/test_*.py` file.
- **Fixed in Phase 31**: all 5 now have test files (65 new tests total,
  all passing). One genuine finding surfaced while writing these:
  `uncertainty_propagation.py`'s Monte Carlo `Qc` confidence band is
  architecturally near-zero-width (std ~1e-13) because cooling capacity
  does not depend on the calibrated loss coefficients in this codebase --
  only `COP_electrical` does. The module's own `Qc_p05`/`Qc_p95` columns
  are therefore not actually informative uncertainty bands; only the
  `COP_electrical_p05`/`p95` columns carry real calibration uncertainty.
  See `tests/test_uncertainty_propagation.py`'s
  `test_qc_is_essentially_invariant_to_calibration_uncertainty_unlike_cop`.

### 3.3 Material x n_layers cross-product was documented as a follow-up, never attempted
`run_layered_optimization()`'s own docstring explicitly scoped this out.
- **Fixed in Phase 31**: implemented as
  `run_layered_optimization_material_family_cross_product()` in
  `core/optimize.py`, reusing the existing per-family and per-n_layers
  Pareto-filtering machinery with no new NSGA-III formulation needed.
  Wired into `main.py` as step "11g.", **off by default** (pass
  `--layered-material-cross-product` to run it -- it multiplies step
  11f.'s already-reduced runtime by ~5x, one per material family).
  Surfaced and fixed a real latent bug in the process:
  `run_layered_optimization()`'s `out_csv` parameter crashed on
  `out_csv=None` (an inconsistency with `per_n_layers_out_dir`, which
  already handled `None` correctly) -- now guarded, with a regression
  test.

### 3.4 Several tests silently overwrite real `results/*.txt` and `results/*.csv` output files
Discovered while preparing this ledger: `tests/test_hysteresis_sensitivity.py` correctly
overrides `run_hysteresis_sensitivity()`'s `out_path` (the human-readable
report) with a `tmp_path`-based scratch file, but does **not** override
that function's `out_csv_on`/`out_csv_off` parameters -- which default to
the real `results/pareto_front_hysteresis_{on,off}.csv` files. Similarly,
`core/plots.py`'s figure-generation code calls
`optimize.run_optimization(out_csv=str(RESULTS_DIR / 'pareto_front.csv'))`
without overriding `per_material_out_dir`, which defaults to the real
`results/pareto_front_by_material/` directory. Running the test suite
(specifically `tests/test_hysteresis_sensitivity.py` and
`tests/test_plots.py`) therefore silently overwrites these real result
files with test-scale (small `pop_size`/`n_gen`) data as a side effect --
discovered directly in Phase 31 when a broad regression run corrupted
several of these files, which then had to be restored from the original
committed content before this phase's own changes could be packaged.
- Status: **open, not fixed in Phase 31** (out of scope for this pass --
  flagged rather than silently worked around). The correct fix is for
  every test that calls a function with a real default output path to
  explicitly override *every* such parameter, not just the most obviously
  named one; a stronger structural fix would be for these functions to
  default to `None` (skip writing) rather than a real repository path,
  matching the convention Phase 31 already applied to
  `run_layered_optimization()`'s `out_csv` in Section 3.3 above.
- Practical consequence: if you need trustworthy, current copies of
  `results/pareto_front.csv`, `results/pareto_front_by_material/*.csv`,
  `results/pareto_front_hysteresis_{on,off}.csv`, or
  `results/hysteresis_sensitivity.txt`, regenerate them with a fresh,
  deliberate call to `main.py` or the specific module function (not by
  running the test suite) after any change to `core/loss_model.py`'s
  calibration or NSGA-III search parameters.

---

## 4. Scope items not attempted in Phase 31 (still open)

- **Full bottom-up manufactured-system BOM.** Cost model still uses an
  order-of-magnitude multiplier for pumps/motors/controls rather than
  individually-priced components. See `core/economics.py` module
  docstring.
- **Reference-book gaps.** Kitanovski et al. (2015) pp.104-109
  (closed-form cycle-topology relations) and Tishin & Spichkin (2003)
  Ch.11 (passive regenerators) remain inaccessible even with the full
  `Papers/` corpus now present (see Section 5, items 5.1-5.2) -- several
  "qualitative ranking only" caveats in
  `core/passive_regenerator_analysis.py` and cycle-type sensitivity work
  trace back to this.
- **India/data-center-specific techno-economics.** `pue_annualized.py`
  and `water_usage.py` exist and are tested, but have not been run with
  actual Indian climate-zone data or commercial electricity tariffs.
- **Regeneration of stale `results/*.txt` diagnostic files** listed in
  Section 2.1 above.
- **Papers/ subfolders not yet mined**: `Reviews/`, `Data center
  cooling/`, and `AMR Theory and Modeling/` were not systematically
  cross-checked against the codebase's citations in Phase 35 (time-
  bounded scope decision, not a finding of any kind) -- a reasonable next
  pass for whoever picks this up next.

---

## 5. Phase 35: `Papers/` corpus verification pass

(Numbering note: this project's own numbering had independently reached
Phase 33/34 -- CORE_PLUS_TUSEK_MULTIPOINT calibration, NTU/utilization
diagnostics, beverage-cooler validation -- by the time this pass was
merged in, so it is labeled Phase 35 here rather than colliding with
those; the underlying work below was done directly against the original
`Papers/` corpus and does not depend on or conflict with Phases 33-34.)

The 41-paper corpus this project's citations had always referenced (many
comments say "PDF present in this repo" or similar) was, until this
phase, missing from the delivered project -- present in the original
developer's own working copy, but not included in earlier hand-offs of
this codebase. This phase added the corpus (`Papers/`, 137 MB, 8
subfolders) and used it to directly re-verify -- not just trust -- the
citations it backs. Headline result: the overwhelming majority of
citations checked out EXACTLY against primary sources (Dan'kov et al.
1998's 2T and 7.5T text-quoted values, Pecharsky & Gschneidner 1997's
Gd5Si2Ge2 peak entropy figure, Giguere et al. 1999's Gd cross-check
values, Bjork et al. 2011's $40/$20 per kg figures, Jacobs et al. 2014's
Astronautics 2502W/11.0K/COP=1.9/1.52kg figures) -- strong evidence the
project's prior "paper-mining pass" work was done with real (if
since-separated) primary-source access, not fabricated. Two genuine
errors were found and fixed:

- **`LITERATURE_DELTA_T_AD[5.0]` was wrong (14.6K, corrected to 12.3K)**.
  A pixel-calibrated re-digitization of Dan'kov et al.'s own Fig. 10 (not
  previously digitized in this project -- the old value was flagged as
  "plausible but not confirmed") gives 12.2-12.5K at 5T, cross-validated
  by reproducing the paper's own prose-stated ~15K figure at 7.5T to
  within 0.5K using the identical method. This flips the model's reported
  5T behavior from an apparent -7.5% underestimate to a +9.8%
  overestimate -- which is actually a MORE physically consistent result
  (mean-field theory should overestimate at every field near Tc, not flip
  sign at exactly one field), and substantially narrows what was
  previously described as a "genuine cross-paper discrepancy" against
  Giguere et al.'s independent Gd measurement. See `core/validation.py`'s
  updated `LITERATURE_DELTA_T_AD` comment for the full digitization
  methodology. This also flipped a knock-on finding in
  `core/inhomogeneous_broadening.py`: what was reported as "narrows the 1T
  error but widens the 5T error -- a genuine trade-off" is, under the
  corrected reference value, "narrows BOTH errors simultaneously" -- no
  code change was needed for this, the module's own conditional trade-off
  detection logic adapted automatically once the input data was corrected.
- **`core/mce_material.py`'s docstring said "-9.5 J/kg/K" for
  Gd5Si2Ge2's peak entropy change at 5T, corrected to "-18 J/kg/K"**.
  Direct inspection of Pecharsky & Gschneidner (1997) Fig. 4 confirms the
  peak is at ~18 J/(kg K), matching `core/giguere_validation.py`'s
  already-correct citation of the same paper/quantity -- this was a
  documentation-only inconsistency (the wrong value was never read by any
  function), not a computational bug.

Two genuinely new, additive citations were incorporated into
`core/economics.py`: Tura & Rowe (2014)'s $42/kg NdFeB and $10-20/kg bulk
Gd figures (triangulating Bjork et al. 2011's $40/$20 per kg), and Bjork,
Bahl & Nielsen (2016)'s specific $100 magnet + $40 MCM / 50W worked
example (a real literature anchor for `AMR_MAGNETIC`'s previously-
undocumented $2200/kW placeholder -- confirms it is reasonable and, if
anything, conservative, without resolving the separate HX/pump/motor/
controls gap, which the cited paper's own cost model shares).

### 5.1 Kitanovski et al. (2015) pp.104-109 remains inaccessible
Checked directly this phase: `Papers/Reference Books/Kitanovski et al.,
Magnetocaloric Energy Conversion (2015).pdf` is only a 30-page front-
matter/preface excerpt (table of contents and introduction, ending partway
into Chapter 1), not the full ~250+ page monograph -- it does not contain
Chapter 4 (Active Magnetic Regeneration) at all, so Sect. 4.1.1-4.1.4 (pp.
104-109), the source this project's `core/amr_cycle.py` cycle-type
factors were always meant to derive from, is still not available. Status
unchanged: qualitative Carnot >= Ericsson >= Brayton ranking only, not the
book's own closed-form relations.

### 5.2 Tishin & Spichkin (2003) remains an image-only, non-extractable PDF
Confirmed directly this phase: `Papers/Reference Books/Tishin & Spichkin,
The Magnetocaloric Effect and its Applications (2003).pdf` extracts zero
characters of real text from its sampled pages (scanned images only, no
OCR text layer) -- matching the exact finding already flagged in
`core/baseline_cooling.py` (Ch. 11, passive regenerators) and
`core/inhomogeneous_broadening.py` (Sect. 2.8, inhomogeneous
ferromagnets). Status unchanged: both remain open; OCR-ing the specific
needed chapters (not attempted this phase, out of scope for this pass)
would be the concrete next step.

Confirmed still-genuinely-open even with the full corpus available (not
just previously unchecked): (Mn,Fe)2(P,Si)/Ga1-xCMn3+x/Mn-Cu-Co-Ge have no
$/kg figure anywhere in the corpus (`Papers/Economics/The Resource Basis
of Magnetic Refrigeration.pdf` discusses these materials only in terms of
rare-earth supply criticality, not unit cost); Kitanovski et al. (2015)
pp.104-109 and all of Tishin & Spichkin (2003) remain inaccessible (see
Sections 5.1-5.2 above) even in this corpus.

---

## How to use this file

If you are preparing a presentation, report, or paper section based on
this codebase and need to state a limitation, check here first for the
authoritative, most-recently-verified wording -- then follow the "Source"
pointer for the full derivation/diagnosis before citing a specific number.
If you resolve one of the "Status: open" items above, update this file in
the same commit; do not let the ledger drift out of sync with the code
the way the pre-Phase-31 scattered version did.
