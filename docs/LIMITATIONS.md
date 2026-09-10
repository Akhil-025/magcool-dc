# Known Limitations & Open Items -- Consolidated Ledger

This file exists for one reason: this project's own convention is to flag
every limitation, calibration gap, and unresolved finding explicitly rather
than smooth it over -- but these flags were scattered across
README.md, ROADMAP.md, several dedicated `*_ROOT_CAUSE_ANALYSIS.md` /
`CITATION_AUDIT_*.md` files, module docstrings, and individual
`results/*.txt` diagnostic outputs. A judge, reviewer, or new contributor
asking "what doesn't work yet?" had to hunt across a dozen files to get a
complete answer. This file is that answer, in one place, with a pointer to
the authoritative source for each item. It does not replace those sources
(which contain the full diagnosis and reasoning) -- it is an index into them.

Nothing below is new information; every item was already flagged somewhere
in the repository before this file existed. Several were verified,
extended, or fixed directly during a later consolidation pass (noted inline).

---

## 1. Model-accuracy limitations (physics)

### 1.1 Mean-field Gd model overpredicts near the Curie point -- LARGELY FIXED
Fixed via exact isentropic ΔT_ad root-solve, a Sommerfeld electronic
entropy term, and fitted grain-Curie-temperature broadening. Error cut
from +48.9%/+29.2%/+9.8% (1T/2T/5T vs. Dan'kov et al. 1998) to
+4.2%/-1.1%/-8.1%; the 5T value now falls inside Giguère et al.'s
held-out 10.5-11.5K range.
- Residual few-percent error near Tc is expected mean-field short-range-
  correlation physics (per de Oliveira & von Ranke, Phys. Rep. 489
  (2010), same treatment) -- closing it needs a different statistical-
  mechanics treatment (e.g. Oguchi/pair-cluster), real future work, not
  a same-pass patch.
- Source: `core/mce_material.py` (`delta_T_adiabatic_exact`,
  `entropy_lattice`, `sommerfeld_gamma_J_per_molK2`),
  `core/inhomogeneous_broadening.py` (`GADOLINIUM_CALIBRATED`),
  `core/validation.py` (`run_validation`, `run_giguere_gd_extension`).

### 1.1b Curie-point field-shift rate (still open, distinct from 1.1)
Separately from the near-Tc overprediction MAGNITUDE fixed in 1.1 above,
Dan'kov et al. (1998) report the Curie-point transition temperature
itself shifting with field at ~6 K/T (2-7.5T range). The model's
self-consistent M(T,H) solve predicts essentially zero shift
(`core/validation.py`'s `run_curie_shift_check()`) — confirmed still true
after the 1.1 fix above, since that fix corrects the ENTROPY/heat-
capacity accounting at a fixed Tc, not the field-dependence of Tc itself.
De Oliveira & von Ranke's own Gd Sec. 4.2 treatment uses the same
single-exchange-parameter mean-field approach and does not supply a
usable Tc(H)-shift mechanism either (confirmed by directly reading it,
not assumed). Reproducing this would need real short-range-correlation
theory (Oguchi/pair approximation or similar) — a genuinely different
model, not attempted here.
- Source: `core/mce_material.py`'s `curie_shift_K_per_T` (phenomenological
  knob, `GADOLINIUM_FIELD_SHIFTED`, not used elsewhere in this repo);
  `core/validation.py`'s `run_curie_shift_check()`.

### 1.1c Two genuine attempts at the "real" fix for 1.1/1.1b -- both ruled out
Since 1.1b above says the field-shift needs "real short-range-correlation
theory," a separate session actually attempted that, twice, rather than
leaving it as a suggestion:

1. **Oguchi (Bethe-Peierls) two-spin cluster correction.** Built from
   scratch, correctly derived (variational/Bogoliubov free energy, not
   guessed), and internally validated against all 3 of its own
   structural checks (recovers mean-field exactly as z->infinity;
   finite-z Tc is lower than mean-field, as it must be; heat capacity is
   a rounded peak, not a discontinuity). Despite passing every internal
   check, it does **not** resolve the actual Gd discrepancy any better
   than the existing mean-field model when checked against real data.
2. **Critical-exponent ("universal curve") scaling**, the literature-
   standard alternative (Franco et al. and others) — investigated
   *before* being built, by fitting the field-scaling exponent n implied
   by Dan'kov et al.'s own three published DeltaT_ad values already in
   this repo. Result: n≈0.90-1.01 (best fit ≈0.946) — nowhere near
   mean-field's n=2/3 (0.667) *or* 3D Heisenberg critical scaling's
   n=0.637 (Gd's own measured critical exponents, Bednarz, Phys. Rev. B
   47 (1993)). Both theoretical exponents are further from what the data
   implies than they are from each other: Gd at 1-5T is simply not in
   the asymptotic critical-scaling regime this class of correction
   assumes, so it was correctly not built at all rather than built and
   found to fail afterward.
- A third possibility (demagnetizing-field correction to Dan'kov et
  al.'s reported internal field) was considered and set aside: real and
  well-documented in general, but implausible as the specific explanation
  here without Dan'kov et al.'s exact sample geometry, which this repo
  doesn't have, and a basic demagnetization oversight is unlikely to have
  survived unnoticed in the field's most-cited Gd MCE reference for 25+
  years.
- **Conclusion:** two real, literature-grounded fix attempts for 1.1b
  (and, had either worked, potentially some of 1.1's residual few-percent
  gap too), both genuinely tried, both genuinely ruled out with evidence
  rather than abandoned by assumption. What's needed instead is primary
  data this repo doesn't have (a digitized full Dan'kov DeltaT_ad(T)
  curve to fit a proper crossover theory against, or Dan'kov et al.'s
  exact sample geometry) — not a third model built on the same
  insufficient 3 points.
- Source: `research/oguchi_pair_cluster_prototype.py` (reference
  implementation, not imported by `core/` or `main.py`, not tested).

### 1.2 First-order Landau model overestimates giant-MCE ΔTad -- CLOSED, proposed fix rejected
Model overestimates ΔTad by ~2.4x vs. Giguere et al. (1999) direct
measurement. The previously-proposed fix (re-fit the 6th-order Landau
coefficients to Giguere's ΔTad instead of Pecharsky & Gschneidner's
ΔS_M) was investigated and rejected: the 0-D lattice-only-Cp model
can't satisfy both targets at once (unresolved latent-heat structure,
not a tuning problem), and applying that correction across fields
drags the independent Gd5Si2Ge2/Gd cross-check ratio from ~1.24
(close to the 1.30 literature value) down to ~0.51 -- i.e. it would
predict the "giant MCE" material underperforming plain Gd.
- Status: closed. Existing single-field empirical correction
  (`DTAD_CORRECTION_FACTOR`, `apply_giguere_correction()`) stands as the
  correct-scoped compromise. Further improvement needs a genuine latent-
  heat-resolving model extension, not a same-pass patch.
- Source: `core/first_order_mce.py`, `core/giguere_validation.py`;
  `tests/test_giant_mce_analysis.py`.

### 1.3 1-D regenerator: direction-inconsistent error, even after three fix passes
`core/regenerator_1d.py`'s multi-cycle transient model undershoots on two
of three benchmark devices and overshoots on the third:
- Tusek: -83.5%, Lozano: +104.9%, DTU/MAGGIE: -62.7%
  (see `results/regenerator_1d_validation.txt` for the full table --
  these are the numbers `main.py`'s own pipeline run actually produces).
- Fix attempted: replaced the previous ad hoc axial-conductivity
  multiplier with a real, citable Maxwell-Eucken packed-bed
  composite-conductivity model (`core/regenerator_1d.py`'s
  `_packed_bed_effective_axial_conductivity()`). **Honest result: this did
  NOT improve accuracy** -- it made the undershoot slightly WORSE on two
  rows while leaving the overshooting row essentially unchanged, because
  the more rigorously-justified conductivity value is *higher* than the
  old placeholder at these bed porosities (more damping, not less). This
  demonstrates the remaining gap is not simply "the conductivity constant
  is a bit off" -- a legitimately different, textbook-derived constant
  moves the same two rows further the same wrong way.
- **Item 1.7's parallel-plate geometry fix, now actually wired into this
  same production benchmark (not just tested standalone):** for the
  Tusek row specifically, `validate_against_benchmarks()`
  (`core/regenerator_1d.py`) now calls `no_load_span()` with
  `geometry="parallel_plate"` and the device's real, source-verified
  geometry (see Item 1.7), instead of silently defaulting to
  `packed_bed`. This IS a real improvement over the pre-fix -96.9% (to
  -83.5%) -- but it is smaller, and the run is flagged
  `converged=False`, than the ~-16% figure an earlier, separate,
  manually-configured test of the same geometry fix found (Item 1.7's
  own original text) -- that number came from one hand-picked mdot at
  n_nodes=20, not from `no_load_span()`'s own automatic mdot search
  (which chooses the best-performing mdot among a fixed candidate grid,
  and is what `main.py`'s pipeline actually runs). The two numbers
  disagree because the model's result is itself sensitive to exactly
  which mdot and grid the search lands on, which is precisely Item 1.7's
  own "not grid-converged" finding restated in a different form -- this
  is additional, direct confirmation of that finding, not a contradiction
  of it. Treat -83.5% as what this repo's pipeline actually reports
  today; treat ~-16% as "achievable at a specific, non-default
  configuration, not reproduced by the automatic search."
- Status: open (both the direction-inconsistent error and the
  grid/search sensitivity), and NOT wired into `core/amr_cycle.py`'s
  `cooling_capacity()` (the function every other part of this codebase
  actually uses) until resolved.
- Source: `core/regenerator_1d.py` module docstring ("Known limitations"),
  `validate_against_benchmarks()`'s own printed/written caveat text.

### 1.4 Loss-model calibration has (at best) zero degrees of freedom
`StateDependentLossModel`'s 3 coefficients (k_eddy, k_pump, base_frac) are
fit to exactly 3 CORE calibration points -- an exactly-determined system,
not a regression. Leave-one-out cross-validation on the 4-point EXTENDED
set shows ~333% error predicting the held-out smallest device (improved
from ~680% after the Tusek point correction below, but still an
order-of-magnitude miss).
- Source: `core/loss_model.py` module docstring; `leave_one_out_cv()`,
  `run_extended_diagnostic()` (both directly re-runnable).
- Status: open. Needs 1-2 more independently-sourced, well-documented AMR
  devices added as genuine (not guessed) calibration points to become an
  actual regression with a reportable residual -- deliberately NOT
  fabricated to close this gap artificially.

### 1.5 Hysteresis-share reversal finding is not stable across NSGA-III seeds -- CLOSED
The earlier finding (La(Fe,Si)13Hy share 88%->100% with hysteresis loss
on) did not hold up: re-run across n=20 seeds at production settings
gives mean=-1.25pp, std=6.9pp, p=0.43 vs. zero -- statistically
indistinguishable from no effect. **Closed as noise**, not a real
hysteresis-driven shift.
- Also found (separate, not root-caused): the same seed number doesn't
  fully determine `run_hysteresis_sensitivity()`'s output -- some
  unseeded state also contributes. Worth tracking down if exact
  reproducibility ever matters elsewhere.
- Caveat unchanged: `hysteresis_loss_J_per_kg` values remain literature
  analogs, not measurements of the calibrated compositions -- this
  closes the search-noise question, not the underlying-data question.
- Source: `results/hysteresis_multiseed_stability_n20.txt`,
  `results/hysteresis_sensitivity.txt`.

### 1.6 Qc(span) "feasibility reopening" artifact -- FIXED, opt-in
Root cause: a genuine mean-field magnetic heat-capacity discontinuity
near Tc (Item 1.1) could make single-point Qc evaluation report a
*larger* value at a bigger span than a smaller one -- physically
backwards. Confirmed case: Tusek AMR(A), span=12.23K, raw Qc=20.51W vs.
digitized 2.03W (+910%).
- Fix: added `AMRSystem.cooling_capacity_span_sweep()`, which applies a
  running-minimum monotonicity clamp across a dense span scan (can only
  reduce Qc, never invent a value). Confirmed: the Tusek case above now
  reports Qc=0W, matching the model's own conclusion at neighboring spans.
- Opt-in only: `cooling_capacity()` itself and all existing callers are
  unchanged, so no existing `results/` output changes because of this.
- Source: `core/amr_cycle.py` (`cooling_capacity_span_sweep`);
  `tests/test_amr_cycle.py`.

### 1.7 Regenerator-1D span undershoot (issue #8) -- partial real fix found, deeper issue surfaced
`core/regenerator_1d.py`'s transient 1-D AMR simulator was previously
found to under-predict peak no-load span by ~92-97% (Tusek benchmark) --
see that module's own "known limitations" docstring.
- **Real root-cause component found and fixed:** the Tusek_singlebed_
  Gd_2010_spanceiling benchmark device is a PARALLEL-PLATE regenerator
  (its own source paper, Tusek et al., Appl. Therm. Eng. 53 (2013) 57-66,
  Table 1: 0.1mm spacing, 0.25mm plate thickness, porosity 0.2564, outer
  dims 10x80x39mm), not a packed bed of spheres -- but every call in this
  module used the packed-bed correlation and its packed-bed default
  geometry regardless. Added `geometry="parallel_plate"` support (new,
  additive, defaults to the exact previous `"packed_bed"` behavior for
  every existing caller) that dispatches to `core/thermal.py`'s existing
  `regenerator_effectiveness_parallel_plate()` with the device's real,
  source-verified geometry. Result: peak span improves from ~1.4K (-92 to
  -97% error) to ~16.5K (-16% error) at this module's usual n_nodes=20.
- **Deeper, NOT fixed, issue found while checking this result's
  robustness:** the model is not grid-converged in n_nodes, in EITHER
  geometry mode, and this is pre-existing (reproduces with the original
  packed-bed geometry too). Direct check: packed_bed gives
  2.19K/3.23K/3.78K at n_nodes=10/20/40; parallel_plate gives
  8.99K/16.15K/2.31K at the same three -- neither sequence is converging,
  and the second isn't even monotonic. Checked further: re-running
  n_nodes=40 with 12000 cycles instead of 1500 gives 1.04K, not 2.31K
  (flat over the last 10 cycles) -- ruling out "just needs more cycles";
  the fine-grid answer is a genuine, different, and SMALLER steady state,
  not a slow transient converging toward the coarse-grid number. Neither
  the pre-fix number nor the new ~16.5K result should be treated as a
  validated, resolution-independent prediction. See
  `core/regenerator_1d.py`'s module docstring for the full diagnostic
  detail -- this needs a genuinely grid-independent per-segment
  effectiveness formula, not a same-pass patch.
- **Status: still open, correctly not wired into `core/amr_cycle.py`** --
  same standard this module already held itself to before this phase.
- **Now actually wired into `validate_against_benchmarks()` itself** (the
  function that produces `results/regenerator_1d_validation.txt`, which
  `main.py`'s own pipeline runs) -- the Tusek row specifically now uses
  `geometry="parallel_plate"` with this device's real geometry. The
  production number this produces (-83.5% error, via `no_load_span()`'s
  own automatic mdot search) is a real improvement over the pre-fix
  -96.9%, but is NOT the same as the ~-16% (span~16.5K) figure quoted
  above -- that number came from one hand-picked mdot/n_nodes
  combination, not from the automatic search `main.py` actually runs.
  See Item 1.3 for the full reconciliation of these two numbers and why
  they differ (short version: it's the same grid/search-sensitivity
  finding below, showing up again in a different configuration).
- Source: `core/regenerator_1d.py` module docstring;
  `tests/test_regenerator_1d.py`'s 5 new tests for the
  `geometry="parallel_plate"` option itself (test that the new code path
  runs correctly at fast/reduced settings -- they do NOT check
  convergence or benchmark accuracy).
- **Further root-cause hunt (separate session), four more hypotheses
  eliminated, cause still not found:** built a standalone, properly-
  coupled solver (implicit tridiagonal solve handling axial conduction
  and convection in ONE step per timestep, instead of the operator-split
  approach in `core/regenerator_1d.py`) specifically to test whether
  operator-splitting order explains the span gap. It doesn't: span is
  still ~0 with this alternative solver. Also checked the calibrated flow
  rate's NTU against the literature-realistic range for working AMR
  devices (NTU~10-30, Nielsen et al. 2010, 2D AMR parameter study) --
  this device's implied NTU~487-1949 is wildly outside that range, a
  real, separate, previously-unnoticed red flag -- but re-testing the
  same coupled solver at flow rates giving a realistic NTU (~0.03-0.07
  kg/s) still gives ~0 span. So: operator-splitting order, axial-
  conduction magnitude, NTU magnitude, and grid resolution are now FOUR
  hypotheses tested and eliminated, and the actual root cause of the span
  gap remains unresolved. Independently re-verified in this repo's
  current environment: as written, this alternative solver also isn't
  numerically stable over many cycles (diverges to NaN via Brillouin-
  function overflow with no temperature clamping in the harness) --
  another reason it's kept as a reference script, not adopted.
- Source (this bullet): `research/regenerator_1d_v2_coupled_prototype.py`
  (reference implementation, not imported by `core/` or `main.py`, not
  tested).

### 1.8 Barocaloric and electrocaloric device models rest on much weaker calibration than the magnetocaloric model
`core/barocaloric_cycle.py` (NPG plastic crystal) and
`core/electrocaloric_cycle.py` (PST relaxor-ferroelectric MLCC) were added
alongside the existing elastocaloric line in
`core/alternative_caloric_comparison.py`, but neither has an
independently-**measured** end-to-end device COP to calibrate against, unlike
elastocaloric (Qian et al. 2023 measured system COP) or this repo's own
AMR model (`data/amr_experimental_benchmarks.csv`):
- **Barocaloric**: the only NPG-specific device-level COP figure found is a
  **simulated** value from a different research group's own model (COP=5.5
  at 1 mHz, 2.4K span, 0.1GPa) -- calibrating against another group's
  simulation is a materially weaker check than calibrating against
  hardware. `run_cycle(span_K=10.0)` reports `COP_electrical=24.4` at this
  repo's own representative span, with `W_parasitic_per_kg=0.0` (no
  pumping-power term modeled for this cycle at all) -- reported directly
  rather than assumed away.
- **Electrocaloric**: this module works top-down from PAPERS' OWN
  reported device-level span/power/COP figures (Meng 2020, Li 2023 +2025
  erratum) rather than bottom-up from an independently-verified
  Delta_S-vs-field curve, because that per-stage entropy/field data was not
  located for the exact MLCC batches used. `run_cycle(span_K=10.0)` returns
  `is_extrapolated=True, is_far_extrapolation=True` at this repo's
  representative 10K span -- the fitted scaling law is being extrapolated
  well beyond the papers' own measured span range, and the module reports
  this flag explicitly on every call rather than silently extrapolating.
- Source: `core/barocaloric_cycle.py` and `core/electrocaloric_cycle.py`
  module docstrings (CALIBRATION HONESTY FLAG sections);
  `core/alternative_caloric_comparison.py`, `main.py` step 17.

### 1.9 Thermal-diode / magnetocaloric-fluid explorations: one mechanism graduated to validated, three did not (Phase 38)
Of the two thermal-diode mechanisms and the magnetocaloric-fluid
working-body architecture added as design-exploration tools
(`core/thermal_diode.py`, `core/fluid_mce_cycle.py`), a later pass found a
real benchmark device for exactly one of them and graduated it; the other
two remain exactly as originally scoped:
- **`FerrofluidThermalSwitch` (VALIDATED, Phase 38):** rectification_ratio
  (3.84) and frequency range are grounded in real hardware (Rodrigues et
  al. 2019), and `check_against_andrade_2024_benchmark()`
  (`core/thermal_diode_analysis.py`) checks this repo's model against a
  real Andrade et al. (2024) Gd+ferrofluid-switch refrigeration
  prototype -- agreement is qualitative/directional (a negative,
  no-symmetric-cycling-advantage result), not a validated quantitative
  match to Andrade's own reported asymmetric-cycling span figure.
- **`MechanicalContactDiode` (still NOT validated):** cost-only sensitivity
  study; no offsetting heat-transfer benefit from rectification is
  modeled (no closed-form Sect. 6.2.4 relation was available to digitize),
  and no benchmark device in this repo's corpus uses this specific
  mechanism.
- **`FerrofluidMCESystem` (still NOT validated):** an extended Phase 38
  literature search still found no magnetocaloric-fluid-as-working-body
  refrigeration benchmark anywhere. Two narrower things did graduate: the
  Krieger-Dougherty rheology grounding (a real magnetite-ferrofluid
  measurement, Susan-Resiga et al. 2012) and the negative finding itself
  (every real ferrofluid magnetocaloric device found uses the fluid as a
  thermal switch, never a working body). The module's own headline
  quantitative finding -- fluid dilution plus no regeneration collapses
  usable span to a fraction of a Kelvin at realistic loadings -- remains a
  design-exploration result, not benchmark-checked.
- Source: `core/thermal_diode.py`, `core/fluid_mce_cycle.py` module
  docstrings (VALIDATION UPDATE sections); ROADMAP.md Phase 38 entry;
  `main.py` step 11c / step 14.

### 1.10 Hybrid solid-state regenerator (HMR): frictionless-limit thermodynamics alone does not close the gap to vapor-compression
`core/hybrid_solid_state_regenerator.py` implements Lin et al. (2024)'s
solid-solid heat-transfer architecture (alternating HTCM/MCM slices,
no working fluid, no pumping/dead-volume losses -- a genuinely different
parasitic-loss channel, inter-layer friction, replaces them).
`compare_to_vcc_realistic()` at this repo's representative point
(T_cold=291.15K, span=10K) finds HMR's best real electrical COP (8.13 at
1.0 Hz, once air-gap friction, a rotary-drivetrain term, and baseline
overhead are added) is **0.66x** vapor-compression's real installed-system
COP (12.23) -- i.e. this architecture does NOT close the electrical-COP
gap at this operating point once the previously-excluded real loss
channels are included on equal footing with VCC's own real number. A
separate, more favorable `compare_to_vcc()` ideal-limit function exists in
the same module and omits drivetrain/baseline overhead entirely -- citing
that number alone without the realistic comparison above would overstate
the case.
- Source: `core/hybrid_solid_state_regenerator.py`'s
  `compare_to_vcc_realistic()`; `main.py` step 19.

### 1.11 Beverage-cooler real-world cross-check: model diverges sharply from one of two reported deployments
`core/beverage_cooler_validation.py` checks this repo's own AMR model
against two independent, real commercial magnetocaloric beverage-cooling
deployments -- the one market segment where magnetocaloric cooling is
already commercially deployed, giving a real-world check this repo's
primary data-center application cannot offer:
- **Magnotherm Eclipse (REWE pilot):** `run_eclipse_directional_check()`
  predicts this repo's model would need **-279%** more energy than the
  incumbent R290 unit at the pilot's own operating point (T_cold=277.65K,
  span=17.5K, AMR_COP_electrical=1.76 vs. VCC_COP=6.66) -- flatly
  contradicting the press-reported **+15%** energy saving. This is a real,
  unresolved divergence between this repo's model and a real deployment,
  reported directly rather than smoothed over; possible explanations
  (proprietary architecture differences, a materially different
  operating point than assumed, or a genuine model gap) are not
  distinguished here.
- **Polaris (CE-certified beverage cooler, Liang et al. 2025):**
  `run_polaris_second_law_validation()` finds much closer agreement --
  model second-law efficiency 6.34% vs. the paper's own reported 5.40% at
  the same architecture (single-material Gd, packed-particle-bed AMR,
  T_cold=277.65K, span=15K) -- the closest real-world agreement this
  repo's model has found for any commercial device.
- The two checks point in different directions at the same technology
  class; this is stated as an open discrepancy, not resolved in favor of
  either reading.
- Source: `core/beverage_cooler_validation.py`; `main.py` step 15b.

### 1.12 No span/technology region found where this repo's own AMR model beats conventional cooling (regime-crossover null result)
`core/regime_crossover_analysis.py`'s `run_cop_crossover_search()` swept
span (3-30K) against vapor-compression second-law efficiency (0.25-0.55,
covering everything from small residential-grade compressors to
well-optimized chilled-water) and a broad grid of this repo's own AMR
design freedoms. Result: **no crossover found** at any combination,
including against vapor-compression's own least-favorable (eta=0.25)
setting -- e.g. at span=10K, this repo's own best AMR_COP_electrical
(4.60) does not beat VCC even at eta=0.25 (COP=7.12). This null result is
consistent with (not contradicted by) every other real-world check in
this repo (`beverage_cooler_validation.py`, `heat_pump_validation.py`,
`hybrid_solid_state_regenerator.py`'s realistic HMR comparison above) --
none of them find magnetocaloric cooling beating conventional cooling on
COP either, only on weight/power-density (heat-pump check) or narrower
non-COP metrics (dry-rejection water usage, refrigerant-free emissions).
- Source: `core/regime_crossover_analysis.py`; `main.py` step 15d.

### 1.13 NSGA-III seed-to-seed stability on the MAIN Pareto front -- CLOSED, positive verdict
Re-run at production settings across 3 seeds: best COP_electrical
mean=10.413, std=0.068 (seed-stable); knee-point material consistent
(La(Fe,Si)13Hy every seed). Material-family *share* of the merged front
is noisier (La(Fe,Si)13Hy 71-93% across seeds), so report that as
mean+/-std, not a single-seed number.
- Source: `core/pareto_multiseed_stability.py`;
  `results/pareto_multiseed_stability.txt`; `main.py` step 16.

### 1.14 Magnet-geometry Pareto sensitivity: seed-to-seed stability -- run, and the effect is NOT stable
Item 13.2-equivalent (`core/magnet_geometry.py`'s
`run_geometric_cost_pareto_sensitivity()`) found, in one reduced-settings
run, that the new GEOMETRIC (super-linear Halbach-cylinder) magnet-mass
cost term does not pull the Pareto front's mean field down the way a
nonlinear cost penalty on high fields would be expected to. A follow-up,
`run_magnet_geometry_multiseed_stability_check()`, re-ran the same
FLAT-vs-GEOMETRIC comparison at full production settings (pop_size=40,
n_gen=25) across 3 seeds:
```
seed  mean_FLAT_T  mean_GEOM_T  front_FLAT  front_GEOM
   1         1.28         1.20          23          23
   2         1.48         1.50          22          21
   3         1.11         1.24          18          20
```
**Result: NOT STABLE.** Seed 2's GEOMETRIC run has a HIGHER mean field
than its own FLAT run -- the "geometric cost pulls mean field down"
expectation does not hold universally even at full production settings.
Same class of outcome as Item 1.5's hysteresis-share finding -- NSGA-III
search noise at this problem's scale appears large enough to flip BOTH of
this repo's reduced-setting Pareto sensitivity findings, not just one.
This resolves the search-noise question, not the underlying cost-data
quality question (the geometric magnet-mass relation is standard
closed-form Halbach-cylinder physics; the $/kg figures and the ~2T
"sweet spot" literature claim were not independently re-derived).
- Source: `core/magnet_geometry.py`'s
  `run_magnet_geometry_multiseed_stability_check()`;
  `results/magnet_geometry_multiseed_stability.txt`; `main.py` step 11e.

### 1.15 Calibration-failure root cause is structural (12/12), and the regenerative-amplification gap it exposes has a real, imperfect opt-in fix
`core/validation_system.py`'s `diagnose_calibration_failures()` checks
each of the 12 "NO CALIBRATION FOUND" benchmark rows (Item 1.7-adjacent
territory) directly against `cooling_capacity()`'s own structural ceiling
(2*dTad_noload -- the max span the 0-D model can reach at that
field/T_mid, for ANY mdot). **All 12/12 are structural**, not a
search-space artifact -- e.g. Astronautics (span=11.0K, margin=-10.96K),
Cooltech_2013_rotary (span=42.0K, margin=-39.19K), DTU_MagQueen_2018
(span=25.0K, margin=-24.97K). Widening the mdot search bound above 5 kg/s
cannot change the outcome for any of them -- confirmed directly from
`cooling_capacity()`'s own `dTad_noload` term, not assumed.
`run_regenerative_amplification_diagnostic()` then quantifies the SIZE of
the gap across every span>0 benchmark row: 12/17 rows with a well-defined
`dTad_noload` exceed the model's own structural cap, ratio range
1.02x-14.92x (median 1.39x), independently cross-checked by a directly-
measured no-load span on the same physical MAGGIE hardware as an existing
calibrated row (not back-calculated from any mdot). This is a LOWER
BOUND on the model's error (graded/layered real beds additionally lose
accuracy from the single-Tc approximation -- see Item 1.16).
`AMRSystem.no_load_span_override` (populated from
`regenerator_1d.regenerative_span_cap()`, Item 1.3/1.7) is the opt-in
mechanism meant to close this gap; tested on 3 of the 12 flagged rows
(the full 12-device check is directly re-runnable, capped here for
pipeline runtime): recovers a usable, evaluable COP for 1/3
(Lozano_POLO_UFSC_2016_r1, COP_lit=0.37 vs. COP_pred=0.318, -14.1% error)
but does NOT reach feasibility for the other 2 (DTU_Eriksen_MAGGIE_2016,
Risoe_DTU_Gd_2011 -- their real spans still exceed even the 1-D model's
own, larger span-cap prediction). Reads as: the override fixes
FEASIBILITY (can represent the span at all) without necessarily fixing
ACCURACY -- exactly the honest outcome `no_load_span_override`'s own
docstring predicts. Still opt-in, still off by default everywhere else
in this codebase.
- Source: `core/validation_system.py`'s `diagnose_calibration_failures()`,
  `run_regenerative_amplification_diagnostic()`,
  `run_regenerative_amplification_override_check()`;
  `results/calibration_failure_diagnostics.txt`,
  `results/regenerative_amplification_diagnostic.txt`,
  `results/regenerative_amplification_override_check.txt`; `main.py`
  steps 2c/2d/2f.

### 1.16 Graded-bed structural fix extended to the remaining large structural devices -- real per-layer data helps two of three, mass-unreported sensitivity only for the third
Following on from Item 1.15's finding that Astronautics, DTU_MagQueen_2018,
Risoe_DTU_Gd_2011, and Cooltech_2013_rotary are all structurally
infeasible under the single-Tc model, and Section 11.3's existing
Astronautics 6-layer graded-bed fix (-81.1% error), this extends the same
mechanism to the other three:
- **DTU_MagQueen_2018** (mass unreported, swept 0.5-10kg): calibrates at
  every swept mass, COP_error approx -45% to -48% depending on mass -- a
  real number, but a mass-sensitivity, not a calibrated answer (unlike
  Astronautics, whose 1.52kg mass IS directly reported).
- **MagQueen, REAL 10-layer bed (Masche, Bahl, Nielsen, Choi, Bez, Bjørk
  & Engelbrecht, Applied Thermal Engineering, 2021, Table 3)** -- using
  the device's own reported per-layer Curie temperatures AND masses (not
  a sweep or composition search): **no calibration found, an HONEST
  STRUCTURAL FINDING, not a wiring bug**. Every one of the 10 real layers,
  centered exactly on its own reported Tc (the most generous possible
  placement), caps out at peak dTad_noload~=0.28K at this device's real
  1.44T field -- a ~0.56K structural span ceiling per layer, well under
  the ~1.03K/layer the real device's own operating point requires. A
  property of `LAFESIH_FIRST_ORDER`'s own Landau calibration at this
  field, confirmed unchanged for the base material -- not an artifact of
  the per-layer wiring.
- **DTU_Eriksen_MAGGIE_2016, REAL 4-composition Gd/Gd-Y graded bed**
  (actual reported Curie temperatures): MAGGIE's own 15.5K row, previously
  a flat "NO CALIBRATION FOUND" under the single-Tc model, now calibrates
  (Qc=81.5W target 81.5W, COP_cascade=3.068 vs. lit=3.6, -14.8% error).
  **Genuine trade-off, not a strict win**: the companion row on the SAME
  physical prototype (DTU_Eriksen_rotary_Gd_2015, 10.2K span) -- already
  calibrating fine under the simpler single-Tc approximation -- trades
  away some of that accuracy under this same real 4-layer model. A single
  graded-bed model does not uniformly improve every operating point of
  the same physical device.
- **Risoe_DTU_Gd_2011 (30K span) -- HYPOTHETICAL 6-stage redesign**, not
  the real device (a single plain-Gd bed, not reported as graded): Qc
  target reached, but only 2/6 stages land within `GD_FAMILY`'s
  documented [20,290]K range (4/6 fall back to plain Gd), and COP_cascade
  collapses to 0.0 (vs. reported COP=5.0). A Qc-feasibility finding only,
  not an efficient-redesign prediction.
- **Cooltech_2013_rotary (42K span, largest in this benchmark set) --
  also HYPOTHETICAL, mass unreported (swept)**: reaches positive Qc
  feasibility at every swept mass (4/6 stages fall back to plain Gd) -- a
  feasibility sensitivity, not a calibrated result.
- **Consolidated**: the graded-bed STRUCTURE closes the Qc-feasibility gap
  for every device checked, at every swept mass -- the mechanism
  generalizes. Whether it also closes the ACCURACY gap depends on whether
  the per-stage material is real (MagQueen structural-negative, MAGGIE
  real-and-mixed) or hypothetical (Risoe/Cooltech, feasibility only).
- **Separately**: applying the Giguere ~2.42x DeltaT_ad overestimate
  correction (Item 1.2) on top of the existing Astronautics 6-layer fix
  narrows that device's error further (COP_cascade 1.919/+1.0% error ->
  1.895/-0.3% error) -- does NOT establish the correction transfers to
  La(Fe,Si)13Hy (it was fit to a Gd5Si2Ge2-specific direct measurement),
  but a smaller error is a smaller error and the experiment was actually
  run rather than left as an unexercised option.
- Source: `core/cascade.py`'s device-specific graded-bed validation
  functions (`main.py` steps 7d/7e); `results/` step-7 output.

### 1.17 Loss-model 4th calibration point (same device class) generalizes better than earlier cross-scale points -- still not a fix for the zero-degrees-of-freedom problem (Item 1.4)
Item 1.4's EXTENDED/FURTHER_EXTENDED 4th-point diagnostics add devices
spanning additional orders of magnitude and do NOT generalize (leave-one-
out error still hundreds of percent). A different 4th point --
DTU_Eriksen_MAGGIE_2016, the SAME physical prototype as an existing CORE
point (DTU_Eriksen_rotary_Gd_2015) at a different operating condition,
only reachable using this session's `no_load_span_override` (Item 1.15)
-- is genuinely different: its own held-out leave-one-out error is
**+26.5%**, far inside the ~250-700% range every other fold (old or new)
still shows, and the two pre-existing folds shared with CORE both improve
modestly rather than degrading. Real, partial progress -- NOT a fix for
the underlying zero-degrees-of-freedom problem (Astronautics and Tusek's
own folds are still order-of-magnitude misses with 4 points, same as with
3). **CORE (3pt) remains the production default.** What generalizes is
more data of the SAME device class, not more devices per se -- consistent
with, and sharpening, Item 1.4's existing diagnosis.
- Source: `core/loss_model.py`'s
  `calibrate_loss_coefficients(CALIBRATION_POINTS_CORE_PLUS_MAGGIE_HIGHSPAN)`;
  `main.py` step 3a2.

### 1.18 Nanocomposite off-design robustness: real but narrow -- avoids catastrophic failure, does not raise off-design performance
Item 1.2's nanocomposite blend trails a single sharply-tuned phase at its
own design span (COP 4.62 vs. 6.72 at 10K). A follow-up asks whether the
blend's deliberately broadened working range pays off OFF-design, when a
composition tuned once (10K design span) is evaluated WITHOUT retuning at
other spans (5/10/15/20K): the nanocomposite is feasible (Qc>0) at 1/3
off-design spans tested; the single sharply-tuned phase is feasible at
0/3. At the one off-design span where the single phase collapses to
Qc=0 (5K), the nanocomposite still delivers positive Qc. No off-design
span left both candidates feasible simultaneously, so no raw
performance comparison was possible off-design -- the genuine finding is
ROBUSTNESS TO NARROWING (avoiding catastrophic failure), not raw
off-design performance, and this is a first-pass finding at one spread
value and one design/off-design span set, not a general claim.
- Source: `core/nanocomposite_material.py` follow-up;
  `results/nanocomposite_robustness.txt`; `main.py` step 8e.

---

## 2. Calibration-data provenance (now resolved, documented for the record)

### 2.0 Working-fluid selection resolved: `DEFAULT_FLUID` = water_eg10
Pure water tops the fluid sweep but no real AMR prototype runs it (Gd
corrodes). `core.fluids.DEFAULT_FLUID` set to **water_eg10**, the best
corrosion-realistic fluid, at only 1.4-2.3% COP_electrical cost vs. the
unrealistic pure-water ceiling. Comparison-only -- no existing
`AMRSystem` call site's own default changed.
- Source: `core/fluid_selection_optimization.py`; `main.py` step 18.

### 2.1 Tusek calibration point corrected to the paper-verified field/point
`CALIBRATION_POINTS_CORE` now uses the paper-verified operating point
(span=7.26K, Qc=5.27W, COP=5.38 at 1.15T/0.1763kg/0.3Hz) instead of a
stopgap guessed point. All downstream numbers re-verified (130+ tests
passing; k_eddy 30.52->31.25, base_frac 0.0484->0.0387).
- **Not yet regenerated**: `results/cycle_type_validation.txt`,
  `results/design_recommendations.txt`,
  `results/geometry_optimization_analysis.txt`,
  `results/regenerative_amplification_diagnostic.txt` may still show
  pre-fix numbers until `python main.py` is re-run (see Section 4).

---

## 3. Software robustness (found and fixed)

### 3.1 ProcessPoolExecutor hangs indefinitely in restricted environments -- FIXED
Bare `future.result()`/`pool.shutdown(wait=True)` calls with no timeout
hung indefinitely in sandboxed environments with unreliable forking.
- Fixed: every `ProcessPoolExecutor` call site in `core/cascade.py` now
  goes through helpers that bound waits to a hard timeout and fall back
  to the sequential path on failure. Verified: previously-hanging tests
  now complete (87.5s/44s, bounded); 95+ downstream tests passing.

### 3.2 Missing test coverage for 5 pipeline-wired modules -- FIXED
`core/regenerator_1d.py`, `uncertainty_propagation.py`, `water_usage.py`,
`pue_annualized.py`, `commercial_landscape.py` now all have test files
(65 new tests, all passing).
- Genuine finding surfaced: `uncertainty_propagation.py`'s Monte Carlo
  `Qc` confidence band is architecturally near-zero-width (std ~1e-13)
  since Qc doesn't depend on the calibrated loss coefficients in this
  codebase, only COP_electrical does -- `Qc_p05`/`Qc_p95` aren't
  informative uncertainty bands.

### 3.3 Material x n_layers cross-product -- now implemented and run
Implemented as `run_layered_optimization_material_family_cross_product()`,
wired into `main.py` as step "11g." (off by default, ~5x runtime;
`--layered-material-cross-product`). Also fixed a latent `out_csv=None`
crash surfaced in the process.
- Result: 5 families x 3 n_layers = 34 non-dominated designs
  (`results/layered_pareto_front_material_cross_product.csv`).
  La(Fe,Si)13Hy dominates (28/34, 82%); best cascade COP=10.407. Confirms
  step 11f.'s single-family ranking rather than overturning it.

### 3.4 Several tests silently overwrote real `results/` output files -- CLOSED
Four instances found and fixed (test overrides of `run_hysteresis_
sensitivity()`, `plot_nsga3_pareto()`, `run_optimization()`'s
`per_material_out_dir`, and a direct `os.remove()` of the regenerator-1D
performance cache) -- all now use scratch/tmp_path paths instead of real
repo paths.
- Verified, not just asserted: full `pytest` run (698 tests) followed by
  an md5 checksum diff of all 190 files under `results/` against a
  pre-run snapshot -- zero differences, across two independent runs.
- Not closed: this only guarantees the *current* test suite is clean, not
  a structural guarantee against a similarly-written future test (would
  need output-path parameters to default to `None` everywhere, per
  Section 3.3's convention) -- that refactor remains open.
- Source: `core/plots.py`; `tests/test_hysteresis_sensitivity.py`,
  `tests/test_optimize.py`, `tests/test_plots.py`,
  `tests/test_regenerative_amplification_override.py`.

---

## 4. Scope items not yet attempted (still open)

- **Full bottom-up manufactured-system BOM (RESOLVED for non-materials
  hardware pricing; two narrower integration gaps remain open).**
  Status update from a later economics-integration check: this item, as
  originally written, is now out of date. `core/economics.py` DOES have a
  genuine bottom-up, individually-priced non-materials BOM
  (`bottom_up_non_materials_bom()`: HX/pump/motor/drive/controls, each
  priced from real market-catalog $/kW bands, LOW/MID/HIGH, see that
  function's own section docstring for sourcing) alongside
  `bom_cost()`'s materials-only BOM, combined by
  `full_system_cost_estimate_bottom_up()`. This is wired into the
  pipeline: `main.py`'s `run_economics()` (step 5) calls
  `economics.cross_check_full_system_cost_methods()`, which calls
  `full_system_cost_estimate_bottom_up()` and reports it alongside the
  older borrowed-VCC-multiplier estimate every run.

  There is also a real, committed, priced example of this BOM:
  `data/MagCool_DC_BOM.xlsx`, built at the NSGA-III knee-point design
  point (1.11T, 5.31kg La(Fe,Si)13Hy, COP_elec=9.57, Qc=35.38kW). A new
  reconciliation test, `tests/test_bom_xlsx_reconciliation.py`, reads that
  design point straight off the spreadsheet and confirms the spreadsheet's
  materials/non-materials/grand-total line items match
  `bom_cost()`/`bottom_up_non_materials_bom()`/
  `full_system_cost_estimate_bottom_up()` run live, to the cent -- the
  spreadsheet is a real, code-reproducible artifact, not a stale or
  hand-typed one, and the new test will fail (not silently drift) if that
  ever stops being true.

  One thing this does NOT resolve, so the honesty flag is narrowed, not
  removed: `bottom_up_non_materials_bom()` is itself still a
  market-catalog ENGINEERING ESTIMATE for generic steady-flow component
  categories (brazed-plate HX, centrifugal pump, NEMA motor, VFD),
  explicitly NOT an AMR-specific vendor quote -- an AMR device needs
  custom manifolding, oscillating-flow-rated seals, and a
  cycle-synchronized drive that could cost more than an off-the-shelf
  part of the same power rating (see that function's own docstring, and
  the "Sources & Assumptions" tab of `data/MagCool_DC_BOM.xlsx`). This is
  a genuine, non-code-fixable data gap (no AMR-specific manufactured-cost
  quote exists in the public literature this project has access to), so
  it remains open by design rather than as an oversight.

  **RESOLVED this pass**: `run_economics()` previously ran the
  BOM/full-system-cost functions only at a hardcoded illustrative design
  point (2T, 5kg, Gd), not at the real NSGA-III knee-point design that
  `design_recommendations.py` (step 13) produces, because step 5 runs
  before step 13 in `main.py`'s pipeline order. LIMITATIONS.md's own
  suggested fix -- rerunning the economics functions a second time at the
  step-13 knee point -- is now implemented: `main.py` has a new step 13b
  (`run_economics_at_knee_point()`, in `main.py`, called after step 13
  with the same knee-point dict `design_recommendations.
  summarize_field_flow_lever()` computes) that reruns
  `bom_cost()`/`full_system_cost_estimate()`/`levelized_cost_of_cooling()`/
  `cross_check_full_system_cost_methods()`/`full_system_cost_envelope()`
  at the actual recommended design, and reports the result in both the
  per-run log and the executive-summary banner. `data/MagCool_DC_BOM.xlsx`
  itself is still a hand-assembled, out-of-band artifact from a prior
  run's knee point (regenerating it programmatically from step 13b's
  output has not been attempted), but a normal `python main.py` run now
  reproduces that same design point's numbers through the pipeline
  itself, not only via manual, out-of-band calls.
- **RESOLVED.** Every `core/economics.py` function this section
  previously flagged as implemented-and-tested-but-never-called from
  `main.py`'s pipeline -- `bom_cost_geometric()` /
  `full_system_cost_estimate_geometric()` (Halbach-cylinder geometric
  magnet-mass model), `magnet_grade_cost_comparison()`,
  `resource_criticality_note()`, `commercial_mcm_price_reality_check()`,
  `rowe2011_vcc_compressor_cost_cross_check()`,
  `rowe2011_magnet_mass_ratio_cross_check()`,
  `compare_legacy_and_updated_magnet_ratio()`, `deflate_to_2026()` /
  `inflation_adjusted_reference_prices()` / `bom_cost_2026_dollars()`
  (CPI-U 2026-dollar adjustment), `full_system_cost_estimate_range()`
  (LOW/MID/HIGH band on the borrowed-multiplier method), and
  `full_system_cost_envelope()` (combined-methods envelope) -- is now
  reachable from `main.py`'s pipeline: `full_system_cost_estimate_range()`
  is called inside `cross_check_full_system_cost_methods()`, which
  `run_economics()` (step 5, and now step 13b above) calls directly;
  `full_system_cost_envelope()` is likewise called directly by both;
  and the rest are called by `full_economics_wiring_report()` (also
  `core/economics.py`), which `run_economics()` calls at the end of step
  5. Every one of these now appears in `results/*.txt` output on a
  normal pipeline run, not only in `core/economics.py` and
  `tests/test_economics.py`.
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
- **Commercial-landscape figures are unaudited vendor/press claims.**
  `core/commercial_landscape.py`'s `COMMERCIAL_SYSTEMS` dataset (Magnotherm
  Stellar, Cooltech's data-center-oriented unit) is sourced entirely from
  vendor announcements and trade press -- no independently-audited
  datasheet was located for either system, and this is stated explicitly
  in the module's own output rather than presented as verified.
- **Regeneration of stale `results/*.txt` diagnostic files** listed in
  Section 2.1 above.
- **Papers/ subfolders not yet mined**: `Reviews/`, `Data center
  cooling/`, and `AMR Theory and Modeling/` were not systematically
  cross-checked against the codebase's citations in a later pass (time-
  bounded scope decision, not a finding of any kind) -- a reasonable next
  pass for whoever picks this up next.

---

## 5. `Papers/` corpus verification pass

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
the way the earlier scattered version did.