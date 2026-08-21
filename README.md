# magcool-dc

Physics-based simulation suite evaluating **magnetocaloric (magnetic) cooling
for data centers**, benchmarked against vapor-compression CRAC/CRAH and
direct liquid cooling.

## Index

This file is a running log, not a reorganized document — sections are
mostly in the order they were written, and (see note below) that is
**not** the same as chronological phase order. Use this index rather
than scrolling.

**Start here:**
- [Status](#status) — current phase, one-paragraph-per-phase changelog (long)
- [What's implemented](#whats-implemented) — per-module table, the fastest way to see what a given file does and its most recent change
- [Quick start](#quick-start) — how to run it
- [Repo layout](#repo-layout) — directory map, file-by-file

**Findings write-ups, in file order (see chronology note below):**
- [Phase 15 findings](#phase-15-findings) — BOM cost model, Hypereg pumping analysis, material+geometry co-optimization
- [Phase 16 findings](#phase-16-findings) — thermal-hysteresis loss, NSGA-III front-composition sensitivity
- [Phase 17-23 findings](#phase-17-23-findings) — cycle-type switch (Brayton/Ericsson/Carnot), thermal diode, Halbach magnet geometry, magnetocaloric fluids, passive/hybrid regenerators, inhomogeneous broadening + nanocomposite materials, elastocaloric reference line
- [1-D regenerator model findings](#1d-regenerator-model-findings) — **most recent work**: the regenerative-amplification gap, the low-mdot fix, the opt-in `no_load_span_override`, and the 4-of-5-device override check
- [Phase 12 findings](#phase-12-findings) — paper-mining pass, part 3
- [Phase 11 findings](#phase-11-findings) — paper-mining pass, part 2
- [Phase 10 findings](#phase-10-findings) — paper-mining pass
- [Phase 7 findings](#phase-7-findings) — (so far)
- [Phase 6 findings](#phase-6-findings)
- [Phase 5 findings](#phase-5-findings)
- [Phase 4 findings](#phase-4-findings)
- [Phase 3 findings](#phase-3-findings)
- [Key Phase 2 finding: ideal vs. electrical COP](#phase-2-finding)
- [Sobol sensitivity: a genuine model-structure finding](#sobol-sensitivity)
- [Validation snapshot (Gd MCE, mean-field vs. Dan'kov et al. 1998)](#validation-snapshot)

**Chronology note:** phase numbers do NOT increase monotonically as you
scroll down this file. Phases 3-12 (paper-mining passes) were written
*before* Phases 15-23 and the 1-D regenerator work, but sit physically
*below* them in this file — they were the earliest, most foundational
passes and were left in place rather than reordered, since reordering
a running log risks breaking existing external links/citations into it.
If you want strict chronological order, read bottom-to-top starting
from [Key Phase 2 finding](#phase-2-finding), then jump up to
[Phase 15](#phase-15-findings) onward. If you just want the current
state of a specific module, [What's implemented](#whats-implemented)
is faster than reading findings sections in any order.

<a id="status"></a>
## Status: Phase 13 (paper-mining pass part 6: traced/corrected the `DTU_rotary_Gd_2016` citation, promoted `DTU_Eriksen_rotary_Gd_2015` into the CORE calibration slot) done, Phase 14 (bug fixes: corrected a mean-field-vs-first-order GD5SI2GE2 material mixup in the cascade comparison/fig20/cascade.py demo; added `core/material_family_comparison.py`, a four-way Gd/Gd5Si2Ge2/GD-family/LAFESIH-family/MNFEPSI-family ranking at the same ASHRAE point, wired into `main.py` step 8d and `plots.py` fig26; documented the `span_fraction` linear-clamp approximation rather than inventing an unsourced smoothing function; confirmed the full-system BOM cost gap, reference-book OCR, and two flagged CSV rows are already correctly left open/flagged, no change needed; Tušek et al. (2013) Figs. 10-11 digitization still open, being done manually in WebPlotDigitizer; 170/170 tests passing) done, Phase 15 (a full-system BOM cost model in `core/economics.py` — soft-magnetic-yoke cost, an order-of-magnitude full-system-cost estimate, and a CRF-based levelized cost of cooling; confirmed the multi-bed-rotary-vs-reciprocating loss question was already answered by existing `RotaryDriveLossModel`/`analyze_parasitic_fraction_scaling` infrastructure, no new term needed; a Hypereg-style parallel-hydraulic pumping-power analysis (`core/hypereg_analysis.py`), grounded in a direct reading of Klinar et al. (2024); and material+geometry co-optimization inside NSGA-III (`core/optimize.py`) — particle diameter as a 7th design variable wired through `core/amr_cycle.py`'s new geometry-aware pumping-power accounting, and material family as a per-family NSGA-III search merged post-hoc into one Pareto front; 153/153 tests passing. **Note**: `core/design_recommendations.py`, referenced in the original Phase 15 plan, does not exist in this project snapshot — see `ROADMAP.md` Phase 15 item 1 for the flagged discrepancy) done, Phase 16 (thermal-hysteresis loss quantified for the first time: new `hysteresis_loss_J_per_kg` field on `FirstOrderMCEMaterial` [core/first_order_mce.py], literature-analog placeholder values for all three first-order families (each heavily honesty-flagged), wired into `AMRSystem.run()` [core/amr_cycle.py] as an additional parasitic-power term; new `core/hysteresis_sensitivity.py` A/B diagnostic asks whether Phase 15's "100% La(Fe,Si)13Hy" merged-front result survives — result: front composition DID shift, but in the opposite direction than naively expected (MORE La(Fe,Si)13Hy-dominant, not less), confirmed to be a genuine NSGA-III search-dynamics effect rather than a bug via an independent fixed-design-point sanity check; open item flagged for a full-resolution rerun; 216/216 tests passing) done, Phase 17 (AMR cycle-topology switch — `cycle_type` on `AMRSystem` [core/amr_cycle.py], `"brayton"` [default, unchanged behavior] / `"ericsson"` / `"carnot"` — with the honest caveat that this project's own copy of Kitanovski et al. (2015) excludes pp. 104-109 [Sect. 4.1.1-4.1.4], so the `CYCLE_TYPE_FACTORS` multipliers encode only the qualitative Carnot >= Ericsson >= Brayton ranking, not the book's own closed-form relations; new `validation_system.run_cycle_type_validation()` re-checks rotary-named benchmark devices as Ericsson-like, finding one genuine improvement [`DTU_Eriksen_rotary_Gd_2015`: COP error -2.1% -> +0.6%] out of two comparable rotary devices; 236/236 tests passing) done, Phase 18 (a scoped-down mechanical-contact active thermal diode — `core/thermal_diode.py`'s `MechanicalContactDiode` and `core/thermal_diode_analysis.py` — after checking the plan's own premise directly and finding `AMRSystem` has no internal frequency ceiling for a diode to "unlock", so the module ships as a parasitic-cost-only sensitivity study [no benchmark device uses thermal diodes] rather than the fuller NSGA-III categorical variant originally proposed; 254/254 tests passing) done, a follow-up pass closing three flagged Phase 16-18 open items (a full-production-settings multiseed hysteresis-reversal check that found the Phase 16 reversal is **not** stable at production NSGA-III settings — `results/hysteresis_multiseed_stability.txt`; a higher-mdot Hypereg re-sweep confirming a >2x larger benefit at 0.3 kg/s than at the 0.08 kg/s baseline; and a `cycle_type` sensitivity check on the Astronautics graded-bed cascade, which found the Ericsson reclassification does **not** narrow that device's much larger -81.1% error, unlike the single-bed DTU case) done, Phase 19 (a closed-form Halbach-cylinder field-vs-mass magnet geometry model — `core/magnet_geometry.py`, `economics.bom_cost_geometric()` — replacing the old flat $/T magnet-mass ratio with a genuinely super-linear-in-field relation, wired into `optimize.cost_index()` as an opt-in `use_geometric_magnet_mass` flag; a citation correction found and documented for the Bjørk et al. field-vs-cost reference; 293/293 tests passing) done, Phase 20 (magnetocaloric fluids as a new working-body class — `core/fluid_mce_cycle.py`'s `FerrofluidMCESystem`, `core/fluid_mce_analysis.py` — built on standard Krieger-Dougherty/Darcy-Weisbach relations since neither source book's fluids chapter was available; headline finding: the lack of regeneration collapses usable span to well under a Kelvin versus solid AMR at the same field/flow, despite an interior COP optimum existing around phi≈0.10-0.20; 313/313 tests passing) done, Phase 21 (passive/hybrid magnetic regenerators — `core/passive_regenerator_analysis.py`, `baseline_cooling.augmented_regenerator_cop()` — recombining this repo's own existing heat-capacity and baseline-COP models rather than digitizing Tishin's unreadable Ch. 11; found Gd's own Curie point sitting inside the ASHRAE operating window gives the largest, still-modest COP gain [+0.24%], shrinking as span widens; 326/326 tests passing) done, **Phase 22 (three sub-items: item 1 — Gaussian inhomogeneous/polycrystalline Tc-broadening [`core/inhomogeneous_broadening.py`], which narrows the model's 1T ΔT_ad overestimate but widens its 5T underestimate, a genuine trade-off rather than a clean fix; item 2 — an engineered multi-phase nanocomposite material family [`core/nanocomposite_material.py`], which underperforms a perfectly-tuned single phase at its own design span but uniquely survives at an off-design span where the sharply-tuned single phase collapses to zero capacity; item 3 — magnetoelastic/anisotropic contributions deliberately left unmodeled [no digitizable data, no plan deliverable] plus a qualitative amorphous-materials cost/performance note in `core/economics.py`; 353/353 tests passing)** done, **Phase 23 (an elastocaloric reference comparison line — `baseline_cooling.elastocaloric_reference_cop()` — added as a static literature-sourced horizontal band on `plots.py` fig08 and in `main.py`'s baseline sweep, exactly as scoped ["a comparison row, not a new simulated device"]; representative COP=5.8 anchor from Qian et al. 2023's simulated multimode elastocaloric system, low end from Wu et al. 2023's measured narrow-span device, since neither of this repo's two source books covers elastocalorics; 356/356 tests passing)** done, and a diagnostic-only follow-up (`core/regenerator_1d.py`, a 1-D transient blow-by-blow AMR simulator checking the 2d regenerative-amplification gap directly: fixed a low-mdot degeneracy by adding axial conduction the model previously lacked, which exposed a still-open, directionally-inconsistent calibration gap; wired into `main.py` as step 2e, additive/diagnostic only, not feeding `cooling_capacity()` or any downstream result by default; PLUS an opt-in `no_load_span_override` on `AMRSystem` + `regenerative_span_cap()`, letting that fix actually be used per design point, and `run_regenerative_amplification_override_check()` (step 2f) honestly showing it recovers usable COP predictions within ~±16% on 4 of 5 previously-structurally-infeasible benchmark devices, and doesn't help on the 5th) done — see `ROADMAP.md`

<a id="whats-implemented"></a>
## What's implemented

| Module | Purpose |
|---|---|
| `core/mce_material.py` | Mean-field (Brillouin/Weiss) model for continuous-transition materials (Gd, La0.7Ca0.3MnO3); Gd validated, Gd5Si2Ge2 flagged as invalid for this framework |
| `core/first_order_mce.py` | Extended (6th-order) Landau model for first-order/giant-MCE materials: Gd5Si2Ge2, La(Fe,Si)13Hy, and (Mn,Fe)2(P,Si) — **Phase 16: added `hysteresis_loss_J_per_kg` field (literature-analog placeholder values, heavily honesty-flagged per family; 0.0 default preserves old behavior). This session: `_equilibrium_m()` gained an exact closed-form fast path for h==0 (roughly half of all calls, since `H_initial=0.0` is the default nearly everywhere), replacing a general degree-5 `np.roots()` eigenvalue solve with a quadratic formula. Profiled as the dominant cost of `core.cascade.run_graded_cascade()` (called heavily by `main.py`'s NSGA-III optimization stages, notably 11f) — this alone cut its per-call runtime by ~30% (verified: 1.02s → 0.72s/call). Verified bit-for-bit equivalent in the only quantity ever used downstream (`m**2`, via `delta_S_isothermal()`) across 1600 test cases spanning all four material families; the two solutions' raw signs can differ (a physically meaningless ±m degeneracy at zero field, spontaneous symmetry breaking — the original `np.roots()` path had no defined sign convention here either) but this is never observed since sign is never used** |
| `core/giant_mce_analysis.py` | Formal Gd vs. giant-MCE comparison → `results/giant_mce_analysis.txt` |
| `core/material_family_comparison.py` | Four-way material family ranking (Gd, Gd5Si2Ge2-fixed, and the three composition-tunable families) at the same ASHRAE point → `results/material_family_comparison.csv`/`.txt`, fig26 (Phase 14) |
| `core/emissions.py` | Refrigerant-free GWP/emissions comparison |
| `core/amr_cycle.py` | 0-D AMR cycle model: cooling capacity, ideal/electrical COP, optional NTU-derived effectiveness, optional blow-fraction asymmetry (Phase 10) — **Phase 15: optional `particle_diameter`/`bed_cross_section_area`/`hypereg_n_parallel` params wire regenerator geometry into both the NTU effectiveness calculation and a geometry-explicit pumping-power term that replaces (not adds to) the generic loss-model `k_pump` term; `None` by default, fully backward-compatible. Phase 16: `_hysteresis_power_W()` adds `hysteresis_loss_J_per_kg * mass_regenerator * frequency` to `W_parasitic` unconditionally (both the loss_model and constant-parasitic_fraction paths); 0.0 for GADOLINIUM. Phase 17: optional `cycle_type` (`"brayton"` default/unchanged, `"ericsson"`, `"carnot"`) applies `CYCLE_TYPE_FACTORS` multipliers to `cooling_capacity()`'s Qc and `magnetic_work()`'s `eta_2nd_law` — illustrative, qualitatively-ordered placeholders, not digitized from Kitanovski et al.'s own Sect. 4.1.1-4.1.4 (pages not in this project's copy of the book). Phase 18: optional `thermal_diode` param (`None` default) adds `_diode_switching_power_W()` to `W_parasitic`, cost-only (no offsetting heat-transfer benefit modeled — see `core/thermal_diode.py`)** |
| `core/hysteresis_sensitivity.py` | **Phase 16 addition**: A/B comparison of `core.optimize.run_optimization()`'s merged Pareto front with material hysteresis loss on vs. forced off, to check whether Phase 15's "100% La(Fe,Si)13Hy" finding is robust to it → `results/hysteresis_sensitivity.txt`. **Phase 16-18 follow-up**: `run_hysteresis_multiseed_stability_check()` reruns the A/B check at full production NSGA-III settings across 3 seeds — the original reduced-setting reversal does **not** hold up → `results/hysteresis_multiseed_stability.txt` |
| `core/thermal.py` | NTU packed-bed regenerator effectiveness model, packed-bed/parallel-plate pumping power (Tušek et al. 2013) — **Phase 15: added `pumping_power_packed_bed_hypereg()`, a parallel-sub-regenerator pressure-drop-reduction variant motivated by Klinar et al. (2024). Phase 21: `regenerator_effectiveness()` gained an optional `cp_solid` override (default `None` = unchanged, fixed Gd-only behavior) so a caller can isolate a material's own magnetic heat-capacity anomaly, used by `passive_regenerator_analysis.py`** |
| `core/hypereg_analysis.py` | **Phase 15 addition**: demonstrates the Hypereg parallel-hydraulic pumping-power benefit at this repo's representative operating point → `results/hypereg_analysis.txt`, `results/hypereg_findings.md`. **Phase 16-18 follow-up**: `sweep_n_parallel_at_higher_mdot()` re-checks the benefit at 0.3 kg/s (~4x baseline), finding it >2x larger there |
| `core/loss_model.py` | State-dependent eddy/pumping/base loss model — **Phase 6: added a 4th benchmark device and found the extended fit is unstable (negative coefficients, leave-one-out errors up to +1639%); production default stays on the stable 3-point CORE fit, instability documented via `run_extended_diagnostic()`. Phase 15: `parasitic_power()` gained a `pumping_power_override` parameter (default `None` = unchanged behavior) so `amr_cycle.py` can substitute a geometry-explicit pumping term without double-counting against the generic `k_pump` term; confirmed the rotary-vs-reciprocating loss-topology question was already closed by the existing `RotaryDriveLossModel`/`analyze_parasitic_fraction_scaling`, documented rather than duplicated. This session: added `CALIBRATION_POINTS_CORE_PLUS_MAGGIE_HIGHSPAN`, a genuine 4th point unlike EXTENDED/FURTHER_EXTENDED's larger-scale-jump devices — it's the SAME physical prototype as an existing CORE point at a different condition, only reachable via `amr_cycle.py`'s new `no_load_span_override`. Honest leave-one-out result: its own held-out fold predicts within +30% (vs. ~250-700% for every other fold, old or new), and it modestly improves rather than degrades the two folds it shares with CORE — real, partial progress, but Astronautics and Tušek's folds remain order-of-magnitude misses either way, so CORE (3pt) stays the production default. `run_core_plus_maggie_highspan_diagnostic()`, wired into `main.py` as step 3a2, additive** |
| `core/thermal_diode.py` | **Phase 18 addition**: `MechanicalContactDiode` dataclass (forward/reverse conductance, `rectification_ratio`, `switching_power_W()`, illustrative `actuation_energy_J_per_cycle` — round-number placeholder, since Kitanovski et al.'s own Sect. 6.2.4 pages aren't in this project's copy of the book) and `cycle_time_reduction_factor()`; rectification ratio (`DEFAULT_MECHANICAL_CONTACT_DIODE`) later grounded against Bywaters & Griffin's cryogenic piezo-actuated gas-gap heat-switch data |
| `core/thermal_diode_analysis.py` | **Phase 18 addition**: `check_frequency_ceiling_claim()` (finds `AMRSystem` has no internal frequency cap for a diode to unlock — NSGA-III's 5 Hz search bound is an unexplained round number, not a physical limit) and `sweep_frequency_with_and_without_diode()` (no benchmark device uses thermal diodes, so this is a within-model cost-only sensitivity study) → `results/thermal_diode_analysis.txt` |
| `core/magnet_geometry.py` | **Phase 19 addition**: closed-form idealized-Halbach-cylinder field-vs-mass relation (`halbach_bore_field_T()`, `halbach_magnet_mass_kg()`, `halbach_field_vs_mass()` — genuinely super-linear in field, unlike the old flat $/T ratio), since Kitanovski et al.'s own Ch. 3 (Magnetic Field Sources) isn't in this project's copy of the book; `bjork_qualitative_check()` honestly reports its own simple fixed-mass cost-per-Kelvin proxy does **not** reproduce the literature's ~2T cost optimum → `results/magnet_geometry_analysis.txt` |
| `core/fluid_mce_cycle.py` | **Phase 20 addition**: `FerrofluidMCESystem`, a sibling to `AMRSystem` for magnetocaloric fluids/ferrofluids as a continuous-flow working body (no packed-bed regenerator) — `krieger_dougherty_viscosity()`, `suspension_delta_T_adiabatic()` (mixture heat-capacity dilution), `pumping_power_pipe_flow()` (Darcy-Weisbach); `eta_2nd_law_fluid` is an explicitly uncalibrated constant (no benchmark device exists for this working-body class) |
| `core/fluid_mce_analysis.py` | **Phase 20 addition**: `volume_fraction_sweep()` (finds an interior COP optimum near phi≈0.10-0.20) and `compare_to_solid_amr_and_liquid_cooling()` (headline finding: no regeneration collapses usable span to well under a Kelvin vs. solid AMR at the same field/flow) → `results/fluid_mce_analysis.txt` |
| `core/passive_regenerator_analysis.py` | **Phase 21 addition**: `compare_candidate_materials()`/`span_sweep()` for passive/hybrid magnetic regenerators — a conventional gas cycle's internal regenerator loaded with an MCE material to exploit its Curie-point heat-capacity anomaly, built by recombining this repo's own existing material and baseline-COP models (Tishin & Spichkin's relevant chapter is an image-only PDF, unreadable here) → `results/passive_regenerator_analysis.txt` |
| `core/inhomogeneous_broadening.py` | **Phase 22 item 1 addition**: `BroadenedMagnetocaloricMaterial` (Gauss-Hermite-quadrature ensemble over a Gaussian-distributed grain-to-grain Curie temperature, `with_Tc()` on `MagnetocaloricMaterial`) — sweeping `sigma_Tc` narrows the mean-field model's 1T ΔT_ad overestimate but widens its 5T underestimate, a genuine trade-off rather than a clean fix → `results/inhomogeneous_broadening.txt` |
| `core/nanocomposite_material.py` | **Phase 22 item 2 addition**: `WeightedMaterialEnsemble`/`nanocomposite_tuned_material()` — a 3-phase triangular-weighted blend of composition-tuned La(Fe,Si)13Hy phases mixed at the ΔT_ad level (not entropy/heat-capacity, since `dTad_correction` is a whole-ratio correction); underperforms a perfectly-tuned single phase at its own design span, but uniquely survives at an off-design span where the single phase collapses to zero → `results/nanocomposite_robustness.txt` |
| `core/optimize.py` | NSGA-III multi-objective optimization — **Phase 15: 7 design variables (field, frequency, flow, mass, effectiveness, blow fraction, + new particle diameter), material now co-optimized as a design choice (Gd + 3 composition-tunable giant-MCE families, each run separately through NSGA-III and merged post-hoc into one globally non-dominated Pareto front — see `ROADMAP.md` Phase 15 for why "merge post-hoc" was chosen over a native mixed-variable formulation); cost objective upgraded to `economics.bom_cost()` (family-specific MCM pricing + SMM yoke term). Phase 19: `cost_index()` gained an opt-in `use_geometric_magnet_mass=False` flag to use the super-linear Halbach-geometry magnet-mass model instead of the flat $/T ratio** |
| `core/cascade.py` | Multi-stage cascade AMR design; Curie-graded beds pluggable across Gd5Si2Ge2, La(Fe,Si)13Hy, and (Mn,Fe)2(P,Si) families — **Phase 16-18 follow-up: `cycle_type` threaded through `run_graded_cascade()`/`validate_astronautics_graded_bed()`; `run_astronautics_cycle_type_sensitivity()` finds the Ericsson reclassification does NOT narrow the Astronautics graded-bed's -81.1% error, unlike the single-bed DTU case** |
| `core/baseline_cooling.py` | Vapor-compression and liquid-cooling COP correlations — **Phase 21: `augmented_regenerator_cop()` models a passive/hybrid magnetic regenerator boosting a conventional cycle's COP (see `core/passive_regenerator_analysis.py`). Phase 23: `elastocaloric_reference_cop()` adds a static, literature-sourced (not simulated) elastocaloric COP reference point/range for comparison** |
| `core/economics.py` | CAPEX/OPEX comparison, materials-cost floor grounded in Bjørk et al. (2011) — **Phase 15: added `bom_cost()` (+ soft-magnetic-material yoke cost, Silva et al. 2017), `full_system_cost_estimate()` (order-of-magnitude full-system multiplier from Russek & Zimm 2006's vapor-compression-AC manufactured-cost benchmark), `levelized_cost_of_cooling()` (CRF-based $/kWh, a second independent cost methodology), and `MCM_COST_PER_KG_BY_FAMILY` (per-material-family MCM pricing, e.g. La(Fe,Si)13Hy at $8/kg). Phase 19: added `geometric_magnet_mass_kg()`/`bom_cost_geometric()`/`full_system_cost_estimate_geometric()`, additive siblings using the Halbach-geometry mass model. Phase 22 item 3: added a qualitative `amorphous_material_cost_performance_note()` (not wired into any priced dict — no sourced $/kg or ΔS_M figure exists in this repo's corpus)** |
| `core/material_family_comparison.py` | Material family ranking at the same ASHRAE point → `results/material_family_comparison.csv`/`.txt`, fig26 (Phase 14) — **Phase 22 item 2: extended from four to six candidates with `NANOCOMPOSITE_FAMILY` added to `TUNABLE_FAMILIES`** |
| `core/regenerator_1d.py` | **Diagnostic addition, follow-up to the 2d gap-quantification step**: a standalone 1-D transient (blow-by-blow, multi-cycle) AMR simulator — the actual regenerator physics `amr_cycle.py`'s 0-D `cooling_capacity()` cannot represent, built to check whether an explicit simulation reproduces the spans 2d shows the 0-D model is structurally capped below. Fixed an initial low-mdot degeneracy (predicted span grew without bound as flow decreased, no interior maximum) by adding axial node-to-node conduction, which the model had lacked entirely. That fix exposed a follow-on problem: the model's error is now directionally inconsistent across the three directly-measured benchmarks (over on one device, under on two others), because the axial-conductivity value used is physically reasoned but not independently calibrated. Its `regenerative_span_cap()` feeds `amr_cycle.py`'s new opt-in `no_load_span_override` (default `None`, every other caller unaffected) — `validation_system.run_regenerative_amplification_override_check()` (step 2f) shows this recovers usable COP predictions within ~±16% on 4 of 5 benchmark devices the old cap makes structurally infeasible, and honestly doesn't help on the 5th. **Neither is used by default anywhere else in this repo — `cooling_capacity()`'s default behavior is unmodified.** Wired into `main.py` as steps 2e/2f, additive/diagnostic-only → `results/regenerator_1d_validation.txt`, `results/regenerative_amplification_override_check.txt` |
| `core/validation.py`, `validation_system.py`, `giguere_validation.py` | Material- and system-level validation against literature/prototypes; `validation.py` also cross-checks Gd at 7T (Giguere et al.) and a Dan'kov et al. Curie-shift held-out prediction (not reproduced by the model, documented); `validation_system.py` also field-sensitivity-checks the Chubu Electric/Toshiba pair and reachability-checks capacity-only rows (Cooltech's 42K stress test); `giguere_validation.py` also cross-checks the Gd5Si2Ge2 dTad correction against Pecharsky & Gschneidner (1997)'s independent peak-ratio figure — **Phase 17: `validation_system.py` gained `infer_cycle_type_for_device()`/`run_cycle_type_validation()`, finding the Ericsson reclassification narrows `DTU_Eriksen_rotary_Gd_2015`'s COP error from -2.1% to +0.6%** |
| `core/sensitivity.py`, `rsm.py` | Sobol sensitivity, RSM surrogate |
| `main.py` | Full comparison across the ASHRAE TC9.9 thermal envelope: 15 numbered stages (plus sub-stages, e.g. 1b/2b/3b-3d/7b-7c/8b-8e/9b/11b-11d, and 2e/2f/3a2 for the 1-D regenerator diagnostics above) and 35 figures (fig35, this session, shows the regenerative-amplification override check visually for the first time — see the "1-D regenerator model findings" section) — see `python main.py --help`-style stage list in the module's own docstring. **`--quick`** skips 2e/2f/3a2 (each a multi-cycle transient regenerator simulation, ~30-90s per device; all three are additive/diagnostic-only, so skipping them changes nothing about any other stage's output). Even without `--quick`, `core/regenerator_1d.py`'s `no_load_span()` disk-caches every simulation (`results/.regenerator_1d_cache.json`, gitignored), so repeated runs only pay the simulation cost once per unique input combination |

<a id="quick-start"></a>
## Quick start

```bash
pip install -r requirements.txt
python -m core.validation            # MCE material model vs. literature
python -m core.validation_system      # AMR system model vs. published prototypes (point-wise + curve-level)
python -c "from core.validation_system import run_cycle_type_validation as f; f()"  # Phase 17: rotary-device Ericsson-like cycle_type sensitivity
python -m core.first_order_mce         # first-order Landau model calibration check
python -m core.giant_mce_analysis       # Gd vs. giant-MCE, formal comparison
python -m core.material_family_comparison # four-way material family ranking at the ASHRAE point
python -m core.emissions                 # refrigerant-free GWP/emissions case
python -m core.loss_model                 # CORE calibration + Phase 6 EXTENDED diagnostic
python -m core.thermal                     # NTU regenerator effectiveness sweeps
python -m core.hypereg_analysis             # Phase 15: Hypereg parallel-hydraulic pumping-power sweep
python -m core.sensitivity                  # Sobol, Phase 2 vs. Phase 3 modes
python -m core.rsm                           # RSM surrogate for cooling capacity
python -m core.optimize                       # NSGA-III Pareto front — Phase 15: material + geometry co-optimized, grounded BOM cost model
python -m core.cascade                         # multi-stage cascade, Gd vs. Gd5Si2Ge2
python -m core.hysteresis_sensitivity          # Phase 16: hysteresis-on-vs-off A/B check of the Phase 15 material finding
python -c "from core.hysteresis_sensitivity import run_hysteresis_multiseed_stability_check as f; f()"  # Phase 16-18 follow-up: multiseed reversal-stability check at production settings
python -m core.thermal_diode_analysis          # Phase 18: mechanical-contact thermal-diode frequency-ceiling + cost sensitivity study
python -m core.magnet_geometry                  # Phase 19: Halbach-cylinder field-vs-mass geometry model + cost-per-Kelvin sweep
python -m core.fluid_mce_analysis                # Phase 20: magnetocaloric-fluid (ferrofluid) working-body sweep + solid-AMR/liquid-cooling comparison
python -m core.passive_regenerator_analysis       # Phase 21: passive/hybrid magnetic regenerator augmentation of a conventional gas cycle
python -m core.inhomogeneous_broadening            # Phase 22 item 1: Gaussian Tc-broadening sensitivity of the mean-field Gd model
python -m core.nanocomposite_material               # Phase 22 item 2: nanocomposite (multi-phase blended) La(Fe,Si)13Hy family + off-design robustness check
python main.py                                       # full pipeline: validation, economics/BOM, cascade, sensitivity, optimization, 35 figures
python main.py --quick                                # same, but skips the slow 1-D regenerator diagnostic stages (2e/2f/3a2) -- see "What's implemented"'s main.py row
```

### Phase 15 module usage notes

```python
# Geometry-aware AMR system (backward-compatible: omit particle_diameter for old behavior)
from core.amr_cycle import AMRSystem
from core.mce_material import GADOLINIUM
from core.loss_model import StateDependentLossModel

sys_ = AMRSystem(GADOLINIUM, mu0H_max=1.5, mass_regenerator=5.0, frequency=1.0,
                  fluid_mdot=0.08, regenerator_effectiveness=0.85,
                  loss_model=StateDependentLossModel(), use_ntu_thermal_model=True,
                  particle_diameter=0.0005,      # 0.5mm -- feeds NTU eps AND replaces
                                                   # the generic k_pump term with a
                                                   # geometry-explicit pumping power
                  hypereg_n_parallel=4)            # optional: Hypereg-style parallel split
result = sys_.run(T_cold=291.0, T_span=10.0)

# Full-system BOM cost model
from core import economics
bom = economics.bom_cost(mu0H_max=1.5, mass_regenerator=5.0, family_name="Gd")
full_system = economics.full_system_cost_estimate(1.5, 5.0, family_name="Gd")
lcoc = economics.levelized_cost_of_cooling(1.5, 5.0, Qc_avg_W=1000, COP_electrical=5.0)
```

Running `python -m core.optimize` (or `main.py`'s step 11) now writes
`results/pareto_front.csv` (merged, globally non-dominated across all
material candidates) AND `results/pareto_front_by_material/<material>.csv`
(each material's own front, before merging, for transparency).

<a id="phase-15-findings"></a>
## Phase 15 findings

Full writeup and honesty flags in `ROADMAP.md`'s Phase 15 section and
`results/hypereg_findings.md`; summary here.

**Material + geometry co-optimization** (`core/optimize.py`): at this
repo's fixed representative operating point (T_cold=291K, 10K span),
running NSGA-III separately per material candidate and merging the
resulting fronts (option (b) from the plan — see `ROADMAP.md` for why
this was chosen over a native mixed-variable formulation) found the
merged, globally non-dominated Pareto front was **100% La(Fe,Si)13Hy**
(23/23 designs) at production settings — Gd and the composition-tunable
Gd-Si-Ge family were all cross-material-dominated at this particular
point. This is a result of the search, not an assumption baked in;
`results/pareto_front_by_material/*.csv` keeps each material's own
front so the conclusion is checkable rather than asserted.
`particle_diameter` (the new geometry design variable) spans roughly
0.3-1.9mm across the merged front — a real, active search dimension.

**Hypereg** (`core/hypereg_analysis.py`, `results/hypereg_findings.md`):
Klinar et al. (2024)'s Hypereg concept is a **hydraulic** idea (parallel
vs. series sub-regenerator flow reducing pressure-drop length), not an
eddy-current one — implemented as a pumping-power-only variant of the
existing packed-bed correlation. At this repo's lab-scale representative
point, splitting into 4 parallel sub-beds (the paper's own illustrative
example) raised COP_electrical by about 0.2% (5.264→5.275), saturating
around n=16 at 5.278 — real but modest in this model, since pumping power
is only one of three loss channels here and not the dominant one at this
scale.

**Full-system BOM** (`core/economics.py`): at the 2T/5kg Gd design point
used throughout the pipeline, the materials BOM (magnet + MCM + new
soft-magnetic-yoke term) comes to $1,375; an order-of-magnitude
full-system estimate (using a vapor-compression-AC manufactured-cost
benchmark as a sanity-check multiplier, NOT an AMR-specific quote) is
about $13,750; and a CRF-based levelized cost of cooling comes to about
$0.034/kWh_cooling. A genuine bottom-up AMR-specific BOM (HX, pump,
motor, controls) remains open — see `ROADMAP.md`'s Phase 16 candidates.

**Rotary AMR topology**: re-examining the plan's "does loss behavior
differ for rotary AMR / multi-bed topologies" question found it was
already answered by existing infrastructure
(`core.loss_model.RotaryDriveLossModel`,
`analyze_parasitic_fraction_scaling()`) — documented explicitly rather
than duplicated; see `core/loss_model.py`'s module docstring.

<a id="phase-16-findings"></a>
## Phase 16 findings

Full writeup and honesty flags in `ROADMAP.md`'s Phase 16 section;
summary here.

**Hysteresis loss, quantified for the first time** (`core/first_order_mce.py`,
`core/amr_cycle.py`): thermal-hysteresis loss — real, irreversible energy
dissipated each cycle by first-order materials, which Gd genuinely does
not pay — was previously a documented but entirely unquantified honesty
flag (prose caveats only). It is now a real number
(`hysteresis_loss_J_per_kg`, literature-analog values for all three
first-order families, each heavily flagged as approximate — see
`core/first_order_mce.py`) that adds `hysteresis_loss_J_per_kg * mass *
frequency` to `W_parasitic` in `AMRSystem.run()`, lowering
`COP_electrical` for every first-order material while leaving Gd
untouched. A direct fixed-design sanity check confirms the mechanism
itself is correct and correctly-signed (COP_electrical dropped from
9.877→9.212 at one representative design point, exactly matching the
expected 218.2W hysteresis penalty).

**Does Phase 15's "100% La(Fe,Si)13Hy" result survive?**
(`core/hysteresis_sensitivity.py`, `results/hysteresis_sensitivity.txt`):
run at reduced NSGA-III settings (`pop_size=32, n_gen=15` — smaller than
`run_optimization()`'s own 40/25 production default, purely for
wall-clock reasons), the merged front's La(Fe,Si)13Hy share went from
88% (hysteresis off) to **100%** (hysteresis on) — the opposite of the
naively expected direction. This is a genuine NSGA-III search-dynamics
effect (the loss term reshapes each material's own per-family Pareto
front rather than uniformly shifting every point downward, since it
scales with `mass * frequency` rather than being a flat penalty), not a
bug — but it was only checked at reduced search resolution. **Open item**:
rerun `core.hysteresis_sensitivity.run_hysteresis_sensitivity()` at full
production settings (and multiple seeds) to check whether this reversal
is stable or is itself search noise at the reduced setting.

**What this phase deliberately left open**: none of the three
`hysteresis_loss_J_per_kg` values are read directly off the actual
calibrated compositions' own measured hysteresis loops (all three are
literature analogs for related-but-different compositions, weakest for
the Mn-Fe-P-Si family); hysteresis is treated as a fixed per-family
constant rather than varying with tuned Tc within a family; and it is
accounted for as an additive parasitic-power term rather than folded
into the ideal-cycle thermodynamics or exergy efficiency. See
`ROADMAP.md`'s Phase 16 section for the full discussion.

<a id="phase-17-23-findings"></a>
## Phase 17-23 findings

Full writeups and honesty flags live in `ROADMAP.md`'s own Phase 17-23
sections; short summary here. All six phases share the same recurring
constraint: this project's copies of both reference books are heavily
gapped (Kitanovski et al. 2015 is a 30-page front-matter/Ch.1-only
excerpt; Tishin & Spichkin 2003 is a scanned, image-only PDF with no
extractable text). Where a plan named a specific book section as its
source, that gap is stated explicitly rather than papered over, and the
module ships on standard textbook physics or independently-sourced
external literature instead.

**Phase 17 — cycle topology** (`core/amr_cycle.py`'s `cycle_type`):
`"ericsson"` narrows one comparable rotary device's COP error from -2.1%
to +0.6%, but this is a single-device result — not enough to confirm
either the qualitative "rotary → Ericsson" mapping or the specific
multiplier values used.

**Phase 18 — thermal diode** (`core/thermal_diode.py`,
`core/thermal_diode_analysis.py`): checking the plan's own premise found
`AMRSystem` has no internal frequency ceiling for a diode to unlock, and
no benchmark device uses one — so this ships as a parasitic-cost-only
sensitivity study, deliberately not paired with an (undigitizable)
heat-transfer benefit from `rectification_ratio`.

**Phase 19 — magnet geometry** (`core/magnet_geometry.py`): a standard
idealized-Halbach-cylinder relation replaces the old flat $/T magnet-mass
ratio with a genuinely super-linear one. A reduced-resolution spot check
(pop_size=12/n_gen=5) initially suggested the geometric cost term pulls
the optimizer's mean field down; the full production-settings rerun
(pop_size=32/n_gen=15) **reversed** that direction instead — mean field
went *up* (1.31T → 1.54T), not down — the same reduced-settings-vs-
production-settings reversal pattern Phase 16's own hysteresis check hit.
No multi-seed confirmation has been run for this one yet, so treat the
direction as **open**, not settled, until that's done (see `ROADMAP.md`'s
Phase 19 follow-up note). Also corrected a mis-citation this project and
its own plan had both been carrying (the Bjørk et al. field-vs-cost
tradeoff paper is arXiv:1410.6248, not :1410.1987).

**Phase 20 — magnetocaloric fluids** (`core/fluid_mce_cycle.py`,
`core/fluid_mce_analysis.py`): the headline finding is that, absent a
regenerator, the mixture-heat-capacity dilution model collapses usable
span to well under a Kelvin — dramatically less than solid AMR at the
same field/flow — even though an interior COP optimum vs. particle volume
fraction (phi≈0.10-0.20) does genuinely exist in this model.

**Phase 21 — passive/hybrid regenerators** (`core/baseline_cooling.py`'s
`augmented_regenerator_cop()`, `core/passive_regenerator_analysis.py`):
reuses this repo's own existing heat-capacity and baseline-COP models
rather than digitizing Tishin's unreadable chapter. Gd's own Curie point
(294K) sitting inside the ASHRAE operating window gives it the largest —
still modest — COP gain (+0.24%), shrinking as span widens.

**Phase 22 — three sub-items**: item 1's Gaussian Tc-broadening
(`core/inhomogeneous_broadening.py`) narrows the mean-field model's 1T
ΔT_ad overestimate (+48.9%→+20.9% at σ_Tc=5K) but widens its 5T
underestimate (-7.5%→-14.2%) — a genuine trade-off, not a clean fix. Item
2's nanocomposite family (`core/nanocomposite_material.py`) underperforms
a perfectly-tuned single La(Fe,Si)13Hy phase at its own design span, but
is the only candidate that still delivers cooling capacity at an
off-design span where the sharply-tuned single phase collapses to zero.
Item 3 deliberately builds no magnetoelastic/anisotropy model (the plan
itself named no deliverable for it) and adds only a qualitative,
un-numbered amorphous-materials cost/performance note to
`core/economics.py`.

**Phase 23 — elastocaloric reference line** (`core/baseline_cooling.py`'s
`elastocaloric_reference_cop()`): a static, literature-sourced comparison
row/band on fig08 and in the baseline sweep — not a simulated device, per
the plan's own scoping. Representative COP=5.8 from Qian et al. (2023)'s
simulated multimode elastocaloric system; low end 3.7 from Wu et al.
(2023)'s measured but narrow-span (~1K) device.

**Follow-up pass (closing three flagged Phase 16-18 items)**: at full
production NSGA-III settings and multiple seeds, Phase 16's own
hysteresis-reversal finding (88%→100% La(Fe,Si)13Hy share) does **not**
hold up — it was largely a search-noise artifact of the smaller
diagnostic setting (`results/hysteresis_multiseed_stability.txt`).
Hypereg's benefit is meaningfully larger (>2x) at a higher flow rate
(0.3 vs. 0.08 kg/s). And the "rotary → Ericsson" cycle-type reclassification
does **not** generalize to the Astronautics graded-bed device the way it
did for the single-bed DTU device — its much larger -81.1% error is
essentially unchanged.

<a id="1d-regenerator-model-findings"></a>
## 1-D regenerator model findings (diagnostic, not wired into the pipeline)

Full writeup in `core/regenerator_1d.py`'s own module docstring;
summary here. The 0-D `cooling_capacity()` model caps achievable span at
`2*dTad_noload` — the material's own single-blow adiabatic ΔT — and has
no way to represent *regenerative amplification*, the mechanism that
lets real beds reach spans several times that (confirmed directly:
DTU_Eriksen_MAGGIE_2016's measured 29.2K no-load span vs. the 0-D
model's own ~14K structural ceiling at the same field). `regenerator_1d.py`
is a real fix attempt: an actual multi-cycle blow-by-blow transient
simulation of the bed, the kind of model the cited literature uses.

**Bugs found and fixed building it**: an explicit-Euler numerical
instability, an incorrect NTU/dt rescaling that silently killed all heat
transfer, a field-unit error (Tesla passed where the material function
expects A/m) that silently zeroed out the magnetocaloric effect, and a
low-mdot degeneracy where predicted span grew without bound as flow
decreased, with no interior maximum. That last one is now fixed: the
model was missing axial (node-to-node) solid conduction entirely, so as
flow dropped its only inter-node coupling (via the fluid) kept getting
*more* effective with nothing to counteract it. Adding axial conduction
gives the search a genuine interior maximum in mdot instead of an
open-ended one.

**What that fix exposed**: the model's accuracy against the three
directly-measured no-load-span benchmarks is now well-defined but
directionally inconsistent — it overshoots one device and undershoots
two others (+112%, −92%, −61%; see `results/regenerator_1d_validation.txt`).
The axial-conductivity value used is physically reasoned (packed spheres
are point-contact, not a continuous rod, so it's dominated by the
stagnant-fluid path rather than bulk metal conductivity) but not
independently calibrated against a literature source for this specific
geometry — tuning it to fit these three points would just be curve-fitting.

**Wiring**: run as step 2e in `main.py`, additive and diagnostic-only.
It does not feed `cooling_capacity()`, the comparison tables, the Pareto
front, or anything else downstream — an earlier attempt at a different
fix for the same underlying gap (reusing `cascade.py`'s series-stage
machinery with one material sliced into layers) was tested and rejected
because it produced a confidently wrong +220% COP error on a device that
previously had no calibration at all. Wiring in a span predictor that is
still directionally inconsistent would repeat that mistake in a new
form. **Open item**: source an independent packed-bed axial-conductivity
correlation (or measurement) for this geometry before relying on this
model's numbers, or wiring it any further into the pipeline.

**Follow-up: an opt-in override, so the fix is usable, not just
diagnosed.** `AMRSystem` gained `no_load_span_override` (default `None`
— every existing caller, test, and Pareto/Sobol run is byte-for-byte
unaffected). When set to a span (K), it replaces `2*dTad_noload` as
where `cooling_capacity()`'s Qc reaches zero, letting the model reach
spans a single-blow evaluation structurally cannot — populated from the
new `core.regenerator_1d.regenerative_span_cap()` (a thin wrapper around
the transient simulation above; ~30-90s per call, so it's meant to be
computed once per design point and cached, never called from inside a
sweep). `validation_system.run_regenerative_amplification_override_check()`
(wired into `main.py` as step 2f, capped at 3 devices there for pipeline
runtime — pass `max_devices=None` to check every flagged device directly)
tests this honestly against every COP-bearing benchmark row the 2d
diagnostic flags as structurally infeasible under the old cap:

| Device | Old cap | 1-D span cap | COP (lit) | COP (predicted) | Error |
|---|---|---|---|---|---|
| DTU_Eriksen_MAGGIE_2016 | 11.9 K | 21.0 K | 3.60 | 3.04 | −15.7% |
| Risoe_DTU_Gd_2011 | 6.2 K | 3.8 K | 5.00 | — | still infeasible |
| Lozano_POLO_UFSC_2016_r1 | 6.3 K | 23.7 K | 0.37 | 0.32 | −13.9% |
| Lozano_POLO_UFSC_2016_r2 | 6.5 K | 29.8 K | 0.44 | 0.47 | +7.4% |
| Lozano_POLO_UFSC_2016_r5 | 6.3 K | 29.8 K | 0.58 | 0.62 | +7.2% |

4 of 5 previously-impossible (hard `Qc=0`) predictions become usable,
evaluable ones within roughly ±16% — a real result, not a rounding
error, and consistent with the "not yet independently calibrated but
qualitatively sound" status above. The fifth (Risoe_DTU_Gd_2011) is the
directional-inconsistency problem cutting the other way: the 1-D model's
own span cap for that specific device comes out *smaller* than the old
2×dTad_noload cap, so the override doesn't help there at all — exactly
the kind of device-to-device inconsistency the calibration gap above
predicts, shown here on COP rather than span. Still opt-in and still off
by default everywhere else in this codebase for the same reason as
before: this is a usable capability, not yet a validated one.

**Figure**: this whole override check was text-only
(`results/regenerative_amplification_override_check.txt`) until this
session's `fig35` (`core/plots.py`'s `plot_regenerative_amplification_
override_check()`, wired into `main.py`'s step 12 and `run_all()`) —
two panels: span reach (old cap vs. 1-D cap vs. the device's real
reported span, so a device the 1-D cap still doesn't reach is visibly a
shorter bar than "actual span," not silently omitted) and COP recovered
vs. literature for whichever devices the override does reach.

<a id="phase-12-findings"></a>
## Phase 12 findings (paper-mining pass, part 3)


Went through everything not yet touched in Parts 1-2 (the economics paper
in full, both "rotary refrigerator" papers, the 1997 discovery paper, the
solid-state caloric cooling review, both reference books) — see
`Paper_Mining_Recommendations_Part3.md`. Full details in `ROADMAP.md`'s
Phase 12 section.

- **Cooltech 2013 (42K span stress test) and DTU MagQueen (LAFESIH
  cross-check)** added to `data/amr_experimental_benchmarks.csv`, both
  numbers confirmed directly against the source PDF. Cooltech's 42K span
  is the largest in this benchmark set; it does not calibrate at any flow
  rate — consistent with the model's existing struggles at large spans.
  MagQueen is genuinely a heat pump (Qh=1500W, COP_h=5, not Qc/COP_c) —
  its CSV row's Qc=1200W/COP=4.0 are derived via the standard Qh=Qc+W
  identity, clearly flagged as derived rather than measured. Added
  `run_capacity_only_calibration_check()` so capacity-only rows (which
  `run_system_validation()` silently skips) get an actual reported result.
- **Cross-checked the Gd5Si2Ge2 ΔT_ad correction factor** against
  Pecharsky & Gschneidner (1997)'s own ~1.30 peak-ratio figure — a second
  independent source and field range. Genuinely interesting result: the
  raw (uncorrected) model's ratio at 5T (~1.24) is close to 1.30, but the
  Giguere-corrected model's ratio (~0.51) predicts Gd5Si2Ge2
  *underperforms* plain Gd — evidence the correction factor (fit at a
  single 7T point) shouldn't be treated as field-independent.
- **Documentation-only findings**: confirmed the economics paper
  (Bjørk et al. 2011) is fully mined; identified which of two similarly-
  named "rotary refrigerator" papers is design-only (no performance
  numbers) and corrected a mis-attributed citation in `Literature_Review.md`
  as a result; flagged (not added) an unconfirmed 18K span number for the
  Astronautics device and two under-documented reference books.

<a id="phase-11-findings"></a>
## Phase 11 findings (paper-mining pass, part 2)

A deeper second pass through the remaining unmined papers (tables, in-text
numeric callouts, secondary reviews' own citation tables) — see
`Paper_Mining_Recommendations_Part2.md`. Full details in `ROADMAP.md`'s
Phase 11 section.

- **Two free, zero-new-sourcing checks from papers already in the repo**:
  extended `validation.py`'s Gd checks to 7T using Giguere et al. (1999)'s
  own pure-Gd cross-check paragraph (the model overestimates relative to
  it, a real cross-paper literature disagreement — reported, not hidden),
  and added a Curie-point field-shift check against Dan'kov et al.
  (1998)'s reported ~6 K/T rate. The shift check does **not** pass: the
  model's own emergent peak-ΔT_ad temperature is pinned at ~294.5K
  regardless of field (fitted ~0 K/T), confirmed with a sub-Kelvin-
  precision optimizer rather than a coarse grid — a genuine, documented
  limitation of this mean-field formulation.
- **Chubu Electric/Toshiba field-sensitivity pair**: added to
  `data/amr_experimental_benchmarks.csv` as a secondary-source entry
  (same caveat as the existing Okamura row), plus a new
  `run_field_sensitivity_check()` (the field-axis analog of the existing
  span-axis `run_curve_validation()`) to actually exercise it. Honest
  result: the 4T anchor point doesn't calibrate at all — this device's
  reported 26K span at 4.856kg Gd/0.167Hz exceeds what any fluid flow
  rate in this single-stage 0-D model can achieve — the same kind of
  finding already on record for `Risoe_DTU_Gd_2011`.
- **Flagged, not built**: a second, independent Riso Lab data point and
  the Institute of Tech./Chubu near-zero-span extreme (both from the same
  comparative prototype table) were identified but not added — neither
  was on the paper-mining pass's own "updated priority list," and the
  Riso point's primary source isn't in this repo's `Papers/`.

<a id="phase-10-findings"></a>
## Phase 10 findings (paper-mining pass)

A cross-reference of this repo's current state against the papers in
`Papers/` (see `Paper_Mining_Recommendations.md`) surfaced two concrete,
numerically-anchored additions and two documentation-only flags. Full
details in `ROADMAP.md`'s Phase 10 section.

- **Blow-fraction asymmetry** (Masche et al. 2022): the AMR cycle model
  previously had no notion of flow-waveform asymmetry at all — Qc and
  second-law efficiency implicitly assumed a symmetric 50/50 blow split.
  A calibrated `blow_fraction` parameter now exists in `amr_cycle.py`
  (default 0.5 reproduces every prior result exactly) and as a 6th NSGA-III
  design variable in `optimize.py`. The optimizer independently converges
  toward blow_fraction≈0.37-0.43 across the Pareto front — close to the
  source paper's own 0.416 "best found" value, though not itself
  independent validation since the same calibration data feeds both.
- **(Mn,Fe)2(P,Si) material family** (Hanggai et al. 2026): a third
  pluggable Curie-graded family (`MNFEPSI_FAMILY`), calibrated to the
  paper's cross-validated peak |ΔS_max|~17.6 J/(kg K) at 2T. Notable
  because its Tc window (295.3-331.2K, directly measured across five real
  compositions) sits almost entirely at or above the ASHRAE data-center
  range — the opposite tension from Gd5(SixGe1-x)4, whose giant-MCE ceiling
  sits just below it.
- **Flagged, not built**: the Tušek 2013 comprehensive-comparison paper is
  the only corpus candidate that could validate `geometry_analysis.py`'s
  parallel-plate effectiveness model against a real device — but the exact
  numbers are in the same undigitized Figs. 10-11 already open since
  Phase 7.

<a id="phase-7-findings"></a>
## Phase 7 findings (so far)

**-1. Lifetime operating cost added — HX/pump/motor/controls CAPEX still open.**
Searched for a published $ breakdown of AMR heat-exchanger/pump/motor-drive/
controls capital cost; found none beyond the two materials-cost studies
already cited, so that gap is not fabricated shut here. What *is* new: a
second paper by the same group, Bjørk, Bahl & Nielsen, "The lifetime cost
of a magnetic refrigerator" (Int. J. Refrig. 63 (2016) 48-62), prices
device *operating* cost (electricity at $0.10/kWh) over a stated device
lifetime — a real cost component `economics.py` previously had nothing
for. New `economics.lifetime_cost()` combines this with the existing
`material_cost()` floor and reports both pieces separately (e.g. for the
Astronautics-scale benchmark device at 15-year lifetime: $293 materials
floor vs. $7,338 lifetime electricity — operating cost dominates at this
scale, consistent with the source paper's own finding that operating cost
rivals or exceeds magnet cost). The 2016 paper explicitly states, like its
2011 predecessor, that manufacturing, transportation, maintenance and
auxiliary systems are excluded — so HX/pump/motor/controls/enclosure
hardware CAPEX remains genuinely unquantified. (`core/economics.py`.)

**0. A small, honest curve-level validation check (not a digitization).**
`validation_system.py` previously validated only a single (span, Qc) point
per device. The roadmap's next step was to digitize full published
characteristic curves (Tušek 2010, Nielsen 2011). *Correction:* both source
papers are actually present in this repository (`Papers/AMR Theory and
Modeling/Development of a rotary magnetic refrigerator.pdf` for Tušek 2010;
`Papers/AMR systems and prototypes/Review on numerical modeling...pdf` for
Nielsen 2011) — an earlier pass mistakenly reported them as missing. What's
still true, though: neither one actually contains the multi-point
experimental characteristic-curve data this item needs (Tušek 2010 is a
device-construction paper; Nielsen 2011 is a numerical-modeling review), so
no digitization was fabricated from them. The genuine curve source turned
out to be a third paper, Tušek et al. (2013) (see Phase 7/ROADMAP.md) — also
in `Papers/`, but not yet digitized either.
Instead, 3 of the 5 benchmark devices (Astronautics, DTU, Risø/DTU) already
have a second, companion data point at a different span for the same
physical device (a zero-span max-capacity reading or a max-span
zero-capacity reading); these are now grouped via a new `device_group`
column in `data/amr_experimental_benchmarks.csv`. `run_curve_validation()`
calibrates fluid mdot at the normal operating point exactly as before, then
uses that *same* calibrated system to predict Qc at the companion span —
a genuinely independent check of the model's Qc(span) curve *shape*, since
the companion point plays no part in calibration. Result: Astronautics'
predicted zero-span Qc is 3577.5W vs. 3042W reported (+17.6%); DTU's
predicted Qc at its reported no-load span is 0.0W, matching the reported
0W exactly; Risø/DTU still fails to calibrate at all, consistent with the
Phase 2/6 finding that this device doesn't fit the model. Tušek and Okamura
remain single-point (no companion data available) and are unaffected.
(`core/validation_system.py`, `results/curve_validation.csv`.)

**1. Better solver, real improvement, still not a fix.** Phase 6 flagged
that the 4-point EXTENDED loss-model fit gave negative (unphysical)
k_eddy/base_frac and leave-one-out errors up to +1639%. The negative
coefficients turned out to be partly a *solver* artifact: the fit used
unconstrained least squares and clipped negative results to zero
afterwards, which doesn't minimize anything under the actual physical
constraint (loss coefficients ≥ 0). Re-solving with non-negative least
squares (NNLS, Lawson & Hanson 1974, `scipy.optimize.nnls`) — the
textbook-correct tool for this constraint — removes the negative
coefficients outright and cuts the worst leave-one-out error from
+1639% to +682% (still predicting the smallest device, Tušek, from the
other three). Real progress, but +682% confirms Phase 6's core claim
that a linear model can't span 6.5W-2502W devices; the CORE 3-point fit
stays the production default. (`core/loss_model.py`, NNLS everywhere the
old code did lstsq-then-clip.)

**2. The size/scale-term hypothesis, tested and not confirmed.** Phase 6
speculated a device-size term would fix this. Sorting the four benchmark
devices by Qc and comparing parasitic fraction (`W_parasitic/Qc`) shows
Tušek (smallest, 6.5W) at 11.7%, DTU (102.8W, corrected — see
Paper-Mining Pass Part 6 below) at 25.5%, Okamura (200W) at 36.7%,
Astronautics (largest, 2502W) at 45.3% — a clean *monotonic increase*
with device scale, the opposite direction from a simple fixed-overhead/
economies-of-scale story, which would predict the fraction *falling* as
devices get bigger. (An earlier pass had a fabricated 818W/17.1% figure
for the DTU point that happened to break monotonicity outright; the
qualitative conclusion — that a fixed-overhead term isn't supported —
is unchanged by the correction, since a rising trend doesn't support it
either, but the specific "non-monotonic" claim no longer holds and has
been corrected here.) Astronautics' own source paper attributes its high
fraction to that device's specific "mediocre" electrical-component
efficiency, not a generic size law. This corrects the Phase 6 roadmap
item: a size/scale term isn't supported by the data actually in hand.
The more plausible driver is heterogeneous drivetrain topology/
component-efficiency class across devices (rotary vs. single-bed,
motor/inverter grade) — the concrete next step is more benchmark devices
with independently reported component efficiencies, not a
size-dependent term. (`analyze_parasitic_fraction_scaling()` in
`core/loss_model.py`.)

<a id="phase-6-findings"></a>
## Phase 6 findings

**1. Adding data revealed model fragility — and that's the finding.** Found
a fully-specified 4th benchmark device (Okamura & Hirano 2013: Qc=200W,
COP=2.5 at 5K span, 1.1T, 1kg Gd) via a targeted literature search and added
it to `loss_model.py`'s calibration set, converting the exactly-determined
3-point fit to an over-determined 4-point one — the Phase 3 "zero held-out
slack" flag directly addressed. Result: least-squares over 4 points gives
**negative k_pump and base_frac** (unphysical), and leave-one-out
cross-validation shows the fit does not generalize — errors up to **+1639%**
when predicting the smallest device (Tušek, 6.5W) from the other three,
which span up to 2502W. Four orders of magnitude in scale is apparently too
much to pool into one linear loss model. **This is reported, not hidden**:
the production `StateDependentLossModel()` default still calibrates on the
stable 3-point CORE set; the unstable 4-point result is exposed separately
via `run_extended_diagnostic()` for transparency. The fix (a size/scale term
in the loss model) is a concrete Phase 7 item, not a data-quantity problem —
more devices of similar heterogeneity would likely make this worse, not
better, without a structural model change first.

**2. Economics grounded in an actual costing study.** Replaced the earlier
loosely-sourced $175/kg placeholder with $40/kg (NdFeB magnet) and $20/kg
(Gd) from Bjørk, Bahl & Smith's dedicated magnetic-refrigerator
cost-minimization study (Int. J. Refrig. 34 (2011) 1805-1816), plus their
reported magnet:MCM mass ratios to approximate magnet mass from
`mass_regenerator`. Re-running `optimize.py` with this shifts Pareto-front
costs to $300-5,000 (materials-only floor), consistent with Bjørk et al.'s
own small-device examples ($7-35) scaled up — down from the previous
placeholder's less-grounded ~$500-50,000 range. Still a materials-only
floor, not full system cost (heat exchangers/pumps/motor/controls excluded)
— reconciling that is Phase 7.

<a id="phase-5-findings"></a>
## Phase 5 findings

**1. The giant-MCE material is now honestly modeled — via a different
physics framework, not a patched version of the old one.** `first_order_mce.py`
implements an extended Landau free-energy model (quadratic-quartic-sextic,
the standard tractable route to a first-order transition, in the same
family as Bean-Rodbell) calibrated to the literature's peak entropy change
(≈-18 J/(kg K) at 5T). It reproduces the physically correct signature that
mean-field theory cannot: a genuine first-order jump, with the transition
temperature shifting with field, peaking at T=286.4 K (not exactly at the
nominal Tc=276 K) and giving ΔT_ad≈23 K at 5T — a legitimately "giant" number.

**2. The Curie-matching principle, demonstrated cleanly.** At the ASHRAE
operating point (291 K), Gd5Si2Ge2 collapses to zero cooling capacity — its
peak-effect window (286.4 K / 13.2°C) sits ~5 K below the ASHRAE range.
Moved to its own favorable window (T_cold=281.4 K), it delivers Qc=5319 W,
COP_electrical=7.76 — strong performance, and Gd fails at that same point for
the mirror-image reason. Gd's own peak (296.5 K / 23.4°C) happens to sit
inside the ASHRAE range, which is presumably why it's the standard
room-temperature reference material. **This does not overturn the Phase 1-4
conclusion** (Gd still trails vapor-compression/liquid cooling on COP within
the ASHRAE range) — but it does identify the concrete, literature-supported
next step: the Gd5(SixGe1-x)4 family has composition-tunable ordering
temperature (documented range ~20-276 K, with pure Gd5Si4 at 335 K), so a
composition tuned nearer 291-300 K, if it retains first-order character, is
the genuinely promising untested direction. That's a materials-synthesis
question this simulation suite can motivate but not answer — Phase 6.

**3. A bigger MCE buys capacity, not efficiency.** Even correctly targeted,
Gd5Si2Ge2's COP_electrical (7.76) is close to Gd's own COP_electrical at its
matched point (7.42) despite ~4x the cooling capacity — consistent with
Phase 3's Sobol finding that COP is governed by frequency/flow/field-dependent
losses, not material choice. Material research buys you smaller/cheaper
hardware for a given cooling load, not a fundamentally better COP.

**4. The refrigerant-free case is real but doesn't rescue the emissions
comparison on its own.** `emissions.py` quantifies both refrigerant-leak and
operational emissions. At representative COPs (AMR 5.0 vs. VCC 12.0 vs.
liquid 20.0), AMR's *lower* COP makes its *operational* emissions the
highest of the three, even with zero refrigerant emissions — the
refrigerant-free story is a real, categorical benefit (leak risk, phase-out
liability) but only becomes an emissions *win* once the COP gap identified
in Phase 1-4 closes. Report both numbers plainly in the paper.

<a id="phase-4-findings"></a>
## Phase 4 findings

**1. The mass gap is closed.** `core/thermal.py` computes regenerator
effectiveness from packed-sphere-bed geometry (Wakao & Kaguei 1982 Nusselt
correlation) and a utilization-factor degradation term, so `mass_regenerator`
now genuinely trades off against `eps`/Qc instead of being pure waste. Re-running
`optimize.py` with this active spreads the Pareto front's mass values across
1.0-14.5 kg (vs. pinned at the floor in Phase 3) — the mass/cost tradeoffs on
the front are now trustworthy. **Caveat**: the utilization-degradation term
`(1 - 0.3*U)` is qualitatively motivated, not independently fit — flagged for
Phase 5/6 replacement with a digitized literature curve.

**2. The giant-MCE material can't be fairly evaluated yet — and that's a
finding, not a dead end.** Running Gd5Si2Ge2 through `cascade.py` at the
data-center operating point gives *zero* cooling capacity everywhere. Root
cause: the mean-field/Brillouin framework in `mce_material.py` is built for
continuous (second-order) transitions, valid for Gd — but Gd5Si2Ge2's actual
"giant" effect comes from a first-order magnetostructural transition that
mean-field theory structurally cannot capture, so the model underpredicts
its ΔT_ad by roughly an order of magnitude. This is now flagged explicitly
in `mce_material.py` rather than silently producing a misleadingly-small
number. **Consequence for the paper**: Gd remains the only material in this
suite that's honestly validated, and it happens to be a good fit anyway —
its T_c = 294 K (21°C) sits almost exactly inside the ASHRAE 18-27°C
recommended supply range, which is presumably no accident (it's why Gd is
the standard reference material for *room-temperature* magnetic
refrigeration in the literature). Fairly testing whether a giant-MCE
material changes the Phase 3 conclusion requires a Bean-Rodbell/Landau
model (Phase 5) — not a difference of degree from what's here, a difference
in the physics needed.

<a id="phase-3-findings"></a>
## Phase 3 findings

**1. Sobol resolved.** Replacing the constant `parasitic_fraction` with
state-dependent eddy (~f²H²), pumping (~ṁ²) and base-overhead terms restores
real sensitivity: frequency dominates (ST=0.68, since eddy loss scales with
f²), flow rate next (ST=0.22), field third (ST=0.09) — see
`results/sobol_results.txt` vs. the archived
`results/sobol_results_phase2_constant.txt`. **Caveat**: the 3 loss
coefficients come from an exactly-determined 3-point fit (zero held-out
data) — treat magnitudes as illustrative pending Phase 4's larger benchmark
set.

**2. Optimization exposed a model gap, not hidden it.** Every NSGA-III
Pareto-optimal design sits at the regenerator-mass lower bound, because
`cooling_capacity()` in `amr_cycle.py` still doesn't depend on
`mass_regenerator` at all — more material is pure cost with zero modeled
benefit. This is flagged in `optimize.py`'s own output and is a
prerequisite fix (via the Phase 4 NTU thermal model) before the
mass/cost tradeoffs on the Pareto front can be trusted for the paper.

**3. Cascade staging fixes feasibility, not competitiveness.** Multi-stage
AMR (2/3/4 stages) does recover the span range the single-stage design lost
above 16 K — but at every span in the ASHRAE 5-20 K range, both single- and
multi-stage AMR (Gd, 2 T, this design point) sit well **below** both
vapor-compression and liquid-cooling electrical COP (e.g. at 10 K span:
1-stage 7.2 / 2-stage 4.9 vs. VCC 12.2 / liquid 19.9). **This is the
honest headline result of Phase 1-3**: at Gd/2T and this level of loss
modeling, magnetic cooling does not out-COP conventional data-center
cooling. The paper's contribution is the validated methodology and the
quantified gap + its drivers (frequency-dependent eddy loss, the
mass/effectiveness coupling gap, single-stage span limits) — not a
"magnetic cooling wins" claim. Whether a higher-field design (Gd5Si2Ge2,
already in the materials library) or the GWP/refrigerant-free case
(Phase 5) changes this conclusion is the natural next question.

<a id="phase-2-finding"></a>
## Key Phase 2 finding: ideal vs. electrical COP

Published AMR COP figures are **electrical** (include pump + magnet-motor-drive
power). The model's *ideal* magnetic-cycle-only COP overpredicts these by
118-619%. Adding a calibrated parasitic-loss fraction (default 0.15,
literature range 0.12-0.45 across three benchmark devices) brings the two
comparable lab-scale devices to single-digit-percent agreement:

| Device | Span | COP (lit / ideal / electrical) | Error (electrical) |
|---|---|---|---|
| DTU rotary Gd (2016) | 10.1 K | 4.20 / 14.88 / 4.60 | +9.6% |
| Tušek single-bed Gd (2010) | 15.0 K | 4.60 / 10.02 / 4.00 | -13.0% |
| Astronautics rotary (2014, naval-cooler scale) | 11.0 K | 1.90 / 13.66 / 4.48 | +135.8% (outlier — paper itself cites "mediocre" electrical-component efficiency at that scale) |

**Using the correct (electrical) COP, single-stage AMR at 2 T does NOT beat
vapor-compression across most of the ASHRAE 5-20 K span range** at this
design point (AMR electrical COP ~4-5.5 vs. VCC ~6-24) — a materially
different conclusion than the naive ideal-COP comparison would suggest, and
the central design-motivation for Phase 3 (multi-stage, higher field,
state-dependent loss modeling).

<a id="sobol-sensitivity"></a>
## Sobol sensitivity: a genuine model-structure finding

Sobol analysis (`results/sobol_results.txt`) found electrical COP is ~99.9%
sensitive to `parasitic_fraction` alone — field, frequency, and flow rate
show ~0 sensitivity to COP (though they do affect *cooling capacity*, Qc).
This is because the current model makes `eta_2nd_law` and
`parasitic_fraction` constants rather than state-dependent functions, so
they cancel out of the COP ratio algebraically. This is flagged as a
required Phase 3 upgrade, not silently patched.

<a id="validation-snapshot"></a>
## Validation snapshot (Gd MCE, mean-field model vs. Dan'kov et al. 1998)

| Field | Literature ΔT_ad | Model ΔT_ad | Error |
|---|---|---|---|
| 1 T | 3.2 K | 4.76 K | +48.9% |
| 2 T | 5.8 K | 7.49 K | +29.2% |
| 5 T | 14.6 K | 13.51 K | -7.5% |

Mean-field theory is known to overpredict ΔT_ad near T_c because it
neglects short-range spin correlations / critical fluctuations (de Oliveira
& von Ranke, *Phys. Rep.* 489 (2010) 89-159).

<a id="repo-layout"></a>
## Repo layout

```
core/                              physics, validation, sensitivity, surrogate, and economics modules
tests/                             pytest suite (452 tests as of this session's regenerative-amplification-override work)
data/                              literature-sourced parameter tables + digitized prototype benchmarks
results/                           generated comparison tables, figures, Sobol results, RSM coefficients
main.py                            top-level comparison driver (15 numbered stages + sub-stages, 35 figures)
diagnose_cop_gap.py                standalone diagnostic script
ROADMAP.md                         phased plan/changelog — the authoritative, detailed record of every phase
phase_plan.md                      forward-looking plan for not-yet-started phases/items
Literature_Review.md               source-by-source literature review and citation notes
magcool-dc_technical_walkthrough.md  narrative walkthrough of the model for a non-code audience
hypereg_findings.md                Phase 15 Hypereg literature-findings note (also written to results/)
hysteresis_multiseed_stability.txt  Phase 16-18 follow-up multiseed stability-check output (also written to results/)
```