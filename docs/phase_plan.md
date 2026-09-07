Here's a full implementation plan, structured the way your own `ROADMAP.md` phases are (goal → data source → model → integration points → calibration/validation → tests → deliverables), in priority order by value/effort.

---

## Phase 16 — Hysteresis as a quantified loss term (highest priority — may change your Phase 15 headline result)

**Why first:** your Pareto front is 100% La(Fe,Si)₁₃Hy, a first-order material, and hysteresis is currently invisible to every objective function. This is the one change most likely to alter a conclusion you've already published in `README.md`.

**Data sources**
- Tishin & Spichkin Ch.2 §2.5 ("MCE at the first-order transitions") and Ch.7/8 tables give hysteresis widths (ΔT_hys, typically field-dependent) for La-Fe-Si-H and Gd₅Si₂Ge₂-type alloys.
- Hanggai et al. (2026), already in your `Literature_Review.md`, is a melt-spun Mn-Fe-P-Si paper — check whether it reports thermal hysteresis width alongside the ΔS values you already extracted (worth a targeted re-read of that PDF section, not the whole thing).
- Kitanovski §2.1.4 ("near-zero hysteresis of the MCE") frames why this is a selection criterion, not just an efficiency correction — some first-order compositions are excluded outright in practice, which your `first_order_mce.composition_tuned_material()` currently never does.

**Model**
1. Add a `thermal_hysteresis_K` field to `FirstOrderMCEMaterial` (default `0.0` for continuous-transition/mean-field materials like Gd, which physically have none — keeps `MagnetocaloricMaterial`/`GADOLINIUM` untouched).
2. Hysteresis loss per AMR cycle ≈ irreversible entropy generation over the loop width: a standard approximation is `q_hys ≈ (ΔT_hys / ΔT_ad) * |latent-heat-equivalent entropy change| * T`, or more simply model it as an *additional* parasitic power term analogous to your existing eddy/pump/base terms:
   `W_hys = k_hys * frequency * mass_regenerator * f(ΔT_hys, mu0H)`
   with `k_hys` calibrated the same way `k_eddy`/`k_pump`/`base_frac` are — via `calibrate_loss_coefficients()` in `loss_model.py`, but this coefficient won't have your existing 3-point CORE benchmark set to fit against (none of your current `amr_experimental_benchmarks.csv` rows are on a hysteretic material at a documented hysteresis width). So the honest path is: implement it as an *estimate* from the digitized M-H or S-T loop width in the literature (a fixed per-material constant, not a fitted-to-benchmarks coefficient), and flag it exactly the way you flagged the DTAD_CORRECTION_FACTOR field-dependence issue in `giguere_validation.py` — i.e., document that it's a first-pass estimate, not a calibrated fit.
3. Add `hysteresis_power_override` parameter to `StateDependentLossModel.parasitic_power()`, following the exact `pumping_power_override` pattern (default `None` = unchanged behavior, so every existing caller is unaffected). This preserves backward compatibility exactly as Phase 15 did for `pumping_power_override`.

**Integration points**
- `core/first_order_mce.py`: add hysteresis field to `FirstOrderMCEMaterial`, populate for `LAFESIH_FIRST_ORDER`, `GD5SI2GE2` (whichever constant name it has), `MNFEPSI_FIRST_ORDER`.
- `core/loss_model.py`: add the override param to both `StateDependentLossModel.parasitic_power()` and `RotaryDriveLossModel.parasitic_power()` (it calls `super()`, so only the base class needs the actual math).
- `core/amr_cycle.py`: `AMRSystem.__init__` reads `material.thermal_hysteresis_K` if present and passes it through to the loss model call inside `cooling_capacity()`/`magnetic_work()`, mirroring exactly how `particle_diameter` currently threads through to `_geometry_pumping_power_W()`.
- `core/optimize.py`: no new design variable needed — hysteresis loss is a material-intrinsic property, not a tunable knob, so it just changes the `f1`/`f2` objective values `AMRDesignProblem._evaluate()` computes per-material. This is the whole point: it should make first-order candidates *look worse* in the merged Pareto front if the effect is real.
- `core/cascade.py`: same hookup, since cascade beds are explicitly first-order-material-heavy.

**Validation**
- Add `run_hysteresis_sensitivity()` to a new or existing diagnostic module: rerun `core/optimize.run_optimization()` with hysteresis on vs. off, report whether the merged front's material composition changes. This is the actual deliverable — a documented answer to "is the Phase 15 La(Fe,Si)₁₃Hy dominance robust to hysteresis?"
- New `tests/test_loss_model.py` cases: `hysteresis_power_override=None` reproduces old numbers exactly (regression-guard, same pattern as your existing pumping override tests); nonzero override adds linearly.

**Deliverables**: `results/hysteresis_sensitivity.txt`, ROADMAP Phase 16 entry documenting the estimate's provenance and honesty-flagging that it's not benchmark-calibrated.

**Effort**: medium — mostly plumbing given your existing override pattern; the literature-value extraction (finding actual ΔT_hys numbers) is the slow part.

---

## Phase 17 — AMR cycle topology: Ericsson-like and Carnot-like variants

**Data source**: Kitanovski §4.1.1 (Ericsson-like), §4.1.2 (hybrid Brayton–Ericsson), §4.1.3 (Carnot-like), §4.1.4 (max specific cooling power per cycle type) — these are closed-form relations, not requiring digitization, since the chapter derives them analytically from the same H-h-T-s framework your Chapter-1-equivalent thermodynamics already uses conceptually.

**Model**
1. Refactor `AMRSystem` so the cycle-shape assumption is a strategy object rather than baked into `cooling_capacity()`/`magnetic_work()`. Add a `cycle_type` parameter (`"brayton"` default = current behavior, `"ericsson"`, `"carnot"`) to `AMRSystem.__init__`.
2. Each cycle type changes how the magnetization/demagnetization and blow-flow phases relate to the effective temperature-entropy loop area — Kitanovski's §4.1.4 gives max specific cooling power formulas per cycle directly, so this is closer to swapping a closed-form multiplier than rewriting the numerical core.
3. Keep `"brayton"` as the literal default so nothing downstream (existing tests, `optimize.py`, `cascade.py`) changes unless explicitly opted in — same backward-compatibility discipline as `particle_diameter=None`.

**Integration points**
- `core/amr_cycle.py`: new `_cycle_shape_factor(cycle_type, ...)` helper, called from `cooling_capacity()`.
- `core/validation_system.py`: check whether any of your 12 benchmarked devices report which cycle topology they use (several rotary devices are effectively Ericsson-like due to continuous rather than stepped field profiles) — this could *improve* an existing validation fit, not just add a feature.
- `core/optimize.py`: optionally add `cycle_type` as a categorical 8th "design choice" analogous to how material family is handled (run separately per cycle type, merge post-hoc — same "option (b)" pattern Phase 15 already established for materials).

**Validation**: rerun `validation_system.py` per device with cycle_type inferred from its drive mechanism (rotary continuous-field devices → Ericsson-like; reciprocating stepped-field devices → Brayton-like) and see if per-device residuals shrink.

**Tests**: `tests/test_amr_cycle.py` — add cycle_type="ericsson"/"carnot" cases checking the specific-cooling-power ordering matches Kitanovski's own qualitative ranking (Carnot-like should bound the other two from above, since it's the ideal reference).

**Effort**: medium-high — this touches the core cycle math, needs care to not perturb existing Brayton-path numbers (170 tests currently passing).

---

## Phase 18 — Thermal diode module (new capability, not a refinement)

**Data source**: Kitanovski Ch.6 in full — §6.2 (active: thermoelectric, thermionic, spincaloritronic, mechanical-contact diodes), §6.3 (passive: bulk and nanoscale rectifiers), §6.6 (three concrete architectures: single-stage, cascade, AMR each combined with thermal diodes).

**Model** — this is genuinely new physics, so scope it narrowly first:
1. Start with the mechanical-contact active thermal diode (§6.2.4), since it's the mechanism actually used in real prototypes cited elsewhere in your corpus, rather than the more exotic thermoelectric/spincaloritronic options — lower risk of building something with no benchmark to check against.
2. New file `core/thermal_diode.py`, structured like `thermal.py`: a `MechanicalContactDiode` class with `forward_conductance`, `reverse_conductance` (the rectification ratio is the key figure of merit — Kitanovski reports typical ratios in §6.5's review), and a `cycle_time_reduction_factor()` that estimates how much a diode-assisted design could raise achievable AMR frequency vs. valve-switched flow (this is the actual payoff — higher frequency → higher power density, directly comparable to your existing `frequency` design variable).
3. Add an `AMRSystem` variant or flag (`thermal_diode_assisted=False` default) that, when true, relaxes the frequency ceiling implicitly assumed elsewhere (check whether `optimize.py`'s `_XU` upper bound on frequency is currently set by a mechanical-switching limit — if so, that's the exact number a diode-assisted design should be allowed to exceed).

**Integration points**
- `core/loss_model.py`: diode switching isn't free — needs its own small parasitic term, same additive pattern.
- `core/optimize.py`: `thermal_diode_assisted` as a categorical variant, run as its own per-configuration NSGA-III search merged post-hoc (same pattern as material families and, per Phase 17, cycle types) — you're building a consistent "run per categorical option, merge fronts" idiom across the codebase, worth actually factoring into a shared helper at this point rather than copy-pasting a third time.

**Validation**: this is the one area where you may genuinely have no benchmark device in `amr_experimental_benchmarks.csv` — flag honestly (same as you did for the parallel-plate validation gap) rather than forcing a fit. Document as a design-exploration tool, not a validated feature, until/unless a diode-equipped prototype paper surfaces.

**Effort**: high — new physics, no benchmark, real risk of building something un-calibratable. Consider scoping down to just "what frequency ceiling would need to be broken for this to matter" as a sensitivity study before building the full diode model.

---

## Phase 19 — Magnetic field source: field-vs-mass geometry model

**Data source**: Bjørk et al. arXiv:1410.1987 (already cited in your `Literature_Review.md` for the cost ratio, but its actual technical contribution — optimized Halbach field-vs-magnet-mass trade-off curves — hasn't been extracted); Kitanovski §3.4 (2D static/rotary Halbach, 3D Halbach assemblies) and §3.5 (comparative evaluation table of magnet assembly designs) for the geometric coefficients.

**Model**
1. New `core/magnet_geometry.py`: a `halbach_field_vs_mass(mu0H_target, air_gap_volume)` function returning required magnet mass, based on the Halbach cylinder analytical relation (field inside a Halbach cylinder is a closed-form function of the ratio of outer/inner radius and remanence — this is standard and doesn't need digitization, it's textbook permanent-magnet theory that Kitanovski §3.2/§3.4 lays out).
2. Replace `economics.py`'s current flat $/kg-with-a-ratio-fudge-factor with `bom_cost()` calling into this new geometric relation for the magnet-mass term specifically, keeping the MCM cost term untouched.
3. This directly changes `optimize.py`'s `cost_index()` — currently `mu0H` and `mass` (of MCM) are independent design variables with field cost entering only linearly; after this change, achieving high `mu0H` should cost *nonlinearly* more magnet mass for a fixed air-gap geometry, which is physically real and currently absent.

**Integration points**
- `core/economics.py`: `bom_cost()`, `full_system_cost_estimate()`.
- `core/optimize.py`: `cost_index()` — this is a real behavior change to the optimizer, not just a new capability, so treat it like Phase 15's BOM cost upgrade: rerun the full Pareto search and diff the new front against the old one, documenting what shifted.

**Validation**: Bjørk et al.'s own reported "2T is the field/cost sweet spot" claim (already summarized in your lit review) is a checkable prediction — see if your new geometric cost curve reproduces that optimum, or explain why it doesn't.

**Effort**: medium — the physics is closed-form and well-established, main work is careful integration into the existing cost pipeline without breaking the Phase 15 BOM work.

---

## Phase 20 — Magnetocaloric fluids as an alternative working-body class

**Data source**: Kitanovski Ch.5 in full (ferrofluid and magnetorheological suspension rheology, §5.4/§5.5 device design notes); Tishin doesn't cover this (its Ch.11 passive-regenerator focus is solid-state).

**Scope decision needed before building**: this is a different working-body topology (fluid flows and is magnetized in place, vs. solid regenerator bed with separate heat-transfer fluid) — it's closer to a new sibling of `AMRSystem` than a parameter on it. Recommend a new `core/fluid_mce_cycle.py` with its own `FerrofluidMCESystem` class rather than trying to force it through `AMRSystem`'s existing packed-bed/parallel-plate geometry assumptions in `thermal.py`.

**Model**
1. Effective magnetocaloric entropy change of a suspension scales with particle volume fraction — Kitanovski §5.4/§5.5 gives the relevant rheology-adjusted relations (viscosity vs. volume fraction competes against MCE intensity vs. volume fraction, an optimization the book flags but Chapter 5 doesn't fully resolve either — genuinely open territory).
2. Pumping power model differs fundamentally: no packed-bed pressure drop, but pumping a more viscous magnetic fluid still costs power — needs `thermal.py`'s pumping-power correlations replaced with fluid-suspension viscosity correlations from §5.2.

**Why this matters for a data-center pitch specifically**: your `README.md` already frames this whole project against liquid cooling as a comparison baseline in `baseline_cooling.py`. A ferrofluid-based magnetocaloric system is structurally a liquid-cooling system with magnetocaloric fluid instead of water/dielectric — worth at least a first-pass comparison against your existing direct-liquid-cooling baseline, since it's the closest architecture your competitor technology already uses.

**Validation**: Kitanovski §5.4.2 reviews prior refrigeration/heat-pumping experimental work with magnetocaloric fluids — check whether any of those studies report Qc/COP numbers usable the way your `amr_experimental_benchmarks.csv` does for solid AMR.

**Effort**: high — new system class, new heat-transfer correlations, uncertain benchmark availability. Lowest priority of the "new capability" items unless the data-center liquid-cooling comparison angle is specifically valuable to you.

---

## Phase 21 — Passive/hybrid regenerator configurations

**Data source**: Tishin §11.1 (passive magnetic regenerators — rare-earth intermetallics/metals used as passive regenerator material in conventional gas-cycle refrigerators), §11.2.3 (magnetically-augmented gas regenerators), §11.2.4 (hybrid magnetic working bodies).

**Model**
1. This is architecturally simpler than Phase 18/20: a passive regenerator doesn't cycle magnetization actively — its magnetic heat capacity anomaly near Tc is exploited passively to boost a conventional gas-cycle regenerator's effective heat capacity at low temperature. Since your `baseline_cooling.py` already has vapor-compression correlations, this is a genuine hybrid: `baseline_cooling`'s gas cycle + a passive-regenerator heat-capacity boost term.
2. New function in `baseline_cooling.py`: `augmented_regenerator_cop(base_cop, passive_regenerator_material, T_range)` — the augmentation factor is a function of how well the material's heat capacity peak (already modeled in `mce_material.py`'s entropy curves — you have this data already, just haven't used it this way) aligns with the gas cycle's cold-end temperature.

**Integration points**
- `core/baseline_cooling.py`: new augmented-COP function.
- `core/mce_material.py`: reuse existing heat-capacity/entropy data — no new material physics needed, just a new *use* of data you already compute.

**Why this is actually cheap**: unlike Phases 18/20, this doesn't need new physics or new benchmark data — it recombines your existing `mce_material.py` entropy curves with your existing `baseline_cooling.py` gas-cycle correlations in a new way. Good candidate to actually do before the harder items above.

**Effort**: low-medium.

---

## Phase 22 — Material model refinements (magnetoelastic/anisotropic, inhomogeneous broadening, superparamagnetic/nanocomposite)

**Data source**: Tishin §2.8 (inhomogeneous ferromagnets), §2.9/Ch.10 (superparamagnetic/nanocomposite/molecular cluster systems), §2.10 (anisotropic/magnetoelastic contributions), §2.12 (MCE–elastocaloric coupling), Ch.9 (amorphous materials).

**Priority ranking within this phase** (don't do all four at once):
1. **Inhomogeneous/polycrystalline broadening (§2.8)** is the highest-value of the four, because it directly explains a discrepancy you've *already documented and left open* — your `ROADMAP.md`/`giguere_validation.py` note the model's sharp mean-field transition is "narrower than the real, hysteresis/inhomogeneity-broadened transition." A first-pass Gaussian-broadening convolution of your existing mean-field entropy curve (a few lines added to `mce_material.py`, calibrated against the same Dan'kov et al. Gd curve you already validate against) could close a gap you've named but not addressed.
2. **Superparamagnetic/nanocomposite (§2.9, Ch.10)** as a genuinely new *candidate material family* alongside Gd/GD5SI2GE2/LAFESIH/MNFEPSI in `material_family_comparison.py` — lower priority since it needs new composition-tuning data your corpus doesn't currently have digitized.
3. **Magnetoelastic/anisotropic (§2.10)** and **amorphous materials (Ch.9)** — lowest priority; interesting physics but no clear near-term payoff for the data-center application specifically (amorphous materials trade lower cost for lower peak ΔS, worth a one-line cost/performance note in `economics.py` rather than a full model).

**Effort**: item 1 is low-medium and closes a documented gap; items 2–3 are exploratory and lower-value for this application.

---

## Phase 23 — Alternative caloric baseline (elastocaloric comparison point)

**Data source**: Kitanovski Ch.10 §10.3 (elastocaloric effect, materials, device concepts).

**Model**: not a new device to simulate in full — just add elastocaloric as a *comparison row* in the same table/plot your `baseline_cooling.py` already uses for vapor-compression and liquid cooling, using published elastocaloric COP/exergy-efficiency figures as a static reference point (similar treatment to how Carnot COP is already just a reference line in your `plots.py` fig08, not a simulated system).

**Integration points**: `core/baseline_cooling.py` (a `elastocaloric_reference_cop()` lookup, not a simulation), `core/plots.py` fig08 gets a fourth reference line.

**Why low effort, still worth doing**: elastocaloric is magnetocaloric's closest competitor technology (better hysteresis in some NiTi-based systems, no rare-earth dependency) and you're already asking the "is magnetocaloric worth it vs. alternatives" question in your baseline comparisons — this just adds one more honest competitor to that comparison at very low implementation cost.

**Effort**: low.

---

## Suggested build order

1. **Phase 16** (hysteresis) — cheapest high-value item, directly interrogates your existing headline result.
2. **Phase 21** (passive/hybrid regenerators) — cheap, reuses existing data.
3. **Phase 23** (elastocaloric reference line) — cheap, one comparison row.
4. **Phase 22 item 1** (inhomogeneous broadening) — closes a gap you've already named in ROADMAP.md.
5. **Phase 19** (magnet field-vs-mass geometry) — medium effort, closed-form physics, meaningfully improves cost realism.
6. **Phase 17** (Ericsson/Carnot cycle variants) — medium-high, touches core cycle math, needs careful regression testing against your 170 passing tests.
7. **Phase 18** (thermal diodes) — high effort, no current benchmark, treat as exploratory.
8. **Phase 20** (magnetocaloric fluids) — highest effort, new system class, do last or skip unless the liquid-cooling-comparison angle is a priority for you.

Each phase should get its own `ROADMAP.md` entry in your existing style — goal, what was found, what was deliberately not done and why, test count before/after — and a new `tests/test_<module>.py` following your existing per-module test file convention.