# magcool-dc

Physics-based simulation suite evaluating **magnetocaloric (magnetic) cooling
for data centers**, benchmarked against vapor-compression CRAC/CRAH and
direct liquid cooling.

## Status: Phase 6 (extended validation stress-test + grounded economics) — see `ROADMAP.md`

## What's implemented

| Module | Purpose |
|---|---|
| `core/mce_material.py` | Mean-field (Brillouin/Weiss) model for continuous-transition materials (Gd, La0.7Ca0.3MnO3); Gd validated, Gd5Si2Ge2 flagged as invalid for this framework |
| `core/first_order_mce.py` | Extended (6th-order) Landau model for first-order/giant-MCE materials |
| `core/giant_mce_analysis.py` | Formal Gd vs. giant-MCE comparison → `results/giant_mce_analysis.txt` |
| `core/emissions.py` | Refrigerant-free GWP/emissions comparison |
| `core/amr_cycle.py` | 0-D AMR cycle model: cooling capacity, ideal/electrical COP, optional NTU-derived effectiveness |
| `core/thermal.py` | NTU packed-bed regenerator effectiveness model |
| `core/loss_model.py` | State-dependent eddy/pumping/base loss model — **Phase 6: added a 4th benchmark device and found the extended fit is unstable (negative coefficients, leave-one-out errors up to +1639%); production default stays on the stable 3-point CORE fit, instability documented via `run_extended_diagnostic()`** |
| `core/optimize.py` | NSGA-III multi-objective optimization — **Phase 6: cost objective now uses `economics.material_cost()`, grounded in Bjørk et al. (2011)'s $40/kg magnet + $20/kg MCM figures** |
| `core/cascade.py` | Multi-stage cascade AMR design, Gd and Gd5Si2Ge2 |
| `core/baseline_cooling.py` | Vapor-compression and liquid-cooling COP correlations |
| `core/economics.py` | CAPEX/OPEX comparison, materials-cost floor grounded in Bjørk et al. (2011) |
| `core/validation.py`, `validation_system.py` | Material- and system-level validation against literature/prototypes |
| `core/sensitivity.py`, `rsm.py` | Sobol sensitivity, RSM surrogate |
| `main.py` | Full comparison across the ASHRAE TC9.9 thermal envelope |

## Quick start

```bash
pip install -r requirements.txt
python -m core.validation            # MCE material model vs. literature
python -m core.validation_system      # AMR system model vs. published prototypes
python -m core.first_order_mce         # first-order Landau model calibration check
python -m core.giant_mce_analysis       # Gd vs. giant-MCE, formal comparison
python -m core.emissions                 # refrigerant-free GWP/emissions case
python -m core.loss_model                 # CORE calibration + Phase 6 EXTENDED diagnostic
python -m core.thermal                     # NTU regenerator effectiveness sweeps
python -m core.sensitivity                  # Sobol, Phase 2 vs. Phase 3 modes
python -m core.rsm                           # RSM surrogate for cooling capacity
python -m core.optimize                       # NSGA-III Pareto front, grounded cost model
python -m core.cascade                         # multi-stage cascade, Gd vs. Gd5Si2Ge2
python main.py                                   # single-stage AMR vs. VCC vs. liquid cooling
```

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

## Sobol sensitivity: a genuine model-structure finding

Sobol analysis (`results/sobol_results.txt`) found electrical COP is ~99.9%
sensitive to `parasitic_fraction` alone — field, frequency, and flow rate
show ~0 sensitivity to COP (though they do affect *cooling capacity*, Qc).
This is because the current model makes `eta_2nd_law` and
`parasitic_fraction` constants rather than state-dependent functions, so
they cancel out of the COP ratio algebraically. This is flagged as a
required Phase 3 upgrade, not silently patched.

## Validation snapshot (Gd MCE, mean-field model vs. Dan'kov et al. 1998)

| Field | Literature ΔT_ad | Model ΔT_ad | Error |
|---|---|---|---|
| 1 T | 3.2 K | 4.43 K | +38.4% |
| 2 T | 6.3 K | 7.16 K | +13.6% |
| 5 T | 14.6 K | 13.18 K | -9.7% |

Mean-field theory is known to overpredict ΔT_ad near T_c because it
neglects short-range spin correlations / critical fluctuations (de Oliveira
& von Ranke, *Phys. Rep.* 489 (2010) 89-159).

## Repo layout

```
core/            physics, validation, sensitivity, surrogate, and economics modules
data/            literature-sourced parameter tables + digitized prototype benchmarks
results/         generated comparison tables, Sobol results, RSM coefficients
main.py          top-level comparison driver
ROADMAP.md       phased plan to reach -suite parity
LITERATURE_REVIEW.md
NOMENCLATURE.md
```
