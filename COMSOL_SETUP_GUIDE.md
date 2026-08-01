# COMSOL 2-D/3-D Regenerator-Bed Setup Guide

**Status: optional, Phase 7 item, not run.** This is a setup guide for
someone with COMSOL access to extend the project beyond its 0-D model — no
COMSOL license was available to actually build or run the model here, so
nothing below is validated against a real COMSOL solve. Every numbered
setting is written to be directly checkable against `core/thermal.py` and
`core/amr_cycle.py` so the person building it can confirm the 2-D model
degenerates to the same effectiveness/Qc numbers this repo's 0-D model
already produces, before trusting anything the 2-D model shows beyond that.

## 1. Why go to 2-D/3-D at all

The 0-D model in this repo (`thermal.py`) treats the regenerator as a single
lumped NTU/effectiveness number:

```
eps = NTU/(NTU + 2) * (1 - 0.3*U)
```

with `NTU` from a Wakao–Kaguei packed-bed correlation and `U` an ad-hoc
utilization-degradation term the README/ROADMAP already flags as
"qualitatively motivated, not independently fit." A 2-D/3-D model replaces
that single closed-form `eps` with a resolved axial (and optionally radial)
temperature profile through the bed over a full magnetization cycle, which:

- replaces the `(1 - 0.3*U)` fudge factor with an actual resolved
  utilization effect,
- resolves axial thermal conduction losses (the dispersion/conduction term
  Nielsen et al. 2011 flags — see `LITERATURE_REVIEW.md` — as a real,
  currently-unmodeled loss mechanism),
- lets you check whether the bed's response time is fast enough at the
  frequencies `loss_model.py` calibrates against (0.155–1.4 Hz across the
  benchmark devices in `data/amr_experimental_benchmarks.csv`).

## 2. Physics setup

Use COMSOL's **Heat Transfer in Porous Media** interface (or, if you want
separate fluid/solid temperatures rather than a local-thermal-equilibrium
assumption, couple **Heat Transfer in Solids** + **Heat Transfer in
Fluids** with a volumetric heat-exchange term — this is the more
defensible choice given AMR beds routinely violate local thermal
equilibrium at high frequency, which is exactly the regime this project's
benchmark devices sit in).

- **Domain**: a 2-D axisymmetric domain is enough to capture axial
  conduction/dispersion and radial wall losses without a full 3-D mesh;
  go 3-D only if you need to resolve genuinely 3-D flow (e.g. a rotary
  wedge-shaped bed rather than a straight packed column).
- **Porous medium properties**: porosity `phi`, particle diameter `d_p`,
  and packed density come straight from `thermal.py`'s existing constants
  — use the same values so the 2-D and 0-D models are geometrically
  identical, which is what makes the degeneracy check in §5 meaningful.
- **Interstitial heat transfer coefficient**: reuse the same Wakao–Kaguei
  correlation (`Nu = 2 + 1.1 Re^0.6 Pr^(1/3)`) as a user-defined
  volumetric heat-transfer coefficient between the two temperature fields,
  rather than switching correlations — otherwise a mismatch between the
  2-D and 0-D `h` values will masquerade as a "geometry effect" that's
  actually just a different correlation.

## 3. The magnetocaloric source term

This is the part that's genuinely new physics, not just re-deriving the
0-D model in more dimensions. The MCE enters as a volumetric heat source
in the solid domain during each half-cycle:

```
Q_MCE = rho_solid * C_p(T) * dT_ad/dt   [W/m^3]
```

or, more directly reusing this repo's material models, as an entropy-based
source using `mce_material.py`'s `entropy_magnetic(T, H)` (Gd,
Brillouin/Weiss) or `first_order_mce.py`'s Landau free energy (Gd5Si2Ge2):

```
Q_MCE = -T * d(S_magnetic)/dt = -T * (dS/dH) * (dH/dt)
```

Two things worth being careful about, both because they're exactly where
Phase 8 found a real bug in the 0-D model:

- **Compute `dS/dH` with the same Newton-solver/small-x-series code path
  `mce_material.py` now uses** (post-Phase-8 fix), not a re-derivation —
  the entropy floor bug that silently zeroed out zero-field entropy near
  Tc (see ROADMAP.md Phase 8) is exactly the kind of thing that's easy to
  reintroduce if this term is rewritten from scratch in COMSOL's
  expression syntax.
- **`dH/dt`** should match the field waveform actually used for
  calibration in `loss_model.py`/`validation_system.py` (discrete on/off
  vs. a ramped field — Nielsen et al. 2011 Fig. 4 shows both conventions;
  the benchmark devices in `amr_experimental_benchmarks.csv` are rotary/
  reciprocating devices, which are closer to a smoothed trapezoidal ramp
  than an instantaneous step).

Implement `Q_MCE` as an **Events**-triggered or time-dependent analytic
expression tied to the same magnetization/demagnetization timing
`amr_cycle.py` assumes, so the cycle-average Qc this model produces is
comparable to the 0-D model's Qc at the same frequency and span.

## 4. Boundary conditions and solver

- **Hot/cold reservoir ends**: fixed-temperature or convective boundary
  conditions at `T_hot`/`T_cold`, matching whatever `ashrae` envelope
  point (`main.py`) you're trying to reproduce.
- **Wall/housing losses**: don't assume adiabatic walls by default —
  Nielsen et al. 2011 §2.11 (see `Papers/AMR systems and prototypes/`)
  explicitly notes most AMR models assume perfect insulation and flags
  this as an under-studied loss mechanism; a lumped convective boundary
  with an estimated housing thermal resistance is a defensible first cut.
- **Time-dependent solver**: run to **cyclic steady state**, not a single
  cycle — periodic AMR models typically need 10-50+ cycles before the
  cold-end temperature profile stops drifting cycle-to-cycle. Use
  COMSOL's built-in periodic/cyclic convergence check rather than a fixed
  cycle count if available.
- **Mesh**: bias toward the hot/cold ends and any wall boundary layer;
  packed-bed NTU correlations already smooth over pore-scale detail, so a
  fine volume-averaged mesh is appropriate — don't try to resolve
  individual particles unless you're deliberately doing pore-scale CFD,
  which is a different (much more expensive) model than this guide
  covers.

## 5. The degeneracy check (do this before trusting anything else)

Before using the 2-D model for anything new, confirm it reproduces the
existing 0-D numbers in the regime where they should agree:

1. Set the 2-D model's porosity, particle diameter, bed mass, frequency
   and field to match one of the calibrated benchmark devices in
   `data/amr_experimental_benchmarks.csv` (Tušek single-bed Gd is the
   smallest/cheapest to mesh).
2. Run to cyclic steady state and extract cycle-averaged Qc and COP.
3. Compare directly against `python -m core.validation_system`'s printed
   Qc(model) for that device.
4. If they disagree by more than the ~10-20% noise band Phase 2/6 already
   found acceptable for lab-scale devices, the discrepancy is telling you
   something about which of `thermal.py`'s simplifying assumptions
   (`(1-0.3U)` term, local thermal equilibrium, adiabatic walls) actually
   matters most for that device — which is the useful output of doing
   this at all, independent of whether the 2-D model is "better."

## 6. What this guide deliberately does not do

- It does not run a COMSOL solve or report any 2-D/3-D numbers — no
  COMSOL license was available in this environment. Anything above is a
  setup specification, not a validated result.
- It does not attempt pore-scale (particle-resolved) CFD — that's a
  different, far more expensive model than "2-D/3-D regenerator bed" as
  scoped in ROADMAP.md.
- It does not extend to a full system model (valves, manifolds, external
  HX) — see the separate, still-open full-system-cost item in
  `ROADMAP.md` for why that's out of scope here too.