# COMSOL 2-D/3-D AMR Regenerator-Bed Setup Guide

## Status

**This is a setup specification, not a validated result.** No COMSOL
license was available in the environment this repo was developed in, so
nothing in this document has actually been built, meshed, or solved in
COMSOL — every geometry, boundary condition, and material property below
is transcribed directly from this repo's own 0-D/1-D correlations
(`core/thermal.py`, `core/mce_material.py`, `core/first_order_mce.py`,
`core/geometry_analysis.py`) so that a 2-D/3-D model, if built, is at
least *consistent* with the physics this repo already validates — it is
not independent confirmation of anything. Treat every number produced by
following this guide as unvalidated until it passes the degeneracy check
in [§8](#8-degeneracy-check-do-this-before-trusting-anything-else).

## 0. Why this exists

`core/amr_cycle.py` and `core/thermal.py` are 0-D/lumped models: the
regenerator bed is treated as a single effectiveness number (`eps`) fed
by an NTU correlation, not as a spatial temperature field. That is
adequate for the system-level COP comparisons this repo is built around
(see `amr_cycle.py`'s module docstring: "good enough for system-level COP
comparison... NOT a replacement for a full 2-D/3-D COMSOL regenerator-bed
solve"), but it cannot show axial temperature gradients, thermal
dispersion, or demagnetization-front propagation within the bed — things
that matter if the next step is optimizing internal bed geometry (channel
shape, layer transitions in a Curie-graded bed, wall conduction losses)
rather than just system-level mass/frequency/field trade-offs.

This guide specifies **one AMR unit cell** (one packed-bed or
parallel-plate channel, periodic in the cycle) at the same operating
point this repo already reports numbers for, so the 2-D/3-D result can be
checked against a real, calibrated benchmark before it's trusted for
anything new.

## 1. Scope

- **In scope:** a single-material (Gd), single-stage regenerator bed,
  packed-sphere OR parallel-plate geometry, one full AMR cycle
  (magnetize → cold-to-hot flow → demagnetize → hot-to-cold flow),
  reproducing the bed-internal temperature field and the resulting Qc/COP
  at a literature-calibrated operating point.
- **Out of scope for this guide:** the 6-layer Curie-graded La(Fe,Si)13Hy
  bed (`core/cascade.py`'s `run_graded_cascade`/`validate_astronautics_
  graded_bed`), multi-stage cascades, magnet-circuit field mapping (this
  guide assumes a prescribed, spatially uniform `mu0H(t)` waveform, not a
  solved magnetic field), and structural/mechanical analysis of the
  regenerator housing. Each of those is a genuinely separate model — do
  not silently extend this geometry to cover them without saying so.

## 2. Governing physics

### 2.1 Fluid domain — porous-media flow

Use COMSOL's **Porous Media and Subsurface Flow** module (Brinkman
Equations interface, or Free and Porous Media Flow if you need explicit
resolution of the void space between spheres in a small representative
volume). For a bed-averaged 1-unit-cell model, Brinkman equations with a
Darcy-Forchheimer drag term calibrated to the same friction-factor
correlations `core/thermal.py` already uses are the right level of
fidelity:

- **Packed bed** (Tušek, Kitanovski, Poredoš, *Int. J. Refrig.* 36 (2013)
  1456-1464, Eq. 5): `f = 23.462 · Re^-0.6716`, valid `10 < Re < 5e5`,
  with `Re = ρ_f u_s d_p / μ_f` (superficial velocity, particle
  diameter). This is exactly `core/thermal.py`'s
  `pressure_drop_packed_bed()` — reproduce the same `f(Re)` as a User
  Defined Function or Analytic function in COMSOL rather than falling
  back to a generic Ergun-equation default, so the 2-D pressure drop is
  directly comparable to the 0-D number.
- **Parallel plate** (same paper, Eq. 6): laminar `f = 24/Re`, valid
  `Re < 2300`, with hydraulic diameter `d_h = 2·(plate spacing)`.
- Hydraulic diameter for the packed bed: `d_h = 4·V_bed·ε/A_total` (Eq.
  7 of the same paper), equivalent to the more common
  `d_h = (2/3)·d_p·ε/(1-ε)`.

### 2.2 Solid domain — MCE heat source

The regenerator solid is NOT a passive conductor: during magnetization
and demagnetization it releases/absorbs heat via the magnetocaloric
effect. Model this as a **volumetric heat source term** in COMSOL's Heat
Transfer in Solids (or Porous Media) interface, active only during the
magnetize/demagnetize steps (see §5), computed from this repo's own
`ΔT_ad(T, H)`:

```
q_MCE(T) = ρ_solid · C_total(T) · ΔT_ad(T, μ0H_max) / t_mag
```

where `t_mag` is the magnetization ramp duration (from the AMR cycle
frequency and blow fraction — see §5) and `ΔT_ad(T, μ0H_max)` and
`C_total(T)` are exactly the functions `core/mce_material.py`'s
`MagnetocaloricMaterial.delta_T_adiabatic()` / `.total_heat_capacity()`
compute for Gd (mean-field/Brillouin model), or
`core/first_order_mce.py`'s Landau-model equivalent for a giant-MCE
material. **Do not re-derive this from scratch in COMSOL** — export a
lookup table instead (see §2.3), so the 2-D model's source term is
provably the same physics as the 0-D model it's being checked against.

### 2.3 Building the source-term lookup table

Run this repo's own material model to generate a `(T, ΔT_ad)` table and
import it into COMSOL as an **Interpolation function** (Global
Definitions → Functions → Interpolation):

```python
import numpy as np
from core.mce_material import GADOLINIUM

mu0H = 1.13  # T -- matches the DTU_Eriksen_rotary_Gd_2015 benchmark, see §8
Ts = np.linspace(270, 310, 401)
dTad = GADOLINIUM.delta_T_adiabatic(Ts, mu0H / (4 * np.pi * 1e-7))
C = GADOLINIUM.total_heat_capacity(Ts, mu0H / (4 * np.pi * 1e-7))
np.savetxt("gd_dTad_vs_T_1p13T.csv",
           np.column_stack([Ts, dTad, C]), delimiter=",",
           header="T_K,dTad_K,C_total_J_per_kgK", comments="")
```

Import the resulting CSV directly (File → Import in the Interpolation
function node) rather than hand-fitting a polynomial — this keeps the
COMSOL source term pinned to the exact numbers `core/mce_material.py`
produces, including its mean-field-model overprediction near Tc that
`core/validation.py` already documents (model overpredicts Dan'kov et
al. 1998 by +29-49% at 1-2T, improves to -7.5% at 5T — see
`results/pipeline.log` step 1). **A 2-D model built on this table
inherits that same known bias; it does not fix it.**

## 3. Geometry

### 3.1 Packed-bed unit cell (recommended starting point)

A 2-D axisymmetric or 3-D representative-volume model of a cylindrical
bed segment:

| Quantity | Value | Source |
|---|---|---|
| Particle diameter `d_p` | 0.5 mm default; sweep 0.05-2 mm per `core/geometry_analysis.py` | Phase 7 geometry sweep found the COP-optimal packed-bed diameter at 0.5mm for this repo's 291K/10K-span/2kg/1Hz/1.5T operating point (see `results/pipeline.log` step 3c) — NOT necessarily optimal at the different (1.13T/10.2K-span/1.7kg/0.75Hz) benchmark point used for the degeneracy check in §8; re-sweep if the geometry itself is the thing being studied. |
| Porosity `ε` | 0.365 | `core/thermal.py`'s `regenerator_effectiveness()` default |
| Bed cross-section | 0.002 m² (~5×4 cm face) | `core/thermal.py` default, representative of lab-scale devices in `data/amr_experimental_benchmarks.csv` |
| Solid density `ρ_Gd` | 7900 kg/m³ | `core/thermal.py`'s `RHO_GD` |

### 3.2 Parallel-plate alternative

| Quantity | Value | Source |
|---|---|---|
| Plate spacing (fluid gap) | 0.1-0.2 mm typical; Phase-7 COP-optimum found at 0.1mm at the 291K/10K operating point | `core/geometry_analysis.py` |
| Plate thickness | 0.5 mm default | `core/thermal.py`'s `regenerator_effectiveness_parallel_plate()` |
| Hydraulic diameter | `d_h = 2 × spacing` | Eq. 7 slot limit |

Build both if internal-geometry trade-offs are the actual research
question — `core/geometry_analysis.py`'s Step 3/4 sweeps (see
`results/pipeline.log`) already show a genuine COP-vs-diameter/spacing
interior optimum in the 0-D model; a 2-D model is the natural next check
on whether that optimum survives once axial/radial gradients are
resolved.

## 4. Material properties

| Property | Value | Source |
|---|---|---|
| Gd solid density | 7900 kg/m³ | `core/thermal.py::RHO_GD` |
| Gd solid specific heat (lattice, near-room-T representative) | 236 J/(kg·K) | `core/thermal.py::CP_SOLID_GD` — see that constant's comment: Dan'kov et al. (1998) report a sharp λ-anomaly peak near 300 J/(kg·K) AT Tc, not a flat value; 236 J/(kg·K) is representative of the broader range away from the peak. **If the model geometry spans temperatures near Tc=294K, use `core/mce_material.py`'s full `total_heat_capacity(T)` (lattice + magnetic λ-anomaly) instead of this flat constant**, imported as a second interpolation table alongside §2.3's ΔT_ad table. |
| Gd Curie temperature Tc | 294.0 K | `core/mce_material.py::GADOLINIUM` |
| Gd Debye temperature θ_D | 169.0 K | same |
| Gd molar mass | 157.25 g/mol | same |
| Water density | 997 kg/m³ | `core/thermal.py::water_properties()` |
| Water specific heat | 4186 J/(kg·K) | same |
| Water dynamic viscosity | 8.9×10⁻⁴ Pa·s | same |
| Water thermal conductivity | 0.606 W/(m·K) | same |

These water properties are constant-value simplifications (the module
docstring for `water_properties()` says so directly: "adequate for this
0-D estimate; a full model would use IAPWS correlations"). A 2-D/3-D
model with real temperature gradients across the bed is exactly the case
where that simplification stops being adequate — replace with
COMSOL's built-in IAPWS-97 water material rather than reusing this
constant table, and note in any writeup that this is a genuine
improvement over the 0-D model, not a discrepancy to explain away.

## 5. AMR cycle timing (COMSOL Events interface)

Four-step cycle, matching `core/amr_cycle.py`'s model structure:

1. **Adiabatic magnetization** (`H: 0 → H_max`): MCE source term
   (§2.2) active, fluid flow OFF.
2. **Cold-to-hot flow**: MCE source term OFF, fluid flow ON in the
   cold→hot direction, duration set by the cycle's **blow fraction**
   (`core/amr_cycle.py::BLOW_FRACTION_MASCHE` — default 0.5, i.e.
   symmetric; if reproducing a specific device, check whether its own
   blow fraction is documented — most benchmark devices in
   `data/amr_experimental_benchmarks.csv` implicitly assume 0.5).
3. **Adiabatic demagnetization** (`H: H_max → 0`): MCE source term
   active with the sign flipped, fluid flow OFF.
4. **Hot-to-cold flow**: MCE source term OFF, fluid flow ON in the
   hot→cold direction.

Total period = `1/frequency`. Use COMSOL's Events interface to switch
between these four sub-models within one time-dependent study, or run
four separate quasi-steady studies chained by their output temperature
field if a full transient solve is too expensive. Either way, run
**multiple full cycles** and check the bed's periodic-steady-state
temperature profile has converged (cycle-to-cycle ΔT at any point below
some tolerance, e.g. 0.01K) before reading off Qc — a single-cycle result
from a bed starting at uniform T is not the AMR's actual operating point.

## 6. Boundary conditions

- **Cold end**: fixed-temperature or fixed-heat-flux reservoir at
  `T_cold`, active only during hot-to-cold flow (step 4) when the fluid
  is drawing heat from this reservoir into the bed.
- **Hot end**: fixed-temperature reservoir at `T_hot = T_cold + span`,
  active only during cold-to-hot flow (step 2).
- **Bed walls**: adiabatic (no radial heat loss) as a first pass,
  matching this repo's 0-D model, which has no wall-loss term at all.
  Adding a wall-conduction loss term here would be a genuine addition
  beyond what `amr_cycle.py`/`loss_model.py` currently capture — if
  added, say so explicitly rather than silently changing what's being
  compared.
- **Fluid inlet**: mass flow rate `mdot` (see §7 for the calibration
  value), OFF during magnetize/demagnetize steps.

## 7. Operating point for the degeneracy check

Use **DTU_Eriksen_rotary_Gd_2015** — the one benchmark device in
`data/amr_experimental_benchmarks.csv` that calibrates cleanly under this
repo's own 0-D model (see `core/loss_model.py::CALIBRATION_POINTS_CORE`
and `results/pipeline.log` step 2):

| Parameter | Value |
|---|---|
| Material | Gd (single-Tc approximation of the real Curie-graded 11-layer bed — same simplification `core/validation_system.py` uses for this row) |
| μ0H_max | 1.13 T |
| Regenerator mass | 1.7 kg |
| Frequency | 0.75 Hz |
| T_cold | 289 K (per `loss_model.py`'s calibration comment) |
| Span | 10.2 K |
| Calibrated mdot | 0.084666 kg/s (reproduces Qc=102.8W exactly under the 0-D model) |

Source: Eriksen, Engelbrecht, Bahl, Bjørk, Nielsen, Insinga, Pryds,
"Design and experimental tests of a rotary active magnetic regenerator
prototype," *Int. J. Refrigeration* (2015),
doi:10.1016/j.ijrefrig.2015.05.004. Reported result: Qc=102.8 W at
span=10.2 K, COP=3.1 ("the COP of 3.1 is 11.3% of the Carnot
efficiency"). This repo's own 0-D model reproduces this to -2.1% error
(`results/pipeline.log`, step 2: `DTU_Eriksen_rotary_Gd_2015 ...
err=-2.1% implied_parasitic=0.255`).

## 8. Degeneracy check (do this before trusting anything else)

**Before drawing any conclusion from the 2-D/3-D model that isn't already
in this repo, reproduce this one number:**

> Qc = 102.8 W at span = 10.2 K, μ0H=1.13T, mass=1.7kg, f=0.75Hz,
> mdot=0.084666 kg/s.

If the COMSOL model's Qc at periodic steady state is not within roughly
the same error band the 0-D model already achieves (-2.1%, or generously
±10-15% to allow for genuine spatial effects the 0-D model can't
capture), **do not trust any new geometry/gradient conclusion from the
model until the discrepancy is understood** — it more likely means a
units error, a wrong blow-fraction assumption, or an under-resolved mesh
than a genuine new physical finding. This mirrors the standard this
repo already holds itself to everywhere else (see e.g.
`core/validation.py`, `core/giguere_validation.py`,
`core/validation_system.py`) — a new, more sophisticated model earns
trust by first reproducing what the simpler, already-checked model gets
right, not by producing an interesting-looking number no one has checked.

## 9. What a working 2-D/3-D model would add beyond this repo's 0-D model

- Real axial temperature gradients within the bed (the 0-D model
  represents the whole bed by one `eps` number).
- Whether the linear `span_fraction = max(0, 1 - T_span/(2·dTad_noload))`
  approximation in `amr_cycle.py::cooling_capacity()` — which that
  function's own docstring flags as producing "a sharper, straight-line
  cutoff... than a real AMR device would show" — is actually a
  reasonable approximation, by directly resolving the spatial
  temperature profile near the no-load span limit instead of assuming a
  shape with no literature source (see that docstring's honesty note).
- A genuine check on whether the packed-bed/parallel-plate COP optima
  `core/geometry_analysis.py` finds at a *fixed representative mdot*
  (documented there as a real methodological simplification, since free
  mdot optimization is degenerate in this repo's 2nd-law work model)
  survive when mdot and geometry are optimized jointly with spatially
  resolved thermal-hydraulics.
- For a Curie-graded bed specifically (out of scope for this guide, but
  the natural next extension): whether `cascade.py`'s treatment of each
  graded layer as independently peak-tuned, ignoring inter-layer axial
  conduction, materially changes the graded-cascade Qc/COP numbers in
  `results/graded_cascade_comparison.csv`.

## 10. Known limitations of this guide itself

- Never built or solved — see §0.
- The MCE source-term lookup table (§2.3) inherits the mean-field
  model's own documented ~+30-50% overprediction of ΔT_ad near Tc at low
  field (`core/validation.py`); a 2-D model built on it will reproduce
  that bias faithfully, not correct it.
- Water properties are treated as temperature-independent constants
  unless explicitly replaced per §4's note.
- No wall-conduction or radiative loss term is specified (adiabatic
  walls, §6) — this matches the 0-D model's scope but is itself a
  simplification worth flagging in any writeup.
- Magnetic field is assumed spatially uniform and prescribed as a time
  waveform, not solved from an actual magnet-circuit model.

## References

- Tušek, Kitanovski, Poredoš, "Geometrical optimization of packed-bed and
  parallel-plate active magnetic regenerators," *Int. J. Refrig.* 36
  (2013) 1456-1464.
- Wakao & Kaguei (1982), packed-bed Nusselt correlation.
- Nickolay & Martin (2002), parallel-plate laminar-entry Nusselt
  correlation (as used in `core/thermal.py`'s
  `regenerator_effectiveness_parallel_plate()`).
- Eriksen, Engelbrecht, Bahl, Bjørk, Nielsen, Insinga, Pryds,
  *Int. J. Refrigeration* (2015), doi:10.1016/j.ijrefrig.2015.05.004
  (§7-8 benchmark device).
- Pecharsky & Gschneidner, *J. Magn. Magn. Mater.* 200 (1999) 44-56;
  Dan'kov et al., *Phys. Rev. B* 57, 3478 (1998) (mean-field model
  calibration and its documented limitations, `core/mce_material.py`,
  `core/validation.py`).
- Kitanovski et al., *Magnetocaloric Energy Conversion*, Springer (2015),
  Ch. 2 (ΔT_ad small-signal approximation) and Ch. 6 (2nd-law AMR work
  decomposition, `core/amr_cycle.py::magnetic_work()`).