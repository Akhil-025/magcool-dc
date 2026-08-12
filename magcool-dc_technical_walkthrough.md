# magcool-dc — Technical Walkthrough (Final / Phase 23)

This traces the model from first-principles magnetism to the final
multi-stage, multi-objective, emissions- and economics-aware comparison —
now extended through alternative working-body architectures, materials
science sensitivity studies, and a static elastocaloric comparison line —
in the order the physics and design logic actually build up. It supersedes
the Phase 6 walkthrough. Every numeric result below was re-verified against
a fresh `python main.py` run (356/356 tests passing); every section states
what's new, what changed, and — carried over from the Phase 6 document's
own convention — what still doesn't work, rather than smoothing that over.

**How to read the honesty flags below:** several late-phase modules (Phases
17–23) were built against two reference books — Kitanovski et al. (2015),
*Magnetocaloric Energy Conversion*, and Tishin & Spichkin (2003), *The
Magnetocaloric Effect and its Applications* — whose specific cited chapters
turned out not to be present in this project's copies (checked directly,
not assumed; see `Literature_Review.md` Section 6 for the full per-chapter
accounting). Where that happened, the affected module is scoped down to a
design-exploration/sensitivity study rather than a validated feature, and
says so in its own docstring and in this walkthrough.

---

## 1. Material physics — `mce_material.py`

### 1.1 Mean-field (Weiss) model

Gd and Gd-based alloys are close to localized-moment ferromagnets, so
magnetization is treated with **mean-field theory**: each spin sees an
effective field made of the applied field *H* plus a molecular field
proportional to bulk magnetization *M*:

```
H_eff = H + λM
λ = 3 kB Tc / [N g² μB² J(J+1) μ0]
```

λ is derived from the known Curie temperature, not fit — computed once per
material.

### 1.2–1.4 Magnetization, entropy, heat capacity

Self-consistent Brillouin-function magnetization, entropy and ΔT_ad from
the Brillouin free energy, and a two-term heat capacity (Debye lattice +
magnetic λ-anomaly).

### 1.5 Gd5Si2Ge2 flagged as structurally invalid here

The mean-field/Brillouin framework is built for **second-order**
(continuous) transitions, valid for Gd. Gd5Si2Ge2's "giant" MCE comes from
a **first-order, magnetostructurally-coupled** transition this framework
cannot capture — running it through `delta_T_adiabatic()` here
*underpredicts* the real effect by roughly an order of magnitude. Gd5Si2Ge2
is kept in this module only as a **materials-library placeholder**; the
credible model lives in `first_order_mce.py` (Section 2).

### 1.6 New: a genuine, documented mean-field limitation — the Curie-point field shift

Dan'kov et al. (1998) also report the Curie-point transition temperature
rising almost linearly with field, ~6 K/T above 2 T. `validation.py`'s
`run_curie_shift_check()` sweeps the model's own peak-ΔT_ad temperature
across 2.0–7.5 T (12 field points, sub-K bounded-Brent precision) and
fits the resulting shift rate. Result:

```
Every field point -> peak at T = 294.5000 K   (no shift at all)
Fitted shift rate: -0.0000 K/T   (Dan'kov et al. 1998 report ~6.0 K/T)
```

The model's peak temperature does **not** move with field — a genuine
mean-field limitation (the molecular-field constant λ is fixed by the
zero-field Tc and never revisited), not a numerical artifact. This is
consistent with de Oliveira & von Ranke (2010)'s review finding that
mean-field theory's blind spot is sharpest right at the critical point.

### 1.7 New (Phase 22 item 1): inhomogeneous / polycrystalline Tc-broadening — `inhomogeneous_broadening.py`

Real polycrystalline Gd is not a single sharp Tc — grain-to-grain
composition/strain variation smears the transition. `BroadenedMagnetocaloricMaterial`
represents this as a Gauss–Hermite-quadrature ensemble over a
Gaussian-distributed Tc (`with_Tc()` on the base material), since Tishin &
Spichkin (2003)'s own Sect. 2.8 ("Inhomogeneous ferromagnets") — the
plan's intended source — is an image-only PDF here (same finding already
flagged for Ch. 11 in Phase 21, Section 8 below), so what's implemented is
the standard literature treatment of Tc-distribution broadening, not book
content.

Sweeping σ_Tc at Dan'kov et al.'s own three fields:

| σ_Tc | Peak ΔT_ad, 1 T | Peak ΔT_ad, 5 T | FWHM, 5 T |
|---|---|---|---|
| 0.0 K (sharp) | 6.464 K | 18.826 K | 31.44 K |
| 1.0 K | 5.619 K | 17.857 K | 35.47 K |
| 3.0 K | 4.781 K | 16.411 K | 40.93 K |
| 5.0 K | 4.177 K | 15.387 K | 45.38 K |

Checked against the pipeline's own material-level validation errors (Step
3 below): broadening genuinely **narrows** the worst-field error against
Dan'kov et al.'s three points, from **+48.9%** (sharp model, σ_Tc=0) to
**+20.9%** at σ_Tc=5.0 K — driven mainly by the 1 T point. But this is a
real trade-off, not a clean win: the 5 T error moves the *other* direction,
−7.5% (sharp) → −14.2% (broadened), consistent with global smoothing of
the transition rather than a field-selective fix. The 1 T error is still
falling monotonically at σ_Tc=5.0 K (the largest value swept), so this
sweep has **not** located an interior optimum — only shown a direction —
and with no real digitized ΔT_ad(T) curve in this repo's corpus to fit
σ_Tc against, no specific σ_Tc value is adopted as calibrated.

---

## 2. First-order materials — `first_order_mce.py`

Replaces the invalid mean-field treatment for first-order/giant-MCE
materials with an **extended Landau free-energy model** (Bean & Rodbell
1962's phenomenological family), expressed as an order-parameter expansion:

```
f(m, τ, h) = (A/2)(τ−1)m² + (B/4)m⁴ + (C/6)m⁶ − hm         (τ = T/Tc, m = M/M_sat)
```

B < 0 produces the discontinuous jump in equilibrium *m* at h=0, τ=1; C>0
keeps the free energy stable. **Equilibrium m(τ,h)** is the real root of
`A(τ−1)m + Bm³ + Cm⁵ = h` that *globally minimizes f* (the reversible
branch through the jump — real hysteresis at the transition is a known,
separately-addressed limitation, see below). **Entropy** via the envelope
theorem: `ΔS_M(τ,h)/(NkB) = −(A/2)·[m(τ,h)² − m(τ,0)²]`.

**Calibration**: (A, B, C) = (10, −4, 8), grid-searched to reproduce the
widely-cited peak |ΔS_M| ≈ 18 J/(kg·K) at 5 T near Tc = 276 K (Pecharsky &
Gschneidner 1997). Checked directly by running the module:

```
1T: dS=-0.25 J/(kg K)   dTad=0.31 K
2T: dS=-0.49 J/(kg K)   dTad=0.61 K
5T: dS=-1.17 J/(kg K)   dTad=1.46 K
Target: dS ~ -18 J/(kg K) at 5T (Pecharsky & Gschneidner 1997 review value)
```

Three pluggable families now exist on this framework: `GD5SI2GE2_FIRST_ORDER`
(fixed composition, Tc=276 K), `GD_FAMILY` (Gd5(SixGe1-x)4(-Ga),
composition-tunable ~20–290 K), `LAFESIH_FAMILY` (La(Fe,Si)13Hy,
composition-tunable ~190–340 K), and `MNFEPSI_FAMILY` (Mn,Fe)2(P,Si),
tunable 295.3–331.2 K, grounded directly in Hanggai et al. 2026's
five-composition Tc-vs-composition data).

### 2.1 Thermal-hysteresis loss (new field, first exercised in Phase 16)

`hysteresis_loss_J_per_kg` is a per-family literature-analog placeholder
(Gd5Si2Ge2: 8.0 J/kg; La(Fe,Si)13Hy: 12.3 J/kg; Mn-Fe-P-Si: 25.0 J/kg;
0.0 for Gd, preserving old behavior everywhere it isn't explicitly turned
on), added to `W_parasitic` in `amr_cycle.py` as
`hysteresis_loss_J_per_kg * mass_regenerator * frequency` — see Section 6.

### 2.2 New (Phase 22 item 2): engineered multi-phase nanocomposite — `nanocomposite_material.py`

`WeightedMaterialEnsemble`/`nanocomposite_tuned_material()` blends three
composition-tuned La(Fe,Si)13Hy phases with triangular weighting, mixing at
the ΔT_ad level (not entropy/heat-capacity, since `dTad_correction` is
already a whole-ratio correction). Motivated by Tishin & Spichkin's Sect.
2.9/Ch. 10 (superparamagnetic/nanocomposite materials) and by the same
Tc-broadening mechanism as Section 1.7 above — both inaccessible in this
project's copy of the book — so this is again the standard literature
treatment of deliberate multi-phase blending, not digitized book content.

**Design-span comparison** (single tuned phase vs. nanocomposite, both
tuned to their own T_mid): the sharply-tuned single phase wins outright at
its own design span (COP=6.726, Qc=4995.8 W at 10 K span vs. the
nanocomposite's COP=4.616, Qc=1753.7 W) — spreading MCE weight across three
phases costs peak performance where the single phase is targeted correctly.

**Off-design robustness follow-up** (both candidates tuned *once* at the
10 K design span's own T_mid=296.15 K, Tc_design=285.40 K, then evaluated
*without retuning* at 5/15/20 K):

| span | nanocomposite COP | nanocomposite Qc | single-phase COP | single-phase Qc |
|---|---|---|---|---|
| 5 K | 3.417 | 863.9 W | 0.000 | 0.0 W |
| 10 K (design) | 4.616 | 1753.7 W | 6.726 | 4995.8 W |
| 15 K | 0.000 | 0.0 W | 0.000 | 0.0 W |
| 20 K | 0.000 | 0.0 W | 0.000 | 0.0 W |

At the one off-design span where either candidate stays feasible (5 K), the
sharply-tuned single phase collapses to Qc=0 (its own no-load ΔT_ad no
longer covers the mismatched span) while the nanocomposite's deliberately
spread working range still delivers positive Qc. **Held to precisely**:
this is a finding about *robustness to narrowing* (avoiding a catastrophic
Qc=0 failure), not about raw performance — and it's a first-pass result at
one spread value and one design/off-design span set, not a general claim.

---

## 3. Material-level validation — `validation.py`, `giguere_validation.py`

Gd's ΔT_ad at 1/2/5 T vs. Dan'kov et al. (1998):

```
mu0H=1.0 T | lit=3.20 K | model=4.76 K | error=+48.9%
mu0H=2.0 T | lit=5.80 K | model=7.49 K | error=+29.2%
mu0H=5.0 T | lit=14.60 K | model=13.51 K | error=-7.5%
```

Systematic high bias at low field, expected from mean-field theory's known
blind spot near Tc. Extended to 7 T against Giguère et al. (1999)'s
independent Gd cross-check: at 5 T the model (13.51 K) sits +22.8% above
the 10.5–11.5 K literature range's midpoint; at 7 T (16.71 K) it sits
+33.7% above the 12.0–13.0 K range's midpoint — **outside both ranges**, a
real disagreement between two published Gd datasets, reported rather than
reconciled.

**Giguère direct-measurement cross-check of the first-order Landau model**
(`giguere_validation.py`), Gd5Si2Ge2 at 7 T: model peak ΔT_ad = 24.17 K
(at T=286.4 K) vs. Giguère et al.'s **direct** measurement of 10.0 K
(9.9 K by an independent Clausius–Clapeyron cross-check) — model
overestimates the direct figure by **+142%** (2.42×). Against Giguère's
own **indirect** (Maxwell-relation) figure of 14.9 K, the model
overestimates by +62% (1.62×) — for reference, Giguère et al.'s *own*
indirect-vs-direct gap is 1.49×. This is not a contradiction of the paper:
it's additive — the paper's own finding is that Maxwell-relation ΔS_M
overestimates ΔT_ad for a first-order transition, and this model's
lattice-only C_lattice denominator (no latent-heat correction, honesty flag
#1 in `first_order_mce.py`) compounds the same direction of error. **This
model should not be refit to match Giguère's direct value** — doing so
would abandon its documented calibration to the peak ΔS_M literature value,
and the 0-D lattice-only-C_p framework can't match both simultaneously.
ΔT_ad predictions from this module are treated as upper-bound-ish, roughly
2–2.5× optimistic, everywhere downstream (cascade capacity/COP numbers
included) via `apply_giguere_correction`.

**Correction-factor field-dependence check**: applying the single-field
(7 T) correction factor at 2 T/5 T *overcorrects* — the raw model's ratio
to Gd (5 T: 1.24) sits closer to Pecharsky & Gschneidner's ~1.30 than the
corrected ratio (0.51), which would wrongly predict Gd5Si2Ge2 underperforms
plain Gd. The correction factor is field-specific, not a universal
constant — documented rather than silently applied everywhere.

---

## 4. From material to machine — `amr_cycle.py`

The AMR 4-step cycle (Barclay 1982) and the core cooling-capacity /
Carnot-work / two-COP equations are the foundation. Every later phase adds
an **optional, backward-compatible plug-in point** to `AMRSystem`, each
`None`/off by default:

```
AMRSystem(..., loss_model=None, use_ntu_thermal_model=False,
          blow_fraction=0.5, particle_diameter=None,
          bed_cross_section_area=None, hypereg_n_parallel=1,
          cycle_type="brayton", thermal_diode=None)
```

- **`loss_model`** (Section 6): state-dependent `W_parasitic`, replacing a
  constant `parasitic_fraction · Qc`.
- **`use_ntu_thermal_model`** (Section 5): regenerator effectiveness from
  bed geometry/NTU instead of a fixed constant.
- **`blow_fraction`** (Masche et al. 2022): asymmetric cold-to-hot vs.
  hot-to-cold flow split, a 7th NSGA-III design variable alongside particle
  diameter.
- **`particle_diameter` / `bed_cross_section_area` / `hypereg_n_parallel`**
  (Section 5): geometry-explicit pumping power, *replacing* (not adding
  to) the generic `k_pump` loss-model term when set.
- **`cycle_type`** (Section 7's Phase 17 sensitivity check): `"brayton"`
  (default, unchanged), `"ericsson"`, `"carnot"` — illustrative,
  qualitatively-ordered multipliers on Qc and η₂, not digitized from
  Kitanovski et al.'s own Sect. 4.1.1–4.1.4 (not in this project's copy).
- **`thermal_diode`** (Section 9): adds `_diode_switching_power_W()` to
  `W_parasitic`, cost-only.
- **Hysteresis power** (Section 2.1): `hysteresis_loss_J_per_kg *
  mass_regenerator * frequency` is added to `W_parasitic` unconditionally
  in both the loss-model and constant-parasitic-fraction paths — 0.0 for
  Gd, so old behavior is preserved unless a first-order material with a
  nonzero placeholder is used.

**Consequence worth being explicit about**: `exergy_eff ≡ η₂(ε)` is still
algebraically degenerate — that identity doesn't change just because ε
can now vary with state or `cycle_type`; it's still not an independent
diagnostic.

---

## 5. Regenerator thermal & geometry model — `thermal.py`, `geometry_analysis.py`, `hypereg_analysis.py`

For a packed sphere bed:

```
V_bed = mass_regenerator / (ρ_Gd · (1−φ))                    porosity φ=0.365 default
a = 6(1−φ)/d_p                        A_total = a · V_bed
Re = ρ_f · u_s · d_p / μ_f             Nu = 2 + 1.1·Re^0.6·Pr^(1/3)     (Wakao & Kaguei 1982)
h = Nu·k_f / d_p
NTU = h·A_total / (ṁ·cp_f)
U = ṁ·cp_f / (2·f·mass_regenerator·cp_solid)                   utilization
ε = [NTU/(NTU+2)] · (1 − 0.3·min(U,1))                          clipped to [0, 0.97]
```

Verified: effectiveness rises from 0.668 at 0.5 kg to 0.970 at 10–15 kg
(mass sweep, f=1 Hz), and from 0.692 at 0.25 Hz to 0.962 at 4 Hz (frequency
sweep, mass=2 kg) — both monotonic and physically sensible.

### 5.1 Geometry-dependent pumping power (Tušek, Kitanovski & Poredoš 2013)

Before this term existed, the pre-existing model had no coupling between
particle/plate geometry and pumping cost at all — effectiveness rose
monotonically as particle diameter shrank toward zero (0.8190 at 2 mm to
0.8936 at 0.001 mm, confirmed directly), which meant no geometry optimum
was even possible in principle. Adding the paper's Ergun/laminar-channel
hydraulic pumping terms (Eqs. 5–7) fixes this. Because free COP-only
optimization over ṁ is itself degenerate in this repo's model (ideal COP
rises monotonically as ṁ→0, driving Qc to zero with it — confirmed by
sweeping ideal COP from 14.149 at 0.5 kg/s to 15.830 at ≤0.005 kg/s), the
geometry sweep below is run at a **fixed representative** ṁ=0.08 kg/s
rather than re-optimizing ṁ per geometry the way the paper's own dynamic
model does — a real, stated methodological difference:

```
Packed-bed sphere diameter:   COP-optimal 0.5 mm (Qc=821.4 W, COP_aug=15.27)
Parallel-plate spacing:       COP-optimal 0.1 mm (Qc=819.4 W, COP_aug=15.27)
```

Both are genuine interior optima — smaller geometry raises effectiveness
only marginally further while pumping power keeps growing. The paper's own
Table 3 optima (0.07/0.17 mm packed-bed; 0.035/0.075 mm parallel-plate) are
not expected to match exactly given the different envelope, flow policy,
and idealized (no pump/motor efficiency) pumping-power accounting used
here — the qualitative confirmation (an interior optimum now exists where
none could before) is the result.

### 5.2 Hypereg parallel-hydraulic pumping power (Klinar et al. 2024)

`pumping_power_packed_bed_hypereg()` splits a series bed into n parallel
sub-regenerators. At the representative operating point (5 kg, mdot=0.08
kg/s, f=1 Hz), splitting into n=16 raises COP_electrical only modestly,
5.264→5.278 (n=1→16) — pumping loss is one of three loss channels, not the
dominant one there. At mdot=0.3 kg/s (~4× baseline), the same split gives
7.034→7.206, a **2.45%** relative gain vs. 0.27% at baseline — consistent
with pumping loss's ṁ²-scaling making it a bigger share of W_parasitic at
higher flow. The qualitative mechanism is confirmed representable in this
repo's own model; no validated optimum n or device-level prediction is
claimed (no benchmark device uses this architecture).

---

## 6. State-dependent parasitic loss model — `loss_model.py`

A constant `parasitic_fraction` makes `COP_electrical` algebraically blind
to field/frequency/flow (Section 9's own Sobol finding). Three
standard-scaling-law terms, **fit, not derived from first principles**:

```
W_eddy = k_eddy · f² · (μ0H)²        (eddy-current loss ~ (dB/dt)²)
W_pump = k_pump · ṁ²                   (Darcy-flow pumping power ~ ṁ²)
W_base = base_frac · Qc                 (controls/inverter/bearing overhead)
```

**CORE calibration** (production default, exactly-determined 3×3 NNLS
solve against Astronautics, DTU_Eriksen_rotary_Gd_2015, and Tušek
single-bed):

```
k_eddy    =  30.519740  W / (Hz^2 * T^2)
k_pump    =   0.000000  W / (kg/s)^2
base_frac =   0.048396  (x Qc)
residual norm = 5.0537 W
```

**Extended (4-point, diagnostic) fit** adding Okamura & Hirano (2013) via
NNLS (non-negative least squares — replacing an earlier unconstrained
least-squares pass whose negative, unphysical coefficients were a known
Phase 6 finding):

```
k_eddy    =  15.087819  W / (Hz^2 * T^2)
base_frac =   0.253013  (x Qc)          residual norm = 12.1329 W
```

Leave-one-out cross-validation on this 4-point set: Astronautics −24.2%,
DTU +193.4%, Okamura −36.5%, and Tušek (smallest device, held out from the
other three which span up to 2502 W) **+472.4%** — an order-of-magnitude
miss. NNLS removes the unphysical negative coefficients and improves the
worst LOO error from an earlier +1639% (unconstrained) to +472.4%, but a
better-behaved solver alone doesn't make one linear model generalize
across four orders of magnitude of device scale.

**Size-effect hypothesis, tested and not adopted**: sorting `W_parasitic /
Qc` by device scale gives Tušek (6.5 W) 0.117, DTU (102.8 W) 0.255, Okamura
(200 W) 0.367, Astronautics (2502 W) 0.453 — **monotonically increasing**
with Qc, the opposite direction from a fixed-overhead/economies-of-scale
story (which predicts the fraction *falling* as devices grow). This is
confirmation of the wrong direction, not of the size-effect hypothesis; a
size/scale term is not adopted. The `CORE` 3-point fit remains the
production default everywhere else in the codebase; the EXTENDED fit is
exposed only for transparency.

---

## 7. System-level validation — `validation_system.py`

Calibration/validation methodology: solve for the one free ṁ that
reproduces reported Qc, then independently check COP_electrical, across a
16-row benchmark set (`data/amr_experimental_benchmarks.csv`, 12–13 device
groups):

| Device | Span | COP (lit / ideal / elec) | Error |
|---|---|---|---|
| Tušek single-bed Gd (2010) | 7.3 K | 5.38 / 20.70 / 1.27 | −76.4% |
| DTU_Eriksen_rotary_Gd (2015) | 10.2 K | 3.10 / 14.73 / 3.03 | −2.1% |
| Lozano/POLO-UFSC r4 (2016) | 6.1 K | 0.58 / 24.64 / 0.66 | +13.5% |
| Lozano/POLO-UFSC r6 (2016) | 5.0 K | 0.65 / 30.06 / 0.63 | −3.5% |
| Lozano/POLO-UFSC r7 (2016) | 3.7 K | 0.76 / 40.62 / 0.85 | +11.6% |
| Lozano/POLO-UFSC r8 (2016) | 3.7 K | 0.83 / 40.62 / 0.91 | +10.2% |
| Okamura & Hirano (2013) | 5.0 K | 2.50 / 30.06 / 3.76 | +50.2% |

Mean |electrical-COP error| across these 7 calibratable rows: **~24%**.
Nine other rows (Astronautics, DTU_Eriksen_MAGGIE_2016, Risø/DTU 2011,
Lozano r1/r2/r3/r5, DTU_MagQueen_2018, ChubuToshiba) return **NO
CALIBRATION FOUND** — their reported Qc is unreachable within a physically
plausible ṁ ∈ [1e-6, 5] kg/s at the device's own field/mass/frequency, or
(ChubuToshiba) at its own anchor field. This is reported directly as a
finding, not hidden: the largest, most data-center-scale devices are
exactly the ones this repo's calibration approach fails on.

**Curve-level check** (Tušek et al. 2013's digitized Figs. 10–11 — see
`Literature_Review.md`'s corrected geometry entry): calibrating AMR(A) at
V*=0.95 against its own anchor point (span=7.26 K, Qc=5.27 W, mdot=0.0074
kg/s) and predicting the intermediate point (span=12.23 K) gives model=
18.01 W vs. lit=2.03 W, **+787.0%** error — a large miss, reported as-is.

### 7.1 Cycle-topology sensitivity (Phase 17) — `run_cycle_type_validation()`

Re-checks rotary-named devices (a naming-convention proxy for cycle
topology, not a literature-confirmed classification per device — none of
the 16 benchmark rows report an explicit AMR cycle-topology classification
in their source papers) as Ericsson-like instead of the default
Brayton-like:

```
DTU_Eriksen_rotary_Gd_2015: COP_err(brayton)=-2.1%  COP_err(ericsson)=+0.6%   [improved]
Astronautics_rotary_2014:   not comparable (did not calibrate under either cycle type)
```

One improvement out of two comparable rotary devices — a directional
sensitivity check, not a validated per-device assignment. Applied to the
6-layer Astronautics graded-bed reproduction (Section 11.2), the same
reclassification does **not** narrow that device's much larger −81.1%
error — the remaining gap there is dominated by other documented issues
(single-Tc-vs-6-real-layers approximation, the ~2.4× ΔT_ad overestimate
from Section 3), not cycle topology.

---

## 8. Baseline technologies — `baseline_cooling.py`

Vapor-compression at `η₂,vcc = 0.42 · Tc/(Th−Tc)`; liquid cooling blended
between free-cooling and mechanical-assist hours.

### 8.1 New (Phase 21): passive/hybrid regenerator augmentation

`augmented_regenerator_cop()` asks whether loading a *conventional* gas
cycle's own regenerator with an MCE material's Curie-point heat-capacity
anomaly can boost that cycle's own COP — motivated by Tishin & Spichkin
(2003) Ch. 11, whose specific content could not be digitized here (confirmed
image-only, no text layer), so this recombines the repo's own existing
heat-capacity and baseline-COP models with an illustrative,
literature-range-anchored effectiveness-to-COP ceiling (capped at +8%),
not a fitted or digitized coefficient.

```
Base VCC COP at Tc=291.15K, Th=301.15K (span=10.0K): 12.2283
Gd (Tc=294.0K, inside window)          eps 0.829 -> 0.860   COP 12.2283 -> 12.2582  (+0.24%)
Gd5Si2Ge2 (Tc=276.0K, outside window)  eps 0.876 -> 0.876   COP unchanged            (+0.00%)
La0.7Ca0.3MnO3 (Tc=267.0K, outside)    eps 0.943 -> 0.943   COP unchanged            (+0.00%)
```

The alignment-vs-mismatch pattern is confirmed directly, not assumed:
candidates whose own Tc falls outside the operating window show ~0% gain
by construction (delta_eps is clipped at 0). Across a 5–20 K span sweep at
fixed T_cold, Gd's gain shrinks from +0.40% (5 K span) to +0.14% (20 K
span) as the window widens and the alignment effect dilutes. Every gain is
capped by construction at 8%, so these numbers read as "this mechanism
could plausibly be worth up to X%, IF the literature-range ceiling holds
for a magnetically-augmented regenerator specifically" — not a validated
device-level COP prediction.

### 8.2 New (Phase 23): elastocaloric reference line

`elastocaloric_reference_cop()` adds a **flat, static literature anchor**
(COP=4.63, the 3.7–5.8 range's midpoint) to the baseline sweep and fig08 —
explicitly *not* a span-dependent simulation like the AMR/VCC/liquid rows.
Sources: Qian et al. (2023), *Science* 380, 722–727 (simulated
steady-state system COP=5.8, up to a 22.5 K span) and Wu et al. (2023),
*Nat. Commun.* 14, 7982 (measured system COP=3.7, at a much narrower ~1 K
span) — neither of this repo's two source books covers elastocalorics at
all, so this comparison row is sourced entirely from these two external
papers, added exactly as scoped: "a comparison row, not a new simulated
device."

---

## 9. Sensitivity analysis — `sensitivity.py`

Sobol/Saltelli (768 samples, N_base=64), two modes:

**Constant-loss mode** (`results/sobol_results_phase2_constant.txt`):

```
parasitic_fraction    ST=1.0052
regen_effectiveness   ST=0.0017
mu0H_max_T, frequency_Hz, fluid_mdot_kgs   ST≈0.0000
```

COP_electrical is algebraically independent of Qc in this mode — both
W_mag and W_parasitic scale linearly with Qc and cancel, so field,
frequency and flow change *how much* cooling you get but not, in this
model, *how efficiently* you get it.

**State-dependent loss mode** (`results/sobol_results.txt`):

```
frequency_Hz          ST=0.8772
fluid_mdot_kgs         ST=0.1346
regen_effectiveness      ST=0.0118
mu0H_max_T                 ST=0.0061
parasitic_fraction            ST=0.0000  (unused in this mode)
```

Frequency and flow now carry **real** sensitivity, because raising them
also raises the eddy-current/pumping loss terms via `loss_model.py`,
producing a genuine efficiency-vs-capacity trade-off. Caveat carried
directly in the output: the loss coefficients themselves come from an
exactly-determined 3-point fit (Section 6), so treat the *magnitude* of
these sensitivities as illustrative, not converged.

---

## 10. Surrogate model — `rsm.py`

Quadratic response-surface for cooling capacity Qc across
[mu0H_max_T, frequency_Hz, fluid_mdot_kgs, regen_eff, span_K]
(300 train / 100 held-out test samples):

```
Held-out R^2 = 0.8620, RMSE = 390.07 W
```

Used to keep future optimizer inner loops fast in principle; `optimize.py`
in fact calls `AMRSystem.run()` directly, since NSGA-III's population
sizes here are cheap enough not to need the surrogate.

---

## 11. Multi-stage cascade design — `cascade.py`

N identical single-stage AMR modules in series, each covering an equal
share of total span; the coldest stage sets the deliverable Qc, and every
downstream stage is re-run at its local span and scaled to match it
(steady-state series constraint): `COP_cascade = Qc_target / Σᵢ W_i`.

### 11.1 Plain-material cascade (Gd, Gd5Si2Ge2)

Single-stage COP collapses toward zero as span approaches the material's
own no-load ΔT_ad ceiling; staging recovers span at a COP cost. The
Gd5Si2Ge2 variant fails at *every* span in the ASHRAE-anchored sweep at all
stage counts — consistent with Section 12's finding that its favorable
window sits below this operating range entirely.

### 11.2 Curie-graded cascade (composition-tuned per stage)

At 10 K span, 3 stages: graded Qc=2388.1 W, COP=2.443 vs. plain-Gd
Qc=1258.0 W, COP=2.415 — grading each stage's composition to track the
local fluid temperature gives more capacity at a comparable COP. Across the
full 64-cell span/stage-count sweep, only 38/64 cells stay fully within the
documented 20–290 K giant-MCE composition range; 25 fall back partially to
plain Gd (reported explicitly per cell, not silently substituted).

### 11.3 Does a 6-layer Curie-graded La(Fe,Si)13Hy bed reproduce Astronautics_rotary_2014?

Six stages tuned to the device's own reported layer temperatures
(294.8–303.6 K), calibrated to reproduce the reported Qc=2502.0 W (mdot=
0.03478 kg/s): predicted COP=0.359 vs. reported COP=1.9, **−81.1% error** —
a large miss, but a genuine improvement in *kind* over Section 3's
single-layer material, which returned a flat "NO CALIBRATION FOUND" for
this same device entirely. The Landau-model Tc-vs-peak-effect offset
(+10.4 K at 2 T, Section 12) is independently confirmed to within ~1 K here
too: the graded bed's own required layer offsets come out +11.1 to +11.5 K
above Jacobs et al. (2014) Table 1's actual reported layer Curie
temperatures. Section 7.1's Ericsson cycle-type reclassification does
**not** narrow this device's error (still −81.1%) — the remaining gap is
dominated by the single-Tc-per-stage approximation of the real 6-layer bed
and the ~2.4× ΔT_ad overestimate (Section 3), not cycle topology.

---

## 12. Giant-MCE analysis — `giant_mce_analysis.py`

Does the giant-MCE material change the COP-competitiveness conclusion
against the ASHRAE 18–27°C (291–300 K) range?

```
Gd peak-effect temperature:        294.5 K (21.4 C) -- INSIDE the ASHRAE range
Gd5Si2Ge2 peak-effect temperature: 286.4 K (13.2 C) -- BELOW the range by ~4.6 K
Landau Tc-vs-peak offset at 2.0T:  +10.4 K above nominal Tc=276.0K (independently
                                     confirmed to ~1K by Section 11.3's graded-bed check)
```

**Test 1** (both at the ASHRAE point, T_cold=291 K, span=10 K): Gd delivers
Qc=1466.1 W, COP_elec=5.12; Gd5Si2Ge2 collapses to Qc=0 (span sits ~9.6 K
from its own peak). **Test 2** (Gd5Si2Ge2 at its own favorable point,
T_cold=281.4 K, span=10 K): Gd5Si2Ge2 delivers Qc=5319.3 W, COP_elec=6.92
(Gd fails at this same point, the mirror-image reason).

**Conclusion, held to precisely**: the giant-MCE effect is real and large
within its own narrow window, but that window is mistargeted for
data-center duty as-is — this does not overturn the Gd-trails-baselines
conclusion. Pecharsky & Gschneidner (1997)'s Gd5(SixGe1-x)4 tunability
(~20–335 K) makes a composition tuned into the ASHRAE range a genuinely
open materials-synthesis question outside what this simulation suite alone
can answer. Also noted: even correctly targeted, Gd5Si2Ge2's COP_electrical
(6.92) sits close to Gd's own (5.12) despite ~4× the Qc — consistent with
Section 9's finding that COP is driven mainly by loss-model terms, not
which material is loaded into the regenerator; a bigger MCE mostly buys
more capacity per kg, not better efficiency.

### 12.1 New: six-way material family comparison — `material_family_comparison.py`

At the representative span (10 K), ranked by 1-stage COP_electrical:

```
1. La(Fe,Si)13Hy (tuned)                        COP=6.72  Qc=4989.1 W
2. Gd (fixed)                                    COP=5.09  Qc=1443.4 W
3. Nanocomposite (LAFESIH 3-phase blend, tuned)  COP=4.62  Qc=1753.7 W
4. Gd5(SixGe1-x)4(-Ga) (tuned)                   COP=4.31  Qc=1352.3 W
   Gd5Si2Ge2 (fixed comp.)                       INFEASIBLE at this span
   (Mn,Fe)2(P,Si) (tuned)                        window doesn't cover this
                                                   point -> falls back to Gd
```

For reference, baselines at this point: VCC COP=12.23, Liquid COP=19.89 —
the best AMR candidate still trails both. This table's value is a fair,
apples-to-apples ranking *among* the giant-MCE options themselves,
including whether each family's documented tunability window can even
reach the ASHRAE point — it does not by itself change the AMR-vs-baseline
conclusion (Section 16).

---

## 13. Multi-objective optimization — `optimize.py`

NSGA-III (Deb & Jain 2014, via pymoo) over field, frequency, ṁ,
regenerator mass, effectiveness, blow fraction, and particle diameter (7
design variables), maximizing COP_electrical and Qc while minimizing cost,
with **material now co-optimized as a categorical choice**: each of Gd,
Gd5(SixGe1-x)4(-Ga) (tuned), and La(Fe,Si)13Hy (tuned) is run through its
own NSGA-III pass and the results merged post-hoc into one globally
non-dominated Pareto front ((Mn,Fe)2(P,Si) is dropped at this operating
point — its required Tc=283.5 K sits outside its own [295.3, 331.2] K
window).

```
Best electrical COP    La(Fe,Si)13Hy  H=1.138T f=0.423Hz mdot=0.4835kg/s
                       mass=1.15kg d_p=0.6309mm bf=0.414 -> COP=9.528, Qc=25480.78W, cost=$175.9
Best cooling capacity  La(Fe,Si)13Hy  H=2.73T  f=2.559Hz mdot=0.4943kg/s
                       mass=10.5kg d_p=0.1332mm bf=0.42  -> COP=6.211, Qc=39282.46W, cost=$3738.8
Lowest cost            (same design as best electrical COP)
Knee point (balanced)  La(Fe,Si)13Hy  H=1.056T f=2.47Hz  mdot=0.4985kg/s
                       mass=6.31kg d_p=0.3416mm bf=0.409 -> COP=8.908, Qc=36287.31W, cost=$900.4
```

21 Pareto-optimal designs found; **100% are La(Fe,Si)13Hy** — Gd and
Gd5(SixGe1-x)4(-Ga) designs are globally dominated at this operating point.
`particle_diameter` spans 0.056–1.835 mm across the front (a real, active
search dimension, not degenerate at either bound); regenerator mass spans
1.15–15.00 kg.

### 13.1 Does thermal-hysteresis loss change the material-selection result? (Phase 16)

A/B rerun with hysteresis loss on vs. forced off (reduced settings,
pop=32/gen=15, for tractability — see the module's own honesty flag on why
this differs from the 40/25 production default):

```
Material                     OFF (pre-Ph16)   ON (Ph16)
Gd                                    1              0
Gd5(SixGe1-x)4(-Ga)                   2              0
La(Fe,Si)13Hy                        21             21
TOTAL front                          24             21
```

La(Fe,Si)13Hy's share rises from 88% (OFF) to 100% (ON) — the *opposite*
direction from the naive expectation (a material with a nonzero hysteresis
placeholder should lose ground when hysteresis is turned on). A follow-up
**full-production-settings, 3-seed stability check**
(`results/hysteresis_multiseed_stability.txt`) subsequently found this
reversal does **not** hold up at full NSGA-III settings — flagged as an
open item, not a settled finding, rather than left uncorrected.

### 13.2 Magnet-geometry (Halbach) cost-term sensitivity (Phase 19)

A/B rerun with the old flat per-Tesla magnet-mass ratio vs. the new
super-linear Halbach-cylinder relation (Section 15):

```
Material                    FLAT   GEOMETRIC
Gd5(SixGe1-x)4(-Ga)            0           1
La(Fe,Si)13Hy                 21          13
TOTAL front                   21          14
mu0H_max_T range:  FLAT 1.06-1.80T (mean 1.31T)  |  GEOMETRIC 1.01-2.88T (mean 1.54T)
```

The expected direction (a nonlinear magnet-mass cost should discourage
very high fields, pulling the GEOMETRIC front's mean field *down*) did
**not** hold in this run — stated plainly rather than assumed, with
plausible causes named (NSGA-III search noise at reduced pop/gen settings,
or the field/COP/Qc trade-off already favoring moderate fields for other
reasons even before the geometric cost term is added).

**Known stale artifact carried from earlier phases**: `optimize.py`'s
module docstring still describes the cost objective in its old form
("~$175/kg Franco et al. 2018 placeholder + flat field² term"); the actual
`cost_index()` implementation calls `economics.bom_cost()`/`bom_cost_geometric()`
(Bjørk et al.-grounded, Section 15) and the two now disagree.

---

## 14. Refrigerant-free emissions comparison — `emissions.py`

At the representative operating point (capacity=1.29 kW, span=10 K, AMR
COP=4.63, VCC COP=12.23, Liquid COP=19.89 — all from Section 16's baseline
sweep):

```
Magnetic (AMR)     refrigerant=0.00 tCO2e/yr  operational=1.21  total=1.21
Vapor-compression  refrigerant=0.04 tCO2e/yr  operational=0.46  total=0.50
Liquid cooling     refrigerant=0.03 tCO2e/yr  operational=0.28  total=0.31
```

**Held to precisely, not oversold**: at this operating point, AMR's lower
COP (the standing Section 16 finding) makes its operational emissions the
highest of the three — the refrigerant-free story is real (zero HFC/HFO
leak-GWP, standard AR5 100-yr GWP values with a 4%/yr leak assumption) but
does not, on its own, flip the total-emissions comparison unless AMR's COP
gap is closed first.

---

## 15. Economics — `economics.py`

### 15.1 Materials-cost floor and full-system estimate (Phase 15)

```
material_cost(mu0H, mass_regenerator):
    magnet_mass = 3.0 * mu0H * mass_regenerator      (flat ratio, Bjork et al. 2011-fit)
    return $40/kg * magnet_mass + $20/kg * mass_regenerator
```

At the representative design point (H=2.0T, 5.0kg Gd): materials BOM =
magnet $1200 (30.00 kg) + MCM $100 + SMM yoke $75 (15.00 kg) = **$1,375**
total. Order-of-magnitude full-system estimate (materials BOM × 10×,
grounded in Russek & Zimm's vapor-compression-AC manufactured-cost
benchmark): **$13,750**. Levelized cost of cooling (CRF-based, Silva et
al. 2017 methodology, 15-yr life, 6% discount): **$0.0341/kWh_cooling**
($0.0126 capital + $0.0216 electricity, materials-only capital basis).

CAPEX/OPEX comparison at the sizing basis (capacity=1.29 kW, span=10 K):

```
Magnetic (AMR)               CAPEX $2,833   OPEX $95/yr
Vapor-compression CRAC/CRAH   CAPEX $451    OPEX $221/yr
Direct liquid cooling         CAPEX $708    OPEX $118/yr
```

AMR CAPEX/OPEX per kW remain explicit pre-commercial placeholders;
`material_cost()` gives a materials-only floor for comparison, not a
device-level bid. **Honesty flag carried directly**: heat exchangers,
pumps, motor/drive, and controls are excluded from the BOM, and Bahl et
al. (2014) note these typically dominate total AMR system cost.

### 15.2 Cost by material family (Phase 15)

```
Gd                                     BOM=$1,375   full-system=$13,750
Gd5(SixGe1-x)4(-Ga) (tuned)            BOM=$1,375   full-system=$13,750
La(Fe,Si)13Hy (tuned)                  BOM=$1,315   full-system=$13,150
```

((Mn,Fe)2(P,Si) dropped at this operating point — same tunability-window
gap as Section 12.1.)

### 15.3 New (Phase 19): geometric Halbach-cylinder magnet mass

`geometric_magnet_mass_kg()`/`bom_cost_geometric()` replace the flat
per-Tesla ratio with a closed-form idealized Halbach-cylinder relation
(Section 2's Bjørk citation correction) — confirmed genuinely super-linear
in field: 1.69× the flat-ratio magnet mass at 1.0 T, rising to **13.98×**
at 3.0 T across `optimize.py`'s own [1.0, 3.0] T search bounds. This is
the missing nonlinearity flagged as absent from the pre-Phase-19 model;
its effect on the NSGA-III material front is Section 13.2 above.

### 15.4 New (Phase 22 item 3): amorphous-material cost/performance note (qualitative only)

Melt-spun/metallic-glass MCM candidates trade *lower* manufacturing cost
(single continuous melt-spinning step, skipping the slow
single-crystal-growth/annealing this repo's priced crystalline families
need) for *lower* peak ΔS_M/ΔT_ad — structural disorder cheap to produce
also broadens and shallows the transition, the same peak-vs-width
trade-off Sections 1.7 and 2.2 already quantify for grain-to-grain
inhomogeneity and deliberate multi-phase blending respectively. No
amorphous-MCM $/kg figure or ΔS_M value is digitized anywhere in this
repo's corpus (Tishin Ch. 9, the natural source, is image-only here), so
this is recorded as a **qualitative note only** — not wired into
`MCM_COST_PER_KG_BY_FAMILY` or any cost/performance calculation, since
adding a numeric placeholder would mean inventing a number this repo has
no basis for.

---

## 16. Baseline comparison sweep & top-level driver — `main.py`, `baseline_cooling.py`

Single-stage AMR vs. VCC vs. liquid cooling vs. Carnot (and, since Phase
23, the static elastocaloric reference) across the ASHRAE 5–20 K span
sweep at Tc=18°C:

| Span (K) | AMR elec. COP | Stages | VCC COP | Liquid COP | Carnot | Elastocaloric ref. |
|---|---|---|---|---|---|---|
| 5 | 5.46 | 1 | 24.46 | 24.78 | 58.23 | 4.63 |
| 10 | 4.63 | 1 | 12.23 | 19.89 | 29.12 | 4.63 |
| 16 | 3.91 | 1 | 7.64 | 18.06 | 18.20 | 4.63 |
| 17 | 2.43 | 2 | 7.19 | 17.88 | 17.13 | 4.63 |
| 20 | 2.33 | 2 | 6.11 | 17.45 | 14.56 | 4.63 |

AMR trails vapor-compression and liquid cooling on electrical COP at every
span in this range, with the gap widening further above 16 K where the
automatic cascade fallback (`core.cascade.staged_baseline_result()`, since
a single stage's own no-load ΔT_ad can no longer cover the span) costs a
further COP step-down. The elastocaloric reference is a flat anchor, not a
simulated curve — its value happens to sit close to AMR's own COP near
10 K span, which is a coincidence of the chosen anchor value rather than a
finding about the two technologies converging.

---

## 17. Figure generation and consolidated design recommendations

`plots.py` generates 34 figures (PNG+PDF) covering every section above —
validation, AMR characteristic curves, cascade/graded staging, sensitivity,
RSM, NSGA-III, economics, emissions, Tc-broadening, nanocomposite
robustness, thermal-diode, fluid-MCE, passive-regenerator, cycle-type,
hysteresis- and magnet-geometry Pareto sensitivity. `design_recommendations.py`
re-organizes the above by demonstrated Sobol sensitivity (not a new
analysis) and prints a starting design point:

```
mu0H_max = 1.056 T,  frequency = 2.47 Hz,  fluid_mdot = 0.4985 kg/s,
mass_regenerator = 6.31 kg,  regen_effectiveness = 0.881
-> COP_electrical = 8.908,  Qc = 36287.31 W,  cost_index = $900.4
   material = La(Fe,Si)13Hy (tuned, Tc=285.3K)
```

Ranked levers: (1) **operating frequency** — the dominant Sobol lever
(ST=0.877); (2) **material/composition choice** — La(Fe,Si)13Hy over plain
Gd is a +32% relative COP gain at the ASHRAE point; (3) **Curie-temperature
grading** — a modest but real Qc/COP gain over uniform-material staging;
(4) **regenerator geometry** — target the interior optima (0.5 mm
packed-bed / 0.1 mm parallel-plate), not minimum particle/channel size; (5)
**field/flow balance** — now co-optimized with material and geometry, so
this lever is best read via the Pareto front rather than in isolation.

---

## 18. Alternative working-body architectures (design-exploration only)

Two directions were explored without any benchmark device in this repo's
corpus to validate against — both are stated as design-exploration tools,
not validated features, consistent with the honesty-flag convention used
throughout this document.

### 18.1 Mechanical-contact thermal diode — `thermal_diode.py`, `thermal_diode_analysis.py` (Phase 18)

Motivated by Kitanovski et al. (2015) Ch. 6, whose pages aren't in this
project's copy of the book (Section 6 of `Literature_Review.md`).
`check_frequency_ceiling_claim()` first checked the plan's own premise
directly: `AMRSystem` has **no internal frequency ceiling** for a diode to
relax — the only bound anywhere is `optimize.py`'s unexplained 5.0 Hz
NSGA-III search bound, not tied to a mechanical-switching limit in any
comment or roadmap entry. The illustrative diode actuation-switching-power
cost then reduces COP_electrical by at most **0.03%** across 0.5–8 Hz — a
small cost-only accounting, since no closed-form relation for how
`rectification_ratio` would improve AMR cycle performance was available to
digitize (so no offsetting heat-transfer benefit is modeled). Net effect is
therefore ≤ the no-diode baseline by construction; this is not a claim that
real thermal diodes are a net negative.

### 18.2 Magnetocaloric fluids (ferrofluid/MR suspension) — `fluid_mce_cycle.py`, `fluid_mce_analysis.py` (Phase 20)

Built on standard Krieger–Dougherty viscosity and Darcy–Weisbach pipe-flow
relations, since neither source book's fluids content is present here.
Volume-fraction sweep (T_cold=291 K, mu0H=1.5 T, mdot=0.05 kg/s, each φ
evaluated at its own favorable span):

```
phi=0.10   dTad=0.3705K   span=0.334K   COP_electrical=108.77   <- best
```

A genuine interior COP optimum exists in φ. But at realistic loadings
(comparison at φ=0.20, span=0.682 K): ferrofluid Qc=8.34 W vs. solid AMR
Qc=723.00 W at the same field/flow — the mixture-heat-capacity dilution
model combined with this architecture's lack of a regenerator collapses
usable span to well under a Kelvin, dramatically less than solid AMR's
regenerator-amplified span. Notably, the ferrofluid system's own
COP_electrical (75.9) is *higher* than solid AMR's (6.5) at this shared,
ultra-narrow span (small W_parasitic relative to W_mag at this flow rate,
plus Qc mostly cancelling out of the COP ratio) — but both trail liquid
cooling (86.6) and vapor-compression (179.1) baselines at this span. **The
real headline is the span, not this specific COP comparison**: it holds
only at the ferrofluid system's own tiny achievable span, which solid AMR
could trivially also hit, with far more Qc.

---

## Known inconsistencies carried into this final state (not fixed as shipped)

Found by executing the modules and checking arithmetic by hand, not
flagged by the code itself — kept unfixed deliberately, as a record of
where documentation drifted from a later code change:

1. `optimize.py`'s module docstring still describes the old ~$175/kg
   Franco-et-al.-based cost objective; the actual `cost_index()`
   implementation calls `economics.bom_cost()`/`bom_cost_geometric()`
   (Bjørk et al.-grounded, Section 15) and the two no longer agree
   (Section 13.2).
2. Phase 16's hysteresis-on-vs-off material-selection reversal (Section
   13.1) does not survive a full-production-settings, multiseed rerun —
   an open item, not a corrected number, since the original reduced-setting
   result is real, just not representative of production settings.
3. Phase 19's magnet-geometry Pareto sensitivity (Section 13.2) did not
   reproduce the expected "geometric cost term lowers the front's mean
   field" direction in this specific run — reported as a negative result
   rather than re-run until it matched expectation.

---

## Equation summary (all in one place)

| Quantity | Equation |
|---|---|
| Molecular field constant | λ = 3kBTc / [Ng²μB²J(J+1)μ0] |
| Effective field | Heff = H + λM |
| Magnetization (Gd, 2nd-order) | M = NgμBJ·B_J(x), x = gμBJμ0(H+λM)/kBT |
| Landau free energy (1st-order families) | f(m,τ,h) = (A/2)(τ−1)m² + (B/4)m⁴ + (C/6)m⁶ − hm |
| Isothermal ΔS (2nd-order) | ΔS_M(T,H) = S_M(T,H) − S_M(T,0) |
| Isothermal ΔS (1st-order) | ΔS_M(τ,h) = −(A/2)[m(τ,h)² − m(τ,0)²]·NkB |
| Adiabatic ΔT | ΔT_ad = −T·ΔS_M / C_total |
| Giguère first-order correction | ΔT_ad,corrected ≈ ΔT_ad,raw / 2.42 (7T-fit factor; overcorrects at 2T/5T) |
| Lattice Cp (Debye) | C_lat = 9nR(T/θD)³∫x⁴eˣ/(eˣ−1)²dx |
| Magnetic Cp (2nd-order only) | C_mag = T·dS_M/dT |
| Gaussian Tc-broadening | Tc ~ N(Tc0, σ_Tc²), Gauss–Hermite quadrature over the ensemble |
| Regenerator eff. (NTU model) | ε = [NTU/(NTU+2)]·(1−0.3·min(U,1)) |
| NTU | NTU = h·A_total/(ṁ·cp_f) |
| Utilization | U = ṁ·cp_f / (2f·m_reg·cp_solid) |
| Packed-bed pumping power (Ergun-type) | ΔP·Q̇, per Tušek, Kitanovski & Poredoš (2013) Eqs. 5-7 |
| Hypereg parallel pumping power | W_pump,n = W_pump,1 / f(n) (n parallel sub-regenerators) |
| Cooling capacity | Qc = ε·ṁcp·ΔT_ad,noload·max(0, 1−span/2ΔT_ad,noload) |
| Carnot work | W_carnot = Qc(Th/Tc−1) |
| 2nd-law efficiency | η₂ = 0.35+0.20ε (× CYCLE_TYPE_FACTORS[cycle_type]) |
| Magnetic work | W_mag = W_carnot/η₂ |
| Parasitic power (state-dependent) | W_parasitic = k_eddy·f²H² + k_pump·ṁ² + base_frac·Qc + hysteresis_loss·m_reg·f + diode_switching_power |
| Ideal COP | COP = Qc/W_mag |
| Electrical COP | COP_e = Qc/(W_mag + W_parasitic) |
| Exergy efficiency | ηₑₓ = COP/COP_carnot ≡ η₂ (still degenerate) |
| Carnot COP | COP_carnot = Tc/(Th−Tc) |
| Cascade COP | COP_cascade = Qc_target / Σᵢ W_i |
| Halbach cylinder bore field | B = Br·ln(Ro/Ri) (Halbach 1980) |
| Geometric magnet mass | m_magnet = ρ_magnet · π(Ro²−Ri²) · L (super-linear in required field) |
| Materials cost floor (flat ratio) | cost = $40/kg·(3·μ0H·m_reg) + $20/kg·m_reg |
| Full-system cost estimate | cost_full ≈ 10 × materials BOM (Russek & Zimm-benchmarked order-of-magnitude) |
| Levelized cost of cooling | LCOC = (CRF·CAPEX + annual electricity cost) / annual kWh_cooling |
| Refrigerant emissions | tCO2e = charge_kg·leak_rate·GWP / 1000 |
| Operational emissions | tCO2e = (capacity·hours·load/COP)·CO2_per_kWh / 1000 |
| Ferrofluid suspension viscosity | Krieger–Dougherty relation (phi-dependent) |
| Ferrofluid pumping power | Darcy–Weisbach pipe-flow pressure drop |
| Elastocaloric reference (static) | COP_elasto ∈ [3.7, 5.8], midpoint 4.63 (Wu et al. 2023 / Qian et al. 2023) |