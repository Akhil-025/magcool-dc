# magcool-dc — Technical Walkthrough (Final / Phase 6)

This traces the model from first-principles magnetism to the final
multi-stage, multi-objective, emissions- and economics-aware comparison, in
the order the physics and design logic actually build up. It supersedes the
Phase 2 walkthrough — every section below reflects the code as it stands at
Phase 6, including what changed, what got added, and what still doesn't
work.

---

## 1. Material physics — `mce_material.py`

### 1.1 Mean-field (Weiss) model

Unchanged from Phase 2. Gd and Gd-based alloys are close to localized-moment
ferromagnets, so magnetization is treated with **mean-field theory**: each
spin sees an effective field made of the applied field *H* plus a molecular
field proportional to bulk magnetization *M*:

```
H_eff = H + λM
λ = 3 kB Tc / [N g² μB² J(J+1) μ0]
```

λ is derived from the known Curie temperature, not fit — computed once in
`__post_init__` for every material in the library.

### 1.2–1.4 Magnetization, entropy, heat capacity

Also unchanged: self-consistent Brillouin-function magnetization, entropy
and ΔT_ad from the Brillouin free energy, and a two-term heat capacity
(Debye lattice + magnetic λ-anomaly). See the Phase 2 walkthrough for the
full derivations — the equations are identical.

### 1.5 New: Gd5Si2Ge2 flagged as structurally invalid here

Phase 4 added an explicit honesty flag directly in this module: the
mean-field/Brillouin framework is built for **second-order** (continuous)
transitions, valid for Gd. Gd5Si2Ge2's "giant" MCE comes from a
**first-order, magnetostructurally-coupled** transition that this framework
cannot capture — running it through `delta_T_adiabatic()` here
*underpredicts* the real effect by roughly an order of magnitude (~1 K at 2
T here vs. several-to-double-digit K reported experimentally). Gd5Si2Ge2
is kept in this module only as a **materials-library placeholder**; the
credible model for it lives in the new module below.

---

## 2. New: First-order materials — `first_order_mce.py` (Phase 5)

Replaces the invalid mean-field treatment of Gd5Si2Ge2 with an **extended
Landau free-energy model**, the standard tractable route to a genuine
first-order transition (same family as Bean & Rodbell 1962), expressed
directly as a phenomenological order-parameter expansion rather than via
self-consistent lattice strain:

```
f(m, τ, h) = (A/2)(τ−1)m² + (B/4)m⁴ + (C/6)m⁶ − hm         (τ = T/Tc, m = M/M_sat)
```

B < 0 (negative quartic term) is what produces the discontinuous jump in
equilibrium *m* at h=0, τ=1; the positive sextic term C keeps the free
energy stable. Reduced field uses the same natural scale as the Brillouin
argument in `mce_material.py`: `h = gJμB·μ0H / (kB·Tc)`.

**Equilibrium m(τ,h)** is found as the real root of
`A(τ−1)m + Bm³ + Cm⁵ = h` that *globally minimizes f*, not just any
stationary point — this correctly selects the reversible branch through the
jump (real hysteresis at the transition is a known limitation this
simplified treatment doesn't capture).

**Entropy** (via the envelope theorem, since ∂f/∂m = 0 at equilibrium):

```
S(τ,h)/(NkB) = −(A/2)·m(τ,h)²
ΔS_M(τ,h) = −(A/2)·[m(τ,h)² − m(τ,0)²]·NkB
```

**Calibration**: (A, B, C) = (10, −4, 8), grid-searched to reproduce the
widely-cited peak |ΔS_M| ≈ 18 J/(kg·K) at 5 T near Tc = 276 K (Pecharsky &
Gschneidner 1997 / Gschneidner & Pecharsky review). The literature itself
spans roughly 10–18.5 J/(kg·K) depending on sample preparation, and this
calibration deliberately targets the upper end — one defensible choice, not
the only one.

**Honesty flags baked into the module** (both real limitations, not
cosmetic caveats):
1. ΔT_ad here divides by the Debye lattice heat capacity *only* (no
   magnetic λ-anomaly term, unlike Gd's treatment) — appropriate for a
   continuous transition, not fully correct for a first-order one where
   latent-heat structure matters. Treat this ΔT_ad as upper-bound-ish, not
   validated.
2. There is **no independent validation dataset** for this model anywhere
   in the repo (unlike Gd's Dan'kov et al. 1998 check) — it's calibrated to
   one literature number and not cross-checked against a second, held-out
   measurement. A concrete next step (flagged, not done): add Giguere et
   al. (PRL 83, 2262 (1999))'s direct ΔT_ad measurement as that check.

---

## 3. Material-level validation — `validation.py`

Unchanged from Phase 2 — Gd's ΔT_ad at 1/2/5 T vs. Dan'kov et al. (1998),
max error 38.4%, systematic high bias at low field expected from mean-field
theory's known blind spot near Tc.

---

## 4. From material to machine — `amr_cycle.py`

The AMR 4-step cycle itself (Barclay 1982) and the core cooling-capacity /
Carnot-work / two-COP equations are unchanged from Phase 2 — see that
walkthrough for the full derivation. What Phase 3–4 added is **two optional
plug-in points**, both off by default so old behavior (and
`results/comparison_table.csv`) is preserved unless explicitly requested:

```
AMRSystem(..., loss_model=None, use_ntu_thermal_model=False)
```

- **`loss_model`** (Section 6 below): if given, `W_parasitic` is computed as
  a function of (frequency, field, mdot, Qc) instead of the constant
  `parasitic_fraction · Qc`. This is the fix for the structural blindness
  found in Phase 2's Sobol run (Section 8).
- **`use_ntu_thermal_model`** (Section 5 below): if set, regenerator
  effectiveness ε is computed from bed geometry/NTU instead of the fixed
  constant, via a new `_effective_eps()` helper used everywhere ε
  previously appeared directly (both `Qc` and `eta_2nd_law`).

**Consequence worth being explicit about**: `exergy_eff ≡ η₂(ε)` is still
algebraically degenerate exactly as in Phase 2 (Section 3.4 of the original
walkthrough) — that identity doesn't change just because ε can now vary
with state; it's still not an independent diagnostic, just a differently-
computed input.

---

## 5. New: Regenerator thermal model — `thermal.py` (Phase 4)

Replaces the fixed `regenerator_effectiveness = 0.85` placeholder that made
`mass_regenerator` cosmetically irrelevant to `Qc` through Phase 1–3 (this
is exactly what Phase 3's `optimize.py` run exposed — every Pareto design
sat at the minimum-mass bound because mass bought nothing). Model chain, for
a packed sphere bed (the geometry used in most published AMR prototypes):

```
V_bed = mass_regenerator / (ρ_Gd · (1−φ))                    porosity φ=0.365 default
a = 6(1−φ)/d_p                        A_total = a · V_bed
Re = ρ_f · u_s · d_p / μ_f             Nu = 2 + 1.1·Re^0.6·Pr^(1/3)     (Wakao & Kaguei 1982)
h = Nu·k_f / d_p
NTU = h·A_total / (ṁ·cp_f)
U = ṁ·cp_f / (2·f·mass_regenerator·cp_solid)                   utilization
ε = [NTU/(NTU+2)] · (1 − 0.3·min(U,1))                          clipped to [0, 0.97]
```

The `NTU/(NTU+2)` term is the standard Kays & London balanced periodic-flow
regenerator handbook formula. **Honesty flag carried in the module**: the
`(1 − 0.3U)` utilization-degradation term is qualitatively motivated by
published ε-NTU-U curves (Nielsen et al. 2011; Trevizoli et al. 2016) but is
*not itself independently fit to data* the way `loss_model.py`'s
coefficients are — a placeholder correction, flagged as a target for
replacement with a digitized literature fit.

Verified by running the module directly: effectiveness rises from 0.668 at
0.5 kg to 0.970 at 10–15 kg (mass sweep, f=1 Hz), and from 0.692 at 0.25 Hz
to 0.962 at 4 Hz (frequency sweep, mass=2 kg) — both monotonic and
physically sensible.

---

## 6. New: State-dependent parasitic loss model — `loss_model.py` (Phase 3, extended Phase 6)

Directly motivated by the Phase 2 Sobol finding (Section 8): a constant
`parasitic_fraction` makes `COP_electrical` algebraically blind to
field/frequency/flow. Three standard-scaling-law loss terms, **fit, not
derived from first principles**:

```
W_eddy = k_eddy · f² · (μ0H)²        (eddy-current loss ~ (dB/dt)²)
W_pump = k_pump · ṁ²                   (Darcy-flow pumping power ~ ṁ²)
W_base = base_frac · Qc                 (controls/inverter/bearing overhead)
W_parasitic = W_eddy + W_pump + W_base
```

**CORE calibration** (production default): exactly-determined 3×3 linear
solve against the three Phase 2 benchmark devices (Astronautics, DTU,
Tušek), using their calibrated ṁ and each device's actual field/frequency.
Verified by running the module:

```
k_eddy    =    1.999  W/(Hz²·T²)
k_pump    =  776.887  W/(kg/s)²
base_frac =    0.061  (× Qc)
```

All three fit residuals are ~0.00 W (exactly-determined system, as
expected).

**Phase 6 attempted extension — and its own negative result, reported
loudly rather than hidden.** A 4th device (Okamura & Hirano 2013, found via
a targeted literature search) was added to make this an over-determined
4-point fit. Running `run_extended_diagnostic()` confirms this makes things
*worse*, not better:

```
k_eddy    =   71.303  W/(Hz²·T²)
k_pump    = -907.604  W/(kg/s)²   <- unphysical, negative
base_frac =   -0.065  (× Qc)       <- unphysical, negative
```

Leave-one-out cross-validation on the 4-point set gives errors from +2.6%
(Astronactics) up to **+1638.9%** (Tušek, the smallest device, predicted
from the other three which span up to 2502 W). **Conclusion, and the actual
Phase 6 production decision**: pooling four orders of magnitude of device
scale into one linear loss model doesn't generalize. The `CORE` 3-point fit
remains the default used everywhere else in the codebase (`amr_cycle.py`,
`cascade.py`, `optimize.py`, `giant_mce_analysis.py`); the unstable 4-point
`EXTENDED` fit is exposed only via `run_extended_diagnostic()` for
transparency. A 5th candidate device (Risø/DTU 2011, 30 K span) couldn't
even be calibrated — its reported 1.1 T field can't reach positive Qc at
that span in this model, independent corroboration of the single-stage
span-collapse finding from Phase 1.

**Known stale artifact**: `validation_system.py` (Section 7 below), a
separate module, still prints a hardcoded "118–619%" range for
ideal-vs-electrical COP overprediction. That range was correct for the
original 3-benchmark set, but the 4th benchmark added to
`data/amr_experimental_benchmarks.csv` in this same phase has its own
ideal-vs-literature overprediction of **1102%** — outside the printed
range. This was checked by hand (not by the code) and has not been
corrected in the source as shipped.

---

## 7. System-level validation — `validation_system.py`

Calibration/validation methodology (solve for the one free ṁ that
reproduces reported Qc, then independently check COP_electrical) is
unchanged from Phase 2. The **benchmark set grew from 3 to 4 devices** with
the addition of Okamura & Hirano (2013) (data file update, not a code
change — the module itself is byte-identical to Phase 2):

| Device | Span | COP (lit / ideal / electrical) | Error |
|---|---|---|---|
| DTU rotary Gd (2016) | 10.1 K | 4.20 / 14.88 / 4.60 | +9.6% |
| Tušek single-bed Gd (2010) | 15.0 K | 4.60 / 10.02 / 4.00 | −13.0% |
| Astronautics rotary (2014) | 11.0 K | 1.90 / 13.66 / 4.48 | +135.8% (outlier, flagged in the paper itself) |
| **Okamura & Hirano (2013)** | 5.0 K | 2.50 / 30.06 / 5.46 | +118.3% |

Mean electrical-COP error moves from 52.8% (3 devices) to **69.2%** (4
devices) — see the stale-range note in Section 6 above for the one place
this expanded dataset wasn't fully propagated.

---

## 8. Baseline technologies — `baseline_cooling.py`

Unchanged from Phase 2: vapor-compression at `η₂,vcc = 0.42 · Tc/(Th−Tc)`,
liquid cooling blended between free-cooling and mechanical-assist hours
(`f_econ = 0.6`, `COP_pump_only = 25`).

---

## 9. Sensitivity analysis — `sensitivity.py` (Phase 3 update)

Same Sobol/Saltelli method as Phase 2, but now runs in **two modes** and
writes both, so the fix is directly demonstrated rather than asserted:

**Phase 2 mode** (constant `parasitic_fraction`, `results/sobol_results_phase2_constant.txt`):
same degenerate result as before — `parasitic_fraction` ST≈0.999,
everything else ST≈0.

**Phase 3 mode** (state-dependent `loss_model`, `results/sobol_results.txt`),
verified by running the module:

```
frequency_Hz          ST=0.684
fluid_mdot_kgs         ST=0.221
mu0H_max_T              ST=0.095
regen_effectiveness      ST=0.031
parasitic_fraction         ST=0.000   (unused in this mode, kept as a no-op slot)
```

Field/frequency/flow now carry **real** sensitivity in COP_electrical,
because raising frequency or field no longer buys more cooling for free —
it also raises the eddy-current loss term via `loss_model.py`. Caveat
carried directly in the output: the loss coefficients themselves come from
an exactly-determined 3-point fit, so treat the *magnitude* of these new
sensitivities as illustrative, not a converged final number.

---

## 10. Surrogate model — `rsm.py`

Unchanged from Phase 2 (byte-identical). Still the Qc-targeted quadratic
response surface, R² = 0.9423 held-out, used to keep future optimizer inner
loops fast — though `optimize.py` (Section 12) in fact calls
`AMRSystem.run()` directly rather than this surrogate, since NSGA-III's
population sizes here were cheap enough not to need it.

---

## 11. New: Multi-stage cascade design — `cascade.py` (Phase 3)

Directly answers the Phase 1 finding that a single AMR stage at μ0H=2T
collapses to zero `Qc` above ~16 K span — inside the project's 5–20 K
ASHRAE target range. N identical single-stage AMR modules in series, each
covering an equal share of the total span; the coldest stage sets the
deliverable `Qc` (the bottleneck), and every downstream stage is re-run at
its local span and scaled to match that same `Qc` (steady-state series
constraint):

```
W_total = Σᵢ W_i(Qc_target, span/N)          COP_cascade = Qc_target / W_total
```

Verified by running the module (Gd, 2 T/stage, 5 kg/stage, NTU thermal
model on): 1-stage COP at 16 K span drops to 2.31 vs. VCC's 7.64, and
1-stage becomes infeasible above 16 K, while 2/3/4-stage designs remain
feasible out to 20 K — staging genuinely recovers the span the single-stage
design loses, though at a COP cost relative to the (still-losing) baseline
comparisons.

**Simplification flagged directly in the module**: every stage uses the
same material (Gd) for tractability. Real cascade designs use graded
Curie-temperature materials per stage — not modeled here, an explicit
deferred item.

Running the Gd5Si2Ge2 variant of this same script shows the giant-MCE
material failing at *every* span in the ASHRAE-anchored sweep (all stage
counts): consistent with Section 12's finding that its favorable window
sits below this operating range entirely, regardless of staging.

---

## 12. New: Giant-MCE analysis — `giant_mce_analysis.py` (Phase 5)

Now that `first_order_mce.py` gives Gd5Si2Ge2 a physically appropriate
model, this answers the question Phase 4 left open: does the giant-MCE
material change the COP-competitiveness conclusion? Verified end to end by
running the module:

- Gd's own peak-effect temperature: **296.5 K** — inside the ASHRAE range.
- Gd5Si2Ge2's own peak-effect temperature: **286.4 K** — ~4.6 K *below* the
  ASHRAE range.
- At the ASHRAE point (291 K, 10 K span): Gd delivers Qc=1241.5 W,
  COP_elec=7.42; Gd5Si2Ge2 collapses to **Qc=0**.
- At Gd5Si2Ge2's own favorable point (281.4 K, 10 K span): Gd5Si2Ge2
  delivers Qc=5319.3 W, COP_elec=7.76 (Gd fails at that same point, for the
  mirror-image reason).

**Conclusion, held to precisely (not oversold)**: the giant-MCE effect is
real and large *within its own narrow transition window*, but that window
is mistargeted for data-center duty as-is. This does not overturn the
earlier Gd-trails-baselines conclusion. What it *does* support: the
Gd5(SixGe1-x)4 family is composition-tunable across ~20–335 K (Pecharsky &
Gschneidner 1997), so a composition tuned into the ASHRAE range — if it
retains first-order/giant character there — is a genuinely open,
materials-synthesis question outside what this simulation suite alone can
answer. Also noted: even correctly targeted, Gd5Si2Ge2's COP_electrical
(7.76) is close to Gd's own (7.42) despite ~4× the Qc — consistent with the
Section 9 finding that COP is driven mainly by loss-model terms, not which
material is loaded into the regenerator; a bigger MCE mostly buys more
capacity per kg, not better efficiency.

---

## 13. New: Multi-objective optimization — `optimize.py` (Phase 3, cost objective updated Phase 6)

NSGA-III (Deb & Jain 2014, via the pymoo implementation) over 5 design
variables (field, frequency, ṁ, regenerator mass, effectiveness), 3
objectives: maximize COP_electrical (using the Phase 3 loss model, not the
Phase 2 constant), maximize Qc, minimize cost. Verified by running the
module (pop=60, gen=40, 34 Pareto-optimal designs found):

```
Best electrical COP    H=2.98T f=0.32Hz mdot=0.057kg/s mass=13.37kg -> COP=7.89, Qc=1566W,  cost=$5047
Best cooling capacity  H=3.00T f=3.17Hz mdot=0.498kg/s mass= 8.54kg -> COP=6.58, Qc=13519W, cost=$3244
Lowest cost            H=2.58T f=0.34Hz mdot=0.024kg/s mass= 1.00kg -> COP=7.55, Qc= 462W,  cost=$331
Knee point (balanced)  H=2.99T f=1.94Hz mdot=0.496kg/s mass= 5.04kg -> COP=6.76, Qc=12178W, cost=$1912
```

`mass_regenerator` now genuinely spans 1.00–14.75 kg across the front
(rather than sitting at the floor as in the Phase 3 constant-ε run) because
`thermal.py`'s NTU model makes more regenerator mass a real Qc/ε lever, not
free waste.

**Phase 6 cost-objective update, confirmed in source**: `cost_index()`
calls `economics.material_cost()` directly (verified: `optimize.py` imports
`material_cost` from `core.economics` and its only use is inside
`cost_index`), which is grounded in Bjørk et al. (2011)'s per-kg costs
(Section 14) rather than the earlier $175/kg Franco et al. placeholder.

**Known stale artifact**: the module's own docstring still describes the
cost objective in its old Phase-3 form ("Gd material cost ~$175/kg (Franco
et al. 2018) + a magnet-cost term scaling with field² × regenerator mass")
— this text was not updated when the Phase 6 code change landed, so the
docstring and the actual `cost_index()` implementation now disagree. Same
class of drift as the stale range noted in Section 6.

---

## 14. New: Refrigerant-free emissions comparison — `emissions.py` (Phase 5)

Quantifies the one categorical, non-COP-dependent advantage AMR has: no
HFC/HFO refrigerant charge, hence zero refrigerant-leak GWP and zero
F-gas-type phase-out liability. Refrigerant emissions use standard AR5
100-year GWP values (R-410A: 2088, R-134a: 1430) with a 4%/year leak-rate
assumption; operational emissions use a representative grid CO2 intensity.
Verified by running the module (100 kW, illustrative COPs AMR=5.0,
VCC=12.0, Liquid=20.0):

```
Magnetic (AMR)      refrigerant=  0.00 tCO2e/yr  operational= 87.07  total= 87.07
Vapor-compression   refrigerant=  3.34 tCO2e/yr  operational= 36.28  total= 39.62
Liquid cooling      refrigerant=  2.29 tCO2e/yr  operational= 21.77  total= 24.06
```

**Held to precisely, not oversold, in the module's own output**: at these
representative COPs, AMR's *lower* COP (the Phase 1–4 finding) makes its
operational emissions the highest of the three — the refrigerant-free
story is real but does not, on its own, flip the total-emissions
comparison unless AMR's COP gap is closed first.

---

## 15. Economics — `economics.py` (Phase 6 grounding update)

The linear CAPEX/OPEX bookkeeping layer from Phase 2 is unchanged, but a
new bottom-up materials-cost function replaces an earlier loosely-sourced
$175/kg placeholder:

```
material_cost(μ0H_max, mass_regenerator):
    magnet_mass = 3.0 · μ0H_max · mass_regenerator      (rough fit to 2 worked
                                                            examples, NOT a
                                                            validated scaling law)
    return $40/kg · magnet_mass + $20/kg · mass_regenerator
```

Unit costs ($40/kg NdFeB magnet, $20/kg Gd) and the magnet:MCM mass-ratio
approximation both come from Bjørk, Bahl & Smith's dedicated
magnetic-refrigerator cost-minimization study (Int. J. Refrig. 34 (2011)
1805–1816). Re-running the Pareto front with this shifts materials-only
cost estimates to roughly $300–5,000, consistent with Bjørk et al.'s own
small-device examples ($7–35) scaled up.

**Honesty flag carried in the module**: this is a **materials-only floor**
— heat exchangers, pumps, motor/drive, and controls are excluded, and Bahl
et al. (2014) note these typically dominate total AMR system cost. The
older `$/kW` placeholder used by `main.py`'s comparison table is explicitly
noted as *not yet reconciled* with this bottom-up figure — a concrete,
named next step, not silently patched over.

---

## 16. Top-level driver — `main.py`

Unchanged from Phase 2: single-stage AMR vs. VCC vs. liquid cooling across
the ASHRAE 5–20 K span sweep at Tc=18°C, writing
`results/comparison_table.csv`. Because all new plug-in points
(`loss_model`, `use_ntu_thermal_model`) default to off, this file's output
is **byte-identical** to the Phase 2 version — the single-stage
conclusion (AMR trails VCC on electrical COP, and its Qc collapses above
~16 K span) stands unchanged; `cascade.py` (Section 11) is the module that
addresses the span limitation, not a rewrite of `main.py` itself.

---

## Known inconsistencies carried into this final state (not fixed as shipped)

Two places where new data/code landed but an adjacent hardcoded description
wasn't updated to match — found by executing the modules and checking the
arithmetic by hand, not flagged by the code itself:

1. `validation_system.py`'s printed "118–619%" ideal-vs-electrical
   overprediction range predates the 4th benchmark device; the true range
   with that device included is **118–1102%**.
2. `optimize.py`'s module docstring still describes the old $175/kg
   Franco-et-al.-based cost objective; the actual `cost_index()`
   implementation was updated to call `economics.material_cost()`
   (Bjørk et al.-grounded) and the two no longer agree.

---

## Equation summary (all in one place)

| Quantity | Equation |
|---|---|
| Molecular field constant | λ = 3kBTc / [Ng²μB²J(J+1)μ0] |
| Effective field | Heff = H + λM |
| Magnetization (Gd, 2nd-order) | M = NgμBJ·B_J(x), x = gμBJμ0(H+λM)/kBT |
| Landau free energy (Gd5Si2Ge2, 1st-order) | f(m,τ,h) = (A/2)(τ−1)m² + (B/4)m⁴ + (C/6)m⁶ − hm |
| Isothermal ΔS (2nd-order) | ΔS_M(T,H) = S_M(T,H) − S_M(T,0) |
| Isothermal ΔS (1st-order) | ΔS_M(τ,h) = −(A/2)[m(τ,h)² − m(τ,0)²]·NkB |
| Adiabatic ΔT | ΔT_ad = −T·ΔS_M / C_total |
| Lattice Cp (Debye) | C_lat = 9nR(T/θD)³∫x⁴eˣ/(eˣ−1)²dx |
| Magnetic Cp (2nd-order only) | C_mag = T·dS_M/dT |
| Regenerator eff. (NTU model) | ε = [NTU/(NTU+2)]·(1−0.3·min(U,1)) |
| NTU | NTU = h·A_total/(ṁ·cp_f) |
| Utilization | U = ṁ·cp_f / (2f·m_reg·cp_solid) |
| Cooling capacity | Qc = ε·ṁcp·ΔT_ad,noload·max(0, 1−span/2ΔT_ad,noload) |
| Carnot work | W_carnot = Qc(Th/Tc−1) |
| 2nd-law efficiency | η₂ = 0.35+0.20ε |
| Magnetic work | W_mag = W_carnot/η₂ |
| Parasitic power (state-dependent) | W_parasitic = k_eddy·f²H² + k_pump·ṁ² + base_frac·Qc |
| Ideal COP | COP = Qc/W_mag |
| Electrical COP | COP_e = Qc/(W_mag + W_parasitic) |
| Exergy efficiency | ηₑₓ = COP/COP_carnot ≡ η₂ (still degenerate) |
| Carnot COP | COP_carnot = Tc/(Th−Tc) |
| Cascade COP | COP_cascade = Qc_target / Σᵢ W_i |
| Materials cost floor | cost = $40/kg·(3·μ0H·m_reg) + $20/kg·m_reg |
| Refrigerant emissions | tCO2e = charge_kg·leak_rate·GWP / 1000 |
| Operational emissions | tCO2e = (capacity·hours·load/COP)·CO2_per_kWh / 1000 |
