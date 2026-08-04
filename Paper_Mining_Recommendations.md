# Paper-Mining Pass: What to Add, What to Validate

Cross-referenced `magcool-dc`'s current state (README/ROADMAP/Literature_Review.md,
`data/amr_experimental_benchmarks.csv`, `core/*.py`) against the papers in `Papers/`
that the codebase's own roadmap had not yet fully mined. Goal: concrete model
additions, each anchored to a specific number from a specific paper — not
general "consider adding X" suggestions.

Four papers below yielded genuinely new, usable content. One "new" file turned
out to be a duplicate. Everything else already in `Papers/` is already reflected
in ROADMAP.md's Phase 7 entries (Tušek 2013 Figs. 10–11 digitization, full HX/
pump/motor BOM cost) — those remain open for the reasons already stated there;
repeating them here would just restate the roadmap.

---

## 1. Add flow-waveform asymmetry (blow fraction) to the AMR cycle model — real, sizeable validation target

**Source:** Masche, Liang, Engelbrecht, Bahl, *"Improving magnetic cooling
efficiency and pulldown by varying flow profiles,"* Applied Thermal
Engineering 215 (2022) 118945 (`Papers/AMR Theory and Modeling/`). DTU rotary
device, 13 trapezoidal beds, 295 g Gd spheres/bed, solenoid-valve-controlled
blow fraction.

**What's missing:** `core/amr_cycle.py` and `core/thermal.py` have no notion
of *blow fraction* — the fraction of the AMR cycle period during which fluid
flows in the cold-to-hot vs. hot-to-cold direction. The model implicitly
assumes a symmetric 50/50 split. This paper shows that's a real, first-order
lever, not a second-order detail:

> At a 16 K span, utilization U = 0.32, 1.4 Hz: increasing the blow fraction
> from 25.0% to 41.6% raised cooling capacity from **70 W to 330 W**
> (a ~4.7× increase) and second-law efficiency from **2.6% to 17.4%**.
> Best blow fraction found across both the 6 K and 16 K spans tested: **~41.6%**.
> Lower blow fractions favor faster temperature pulldown — the system reached
> ~14 K span **~30% faster** at a lower blow fraction, trading steady-state
> capacity for pulldown speed.

**Concrete addition:** parameterize `amr_cycle.py`'s cycle model with a
`blow_fraction` (or asymmetric cold-blow/hot-blow duration split) rather than
a fixed symmetric blow, and expose it in `optimize.py` as another
NSGA-III decision variable alongside frequency/field/mdot. This is a genuine
new degree of freedom the model currently can't represent at all, and it's
the single largest reported model gap I found — a 4.7× Qc swing from one
operating parameter, at a device geometry already close to the ones in the
current benchmark set.

**Validation target:** the (70 W, 2.6%) → (330 W, 17.4%) pair at fixed
16 K span/U=0.32/1.4 Hz is a clean two-point curve to calibrate a
blow-fraction-dependent effectiveness or loss term against — same
calibrate-then-validate pattern as `validation_system.py` already uses.

---

## 2. No benchmark device currently validates the parallel-plate geometry model

**Source:** Tušek, Kitanovski, Zupan, Prebil, Poredoš, *"A comprehensive
experimental analysis of gadolinium active magnetic regenerators,"* Applied
Thermal Engineering 53 (2013) 57–66 — the *other* Tušek 2013 paper already
flagged in ROADMAP.md for its unfinished Figs. 10–11 digitization.

**What I checked that the roadmap notes didn't emphasize:** every row in
`data/amr_experimental_benchmarks.csv` is a packed-bed (or layered packed-bed)
device. `core/geometry_analysis.py` (Phase 7) added a
`regenerator_effectiveness_parallel_plate()` model, but nothing in the
benchmark set actually exercises it against a real parallel-plate device.
This paper is the natural candidate: its highlight finding, confirmed in the
abstract text (not requiring Fig. 10/11 digitization), is

> parallel-plate AMR with ~25% porosity, plates oriented parallel to the
> field, 1.15 T: **20 K temperature span**, the largest reported for a
> parallel-plate AMR at this field at the time of publication. Parallel-plate
> AMRs showed higher COP than packed-bed AMRs across the comparison.

The exact (Qc, COP) pair at that operating point is only in Figs. 10–11
(same digitization blocker ROADMAP.md already documents — I confirmed this
directly rather than re-attempting the extraction), so this can't yet become
a full calibrate-then-validate CSV row. But it's worth flagging as a
**named, specific target** for whoever does pick up the Fig. 10/11
digitization: it isn't just "more curve data," it's the only available
opportunity to validate `geometry_analysis.py`'s parallel-plate branch at
all, which currently ships with zero device-level validation.

---

## 3. A second giant-MCE material family is available and better field-effective than Gd5Si2Ge2 — worth adding alongside LaFeSiH

**Source:** Hanggai, Yibole, Guillou, Kwakernaak, van Dijk, Brück,
*"Preparation of Fe-rich giant magnetocaloric (Mn,Fe)2(P,Si) ribbons and
calorimetric analysis of the first-order magnetic transition,"* Acta
Materialia 302 (2026) 121677 (`Papers/Magnetocaloric effect and materials
physics/Impact of F and S Doping...pdf` — the file's actual title differs
from its filename; confirmed via the PDF's own header).

**Numbers extracted directly from the text** (no digitization needed):

| Composition | \|ΔS_max\| | Field change | Note |
|---|---|---|---|
| Mn0.60Fe1.30P0.66Si0.34 (parent) | 12 J kg⁻¹ K⁻¹ | — | baseline |
| Mn0.62Fe1.28... | 16.66 J kg⁻¹ K⁻¹ | 2 T | +39% vs. parent |
| highest-Mn/Si composition tested | 17.61 J kg⁻¹ K⁻¹ | 2 T | "40% enhancement... compared to parent compound" (paper's own framing) |

Reconstructed latent heat ≈ 19.97 kJ/kg for the MnFe2P-based system,
compared against an independently-sourced elastic strain energy of
30.89 kJ/kg — the paper's own cross-check, included here only as provenance
context, not something to re-derive.

**Why this is worth adding, not just noting:** (Mn,Fe)2(P,Si) is the other
major room-temperature first-order giant-MCE family besides Gd5Si2Ge2 and
La(Fe,Si)13Hy, and it's Curie-tunable by composition (same "graded cascade"
mechanism `cascade.py`'s `GradedFamily` abstraction already supports for
`GD_FAMILY`/`LAFESIH_FAMILY`). Given `cascade.py` was explicitly generalized
in Phase 9's addendum to a pluggable-family design specifically so a third
family could be added without re-deriving the staging logic, this is close
to a drop-in: a `MNFEPSI_FAMILY` alongside the existing two, calibrated to
the ~17.6 J/(kg·K) at 2 T peak entropy change figure above (same grid-search
calibration method already used for `GD5SI2GE2_FIRST_ORDER` and
`LAFESIH_FIRST_ORDER`).

**Honesty flag to carry over if this is added:** this paper reports
ΔS_max (Maxwell-relation entropy change), not a direct ΔT_ad measurement —
exactly the indirect-vs-direct gap that `giguere_validation.py` already
found overstates the model's ΔT_ad by ~2.4× for Gd5Si2Ge2. No equivalent
direct-measurement cross-check paper for (Mn,Fe)2(P,Si) is in this corpus,
so the same `dTad_correction`-style caveat would need to be carried
forward rather than assumed away.

---

## 4. Confirmed duplicate: "nine-layer active regenerator" is not a new paper

`Papers/AMR Theory and Modeling/Performance evaluation of a nine-layer active
regenerator.pdf` and `Papers/AMR systems and prototypes/The performance of a
large-scale rotary magnetic refrigerator.pdf` are **the same PDF** (Jacobs,
Auringer, Boeder, Chell, Komorowski, Leonard, Russek, Zimm, Int. J. Refrig.,
DOI 10.1016/j.ijrefrig.2013.09.025 — the Accepted-Manuscript version of the
paper already cited in `amr_experimental_benchmarks.csv` as
`Astronautics_rotary_2014`). Confirmed by opening both and comparing text —
identical title, DOI, and abstract (3042 W zero-span / 2502 W at 11 K span,
the exact numbers already in the CSV). Worth noting only so nobody spends
time later searching this file for "nine-layer" content that isn't there —
the six-layer Astronautics bed already in the benchmark set is the real
device this file describes.

---

## 5. Lower-priority, review-level directions (context, not action items)

From Zhang, Wu, He, Wang, Yu, *"Solutions to obstacles in the
commercialization of room-temperature magnetic refrigeration,"* Renewable
and Sustainable Energy Reviews 143 (2021) 110933 — a survey of what's
limiting AMR commercially, useful as *framing* for the project's discussion
section rather than as a model addition:
- Heat-transfer-fluid enhancement (nanofluids, liquid metals) and MCM
  shaping/high-conductivity inserts are identified as the two most
  practical near-term heat-regeneration improvements — relevant context if
  `thermal.py`'s NTU model is ever extended beyond water/Wakao-Kaguei, but
  no fluid-property numbers were extracted since nothing in the current
  model needs them yet.
- Fully-solid-state and multi-caloric (combined field) cycles are flagged as
  the likely long-term direction — out of scope for this repo's AMR-cycle
  architecture, noted only for the paper's own discussion section if it
  addresses future work.

No new numeric targets came out of this one; it's a discussion-section
citation, not a validation source.

---

## Summary — priority order if picking one thing to build next

1. **Blow-fraction asymmetry in `amr_cycle.py`** (§1) — largest reported
   effect size (4.7× Qc), two clean numbers to calibrate against, fits the
   existing calibrate-then-validate pattern directly.
2. **(Mn,Fe)2(P,Si) as `MNFEPSI_FAMILY`** (§3) — genuinely new material,
   near-drop-in given the Phase 9 `GradedFamily` refactor, real ΔS_max
   numbers in hand.
3. **Flag the parallel-plate validation gap** (§2) for whoever eventually
   does the Tušek 2013 Fig. 10/11 digitization already on the roadmap —
   redirect that effort specifically at closing this gap rather than
   generic "more curve data."
