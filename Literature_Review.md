## Literature Review — Magnetocaloric Cooling for Data Centers

*Updated through Phase 23 (elastocaloric baseline comparison). Corrections
to earlier drafts are marked **CORRECTED**; new entries added to cover
Phases 15–23 are marked **NEW**. Where this project's own copies of the
two reference books (Kitanovski et al. 2015; Tishin & Spichkin 2003) could
not supply a cited section, that is stated explicitly rather than silently
worked around — see the consolidated book-access note in Section 6.*

### 1. Magnetocaloric Effect and Magnetocaloric Materials

#### Giant Magnetocaloric Effect

* **Pecharsky & Gschneidner**, *Physical Review Letters*, **78**, 4494–4497 (1997).

  * First reported the **giant magnetocaloric effect (GMCE)** in **Gd₅Si₂Ge₂**.
  * Demonstrated that first-order magnetostructural phase transitions can produce substantially larger entropy and temperature changes than conventional second-order ferromagnets.
  * Established the foundation for modern room-temperature magnetic refrigeration research.
  * Also reports composition tunability across the Gd₅(SixGe1-x)₄ pseudo-binary system — ordering temperature spans roughly 30 K to 335 K (Gd₅Si₄ itself) as Si:Ge ratio is varied. Used directly by `core/first_order_mce.py`'s `GD_FAMILY` composition-tuning and, in `core/giant_mce_analysis.py`, to frame the "untested direction" for closing the ASHRAE-range COP gap once Gd₅Si₂Ge₂'s own peak is found to sit ~5 K below the data-center supply range (see Section 1's Gadolinium-vs-giant-MCE finding, reproduced in `results/giant_mce_analysis.txt`).

#### Gadolinium as the Benchmark Material

* **Pecharsky & Gschneidner**, *Journal of Magnetism and Magnetic Materials*, **200**, 44–56 (1999).

  * Comprehensive review of the magnetic, thermodynamic and magnetocaloric properties of gadolinium.
  * Identified Gd as the standard reference material for room-temperature magnetic refrigeration owing to its Curie temperature (~294 K), making it directly relevant to data-center operating temperatures.

#### Experimental Validation of Magnetocaloric Properties

* **Dan'kov, Tishin, Pecharsky & Gschneidner**, *Physical Review B*, **57**, 3478–3490 (1998).

  * Measured the adiabatic temperature change (ΔTₐd) of gadolinium under magnetic fields of **1 T, 2 T and 5 T** near its Curie temperature.
  * These measurements remain one of the most widely used experimental datasets for validating theoretical and numerical MCE models.
  * Also reports the Curie-point transition temperature increasing almost linearly with field, ~6 K/T above 2 T up to 7.5 T (zero-field TC = 294 ± 1 K, confirmed by four independent techniques). Used in `core/validation.py`'s `run_curie_shift_check()` as a held-out check of the mean-field model's own emergent peak-shift-with-field — the model does NOT reproduce this rate (fitted ~0.0 K/T vs. the reported ~6 K/T across 2.0–7.5 T), a genuine, documented mean-field limitation rather than a re-confirmed match (`main.py` step 1).
  * **NEW (Phase 22 item 1):** the same three field points anchor `core/inhomogeneous_broadening.py`'s Gaussian polycrystalline Tc-broadening sensitivity check. Broadening the model's single sharp Tc into a Gaussian ensemble (σ_Tc up to 5 K) narrows the worst-field error against these three points from +48.9% (sharp model) to +20.9% (σ_Tc=5 K) — but the improvement is one-sided: the 5 T point moves the *other* direction, from −7.5% to −14.2%, consistent with global smoothing of the transition rather than a field-selective fix. The 1 T error is still falling monotonically at the edge of the swept range, so no interior optimum was located, and no real digitized ΔTₐd(T) curve exists in this repo's corpus to fit σ_Tc against — a sensitivity finding, not a calibration.

* **Giguère, Foldeaki, Ravi Gopal, Chahine, Bose, Frydman & Barclay**, *Physical Review Letters*, **83**, 2262 (1999).

  * Primarily a direct-measurement Gd₅Si₂Ge₂ paper (already used by `core/giguere_validation.py`), but its methods section also independently cross-checks pure Gd: ΔTₐd ≈ 10.5–11.5 K at 5 T (this paper's high-purity sample vs. AMES laboratory's measurement, agreeing within 1 K) and ≈ 12–13 K at 7 T (industrial- vs. high-purity Gd respectively), agreeing with Brown (1976)'s independently reported 14 K.
  * Used in `core/validation.py`'s `run_giguere_gd_extension()` as a second, independent Gd dataset extending validation to 7 T. The model (calibrated to Dan'kov et al.'s higher 5 T value) sits outside both ranges: +22.8% vs. the 5 T midpoint, +33.7% vs. the 7 T midpoint — a real disagreement between two published Gd measurements, reported rather than reconciled.
  * Its own **direct** ΔTₐd measurement for Gd₅Si₂Ge₂ (10.0 K at 7 T; 9.9 K by an independent Clausius–Clapeyron cross-check) is the primary validation target for `core/giguere_validation.py`'s check of the first-order Landau model (see next entry). Its own indirect (Maxwell-relation) figure for the same point is 14.9 K — a 1.49× indirect-vs-direct gap the paper itself reports, used as the point of comparison for this repo's own, larger 2.42× gap.

* **Pecharsky & Gschneidner**, *Physical Review Letters*, **78**, 4494 (1997).

  * The original Gd₅Si₂Ge₂ "giant" magnetocaloric effect discovery paper — already cited above for `GD5SI2GE2_FIRST_ORDER`'s Tc and J, but its own text also reports a heat-capacity-derived comparison not previously extracted: "the ΔTₐd values of Gd₅Si₂Ge₂ are larger than the corresponding ΔTₐd values for Gd by about 30%, comparing the peak values, regardless of the temperature" (Fig. 6, field changes of 0→2T and 0→5T).
  * Used in `core/giguere_validation.py`'s `run_pecharsky_ratio_check()` as a second, independent primary source (heat-capacity-based, not the pulse-field-thermometry method Giguère et al. used) and a second field range (2 T/5 T, vs. the single 7 T point `DTAD_CORRECTION_FACTOR` is fit to). Result: the *raw* (uncorrected) model's peak ratio at 5 T (~1.24) sits close to this paper's ~1.30 — but applying the Giguère-derived correction factor overcorrects it to ~0.51, predicting Gd₅Si₂Ge₂ *underperforms* plain Gd, which contradicts the giant-MCE premise. Documented as evidence the 7 T correction factor should not be extrapolated as field-independent.

#### Magnetocaloric Theory

* **de Oliveira & von Ranke**, *Physics Reports*, **489**, 89–159 (2010).

  * Comprehensive review of thermodynamic models describing the magnetocaloric effect.
  * Highlights the limitations of mean-field theory near the Curie temperature due to critical magnetic fluctuations and discusses improved theoretical formulations — the same limitation `core/validation.py`'s Curie-shift check (above) and `core/mce_material.py`'s Section 1.5 honesty flag (Gd₅Si₂Ge₂ is structurally invalid under this framework) independently surface in this repo's own model.

* **Bean & Rodbell**, *Physical Review*, **126**, 104 (1962).

  * The classical magnetoelastic Landau-expansion treatment of a first-order ferromagnetic transition (MnAs) — the same phenomenological family (order-parameter free energy with a negative quartic term stabilized by a positive sextic term) `core/first_order_mce.py`'s extended Landau model for Gd₅Si₂Ge₂/La(Fe,Si)₁₃Hy/(Mn,Fe)₂(P,Si) follows, expressed directly as an order-parameter expansion rather than via self-consistent lattice strain.

#### Alternative Magnetocaloric Materials

* **Guo et al.**, *Applied Physics Letters*, **78**, 1142 (1997).
* **Phan & Yu**, *Journal of Magnetism and Magnetic Materials*, **308**, 325–340 (2007).

  * Investigated perovskite manganites and other rare-earth alloys exhibiting significant room-temperature magnetocaloric effects.
  * Discuss material selection criteria including entropy change, hysteresis, thermal conductivity and operating temperature. La₀.₇Ca₀.₃MnO₃, from this materials class, is the third candidate `core/passive_regenerator_analysis.py`'s Phase 21 passive-regenerator comparison uses precisely because its Curie point (267 K) sits *outside* the ASHRAE operating window — the "zero-gain control case" that confirms the passive-regenerator benefit is a genuine Tc-alignment effect and not a fixed per-material bonus (see Section 2's Passive/Hybrid Regenerators entry).

#### (Mn,Fe)₂(P,Si) Giant Magnetocaloric Materials

* **Hanggai, Yibole, Guillou, Kwakernaak, van Dijk & Brück**, *Acta Materialia*, **302**, 121677 (2026).

  * Preparation and calorimetric analysis of Fe-rich melt-spun Mn₀.₆₀₊ₓFe₁.₃₋ₓP₀.₆₆₋ᵧSi₀.₃₄₊ᵧ (0≤x≤0.08, x=2y) ribbons.
  * Reports Curie temperature increasing linearly from 295.3 K (parent) to 331.2 K (highest Mn/Si tested) across five measured compositions — a directly-measured, room-temperature-adjacent tunability window.
  * Reports peak |ΔS_max| ~17.6 J/(kg K) at a 2T field change for the highest-Mn/Si composition, cross-validated by two independent methods (16.66 J/(kg K) calorimetric, 17.61 J/(kg K) magnetization) — a ~40% enhancement over the 12 J/(kg K) parent-compound value at the same field.
  * Used to add `MNFEPSI_FAMILY`, a third pluggable Curie-graded material family alongside the Gd₅(SixGe1-x)₄ and La(Fe,Si)₁₃Hy families — notable because its Tc window [295.3, 331.2] K sits almost entirely *at or above* the ASHRAE data-center supply range, the opposite tension from the Gd₅(SixGe1-x)₄ family's ceiling sitting just below it. In practice, the six-way material comparison (`core/material_family_comparison.py`) finds this family's window does not cover the representative T_mid≈296 K operating point at any tested span, so it falls back to plain Gd there rather than being independently ranked — reported explicitly in `results/material_family_comparison.csv` rather than silently substituted.

---

## 2. Active Magnetic Regenerator (AMR) Systems

#### Active Magnetic Regenerator Concept

* **Barclay**, US Patent **4332135** (1982).

  * Introduced the Active Magnetic Regenerator (AMR) cycle.
  * Demonstrated how cyclic magnetization and demagnetization combined with regenerative heat exchange could achieve practical magnetic refrigeration.

#### AMR Prototype Comparisons (Secondary Sources)

* **Kamran, Ahmad & Wang**, *Renewable and Sustainable Energy Reviews*, **133**, 110247 (2020).

  * Review with a comparative summary table (Table 2) of ~12 AMR prototypes (Chubu Electric/Toshiba, Institute of Tech. Chubu, Nanjing University, Riso Lab, POLO, U. Salerno, G2E Grenoble, U. Tokyo, Wroclaw, Teyber et al., among others), citing each device's own primary paper rather than reporting original measurements.
  * Table did not extract cleanly as text (merged/multi-value cells) — read from a rendered page image instead of `pdftotext`/`pdfplumber` output.
  * Used to add the Chubu Electric/Toshiba two-field-point row (`ChubuToshiba_Gd_2016_4T`/`_2T`, original device ref [69] in this review, not itself in this repo's `Papers/`) — SECONDARY SOURCE, same caveat as the `Okamura_Hirano_2013` row. Several other rows in the same table (Institute of Tech. Chubu's 540W near-zero-span extreme; a second, independent Riso Lab data point with different regenerator count/frequency than the existing `Risoe_DTU_Gd_2011` row; Teyber et al.'s 100K-span superconducting-magnet outlier) were identified but left un-added — flagged in `ROADMAP.md` as a scope decision pending either a clearer boundary or the primary source papers.

* **Greco, Aprea, Maiorino & Masselli**, *International Journal of Refrigeration* (2019).

  * "A review of the state of the art of solid-state caloric cooling processes at room-temperature before 2019" — Table 2 lists ~25 AMR prototypes built 2009–2018 (Q̇_ref,max, ΔT_span,max, field, MCM/HTF), larger and cleaner than the Kamran/Ahmad/Wang table above, though the PDF carries a diagonal "ACCEPTED MANUSCRIPT" watermark that bleeds letter fragments into extracted text on several rows — checked against surrounding column structure, but not a rendered-image read like the Kamran table.
  * Used to add `Cooltech_2013_rotary` (42 K span — the largest span in this benchmark set — confirmed directly in the review's own body text, not table-only, and used as a capacity-only feasibility stress test) and `DTU_MagQueen_2018` (a La(Fe,Mn,Si)₁₃Hz-based heat pump, giving this repo's first LAFESIH-material benchmark point independent of Astronautics; Qc/COP values are DERIVED from the paper's own reported heating power/COP, flagged in the CSV row). The table's Astronautics 2014 entry lists ΔT_span,max=18K, unconfirmed elsewhere in this corpus and NOT added as a third Astronautics data point — flagged in `ROADMAP.md` as worth checking against the primary Jacobs et al. paper before use.

#### Performance Modelling

* **Tušek, Zupan, Šarlah, Prebil & Poredoš**, *International Journal of Refrigeration*, **33**, 294–300 (2010).

  * **CORRECTED:** this citation previously shared a "characteristic-curve models" description with Nielsen et al. (2011) below. Identifying the underlying PDF (`Development of a rotary magnetic refrigerator.pdf`) confirmed this is actually a mechanical/magnet-design paper for an earlier Ljubljana AMR prototype, NOT a performance-testing or characteristic-curve paper. It reports zero Qc/span numbers.
  * Content is the permanent-magnet assembly design and a "pros and cons" table of the mechanical build (shaft-seal leakage causing rotational friction/heat generation, large magnet-structure weight, assembly/disassembly complexity).
  * Useful only as qualitative engineering-realism context — e.g. for discussing real mechanical parasitic losses that `core/amr_cycle.py`'s idealized 0-D cycle doesn't capture — not as a numeric validation source. (Not to be confused with `Development of a novel rotary magnetic refrigerator.pdf` = Lozano et al. 2016, already the primary source behind the `Lozano_POLO_UFSC_2016` CSV rows.)

* **Nielsen et al.**, *International Journal of Refrigeration*, **34**, 603–616 (2011).

  * Developed characteristic-curve models relating cooling capacity, utilization factor, operating frequency, pressure drop and regenerator effectiveness.
  * These models provide widely accepted benchmarks for predicting AMR performance and validating numerical simulations.

#### Flow-Profile / Blow-Fraction Effects

* **Masche, Liang, Engelbrecht & Bahl**, *Applied Thermal Engineering*, **215**, 118945 (2022).

  * Experimentally varied flow profiles (blow fraction — the fraction of the cycle period spent in cold-to-hot vs. hot-to-cold flow) on a DTU rotary AMR device (13 trapezoidal beds, 295 g Gd spheres/bed, solenoid-valve-controlled blow fraction).
  * At a fixed 16 K span/U=0.32/1.4 Hz operating point, increasing blow fraction from 25.0% to 41.6% raised cooling capacity from 70 W to 330 W (~4.7×) and second-law efficiency from 2.6% to 17.4%; lower blow fractions instead favored faster temperature pulldown.
  * Used in `core/amr_cycle.py` to add `blow_fraction` as a new AMR cycle degree of freedom (`BLOW_FRACTION_MASCHE`, `_blow_fraction_multiplier()`) and as a 6th (later 7th, alongside particle diameter — see the Regenerator Geometry entries below) NSGA-III design variable in `core/optimize.py` — a real, sizeable model gap this repo had no prior representation of at all (Qc/second-law efficiency previously assumed a fixed symmetric 50/50 blow split implicitly).

#### Regenerator Geometry

* **Tušek, Kitanovski, Zupan, Prebil & Poredoš**, *Applied Thermal Engineering*, **53**, 57–66 (2013).

  * Experimentally compared six AMR configurations including packed-bed and parallel-plate regenerators (three sub-devices, "AMR A/B/F", each swept across 3 flow-utilization ratios V*=0.16/0.42/0.95).
  * Achieved approximately **20 K temperature span** under a **1.15 T** magnetic field.
  * Reported that optimized parallel-plate regenerators produced the highest COP because of reduced flow resistance and improved heat transfer.
  * Demonstrated that magnetic refrigeration can outperform vapor-compression systems under certain operating conditions, although performance depends strongly on magnetic field strength and temperature span.
  * **CORRECTED — digitization gap now closed.** Earlier drafts of this review flagged Figs. 10–11 (the Qc(span) and COP(span) curves needed for a real curve-shape validation of `core/geometry_analysis.py`'s packed-bed/parallel-plate models) as undigitized. They have since been pixel-calibrated directly from the source PDF (`data/tusek_ate2013_figs/{fig10_data.csv,fig11_data.csv,tusek_ate2013_figs_notes.md}` — automated marker detection with by-eye series disambiguation at curve crossings; full method and residual-uncertainty notes in that folder's own notes file) and are used by `core/validation_system.py`'s curve-level check. Calibrating AMR(A) at V*=0.95 against its own anchor point (span=7.26 K, Qc=5.27 W) and predicting the intermediate point (span=12.23 K) gives a **+787.0%** error against the digitized literature value (2.03 W) — a large miss, reported directly (`main.py` step 2, "Tusek AMR(A) V*=0.95" lines) rather than treated as a clean validation.

* **Tušek, Kitanovski & Poredoš**, *International Journal of Refrigeration*, **36**, 1456–1464 (2013). (A different Tušek 2013 paper from the one above — numerical geometry optimization, not an experimental comparison.)

  * Numerically optimizes packed-bed sphere diameter and parallel-plate spacing/thickness against cooling capacity and COP for a fixed AMR outer envelope (Gd, water, 0-1T, 15K span). Reports its own operating-point optima in Table 3: packed-bed sphere diameter 0.07 mm (Qc-optimal)/0.17 mm (COP-optimal); parallel-plate spacing 0.035 mm (Qc-optimal)/0.075 mm (COP-optimal, independent of plate thickness).
  * Used in `core/thermal.py`/`core/geometry_analysis.py` to add a geometry-dependent hydraulic pumping-power term (packed-bed Ergun-type and parallel-plate laminar-channel relations, this paper's Eqs. 5–7) coupled to the NTU effectiveness gain — closing a gap where the pre-existing regenerator model had no coupling between particle/plate geometry and pumping cost at all, and so could not show an interior COP optimum vs. geometry even in principle. Re-run at this repo's own fixed representative flow (0.08 kg/s, since Section 9's own finding is that free COP-only optimization over ṁ is itself degenerate — see the "From Sensitivity to Optimization" entry below), the coupled model now *does* reproduce a genuine interior optimum: packed-bed sphere diameter 0.5 mm and parallel-plate spacing 0.1 mm both maximize COP_aug in `core/geometry_analysis.py`'s own sweep (`main.py` step 3c). These specific values are not expected to match the paper's own Table 3 exactly — the two studies differ in operating point, envelope geometry, and (crucially) whether ṁ is re-optimized per geometry — but the qualitative finding (an interior optimum now exists where none could before) is the confirmed result.

#### Hydraulic / Pumping-Power Architecture **[NEW]**

* **Klinar, Muhič, Tušek & Nielsen**, *Advanced Energy Materials*, **14**, 2401739 (2024).

  * Review covering, among other topics, "Hypereg"-style parallel-hydraulic AMR architectures — splitting a single series regenerator bed into several parallel sub-regenerator channels to reduce Darcy-flow pumping-power loss for the same total mass flow, illustrated with a worked example (the review's own Figs. 18–21).
  * Used directly in `core/hypereg_analysis.py` and `core/thermal.py`'s `pumping_power_packed_bed_hypereg()`. At this repo's own representative operating point (5 kg regenerator, mdot=0.08 kg/s, f=1 Hz), splitting into n=16 parallel sub-beds raises COP_electrical only modestly (5.264→5.278, n=1→n=16) because pumping loss is one of three loss channels and not the dominant one there; at a 4× higher flow rate (0.3 kg/s), the same 16-way split gives a 2.45% relative COP gain vs. 0.27% at baseline — consistent with the paper's own framing that the benefit scales with how pumping-dominated the operating point already is (`main.py` step 3d, `results/hypereg_findings.md`). This confirms the review's qualitative mechanism is representable in this repo's own model without claiming a validated optimum n or a device-level performance prediction for an as-yet-unbuilt concept.

#### Thermal Diodes / Rectified Heat Switching **[NEW]**

* **Kitanovski, Tušek, Tomc, Plaznik, Ožbolt & Poredoš**, *Magnetocaloric Energy Conversion*, Springer (2015), Ch. 6 ("Special Heat Transfer Techniques").

  * The plan-of-record source for mechanical-contact active thermal diode design in AMR devices. **Book-access honesty flag**: this project's own copy of Kitanovski et al. (2015) is a ~30-page front-matter/Ch. 1 excerpt only — Ch. 6 (pp. 211–268) is not present. `core/thermal_diode.py`'s `MechanicalContactDiode` therefore ships with illustrative, round-number forward/reverse-conductance and actuation-energy placeholders rather than values derived from this book, and the module is scoped down accordingly (a parasitic-cost-only sensitivity study, not the fuller NSGA-III categorical design variable originally planned).

* **Bywaters & Griffin** (cryogenic piezo-actuated gas-gap heat-switch characterization; used to ground the diode's rectification ratio).

  * Provides an independent, non-AMR-specific data point for mechanical-contact heat-switch forward/reverse conductance ratios, used to anchor `DEFAULT_MECHANICAL_CONTACT_DIODE`'s `rectification_ratio` after the Kitanovski Ch. 6 pages proved inaccessible — a substitute grounding, not a room-temperature-AMR-specific measurement, and flagged as such in the module.
  * **Finding, held to precisely:** `core/thermal_diode_analysis.py` first checked the plan's own premise directly — does this repo's model even have a mechanical-switching frequency ceiling for a diode to relax? It does not: `AMRSystem` has no internal frequency cap on `cooling_capacity()`/`magnetic_work()` (frequency enters only monotonic, uncapped loss terms); the only frequency bound anywhere is `core/optimize.py`'s unexplained 5.0 Hz NSGA-III search-space bound, not a physical constraint tied to mechanical valve switching in any comment or roadmap entry. The illustrative diode actuation cost then reduces COP_electrical by at most 0.03% across 0.5–8 Hz — a small, cost-only accounting (no offsetting heat-transfer benefit from `rectification_ratio` is modeled, since no closed-form relation for how rectification ratio improves AMR cycle performance was available to digitize). No benchmark device in this repo's corpus uses thermal diodes, so this module is a design-exploration tool, not a validated feature (`main.py` step 11c).

#### Permanent Magnet Design

* **Halbach**, K., *Nuclear Instruments and Methods*, **169**, 1 (1980).

  * The original closed-form treatment of the ideal cylindrical Halbach-array field, `B = Br·ln(Ro/Ri)`, used directly by `core/magnet_geometry.py`'s `halbach_bore_field_T()` for an idealized 2D Halbach cylinder (infinitely long, remanence Br, inner/outer radii Ri/Ro).

* **Bjørk, Bahl, Smith, Christensen & Pryds**, *International Journal of Refrigeration*, **34**, 1805–1816 (2011). ("Determining the minimum mass and cost of a magnetic refrigerator")

  * The primary source behind `core/economics.py`'s `COST_MCM_PER_KG` ($20/kg) and `COST_MAGNET_PER_KG` ($40/kg, NdFeB N42) constants and the magnet-to-MCM mass-ratio approximation used in the flat-ratio (pre-Phase-19) cost model.
  * CONFIRMED FULLY MINED: checked the full text, not just the abstract. The only content not already reflected in `economics.py`'s constants is Fig. 9 (minimum system cost vs. operating frequency, at fixed 20K span/100W) — a qualitative trend ("increasing frequency reduces cost," no simple optimum found) presented only as a figure, not a digitizable value, and it's the same magnet+MCM-only cost model already used here.

* **Bjørk, Bahl & Nielsen**, "The Halbach Cylinder — Design and Physical Limitations," arXiv:1410.1987. **[CITATION CORRECTED]**

  * An earlier draft of this review cited "Bjørk et al., arXiv:1410.1987" as the source for the "fields near 2 T give an effective cost-vs-performance compromise" claim under the heading *Permanent Magnet Design*, alongside language ("investigated optimized Halbach permanent magnet assemblies… demonstrated that fields near 2 T provide an effective compromise") that describes a field-vs-magnet-mass/design-limits paper, not the field-vs-*cost* claim it was attached to. Checking the arXiv identifier directly (`core/magnet_geometry.py`'s own module docstring, HONESTY FLAG #2) confirms this is a genuinely different paper by an overlapping author list from the Bjørk, Bahl, Smith, Christensen & Pryds (2011) cost paper above — this one covers Halbach-cylinder field-vs-geometry physical limits, not a cost-minimization study. The ~2 T "sweet spot" claim as originally paraphrased here remains sourced only to this (now correctly attributed) reference; it is **not** independently re-derived by this repo's own models. `core/magnet_geometry.py`'s Phase 19 fixed-MCM-mass, dollars-per-Kelvin proxy in fact finds its *own* cost-per-K minimum at 0.5 T, not ~2 T — explicitly *not* claimed as a refutation of the literature figure, since that simple proxy holds MCM mass fixed while sweeping only field, ignoring the system-level cooling-power/device-size trade-offs a real design would also make (plausibly the actual driver of the literature's own reported optimum). What the new Halbach-cylinder relation *does* confirm directly (not merely assert from closed-form algebra) is that magnet mass grows super-linearly with field — 1.69× at 1.0 T rising to 13.98× at 3.0 T relative to the old flat per-Tesla ratio across `core/optimize.py`'s own [1.0, 3.0] T search bounds — the specific nonlinearity this review's earlier language implied without the code actually having it.

#### Thermodynamic and Exergy Analysis

* **Kitanovski et al.**, *Magnetocaloric Energy Conversion*, Springer (2015).

  * Comprehensive reference covering AMR thermodynamics, second-law analysis, exergy efficiency, heat transfer, permanent magnet design, numerical modelling, and system optimization.
  * **Consolidated book-access note** (see Section 6 for the full accounting): this project's own copy is a ~30-page front-matter/Ch. 1 excerpt. Chapters/sections cited by the phase plan but not present in this copy, each separately confirmed rather than assumed: Ch. 3 (Magnetic Field Sources, incl. Sect. 3.4/3.5 on Halbach cylinders — Phase 19), Ch. 4 (AMR performance), Ch. 6 (thermal diodes/heat switches — Phase 18), Ch. 7 (prototypes by country), Ch. 9 (system costs), and the fluids/ferrofluid chapter referenced by the Phase 20 plan (Phase 20).

#### Recent Developments

* **Recent review papers (2020 onwards).**

  * Summarize advances in magnetocaloric materials, regenerator geometries, permanent magnet systems, prototype refrigerators, modelling techniques, and performance metrics.

---

## 3. Data Center Cooling

#### Thermal Operating Guidelines

* **ASHRAE TC9.9**, *Thermal Guidelines for Data Processing Environments*, 5th Edition (2021).

  * Defines recommended inlet air temperatures for Class A1 and A2 data centers.
  * Establishes operating conditions commonly used when evaluating alternative cooling technologies. This repo's representative operating point (T_cold=18°C/291 K) and 5–20 K span sweep are anchored to this envelope throughout `main.py`.

#### Liquid Cooling Standards

* **ASHRAE TC9.9**, *Liquid Cooling Guidelines for Datacom Equipment Centers*, 2nd Edition (2021).

  * Provides recommended facility-water temperature ranges, liquid cooling architectures and operational guidelines for modern high-density servers.

#### Conventional Cooling Technologies

* **Ebrahimi, Jones & Fleischer**, *Renewable and Sustainable Energy Reviews*, **31**, 622–638 (2014).

  * Reviews air cooling, liquid cooling, chilled-water systems and economizers used in data centers.
  * Includes energy efficiency and second-law performance comparisons that serve as useful baselines for emerging cooling technologies, and directly informs `core/baseline_cooling.py`'s liquid-cooling free-cooling/mechanical-assist blend (`f_econ` default 0.50–0.55).

#### Chip-Level Thermal Management

* **Shah, Bash & Patel** (2004).

  * Discuss cooling requirements and power density trends for high-performance processors.
  * Highlights the growing need for efficient cooling solutions capable of supporting increasing rack power densities.

#### Elastocaloric Cooling — Comparison Reference **[NEW, Phase 23]**

* **Qian, Catalini, Muehlbauer, Liu, Mevada, Hou, Hwang, Radermacher & Takeuchi**, *Science*, **380**, 722–727 (2023).

  * Reports a simulated steady-state multimode elastocaloric cooling system with system-level COP=5.8, up to a 22.5 K span.
  * Used as the *high* end of `core/baseline_cooling.py`'s `elastocaloric_reference_cop()` (`ELASTOCALORIC_COP_HIGH = 5.8`).

* **Wu et al.**, "Continuous and efficient elastocaloric air cooling by coil-bending," *Nature Communications*, **14**, 7982 (2023).

  * Reports a *measured* device-level COP=3.7 at a much narrower (~1 K) span.
  * Used as the *low* end (`ELASTOCALORIC_COP_LOW = 3.7`).
  * **Honesty flag carried in the module and reproduced here**: the elastocaloric comparison line in `results/comparison_table.csv`/fig08 is a **flat, static literature anchor (COP=4.63, midpoint of the 3.7–5.8 range) repeated across every span** — not a span-dependent simulation the way the AMR/VCC/liquid-cooling rows are. Neither of this repo's two source books (Kitanovski et al. 2015; Tishin & Spichkin 2003) covers elastocalorics at all, so this is a comparison row sourced entirely from these two external papers, added exactly as scoped ("a comparison row, not a new simulated device").

---

## 4. Techno-Economic Analysis

#### Magnetic Refrigeration Cost

* **Bjørk, Bahl & Smith**, *International Journal of Refrigeration*, **34**, 1805–1816 (2011).

  * See the full entry and correction note under Section 2's *Permanent Magnet Design* — this is the primary source behind `economics.py`'s per-kg cost constants and, as of Phase 15, its `bom_cost()` full-system BOM extension (soft-magnetic yoke term added per Silva et al. 2017's methodology, below).

* **Bahl, Engelbrecht et al.**, *International Journal of Refrigeration*, **37**, 78–83 (2014).

  * Presents cost breakdowns for AMR systems.
  * Identifies permanent magnets and magnetocaloric materials as the dominant contributors to capital cost. Used as the qualitative basis for treating the materials-only BOM figure in `core/economics.py` as a cost *floor*, not a full-system estimate.

#### Full-System Cost Estimation **[NEW, Phase 15]**

* **Russek & Zimm**, (manufactured vapor-compression-AC cost benchmark; NIST/Astronautics collaboration literature).

  * Used as the source for the order-of-magnitude "materials BOM × 10×" full-system cost multiplier in `core/economics.py`'s `full_system_cost_estimate()` — an explicitly order-of-magnitude-only estimate, not a bottom-up HX/pump/motor/controls BOM, which remains an open item.

* **Silva, Paiva, Dutra et al.** (2017) (CRF-based levelized-cost-of-cooling methodology).

  * The methodology basis for `core/economics.py`'s `levelized_cost_of_cooling()` (capital-recovery-factor annualization, 15-yr life, 6% discount rate) — a second, independent cost accounting alongside the CAPEX/OPEX table, giving $0.0341/kWh_cooling at the representative operating point ($0.0126 capital + $0.0216 electricity, materials-only capital basis).

#### Magnetocaloric Material Economics

* **Franco, Blázquez et al.**, *International Journal of Refrigeration*, **57**, 288–298 (2018).

  * Reviews manufacturing processes, rare-earth availability, material cost and commercialization challenges for magnetocaloric materials. An earlier, less-grounded $175/kg Gd-cost placeholder in this repo's cost objective traced to this review before being replaced by the Bjørk et al. (2011)-grounded figures above.

#### Data Center Cooling Economics

* **Lawrence Berkeley National Laboratory (LBNL)** publications.

  * Provide benchmarks for cooling energy consumption, operating cost and energy-efficiency metrics in modern data centers.
  * Offer reference values for comparing the techno-economic feasibility of magnetic refrigeration with conventional cooling systems.

---

## 5. Alternative Working-Body Architectures **[NEW]**

Two exploratory directions were added late in this project's timeline as
*design-exploration* studies rather than validated results — no benchmark
device exists in this repo's corpus for either, a limitation stated
directly rather than papered over.

#### Magnetocaloric Fluids (Phase 20)

* Standard rheology and pipe-flow relations were used in place of a
  digitized fluids-chapter source, since neither Kitanovski et al. (2015)'s
  nor Tishin & Spichkin (2003)'s fluids/ferrofluid content is present in
  this repo's copies (checked directly for both, not assumed):
  the **Krieger–Dougherty** viscosity relation for the phi-dependent
  suspension viscosity, and **Darcy–Weisbach** pipe-flow pressure drop for
  pumping power. `core/fluid_mce_cycle.py`'s `FerrofluidMCESystem` combines
  these with a mixture-heat-capacity dilution model for suspension ΔTₐd.
  `core/fluid_mce_analysis.py`'s volume-fraction sweep finds a genuine
  interior COP_electrical optimum near φ≈0.10 (108.8), but the headline
  finding is architectural, not a COP number: the mixture-dilution model
  combined with this architecture's lack of a regenerator collapses usable
  span to under 1 K at realistic loadings (0.68 K at φ=0.20) — dramatically
  narrower than solid AMR achieves at the same field/flow, whose
  regenerator bed amplifies achievable span well beyond a single stage's
  own ΔTₐd (`main.py` step 14, `results/fluid_mce_analysis.txt`).

#### Passive / Hybrid Magnetic Regeneration (Phase 21)

* **Tishin & Spichkin**, *The Magnetocaloric Effect and its Applications*, IOP Publishing (2003), Ch. 11.

  * The plan-of-record source for passive/hybrid regenerator augmentation of a conventional gas cycle — recombining an MCE material's Curie-point heat-capacity anomaly with a conventional (vapor-compression) regenerator's own effectiveness. **Book-access honesty flag, confirmed directly**: this project's copy of Tishin & Spichkin (2003) is a 486-page scanned, image-only PDF with zero pages returning extractable text via `pdfplumber` — Ch. 11 could not be digitized. `core/baseline_cooling.py`'s `augmented_regenerator_cop()` therefore recombines this repo's *own already-existing* heat-capacity and baseline-COP models instead, with an illustrative, literature-range-anchored effectiveness-to-COP ceiling (capped at +8%) rather than a fitted or digitized coefficient.
  * **Finding, confirmed rather than assumed:** at the representative ASHRAE point, Gd (Curie point 294 K, inside the [291.1, 301.1] K operating window) gives the largest boost — regenerator effectiveness 0.829→0.860, base vapor-compression COP 12.228→12.258 (+0.24%) — while La₀.₇Ca₀.₃MnO₃ (Curie point 267 K, far outside the window) gives +0.00% gain, and this alignment-vs.-mismatch pattern holds across a 5–20 K span sweep, with the gain shrinking as span widens (`main.py` step 15, `results/passive_regenerator_analysis.txt`).

---

## 6. Materials-Property Compendia — Consolidated Book-Access Accounting **[NEW]**

Two reference books recur across the phase plan as the intended primary
source for several late-stage items. Both were checked directly for
extractability rather than assumed accessible, and both turned out to be
only partially usable:

* **Kitanovski, Tušek, Tomc, Plaznik, Ožbolt & Poredoš**, *Magnetocaloric
  Energy Conversion*, Springer (2015). This project's copy is front matter
  plus ~9 pages of Ch. 1 only. Sections cited by the phase plan but absent
  from this copy: Ch. 3 (incl. Sect. 3.4/3.5, Halbach cylinders — Phase
  19), Ch. 4 (AMR performance, Sect. 4.1.1–4.1.4 cycle-topology relations
  — Phase 17), Ch. 6 (special heat-transfer techniques / thermal diodes —
  Phase 18), Ch. 7 (prototypes by country), Ch. 9 (system costs), and its
  fluids-chapter content (Phase 20).

* **Tishin & Spichkin**, *The Magnetocaloric Effect and its Applications*,
  IOP Publishing (2003). 486 pages, entirely scanned page images with no
  OCR text layer — confirmed via `pdfplumber` (zero pages return
  extractable text across every page sampled), not merely assumed from a
  quick check. Likely the deepest materials-property compendium in the
  corpus, but full-book OCR was not attempted given the effort/yield
  trade-off against primary sources already mined. Sections cited by the
  phase plan but inaccessible here: Sect. 2.8 (inhomogeneous ferromagnets
  / Tc-distribution broadening — Phase 22 item 1), Sect. 2.9 and Ch. 10
  (superparamagnetic / nanocomposite materials — Phase 22 item 2), Ch. 9
  (amorphous / melt-spun-ribbon materials — Phase 22 item 3), and Ch. 11
  (regenerator/hybrid-cycle content — Phase 21).

Where an item's plan-of-record source turned out to be inaccessible, the
corresponding module says so explicitly and substitutes either (a) the
standard textbook/literature treatment of the same physical mechanism
(inhomogeneous broadening: standard Gaussian-Tc-ensemble treatment; fluids:
Krieger–Dougherty/Darcy–Weisbach), or (b) a qualitative-only note with no
invented numeric placeholder where no defensible number exists anywhere in
this repo's corpus (Phase 22 item 3's amorphous-material cost/performance
note is the clearest example: it is *not* wired into any priced dict,
because doing so would mean inventing a $/kg or ΔS_M figure this repo has
no basis for). If specific data tables are wanted from either book,
targeted OCR on a stated page range/topic is a cheaper path than either
the whole book or waiting for a more complete Kitanovski copy.

---

## 7. Research Gaps

The literature identifies several areas requiring further investigation.
Some gaps identified in earlier drafts of this review have since been
substantially (though not completely) addressed by this project's own
work, noted inline below.

* Most experimental magnetic refrigeration systems have been developed for
  domestic refrigeration rather than continuous, high-load data-center
  cooling. *(Still true of the underlying literature; this repo's own
  16-device-row benchmark set spans domestic-, lab- and a few
  larger-scale prototypes, none purpose-built for data-center duty.)*

* System-level validation data for AMRs operating within the **20–40°C**
  temperature range relevant to data centers remain limited. *(Partially
  addressed: this repo's benchmark set now spans 16 device rows across
  ~13 device groups, but the largest and most data-center-relevant devices
  — Astronautics, DTU MagQueen, Cooltech, the Risø/DTU 30 K-span device —
  are exactly the ones this repo's own model could NOT calibrate against
  within a physically plausible flow-rate range, an unresolved gap
  documented rather than hidden — see `main.py` step 2's "NO CALIBRATION
  FOUND" rows.)*

* Existing studies often validate either **material properties** or
  **overall system performance**, while comparatively few provide
  comprehensive validation across both levels. *(This repo's own two-tier
  validation strategy — Dan'kov/Giguère at the material level,
  `amr_experimental_benchmarks.csv` at the system level — directly targets
  this gap; see the Validation Strategy section below.)*

* Although first-order magnetocaloric materials exhibit larger entropy
  changes, many suffer from hysteresis, limited operating temperature
  windows and material stability issues, reducing their suitability for
  data-center applications. *(Now partially quantified rather than only
  asserted: Phase 16's literature-analog hysteresis-loss placeholders,
  heavily honesty-flagged, shift the NSGA-III material-selection front's
  composition — though a full-production multiseed check subsequently
  found the original reduced-setting reversal is NOT stable, an open item
  rather than a settled finding.)*

* Comparisons between magnetic refrigeration and modern liquid-cooling
  technologies remain limited, particularly under realistic server
  operating conditions and facility water temperatures. *(Addressed by
  this repo's own `main.py` baseline sweep across the ASHRAE 5–20 K span
  range — AMR trails both vapor-compression and liquid cooling on
  electrical COP throughout, with the gap widening at wider spans, per
  `results/comparison_table.csv`.)*

* Few studies integrate **material selection, AMR thermodynamics,
  permanent magnet optimization and techno-economic analysis** into a
  unified framework for evaluating data-center cooling applications.
  *(Now largely addressed by this repo's own Phase 15 NSGA-III
  material+geometry co-optimization, extended in Phase 19 with a
  physically nonlinear magnet-mass-vs-field cost term — though the
  underlying loss-model calibration remains a small, exactly-determined
  3-point fit, so the *specific* optimum found should be read as
  illustrative of the methodology, not a converged final design.)*

* **New gaps opened by this project's own late-phase exploration, not
  present in the original literature-gap list:** no benchmark device
  exists anywhere in this repo's corpus for either magnetocaloric-fluid
  working bodies (Phase 20) or mechanical-contact thermal diodes (Phase
  18) — both remain design-exploration tools rather than validated
  features, and closing that gap would require either a purpose-built
  literature search beyond this project's own passes or new experimental
  data neither source book supplies.

* **Reference books in this corpus remain largely untapped** — see
  Section 6 above for the full, per-section accounting of what Kitanovski
  et al. (2015) and Tishin & Spichkin (2003) could and could not supply.

---

### Validation Strategy for This Work

The numerical model developed in this study is validated at two
complementary levels:

1. **Material-level validation** using the experimental magnetocaloric
   data of **Dan'kov et al. (1998)** for gadolinium (adiabatic temperature
   change and Curie-shift-with-field), extended to 7 T by **Giguère et
   al. (1999)**'s independent cross-check, and separately for the
   first-order Landau model of Gd₅Si₂Ge₂ against Giguère et al.'s own
   **direct** ΔTₐd measurement and Pecharsky & Gschneidner (1997)'s
   heat-capacity-derived peak ratio (Sections 1 above). A Phase 22
   Gaussian Tc-broadening sensitivity check additionally probes how
   grain-to-grain inhomogeneity would move the material-level error
   pattern, using the standard literature treatment of the mechanism since
   Tishin & Spichkin's own Sect. 2.8 could not be digitized (Section 6).

2. **System-level validation** against an experimental benchmark set
   (`data/amr_experimental_benchmarks.csv`, checked in
   `core/validation_system.py`), which spans 12–13 device groups / 22
   rows: Jacobs et al.'s Astronautics rotary device (2014); the Risø/DTU
   rotary Gd device (Engelbrecht et al., Purdue ICR 2010/2016); the DTU
   "MAGGIE" rotary Gd device, represented by two separate papers on the
   same physical prototype — Eriksen et al., *Int. J. Refrigeration*
   (2015) and Eriksen/Engelbrecht/Bahl/Bjørk, *Sci. Technol. Built
   Environ.* 22(5) (2016); the Tušek single-bed Gd device (2010/2013,
   including the Phase-and-figure-digitized curve check above); Okamura &
   Hirano's device (2013, as reported secondhand in a Trevizoli & Barbosa
   review); the Lozano/POLO-UFSC rotary device (2016, 8 independent
   operating points); and the ChubuToshiba, Cooltech, and DTU MagQueen
   devices (all secondary-sourced, see the CSV's own notes) — comparing
   cooling capacity, COP, utilization factor and temperature span at each
   device's reported operating point, a curve-shape check (Tušek et al.
   2013's digitized Figs. 10–11) for the device that has one, and (Phase
   17) a directional cycle-topology (Ericsson- vs. Brayton-like)
   re-classification sensitivity check for rotary-named devices.

This two-stage validation approach improves confidence in both the
underlying magnetocaloric material model and the overall AMR system
simulation before comparison with conventional data-center cooling
technologies — while several honestly-reported gaps (large-device
calibration failures, an order-of-magnitude curve-shape miss, an
unvalidated loss model beyond 4 points) remain open rather than smoothed
over, consistent with this project's own stated preference for reporting
findings directly over reporting only favorable ones.