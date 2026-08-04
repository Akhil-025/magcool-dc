## Literature Review — Magnetocaloric Cooling for Data Centers

### 1. Magnetocaloric Effect and Magnetocaloric Materials

#### Giant Magnetocaloric Effect

* **Pecharsky & Gschneidner**, *Physical Review Letters*, **78**, 4494–4497 (1997).

  * First reported the **giant magnetocaloric effect (GMCE)** in **Gd₅Si₂Ge₂**.
  * Demonstrated that first-order magnetostructural phase transitions can produce substantially larger entropy and temperature changes than conventional second-order ferromagnets.
  * Established the foundation for modern room-temperature magnetic refrigeration research.

#### Gadolinium as the Benchmark Material

* **Pecharsky & Gschneidner**, *Journal of Magnetism and Magnetic Materials*, **200**, 44–56 (1999).

  * Comprehensive review of the magnetic, thermodynamic and magnetocaloric properties of gadolinium.
  * Identified Gd as the standard reference material for room-temperature magnetic refrigeration owing to its Curie temperature (~294 K), making it directly relevant to data-center operating temperatures.

#### Experimental Validation of Magnetocaloric Properties

* **Dan'kov, Tishin, Pecharsky & Gschneidner**, *Physical Review B*, **57**, 3478–3490 (1998).

  * Measured the adiabatic temperature change (ΔTₐd) of gadolinium under magnetic fields of **1 T, 2 T and 5 T** near its Curie temperature.
  * These measurements remain one of the most widely used experimental datasets for validating theoretical and numerical MCE models.
  * Also reports the Curie-point transition temperature increasing almost linearly with field, ~6 K/T above 2 T up to 7.5 T (zero-field TC = 294 ± 1 K, confirmed by four independent techniques). Used in Phase 11 (`core/validation.py`, `run_curie_shift_check()`) as a held-out check of the mean-field model's own emergent peak-shift-with-field — the model does NOT reproduce this rate (fitted ~0 K/T vs. the reported ~6 K/T), a genuine, documented limitation rather than a re-confirmed match.

* **Giguère, Foldeaki, Ravi Gopal, Chahine, Bose, Frydman & Barclay**, *Physical Review Letters*, **83**, 2262 (1999).

  * Primarily a direct-measurement Gd₅Si₂Ge₂ paper (already used by `core/giguere_validation.py`), but its methods section also independently cross-checks pure Gd: ΔTₐd ≈ 10.5–11.5 K at 5 T (this paper's high-purity sample vs. AMES laboratory's measurement, agreeing within 1 K) and ≈ 12–13 K at 7 T (industrial- vs. high-purity Gd respectively), agreeing with Brown (1976)'s independently reported 14 K.
  * Used in Phase 11 (`core/validation.py`, `run_giguere_gd_extension()`) as a second, independent Gd dataset extending validation to 7 T. The model (calibrated to Dan'kov et al.'s higher 5T value) overestimates relative to this range at both fields — a real disagreement between two published Gd measurements, reported rather than reconciled.

* **Pecharsky & Gschneidner**, *Physical Review Letters*, **78**, 4494 (1997).

  * The original Gd₅Si₂Ge₂ "giant" magnetocaloric effect discovery paper — already cited in this repo for `GD5SI2GE2_FIRST_ORDER`'s Tc and J, but its own text also reports a heat-capacity-derived comparison not previously extracted: "the ΔTₐd values of Gd₅Si₂Ge₂ are larger than the corresponding ΔTₐd values for Gd by about 30%, comparing the peak values, regardless of the temperature" (Fig. 6, field changes of 0→2T and 0→5T).
  * Used in Phase 12 (`core/giguere_validation.py`, `run_pecharsky_ratio_check()`) as a second, independent primary source (heat-capacity-based, not the pulse-field-thermometry method Giguère et al. used) and a second field range (2T/5T, vs. the single 7T point `DTAD_CORRECTION_FACTOR` is fit to). Result: the raw (uncorrected) model's ratio at 5T (~1.24) is close to this paper's ~1.30 — but applying the Giguère-derived correction overcorrects it to ~0.51, predicting Gd₅Si₂Ge₂ underperforms plain Gd. Documented as evidence the correction should not be treated as field-independent.

#### Magnetocaloric Theory

* **de Oliveira & von Ranke**, *Physics Reports*, **489**, 89–159 (2010).

  * Comprehensive review of thermodynamic models describing the magnetocaloric effect.
  * Highlights the limitations of mean-field theory near the Curie temperature due to critical magnetic fluctuations and discusses improved theoretical formulations.

#### Alternative Magnetocaloric Materials

* **Guo et al.**, *Applied Physics Letters*, **78**, 1142 (1997).
* **Phan & Yu**, *Journal of Magnetism and Magnetic Materials*, **308**, 325–340 (2007).

  * Investigated perovskite manganites and other rare-earth alloys exhibiting significant room-temperature magnetocaloric effects.
  * Discuss material selection criteria including entropy change, hysteresis, thermal conductivity and operating temperature.

#### (Mn,Fe)₂(P,Si) Giant Magnetocaloric Materials

* **Hanggai, Yibole, Guillou, Kwakernaak, van Dijk & Brück**, *Acta Materialia*, **302**, 121677 (2026).

  * Preparation and calorimetric analysis of Fe-rich melt-spun Mn₀.₆₀₊ₓFe₁.₃₋ₓP₀.₆₆₋ᵧSi₀.₃₄₊ᵧ (0≤x≤0.08, x=2y) ribbons.
  * Reports Curie temperature increasing linearly from 295.3 K (parent) to 331.2 K (highest Mn/Si tested) across five measured compositions — a directly-measured, room-temperature-adjacent tunability window.
  * Reports peak |ΔS_max| ~17.6 J/(kg K) at a 2T field change for the highest-Mn/Si composition, cross-validated by two independent methods (16.66 J/(kg K) calorimetric, 17.61 J/(kg K) magnetization) — a ~40% enhancement over the 12 J/(kg K) parent-compound value at the same field.
  * Used in Phase 10 (`core/first_order_mce.py`, `core/cascade.py`) to add `MNFEPSI_FAMILY`, a third pluggable Curie-graded material family alongside the Gd₅(SixGe1-x)₄ and La(Fe,Si)₁₃Hy families already in this repo — notable because its Tc window sits almost entirely at or above the ASHRAE data-center supply range, the opposite tension from the Gd₅(SixGe1-x)₄ family's ceiling sitting just below it.

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
  * Used in Phase 11 (`data/amr_experimental_benchmarks.csv`) to add the Chubu Electric/Toshiba two-field-point row (`ChubuToshiba_Gd_2016_4T`/`_2T`, original device ref [69] in this review, not itself in this repo's `Papers/`) — SECONDARY SOURCE, same caveat as the existing `Okamura_Hirano_2013` row. Several other rows in the same table (Institute of Tech. Chubu's 540W near-zero-span extreme; a second, independent Riso Lab data point with different regenerator count/frequency than the existing `Risoe_DTU_Gd_2011` row; Teyber et al.'s 100K-span superconducting-magnet outlier) were identified but left un-added — flagged in `ROADMAP.md` Phase 11 rather than built, pending either a clearer scope decision or the primary source papers.

* **Greco, Aprea, Maiorino & Masselli**, *International Journal of Refrigeration* (2019).

  * "A review of the state of the art of solid-state caloric cooling processes at room-temperature before 2019" — Table 2 lists ~25 AMR prototypes built 2009–2018 (Q̇_ref,max, ΔT_span,max, field, MCM/HTF), larger and cleaner than the Kamran/Ahmad/Wang table above, though the PDF carries a diagonal "ACCEPTED MANUSCRIPT" watermark that bleeds letter fragments into extracted text on several rows — checked against surrounding column structure, but not a rendered-image read like the Kamran table.
  * Used in Phase 12 (`data/amr_experimental_benchmarks.csv`) to add `Cooltech_2013_rotary` (42K span — the largest in this benchmark set — confirmed directly in the review's own body text, not table-only) and `DTU_MagQueen_2018` (a La(Fe,Mn,Si)₁₃Hz-based heat pump, giving this repo's first LAFESIH-material benchmark point independent of Astronautics; Qc/COP values are DERIVED from the paper's own reported heating power/COP, flagged in the CSV row). The table's Astronautics 2014 entry lists ΔT_span,max=18K, unconfirmed elsewhere in this corpus and NOT added as a third Astronautics data point — flagged in `ROADMAP.md` Phase 12 as worth checking against the primary Jacobs et al. paper before use.

#### Performance Modelling

* **Tušek, Zupan, Šarlah, Prebil & Poredoš**, *International Journal of Refrigeration*, **33**, 294–300 (2010).

  * CORRECTED (Paper-Mining Pass Part 3, §3): this citation previously shared a "characteristic-curve models" description with Nielsen et al. (2011) below — confirmed by identifying the underlying PDF (`Development of a rotary magnetic refrigerator.pdf`) that this is actually a mechanical/magnet-design paper for an earlier Ljubljana AMR prototype, NOT a performance-testing or characteristic-curve paper. It reports zero Qc/span numbers.
  * Content is the permanent-magnet assembly design and a "pros and cons" table of the mechanical build (shaft-seal leakage causing rotational friction/heat generation, large magnet-structure weight, assembly/disassembly complexity).
  * Useful only as qualitative engineering-realism context — e.g. for discussing real mechanical parasitic losses that `core/amr_cycle.py`'s idealized 0-D cycle doesn't capture — not as a numeric validation source. (Not to be confused with `Development of a novel rotary magnetic refrigerator.pdf` = Lozano et al. 2016, already the primary source behind the `Lozano_POLO_UFSC_2016` CSV rows.)

* **Nielsen et al.**, *International Journal of Refrigeration*, **34**, 603–616 (2011).

  * Developed characteristic-curve models relating cooling capacity, utilization factor, operating frequency, pressure drop and regenerator effectiveness.
  * These models provide widely accepted benchmarks for predicting AMR performance and validating numerical simulations.

#### Flow-Profile / Blow-Fraction Effects

* **Masche, Liang, Engelbrecht & Bahl**, *Applied Thermal Engineering*, **215**, 118945 (2022).

  * Experimentally varied flow profiles (blow fraction — the fraction of the cycle period spent in cold-to-hot vs. hot-to-cold flow) on a DTU rotary AMR device (13 trapezoidal beds, 295 g Gd spheres/bed, solenoid-valve-controlled blow fraction).
  * At a fixed 16 K span/U=0.32/1.4 Hz operating point, increasing blow fraction from 25.0% to 41.6% raised cooling capacity from 70 W to 330 W (~4.7×) and second-law efficiency from 2.6% to 17.4%; lower blow fractions instead favored faster temperature pulldown.
  * Used in Phase 10 (`core/amr_cycle.py`) to add `blow_fraction` as a new AMR cycle degree of freedom (`BLOW_FRACTION_MASCHE`, `_blow_fraction_multiplier()`) and as a 6th NSGA-III design variable in `core/optimize.py` — a real, sizeable model gap this repo had no prior representation of at all (Qc/second-law efficiency previously assumed a fixed symmetric 50/50 blow split implicitly).

#### Regenerator Geometry

* **Tušek, Kitanovski, Zupan, Prebil & Poredoš**, *Applied Thermal Engineering*, **53**, 57–66 (2013).

  * Experimentally compared six AMR configurations including packed-bed and parallel-plate regenerators.
  * Achieved approximately **20 K temperature span** under a **1.15 T** magnetic field.
  * Reported that optimized parallel-plate regenerators produced the highest COP because of reduced flow resistance and improved heat transfer.
  * Demonstrated that magnetic refrigeration can outperform vapor-compression systems under certain operating conditions, although performance depends strongly on magnetic field strength and temperature span.
  * **Open validation gap (flagged in ROADMAP.md Phase 10):** this is the only paper in the corpus that could validate `core/geometry_analysis.py`'s `regenerator_effectiveness_parallel_plate()` model against a real device — every row in `data/amr_experimental_benchmarks.csv` is currently packed-bed or layered packed-bed. The exact (Qc, COP) pair needed sits only in Figs. 10–11, which remain undigitized (Phase 7 open item).

* **Tušek, Kitanovski & Poredoš**, *International Journal of Refrigeration*, **36**, 1456–1464 (2013). (A different Tušek 2013 paper from the one above — numerical geometry optimization, not an experimental comparison.)

  * Numerically optimizes packed-bed sphere diameter and parallel-plate spacing/thickness against cooling capacity and COP for a fixed AMR outer envelope (Gd, water, 0-1T, 15K span).
  * Reports genuine trade-off optima driven by the competing effects of heat-transfer coefficient and viscous pressure drop: packed-bed sphere diameter 0.07mm (Qc-optimal)/0.17mm (COP-optimal); parallel-plate spacing 0.035mm/0.075mm.
  * Used in Phase 7 (`core/thermal.py`, `core/geometry_analysis.py`) to add a geometry-dependent pumping-power term and a parallel-plate effectiveness model, closing a gap where this repo's regenerator model had no coupling between particle/plate geometry and pumping cost at all.

#### Thermodynamic and Exergy Analysis

* **Kitanovski et al.**, *Magnetocaloric Energy Conversion*, Springer (2015).

  * Comprehensive reference covering:

    * AMR thermodynamics
    * Second-law analysis
    * Exergy efficiency
    * Heat transfer
    * Permanent magnet design
    * Numerical modelling
    * System optimization

#### Experimental AMR Performance

* **Experimental packed-bed AMR studies** (2016).

  * Reported maximum second-law efficiency at temperature spans of approximately **15–20 K**.
  * These operating conditions are representative of many electronics and liquid-cooling applications.

#### Permanent Magnet Design

* **Bjørk et al.**, arXiv:1410.1987.

  * Investigated optimized Halbach permanent magnet assemblies.
  * Quantified trade-offs among magnetic field strength, magnet mass and cooling performance.
  * Demonstrated that magnetic fields near **2 T** provide an effective compromise between performance and permanent magnet cost.

#### Recent Developments

* **Recent review papers (2020 onwards).**

  * Summarize advances in:

    * magnetocaloric materials
    * regenerator geometries
    * permanent magnet systems
    * prototype refrigerators
    * modelling techniques
    * performance metrics

---

## 3. Data Center Cooling

#### Thermal Operating Guidelines

* **ASHRAE TC9.9**, *Thermal Guidelines for Data Processing Environments*, 5th Edition (2021).

  * Defines recommended inlet air temperatures for Class A1 and A2 data centers.
  * Establishes operating conditions commonly used when evaluating alternative cooling technologies.

#### Liquid Cooling Standards

* **ASHRAE TC9.9**, *Liquid Cooling Guidelines for Datacom Equipment Centers*, 2nd Edition (2021).

  * Provides recommended facility-water temperature ranges, liquid cooling architectures and operational guidelines for modern high-density servers.

#### Conventional Cooling Technologies

* **Ebrahimi, Jones & Fleischer**, *Renewable and Sustainable Energy Reviews*, **31**, 622–638 (2014).

  * Reviews air cooling, liquid cooling, chilled-water systems and economizers used in data centers.
  * Includes energy efficiency and second-law performance comparisons that serve as useful baselines for emerging cooling technologies.

#### Chip-Level Thermal Management

* **Shah, Bash & Patel** (2004).

  * Discuss cooling requirements and power density trends for high-performance processors.
  * Highlights the growing need for efficient cooling solutions capable of supporting increasing rack power densities.

---

## 4. Techno-Economic Analysis

#### Magnetic Refrigeration Cost

* **Bjørk, Bahl & Smith**, *International Journal of Refrigeration*, **34**, 1805–1816 (2011).

  * "Determining the minimum mass and cost of a magnetic refrigerator" — the primary source behind `core/economics.py`'s `COST_MCM_PER_KG` ($20/kg) and `COST_MAGNET_PER_KG` ($40/kg, NdFeB N42) constants, and the magnet-to-MCM mass ratio fit.
  * CONFIRMED FULLY MINED (Paper-Mining Pass Part 3, §4): checked the full text, not just the abstract already used. The only content not already reflected in `economics.py`'s constants is Fig. 9 (minimum system cost vs. operating frequency, at fixed 20K span/100W) — a qualitative trend ("increasing frequency reduces cost," no simple optimum found) presented only as a figure, not a digitizable value, and it's the same magnet+MCM-only cost model already used here. Given `ROADMAP.md` already treats full BOM cost (HX/pump/motor) as open pending real data, digitizing this figure wouldn't move that needle. No action taken.

* **Bahl, Engelbrecht et al.**, *International Journal of Refrigeration*, **37**, 78–83 (2014).

  * Presents cost breakdowns for AMR systems.
  * Identifies permanent magnets and magnetocaloric materials as the dominant contributors to capital cost.

#### Magnetocaloric Material Economics

* **Franco, Blázquez et al.**, *International Journal of Refrigeration*, **57**, 288–298 (2018).

  * Reviews manufacturing processes, rare-earth availability, material cost and commercialization challenges for magnetocaloric materials.

#### Data Center Cooling Economics

* **Lawrence Berkeley National Laboratory (LBNL)** publications.

  * Provide benchmarks for cooling energy consumption, operating cost and energy-efficiency metrics in modern data centers.
  * Offer reference values for comparing the techno-economic feasibility of magnetic refrigeration with conventional cooling systems.

---

## 5. Research Gaps

The literature identifies several areas requiring further investigation:

* Most experimental magnetic refrigeration systems have been developed for domestic refrigeration rather than continuous, high-load data-center cooling.

* System-level validation data for AMRs operating within the **20–40°C** temperature range relevant to data centers remain limited.

* Existing studies often validate either **material properties** or **overall system performance**, while comparatively few provide comprehensive validation across both levels.

* Although first-order magnetocaloric materials exhibit larger entropy changes, many suffer from hysteresis, limited operating temperature windows and material stability issues, reducing their suitability for data-center applications.

* Comparisons between magnetic refrigeration and modern liquid-cooling technologies remain limited, particularly under realistic server operating conditions and facility water temperatures.

* Few studies integrate **material selection, AMR thermodynamics, permanent magnet optimization and techno-economic analysis** into a unified framework for evaluating data-center cooling applications.

* **Reference books in this corpus remain largely untapped (Paper-Mining Pass Part 3, §5)**: Kitanovski et al., *Magnetocaloric Energy Conversion* (2015) — the corpus copy is front-matter + ~9 pages of Ch.1 only, NOT the full book (Chapters 4/7/9 on AMR performance, prototypes-by-country and costs are listed in the table of contents but not present in this file). Tishin & Spichkin, *The Magnetocaloric Effect and its Applications* (2003) — 486 pages, but entirely scanned page images with no OCR text layer (confirmed via `pdfplumber`, zero pages return extractable text); likely the deepest materials-property compendium in the corpus, but OCR'ing all 486 pages wasn't attempted given the effort/yield tradeoff against primary sources already mined. If specific data tables are wanted from either book, targeted OCR on a stated page range/topic is a cheaper path than either the whole book or waiting for a more complete Kitanovski copy.

---

### Validation Strategy for This Work (Optional Section)

The numerical model developed in this study is validated at two complementary levels:

1. **Material-level validation** using the experimental magnetocaloric data of **Dan'kov et al. (1998)** for gadolinium, including adiabatic temperature change and entropy variation.

2. **System-level validation** against a five-device experimental benchmark set (`data/amr_experimental_benchmarks.csv`, checked in `core/validation_system.py`): Jacobs et al.'s Astronautics rotary device (2014), the DTU/Risø rotary Gd device (Bahl/Eriksen/Engelbrecht, 2016; Engelbrecht et al., Purdue ICR 2010/2016), the Tušek single-bed Gd device (2010), Okamura & Hirano's device (2013, as reported secondhand in a Trevizoli & Barbosa review), and the Lozano/POLO-UFSC rotary device (2016) — comparing cooling capacity, COP, utilization factor and temperature span at each device's reported operating point, plus a 2-point curve-shape check (companion zero-span/max-span readings) for the three devices that have one.

    *Correction:* an earlier draft of this section stated that system-level validation used "characteristic curves ... reported by Tušek et al. (2010, 2011) and Nielsen et al. (2011)." That was aspirational rather than accurate — see ROADMAP.md Phase 7 for the full account. In brief: both papers are physically present in this repository, but neither contains the multi-point experimental characteristic-curve data this item needs (Tušek 2010 is a device-construction paper; Nielsen 2011 is a numerical-modeling review). The genuine curve source turned out to be a third paper, Tušek et al. (2013) (Figs. 10–11, 9 overlapping series) — present in `Papers/`, but not yet digitized. The benchmark set above is what the validation code actually uses today.

This two-stage validation approach improves confidence in both the underlying magnetocaloric material model and the overall AMR system simulation before comparison with conventional data-center cooling technologies.