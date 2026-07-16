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

#### Magnetocaloric Theory

* **de Oliveira & von Ranke**, *Physics Reports*, **489**, 89–159 (2010).

  * Comprehensive review of thermodynamic models describing the magnetocaloric effect.
  * Highlights the limitations of mean-field theory near the Curie temperature due to critical magnetic fluctuations and discusses improved theoretical formulations.

#### Alternative Magnetocaloric Materials

* **Guo et al.**, *Applied Physics Letters*, **78**, 1142 (1997).
* **Phan & Yu**, *Journal of Magnetism and Magnetic Materials*, **308**, 325–340 (2007).

  * Investigated perovskite manganites and other rare-earth alloys exhibiting significant room-temperature magnetocaloric effects.
  * Discuss material selection criteria including entropy change, hysteresis, thermal conductivity and operating temperature.

---

## 2. Active Magnetic Regenerator (AMR) Systems

#### Active Magnetic Regenerator Concept

* **Barclay**, US Patent **4332135** (1982).

  * Introduced the Active Magnetic Regenerator (AMR) cycle.
  * Demonstrated how cyclic magnetization and demagnetization combined with regenerative heat exchange could achieve practical magnetic refrigeration.

#### Performance Modelling

* **Tušek et al.**, *International Journal of Refrigeration*, **33** (2010).
* **Nielsen et al.**, *International Journal of Refrigeration*, **34**, 603–616 (2011).

  * Developed characteristic-curve models relating cooling capacity, utilization factor, operating frequency, pressure drop and regenerator effectiveness.
  * These models provide widely accepted benchmarks for predicting AMR performance and validating numerical simulations.

#### Regenerator Geometry

* **Tušek, Kitanovski, Prebil & Poredoš**, *Applied Thermal Engineering* (2011).

  * Experimentally compared six AMR configurations including packed-bed and parallel-plate regenerators.
  * Achieved approximately **20 K temperature span** under a **1.15 T** magnetic field.
  * Reported that optimized parallel-plate regenerators produced the highest COP because of reduced flow resistance and improved heat transfer.
  * Demonstrated that magnetic refrigeration can outperform vapor-compression systems under certain operating conditions, although performance depends strongly on magnetic field strength and temperature span.

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

---

### Validation Strategy for This Work (Optional Section)

The numerical model developed in this study is validated at two complementary levels:

1. **Material-level validation** using the experimental magnetocaloric data of **Dan'kov et al. (1998)** for gadolinium, including adiabatic temperature change and entropy variation.

2. **System-level validation** using the characteristic curves and experimental AMR performance reported by **Tušek et al. (2010, 2011)** and **Nielsen et al. (2011)**, including cooling capacity, coefficient of performance (COP), utilization factor and temperature span.

This two-stage validation approach improves confidence in both the underlying magnetocaloric material model and the overall AMR system simulation before comparison with conventional data-center cooling technologies.
