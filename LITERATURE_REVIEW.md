# Literature review — magnetocaloric cooling for data centers

## 1. Magnetocaloric effect & materials
- **Discovery of giant MCE**: Pecharsky & Gschneidner, *Phys. Rev. Lett.* 78, 4494 (1997) — Gd5(Si2Ge2), first-order transition giant magnetocaloric effect.
- **Gd baseline data**: Pecharsky & Gschneidner, *J. Magn. Magn. Mater.* 200 (1999) 44–56 — review/compilation of Gd magnetocaloric properties.
- **Direct ΔT_ad measurement (validation benchmark used here)**: Dan'kov, Tishin, Pecharsky & Gschneidner, *Phys. Rev. B* 57 (1998) 3478 — Gd adiabatic temperature change at 1, 2, 5 T near T_c = 294 K.
- **Mean-field theory limitations near T_c**: de Oliveira & von Ranke, *Phys. Rep.* 489 (2010) 89–159 — critical-fluctuation corrections to mean-field MCE models.
- **Perovskite manganite MCE**: Guo et al., *Appl. Phys. Lett.* 78, 1142 (1997); Phan & Yu, *J. Magn. Magn. Mater.* 308 (2007) 325–340 (review).

## 2. AMR cycle & system performance
- **AMR cycle concept**: Barclay, *US Patent 4332135* (1982) — original active magnetic regenerator concept.
- **Characteristic curve / performance modeling**: Tušek et al., *Int. J. Refrig.* 33 (2010); Nielsen et al., *Int. J. Refrig.* 34 (2011) 603–616.
- **Experimental AMR geometry comparison (parallel-plate vs. packed-bed)**: Tušek, Kitanovski, Prebil, Poredoš, *Appl. Therm. Eng.* (2011) — six-AMR comparison, 20 K span at 1.15 T, parallel-plate lowest porosity best COP.
- **AMR vs. vapor-compression COP at matched thermal levels**: same source — magnetocaloric COP exceeds classical vapor-compression COP between the same thermal reservoirs (cited claim, used to motivate the comparison in this study; the present model reproduces this qualitatively but at data-center-relevant spans the crossover is field- and span-dependent, see `results/comparison_table.csv`).
- **Second-law / exergetic performance**: Bahl & Nielsen, review chapters in Kitanovski et al., *Magnetocaloric Energy Conversion*, Springer (2015), Ch. 4 & 6.
- **Experimental AMR thermodynamic performance, packed Gd spheres**: (ScienceDirect 2016) — second-law efficiency peak at 15–20 K span, motivating the span range swept in `main.py`.
- **Review of AMR refrigerator developments**: (ScienceDirect 2020) — chronological review of materials, field sources, system design, performance indicators.
- **Optimized permanent-magnet (Halbach) design**: arXiv:1410.1987 — field/magnet design tradeoffs used to set the mu0*H_max = 2 T default in `main.py`.

## 3. Data center cooling baselines
- **Thermal envelope / operating points**: ASHRAE TC9.9, *Thermal Guidelines for Data Processing Environments*, 5th ed. (2021) — Class A1/A2 recommended supply range used to set T_cold in `main.py`.
- **Liquid cooling guidelines**: ASHRAE TC9.9, *Liquid Cooling Guidelines for Datacom Equipment Centers*, 2nd ed. (2021) — W-class facility water temperature bands, economizer-hour assumptions used in `baseline_cooling.py`.
- **Data-center cooling technology review**: Ebrahimi, Jones & Fleischer, *Renew. Sustain. Energy Rev.* 31 (2014) 622–638 — second-law efficiency benchmarks for DX/chiller plants.
- **Chip/system cooling power considerations**: Shah, Bash & Patel, *Cooling and Power Considerations for Chips*, ASME (2004).

## 4. Techno-economics
- **AMR system cost breakdown**: Bahl, Engelbrecht et al., *Int. J. Refrig.* 37 (2014) 78–83.
- **Magnetocaloric material cost review**: Franco, Blázquez et al., *Int. J. Refrig.* 57 (2018) 288–298.
- **Data center cooling cost benchmarks**: Lawrence Berkeley National Laboratory, Data Center Energy Efficiency program publications.

## Open items for Phase 2 (see ROADMAP.md)
- Digitize/tabulate the full Tušek et al. and Nielsen et al. characteristic curves as validation targets for `amr_cycle.py` (currently validated only at the material-property level via `validation.py`, not yet at the system/COP level).
- Source a rare-earth-magnet cost quote or NREL/DOE critical-materials report for a defensible CAPEX sensitivity range in `economics.py`.
- Confirm exact ASHRAE Region XV CRC 2026 (SAMUDRA, Chennai) paper format/deadline directly with the chapter — could not be confirmed via public search as of July 2026.
