# Roadmap to pemfc-suite parity

Mirroring the `pemfc` repo's structure realistically takes several build
sessions. This is the phased plan; **Phase 1 is done** as of this session.

## Phase 1 — done this session
- [x] Mean-field MCE material model (`mce_material.py`), validated against
      Dan'kov et al. (1998) Gd data (`validation.py`)
- [x] 0-D AMR cycle model (`amr_cycle.py`)
- [x] Vapor-compression & liquid-cooling baseline COP models (`baseline_cooling.py`)
- [x] Order-of-magnitude economics comparison (`economics.py`)
- [x] Comparison driver across the ASHRAE thermal envelope (`main.py`)
- [x] Literature review with citations (`LITERATURE_REVIEW.md`)

## Phase 2 — system-level validation & sensitivity
- [ ] Digitize published AMR characteristic curves (Tušek 2010, Nielsen 2011)
      as COP/span validation targets — currently only material-level ΔT_ad
      is validated, not full-system COP
- [ ] Sobol/variance-based sensitivity analysis on (field, frequency,
      utilization, regenerator effectiveness) → COP, mirroring
      `sobol_results.txt` in pemfc
- [ ] Response-surface (RSM) surrogate for fast optimization, mirroring
      `rsm.py` / `rsm_coefficients.txt`

## Phase 3 — optimization & multi-stage design
- [ ] Multi-objective optimization (NSGA-III or similar) trading off COP,
      cooling capacity, magnet mass/cost — mirrors `optimize.py`
- [ ] Cascade/multi-stage AMR design to cover the full 5–20 K ASHRAE span
      (Phase 1 model shows single-stage collapse above ~16 K at 2 T — this
      is the key design finding motivating Phase 3)

## Phase 4 — thermal/geometric detail
- [ ] NTU-based regenerator effectiveness model replacing the fixed 0.85
      placeholder (`thermal.py`)
- [ ] Optional COMSOL 2-D/3-D regenerator-bed model + setup guide, mirroring
      `COMSOL_setup.md`, for anyone wanting to extend beyond the 0-D model

## Phase 5 — economics, emissions, and report/paper assembly
- [ ] Refine CAPEX sensitivity with a real rare-earth-magnet cost source
- [ ] Refrigerant-free GWP/emissions comparison (magnetic cooling uses no
      HFC/HFO refrigerant — quantify the avoided-emissions case)
- [ ] Assemble the ASHRAE Region XV CRC paper (LaTeX/Word) + slide deck from
      the validated results, once the chapter confirms the exact 2026
      SAMUDRA submission format and deadline

## Immediate next step
Confirm the ASHRAE Region XV CRC 2026 (SAMUDRA, Chennai) submission
guidelines (page limit, format, deadline) directly with the chapter —
this could not be found via public search as of July 2026, and it
determines how the Phase 5 write-up should be structured.
