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

## Phase 2 — system-level validation & sensitivity — done this session
- [x] Digitized published AMR prototype data (`data/amr_experimental_benchmarks.csv`:
      Astronautics/Jacobs 2014, DTU rotary Gd, Tušek single-bed Gd) as
      system-level COP/span validation targets
- [x] Calibrate-then-validate methodology (`validation_system.py`) —
      revealed that published COP figures are *electrical* (pump + motor
      overhead included), not thermodynamic-cycle-only; added a calibrated
      `parasitic_fraction` to `amr_cycle.py` and got the two comparable
      lab devices to ~10% agreement
- [x] Sobol/variance-based sensitivity analysis (`sensitivity.py` →
      `results/sobol_results.txt`) — found electrical COP is ~99.9%
      sensitive to `parasitic_fraction` alone in the current model
      structure, a genuine finding that motivates Phase 3
- [x] Response-surface (RSM) surrogate for cooling capacity Qc
      (`rsm.py` → `results/rsm_coefficients.txt`, R²=0.94 held-out)

## Phase 3 — state-dependent losses, optimization & multi-stage design
- [ ] **New, motivated by Phase 2's Sobol finding:** make `eta_2nd_law` and
      `parasitic_fraction` state-dependent functions instead of constants —
      eddy-current losses ~ frequency², viscous dissipation ~ mdot²,
      magnet/motor sizing ~ field — using correlations from Tušek et al.
      and Eriksen et al. (2015)
- [ ] Multi-objective optimization (NSGA-III or similar) trading off
      electrical COP, cooling capacity, magnet mass/cost, using the RSM
      surrogate as the fast inner-loop evaluator — mirrors `optimize.py`
- [ ] Cascade/multi-stage AMR design to cover the full 5-20 K ASHRAE span
      (Phase 1/2 model shows single-stage Qc collapse above ~16 K at 2 T,
      AND that single-stage electrical COP trails vapor-compression across
      most of the range — this is the key design finding motivating Phase 3)

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
