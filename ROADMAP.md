# Roadmap to pemfc-suite parity

Mirroring the `pemfc` repo's structure realistically takes several build
sessions. This is the phased plan.

## Phase 1 — done
- [x] Mean-field MCE material model (`mce_material.py`), validated against
      Dan'kov et al. (1998) Gd data (`validation.py`)
- [x] 0-D AMR cycle model (`amr_cycle.py`)
- [x] Vapor-compression & liquid-cooling baseline COP models (`baseline_cooling.py`)
- [x] Order-of-magnitude economics comparison (`economics.py`)
- [x] Comparison driver across the ASHRAE thermal envelope (`main.py`)
- [x] Literature review with citations (`LITERATURE_REVIEW.md`)

## Phase 2 — system-level validation & sensitivity — done
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

## Phase 3 — state-dependent losses, optimization & multi-stage design — done
- [x] State-dependent loss model (`loss_model.py`): eddy ~ f²H², pumping ~ ṁ²,
      base overhead ~ Qc, calibrated against the 3 Phase 2 benchmark devices
      (exactly-determined 3-point fit — flagged as needing more data)
- [x] Re-ran Sobol with the new model: frequency (ST=0.68), flow (ST=0.22)
      and field (ST=0.09) now carry real sensitivity, resolving the Phase 2
      diagnostic finding
- [x] NSGA-III multi-objective optimization (`optimize.py`, COP vs. Qc vs.
      cost) → `results/pareto_front.csv` — exposed that `mass_regenerator`
      has no effect on Qc in the current model, so every Pareto design sits
      at the mass floor (a real finding, not hidden)
- [x] Multi-stage cascade design (`cascade.py`) → `results/cascade_comparison.csv`
      — confirms staging recovers the span range lost above 16 K, but ALSO
      shows both single- and multi-stage Gd/2T AMR trail vapor-compression
      and liquid cooling on electrical COP across the whole ASHRAE range at
      this design point. This is the project's current honest conclusion.

## Phase 4 — close the model gaps Phase 3 exposed — done
- [x] **NTU-based `thermal.py`**: packed-sphere-bed regenerator effectiveness
      from Wakao & Kaguei (1982) Nusselt correlation + utilization-factor
      degradation, wired into `amr_cycle.py` via `use_ntu_thermal_model`.
      Confirmed the fix: NSGA-III Pareto front mass values now spread
      1.0-14.5 kg instead of pinning at the floor
- [x] **Attempted the higher-field/giant-MCE re-run** — and found the real
      blocker: `mce_material.py`'s mean-field/Brillouin framework can't
      capture Gd5Si2Ge2's first-order-transition giant MCE (underpredicts
      ΔT_ad by ~10x), so `cascade.py` correctly returns zero capacity for it
      at the DC operating point rather than a misleadingly small nonzero
      number. Flagged explicitly in `mce_material.py`. Gd's own T_c=294K
      (21°C) already sits inside the ASHRAE 18-27°C range, which is worth
      noting in the paper as a reason Gd, not an exotic alloy, is the
      practical near-term candidate.

## Phase 5 — first-order-transition materials, economics/emissions — done
- [x] **Bean-Rodbell-family model**: `first_order_mce.py` implements an
      extended (6th-order) Landau free-energy model, calibrated to Gd5Si2Ge2's
      literature peak entropy change (~-18 J/(kg K) at 5T) — resolves the
      Phase 4 blocker with a physically appropriate framework, not a patch
- [x] **Giant-MCE vs. Gd formal comparison** (`giant_mce_analysis.py` →
      `results/giant_mce_analysis.txt`): confirmed the Curie-matching
      principle cleanly — Gd5Si2Ge2's peak-effect window (286.4K) sits ~5K
      below the ASHRAE range and collapses to zero there, but performs
      strongly (Qc=5319W, COP_elec=7.76) at its own matched point. Does NOT
      overturn the Phase 1-4 conclusion; identifies the concrete next step
      (composition-tuned Gd5(SixGe1-x)4, literature-documented as tunable)
- [x] **Refrigerant-free GWP/emissions comparison** (`emissions.py`):
      quantified honestly — real categorical benefit (zero leak/phase-out
      risk) but does not overturn the emissions comparison until AMR's COP
      gap closes; report both numbers plainly in the paper

## Phase 6 — extended validation stress-test + grounded economics — done
- [x] **Expanded `data/amr_experimental_benchmarks.csv`**: found a fully
      usable 4th device (Okamura & Hirano 2013, Qc=200W/COP=2.5/5K span/
      1.1T/1kg Gd) via targeted search. Added to `loss_model.py` as an
      EXTENDED calibration set.
- [x] **Stress-tested the loss model — and it failed the stress test.**
      The 4-point least-squares fit gives negative (unphysical) k_pump and
      base_frac, and leave-one-out CV shows errors up to +1639% — four
      orders of magnitude of device scale (6.5W to 2502W) apparently can't
      be pooled into one linear loss model. Reported via
      `run_extended_diagnostic()`, NOT silently absorbed into the
      production default, which stays on the stable 3-point CORE fit.
- [x] Attempted to add the Risoe/DTU 2011 (30K span) device too — could not
      calibrate; consistent with, and further evidence for, the Phase 1
      single-stage span-collapse finding.
- [x] **Grounded `economics.py`** in an actual costing study: Bjørk, Bahl &
      Smith (Int. J. Refrig. 34 (2011) 1805-1816), $40/kg magnet + $20/kg
      MCM, replacing the earlier loosely-sourced placeholder. Wired into
      `optimize.py`'s cost objective via `economics.material_cost()`.

## Phase 7 — remaining open items before the paper is fully evidence-complete
- [ ] **Loss-model scale term**: add a size/capacity-dependent term (or
      separate small-device vs. large-device regimes) to `loss_model.py` —
      Phase 6 showed the current linear form can't span 6.5W-2502W devices.
      Needed before the EXTENDED calibration set can be trusted.
- [ ] Curie-graded cascade once/if a composition-tuned Gd5(SixGe1-x)4 (or
      La(Fe,Si)13Hy — used in the Astronautics benchmark device, also not
      yet modeled here) is added with independent literature validation
      (unlike Gd5Si2Ge2's single-point calibration, see `first_order_mce.py`
      honesty flag #2 — cross-check against Giguere et al. 1999's direct
      ΔT_ad measurement before trusting this material further)
- [ ] Digitize published AMR characteristic curves (Tušek 2010, Nielsen 2011)
      for a true system-level COP-vs-span validation (current validation is
      point-wise, see `validation_system.py`)
- [ ] Full-system cost (heat exchangers, pumps, motor/drive, controls) to
      complement `economics.material_cost()`'s materials-only floor
- [ ] Assemble the ASHRAE Region XV CRC paper (LaTeX/Word) + slide deck —
      blocked on the chapter confirming the exact 2026 SAMUDRA submission
      format/deadline (see LITERATURE_REVIEW.md open items)
- [ ] Optional COMSOL 2-D/3-D regenerator-bed model + setup guide for anyone
      wanting to extend beyond the 0-D model

## Immediate next step
Confirm the ASHRAE Region XV CRC 2026 (SAMUDRA, Chennai) submission
guidelines (page limit, format, deadline) directly with the chapter —
this could not be found via public search as of July 2026, and it
determines how the paper write-up should be structured.
