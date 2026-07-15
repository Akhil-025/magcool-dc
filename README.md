# magcool-dc

Physics-based simulation suite evaluating **magnetocaloric (magnetic) cooling
for data centers**, benchmarked against vapor-compression CRAC/CRAH and
direct liquid cooling — built as the magnetic-cooling counterpart to the
[`pemfc`](https://github.com/Akhil-025/pemfc) PEM fuel cell simulation suite,
for an ASHRAE Region XV student paper submission.

## Status: Phase 2 (system-level validation + sensitivity + surrogate) — see `ROADMAP.md`

## What's implemented

| Module | Purpose |
|---|---|
| `core/mce_material.py` | Mean-field (Brillouin/Weiss) model of the magnetocaloric effect; materials library (Gd, Gd5Si2Ge2, La0.7Ca0.3MnO3) parameterized from published Tc, J, g, Debye temperature |
| `core/amr_cycle.py` | 0-D Active Magnetic Regenerator (AMR) cycle model: cooling capacity, ideal magnetic-cycle COP, and parasitic-inclusive electrical COP |
| `core/baseline_cooling.py` | Vapor-compression and liquid-cooling COP correlations for data centers, referenced to published second-law efficiencies |
| `core/economics.py` | Order-of-magnitude CAPEX/OPEX comparison |
| `core/validation.py` | Validates the MCE material model against Dan'kov et al. (1998) experimental Gd ΔT_ad data |
| `core/validation_system.py` | **Phase 2:** calibrates + validates the AMR *system* model against three published prototypes (Astronautics, DTU, Tušek) |
| `core/sensitivity.py` | **Phase 2:** Sobol global sensitivity analysis of electrical COP over 5 design parameters |
| `core/rsm.py` | **Phase 2:** quadratic response-surface surrogate for cooling capacity, R²=0.94 held-out, for use in Phase 3 optimization |
| `main.py` | Runs the full comparison across the ASHRAE TC9.9 thermal envelope, writes `results/comparison_table.csv` |

## Quick start

```bash
pip install -r requirements.txt
python -m core.validation          # MCE material model vs. literature
python -m core.validation_system    # AMR system model vs. published prototypes
python -m core.sensitivity           # Sobol sensitivity -> results/sobol_results.txt
python -m core.rsm                    # RSM surrogate -> results/rsm_coefficients.txt
python main.py                         # AMR vs. VCC vs. liquid-cooling comparison
```

## Key Phase 2 finding: ideal vs. electrical COP

Published AMR COP figures are **electrical** (include pump + magnet-motor-drive
power). The model's *ideal* magnetic-cycle-only COP overpredicts these by
118-619%. Adding a calibrated parasitic-loss fraction (default 0.15,
literature range 0.12-0.45 across three benchmark devices) brings the two
comparable lab-scale devices to single-digit-percent agreement:

| Device | Span | COP (lit / ideal / electrical) | Error (electrical) |
|---|---|---|---|
| DTU rotary Gd (2016) | 10.1 K | 4.20 / 14.88 / 4.60 | +9.6% |
| Tušek single-bed Gd (2010) | 15.0 K | 4.60 / 10.02 / 4.00 | -13.0% |
| Astronautics rotary (2014, naval-cooler scale) | 11.0 K | 1.90 / 13.66 / 4.48 | +135.8% (outlier — paper itself cites "mediocre" electrical-component efficiency at that scale) |

**Using the correct (electrical) COP, single-stage AMR at 2 T does NOT beat
vapor-compression across most of the ASHRAE 5-20 K span range** at this
design point (AMR electrical COP ~4-5.5 vs. VCC ~6-24) — a materially
different conclusion than the naive ideal-COP comparison would suggest, and
the central design-motivation for Phase 3 (multi-stage, higher field,
state-dependent loss modeling).

## Sobol sensitivity: a genuine model-structure finding

Sobol analysis (`results/sobol_results.txt`) found electrical COP is ~99.9%
sensitive to `parasitic_fraction` alone — field, frequency, and flow rate
show ~0 sensitivity to COP (though they do affect *cooling capacity*, Qc).
This is because the current model makes `eta_2nd_law` and
`parasitic_fraction` constants rather than state-dependent functions, so
they cancel out of the COP ratio algebraically. This is flagged as a
required Phase 3 upgrade, not silently patched.

## Validation snapshot (Gd MCE, mean-field model vs. Dan'kov et al. 1998)

| Field | Literature ΔT_ad | Model ΔT_ad | Error |
|---|---|---|---|
| 1 T | 3.2 K | 4.43 K | +38.4% |
| 2 T | 6.3 K | 7.16 K | +13.6% |
| 5 T | 14.6 K | 13.18 K | -9.7% |

Mean-field theory is known to overpredict ΔT_ad near T_c because it
neglects short-range spin correlations / critical fluctuations (de Oliveira
& von Ranke, *Phys. Rep.* 489 (2010) 89-159).

## Repo layout

```
core/            physics, validation, sensitivity, surrogate, and economics modules
data/            literature-sourced parameter tables + digitized prototype benchmarks
results/         generated comparison tables, Sobol results, RSM coefficients
main.py          top-level comparison driver
ROADMAP.md       phased plan to reach pemfc-suite parity
LITERATURE_REVIEW.md
NOMENCLATURE.md
```
