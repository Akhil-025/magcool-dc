# magcool-dc

Physics-based simulation suite evaluating **magnetocaloric (magnetic) cooling
for data centers**, benchmarked against vapor-compression CRAC/CRAH and
direct liquid cooling — built as the magnetic-cooling counterpart to the
[`pemfc`](https://github.com/Akhil-025/pemfc) PEM fuel cell simulation suite,
for an ASHRAE Region XV student paper submission.

## Status: Phase 1 (core physics + baseline comparison) — see `ROADMAP.md`

## What's implemented

| Module | Purpose |
|---|---|
| `core/mce_material.py` | Mean-field (Brillouin/Weiss) model of the magnetocaloric effect; materials library (Gd, Gd5Si2Ge2, La0.7Ca0.3MnO3) parameterized from published Tc, J, g, Debye temperature |
| `core/amr_cycle.py` | 0-D Active Magnetic Regenerator (AMR) cycle model: cooling capacity, COP, second-law efficiency vs. temperature span, field, frequency |
| `core/baseline_cooling.py` | Vapor-compression and liquid-cooling COP correlations for data centers, referenced to published second-law efficiencies |
| `core/economics.py` | Order-of-magnitude CAPEX/OPEX comparison |
| `core/validation.py` | Validates the MCE model against Dan'kov et al. (1998) experimental Gd ΔT_ad data |
| `main.py` | Runs the full comparison across the ASHRAE TC9.9 thermal envelope, writes `results/comparison_table.csv` |

## Quick start

```bash
pip install -r requirements.txt
python -m core.validation     # sanity-check the MCE model against literature
python main.py                  # run the AMR vs. VCC vs. liquid-cooling comparison
```

## Validation snapshot (Gd, mean-field model vs. Dan'kov et al. 1998)

| Field | Literature ΔT_ad | Model ΔT_ad | Error |
|---|---|---|---|
| 1 T | 3.2 K | 4.43 K | +38.4% |
| 2 T | 6.3 K | 7.16 K | +13.6% |
| 5 T | 14.6 K | 13.18 K | -9.7% |

Mean-field theory is known to overpredict ΔT_ad near T_c because it
neglects short-range spin correlations / critical fluctuations (de Oliveira
& von Ranke, *Phys. Rep.* 489 (2010) 89-159). This is documented, not
hidden — the roadmap's Phase 2 replaces the low-field region with a
Landau-expansion or tabulated-experimental-data correction.

## Repo layout

```
core/            physics + economics modules
data/            literature-sourced parameter tables
results/         generated comparison tables / figures
main.py          top-level comparison driver
ROADMAP.md       phased plan to reach pemfc-suite parity
LITERATURE_REVIEW.md
NOMENCLATURE.md
```
