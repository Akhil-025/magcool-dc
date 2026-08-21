"""
pareto_multiseed_stability.py
==============================
Phase 30 addition (statistical rigor pass, Tier 1 item 5 of the original
review).

`core/optimize.py`'s `run_optimization()` (the main material+geometry
NSGA-III co-optimization behind `results/pareto_front.csv`) has, until
now, only ever been run at a single fixed `seed=1` for its PRODUCTION
output -- unlike `core/hysteresis_sensitivity.py`'s
`run_hysteresis_multiseed_stability_check()` (Phase 16-18 follow-up) and
`core/magnet_geometry.py`'s multiseed check (Paper-Mining Pass review item
4), which already established the pattern this module extends to the
MAIN Pareto front itself.

This module reruns `run_optimization()` across multiple independent NSGA-
III seeds (same pop_size/n_gen as production, only `seed` varies) and
reports variance/spread on the specific headline numbers a paper would
cite:
  - best (max) COP_electrical on the merged front
  - the knee-point (balanced) design's COP_electrical, Qc, and cost
  - each material family's SHARE (%) of the merged non-dominated front
    (the Phase 15/16 "100% La(Fe,Si)13Hy"-style claims this repo has
    made before, and the exact kind of claim `hysteresis_sensitivity.py`
    already found could reverse under different search settings)

This does NOT replace `run_hysteresis_multiseed_stability_check()` (which
tests hysteresis on/off, a different axis) -- it tests plain seed-to-seed
NSGA-III search variance at fixed production settings, the more basic and
previously entirely unchecked question of "if I reran this search with a
different random seed, would the headline number substantially move."
"""

import numpy as np

from core import optimize as optimize_module


def run_pareto_multiseed_stability_check(seeds=(1, 2, 3, 4, 5),
                                          pop_size=40, n_gen=25,
                                          verbose=True):
    """Reruns `optimize.run_optimization()` at production pop_size/n_gen
    across `seeds`, collecting the same headline numbers a paper would
    report, and summarizes their seed-to-seed spread."""
    per_seed = []
    for seed in seeds:
        if verbose:
            print(f"\n--- NSGA-III seed={seed} ---")
        rows = optimize_module.run_optimization(
            pop_size=pop_size, n_gen=n_gen, seed=seed,
            out_csv=f"results/pareto_front_seed{seed}.csv",
            per_material_out_dir=None)  # skip per-material CSVs for the sweep
        if not rows:
            per_seed.append(None)
            continue
        F_cop = [r["COP_electrical"] for r in rows]
        best_cop = max(F_cop)
        # knee point, same normalized-distance definition run_optimization() uses
        F = np.array([[-r["COP_electrical"], -r["Qc_W"], r["cost_index_USD"]] for r in rows])
        Fn = (F - F.min(axis=0)) / (F.max(axis=0) - F.min(axis=0) + 1e-12)
        knee_idx = int(np.argmin(np.linalg.norm(Fn, axis=1)))
        knee = rows[knee_idx]
        material_counts = {}
        for r in rows:
            material_counts[r["material"]] = material_counts.get(r["material"], 0) + 1
        material_share = {m: 100 * c / len(rows) for m, c in material_counts.items()}
        per_seed.append({
            "seed": seed,
            "n_front": len(rows),
            "best_COP_electrical": best_cop,
            "knee_COP_electrical": knee["COP_electrical"],
            "knee_Qc_W": knee["Qc_W"],
            "knee_cost_USD": knee["cost_index_USD"],
            "knee_material": knee["material"],
            "material_share_pct": material_share,
        })

    valid = [s for s in per_seed if s is not None]
    if not valid:
        if verbose:
            print("All seeds returned an empty front -- cannot summarize.")
        return {"per_seed": per_seed, "summary": None}

    best_cops = [s["best_COP_electrical"] for s in valid]
    knee_cops = [s["knee_COP_electrical"] for s in valid]
    knee_qcs = [s["knee_Qc_W"] for s in valid]
    knee_costs = [s["knee_cost_USD"] for s in valid]
    all_materials = sorted({m for s in valid for m in s["material_share_pct"]})
    material_share_stats = {}
    for m in all_materials:
        shares = [s["material_share_pct"].get(m, 0.0) for s in valid]
        material_share_stats[m] = {
            "mean_pct": round(float(np.mean(shares)), 1),
            "std_pct": round(float(np.std(shares)), 1),
            "min_pct": round(float(np.min(shares)), 1),
            "max_pct": round(float(np.max(shares)), 1),
        }
    knee_materials = {s["knee_material"] for s in valid}

    summary = {
        "n_seeds": len(valid),
        "best_COP_electrical_mean": round(float(np.mean(best_cops)), 3),
        "best_COP_electrical_std": round(float(np.std(best_cops)), 3),
        "best_COP_electrical_range": (round(min(best_cops), 3), round(max(best_cops), 3)),
        "knee_COP_electrical_mean": round(float(np.mean(knee_cops)), 3),
        "knee_COP_electrical_std": round(float(np.std(knee_cops)), 3),
        "knee_Qc_W_mean": round(float(np.mean(knee_qcs)), 1),
        "knee_Qc_W_std": round(float(np.std(knee_qcs)), 1),
        "knee_cost_USD_mean": round(float(np.mean(knee_costs)), 1),
        "knee_cost_USD_std": round(float(np.std(knee_costs)), 1),
        "knee_material_consistent_across_seeds": len(knee_materials) == 1,
        "knee_materials_seen": sorted(knee_materials),
        "material_share_stats": material_share_stats,
    }

    if verbose:
        print(f"\n=== Multiseed stability summary ({summary['n_seeds']} seeds) ===")
        print(f"best COP_electrical: mean={summary['best_COP_electrical_mean']:.3f}  "
              f"std={summary['best_COP_electrical_std']:.3f}  "
              f"range={summary['best_COP_electrical_range']}")
        print(f"knee-point COP_electrical: mean={summary['knee_COP_electrical_mean']:.3f}  "
              f"std={summary['knee_COP_electrical_std']:.3f}")
        print(f"knee-point Qc: mean={summary['knee_Qc_W_mean']:.1f}W  "
              f"std={summary['knee_Qc_W_std']:.1f}W")
        print(f"knee-point cost: mean=${summary['knee_cost_USD_mean']:.1f}  "
              f"std=${summary['knee_cost_USD_std']:.1f}")
        print(f"knee-point material consistent across seeds: "
              f"{summary['knee_material_consistent_across_seeds']} "
              f"(materials seen: {summary['knee_materials_seen']})")
        print("material share of the merged Pareto front, mean +/- std across seeds:")
        for m, stats in sorted(material_share_stats.items(), key=lambda kv: -kv[1]["mean_pct"]):
            print(f"  {m:<40} {stats['mean_pct']:5.1f}% +/- {stats['std_pct']:4.1f}%  "
                  f"(range {stats['min_pct']:.0f}-{stats['max_pct']:.0f}%)")
        print()
        print("HONEST FRAMING FOR THE PAPER: report every headline Pareto "
              "number (best COP, knee-point design, material-family share) "
              "as mean +/- std across these seeds, not a single-seed point "
              "value -- exactly the discipline this repo's own "
              "hysteresis_multiseed_stability.txt already applied to the "
              "hysteresis on/off question, extended here to plain seed-to-"
              "seed NSGA-III variance at fixed production settings. If "
              "knee_material_consistent_across_seeds is False, DO NOT "
              "report a single 'best material' headline claim without "
              "this caveat -- state the seed-to-seed material split "
              "explicitly instead.")

    return {"per_seed": per_seed, "summary": summary}


def write_pareto_multiseed_stability_report(
        path="results/pareto_multiseed_stability.txt",
        seeds=(1, 2, 3, 4, 5), pop_size=40, n_gen=25):
    import io, contextlib
    buf = io.StringIO()
    with contextlib.redirect_stdout(buf):
        run_pareto_multiseed_stability_check(seeds=seeds, pop_size=pop_size, n_gen=n_gen)
    with open(path, "w") as f:
        f.write(buf.getvalue())
    print(buf.getvalue())
    print(f"Wrote {path}")


if __name__ == "__main__":
    write_pareto_multiseed_stability_report()
