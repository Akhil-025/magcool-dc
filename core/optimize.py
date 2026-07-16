"""
optimize.py
============
Phase 3: multi-objective design optimization of the AMR system using
NSGA-III (Deb & Jain, IEEE Trans. Evol. Comput. 18(4), 577-601 (2014)),
via the pymoo implementation (Blank & Deb, IEEE Access 8, 89497-89509 (2020)).

Design variables:
    mu0H_max              [1.0, 3.0]   T
    frequency               [0.3, 5.0]   Hz
    fluid_mdot                [0.02, 0.5] kg/s
    mass_regenerator            [1.0, 15.0] kg
    regenerator_effectiveness    [0.6, 0.95]  -

Objectives (all converted to minimization for pymoo):
    f1 = -COP_electrical        (maximize electrical COP; uses the Phase 3
                                   state-dependent loss_model, NOT the Phase
                                   2 constant, so the field/frequency/mdot
                                   choices carry a genuine efficiency cost)
    f2 = -Qc                     (maximize cooling capacity)
    f3 = cost_index              (minimize a literature-informed CAPEX proxy:
                                   Gd material cost ~$175/kg (Franco et al.
                                   2018) + a magnet-cost term scaling with
                                   field^2 * regenerator mass as a rough
                                   permanent-magnet-volume proxy, consistent
                                   with the order-of-magnitude figures in
                                   economics.py / Bahl et al. 2014)

Fixed operating point: T_cold = 291 K, span = 10 K (same as sensitivity.py,
for direct comparability of results).

Output: Pareto front written to results/pareto_front.csv, plus a short
text summary of the extreme and knee-point designs.
"""

import numpy as np
import csv
from pymoo.core.problem import ElementwiseProblem
from pymoo.algorithms.moo.nsga3 import NSGA3
from pymoo.util.ref_dirs import get_reference_directions
from pymoo.optimize import minimize as pymoo_minimize
from core.mce_material import GADOLINIUM
from core.amr_cycle import AMRSystem
from core.loss_model import StateDependentLossModel
from core.economics import material_cost

T_COLD_K = 291.0
SPAN_K = 10.0

_LOSS_MODEL = StateDependentLossModel()
USE_NTU_THERMAL_MODEL = True   # Phase 4: mass_regenerator now has a real
                                 # effect on Qc via core/thermal.py; set False
                                 # to reproduce the Phase 3 mass-floor result


def cost_index(mu0H, mass_regenerator):
    return material_cost(mu0H, mass_regenerator)


class AMRDesignProblem(ElementwiseProblem):
    def __init__(self):
        super().__init__(
            n_var=5, n_obj=3, n_constr=0,
            xl=np.array([1.0, 0.3, 0.02, 1.0, 0.6]),
            xu=np.array([3.0, 5.0, 0.5, 15.0, 0.95]),
        )

    def _evaluate(self, x, out, *args, **kwargs):
        mu0H, freq, mdot, mass, eps = x
        sys_ = AMRSystem(material=GADOLINIUM, mu0H_max=mu0H, mass_regenerator=mass,
                          frequency=freq, fluid_mdot=mdot, regenerator_effectiveness=eps,
                          loss_model=_LOSS_MODEL, use_ntu_thermal_model=USE_NTU_THERMAL_MODEL)
        result = sys_.run(T_COLD_K, SPAN_K)
        f1 = -result.COP_electrical
        f2 = -result.Qc
        f3 = cost_index(mu0H, mass)
        out["F"] = [f1, f2, f3]


def run_optimization(pop_size=60, n_gen=40, seed=1,
                      out_csv="results/pareto_front.csv"):
    ref_dirs = get_reference_directions("das-dennis", 3, n_partitions=8)
    algorithm = NSGA3(pop_size=pop_size, ref_dirs=ref_dirs)
    problem = AMRDesignProblem()
    res = pymoo_minimize(problem, algorithm, ("n_gen", n_gen), seed=seed, verbose=False)

    X, F = res.X, res.F
    rows = []
    for x, f in zip(X, F):
        rows.append({
            "mu0H_max_T": round(x[0], 3), "frequency_Hz": round(x[1], 3),
            "fluid_mdot_kgs": round(x[2], 4), "mass_regenerator_kg": round(x[3], 2),
            "regen_effectiveness": round(x[4], 3),
            "COP_electrical": round(-f[0], 3), "Qc_W": round(-f[1], 2),
            "cost_index_USD": round(f[2], 1),
        })

    with open(out_csv, "w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=rows[0].keys())
        writer.writeheader()
        writer.writerows(rows)

    # Highlight three extreme designs + a naive "knee point" (min L2 distance
    # to the per-objective ideal point in normalized objective space)
    Fn = (F - F.min(axis=0)) / (F.max(axis=0) - F.min(axis=0) + 1e-12)
    knee_idx = np.argmin(np.linalg.norm(Fn, axis=1))
    best_cop_idx = np.argmax(-F[:, 0])
    best_qc_idx = np.argmax(-F[:, 1])
    best_cost_idx = np.argmin(F[:, 2])

    print(f"NSGA-III optimization: {len(rows)} Pareto-optimal designs found "
          f"(pop={pop_size}, gen={n_gen})")
    print(f"Wrote {out_csv}\n")
    for label, idx in [("Best electrical COP", best_cop_idx),
                        ("Best cooling capacity", best_qc_idx),
                        ("Lowest cost", best_cost_idx),
                        ("Knee point (balanced)", knee_idx)]:
        r = rows[idx]
        print(f"{label:<24} H={r['mu0H_max_T']}T  f={r['frequency_Hz']}Hz  "
              f"mdot={r['fluid_mdot_kgs']}kg/s  mass={r['mass_regenerator_kg']}kg  "
              f"eps={r['regen_effectiveness']}  -> COP_elec={r['COP_electrical']}, "
              f"Qc={r['Qc_W']}W, cost=${r['cost_index_USD']}")

    masses = [r["mass_regenerator_kg"] for r in rows]
    if max(masses) - min(masses) < 0.5 and not USE_NTU_THERMAL_MODEL:
        print("\nDIAGNOSTIC FINDING (Phase 3, constant-eps mode): every Pareto-optimal "
              "design sits at the regenerator-mass lower bound because mass_regenerator "
              "doesn't affect Qc without the NTU thermal model. Set USE_NTU_THERMAL_MODEL "
              "= True (default since Phase 4) to see the real mass/cost/performance tradeoff.")
    elif USE_NTU_THERMAL_MODEL:
        print(f"\nPHASE 4: mass_regenerator now spans {min(masses):.2f}-{max(masses):.2f} kg "
              "across the Pareto front (core/thermal.py's NTU model makes more regenerator "
              "material genuinely improve eps/Qc, trading off against the cost objective, "
              "rather than being pure waste as in the Phase 3 constant-eps result).")
    return rows


if __name__ == "__main__":
    run_optimization()
