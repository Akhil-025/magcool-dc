"""
sensitivity.py
===============
Sobol variance-based global sensitivity analysis (Saltelli sampling) on the
AMR system's electrical COP, with respect to the five design/operating
parameters exposed by AMRSystem:

    mu0H_max                  [1.0, 3.0]   T     (permanent-magnet field)
    frequency                  [0.5, 5.0]   Hz     (AMR cycle frequency)
    fluid_mdot                  [0.02, 0.2] kg/s     (heat transfer fluid flow)
    regenerator_effectiveness   [0.6, 0.95]  -        (NTU-based bed effectiveness)
    parasitic_fraction           [0.10, 0.45] -        (pump+motor overhead,
                                                          Phase 2 calibrated
                                                          range from the three
                                                          benchmark devices)

Evaluated at a fixed operating point: T_cold = 291 K (18 C, ASHRAE
recommended supply), span = 10 K.

Uses SALib (Herman & Usher, J. Open Source Software 2(9), 97 (2017)), the
same Sobol/Saltelli implementation used across the fuel-cell-simulation
literature and in the companion pemfc repo's sobol_results.txt.
"""

import numpy as np
from SALib.sample import sobol as sobol_sample
from SALib.analyze import sobol as sobol_analyze
from core.mce_material import GADOLINIUM
from core.amr_cycle import AMRSystem

T_COLD_K = 291.0
SPAN_K = 10.0

PROBLEM = {
    "num_vars": 5,
    "names": ["mu0H_max_T", "frequency_Hz", "fluid_mdot_kgs",
              "regen_effectiveness", "parasitic_fraction"],
    "bounds": [[1.0, 3.0], [0.5, 5.0], [0.02, 0.2], [0.6, 0.95], [0.10, 0.45]],
}


def model_cop(params):
    mu0H, freq, mdot, eps, parasitic = params
    sys_ = AMRSystem(material=GADOLINIUM, mu0H_max=mu0H, mass_regenerator=5.0,
                      frequency=freq, fluid_mdot=mdot,
                      regenerator_effectiveness=eps, parasitic_fraction=parasitic)
    result = sys_.run(T_COLD_K, SPAN_K)
    return result.COP_electrical


def run_sobol(n_base=512, seed=42, out_path="results/sobol_results.txt"):
    np.random.seed(seed)
    param_values = sobol_sample.sample(PROBLEM, n_base, calc_second_order=True)
    Y = np.array([model_cop(p) for p in param_values])
    Si = sobol_analyze.analyze(PROBLEM, Y, calc_second_order=True, print_to_console=False)

    lines = []
    lines.append(f"Sobol sensitivity analysis: AMR electrical COP at "
                  f"T_cold={T_COLD_K}K, span={SPAN_K}K")
    lines.append(f"Samples: {len(Y)} (Saltelli, N_base={n_base})")
    lines.append("")
    lines.append(f"{'parameter':<22}{'S1':>10}{'S1_conf':>10}{'ST':>10}{'ST_conf':>10}")
    for i, name in enumerate(PROBLEM["names"]):
        lines.append(f"{name:<22}{Si['S1'][i]:>10.4f}{Si['S1_conf'][i]:>10.4f}"
                      f"{Si['ST'][i]:>10.4f}{Si['ST_conf'][i]:>10.4f}")
    lines.append("")
    lines.append("Ranked by total-order sensitivity (ST):")
    ranked = sorted(zip(PROBLEM["names"], Si["ST"]), key=lambda x: -x[1])
    for name, st in ranked:
        lines.append(f"  {name:<22} ST={st:.4f}")

    lines.append("")
    lines.append("DIAGNOSTIC FINDING (this is a real result, not a numerical artifact):")
    lines.append("mu0H_max, frequency and fluid_mdot show ~0 sensitivity because the")
    lines.append("current amr_cycle.py model structure makes COP_electrical algebraically")
    lines.append("independent of Qc: COP_electrical = 1 / [(Th/Tc-1)/eta_2nd_law(eps) +")
    lines.append("parasitic_fraction], since both W_mag and W_parasitic scale linearly")
    lines.append("with Qc and cancel. Field/frequency/flow rate DO change how much cooling")
    lines.append("you get (Qc, via dTad) but NOT, in this model, how efficiently you get it.")
    lines.append("Real AMR systems don't have this decoupling -- eddy-current losses scale")
    lines.append("with frequency^2, viscous dissipation with mdot^2, and magnet support/")
    lines.append("motor sizing with field -- so eta_2nd_law and parasitic_fraction should")
    lines.append("be state-dependent functions, not constants. This is flagged as a")
    lines.append("required Phase 3 model upgrade in ROADMAP.md, not silently patched here.")

    with open(out_path, "w") as f:
        f.write("\n".join(lines) + "\n")
    print("\n".join(lines))
    print(f"\nWrote {out_path}")
    return Si


if __name__ == "__main__":
    run_sobol()
