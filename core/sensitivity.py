"""
sensitivity.py
===============
Sobol variance-based global sensitivity analysis (Saltelli sampling) on the
AMR system's electrical COP, with respect to the five design/operating
parameters exposed by AMRSystem:

    mu0H_max                     [1.0, 3.0]   T        (permanent-magnet field)
    frequency                    [0.5, 5.0]   Hz       (AMR cycle frequency)
    fluid_mdot                   [0.02, 0.2] kg/s      (heat transfer fluid flow)
    regenerator_effectiveness    [0.6, 0.95]  -        (NTU-based bed effectiveness)
    parasitic_fraction           [0.10, 0.45] -        (constant parasitic
                                                        power fraction used
                                                        for comparison with
                                                        the state-dependent
                                                        loss model)

Evaluated at a fixed operating point: T_cold = 291 K (18 C, ASHRAE
recommended supply), span = 10 K.

Uses SALib (Herman & Usher, J. Open Source Software 2(9), 97 (2017)) for
Sobol variance-based global sensitivity analysis.
"""

import numpy as np
from SALib.sample import sobol as sobol_sample
from SALib.analyze import sobol as sobol_analyze
from core.mce_material import GADOLINIUM
from core.amr_cycle import AMRSystem
from core.loss_model import StateDependentLossModel

T_COLD_K = 291.0
SPAN_K = 10.0

PROBLEM = {
    "num_vars": 5,
    "names": ["mu0H_max_T", "frequency_Hz", "fluid_mdot_kgs",
              "regen_effectiveness", "parasitic_fraction"],
    "bounds": [[1.0, 3.0], [0.5, 5.0], [0.02, 0.2], [0.6, 0.95], [0.10, 0.45]],
}

# The fifth parameter is used only with the constant-loss formulation.
# When the state-dependent loss model is enabled, it is ignored so that
# both analyses use the same Sobol sampling design and remain directly
# comparable.
_LOSS_MODEL = StateDependentLossModel()


def model_cop(params, use_state_dependent_losses=False):
    mu0H, freq, mdot, eps, parasitic = params
    kwargs = dict(material=GADOLINIUM, mu0H_max=mu0H, mass_regenerator=5.0,
                  frequency=freq, fluid_mdot=mdot, regenerator_effectiveness=eps)
    if use_state_dependent_losses:
        sys_ = AMRSystem(loss_model=_LOSS_MODEL, **kwargs)
    else:
        sys_ = AMRSystem(parasitic_fraction=parasitic, **kwargs)
    result = sys_.run(T_COLD_K, SPAN_K)
    return result.COP_electrical


def run_sobol(n_base=64, seed=42, out_path="results/sobol_results.txt",
              use_state_dependent_losses=False):
    np.random.seed(seed)
    param_values = sobol_sample.sample(PROBLEM, n_base, calc_second_order=True)
    Y = np.array([model_cop(p, use_state_dependent_losses) for p in param_values])
    Si = sobol_analyze.analyze(PROBLEM, Y, calc_second_order=True, print_to_console=False)

    mode = ("State-dependent loss model"
            if use_state_dependent_losses
            else "Constant parasitic-loss model")
    lines = []
    lines.append(f"Sobol sensitivity analysis: AMR electrical COP at "
                  f"T_cold={T_COLD_K}K, span={SPAN_K}K -- {mode}")
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
    if not use_state_dependent_losses:
        lines.append("Interpretation:")
        lines.append("mu0H_max, frequency and fluid_mdot show ~0 sensitivity because the")
        lines.append("current amr_cycle.py model structure makes COP_electrical algebraically")
        lines.append("independent of Qc: COP_electrical = 1 / [(Th/Tc-1)/eta_2nd_law(eps) +")
        lines.append("parasitic_fraction], since both W_mag and W_parasitic scale linearly")
        lines.append("with Qc and cancel. Field/frequency/flow rate DO change how much cooling")
        lines.append("you get (Qc, via dTad) but NOT, in this model, how efficiently you get it.")
        lines.append("Real AMR systems don't have this decoupling -- eddy-current losses scale")
        lines.append("with frequency^2, viscous dissipation with mdot^2, and magnet support/")
        lines.append("motor sizing with field -- so eta_2nd_law and parasitic_fraction should")
        lines.append("be state-dependent functions rather than constants. The")
        lines.append("state-dependent formulation implemented in loss_model.py")
        lines.append("addresses this limitation.")
    else:
        lines.append("Interpretation: with the state-dependent loss model (eddy ~ f²H²,")
        lines.append("pump ~ mdot^2, base overhead ~ Qc) in place of the constant")
        lines.append("parasitic_fraction, mu0H_max/frequency/fluid_mdot now carry real")
        lines.append("sensitivity in COP_electrical, because raising frequency or field no")
        lines.append("longer only raises Qc for free -- it also raises the eddy-current loss")
        lines.append("term, producing a genuine efficiency-versus-capacity trade-off")
        lines.append("that can be explored during multi-objective optimization.")
        lines.append("Caveat: the loss coefficients are calibrated from a small set of")
        lines.append("experimental systems (see loss_model.py). The qualitative trends")
        lines.append("are physically motivated, while the exact sensitivity magnitudes")
        lines.append("should be interpreted in light of the available calibration data.")

    with open(out_path, "w") as f:
        f.write("\n".join(lines) + "\n")
    print("\n".join(lines))
    print(f"\nWrote {out_path}")
    return Si


if __name__ == "__main__":
    print("=" * 100)
    print("Constant parasitic-loss model")
    print("=" * 100)
    run_sobol(out_path="results/sobol_results_phase2_constant.txt",
              use_state_dependent_losses=False)
    print("\n" + "=" * 100)
    print("State-dependent loss model")
    print("=" * 100)
    run_sobol(out_path="results/sobol_results.txt",
              use_state_dependent_losses=True)
