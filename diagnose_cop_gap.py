"""
diagnose_cop_gap.py
====================
Isolates WHY the fixed-condition baseline sweep (stage 4, COP_elec=4.63
at span=10K) differs from the NSGA-III Pareto-optimal design at the same
operating point (COP_elec=9.49). Toggles one factor at a time so each
factor's individual contribution to the ~2.05x gap can be read off.
"""
import numpy as np
from core.mce_material import GADOLINIUM
from core.amr_cycle import AMRSystem
from core.loss_model import StateDependentLossModel

T_COLD_K = 291.15   # stage 4's actual T_cold (18C), close enough to optimizer's 291.0K
SPAN_K = 10.0

loss_model = StateDependentLossModel()

# Stage 4 baseline sweep's fixed point
BASE = dict(mu0H_max=2.0, mass_regenerator=5.0, frequency=2.0,
            fluid_mdot=0.08, regenerator_effectiveness=0.85)

# NSGA-III's best-electrical-COP Pareto point
OPT = dict(mu0H_max=2.993, mass_regenerator=11.55, frequency=0.322,
           fluid_mdot=0.4999, regenerator_effectiveness=0.83, blow_fraction=0.414)


def run(label, params, loss=None, ntu=False, bf=0.5):
    p = dict(params)
    p.pop("regenerator_effectiveness", None) if ntu else None
    sys_ = AMRSystem(material=GADOLINIUM,
                      mu0H_max=p["mu0H_max"], mass_regenerator=p["mass_regenerator"],
                      frequency=p["frequency"], fluid_mdot=p["fluid_mdot"],
                      regenerator_effectiveness=params.get("regenerator_effectiveness", 0.85),
                      loss_model=loss, use_ntu_thermal_model=ntu, blow_fraction=bf)
    r = sys_.run(T_COLD_K, SPAN_K)
    print(f"{label:60s} COP_elec={r.COP_electrical:6.3f}  Qc={r.Qc:9.1f}W  "
          f"W_mag={r.W_mag:7.1f}W  W_parasitic={r.W_parasitic:7.1f}W  eps_used={sys_._effective_eps():.3f}")
    return r.COP_electrical


print(f"Operating point: T_cold={T_COLD_K}K, span={SPAN_K}K\n")
print("--- Starting point (stage 4 baseline sweep, exact) ---")
c0 = run("BASE: f=2Hz mdot=0.08 mass=5kg H=2T eps=0.85(const), constant parasitic=0.15*Qc, bf=0.5",
          BASE, loss=None, ntu=False, bf=0.5)

print("\n--- Step 1: swap ONLY the loss model (constant 0.15*Qc -> calibrated state-dependent) ---")
c1 = run("BASE params, but state-dependent loss model", BASE, loss=loss_model, ntu=False, bf=0.5)

print("\n--- Step 2: also swap constant eps=0.85 -> NTU-computed eps (still BASE mass/f/mdot) ---")
c2 = run("BASE params, state-dependent loss + NTU eps", BASE, loss=loss_model, ntu=True, bf=0.5)

print("\n--- Step 3: also move to the optimizer's field/frequency/mdot/mass ---")
c3 = run("OPT params (H,f,mdot,mass), state-dependent loss + NTU eps, bf=0.5",
          OPT, loss=loss_model, ntu=True, bf=0.5)

print("\n--- Step 4: also apply the optimizer's blow_fraction=0.414 (full reproduction) ---")
c4 = run("OPT params, state-dependent loss + NTU eps + bf=0.414 (should match NSGA-III's 9.494)",
          OPT, loss=loss_model, ntu=True, bf=0.414)

print(f"\nSummary of cumulative COP_elec: {c0:.2f} -> {c1:.2f} -> {c2:.2f} -> {c3:.2f} -> {c4:.2f}")
print(f"Total gap reproduced: {c4/c0:.2f}x (NSGA-III reported 9.494/4.63 = {9.494/4.63:.2f}x)")
print("\nPer-factor multipliers:")
print(f"  loss model swap alone:        {c1/c0:.2f}x")
print(f"  + NTU effectiveness:          {c2/c1:.2f}x")
print(f"  + move to optimizer's H/f/mdot/mass: {c3/c2:.2f}x")
print(f"  + optimizer's blow_fraction:  {c4/c3:.2f}x")