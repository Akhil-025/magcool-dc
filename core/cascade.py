"""
cascade.py
===========
Multi-stage (cascade) AMR model for extending the operating temperature
span beyond the limits of a single-stage system.

Cascade concept: N single-stage AMR modules in series, each handling an
equal share of the total span (Th_total - Tc_total)/N, analogous to
cascade vapor-compression refrigeration. Stage 1 (coldest) absorbs the
data-center heat load Qc at T_cold; its heat rejection becomes the input
to Stage 2, and so on until the final stage rejects heat at T_hot.
Because heat is transferred through the stages in series, each stage must
reject the same cooling load (steady state, neglecting inter-stage losses):

    W_total = Σ W_i(Qc, span/N)
    COP_cascade = Qc / W_total

This simplified model assumes identical gadolinium regenerator stages.
More advanced cascade systems may employ graded Curie-temperature
materials (e.g. Gd alloys) in different stages, but those effects are not
included here.
"""

import numpy as np
import csv
from core.mce_material import GADOLINIUM
from core.amr_cycle import AMRSystem
from core.loss_model import StateDependentLossModel
from core.baseline_cooling import vapor_compression_cop, liquid_cooling_cop

_LOSS_MODEL = StateDependentLossModel()
USE_NTU_THERMAL_MODEL = True


def run_cascade(T_cold_K, total_span_K, n_stages, material=None, mu0H_max=2.0,
                 mass_per_stage=2.0, frequency=1.0, fluid_mdot=0.08,
                 regenerator_effectiveness=0.85):
    """Runs n_stages identical AMR modules in series, each covering
    total_span_K/n_stages, all passing the same Qc through in steady state
    (Qc is set by the coldest/first stage's capacity at its local span)."""
    if material is None:
        material = GADOLINIUM
    span_per_stage = total_span_K / n_stages
    T_local = T_cold_K
    # First stage sets the deliverable Qc (bottleneck of the chain)
    stage1 = AMRSystem(material=material, mu0H_max=mu0H_max,
                        mass_regenerator=mass_per_stage, frequency=frequency,
                        fluid_mdot=fluid_mdot, regenerator_effectiveness=regenerator_effectiveness,
                        loss_model=_LOSS_MODEL, use_ntu_thermal_model=USE_NTU_THERMAL_MODEL)
    r1 = stage1.run(T_local, span_per_stage)
    Qc_target = r1.Qc
    if Qc_target <= 0:
        return {"n_stages": n_stages, "span_per_stage_K": span_per_stage,
                "Qc_W": 0.0, "W_total_W": np.nan, "COP_cascade": 0.0,
                "feasible": False}

    W_total = 0.0
    for i in range(n_stages):
        stage = AMRSystem(material=material, mu0H_max=mu0H_max,
                           mass_regenerator=mass_per_stage, frequency=frequency,
                           fluid_mdot=fluid_mdot, regenerator_effectiveness=regenerator_effectiveness,
                           loss_model=_LOSS_MODEL, use_ntu_thermal_model=USE_NTU_THERMAL_MODEL)
        # each stage handles the same Qc_target at its local span; back out
        # the required work by re-running at span_per_stage and scaling mdot
        # if needed so Qc matches Qc_target (steady-state series constraint)
        r_i = stage.run(T_local, span_per_stage)
        if r_i.Qc > 0:
            scale = Qc_target / r_i.Qc
            W_i = (r_i.W_mag + r_i.W_parasitic) * scale
        else:
            W_i = np.inf
        W_total += W_i
        T_local += span_per_stage

    COP_cascade = Qc_target / W_total if W_total > 0 else 0.0
    return {"n_stages": n_stages, "span_per_stage_K": span_per_stage,
            "Qc_W": round(Qc_target, 1), "W_total_W": round(W_total, 1),
            "COP_cascade": round(COP_cascade, 3), "feasible": True}


def compare_staging(T_cold_C=18.0, spans=range(5, 21), stage_counts=(1, 2, 3, 4),
                     material=None, mass_per_stage=2.0,
                     out_csv="results/cascade_comparison.csv"):
    T_cold_K = T_cold_C + 273.15
    rows = []
    for span in spans:
        T_hot_K = T_cold_K + span
        vcc = vapor_compression_cop(T_cold_K, T_hot_K)
        liq = liquid_cooling_cop(T_cold_K, T_hot_K)
        row = {"span_K": span, "VaporCompression_COP": round(vcc.COP, 2),
               "LiquidCooling_COP": round(liq.COP, 2)}
        for n in stage_counts:
            res = run_cascade(T_cold_K, span, n, material=material, mass_per_stage=mass_per_stage)
            row[f"AMR_{n}stage_COP"] = res["COP_cascade"] if res["feasible"] else None
            row[f"AMR_{n}stage_Qc_W"] = res["Qc_W"] if res["feasible"] else None
        rows.append(row)

    with open(out_csv, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=rows[0].keys())
        writer.writeheader()
        writer.writerows(rows)
    return rows


if __name__ == "__main__":
    from core.mce_material import GD5SI2GE2

    print("Cascade AMR staging vs. baselines, ASHRAE 5-20K span sweep")
    print("(mu0H=2T per stage, 5kg regenerator per stage, f=1Hz, mdot=0.08kg/s, NTU thermal model on)")
    print("=" * 100)
    print("\n--- Material: Gd (baseline) ---")
    rows_gd = compare_staging(material=GADOLINIUM, mass_per_stage=5.0,
                                out_csv="results/cascade_comparison.csv")
    header = f"{'span':>5} {'1-stage':>9} {'2-stage':>9} {'3-stage':>9} {'4-stage':>9} {'VCC':>7} {'Liquid':>7}"
    print(header)
    for r in rows_gd:
        def fmt(v):
            return f"{v:9.2f}" if v is not None else f"{'--':>9}"
        print(f"{r['span_K']:>5} {fmt(r['AMR_1stage_COP'])} {fmt(r['AMR_2stage_COP'])} "
              f"{fmt(r['AMR_3stage_COP'])} {fmt(r['AMR_4stage_COP'])} "
              f"{r['VaporCompression_COP']:>7} {r['LiquidCooling_COP']:>7}")
    print(f"Wrote results/cascade_comparison.csv")

    print("\n--- Material: Gd5Si2Ge2 (giant MCE) ---")
    rows_giant = compare_staging(material=GD5SI2GE2, mass_per_stage=5.0,
                                   out_csv="results/cascade_comparison_giant_mce.csv")
    print(header)
    for r in rows_giant:
        print(f"{r['span_K']:>5} {fmt(r['AMR_1stage_COP'])} {fmt(r['AMR_2stage_COP'])} "
              f"{fmt(r['AMR_3stage_COP'])} {fmt(r['AMR_4stage_COP'])} "
              f"{r['VaporCompression_COP']:>7} {r['LiquidCooling_COP']:>7}")
    print(f"Wrote results/cascade_comparison_giant_mce.csv")

    gd_10K = next(r for r in rows_gd if r["span_K"] == 10)
    giant_10K = next(r for r in rows_giant if r["span_K"] == 10)
    print(f"\nAt 10K span: Gd 1-stage COP={gd_10K['AMR_1stage_COP']} vs. "
          f"Gd5Si2Ge2 1-stage COP={giant_10K['AMR_1stage_COP']} "
          f"(VCC={gd_10K['VaporCompression_COP']}, Liquid={gd_10K['LiquidCooling_COP']})")
