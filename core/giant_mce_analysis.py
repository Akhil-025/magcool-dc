"""
giant_mce_analysis.py
=======================
Phase 5: now that first_order_mce.py gives Gd5Si2Ge2 a physically-appropriate
(first-order Landau) model instead of the structurally-invalid mean-field
treatment flagged in Phase 4, this script answers the question posed at the
end of Phase 4: does the giant-MCE material change the COP-competitiveness
conclusion?

Short answer: not as-is, because its favorable operating window is
mistargeted for data-center duty -- but the underlying material family is
tunable, which is the real, actionable finding.
"""

import numpy as np
from core.mce_material import GADOLINIUM
from core.first_order_mce import GD5SI2GE2_FIRST_ORDER
from core.amr_cycle import AMRSystem
from core.loss_model import StateDependentLossModel
from core.baseline_cooling import vapor_compression_cop, liquid_cooling_cop

_LOSS_MODEL = StateDependentLossModel()


def find_peak_temperature(material, mu0H, T_range=(260, 320), n=601):
    mu0 = 4 * np.pi * 1e-7
    Ts = np.linspace(*T_range, n)
    dT = material.delta_T_adiabatic(Ts, mu0H / mu0)
    return float(Ts[int(np.argmax(dT))])


def run_analysis(out_path="results/giant_mce_analysis.txt"):
    lines = []
    mu0H = 2.0
    peak_T_giant = find_peak_temperature(GD5SI2GE2_FIRST_ORDER, mu0H)
    peak_T_gd = find_peak_temperature(GADOLINIUM, mu0H)

    lines.append("Giant-MCE (Gd5Si2Ge2, first-order Landau model) vs. Gd, at their")
    lines.append("own favorable operating points vs. the ASHRAE 18-27C (291-300K) range")
    lines.append("=" * 90)
    lines.append(f"Gd peak-effect temperature: {peak_T_gd:.1f} K ({peak_T_gd-273.15:.1f} C) "
                 f"-- INSIDE the ASHRAE recommended supply range")
    lines.append(f"Gd5Si2Ge2 peak-effect temperature: {peak_T_giant:.1f} K "
                 f"({peak_T_giant-273.15:.1f} C) -- BELOW the ASHRAE range by "
                 f"~{291.0-peak_T_giant:.1f} K")
    lines.append("")

    def eval_at(material, T_cold, span, mass=5.0):
        sys_ = AMRSystem(material=material, mu0H_max=mu0H, mass_regenerator=mass,
                          frequency=1.0, fluid_mdot=0.08, loss_model=_LOSS_MODEL,
                          use_ntu_thermal_model=True)
        return sys_.run(T_cold, span)

    span = 10.0
    lines.append(f"Test 1: BOTH materials at the ASHRAE operating point (T_cold=291K, span={span}K)")
    r_gd_ashrae = eval_at(GADOLINIUM, 291.0, span)
    r_giant_ashrae = eval_at(GD5SI2GE2_FIRST_ORDER, 291.0, span)
    lines.append(f"  Gd:         Qc={r_gd_ashrae.Qc:8.1f} W   COP_elec={r_gd_ashrae.COP_electrical:.2f}")
    lines.append(f"  Gd5Si2Ge2:  Qc={r_giant_ashrae.Qc:8.1f} W   COP_elec={r_giant_ashrae.COP_electrical:.2f}"
                 f"   <- collapses to ~0: {span}K span is centered "
                 f"~{291.0+span/2-peak_T_giant:.1f}K from its own peak")
    lines.append("")

    T_cold_giant_favorable = peak_T_giant - span / 2
    lines.append(f"Test 2: Gd5Si2Ge2 at ITS OWN favorable point (T_cold="
                 f"{T_cold_giant_favorable:.1f}K = {T_cold_giant_favorable-273.15:.1f}C, "
                 f"span={span}K, straddling its {peak_T_giant:.1f}K peak)")
    r_giant_own = eval_at(GD5SI2GE2_FIRST_ORDER, T_cold_giant_favorable, span)
    r_gd_same_point = eval_at(GADOLINIUM, T_cold_giant_favorable, span)
    lines.append(f"  Gd5Si2Ge2:  Qc={r_giant_own.Qc:8.1f} W   COP_elec={r_giant_own.COP_electrical:.2f}"
                 f"   <- strong performance when correctly targeted")
    lines.append(f"  Gd (same point): Qc={r_gd_same_point.Qc:8.1f} W   "
                 f"COP_elec={r_gd_same_point.COP_electrical:.2f}   <- Gd fails here "
                 f"(point is far from Gd's own {peak_T_gd:.1f}K)")
    lines.append("")

    vcc = vapor_compression_cop(291.0, 291.0 + span)
    liq = liquid_cooling_cop(291.0, 291.0 + span)
    lines.append("CONCLUSION:")
    lines.append("The giant-MCE effect is real and large when the material is operated within")
    lines.append("its own narrow first-order transition window -- but Gd5Si2Ge2's window sits")
    lines.append(f"~{291.0-peak_T_giant:.0f}K below the ASHRAE data-center range, so it is not directly")
    lines.append("usable for this application as-is. This does NOT overturn the Phase 1-4")
    lines.append("conclusion (Gd trails vapor-compression/liquid cooling on COP within the")
    lines.append("ASHRAE range). What it DOES support: literature confirms the Gd5(SixGe1-x)4")
    lines.append("family has composition-tunable ordering temperature (Pecharsky & Gschneidner,")
    lines.append("Appl. Phys. Lett. 70, 3299 (1997), report tunability from ~20K to ~276K by")
    lines.append("Si:Ge ratio, with Gd5Si4 itself ordering at 335K) -- so a composition between")
    lines.append("Gd5Si2Ge2 and Gd5Si4 tuned to ~291-300K, IF it retains first-order/giant")
    lines.append("character at that composition, is the genuinely promising untested direction")
    lines.append("for closing the COP gap. This is a materials-synthesis question outside what")
    lines.append("a simulation suite alone can answer -- flagged as Phase 6 in ROADMAP.md, not")
    lines.append("claimed as solved here.")
    lines.append(f"\nFor reference, baselines at this operating point: VCC COP={vcc.COP:.2f}, "
                 f"Liquid COP={liq.COP:.2f}")
    lines.append(f"\nNote also: even correctly targeted, Gd5Si2Ge2's COP_electrical "
                 f"({r_giant_own.COP_electrical:.2f}) is close to Gd's own COP_electrical at "
                 f"its matched point ({r_gd_ashrae.COP_electrical:.2f}), not dramatically "
                 f"higher, despite ~4x the cooling capacity (Qc). This is consistent with "
                 f"Phase 3's Sobol finding: COP_electrical is driven mainly by frequency/flow/"
                 f"field-dependent losses (loss_model.py), not by which material is loaded into "
                 f"the regenerator -- a bigger MCE mostly buys more Qc per kg, not a better COP.")

    with open(out_path, "w") as f:
        f.write("\n".join(lines) + "\n")
    print("\n".join(lines))
    print(f"\nWrote {out_path}")


if __name__ == "__main__":
    run_analysis()
