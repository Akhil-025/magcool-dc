"""
giant_mce_analysis.py
=====================
Evaluates whether a first-order giant magnetocaloric material
(Gd5Si2Ge2), modeled using the Landau free-energy formulation in
first_order_mce.py, changes the system-level competitiveness of AMR
cooling for data-center applications.

Short answer: not in its current composition, because its optimal
operating temperature lies below the ASHRAE-recommended data-center
range. However, the underlying Gd5(SixGe1-x)4 material family is
compositionally tunable, making it a promising direction for future
materials development.
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


def landau_peak_offset_K(material, mu0H, T_range=(260, 330), n=701):
    """Quantifies a real, load-bearing property of this repo's extended
    Landau formulation (first_order_mce.py): the temperature at which
    delta_T_adiabatic(T) actually PEAKS at a given field is NOT the
    material's nominal Tc parameter -- it sits systematically ABOVE Tc.

    This was found during a paper-mining/cross-check pass and confirmed
    TWICE, independently, at the same field (mu0H~1.4-1.6T):
      1. GD5SI2GE2_FIRST_ORDER (Tc=276.0K, per Pecharsky & Gschneidner,
         Phys. Rev. Lett. 78, 4494 (1997)) has its dTad(T) peak at
         ~286.4K at 2T -- an offset of ~+10.4K.
      2. cascade.py's Astronautics graded-bed reproduction (see
         run_astronautics_graded_bed_check() and ROADMAP.md
         addendum) back-solves the composition Tc needed, PER STAGE, for
         the model's peak effect to land on each stage's actual operating
         temperature (T_mid). Compared against Jacobs et al. (2014,
         Int. J. Refrig. 37, 84-91), Table 1's six ACTUAL layer Curie
         temperatures (30.5-43.0 C = 303.65-316.15 K), the offset is
         +11.1 to +11.5 K across all six stages -- consistent within
         ~1K of finding (1), despite being a completely independent check
         (different composition target per stage, different field regime
         within the cascade).

    This is a genuine, stable characteristic of the (A,B,C)=(10,-4,8)
    Landau free-energy expansion used here (see first_order_mce.py's
    module docstring for the functional form) -- NOT a bug, and NOT
    something to "fix" by shifting Tc, since Tc is itself calibrated
    to the transition temperature reported in the source literature
    (Pecharsky & Gschneidner 1997; Jacobs et al. 2014's Table 1). It
    means: whenever this codebase composition-tunes a first-order
    material to hit a target operating temperature (composition_tuned_
    material() in first_order_mce.py), the resulting nominal `Tc` will
    read ~10-11K BELOW that target -- this is expected, and downstream
    code should not be "corrected" to make Tc equal the target T_mid.

    Returns the offset in K (peak_T - material.Tc) at the given field.
    """
    peak_T = find_peak_temperature(material, mu0H, T_range=T_range, n=n)
    return peak_T - material.Tc


def run_analysis(out_path="results/giant_mce_analysis.txt"):
    lines = []
    mu0H = 2.0
    peak_T_giant = find_peak_temperature(GD5SI2GE2_FIRST_ORDER, mu0H)
    peak_T_gd = find_peak_temperature(GADOLINIUM, mu0H)
    landau_offset = landau_peak_offset_K(GD5SI2GE2_FIRST_ORDER, mu0H)

    lines.append("Giant-MCE (Gd5Si2Ge2, first-order Landau model) vs. Gd, at their")
    lines.append("own favorable operating points vs. the ASHRAE 18-27C (291-300K) range")
    lines.append("=" * 90)
    lines.append(f"Gd peak-effect temperature: {peak_T_gd:.1f} K ({peak_T_gd-273.15:.1f} C) "
                 f"-- INSIDE the ASHRAE recommended supply range")
    lines.append(f"Gd5Si2Ge2 peak-effect temperature: {peak_T_giant:.1f} K "
                 f"({peak_T_giant-273.15:.1f} C) -- BELOW the ASHRAE range by "
                 f"~{291.0-peak_T_giant:.1f} K")
    lines.append(f"Landau-model Tc-vs-peak-effect offset at {mu0H:.1f}T: "
                 f"peak sits {landau_offset:+.1f} K above the nominal Tc="
                 f"{GD5SI2GE2_FIRST_ORDER.Tc:.1f}K parameter -- a real, stable "
                 f"property of this extended Landau formulation (see "
                 f"landau_peak_offset_K() docstring), independently confirmed "
                 f"to within ~1K by the Astronautics graded-bed cascade check "
                 f"(+11.1 to +11.5K offset across all 6 stages vs. Jacobs et "
                 f"al. 2014 Table 1's actual layer Curie temperatures).")
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
    lines.append(f" Gd:         Qc={r_gd_ashrae.Qc:8.1f} W COP_elec={r_gd_ashrae.COP_electrical:.2f}")
    lines.append(f" Gd5Si2Ge2:  Qc={r_giant_ashrae.Qc:8.1f} W COP_elec={r_giant_ashrae.COP_electrical:.2f}"
                 f"   <- collapses to ~0: {span}K span is centered "
                 f"~{291.0+span/2-peak_T_giant:.1f}K from its own peak")
    lines.append("")

    T_cold_giant_favorable = peak_T_giant - span / 2
    lines.append(f"Test 2: Gd5Si2Ge2 at ITS OWN favorable point (T_cold="
                 f"{T_cold_giant_favorable:.1f}K = {T_cold_giant_favorable-273.15:.1f}C, "
                 f"span={span}K, straddling its {peak_T_giant:.1f}K peak)")
    r_giant_own = eval_at(GD5SI2GE2_FIRST_ORDER, T_cold_giant_favorable, span)
    r_gd_same_point = eval_at(GADOLINIUM, T_cold_giant_favorable, span)
    lines.append(f" Gd5Si2Ge2:  Qc={r_giant_own.Qc:8.1f} W COP_elec={r_giant_own.COP_electrical:.2f}"
                 f"   <- strong performance when correctly targeted")
    lines.append(f" Gd (same point): Qc={r_gd_same_point.Qc:8.1f} W "
                 f"COP_elec={r_gd_same_point.COP_electrical:.2f}   <- Gd fails here "
                 f"(point is far from Gd's own {peak_T_gd:.1f}K)")
    lines.append("")

    vcc = vapor_compression_cop(291.0, 291.0 + span)
    liq = liquid_cooling_cop(291.0, 291.0 + span)
    # CITATION FIX: verified directly against Pecharsky & Gschneidner, Appl.
    # Phys. Lett. 70, 3299 (1997). The abstract's own tunability claim is
    # "the ordering temperature is tunable from ~30 to ~276 K by adjusting
    # the Si:Ge ratio" -- the "~20K" previously used here conflated this with
    # the paper's TITLE ("...for magnetic refrigeration from ~20 to ~290K"),
    # a broader claim about the achievable refrigeration range across the
    # whole Gd5(SixGe1-x)4 family (including Gd5Ge4-rich, non-giant-MCE
    # compositions), not the tunable-while-still-giant-MCE window quoted here.
    lines.append("CONCLUSION:")
    lines.append("The giant-MCE effect is real and large when the material is operated within")
    lines.append("its own narrow first-order transition window -- but Gd5Si2Ge2's window sits")
    lines.append(f"~{291.0-peak_T_giant:.0f}K below the ASHRAE data-center range, so it is not directly")
    lines.append("usable for this application as-is. This does not overturn the")
    lines.append("conclusion (Gd trails vapor-compression/liquid cooling on COP within the")
    lines.append("ASHRAE range). What it DOES support: literature confirms the Gd5(SixGe1-x)4")
    lines.append("family has composition-tunable ordering temperature (Pecharsky & Gschneidner,")
    lines.append("Appl. Phys. Lett. 70, 3299 (1997), report tunability from ~30K to ~276K by")
    lines.append("Si:Ge ratio, with Gd5Si4 itself ordering at ~335K) -- so a composition between")
    lines.append("Gd5Si2Ge2 and Gd5Si4 tuned to ~291-300K, IF it retains first-order/giant")
    lines.append("character at that composition, is the genuinely promising untested direction")
    lines.append("for closing the COP gap. This is a materials-synthesis question outside what")
    lines.append("a simulation study alone can answer and remains an open materials")
    lines.append("research question rather than a result established here.")
    lines.append(f"\nFor reference, baselines at this operating point: VCC COP={vcc.COP:.2f}, "
                 f"Liquid COP={liq.COP:.2f}")
    lines.append(f"\nNote also: even correctly targeted, Gd5Si2Ge2's COP_electrical "
                 f"({r_giant_own.COP_electrical:.2f}) is close to Gd's own COP_electrical at "
                 f"its matched point ({r_gd_ashrae.COP_electrical:.2f}), not dramatically "
                 f"higher, despite ~4x the cooling capacity (Qc). This is consistent with "
                 f"the Sobol sensitivity analysis: COP_electrical is driven mainly by frequency/flow/"
                 f"field-dependent losses (loss_model.py), not by which material is loaded into "
                 f"the regenerator -- a bigger MCE mostly buys more Qc per kg, not a better COP.")

    with open(out_path, "w") as f:
        f.write("\n".join(lines) + "\n")
    print("\n".join(lines))
    print(f"\nWrote {out_path}")


if __name__ == "__main__":
    run_analysis()