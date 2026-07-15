"""
main.py
=======
Runs the magnetic (AMR) cooling model against vapor-compression and liquid
cooling baselines across the ASHRAE TC9.9 data-center thermal guideline
envelope, and writes a comparison table + plot to results/.

Operating points: server return air / facility water at Tc = 18-27 C
(ASHRAE Class A1/A2 recommended range), rejecting to Th = Tc + span, where
span is swept 5-20 K to cover close-coupled to chilled-water-plant duty.
Reference: ASHRAE TC9.9, "Thermal Guidelines for Data Processing
Environments", 5th ed. (2021).
"""

import numpy as np
import csv
from core.mce_material import GADOLINIUM
from core.amr_cycle import AMRSystem
from core.baseline_cooling import vapor_compression_cop, liquid_cooling_cop

RESULTS_CSV = "results/comparison_table.csv"


def main():
    T_cold_C = 18.0          # ASHRAE recommended supply, deg C
    T_cold_K = T_cold_C + 273.15
    spans = np.arange(5, 21, 1)  # K

    amr = AMRSystem(
        material=GADOLINIUM,
        mu0H_max=2.0,               # T, representative Halbach permanent-magnet array
        mass_regenerator=5.0,        # kg
        frequency=2.0,                # Hz, typical lab/pilot AMR
        fluid_cp=4186.0,
        fluid_mdot=0.08,
        regenerator_effectiveness=0.85,
    )

    rows = []
    for span in spans:
        T_hot_K = T_cold_K + span
        amr_res = amr.run(T_cold_K, span)
        vcc = vapor_compression_cop(T_cold_K, T_hot_K)
        liq = liquid_cooling_cop(T_cold_K, T_hot_K)
        rows.append({
            "span_K": span,
            "Tc_K": T_cold_K, "Th_K": T_hot_K,
            "AMR_COP": round(amr_res.COP, 2),
            "AMR_Qc_W": round(amr_res.Qc, 1),
            "AMR_2ndlaw_eff": round(amr_res.exergy_eff, 3),
            "VaporCompression_COP": round(vcc.COP, 2),
            "LiquidCooling_COP": round(liq.COP, 2),
            "Carnot_COP": round(vcc.COP_carnot, 2),
        })

    with open(RESULTS_CSV, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=rows[0].keys())
        writer.writeheader()
        writer.writerows(rows)

    print(f"{'span(K)':>8} {'AMR COP':>9} {'VCC COP':>9} {'Liquid COP':>11} {'Carnot':>8}")
    for r in rows:
        print(f"{r['span_K']:>8} {r['AMR_COP']:>9} {r['VaporCompression_COP']:>9} "
              f"{r['LiquidCooling_COP']:>11} {r['Carnot_COP']:>8}")
    print(f"\nWrote {RESULTS_CSV}")


if __name__ == "__main__":
    main()
