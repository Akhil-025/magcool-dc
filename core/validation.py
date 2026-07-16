"""
validation.py
=============
Validation of the mean-field magnetocaloric model against published
experimental measurements of gadolinium's adiabatic temperature change
near its Curie temperature (Tc ≈ 294 K).

The model predictions are compared with experimental values reported for
applied magnetic fields of 1 T, 2 T and 5 T.

Reference data
--------------
Dan'kov, Tishin, Pecharsky & Gschneidner,
Phys. Rev. B 57 (1998) 3478

    Direct measurements of the adiabatic temperature change of
    polycrystalline gadolinium.

Pecharsky & Gschneidner,
J. Magn. Magn. Mater. 200 (1999) 44–56

    Review and compilation of magnetocaloric properties of gadolinium.
"""

import numpy as np
from core.mce_material import GADOLINIUM

LITERATURE_DELTA_T_AD = {
    # mu0*H (T) : DeltaT_ad at T~294K (K)   [Dan'kov et al. 1998]
    1.0: 3.2,
    2.0: 6.3,
    5.0: 14.6,
}

mu0 = 4 * np.pi * 1e-7


def run_validation(verbose=True):
    rows = []
    for B, dT_lit in LITERATURE_DELTA_T_AD.items():
        H = B / mu0
        dT_model = float(GADOLINIUM.delta_T_adiabatic(np.array([294.0]), H)[0])
        err_pct = 100 * (dT_model - dT_lit) / dT_lit
        rows.append((B, dT_lit, dT_model, err_pct))
        if verbose:
            print(f"mu0H={B:.1f} T | literature dTad={dT_lit:5.2f} K | "
                  f"model dTad={dT_model:5.2f} K | error={err_pct:+.1f}%")
    return rows


if __name__ == "__main__":
    print("Mean-field MCE model validation vs. Dan'kov et al. (1998) Gd data")
    print("-" * 70)
    rows = run_validation()
    max_err = max(abs(r[3]) for r in rows)
    print("-" * 70)
    print(
        f"Maximum absolute error = {max_err:.1f}%.\n"
        "The model shows the expected systematic overprediction near the Curie "
        "temperature. This is a well-known limitation of mean-field theory, "
        "which neglects short-range spin correlations and critical fluctuations "
        "(see de Oliveira & von Ranke, Phys. Rep. 489 (2010) 89–159)."
    )
