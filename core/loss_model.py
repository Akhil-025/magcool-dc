"""
loss_model.py
==============
Phase 3: replaces the constant `eta_2nd_law` / `parasitic_fraction` in
amr_cycle.py (Phase 1-2) with physically-structured, state-dependent loss
terms, motivated directly by the Sobol sensitivity finding in
results/sobol_results.txt (COP was ~99.9% sensitive to a single constant
because field/frequency/flow had no mechanism to affect efficiency).

Loss mechanisms modeled (functional forms are standard scaling laws; the
three coefficients are FIT, not derived from first principles):

    W_eddy  = k_eddy  * f^2 * (mu0*H)^2     (eddy-current loss in the
                                                magnet/regenerator support
                                                structure scales with dB/dt
                                                squared, i.e. f^2*H^2 -
                                                standard result, e.g.
                                                Kitanovski et al. 2015 Ch. 6)
    W_pump  = k_pump  * mdot^2                 (pumping power for regenerator
                                                bed flow; DeltaP ~ mdot for
                                                porous-media/Darcy flow at
                                                these Reynolds numbers, so
                                                pump power = mdot*DeltaP/rho
                                                ~ mdot^2)
    W_base  = base_frac * Qc                   (baseline controls/inverter/
                                                bearing overhead, taken as
                                                proportional to duty)

    W_parasitic = W_eddy + W_pump + W_base

CALIBRATION: the three coefficients (k_eddy, k_pump, base_frac) are fit by
solving an exactly-determined 3x3 linear system against the three Phase 2
benchmark devices (Astronautics, DTU, Tusek), using their calibrated mdot
(from validation_system.py) and the field/frequency each device actually
reports.

**Phase 6 update, attempted**: added a 4th device (Okamura & Hirano 2013),
found via a targeted literature search this session, to convert this to an
over-determined 4-point fit. **Result: this made things WORSE, not better,
and that is itself the Phase 6 finding.** Least-squares over the 4 points
gives NEGATIVE k_pump and base_frac (unphysical), and leave-one-out
cross-validation shows the fit does not generalize (errors up to +1639%
when Tusek's low-power lab device is predicted from the other three, which
span 6.5W to 2502W of cooling capacity -- four orders of magnitude
different scales/designs/vintages being pooled into one linear model is
apparently not a safe assumption). This is flagged loudly, not hidden:
`StateDependentLossModel()`'s default still calibrates on the original
3-point CORE set (which IS well-behaved, see Phase 3), and the 4-point
EXTENDED fit is exposed separately via `run_extended_diagnostic()` for
transparency, not used as the production default. A 5th candidate
(Risoe/DTU 2011, 30K span) couldn't even be calibrated: at its reported
1.1T field, cooling capacity formula cannot reach positive Qc at a 30K
span -- consistent with, and further evidence for, the single-stage
span-collapse finding from Phase 1.
"""

import numpy as np

# (device, f_Hz, mu0H_T, mdot_kg_s, Qc_W, W_parasitic_required_W)
# CORE: the stable, well-behaved 3-point set used as the production default.
CALIBRATION_POINTS_CORE = [
    ("Astronautics_rotary_2014", 4.0, 1.44, 1.0854, 2502.0, 0.453 * 2502.0),
    ("DTU_rotary_Gd_2016", 1.4, 1.44, 0.3251, 818.0, 0.171 * 818.0),
    ("Tusek_singlebed_Gd_2010", 0.25, 1.69, 0.0045, 6.5, 0.118 * 6.5),
]
# EXTENDED: CORE + Okamura & Hirano (2013). Diagnostic only -- see the
# module docstring's Phase 6 finding. NOT used by default.
CALIBRATION_POINTS_EXTENDED = CALIBRATION_POINTS_CORE + [
    # frequency not reported in the secondary source for this device --
    # 1.0 Hz placeholder (see data/amr_experimental_benchmarks.csv note)
    ("Okamura_Hirano_2013", 1.0, 1.1, 0.0502, 200.0, 0.367 * 200.0),
]
CALIBRATION_POINTS = CALIBRATION_POINTS_CORE  # backward-compat alias


def _build_system(points):
    A = np.zeros((len(points), 3))
    b = np.zeros(len(points))
    for i, (name, f, H, mdot, Qc, Wp) in enumerate(points):
        A[i, 0] = f ** 2 * H ** 2
        A[i, 1] = mdot ** 2
        A[i, 2] = Qc
        b[i] = Wp
    return A, b


def leave_one_out_cv(points=None, verbose=True):
    """Fit on N-1 points, predict the held-out point's W_parasitic, repeat
    for each point. Reports absolute and percent error per held-out device."""
    points = points if points is not None else CALIBRATION_POINTS_CORE
    results = []
    for i in range(len(points)):
        train = [p for j, p in enumerate(points) if j != i]
        test = points[i]
        A, b = _build_system(train)
        coeffs, *_ = np.linalg.lstsq(A, b, rcond=None)
        k_eddy, k_pump, base_frac = np.clip(coeffs, 0, None)
        name, f, H, mdot, Qc, Wp_true = test
        Wp_pred = k_eddy * f ** 2 * H ** 2 + k_pump * mdot ** 2 + base_frac * Qc
        err_pct = 100 * (Wp_pred - Wp_true) / Wp_true if Wp_true != 0 else float("nan")
        results.append((name, Wp_true, Wp_pred, err_pct))
        if verbose:
            print(f"  held out {name:<28} W_parasitic true={Wp_true:8.2f}W  "
                  f"predicted={Wp_pred:8.2f}W  error={err_pct:+7.1f}%")
    return results


def calibrate_loss_coefficients(points=None, verbose=True, label="CORE (production default)"):
    points = points if points is not None else CALIBRATION_POINTS_CORE
    A, b = _build_system(points)
    coeffs, residuals, rank, sv = np.linalg.lstsq(A, b, rcond=None)
    k_eddy, k_pump, base_frac = coeffs
    pred = A @ coeffs
    if verbose:
        print(f"Calibrated loss-model coefficients [{label}], "
              f"{len(points)} points, 3 unknowns "
              f"({'exactly-determined' if len(points) == 3 else 'over-determined'}):")
        print(f"  k_eddy    = {k_eddy: .6f}  W / (Hz^2 * T^2)")
        print(f"  k_pump    = {k_pump: .6f}  W / (kg/s)^2")
        print(f"  base_frac = {base_frac: .6f}  (dimensionless, x Qc)")
        print("  Fit residuals (predicted - required W_parasitic):")
        for (name, *_, Wp_true), Wp_pred in zip(points, pred):
            print(f"    {name:<28} true={Wp_true:8.2f}W  fit={Wp_pred:8.2f}W  "
                  f"resid={Wp_pred - Wp_true:+7.2f}W")
        for c, name in zip(coeffs, ["k_eddy", "k_pump", "base_frac"]):
            if c < 0:
                print(f"  WARNING: {name} came out negative ({c:.4f}) - "
                      "unphysical for a loss coefficient; clip to 0 in production use.")
    return {"k_eddy": max(k_eddy, 0.0), "k_pump": max(k_pump, 0.0),
            "base_frac": max(base_frac, 0.0), "raw": coeffs}


def run_extended_diagnostic():
    """Phase 6: shows the 4-point EXTENDED fit's instability explicitly,
    for the record -- NOT used to set the production default."""
    print("=" * 90)
    print("PHASE 6 DIAGNOSTIC: adding Okamura & Hirano (2013) as a 4th point")
    print("=" * 90)
    calibrate_loss_coefficients(CALIBRATION_POINTS_EXTENDED, verbose=True,
                                  label="EXTENDED (4pt, diagnostic only)")
    print("\n  Leave-one-out cross-validation on the EXTENDED set:")
    leave_one_out_cv(CALIBRATION_POINTS_EXTENDED, verbose=True)
    print("\n  CONCLUSION: negative coefficients + leave-one-out errors up to "
          "+1639% show a single linear loss model does not generalize across "
          "devices spanning 6.5W to 2502W of cooling capacity. The CORE "
          "3-point fit remains the production default. Pooling heterogeneous")
    print("  device scales into one fit is the wrong move without a size/scale "
          "term in the model -- a concrete Phase 7 item, not a data problem.")


class StateDependentLossModel:
    def __init__(self, k_eddy=None, k_pump=None, base_frac=None):
        if k_eddy is None or k_pump is None or base_frac is None:
            cal = calibrate_loss_coefficients(verbose=False)
            k_eddy = k_eddy if k_eddy is not None else cal["k_eddy"]
            k_pump = k_pump if k_pump is not None else cal["k_pump"]
            base_frac = base_frac if base_frac is not None else cal["base_frac"]
        self.k_eddy = k_eddy
        self.k_pump = k_pump
        self.base_frac = base_frac

    def parasitic_power(self, frequency, mu0H, mdot, Qc):
        W_eddy = self.k_eddy * frequency ** 2 * mu0H ** 2
        W_pump = self.k_pump * mdot ** 2
        W_base = self.base_frac * Qc
        return W_eddy + W_pump + W_base


if __name__ == "__main__":
    calibrate_loss_coefficients()
    print()
    run_extended_diagnostic()
