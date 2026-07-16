"""
loss_model.py
=============
Implements a state-dependent parasitic loss model for the AMR system,
replacing the constant parasitic-power assumption with losses that depend
on operating conditions. The model was motivated by sensitivity analysis,
which showed that a constant parasitic term made electrical COP almost
entirely insensitive to magnetic field, frequency, and fluid flow rate.

Loss mechanisms modeled (functional forms are standard scaling laws; the
three coefficients are fitted, not derived from first principles):

    W_eddy  = k_eddy * f² * (mu0H)²
        Eddy-current losses in the magnet/regenerator support structure
        scale with (dB/dt)² and therefore approximately with f²H²
        (Kitanovski et al., 2015, Ch. 6).

    W_pump  = k_pump * mdot²
        Pumping power for porous-media flow through the regenerator.
        Assuming Darcy-flow behaviour, ΔP ∝ mdot and therefore
        pump power ∝ mdot².

    W_base  = base_frac * Qc
        Baseline electrical overhead (controls, inverter, bearings),
        assumed proportional to cooling duty.

    W_parasitic = W_eddy + W_pump + W_base

Calibration
-----------
The coefficients (k_eddy, k_pump, base_frac) are obtained by solving a
3×3 linear system using three experimentally reported AMR devices with
published operating conditions and calibrated flow rates.

An additional four-point least-squares fit including the Okamura &
Hirano (2013) device was investigated. Although adding more data might
normally improve robustness, the resulting fit produced negative loss
coefficients and poor leave-one-out cross-validation performance,
indicating that a single linear model does not generalize across devices
with widely different sizes and operating conditions. For this reason,
the original three-point calibration remains the default, while the
extended fit is retained only as a diagnostic example.

A fifth candidate dataset (Risø/DTU 2011, 30 K span) could not be
calibrated because the corresponding operating point does not yield
positive cooling capacity under the present AMR model.
"""

import numpy as np

# (device, f_Hz, mu0H_T, mdot_kg_s, Qc_W, W_parasitic_required_W)
# CORE: the stable, well-behaved 3-point set used as the production default.
CALIBRATION_POINTS_CORE = [
    ("Astronautics_rotary_2014", 4.0, 1.44, 1.0854, 2502.0, 0.453 * 2502.0),
    ("DTU_rotary_Gd_2016", 1.4, 1.44, 0.3251, 818.0, 0.171 * 818.0),
    ("Tusek_singlebed_Gd_2010", 0.25, 1.69, 0.0045, 6.5, 0.118 * 6.5),
]
# EXTENDED: CORE + Okamura & Hirano (2013). Retained only as a
# diagnostic comparison; not used as the default calibration.
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
    """Demonstrates the instability of the four-point extended fit.
    The diagnostic is provided for transparency and is not used as the
    production calibration."""
    print("=" * 90)
    print("DIAGNOSTIC: adding Okamura & Hirano (2013) as a fourth calibration point")
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
          "term in the model. A size-dependent loss model would be a logical next step.")


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
