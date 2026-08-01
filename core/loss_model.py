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
The coefficients (k_eddy, k_pump, base_frac) are obtained by fitting a
3-parameter linear model to experimentally reported AMR devices with
published operating conditions and calibrated flow rates. The CORE set
(3 devices, 3 unknowns) is exactly determined.

Phase 7 update: the fit is solved with non-negative least squares (NNLS,
Lawson & Hanson 1974, `scipy.optimize.nnls`) rather than unconstrained
least squares with post-hoc clipping of negative coefficients to zero.
NNLS is the correct tool once the physical constraint is "loss
coefficients cannot be negative": it solves the *constrained* problem
directly (the best non-negative fit), whereas clipping an unconstrained
solution afterwards does not minimize any well-defined objective under
that constraint and can silently mask an unstable fit. For the CORE
3-point set the unconstrained optimum is already non-negative, so NNLS
and the old lstsq-then-clip approach agree exactly (checked by
`tests/test_supporting_modules.py`). For the EXTENDED 4-point set, NNLS
removes the negative k_eddy/base_frac Phase 6 reported (`k_eddy=0` is
now the constrained optimum instead of a clipped negative value), but
leave-one-out error predicting the smallest device (Tušek, 6.5 W) from
the other three is still ~680% (down from +1639% with plain lstsq, but
still an order-of-magnitude miss) -- confirming Phase 6's conclusion
that pooling four orders of magnitude of device scale needs a structural
change to the model, not just a better-behaved solver.

`analyze_parasitic_fraction_scaling()` checks the natural next
hypothesis -- that this is a simple device-size effect, i.e. that
smaller devices carry proportionally more fixed (non-Qc-scaling)
overhead. Sorting the four devices by Qc shows this does NOT hold:
the smallest device (Tušek, 6.5 W) has one of the *lowest* parasitic
fractions (11.8%), and the largest (Astronautics, 2502 W) has the
*highest* (45.3%) -- the opposite of what a fixed-overhead/economies-
of-scale story predicts, and non-monotonic in between (Okamura,
200 W, at 36.7%; DTU, 818 W, at 17.1%). The Astronautics figure is
independently flagged by its own source paper as reflecting "mediocre"
electrical-component efficiency at that scale (a device-specific
engineering choice), not a generic size law. So a size/scale *term* --
as the Phase 6 write-up speculated -- is not supported by the data
actually in hand; what varies device-to-device looks more like
motor/inverter efficiency class and drivetrain topology (rotary vs.
reciprocating single-bed) than raw cooling capacity. Flagged as a
correction to the Phase 7 roadmap: more benchmark devices with
independently reported component efficiencies, not a size term, is the
concrete next step.

A fifth candidate dataset (Risø/DTU 2011, 30 K span) could not be
calibrated because the corresponding operating point does not yield
positive cooling capacity under the present AMR model.
"""

import numpy as np
from scipy.optimize import nnls

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
# FURTHER_EXTENDED: EXTENDED + 4 points from Lozano et al. (2016), the
# POLO/UFSC rotary device (see data/amr_experimental_benchmarks.csv and
# validation_system.py). Of the 8 reported (frequency, flow, span, Qc, COP)
# operating points, only these 4 calibrate at all within mdot in [1e-6,5]
# kg/s -- the other 4 (r1, r2, r3, r5) reported a Qc the model cannot reach
# at that span/field/mass under any flow rate, and are reported as
# calibration failures, not silently dropped (see
# validation_system.run_system_validation() output). Each Wp_required here
# comes from that device's OWN calibrated mdot: Wp = Qc * (1/COP_lit -
# 1/COP_ideal), same formula used to build CORE/EXTENDED above.
CALIBRATION_POINTS_FURTHER_EXTENDED = CALIBRATION_POINTS_EXTENDED + [
    ("Lozano_POLO_UFSC_2016_r4", 0.4, 0.88, 0.3244, 62.5, 1.684 * 62.5),
    ("Lozano_POLO_UFSC_2016_r6", 0.8, 0.88, 0.0430, 81.2, 1.505 * 81.2),
    ("Lozano_POLO_UFSC_2016_r7", 0.4, 0.88, 0.0207, 80.8, 1.291 * 80.8),
    ("Lozano_POLO_UFSC_2016_r8", 0.8, 0.88, 0.0308, 120.4, 1.180 * 120.4),
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
        coeffs, _resid_norm = nnls(A, b)
        k_eddy, k_pump, base_frac = coeffs
        name, f, H, mdot, Qc, Wp_true = test
        Wp_pred = k_eddy * f ** 2 * H ** 2 + k_pump * mdot ** 2 + base_frac * Qc
        err_pct = 100 * (Wp_pred - Wp_true) / Wp_true if Wp_true != 0 else float("nan")
        results.append((name, Wp_true, Wp_pred, err_pct))
        if verbose:
            print(f"  held out {name:<28} W_parasitic true={Wp_true:8.2f}W  "
                  f"predicted={Wp_pred:8.2f}W  error={err_pct:+7.1f}%")
    return results


def calibrate_loss_coefficients(points=None, verbose=True, label="CORE (production default)"):
    """Fit (k_eddy, k_pump, base_frac) via non-negative least squares (NNLS,
    Lawson & Hanson 1974). NNLS solves the constrained problem directly --
    the best-fitting non-negative coefficients -- rather than solving the
    unconstrained problem and clipping negative results afterwards, which
    does not minimize any well-defined objective under the non-negativity
    constraint. When the unconstrained optimum is already non-negative (as
    for the exactly-determined CORE 3-point set) NNLS returns the same
    solution as unconstrained least squares."""
    points = points if points is not None else CALIBRATION_POINTS_CORE
    A, b = _build_system(points)
    coeffs, resid_norm = nnls(A, b)
    k_eddy, k_pump, base_frac = coeffs
    pred = A @ coeffs
    if verbose:
        print(f"Calibrated loss-model coefficients [{label}], "
              f"{len(points)} points, 3 unknowns "
              f"({'exactly-determined' if len(points) == 3 else 'over-determined'}, "
              f"NNLS, residual norm={resid_norm:.4f}):")
        print(f"  k_eddy    = {k_eddy: .6f}  W / (Hz^2 * T^2)")
        print(f"  k_pump    = {k_pump: .6f}  W / (kg/s)^2")
        print(f"  base_frac = {base_frac: .6f}  (dimensionless, x Qc)")
        print("  Fit residuals (predicted - required W_parasitic):")
        for (name, *_, Wp_true), Wp_pred in zip(points, pred):
            print(f"    {name:<28} true={Wp_true:8.2f}W  fit={Wp_pred:8.2f}W  "
                  f"resid={Wp_pred - Wp_true:+7.2f}W")
        for c, name in zip(coeffs, ["k_eddy", "k_pump", "base_frac"]):
            if c == 0.0:
                print(f"  NOTE: {name} pinned to 0 by the non-negativity "
                      "constraint (the unconstrained optimum was negative here).")
    return {"k_eddy": k_eddy, "k_pump": k_pump,
            "base_frac": base_frac, "raw": coeffs}


def analyze_parasitic_fraction_scaling(points=None, verbose=True):
    """Tests the hypothesis (raised in the Phase 6 write-up) that the loss
    model's cross-device instability is a simple device-*size* effect --
    e.g. small devices carrying proportionally more fixed overhead. Sorts
    the benchmark devices by Qc and reports the parasitic fraction
    (W_parasitic / Qc) for each. Returns the sorted (name, Qc, fraction)
    list so the caller/tests can check the (non-)monotonicity claim."""
    points = points if points is not None else CALIBRATION_POINTS_EXTENDED
    rows = [(name, Qc, Wp / Qc) for (name, f, H, mdot, Qc, Wp) in points]
    rows.sort(key=lambda r: r[1])
    if verbose:
        print("Parasitic fraction (W_parasitic / Qc) sorted by device scale:")
        for name, Qc, frac in rows:
            print(f"    {name:<28} Qc={Qc:8.1f} W   fraction={frac:.3f}")
        fracs = [r[2] for r in rows]
        monotonic = all(fracs[i] <= fracs[i + 1] for i in range(len(fracs) - 1)) or \
                    all(fracs[i] >= fracs[i + 1] for i in range(len(fracs) - 1))
        print(f"  Monotonic in device scale? {monotonic}")
        if not monotonic:
            print(f"  CONCLUSION: no monotonic size trend in this {len(rows)}-point "
                  "set -- the smallest device (Tusek, 6.5W) does NOT have the "
                  "highest overhead fraction, and the largest (Astronautics, "
                  "2502W) does NOT have the lowest. A simple size/scale term is "
                  "not supported by the data in hand; the Astronautics outlier is "
                  "independently attributed by its own source paper to 'mediocre' "
                  "electrical-component efficiency at that scale, i.e. a "
                  "device-specific engineering choice, not a generic size law. "
                  "With the Lozano points included, the picture sharpens further: "
                  "the four highest fractions in the whole set (1.18-1.68) are "
                  "ALL Lozano points clustered in the low-to-mid Qc range, sitting "
                  "well above Okamura (200W, 0.367) and Astronautics (2502W, "
                  "0.453) -- i.e. grouped by device/paper, not ordered by scale.")
    return rows


def run_extended_diagnostic():
    """Demonstrates the instability of the four-point extended fit.
    The diagnostic is provided for transparency and is not used as the
    production calibration."""
    print("=" * 90)
    print("DIAGNOSTIC: adding Okamura & Hirano (2013) as a fourth calibration point")
    print("=" * 90)
    calibrate_loss_coefficients(CALIBRATION_POINTS_EXTENDED, verbose=True,
                                  label="EXTENDED (4pt, diagnostic only, NNLS)")
    print("\n  Leave-one-out cross-validation on the EXTENDED set (NNLS per fold):")
    loo = leave_one_out_cv(CALIBRATION_POINTS_EXTENDED, verbose=True)
    worst = max(loo, key=lambda r: abs(r[3]))
    print(f"\n  CONCLUSION: switching from unconstrained lstsq to NNLS removes the "
          f"negative (unphysical) coefficients Phase 6 found, and improves the "
          f"worst leave-one-out error from +1639% to {worst[3]:+.0f}% "
          f"(held-out device: {worst[0]}) -- but that is still an order-of-"
          f"magnitude miss. A better-behaved solver alone does not make a "
          f"single linear model generalize across devices spanning 6.5W to "
          f"2502W of cooling capacity. The CORE 3-point fit remains the "
          f"production default.")
    print("\n  Testing the natural next hypothesis -- that this is a simple "
          "device-size effect:")
    analyze_parasitic_fraction_scaling(CALIBRATION_POINTS_EXTENDED, verbose=True)


def run_further_extended_diagnostic():
    """Adds the 4 calibratable Lozano et al. (2016) POLO/UFSC points to the
    EXTENDED set (8 points total spanning 6.5-2502W) and re-checks NNLS
    fit stability and the scale-monotonicity question. This device class
    is qualitatively different from CORE/EXTENDED: Lozano's own paper
    reports COP of 0.37-0.83 (vs. 1.9-4.6 for every other benchmark
    device) because its motor power (87-145W) is comparable to or exceeds
    its own Qc (61-120W) -- an early-generation, unoptimized rotary
    prototype the paper itself calls "modest ... in comparison with
    established cooling technologies". The implied parasitic fractions
    here (1.18-1.68) exceed 1.0, i.e. parasitic power alone exceeds Qc --
    outside the range of any CORE/EXTENDED device (max 0.453)."""
    print("=" * 90)
    print("DIAGNOSTIC: adding 4 calibratable Lozano POLO/UFSC (2016) points "
          "(8pt FURTHER_EXTENDED)")
    print("=" * 90)
    calibrate_loss_coefficients(CALIBRATION_POINTS_FURTHER_EXTENDED, verbose=True,
                                  label="FURTHER_EXTENDED (8pt, diagnostic only, NNLS)")
    print("\n  Leave-one-out cross-validation on the FURTHER_EXTENDED set (NNLS per fold):")
    loo = leave_one_out_cv(CALIBRATION_POINTS_FURTHER_EXTENDED, verbose=True)
    worst = max(loo, key=lambda r: abs(r[3]))
    lozano_loo = [r for r in loo if r[0].startswith("Lozano")]
    lozano_mean_err = np.mean([r[3] for r in lozano_loo])
    print(f"\n  CONCLUSION: the naive 'worst leave-one-out error' metric actually "
          f"IMPROVES with 8 points ({worst[3]:+.0f}%, held-out device: {worst[0]}) "
          f"versus the 4-point EXTENDED set's +682% -- but this is misleading, "
          f"not genuine progress: it happens because Tusek's small W_parasitic "
          f"becomes easier to hit once more low/mid-Qc points anchor the fit, "
          f"not because the model now generalizes to the Lozano device class. "
          f"The more informative number is that all 4 Lozano points, when held "
          f"out, are underpredicted by a similar amount (mean "
          f"{lozano_mean_err:+.0f}%, individually {[f'{r[3]:+.0f}%' for r in lozano_loo]}) "
          f"-- a CONSISTENT, one-directional miss, not noise. That consistency is "
          f"itself evidence for the standing hypothesis: Lozano's device sits in a "
          f"different, definably-worse motor/inverter efficiency class (its own "
          f"paper: COP 'modest ... in comparison with established cooling "
          f"technologies', motor power comparable to or exceeding Qc) that a "
          f"state-variable-only loss model (frequency, field, mdot, Qc) cannot "
          f"represent, regardless of solver. The CORE 3-point fit remains the "
          f"production default.")
    print("\n  Re-testing scale monotonicity with the enlarged 8-device set:")
    analyze_parasitic_fraction_scaling(CALIBRATION_POINTS_FURTHER_EXTENDED, verbose=True)


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
    print()
    run_further_extended_diagnostic()