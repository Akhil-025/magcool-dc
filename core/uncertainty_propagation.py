"""
uncertainty_propagation.py
===========================
Phase 30 addition.

`core/loss_model.py`'s CORE calibration is an exactly-determined 3-point
NNLS fit (`leave_one_out_cv()` already documents leave-one-out errors up to
+1639% for the unstable EXTENDED 4-point fit, and the CORE 3-point fit
itself has no closed-form parameter covariance since it is exactly
determined). `results/comparison_table.csv` and `results/pareto_front.csv`
currently report single point-value COP/Qc predictions with no uncertainty
band at all, which is a real gap for a paper making quantitative COP
claims.

This module propagates that calibration uncertainty forward via Monte
Carlo (Metropolis-free, direct resampling), NOT by inventing a parameter
covariance matrix the exactly-determined 3-point fit does not have.
Instead, it perturbs the three calibration *data points themselves*
(the measured Wp/Qc/f/H/mdot benchmark values `CALIBRATION_POINTS_CORE`
is built from) within a literature-reasonable measurement-uncertainty
band, re-fits the NNLS loss coefficients for each perturbed draw, and runs
`AMRSystem.run()` at a fixed representative operating point for each
resulting loss model -- giving an empirical distribution over
COP_electrical and Qc that reflects genuine calibration-input
uncertainty, not just numerical fit noise.

This is a pragmatic choice, not a claim of rigor beyond what the input
data supports: with only 3 calibration points and no reported
device-level measurement uncertainties in the source papers, an assumed
+/-X% relative perturbation is itself an assumption, stated explicitly
via `measurement_uncertainty_frac` (default 15%, a round, conservative
number reflecting typical prototype-scale instrumentation uncertainty in
this class of experiment; NOT sourced from a specific paper's stated
error bars, since none of the three CORE benchmark papers report them).
"""

from dataclasses import dataclass
from typing import List

import numpy as np

from core.loss_model import CALIBRATION_POINTS_CORE, _build_system, StateDependentLossModel
from core.amr_cycle import AMRSystem
from core.mce_material import GADOLINIUM
from scipy.optimize import nnls


@dataclass
class MCResult:
    n_draws: int
    measurement_uncertainty_frac: float
    COP_electrical_mean: float
    COP_electrical_std: float
    COP_electrical_p05: float
    COP_electrical_p95: float
    Qc_mean: float
    Qc_std: float
    Qc_p05: float
    Qc_p95: float
    fraction_failed_draws: float


def _perturbed_points(points, frac, rng):
    """Return a copy of CALIBRATION_POINTS_CORE-style tuples with Wp and Qc
    (the two measured quantities the NNLS fit is regressed against/from)
    perturbed by independent Gaussian noise of std = frac * value. f, H,
    mdot are treated as accurately known (frequency/field/flow are set
    points in these experiments, not measured outputs) -- only the derived
    thermal/electrical quantities Qc and Wp carry the assumed measurement
    uncertainty."""
    out = []
    for (name, f, H, mdot, Qc, Wp) in points:
        Qc_p = max(1e-6, rng.normal(Qc, frac * abs(Qc)))
        Wp_p = max(1e-6, rng.normal(Wp, frac * abs(Wp)))
        out.append((name, f, H, mdot, Qc_p, Wp_p))
    return out


def monte_carlo_cop_uncertainty(
    n_draws=2000,
    measurement_uncertainty_frac=0.15,
    T_cold_K=291.15, span_K=10.0,
    mu0H_max=2.0, mass_regenerator=5.0, frequency=2.0, fluid_mdot=0.08,
    regenerator_effectiveness=0.85,
    seed=0, verbose=True,
) -> MCResult:
    """Monte Carlo propagation of CORE calibration-input uncertainty
    through to COP_electrical/Qc at a fixed representative operating point
    (defaults match main.py's own 10K-span baseline row)."""
    rng = np.random.default_rng(seed)
    cops, qcs, n_failed = [], [], 0
    for _ in range(n_draws):
        pts = _perturbed_points(CALIBRATION_POINTS_CORE, measurement_uncertainty_frac, rng)
        A, b = _build_system(pts)
        coeffs, _ = nnls(A, b)
        k_eddy, k_pump, base_frac = coeffs
        loss_model = StateDependentLossModel(k_eddy=k_eddy, k_pump=k_pump, base_frac=base_frac)
        sys_ = AMRSystem(GADOLINIUM, mu0H_max=mu0H_max, mass_regenerator=mass_regenerator,
                          frequency=frequency, fluid_cp=4186.0, fluid_mdot=fluid_mdot,
                          regenerator_effectiveness=regenerator_effectiveness,
                          loss_model=loss_model, use_ntu_thermal_model=True)
        try:
            res = sys_.run(T_cold_K, span_K)
            if res.COP_electrical > 0 and np.isfinite(res.COP_electrical):
                cops.append(res.COP_electrical)
                qcs.append(res.Qc)
            else:
                n_failed += 1
        except Exception:
            n_failed += 1

    cops = np.array(cops)
    qcs = np.array(qcs)
    if len(cops) == 0:
        if verbose:
            print(f"  WARNING: all {n_draws} draws infeasible at span={span_K}K "
                  f"(0-D model structural span cap, see regenerator_1d.py) -- "
                  f"no COP/Qc distribution computable at this point.")
        return MCResult(n_draws=n_draws, measurement_uncertainty_frac=measurement_uncertainty_frac,
                         COP_electrical_mean=float("nan"), COP_electrical_std=float("nan"),
                         COP_electrical_p05=float("nan"), COP_electrical_p95=float("nan"),
                         Qc_mean=float("nan"), Qc_std=float("nan"),
                         Qc_p05=float("nan"), Qc_p95=float("nan"),
                         fraction_failed_draws=1.0)
    result = MCResult(
        n_draws=n_draws,
        measurement_uncertainty_frac=measurement_uncertainty_frac,
        COP_electrical_mean=float(np.mean(cops)),
        COP_electrical_std=float(np.std(cops)),
        COP_electrical_p05=float(np.percentile(cops, 5)),
        COP_electrical_p95=float(np.percentile(cops, 95)),
        Qc_mean=float(np.mean(qcs)),
        Qc_std=float(np.std(qcs)),
        Qc_p05=float(np.percentile(qcs, 5)),
        Qc_p95=float(np.percentile(qcs, 95)),
        fraction_failed_draws=n_failed / n_draws,
    )
    if verbose:
        print(f"Monte Carlo COP/Qc uncertainty, {n_draws} draws, "
              f"+/-{measurement_uncertainty_frac*100:.0f}% assumed calibration-input "
              f"noise, at span={span_K}K, mu0H={mu0H_max}T:")
        print(f"  COP_electrical: mean={result.COP_electrical_mean:.3f}  "
              f"std={result.COP_electrical_std:.3f}  "
              f"90% CI=[{result.COP_electrical_p05:.3f}, {result.COP_electrical_p95:.3f}]")
        print(f"  Qc (W):         mean={result.Qc_mean:.1f}  std={result.Qc_std:.1f}  "
              f"90% CI=[{result.Qc_p05:.1f}, {result.Qc_p95:.1f}]")
        print(f"  fraction of draws with infeasible/failed COP: "
              f"{result.fraction_failed_draws*100:.1f}%")
        print("  HONEST FRAMING FOR THE PAPER: the +/-15% default is a "
              "reasoned assumption, not a source-derived measurement "
              "uncertainty (none of the 3 CORE calibration papers report "
              "device-level error bars). Report this CI as 'uncertainty "
              "under an assumed 15% calibration-input measurement noise', "
              "not as a rigorously derived confidence interval, and run a "
              "sensitivity check at e.g. 10%/25% to show how the band scales.")
    return result


def uncertainty_band_across_spans(spans_K=range(5, 21), n_draws=500,
                                   measurement_uncertainty_frac=0.15,
                                   verbose=True) -> List[dict]:
    """Same Monte Carlo propagation repeated across the full 5-20K ASHRAE
    span sweep, for a confidence-band overlay on comparison_table.csv's
    AMR_COP_electrical curve (Tier 1 item 5 in the original review)."""
    rows = []
    for span in spans_K:
        mc = monte_carlo_cop_uncertainty(n_draws=n_draws,
                                          measurement_uncertainty_frac=measurement_uncertainty_frac,
                                          span_K=float(span), verbose=False)
        rows.append({
            "span_K": span,
            "COP_electrical_mean": round(mc.COP_electrical_mean, 3),
            "COP_electrical_p05": round(mc.COP_electrical_p05, 3),
            "COP_electrical_p95": round(mc.COP_electrical_p95, 3),
            "Qc_mean_W": round(mc.Qc_mean, 1),
            "Qc_p05_W": round(mc.Qc_p05, 1),
            "Qc_p95_W": round(mc.Qc_p95, 1),
        })
        if verbose:
            print(f"  span={span:>2}K  COP_electrical mean={rows[-1]['COP_electrical_mean']:.2f}  "
                  f"90% CI=[{rows[-1]['COP_electrical_p05']:.2f}, "
                  f"{rows[-1]['COP_electrical_p95']:.2f}]")
    return rows


def write_uncertainty_report(path="results/uncertainty_propagation.txt",
                              csv_path="results/comparison_table_with_uncertainty.csv"):
    import io, contextlib, csv as csv_mod
    buf = io.StringIO()
    with contextlib.redirect_stdout(buf):
        monte_carlo_cop_uncertainty(verbose=True)
        print()
        print("Full 5-20K span sweep, 90% CI band:")
        rows = uncertainty_band_across_spans(verbose=True)
    with open(path, "w") as f:
        f.write(buf.getvalue())
    with open(csv_path, "w", newline="") as f:
        w = csv_mod.DictWriter(f, fieldnames=rows[0].keys())
        w.writeheader()
        w.writerows(rows)
    print(buf.getvalue())
    print(f"Wrote {path}")
    print(f"Wrote {csv_path}")


if __name__ == "__main__":
    write_uncertainty_report()
