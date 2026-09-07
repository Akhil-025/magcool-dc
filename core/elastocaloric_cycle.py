"""
elastocaloric_cycle.py
=======================
0-D device-level electrical-COP model for an elastocaloric (NiTi,
stress-driven) cooling cycle, spanning 5-20 K to match this repo's own
comparison_table.csv span range for VCC/AMR.

Structurally mirrors core/amr_cycle.py's own COP_electrical pattern
(Qc / (W_ideal + W_parasitic)), NOT a copy of its magnetic-specific
mechanics:

    1. Required stress: the stress sigma_needed such that the material's
       own Delta_T_ad(sigma_needed) exceeds the target span by a headroom
       margin (span_margin, default 1.3x) -- the SAME structural idea
       amr_cycle.py uses (a regenerator needs a no-load Delta_T_ad bigger
       than the target span; see this repo's own
       regenerative_amplification_override_check.txt for exactly this
       concept in the magnetic case). This repo does NOT have a 1-D
       elastocaloric regenerator model (the magnetic equivalent of
       regenerator_1d.py), so the margin is a documented, adjustable
       placeholder, not a fitted or literature-derived number.
    2. Ideal cooling capacity per unit mass per cycle:
           Qc = c_p * (Delta_T_ad(sigma_needed) - span)
       i.e. the material's adiabatic swing minus what regeneration must
       "spend" reaching the target span -- again, the same structural
       relationship amr_cycle.py uses for cooling_capacity(), not a new
       invented relation.
    3. Work input per unit mass per cycle:
           W = hysteresis_loss_J_per_kg(sigma_needed) + actuator parasitic
       The actuator parasitic term is calibrated (see calibrate_parasitic_
       fraction_to_literature() below) against MEASURED system-level COPs
       from the literature, exactly the way amr_cycle.py's own
       parasitic_fraction was calibrated against measured benchmark
       devices (README.md "Key finding: ideal vs. electrical COP").

Calibration targets (measured, not simulated, system-level COP)
------------------------------------------------------------------
  - Zhang et al., "Continuous and efficient elastocaloric air cooling by
    coil-bending", Nature Communications 14 (2023) 8323: measured system
    COP = 3.7 (cooling power / rotational mechanical power), 10.6 K span.
  - "Real-time AI-optimized elastocaloric cooling... compression-mode
    Ni-Ti systems", ScienceDirect (2025): measured COP = 2.8-3.1 at a span
    implied by the paper's own stated operating point (this repo uses the
    midpoint COP=2.95 at an assumed representative ~15 K span, since the
    exact span was not confirmed from the available excerpt -- FLAGGED as
    an approximation, not a directly-read number).
  - These are BOTH well below the 6.85 "projected, not yet demonstrated"
    figure and the qualitative "up to 22" theoretical-ceiling figure that
    earlier literature searches surfaced -- this module deliberately
    calibrates against the lower, MEASURED numbers, following this
    repo's own stated preference (see README.md's Sobol/COP sections)
    for electrical/measured COP over ideal/projected COP.
"""

import numpy as np
from dataclasses import dataclass

from core.elastocaloric_material import NiTi_binary


@dataclass
class ElastocaloricCycleResult:
    span_K: float
    sigma_needed_MPa: float
    delta_T_ad_at_sigma: float
    Qc_per_kg: float
    W_hysteresis_per_kg: float
    W_parasitic_per_kg: float
    COP_ideal: float          # Qc / W_hysteresis alone (no actuator parasitic)
    COP_electrical: float     # Qc / (W_hysteresis + W_parasitic)
    feasible: bool            # False if material can't reach this span at all


MEASURED_BENCHMARKS = [
    # (span_K, measured_system_COP, source)
    (10.6, 3.7, "Zhang et al., Nat. Commun. 14 (2023) 8323 (coil-bending, measured)"),
    (15.0, 2.95, "'Real-time AI-optimized elastocaloric cooling', ScienceDirect (2025), "
                 "measured COP=2.8-3.1 midpoint, span approximated -- see docstring"),
]


def _sigma_needed_for_span(material, span_K, span_margin=1.3):
    """Smallest stress such that Delta_T_ad(sigma) >= span_margin * span_K.
    Closed-form: transformed_fraction() is linear in sigma up to
    sigma_sat_MPa (see elastocaloric_material.py docstring), so
    Delta_T_ad(sigma) = (T*Delta_S/c_p) * (sigma/sigma_sat) for
    sigma <= sigma_sat, and constant beyond. Solving directly (rather than
    scanning) avoids an expensive nested search inside the parasitic-
    fraction calibration's own search loop below."""
    target = span_margin * span_K
    dT_max = material.delta_T_ad(material.sigma_sat_MPa)
    if target > dT_max:
        return material.sigma_sat_MPa, False
    if target <= 0:
        return 0.0, True
    sigma_needed = target / dT_max * material.sigma_sat_MPa
    return float(sigma_needed), True


def run_cycle(span_K, T_cold=291.15, material=NiTi_binary, span_margin=1.3,
              parasitic_fraction_of_Qc=0.0):
    """One span-point evaluation. parasitic_fraction_of_Qc mirrors
    amr_cycle.py's own W_parasitic = parasitic_fraction * Qc model exactly."""
    T = T_cold + span_K / 2.0
    sigma_needed, feasible = _sigma_needed_for_span(material, span_K, span_margin)
    dT_ad = material.delta_T_ad(sigma_needed, T=T)

    if not feasible:
        return ElastocaloricCycleResult(span_K, sigma_needed, dT_ad, 0.0, 0.0, 0.0,
                                         0.0, 0.0, feasible=False)

    Qc = max(material.c_p_J_per_kgK * (dT_ad - span_K), 0.0)
    W_hys = material.hysteresis_loss_J_per_kg(sigma_needed)
    W_par = parasitic_fraction_of_Qc * Qc

    COP_ideal = Qc / W_hys if W_hys > 0 else 0.0
    COP_elec = Qc / (W_hys + W_par) if (W_hys + W_par) > 0 else 0.0

    return ElastocaloricCycleResult(span_K, sigma_needed, dT_ad, Qc, W_hys, W_par,
                                     COP_ideal, COP_elec, feasible=True)


def calibrate_parasitic_fraction_to_literature(material=NiTi_binary, span_margin=1.3,
                                                 verbose=True):
    """Back-solves parasitic_fraction_of_Qc (a single, constant, repo-wide
    value, matching amr_cycle.py's own constant-parasitic_fraction
    convention) so that this model's COP_electrical matches the MEASURED
    benchmarks in MEASURED_BENCHMARKS as closely as possible in a
    least-squares sense. This is the elastocaloric analogue of this
    repo's own "Key finding: ideal vs. electrical COP" calibration for
    magnetocaloric materials."""
    def cop_electrical_given_fraction(span, frac):
        r = run_cycle(span, material=material, span_margin=span_margin,
                       parasitic_fraction_of_Qc=frac)
        return r.COP_electrical

    def total_sq_err(frac):
        err = 0.0
        for span, cop_lit, _src in MEASURED_BENCHMARKS:
            cop_model = cop_electrical_given_fraction(span, frac)
            err += (cop_model - cop_lit) ** 2
        return err

    # UPGRADE (this pass): the previous 500-point grid over [0, 5]
    # (spacing ~0.01) was a source of avoidable calibration error --
    # the barocaloric module's single-point version of this same pattern
    # was caught by its own test failing outside a 0.05 COP tolerance
    # purely from grid coarseness (see barocaloric_cycle.py's fix).
    # With several benchmark points there's no closed form here, so use
    # scipy's bounded scalar minimizer instead of a fixed grid -- same
    # least-squares objective, much finer precision, no extra assumptions.
    from scipy.optimize import minimize_scalar
    res = minimize_scalar(total_sq_err, bounds=(0.0, 5.0), method="bounded",
                           options={"xatol": 1e-8})
    best_frac, best_err = float(res.x), float(res.fun)

    if verbose:
        print("Elastocaloric parasitic_fraction_of_Qc calibration "
              "(least-squares against MEASURED literature benchmarks):")
        for span, cop_lit, src in MEASURED_BENCHMARKS:
            cop_model = cop_electrical_given_fraction(span, best_frac)
            print(f"  span={span:.1f}K  lit_COP={cop_lit:.2f}  "
                  f"model_COP={cop_model:.2f}  ({src})")
        print(f"  -> calibrated parasitic_fraction_of_Qc = {best_frac:.4f} "
              f"(RMS error = {np.sqrt(best_err / len(MEASURED_BENCHMARKS)):.3f})")

    return best_frac


def run_span_sweep(spans_K=(5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20),
                    material=NiTi_binary, span_margin=1.3, parasitic_fraction=None,
                    verbose=True):
    if parasitic_fraction is None:
        parasitic_fraction = calibrate_parasitic_fraction_to_literature(
            material=material, span_margin=span_margin, verbose=verbose)

    rows = []
    for span in spans_K:
        r = run_cycle(span, material=material, span_margin=span_margin,
                       parasitic_fraction_of_Qc=parasitic_fraction)
        rows.append(r)
    return rows, parasitic_fraction


if __name__ == "__main__":
    rows, frac = run_span_sweep()
    print("\nElastocaloric electrical COP vs. span (calibrated model):")
    for r in rows:
        print(f"  span={r.span_K:.1f}K  feasible={r.feasible}  "
              f"sigma_needed={r.sigma_needed_MPa:.0f}MPa  "
              f"COP_ideal={r.COP_ideal:.2f}  COP_electrical={r.COP_electrical:.2f}")