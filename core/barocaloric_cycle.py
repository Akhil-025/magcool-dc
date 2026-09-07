"""
barocaloric_cycle.py
=====================
0-D device-level electrical-COP model for a barocaloric (NPG plastic
crystal, pressure-driven) cooling cycle, structured identically to
elastocaloric_cycle.py -- see that module's docstring for the shared
Qc/(W_hysteresis + W_parasitic) structure, which itself mirrors
core/amr_cycle.py's own COP_electrical pattern.

CALIBRATION HONESTY FLAG (read before trusting any number here): unlike
elastocaloric_cycle.py, this project could NOT locate an independently
MEASURED end-to-end barocaloric device COP anywhere in the sources
checked. The only NPG-specific device-level COP figure found is a
SIMULATION: "a BC refrigerator operating with a temperature span of 2.4 K
and 0.1 GPa applied pressure can achieve ... a COP as high as 5.5 at
1 mHz cycle frequency" ('Modeling of an Elastocaloric Cooling System for
Determining Efficiency', ResearchGate/arXiv preprint). Two things about
that figure matter before using it as a calibration target:
  1. It is a MODELED figure from a different research group's model, not
     a measurement -- calibrating this repo's own model against another
     group's model, rather than against hardware, is a materially weaker
     form of calibration than elastocaloric_cycle.py's benchmarks. This
     module reports COP against that single simulated point and states
     the weaker confidence explicitly wherever the number is used
     downstream (see alternative_caloric_comparison.py).
  2. 1 mHz is an extremely slow cycle frequency (one full cycle every
     ~17 minutes) -- the cited COP is achieved by trading almost all
     cooling POWER away (a real device running that slowly would deliver
     negligible Qc in watts, even though the per-cycle COP ratio looks
     good). This repo does NOT model frequency-dependent parasitic losses
     for barocaloric cycles (no analogue of amr_cycle.py's W_eddy~f^2 term
     exists here, since no comparable loss-scaling literature for
     pressure-driven cycles was located) -- so this module's COP numbers
     should be read as "per-cycle, ignoring rate-dependent losses",
     NOT as "achievable at data-center-relevant cooling power", which is
     a materially different and unverified claim.
"""

import numpy as np
from dataclasses import dataclass

from core.barocaloric_material import NPG_plastic_crystal


@dataclass
class BarocaloricCycleResult:
    span_K: float
    P_needed_MPa: float
    delta_T_ad_at_P: float
    Qc_per_kg: float
    W_hysteresis_per_kg: float
    W_parasitic_per_kg: float
    COP_ideal: float
    COP_electrical: float
    feasible: bool


# See module docstring: SIMULATED, not measured; single source; weak
# calibration target, used because no measured alternative exists.
SIMULATED_CALIBRATION_TARGET = (
    2.4,   # span_K
    5.5,   # COP at 1 mHz (rate-dependent losses not modeled here -- see docstring)
    "'Modeling of an Elastocaloric Cooling System for Determining "
    "Efficiency' (ResearchGate, PDF), BC refrigerator sub-case, "
    "SIMULATED not measured, 1 mHz -- extremely low cycle rate, see docstring",
)


def _pressure_needed_for_span(material, span_K, span_margin=1.3):
    """Closed-form inversion -- see elastocaloric_cycle.py's
    _sigma_needed_for_span() for why (same linear-clamp material model,
    same rationale for avoiding a scan inside the calibration loop)."""
    target = span_margin * span_K
    dT_max = material.delta_T_ad(material.P_sat_MPa)
    if target > dT_max:
        return material.P_sat_MPa, False
    if target <= 0:
        return 0.0, True
    P_needed = target / dT_max * material.P_sat_MPa
    return float(P_needed), True


def run_cycle(span_K, T_cold=291.15, material=NPG_plastic_crystal, span_margin=1.3,
              parasitic_fraction_of_Qc=0.0):
    T = T_cold + span_K / 2.0
    P_needed, feasible = _pressure_needed_for_span(material, span_K, span_margin)
    dT_ad = material.delta_T_ad(P_needed, T=T)

    if not feasible:
        return BarocaloricCycleResult(span_K, P_needed, dT_ad, 0.0, 0.0, 0.0,
                                       0.0, 0.0, feasible=False)

    Qc = max(material.c_p_J_per_kgK * (dT_ad - span_K), 0.0)
    W_hys = material.hysteresis_loss_J_per_kg(P_needed)
    W_par = parasitic_fraction_of_Qc * Qc

    COP_ideal = Qc / W_hys if W_hys > 0 else 0.0
    COP_elec = Qc / (W_hys + W_par) if (W_hys + W_par) > 0 else 0.0

    return BarocaloricCycleResult(span_K, P_needed, dT_ad, Qc, W_hys, W_par,
                                   COP_ideal, COP_elec, feasible=True)


def calibrate_parasitic_fraction_to_simulation(material=NPG_plastic_crystal,
                                                 span_margin=1.3, verbose=True):
    """Solves for parasitic_fraction_of_Qc against the single SIMULATED
    target above -- see module docstring for why this is a materially
    weaker calibration than the elastocaloric case.

    UPGRADE (this pass): with only one calibration point, COP_electrical =
    Qc / (W_hys + frac*Qc) is monotonically decreasing in frac at fixed
    span, so the exact solution is closed-form:
        frac = 1/COP_target - W_hys/Qc
    The previous implementation instead did a coarse 500-point grid search
    over frac in [0, 5] (spacing ~0.01), which is a leftover from
    elastocaloric_cycle.py's multi-point least-squares pattern (see that
    module's calibrate_parasitic_fraction_to_literature(), where a closed
    form doesn't exist because there are several benchmark spans at once).
    For a single point that grid was needlessly imprecise -- it landed
    >0.05 COP off target (test_calibration_matches_simulated_target_closely
    caught this) purely from grid spacing, not from any physical
    limitation. Solving in closed form removes that error entirely."""
    span, cop_target, src = SIMULATED_CALIBRATION_TARGET

    def cop_given_fraction(frac):
        r = run_cycle(span, material=material, span_margin=span_margin,
                       parasitic_fraction_of_Qc=frac)
        return r.COP_electrical

    r0 = run_cycle(span, material=material, span_margin=span_margin,
                    parasitic_fraction_of_Qc=0.0)
    if r0.Qc_per_kg > 0:
        best_frac = max(0.0, 1.0 / cop_target - r0.W_hysteresis_per_kg / r0.Qc_per_kg)
    else:
        best_frac = 0.0

    if verbose:
        cop_model = cop_given_fraction(best_frac)
        print("Barocaloric parasitic_fraction_of_Qc calibration "
              "(against a SIMULATED, not measured, target -- see docstring):")
        print(f"  span={span:.1f}K  target_COP={cop_target:.2f}  "
              f"model_COP={cop_model:.2f}  ({src})")
        print(f"  -> calibrated parasitic_fraction_of_Qc = {best_frac:.4f}")

    return best_frac


def run_span_sweep(spans_K=(5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20),
                    material=NPG_plastic_crystal, span_margin=1.3,
                    parasitic_fraction=None, verbose=True):
    if parasitic_fraction is None:
        parasitic_fraction = calibrate_parasitic_fraction_to_simulation(
            material=material, span_margin=span_margin, verbose=verbose)

    rows = []
    for span in spans_K:
        r = run_cycle(span, material=material, span_margin=span_margin,
                       parasitic_fraction_of_Qc=parasitic_fraction)
        rows.append(r)
    return rows, parasitic_fraction


if __name__ == "__main__":
    rows, frac = run_span_sweep()
    print("\nBarocaloric electrical COP vs. span (calibrated against a "
          "SIMULATED target -- treat with the weaker confidence the "
          "docstring describes):")
    for r in rows:
        print(f"  span={r.span_K:.1f}K  feasible={r.feasible}  "
              f"P_needed={r.P_needed_MPa:.1f}MPa  "
              f"COP_ideal={r.COP_ideal:.2f}  COP_electrical={r.COP_electrical:.2f}")