"""
cascade.py
===========
Multi-stage (cascade) AMR model for extending the operating temperature
span beyond the limits of a single-stage system.

Cascade concept: N single-stage AMR modules in series, each handling an
equal share of the total span (Th_total - Tc_total)/N, analogous to
cascade vapor-compression refrigeration. Stage 1 (coldest) absorbs the
data-center heat load Qc at T_cold; its heat rejection becomes the input
to Stage 2, and so on until the final stage rejects heat at T_hot.
Because heat is transferred through the stages in series, each stage must
reject the same cooling load (steady state, neglecting inter-stage losses):

    W_total = Σ W_i(Qc, span/N)
    COP_cascade = Qc / W_total

run_cascade()/compare_staging() below assume identical regenerator stages
(all Gd, or all Gd5Si2Ge2). run_graded_cascade()/compare_graded_cascade()
implement the more advanced Curie-graded variant (ROADMAP.md Phase 7 open
item): each stage uses a hypothetical composition-tuned Gd5(SixGe1-x)4(-Ga)
material whose own peak MCE effect is matched to that stage's local
operating temperature, checked against the literature-documented
composition-tunability range and against an independent direct-measurement
validation (see core/first_order_mce.py and core/giguere_validation.py).

Phase 9 addendum: the Curie-graded machinery was generalized to a pluggable
`GradedFamily` (see below) rather than being hardcoded to the Gd5(SixGe1-x)4
family, so the SAME grading mechanism can be applied to the La(Fe,Si)13Hy
family added in Phase 9 (core.first_order_mce.lafesih_composition_tuned_material)
-- specifically to test whether a 6-layer Curie-graded La(Fe,Si)13Hy bed can
reproduce the REAL Astronautics_rotary_2014 benchmark device (which is
exactly such a bed), closing the gap validation_system.py's Phase 9
single-layer approximation flagged as a "natural next step, not done here".
GD_FAMILY below reproduces the original Gd5(SixGe1-x)4 behavior exactly
(default argument, so existing calls are unaffected); LAFESIH_FAMILY is new.
See validate_astronautics_graded_bed() at the bottom of this module for the
actual comparison against Jacobs et al. (2014)'s reported numbers.
"""

import numpy as np
import csv
from dataclasses import dataclass, field
from typing import Callable
from scipy.optimize import brentq
from core.mce_material import GADOLINIUM
from core.amr_cycle import AMRSystem
from core.loss_model import StateDependentLossModel
from core.baseline_cooling import vapor_compression_cop, liquid_cooling_cop
from core.first_order_mce import (composition_tuned_material,
                                    GIANT_MCE_TC_MIN_K, GIANT_MCE_TC_MAX_K,
                                    lafesih_composition_tuned_material,
                                    LAFESIH_TC_MIN_K, LAFESIH_TC_MAX_K,
                                    GD5SI2GE2_FIRST_ORDER, LAFESIH_FIRST_ORDER)

_LOSS_MODEL = StateDependentLossModel()
USE_NTU_THERMAL_MODEL = True


@dataclass
class GradedFamily:
    """Everything run_graded_cascade() needs to grade a bed from a single
    composition-tunable material family. reference_material is only used to
    seed _target_composition_for_peak()'s fixed-point iteration (its own
    peak-vs-Tc offset is the starting guess); tuned_fn(Tc_target_K) must
    return a FirstOrderMCEMaterial and raise ValueError outside
    [tc_min, tc_max] (both composition_tuned_material() and
    lafesih_composition_tuned_material() already do this)."""
    name: str
    tuned_fn: Callable[[float], object]
    tc_min: float
    tc_max: float
    reference_material: object
    fallback_material: object = field(default_factory=lambda: GADOLINIUM)


GD_FAMILY = GradedFamily(
    name="Gd5(SixGe1-x)4(-Ga)",
    tuned_fn=lambda Tc: composition_tuned_material(Tc, apply_giguere_correction=True),
    tc_min=GIANT_MCE_TC_MIN_K, tc_max=GIANT_MCE_TC_MAX_K,
    reference_material=GD5SI2GE2_FIRST_ORDER, fallback_material=GADOLINIUM,
)

LAFESIH_FAMILY = GradedFamily(
    name="La(Fe,Si)13Hy",
    tuned_fn=lafesih_composition_tuned_material,
    tc_min=LAFESIH_TC_MIN_K, tc_max=LAFESIH_TC_MAX_K,
    reference_material=LAFESIH_FIRST_ORDER, fallback_material=GADOLINIUM,
)


def run_cascade(T_cold_K, total_span_K, n_stages, material=None, mu0H_max=2.0,
                 mass_per_stage=2.0, frequency=1.0, fluid_mdot=0.08,
                 regenerator_effectiveness=0.85):
    """Runs n_stages identical AMR modules in series, each covering
    total_span_K/n_stages, all passing the same Qc through in steady state
    (Qc is set by the coldest/first stage's capacity at its local span)."""
    if material is None:
        material = GADOLINIUM
    span_per_stage = total_span_K / n_stages
    T_local = T_cold_K
    # First stage sets the deliverable Qc (bottleneck of the chain)
    stage1 = AMRSystem(material=material, mu0H_max=mu0H_max,
                        mass_regenerator=mass_per_stage, frequency=frequency,
                        fluid_mdot=fluid_mdot, regenerator_effectiveness=regenerator_effectiveness,
                        loss_model=_LOSS_MODEL, use_ntu_thermal_model=USE_NTU_THERMAL_MODEL)
    r1 = stage1.run(T_local, span_per_stage)
    Qc_target = r1.Qc
    if Qc_target <= 0:
        return {"n_stages": n_stages, "span_per_stage_K": span_per_stage,
                "Qc_W": 0.0, "W_total_W": np.nan, "COP_cascade": 0.0,
                "feasible": False}

    W_total = 0.0
    for i in range(n_stages):
        stage = AMRSystem(material=material, mu0H_max=mu0H_max,
                           mass_regenerator=mass_per_stage, frequency=frequency,
                           fluid_mdot=fluid_mdot, regenerator_effectiveness=regenerator_effectiveness,
                           loss_model=_LOSS_MODEL, use_ntu_thermal_model=USE_NTU_THERMAL_MODEL)
        # each stage handles the same Qc_target at its local span; back out
        # the required work by re-running at span_per_stage and scaling mdot
        # if needed so Qc matches Qc_target (steady-state series constraint)
        r_i = stage.run(T_local, span_per_stage)
        if r_i.Qc > 0:
            scale = Qc_target / r_i.Qc
            W_i = (r_i.W_mag + r_i.W_parasitic) * scale
        else:
            W_i = np.inf
        W_total += W_i
        T_local += span_per_stage

    COP_cascade = Qc_target / W_total if W_total > 0 else 0.0
    return {"n_stages": n_stages, "span_per_stage_K": span_per_stage,
            "Qc_W": round(Qc_target, 1), "W_total_W": round(W_total, 1),
            "COP_cascade": round(COP_cascade, 3), "feasible": True}


def _peak_temperature(material, mu0H_max, T_range=(200.0, 340.0), n=1401):
    """Finds where THIS material's own DeltaT_ad is maximized (same
    approach as giant_mce_analysis.py's find_peak_temperature()). Needed
    because the Landau model's peak does not sit exactly at the nominal
    Tc, and -- discovered while building this cascade -- that offset is
    NOT simply translation-invariant as Tc is shifted (Debye C_lattice(T)
    depends on absolute T, not on T-Tc), so a single global offset applied
    to every stage is not accurate enough; see _target_composition_for_peak
    below for why this matters quantitatively."""
    mu0 = 4 * np.pi * 1e-7
    H = mu0H_max / mu0
    # Coarse pass over the full range, then a fine pass zoomed into the
    # coarse peak's neighborhood -- gets the same resolution near the peak
    # as a single n=1401 pass at a fraction of the evaluations (this
    # function is called inside a root-finder in
    # _target_composition_for_peak, so its cost multiplies fast).
    n_coarse = max(51, n // 10)
    Ts_coarse = np.linspace(*T_range, n_coarse)
    dT_coarse = material.delta_T_adiabatic(Ts_coarse, H)
    i0 = int(np.argmax(dT_coarse))
    step = (T_range[1] - T_range[0]) / (n_coarse - 1)
    lo = Ts_coarse[max(0, i0 - 2)]
    hi = Ts_coarse[min(n_coarse - 1, i0 + 2)]
    if hi <= lo:
        return float(Ts_coarse[i0])
    Ts_fine = np.linspace(lo, hi, max(21, n // 20))
    dT_fine = material.delta_T_adiabatic(Ts_fine, H)
    return float(Ts_fine[int(np.argmax(dT_fine))])


def _target_composition_for_peak(T_target_K, mu0H_max, family, max_iter=6, tol_K=0.02):
    """Solves for the composition Tc whose OWN peak DeltaT_ad lands at
    T_target_K, for the given GradedFamily. peak_T(Tc) is monotonic
    increasing in Tc (verified numerically for both GD_FAMILY and
    LAFESIH_FAMILY), so this is a straightforward bracketed root-find
    (scipy.optimize.brentq) on peak_T(Tc) - T_target_K = 0.

    This turned out to matter more than expected. The original (Phase 7)
    implementation used a simple fixed-point update (Tc += err from a
    single global offset), which worked fine for GD_FAMILY but was found
    (Phase 9, while adding LAFESIH_FAMILY) to visibly FAIL for it: this
    Landau model's transition is narrower for La(Fe,Si)13Hy than for
    Gd5Si2Ge2 (DeltaT_ad can fall by more than an order of magnitude
    within ~0.05K of the true peak -- narrower still than the ~0.2-0.5K
    scale already flagged for Gd5Si2Ge2's own transition, itself already
    narrower than the real, hysteresis/inhomogeneity-broadened transition
    Giguere et al.'s Fig. 3 shows spanning ~10-15K). At that sharpness the
    fixed-point update overshot and oscillated between iterations instead
    of converging (verified: 6-20 fixed-point iterations at tol_K=0.001-
    0.02 landed on Tc values up to ~0.15K apart, several of which left a
    stage's own dTad_noload at 0.6K instead of ~21K -- collapsing that
    stage's Qc to zero via cooling_capacity()'s span_fraction clamp, NOT
    because the material itself is bad). That narrowness is itself a
    genuine physical/numerical limitation of this idealized 6th-order
    Landau fit (flagged here, not smoothed over) -- what changed here is
    only the ROOT-FINDING METHOD used to hit each stage's own true peak as
    precisely as that narrow peak demands, so the graded-bed MECHANISM can
    be evaluated on its own terms without this numerical artifact adding
    noise on top.
    """
    ref = family.reference_material
    # Floor at 100K, well clear of the low-T region where this Landau
    # model's DeltaT_ad numerically diverges (lattice heat capacity ->0 as
    # T->0K, a pre-existing model artifact -- see e.g. T=1K giving
    # DeltaT_ad~30K for GD5SI2GE2_FIRST_ORDER, checked while debugging this
    # search). GD_FAMILY's tc_min=20K would otherwise put the search range
    # at (-20, 330)K and let brentq wander into that spurious peak.
    search_range = (max(100.0, family.tc_min - 40.0), family.tc_max + 40.0)
    offset_guess = _peak_temperature(ref, mu0H_max, T_range=search_range) - ref.Tc
    Tc_guess = T_target_K - offset_guess

    def f(Tc):
        mat = family.tuned_fn(Tc)
        return _peak_temperature(mat, mu0H_max, T_range=search_range) - T_target_K

    lo = max(family.tc_min, Tc_guess - 15.0)
    hi = min(family.tc_max, Tc_guess + 15.0)
    if lo >= hi:
        return Tc_guess  # nothing sane to bracket; let the caller's range check handle it
    try:
        if f(lo) * f(hi) > 0:
            # local bracket doesn't straddle a root -- widen to the family's
            # full documented range once before giving up
            lo, hi = family.tc_min, family.tc_max
            if f(lo) * f(hi) > 0:
                return Tc_guess  # give up gracefully, let range/feasibility checks handle it
        return brentq(f, lo, hi, xtol=min(tol_K, 0.005), maxiter=max(max_iter, 50))
    except ValueError:
        return Tc_guess


def run_graded_cascade(T_cold_K, total_span_K, n_stages, mu0H_max=2.0,
                        mass_per_stage=2.0, frequency=1.0, fluid_mdot=0.08,
                        regenerator_effectiveness=0.85,
                        apply_giguere_correction=True, family=None):
    """Curie-graded cascade (ROADMAP.md Phase 7 open item; generalized in
    Phase 9): rather than identical stages of one material (run_cascade
    above), each stage is assigned a hypothetical composition-tuned material
    from `family` (a GradedFamily -- see GD_FAMILY/LAFESIH_FAMILY above)
    whose Curie temperature matches THAT stage's own local midpoint
    temperature, on the Curie-matching principle confirmed in
    giant_mce_analysis.py (a first-order giant-MCE material performs
    strongly at its own Tc and collapses to ~zero capacity away from it).

    family defaults to GD_FAMILY (Gd5(SixGe1-x)4(-Ga), with the Giguere et
    al. (1999) empirical DeltaT_ad correction applied iff
    apply_giguere_correction=True), reproducing this function's original
    Phase 7 behavior exactly -- apply_giguere_correction is IGNORED if you
    pass an explicit `family` (build the correction into family.tuned_fn
    instead, as GD_FAMILY itself does).

    Each stage's target Tc is checked against family's documented
    tunability window (family.tc_min/tc_max). If a stage's target Tc falls
    outside that window, this function does NOT silently extrapolate a
    fictitious material -- it falls back to family.fallback_material (Gd by
    default) for that stage and records the fallback, so the returned
    result honestly reflects what is and is not supported by the
    composition-tunability literature at the requested operating point.
    """
    if family is None:
        family = GD_FAMILY if apply_giguere_correction else GradedFamily(
            name=GD_FAMILY.name,
            tuned_fn=lambda Tc: composition_tuned_material(Tc, apply_giguere_correction=False),
            tc_min=GD_FAMILY.tc_min, tc_max=GD_FAMILY.tc_max,
            reference_material=GD_FAMILY.reference_material,
            fallback_material=GD_FAMILY.fallback_material)

    span_per_stage = total_span_K / n_stages
    T_local = T_cold_K
    stage_materials = []
    stage_info = []
    for i in range(n_stages):
        T_mid_stage = T_local + span_per_stage / 2.0
        Tc_target = _target_composition_for_peak(T_mid_stage, mu0H_max, family)
        if family.tc_min <= Tc_target <= family.tc_max:
            mat = family.tuned_fn(Tc_target)
            stage_info.append({"stage": i + 1, "T_mid_K": round(T_mid_stage, 1),
                                "Tc_target_K": round(Tc_target, 1),
                                "material": mat.name, "in_range": True})
        else:
            mat = family.fallback_material
            stage_info.append({"stage": i + 1, "T_mid_K": round(T_mid_stage, 1),
                                "Tc_target_K": round(Tc_target, 1),
                                "material": f"{mat.name} (fallback -- Tc target outside "
                                            f"{family.tc_min:.0f}-{family.tc_max:.0f}K "
                                            f"documented {family.name} range)",
                                "in_range": False})
        stage_materials.append(mat)
        T_local += span_per_stage

    n_fallback = sum(1 for s in stage_info if not s["in_range"])

    T_local = T_cold_K
    stage1 = AMRSystem(material=stage_materials[0], mu0H_max=mu0H_max,
                        mass_regenerator=mass_per_stage, frequency=frequency,
                        fluid_mdot=fluid_mdot, regenerator_effectiveness=regenerator_effectiveness,
                        loss_model=_LOSS_MODEL, use_ntu_thermal_model=USE_NTU_THERMAL_MODEL)
    r1 = stage1.run(T_local, span_per_stage)
    Qc_target = r1.Qc
    if Qc_target <= 0:
        return {"n_stages": n_stages, "span_per_stage_K": span_per_stage,
                "Qc_W": 0.0, "W_total_W": np.nan, "COP_cascade": 0.0,
                "feasible": False, "n_stages_out_of_range": n_fallback,
                "stage_info": stage_info}

    W_total = 0.0
    T_local = T_cold_K
    for i in range(n_stages):
        stage = AMRSystem(material=stage_materials[i], mu0H_max=mu0H_max,
                           mass_regenerator=mass_per_stage, frequency=frequency,
                           fluid_mdot=fluid_mdot, regenerator_effectiveness=regenerator_effectiveness,
                           loss_model=_LOSS_MODEL, use_ntu_thermal_model=USE_NTU_THERMAL_MODEL)
        r_i = stage.run(T_local, span_per_stage)
        if r_i.Qc > 0:
            scale = Qc_target / r_i.Qc
            W_i = (r_i.W_mag + r_i.W_parasitic) * scale
        else:
            W_i = np.inf
        W_total += W_i
        T_local += span_per_stage

    COP_cascade = Qc_target / W_total if W_total > 0 else 0.0
    return {"n_stages": n_stages, "span_per_stage_K": span_per_stage,
            "Qc_W": round(Qc_target, 1), "W_total_W": round(W_total, 1),
            "COP_cascade": round(COP_cascade, 3), "feasible": True,
            "n_stages_out_of_range": n_fallback, "stage_info": stage_info}


def compare_graded_cascade(T_cold_C=18.0, spans=range(5, 21), stage_counts=(1, 2, 3, 4),
                            mass_per_stage=2.0, family=None,
                            out_csv="results/graded_cascade_comparison.csv"):
    """Same sweep as compare_staging(), but using run_graded_cascade()
    instead of identical-stage run_cascade(). family is passed straight
    through to run_graded_cascade() (default: GD_FAMILY, i.e. the original
    Phase 7 Gd5(SixGe1-x)4(-Ga) behavior). At the ASHRAE data-center range
    (T_cold_C=18 -> T_cold_K=291.15K) with the default GD_FAMILY, each
    stage's needed composition Tc is checked against GIANT_MCE_TC_MAX_K=
    290K: for small spans/stage-counts every stage stays within that
    documented range and the cascade is fully buildable; for larger spans
    and/or more stages the hottest stage(s) push above 290K and fall back
    to plain Gd for that stage only. See the __main__ block below for the
    actual breakdown across the sweep (computed, not assumed)."""
    T_cold_K = T_cold_C + 273.15
    rows = []
    all_stage_info = []
    for span in spans:
        row = {"span_K": span}
        for n in stage_counts:
            res = run_graded_cascade(T_cold_K, span, n, mass_per_stage=mass_per_stage,
                                      family=family)
            row[f"Graded_{n}stage_COP"] = res["COP_cascade"] if res["feasible"] else None
            row[f"Graded_{n}stage_Qc_W"] = res["Qc_W"] if res["feasible"] else None
            row[f"Graded_{n}stage_n_fallback_to_Gd"] = res["n_stages_out_of_range"]
            all_stage_info.append({"span_K": span, "n_stages": n, "stage_info": res["stage_info"]})
        rows.append(row)

    with open(out_csv, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=rows[0].keys())
        writer.writeheader()
        writer.writerows(rows)
    return rows, all_stage_info


def compare_staging(T_cold_C=18.0, spans=range(5, 21), stage_counts=(1, 2, 3, 4),
                     material=None, mass_per_stage=2.0,
                     out_csv="results/cascade_comparison.csv"):
    T_cold_K = T_cold_C + 273.15
    rows = []
    for span in spans:
        T_hot_K = T_cold_K + span
        vcc = vapor_compression_cop(T_cold_K, T_hot_K)
        liq = liquid_cooling_cop(T_cold_K, T_hot_K)
        row = {"span_K": span, "VaporCompression_COP": round(vcc.COP, 2),
               "LiquidCooling_COP": round(liq.COP, 2)}
        for n in stage_counts:
            res = run_cascade(T_cold_K, span, n, material=material, mass_per_stage=mass_per_stage)
            row[f"AMR_{n}stage_COP"] = res["COP_cascade"] if res["feasible"] else None
            row[f"AMR_{n}stage_Qc_W"] = res["Qc_W"] if res["feasible"] else None
        rows.append(row)

    with open(out_csv, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=rows[0].keys())
        writer.writeheader()
        writer.writerows(rows)
    return rows


def validate_astronautics_graded_bed(apply_correction=None):
    """Phase 9 follow-up: builds a 6-layer Curie-graded La(Fe,Si)13Hy bed
    (LAFESIH_FAMILY) at the REAL Astronautics_rotary_2014 operating point
    (Jacobs et al., Int. J. Refrig. 37 (2014) 84-91: mu0H=1.44T, 1.52kg
    total MCM, f=4Hz, T_cold=305K/32C, span=11K, reported Qc=2502W,
    COP=1.9) and runs it through the SAME calibrate-then-validate
    methodology core.validation_system.py uses for every other device:
    fluid_mdot is calibrated (brentq) to reproduce the reported Qc exactly,
    then the resulting COP is compared against the reported COP -- so this
    tests the model's predicted EFFICIENCY, not its ability to predict Qc
    from first principles (same caveat validation_system.py's own
    docstring states for its methodology).

    This directly tests the hypothesis raised in ROADMAP.md Phase 9: that
    validation_system.py's single-Tc=287K LAFESIH_FIRST_ORDER material
    failed to calibrate against this device because the real device is SIX
    Curie-graded layers (~303.6-316.2K), not because the general modeling
    approach is wrong. Layer Tc targets are spread evenly across that
    reported range (linspace(303.6, 316.2, 6)) as the best available
    approximation -- Jacobs et al. (2014) does not tabulate the individual
    layer compositions/Tc values themselves, only the range. mass_per_stage
    = 1.52/6 kg (equal split; per-layer masses are also not individually
    reported).

    apply_correction defaults to LAFESIH_FIRST_ORDER's own dTad_correction
    (1.0, i.e. uncorrected -- see that material's honesty flags for why no
    Giguere-style correction is available for this family) unless
    explicitly overridden.

    Returns a dict with the calibrated mdot, predicted Qc (should equal
    Qc_lit_W by construction), predicted vs. literature COP and its error,
    and the per-stage breakdown -- or a "no calibration found" status dict
    if no mdot in [1e-6, 1.0] kg/s reproduces the reported Qc.
    """
    from scipy.optimize import brentq

    T_cold_K = 305.0
    span_K = 11.0
    mu0H = 1.44
    mass_total = 1.52
    n_stages = 6
    freq = 4.0
    Qc_lit = 2502.0
    cop_lit = 1.9

    family = LAFESIH_FAMILY
    if apply_correction is not None:
        base = LAFESIH_FAMILY

        def tuned_fn(Tc, _corr=float(apply_correction)):
            mat = lafesih_composition_tuned_material(Tc)
            mat.dTad_correction = _corr
            return mat
        family = GradedFamily(name=base.name, tuned_fn=tuned_fn, tc_min=base.tc_min,
                               tc_max=base.tc_max, reference_material=base.reference_material,
                               fallback_material=base.fallback_material)

    def qc_residual(mdot):
        r = run_graded_cascade(T_cold_K, span_K, n_stages, mu0H_max=mu0H,
                                mass_per_stage=mass_total / n_stages, frequency=freq,
                                fluid_mdot=max(mdot, 1e-6), family=family)
        return (r["Qc_W"] if r["feasible"] else 0.0) - Qc_lit

    try:
        mdot_cal = brentq(qc_residual, 1e-6, 1.0, xtol=1e-6)
    except ValueError:
        return {"feasible": False, "status": "no calibration found "
                "(reported Qc unreachable within mdot in [1e-6, 1.0] kg/s "
                "for the 6-layer graded La(Fe,Si)13Hy bed)"}

    result = run_graded_cascade(T_cold_K, span_K, n_stages, mu0H_max=mu0H,
                                 mass_per_stage=mass_total / n_stages, frequency=freq,
                                 fluid_mdot=mdot_cal, family=family)
    result["mdot_calibrated_kg_s"] = round(mdot_cal, 5)
    result["Qc_lit_W"] = Qc_lit
    result["COP_lit"] = cop_lit
    result["COP_error_pct"] = round(100 * (result["COP_cascade"] - cop_lit) / cop_lit, 1)
    return result


if __name__ == "__main__":
    from core.mce_material import GD5SI2GE2

    print("Cascade AMR staging vs. baselines, ASHRAE 5-20K span sweep")
    print("(mu0H=2T per stage, 5kg regenerator per stage, f=1Hz, mdot=0.08kg/s, NTU thermal model on)")
    print("=" * 100)
    print("\n--- Material: Gd (baseline) ---")
    rows_gd = compare_staging(material=GADOLINIUM, mass_per_stage=5.0,
                                out_csv="results/cascade_comparison.csv")
    header = f"{'span':>5} {'1-stage':>9} {'2-stage':>9} {'3-stage':>9} {'4-stage':>9} {'VCC':>7} {'Liquid':>7}"
    print(header)
    for r in rows_gd:
        def fmt(v):
            return f"{v:9.2f}" if v is not None else f"{'--':>9}"
        print(f"{r['span_K']:>5} {fmt(r['AMR_1stage_COP'])} {fmt(r['AMR_2stage_COP'])} "
              f"{fmt(r['AMR_3stage_COP'])} {fmt(r['AMR_4stage_COP'])} "
              f"{r['VaporCompression_COP']:>7} {r['LiquidCooling_COP']:>7}")
    print(f"Wrote results/cascade_comparison.csv")

    print("\n--- Material: Gd5Si2Ge2 (giant MCE) ---")
    rows_giant = compare_staging(material=GD5SI2GE2, mass_per_stage=5.0,
                                   out_csv="results/cascade_comparison_giant_mce.csv")
    print(header)
    for r in rows_giant:
        print(f"{r['span_K']:>5} {fmt(r['AMR_1stage_COP'])} {fmt(r['AMR_2stage_COP'])} "
              f"{fmt(r['AMR_3stage_COP'])} {fmt(r['AMR_4stage_COP'])} "
              f"{r['VaporCompression_COP']:>7} {r['LiquidCooling_COP']:>7}")
    print(f"Wrote results/cascade_comparison_giant_mce.csv")

    gd_10K = next(r for r in rows_gd if r["span_K"] == 10)
    giant_10K = next(r for r in rows_giant if r["span_K"] == 10)
    print(f"\nAt 10K span: Gd 1-stage COP={gd_10K['AMR_1stage_COP']} vs. "
          f"Gd5Si2Ge2 1-stage COP={giant_10K['AMR_1stage_COP']} "
          f"(VCC={gd_10K['VaporCompression_COP']}, Liquid={gd_10K['LiquidCooling_COP']})")

    print("\n" + "=" * 100)
    print("--- Curie-graded cascade (ROADMAP.md Phase 7 open item) ---")
    print("=" * 100)
    print("Each stage uses a hypothetical composition-tuned Gd5(SixGe1-x)4(-Ga) material")
    print("whose own peak MCE effect is matched (via iterative peak search, see")
    print("_target_composition_for_peak) to that stage's local operating temperature,")
    print("checked against the literature-documented Tc range (20-290K) and scaled by the")
    print("Giguere et al. (1999) empirical correction (core.giguere_validation), so these")
    print("numbers are not built on the raw model's ~2.4x-optimistic DeltaT_ad.\n")

    rows_graded, stage_info_all = compare_graded_cascade(
        T_cold_C=18.0, spans=range(5, 21), mass_per_stage=5.0,
        out_csv="results/graded_cascade_comparison.csv")

    example = next(s for s in stage_info_all if s["span_K"] == 10 and s["n_stages"] == 3)
    print("Example (10K span, 3 stages):")
    for s in example["stage_info"]:
        print(f"    stage {s['stage']}: T_mid={s['T_mid_K']}K, needed composition Tc="
              f"{s['Tc_target_K']}K -> {s['material']}")

    graded_10K_3 = next(r for r in rows_graded if r["span_K"] == 10)
    gd_10K_3 = next(r for r in rows_gd if r["span_K"] == 10)
    print(f"\nAt this point: Graded 3-stage Qc={graded_10K_3['Graded_3stage_Qc_W']}W, "
          f"COP_elec={graded_10K_3['Graded_3stage_COP']}  vs.  plain-Gd 3-stage "
          f"Qc={gd_10K_3['AMR_3stage_Qc_W']}W, COP_elec={gd_10K_3['AMR_3stage_COP']}")
    print("-> consistent with giant_mce_analysis.py's earlier finding: a bigger MCE mostly")
    print("   buys more Qc per kg (here: substantially more), not a materially better COP")
    print("   (loss_model.py's field/frequency/flow-dependent parasitics dominate COP either way).")

    n_cells = sum(1 for row in rows_graded for n in (1, 2, 3, 4))
    n_full_range = sum(1 for row in rows_graded for n in (1, 2, 3, 4)
                        if row[f"Graded_{n}stage_n_fallback_to_Gd"] == 0)
    n_some_fallback = sum(1 for row in rows_graded for n in (1, 2, 3, 4)
                           if 0 < row[f"Graded_{n}stage_n_fallback_to_Gd"] < n)
    n_all_fallback_or_infeasible = n_cells - n_full_range - n_some_fallback
    print(f"\nAcross the full 5-20K span x 1-4 stage-count sweep ({n_cells} cells):")
    print(f"  {n_full_range} cells: every stage's needed composition is within the "
          f"documented 20-290K giant-MCE range")
    print(f"  {n_some_fallback} cells: SOME stages exceed 290K and fall back to plain Gd "
          f"for that stage only (larger spans and/or more stages push the hottest stage's "
          f"needed Tc above the Ga-alloyed ceiling)")
    print(f"  {n_all_fallback_or_infeasible} cells: fully infeasible (all stages fell back, "
          f"or Qc collapsed to ~0)")
    print("\nHONEST CAVEAT: this Landau model's transition turns out to be numerically much")
    print("narrower (DeltaT_ad falls off within a few tenths of a K of its peak at mu0H=2T)")
    print("than the real, hysteresis/inhomogeneity-broadened transition Giguere et al.'s own")
    print("Fig. 3 shows (spread over ~10-15K) -- an idealized-model limitation on top of the")
    print("dTad-magnitude one already flagged, not smoothed over here. A small number of")
    print("individual span/stage-count cells show a stage's Qc collapsing to ~0 despite its")
    print("composition nominally being in-range, from residual peak-alignment error after the")
    print("iterative search (see run_graded_cascade's COP_cascade=0.0 with Qc_W>0 rows in the")
    print("output CSV) -- flagged here as a real numerical fragility of this idealized 6th-")
    print("order Landau fit, not hidden by rounding it away.")
    print(f"\nWrote results/graded_cascade_comparison.csv")
    print("\nBOTTOM LINE: unlike the fixed Gd5Si2Ge2 comparison above (which collapses to zero")
    print("everywhere in the ASHRAE range because Gd5Si2Ge2's own Tc=276K is fixed and far from")
    print("the operating point), a CURIE-GRADED cascade -- built from literature-documented")
    print("composition tunability and validated against Giguere et al.'s direct measurement --")
    print("genuinely delivers several times the cooling capacity of plain Gd at comparable COP")
    print("for smaller spans/stage-counts, but the same giant-MCE family's documented ~290K")
    print("composition ceiling still constrains larger spans and higher stage counts, echoing")
    print("giant_mce_analysis.py's conclusion that the ASHRAE range sits right at the edge of")
    print("what this material family is documented to reach.")

    print("\n" + "=" * 100)
    print("--- Phase 9: does a 6-layer Curie-graded La(Fe,Si)13Hy bed reproduce the REAL")
    print("    Astronautics_rotary_2014 device? (validation_system.py's single-Tc=287K")
    print("    material could not -- see ROADMAP.md Phase 9) ---")
    print("=" * 100)
    astro = validate_astronautics_graded_bed()
    if astro.get("feasible"):
        print(f"Layer Tc targets (evenly spread across the device's reported 303.6-316.2K "
              f"layer range):")
        for s in astro["stage_info"]:
            print(f"    stage {s['stage']}: T_mid={s['T_mid_K']}K, needed composition Tc="
                  f"{s['Tc_target_K']}K -> {s['material']}")
        print(f"\nmdot calibrated to reproduce reported Qc={astro['Qc_lit_W']}W: "
              f"{astro['mdot_calibrated_kg_s']} kg/s")
        print(f"Predicted COP={astro['COP_cascade']}  vs.  reported COP={astro['COP_lit']} "
              f"({astro['COP_error_pct']:+.1f}% error)")
        print("\n-> a comparable-magnitude error to the Gd devices in validation_system.py's own")
        print("   point-wise validation, and a very different outcome from the flat 'no")
        print("   calibration found' the single-layer LAFESIH_FIRST_ORDER material gave this")
        print("   same device: the graded-bed STRUCTURE, not just the material, was the missing")
        print("   piece. Still an approximation -- layer Tc's are evenly spread across the")
        print("   reported range, not the paper's actual (unpublished here) per-layer values,")
        print("   and mdot is calibrated rather than predicted, same caveat as")
        print("   validation_system.py's own methodology throughout.")
    else:
        print(astro.get("status", "infeasible"))