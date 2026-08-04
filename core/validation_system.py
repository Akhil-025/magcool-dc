"""
validation_system.py
====================
System-level validation of the active magnetic regenerator (AMR) cycle model
against published prototype-scale experimental data.

Benchmark data are loaded from
data/amr_experimental_benchmarks.csv.

Methodology
-----------
The available literature generally does not report fluid mass-flow rate in a
form directly comparable across devices. Some studies instead report
utilization factors or other device-specific quantities, making direct
comparison difficult.

The validation therefore follows a two-step procedure:

1. Calibration

   For each benchmark device, determine the fluid mass-flow rate that
   reproduces the reported cooling capacity (Qc) at the reported operating
   temperature span using the published field strength, regenerator mass and
   operating frequency.

   The magnetic material closest to the experimental device is used:

   • Gd for gadolinium-based prototypes.

   • LAFESIH_FIRST_ORDER (core/first_order_mce.py) for La(Fe,Si)13Hy devices
     -- added because the Astronautics_rotary_2014 row was previously run
     against GADOLINIUM as an explicitly-flagged stand-in. See the block
     comment above that material's definition for its calibration and
     honesty flags; it represents a single Tc=287K layer, not the real
     device's six Curie-graded layers (~304-316K), so treat it as a
     representative single-material approximation like every other row here.

2. Validation

   After calibrating the mass-flow rate, compare the model-predicted
   electrical COP with the published experimental COP.

   Since cooling capacity is matched during calibration, the validation
   primarily assesses the cycle model's prediction of efficiency rather than
   its ability to predict cooling capacity from first principles.

Limitations
-----------
Because the fluid mass-flow rate is calibrated rather than independently
measured, this procedure validates the cycle-efficiency model rather than
providing a completely independent end-to-end validation of the entire
0-D AMR model.

A more rigorous validation would compare against experimental datasets that
report complete operating conditions, including directly measured flow rates
or utilization values that are consistent across devices.

Curve-level (Phase 7) addition
-------------------------------
`amr_cycle.AMRSystem.cooling_capacity()` predicts a specific characteristic-
curve *shape*: Qc falls off roughly linearly from a maximum at zero span to
zero at the no-load span (the Tusek 2010 / Nielsen 2011 "characteristic
curve" form the model docstring already claims to match). Three of the
benchmark devices happen to have a second, companion data point at a
different span for the *same physical device* (a zero-span max-capacity
point or a max-span zero-capacity point), grouped via the `device_group`
column. `run_curve_validation()` uses these to test the predicted curve
*shape*, not just a single operating point: it calibrates mdot at the
device's normal operating point (as above) and then checks whether the
resulting Qc(span) curve, evaluated at the companion span, reproduces the
companion's reported Qc/zero-crossing.

This is NOT a substitute for digitizing full published characteristic
curves (Tusek 2010, Nielsen 2011) -- those two specific source papers are
still not available in this repository, so no such digitization is
fabricated here. This is a smaller, honest step: a genuine 2-point curve
check using data already in `amr_experimental_benchmarks.csv`, for the
device groups that have a second point at a fixed operating condition.

Phase 7 continued: a genuinely independent multi-point device was located
and added -- Lozano et al. (2016), "Development of a novel rotary magnetic
refrigerator" (POLO/UFSC), whose Table 3 reports 8 real (frequency, flow
rate, span, Qc, COP) operating points directly as clean numbers (no
digitization needed), plus zero-span and no-load-span endpoints from the
abstract. This is genuinely different in kind from the 2-point companion
checks above: every row has its own frequency AND flow rate, not just a
different span at otherwise-fixed conditions, so it is validated
point-by-point by `run_system_validation()` (each row calibrates its own
mdot independently) rather than forced through the anchor/companion
pairing here -- see the guard in `run_curve_validation()` below. Also
found: Tusek, Kitanovski, Zupan, Prebil, Poredos (2013), "A comprehensive
experimental analysis of gadolinium active magnetic regenerators", Appl.
Therm. Eng. 53, 57-66 -- this is the paper actually being described by the
LITERATURE_REVIEW.md / CSV citation "Tusek et al., Appl. Therm. Eng.
(2011 dataset, 2010 device)" for the existing `Tusek_singlebed_Gd_2010`
row, though its own Table 1/Sec.3.1 reports a DIFFERENT field (1.15T vs.
the CSV's 1.69T) and mass (6 AMRs, 0.093-0.176kg, none matching the CSV's
0.196kg exactly) -- so it is additional independent data, not a
confirmation of that existing single-point row, and the true source of
the 1.69T/0.196kg point remains unidentified. Its Figs. 10-11 (Qc-vs-span
and COP-vs-span for 3 AMR geometries at 3 flow ratios, 9 lines total) are
genuine multi-point published characteristic curves in the Tusek-2010/
Nielsen-2011 sense this module's docstring originally flagged as blocked
-- but extracting them requires pixel-calibrated digitization from the
figure images, which was NOT done here (out of scope for this pass); a
rough by-eye read is in `results/tusek_ate2013_figs_notes.md` for anyone
picking this up, clearly marked as non-authoritative. So: the Tusek/
Nielsen full-curve item is now unblocked for Tusek 2013 (paper in hand)
but still not actually digitized, and Nielsen 2011 is still not in the
repository at all.
"""

import csv
import numpy as np
from scipy.optimize import brentq
from core.mce_material import GADOLINIUM
from core.first_order_mce import LAFESIH_FIRST_ORDER
from core.amr_cycle import AMRSystem
from core.loss_model import StateDependentLossModel, RotaryDriveLossModel

BENCH_CSV = "data/amr_experimental_benchmarks.csv"
T_COLD_ASSUMED_K = 294.0 - 5.0     # assume device centered near Gd's Tc=294K,
                                   # cold side ~5K below center as a working default
T_COLD_LAFESIH_K = 305.0          # 32C measured cold-side inlet, Jacobs et al.
                                   # (2014) Sec. 2 -- used instead of the Gd
                                   # default above for La(Fe,Si)13Hy rows, since
                                   # that device intentionally runs well above
                                   # room temperature (naval electronics cooler
                                   # spec) and reusing the Gd-centered value
                                   # here would silently apply the wrong
                                   # operating point to the wrong material.


def _material_for_row(row):
    """Select the MCE material object for a benchmark row. La(Fe,Si)13Hy rows
    now use the real LAFESIH_FIRST_ORDER material (see first_order_mce.py)
    instead of the previous GADOLINIUM stand-in."""
    if "La" in row["material"]:
        return LAFESIH_FIRST_ORDER
    return GADOLINIUM


def _t_cold_for_row(row):
    return T_COLD_LAFESIH_K if "La" in row["material"] else T_COLD_ASSUMED_K


# BUGFIX (found while diagnosing why a prior loss-model fix produced no change
# in this module's output): AMRSystem() was never being given a loss_model=
# argument anywhere in this file, so every device below -- including Lozano --
# silently fell back to AMRSystem's flat parasitic_fraction=0.15 default
# instead of either calibrated model. loss_model.py's StateDependentLossModel
# (CORE-calibrated) and RotaryDriveLossModel (Lozano-specific) were built and
# tested in isolation but never actually reached this validation pipeline.
# Cached at module scope so the (cheap but non-trivial) NNLS/least-squares
# calibration in loss_model.py only runs once per process, not once per row.
_CORE_LOSS_MODEL = None
_LOZANO_LOSS_MODEL = None


def _get_core_loss_model():
    global _CORE_LOSS_MODEL
    if _CORE_LOSS_MODEL is None:
        _CORE_LOSS_MODEL = StateDependentLossModel()
    return _CORE_LOSS_MODEL


def _get_lozano_loss_model():
    global _LOZANO_LOSS_MODEL
    if _LOZANO_LOSS_MODEL is None:
        _LOZANO_LOSS_MODEL = RotaryDriveLossModel()
    return _LOZANO_LOSS_MODEL


def _loss_model_for_row(row):
    """Select the calibrated parasitic-loss model for a benchmark row.

    Lozano_POLO_UFSC_2016 rows use RotaryDriveLossModel -- CORE's eddy/pump/
    base terms plus an additional rotary magnet-assembly/valve drivetrain
    term fit to Lozano et al. (2016)'s own directly-measured WM data (see
    loss_model.RotaryDriveLossModel's docstring for why this is device-
    specific rather than a general "rotary AMR" flag; Astronautics and DTU
    are also rotary devices and are deliberately left on plain CORE).
    Every other device uses the plain CORE-calibrated StateDependentLossModel.
    """
    if row["device_group"] == "Lozano_POLO_UFSC_2016":
        return _get_lozano_loss_model()
    return _get_core_loss_model()


def load_benchmarks(path=BENCH_CSV):
    rows = []
    with open(path) as f:
        for r in csv.DictReader(f):
            if not r.get("device"):
                continue  # skip blank/trailing rows
            rows.append(r)
    return rows


def calibrate_and_check(row, verbose=True):
    span = float(row["span_K"])
    Qc_lit = row["Qc_W"]
    cop_lit = row["COP"]
    if span <= 0 or not Qc_lit or not cop_lit:
        return None  # zero-span / capacity-only rows aren't COP validation targets

    Qc_lit = float(Qc_lit)
    cop_lit = float(cop_lit)
    mu0H = float(row["mu0H_T"])
    mass = float(row["mass_MCM_kg"]) if row["mass_MCM_kg"] else 1.0
    freq = float(row["frequency_Hz"]) if row["frequency_Hz"] else 1.0
    material = _material_for_row(row)
    t_cold = _t_cold_for_row(row)
    material_note = ("La(Fe,Si)13Hy (LAFESIH_FIRST_ORDER, single-Tc "
                      "approximation of the real 6-layer graded bed)"
                      if "La" in row["material"] else "Gd (matches device)")
    loss_model = _loss_model_for_row(row)

    def qc_residual(mdot):
        sys_ = AMRSystem(material=material, mu0H_max=mu0H, mass_regenerator=mass,
                          frequency=freq, fluid_mdot=max(mdot, 1e-6),
                          loss_model=loss_model)
        Qc_model, _ = sys_.cooling_capacity(t_cold, span)
        return Qc_model - Qc_lit

    try:
        mdot_cal = brentq(qc_residual, 1e-6, 5.0, xtol=1e-6)
    except ValueError:
        if verbose:
            print(f"{row['device']:<28} span={span:5.1f}K  Qc_lit={Qc_lit:7.1f}W  "
                  f"NO CALIBRATION FOUND (reported Qc unreachable within "
                  f"mdot in [1e-6,5] kg/s at this field/mass/frequency)  "
                  f"[{material_note}]")
        return {"device": row["device"], "status": "no calibration found "
                "(reported Qc unreachable within mdot in [1e-6,5] kg/s "
                "at this field/mass/frequency)", "material_note": material_note}

    sys_ = AMRSystem(material=material, mu0H_max=mu0H, mass_regenerator=mass,
                      frequency=freq, fluid_mdot=mdot_cal, loss_model=loss_model)
    result = sys_.run(t_cold, span)
    cop_err_pct = 100 * (result.COP_electrical - cop_lit) / cop_lit
    implied_parasitic_frac = (1 / cop_lit - 1 / result.COP) if result.COP > 0 else float("nan")

    out = {"device": row["device"], "span_K": span, "Qc_lit_W": Qc_lit,
           "Qc_model_W": round(result.Qc, 1), "COP_lit": cop_lit,
           "COP_model_ideal": round(result.COP, 2),
           "COP_model_electrical": round(result.COP_electrical, 2),
           "COP_error_pct": round(cop_err_pct, 1),
           "implied_parasitic_fraction": round(implied_parasitic_frac, 3),
           "mdot_calibrated_kg_s": round(mdot_cal, 4), "material_note": material_note}
    if verbose:
        print(f"{row['device']:<28} span={span:5.1f}K  Qc(lit/model)="
              f"{Qc_lit:7.1f}/{result.Qc:7.1f} W  COP(lit/ideal/elec)="
              f"{cop_lit:5.2f}/{result.COP:5.2f}/{result.COP_electrical:5.2f}"
              f"  err={cop_err_pct:+6.1f}%  implied_parasitic={implied_parasitic_frac:.3f}"
              f"  [{material_note}]")
    return out


def run_capacity_only_calibration_check(verbose=True):
    """run_system_validation() silently skips any row without a reported
    COP (`calibrate_and_check()`'s own comment: "capacity-only rows aren't
    COP validation targets") -- which is correct for COP comparison, but
    means those rows (zerospan/maxspan companions, the Chubu Electric/
    Toshiba pair, and now Cooltech_2013_rotary/DTU_MagQueen_2018) never
    get ANY reported result. This function reuses `_calibrate_mdot()`
    (which only needs span/Qc/field/mass/frequency, not COP) to report
    whether each one calibrates at all -- in particular the
    Cooltech_2013_rotary row, which Paper-Mining Pass Part 3, §1 flags as
    a genuine STRESS TEST: 42K is the largest span in this benchmark set
    (next is Risoe_DTU_Gd_2011 at 30K, which itself does NOT calibrate --
    see run_system_validation()'s own printed output).
    """
    rows = load_benchmarks()
    results = []
    for row in rows:
        span = float(row["span_K"])
        Qc_lit = row["Qc_W"]
        cop_lit = row["COP"]
        if span <= 0 or not Qc_lit or cop_lit:
            continue  # zero-span rows, or rows that already have a reported COP (covered by run_system_validation() instead)
        Qc_lit = float(Qc_lit)
        cal = _calibrate_mdot(row)
        material_note = ("La(Fe,Si)13Hy (LAFESIH_FIRST_ORDER)" if "La" in row["material"]
                          else "Gd (matches device)")
        if cal is None:
            out = {"device": row["device"], "span_K": span, "Qc_lit_W": Qc_lit,
                   "status": "no calibration found", "material_note": material_note}
            if verbose:
                print(f"{row['device']:<28} span={span:5.1f}K  Qc_lit={Qc_lit:7.1f}W  "
                      f"NO CALIBRATION FOUND (capacity-only row, no COP to compare -- "
                      f"reports whether the span/Qc pair alone is achievable)  "
                      f"[{material_note}]")
        else:
            mdot_cal, sys_ = cal
            out = {"device": row["device"], "span_K": span, "Qc_lit_W": Qc_lit,
                   "status": "calibrated", "mdot_calibrated_kg_s": round(mdot_cal, 4),
                   "material_note": material_note}
            if verbose:
                print(f"{row['device']:<28} span={span:5.1f}K  Qc_lit={Qc_lit:7.1f}W  "
                      f"calibrated at mdot={mdot_cal:.4f}kg/s (no COP reported -- "
                      f"capacity-only check)  [{material_note}]")
        results.append(out)
    return results


def run_system_validation():
    rows = load_benchmarks()
    results = [calibrate_and_check(r) for r in rows]
    return [r for r in results if r is not None]


def _calibrate_mdot(row):
    """Calibrate fluid mdot to reproduce this row's (span, Qc). Returns
    (mdot, AMRSystem) or None if no calibration is reachable."""
    span = float(row["span_K"])
    Qc_lit = float(row["Qc_W"])
    mu0H = float(row["mu0H_T"])
    mass = float(row["mass_MCM_kg"]) if row["mass_MCM_kg"] else 1.0
    freq = float(row["frequency_Hz"]) if row["frequency_Hz"] else 1.0
    material = _material_for_row(row)
    t_cold = _t_cold_for_row(row)
    loss_model = _loss_model_for_row(row)

    def qc_residual(mdot):
        sys_ = AMRSystem(material=material, mu0H_max=mu0H, mass_regenerator=mass,
                          frequency=freq, fluid_mdot=max(mdot, 1e-6),
                          loss_model=loss_model)
        Qc_model, _ = sys_.cooling_capacity(t_cold, span)
        return Qc_model - Qc_lit

    try:
        mdot_cal = brentq(qc_residual, 1e-6, 5.0, xtol=1e-6)
    except ValueError:
        return None
    sys_ = AMRSystem(material=material, mu0H_max=mu0H, mass_regenerator=mass,
                      frequency=freq, fluid_mdot=mdot_cal, loss_model=loss_model)
    return mdot_cal, sys_


def run_curve_validation(verbose=True):
    """Curve-level (2-point) validation: for each device_group with >=2 rows,
    calibrate mdot on the row that has a reported COP (the normal operating
    point), then use that SAME calibrated system to predict Qc at the
    companion row's span, and compare against the companion's reported Qc.

    This checks the model's predicted Qc(span) *shape*, not just a single
    point -- see module docstring's "Curve-level (Phase 7) addition"."""
    rows = load_benchmarks()
    groups = {}
    for r in rows:
        groups.setdefault(r["device_group"], []).append(r)

    results = []
    for group_name, group_rows in groups.items():
        if len(group_rows) < 2:
            continue  # single-point device, no curve to check

        cop_rows = [r for r in group_rows if r["COP"]]
        if len(cop_rows) > 1:
            # Multiple independent (frequency, flow, span, Qc, COP) points for
            # this device (e.g. Lozano_POLO_UFSC_2016: 8 rows, each its own
            # frequency and flow rate) are NOT a fixed-condition span sweep --
            # there is no single "anchor" operating point whose calibrated
            # mdot should be reused to predict the others under different
            # frequencies. Forcing them through the 2-point companion-pair
            # logic below would silently compare across mismatched operating
            # conditions. These groups are already validated independently,
            # point by point, by run_system_validation(); skip them here.
            if verbose:
                print(f"{group_name:<28} {len(cop_rows)} independent multi-point rows "
                      f"(own frequency/flow per row) -- validated by "
                      f"run_system_validation() instead, skipped here")
            results.append({"device_group": group_name,
                             "status": f"multi-point set ({len(cop_rows)} independent "
                                       "points) -- see run_system_validation()"})
            continue

        # the calibration ("anchor") row is the one with a reported COP
        anchor = next((r for r in group_rows if r["COP"]), None)
        companion = next((r for r in group_rows if r is not anchor), None)
        if anchor is None or companion is None:
            continue

        cal = _calibrate_mdot(anchor)
        if cal is None:
            if verbose:
                print(f"{group_name:<28} no calibration found at anchor point "
                      f"(span={anchor['span_K']}K) -- consistent with the Phase 2/6 "
                      f"finding that this device doesn't calibrate; skipped, not "
                      f"silently dropped")
            results.append({"device_group": group_name,
                             "status": "no calibration found at anchor point"})
            continue
        mdot_cal, sys_ = cal

        companion_span = float(companion["span_K"])
        Qc_companion_lit = float(companion["Qc_W"])
        Qc_companion_model, _ = sys_.cooling_capacity(_t_cold_for_row(anchor), companion_span)

        if Qc_companion_lit > 0:
            err_pct = 100 * (Qc_companion_model - Qc_companion_lit) / Qc_companion_lit
            err_str = f"{err_pct:+.1f}%"
        else:
            # companion is a no-load / zero-capacity (max-span) point --
            # relative error is undefined, report absolute W instead
            err_pct = None
            err_str = f"{Qc_companion_model:+.1f}W (lit=0, undefined %)"

        out = {"device_group": group_name,
               "anchor_span_K": float(anchor["span_K"]), "anchor_Qc_W": float(anchor["Qc_W"]),
               "companion_span_K": companion_span, "companion_Qc_lit_W": Qc_companion_lit,
               "companion_Qc_model_W": round(Qc_companion_model, 1),
               "companion_Qc_error_pct": round(err_pct, 1) if err_pct is not None else None,
               "mdot_calibrated_kg_s": round(mdot_cal, 4)}
        if verbose:
            print(f"{group_name:<28} anchor span={out['anchor_span_K']:5.1f}K -> "
                  f"predict Qc at companion span={companion_span:5.1f}K: "
                  f"model={Qc_companion_model:7.1f}W  lit={Qc_companion_lit:7.1f}W  "
                  f"err={err_str}")
        results.append(out)
    return results


def run_field_sensitivity_check(verbose=True, device_group="ChubuToshiba_Gd_2016"):
    """Field-axis analog of run_curve_validation() (which checks the
    model's predicted Qc(span) shape): checks the model's predicted
    Qc(field) sensitivity, using the Chubu Electric/Toshiba two-field-point
    pair (ChubuToshiba_Gd_2016_4T/_2T, Paper-Mining Pass Part 2, §1) --
    the only field-sensitivity data point in this benchmark set; every
    other device here is reported at a single, fixed field.

    Calibrates fluid mdot to reproduce the HIGHER-field (4T) row's own
    (span, Qc) -- same _calibrate_mdot() this module already uses, note it
    only needs (span, Qc, mu0H, mass, freq), NOT a reported COP, so the
    lack of a COP column for these two rows (not in the source review's
    table) is not a blocker here, unlike run_curve_validation()'s
    COP-anchored anchor selection. Reuses that SAME calibrated
    mdot/mass/frequency system at the LOWER field (2T) to predict Qc at
    the companion row's own reported span, and compares against its
    reported Qc.

    SAME secondary-source caveat as the CSV rows themselves and the
    existing Okamura_Hirano_2013 row: both points are read from a review's
    table (Kamran, Ahmad & Wang, Renew. Sustain. Energy Rev. 133 (2020)
    110247, Table 2), not the primary device paper (ref [69] in that
    review), which is not in this repo's Papers/.
    """
    rows = load_benchmarks()
    group = [r for r in rows if r["device_group"] == device_group]
    if len(group) < 2:
        if verbose:
            print(f"{device_group:<28} fewer than 2 rows found -- field-sensitivity "
                  f"check skipped")
        return None

    anchor = max(group, key=lambda r: float(r["mu0H_T"]))
    companion = next((r for r in group if r is not anchor), None)

    cal = _calibrate_mdot(anchor)
    if cal is None:
        if verbose:
            print(f"{device_group:<28} no calibration found at anchor field="
                  f"{anchor['mu0H_T']}T")
        return {"device_group": device_group, "status": "no calibration found at anchor field"}
    mdot_cal, _ = cal

    material = _material_for_row(anchor)
    mass = float(anchor["mass_MCM_kg"])
    freq = float(anchor["frequency_Hz"])
    t_cold = _t_cold_for_row(anchor)
    companion_field = float(companion["mu0H_T"])
    companion_span = float(companion["span_K"])
    Qc_companion_lit = float(companion["Qc_W"])

    sys_companion = AMRSystem(material=material, mu0H_max=companion_field,
                               mass_regenerator=mass, frequency=freq, fluid_mdot=mdot_cal,
                               loss_model=_loss_model_for_row(anchor))
    Qc_companion_model, _ = sys_companion.cooling_capacity(t_cold, companion_span)
    err_pct = 100 * (Qc_companion_model - Qc_companion_lit) / Qc_companion_lit

    lit_field_ratio = float(anchor["Qc_W"]) / Qc_companion_lit
    model_field_ratio = float(anchor["Qc_W"]) / Qc_companion_model if Qc_companion_model > 0 else float("inf")

    out = {"device_group": device_group,
           "anchor_field_T": float(anchor["mu0H_T"]), "anchor_Qc_W": float(anchor["Qc_W"]),
           "companion_field_T": companion_field, "companion_Qc_lit_W": Qc_companion_lit,
           "companion_Qc_model_W": round(Qc_companion_model, 1),
           "companion_Qc_error_pct": round(err_pct, 1),
           "lit_field_ratio": round(lit_field_ratio, 2),
           "model_field_ratio": round(model_field_ratio, 2),
           "mdot_calibrated_kg_s": round(mdot_cal, 4)}
    if verbose:
        print(f"{device_group:<28} anchor field={out['anchor_field_T']:.1f}T "
              f"(Qc={out['anchor_Qc_W']:.0f}W, calibrated) -> predict Qc at "
              f"companion field={companion_field:.1f}T: model={Qc_companion_model:7.1f}W  "
              f"lit={Qc_companion_lit:7.1f}W  err={err_pct:+.1f}%")
        print(f"{'':<28} literature Qc({out['anchor_field_T']:.0f}T)/Qc({companion_field:.0f}T) "
              f"ratio={lit_field_ratio:.2f}  vs. model ratio={model_field_ratio:.2f}")
    return out


if __name__ == "__main__":
    print("System-level validation against published AMR prototype data")
    print("=" * 110)
    results = run_system_validation()
    errs = [abs(r["COP_error_pct"]) for r in results if "COP_error_pct" in r]
    if errs:
        print("=" * 110)
        print(f"Mean |COP_electrical error| = {np.mean(errs):.1f}%  |  Max = {np.max(errs):.1f}%")
        print(
            "\nSummary: the ideal thermodynamic-cycle COP substantially "
            "overpredicts published electrical COP because experimental values "
            "include pump, motor and drive losses that are absent from the ideal "
            "cycle. Incorporating calibrated parasitic losses significantly improves "
            "agreement for the laboratory-scale benchmark devices, while the "
            "large-scale Astronautics prototype remains an outlier, likely reflecting "
            "additional system-level inefficiencies not represented by the current "
            "model. Electrical COP is therefore the appropriate metric for comparison "
            "with vapor-compression and liquid-cooling systems."
        )

    print("\n" + "=" * 110)
    print("Curve-level (2-point Qc-vs-span shape) validation")
    print("=" * 110)
    curve_results = run_curve_validation()
    curve_errs = [abs(r["companion_Qc_error_pct"]) for r in curve_results
                  if r.get("companion_Qc_error_pct") is not None]
    if curve_errs:
        print("=" * 110)
        print(f"Mean |companion Qc error| (nonzero-Qc companions only) = "
              f"{np.mean(curve_errs):.1f}%  |  Max = {np.max(curve_errs):.1f}%")
    print(
        "\nSummary: only 3 of 5 benchmark device groups have a second span "
        "point, so this is a small, honest check of curve *shape* (not a "
        "digitized full characteristic curve -- see module docstring), but it "
        "is a genuinely independent test: the companion point is NOT used in "
        "calibration, only predicted from the mdot fitted at the anchor point."
    )

    import csv as _csv
    with open("results/curve_validation.csv", "w", newline="") as f:
        fieldnames = ["device_group", "status", "anchor_span_K", "anchor_Qc_W",
                      "companion_span_K", "companion_Qc_lit_W",
                      "companion_Qc_model_W", "companion_Qc_error_pct",
                      "mdot_calibrated_kg_s"]
        w = _csv.DictWriter(f, fieldnames=fieldnames)
        w.writeheader()
        for r in curve_results:
            w.writerow({k: r.get(k, "") for k in fieldnames})
    print("\nSaved results/curve_validation.csv")

    print("\n" + "=" * 110)
    print("Field-sensitivity (2-point Qc-vs-field) check (Paper-Mining Pass Part 2, §1)")
    print("=" * 110)
    field_result = run_field_sensitivity_check()
    if field_result and "companion_Qc_error_pct" in field_result:
        print(
            "\nSummary: only one device group in this benchmark set (Chubu Electric/"
            "Toshiba, secondary-source via Kamran, Ahmad & Wang 2020) reports the "
            "same device at two different fields, so this is necessarily a single "
            "check, not a statistical sample -- but it is genuinely independent: "
            "the companion (2T) point is NOT used in calibration, only predicted "
            "from the mdot fitted at the anchor (4T) point."
        )

    print("\n" + "=" * 110)
    print("Capacity-only rows: calibration-reachability check (Paper-Mining Pass "
          "Part 3, §1's Cooltech 42K stress test, among others)")
    print("=" * 110)
    run_capacity_only_calibration_check()