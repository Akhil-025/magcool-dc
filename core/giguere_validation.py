"""
giguere_validation.py
======================
Cross-checks the first-order Landau model (`first_order_mce.py`) against an
INDEPENDENT experimental dataset, as required by that module's honesty flag
#2 and by ROADMAP.md the earlier open item on the Curie-graded cascade.

Reference data
--------------
Giguère, Foldeaki, Ravi Gopal, Chahine, Bose, Frydman & Barclay,
"Direct Measurement of the 'Giant' Adiabatic Temperature Change in
Gd5Si2Ge2", Phys. Rev. Lett. 83, 2262 (1999).
(PDF present in this repo: Papers/Magnetocaloric effect and materials
physics/Direct Measurement of the "Giant" Adiabatic Temperature Change in
Gd5Si2Ge2.pdf)

What that paper actually reports (read directly from the PDF, not assumed):
  - Two Gd5Si2Ge2 samples (rod + button, both AMES-laboratory arc-melted)
    were DIRECTLY measured (not calculated from magnetization/Maxwell
    relations). At mu0*H = 7 T, the direct measurement gives
    DeltaT_ad ~ 10 K (their Fig. 3), cross-checked independently via a
    Clausius-Clapeyron latent-heat analysis of the same first-order
    transition, which gives DeltaT_ad = 9.9 K -- "excellent agreement".
  - The zero-field transition temperature, from H_crit(T) extrapolation,
    is 272.2 K (increasing field) / 278.2 K (decreasing field) -- i.e.
    straddling the Tc=276 K used for GD5SI2GE2_FIRST_ORDER in this
    codebase, which is therefore a reasonable choice, independently of
    this validation.
  - Critically, the paper ALSO computed the "indirect" value the same way
    the widely-cited "giant" MCE literature does -- integrating Maxwell
    relations over the measured magnetization curves -- and got
    DeltaT_max = 14.9 K, "about the value reported in Ref. [5]" (Pecharsky
    & Gschneidner, PRL 78, 4494 (1997), the source of the ~18 J/(kg K)
    peak DeltaS_M this codebase's Landau model is calibrated to). The
    paper's own conclusion is that this indirect/Maxwell-relation method
    OVERESTIMATES DeltaT_ad for a first-order transition (their words:
    "using the maximum DeltaS value obtained from Maxwell relations
    overestimates DeltaT_ad in case of first order transitions"), by
    about 14.9/10 = 1.49x at 7 T in their own data.

Why this matters here: `first_order_mce.py`'s (A, B, C) Landau coefficients
were calibrated to that SAME indirect/Maxwell-relation "giant" DeltaS_M
literature value (~18 J/(kg K) at 5 T, the "widely-quoted" end of the
Pecharsky & Gschneidner review range). Giguere et al.'s own finding is
exactly the mechanism honesty flag #1 already worried about in the abstract
(lattice-only C_lattice denominator, no first-order latent-heat structure)
-- this module turns that worry into an actual number, computed once here.

Result of this cross-check (see run_validation() below): the model
predicts a peak DeltaT_ad at 7 T of ~24 K, vs. the 10.0 K DIRECTLY measured
by Giguere et al. -- roughly a 2.4x overestimate, WORSE than the 1.49x
overestimate the same paper found for the "ordinary" Maxwell-relation
method it was comparing against. This is consistent with (not merely
duplicating) honesty flag #1's concern: this 0-D model's lattice-only
C_lattice denominator is uncorrected for the first-order transition's own
latent heat, which independently inflates the DeltaT_ad estimate on top of
the Maxwell-relation-vs-direct gap Giguere et al. already identified.

Note on scope: this compares the model's OWN peak (found by scanning T,
same as the module's __main__ block does, since the peak shifts off the
nominal Tc) against Giguere's directly measured 7 T value. It is a
single-field, single-composition check -- exactly the kind of thing
honesty flag #2 asked for, not a full re-validation across field and
composition (no such multi-field direct dataset exists in this repo).

 verification note: this module's docstring already said "PDF
present in this repo" before  -- re-confirmed directly again this
phase (Papers/ is now actually included in the delivered project, closing
a gap where these citations referenced a Papers/ folder that had been
omitted from earlier deliveries of this codebase). All four numbers above
(10.0 K direct, 9.9 K Clausius-Clapeyron, 14.9 K indirect Maxwell, 272.2/
278.2 K zero-field Tc) were independently re-read from the primary-source
PDF this phase and match exactly.
"""

import numpy as np
from core.first_order_mce import GD5SI2GE2_FIRST_ORDER

mu0 = 4 * np.pi * 1e-7

# --- Giguere et al. (1999) reported values, read directly from the PDF ---
GIGUERE_DIRECT_DTAD_7T = 10.0     # K, Fig. 3 direct measurement at 7 T
GIGUERE_CLAUSIUS_CLAPEYRON_DTAD_7T = 9.9   # K, independent latent-heat cross-check
GIGUERE_INDIRECT_MAXWELL_DTAD_7T = 14.9    # K, their own Maxwell-relation calc (~= Ref. [5]'s "giant" value)
GIGUERE_ZERO_FIELD_TC_RANGE = (272.2, 278.2)  # K, increasing/decreasing-field extrapolation


def _model_peak_dTad(mu0H_tesla, T_range=(260.0, 300.0), n=801):
    """Scans T for the model's own peak DeltaT_ad at a given field, same
    approach as first_order_mce.py's __main__ block (the peak is not at
    the nominal Tc since the transition shifts with field)."""
    H = mu0H_tesla / mu0
    Ts = np.linspace(*T_range, n)
    dT = GD5SI2GE2_FIRST_ORDER.delta_T_adiabatic(Ts, H)
    i = int(np.argmax(dT))
    return float(Ts[i]), float(dT[i])


def run_validation(out_path="results/giguere_validation.txt", verbose=True):
    lines = []
    lines.append("First-order Landau model (first_order_mce.py) vs. Giguere et al. (1999)")
    lines.append("DIRECT adiabatic temperature change measurement, Gd5Si2Ge2, 7 T")
    lines.append("=" * 78)

    T_peak, dT_model = _model_peak_dTad(7.0)
    err_vs_direct = 100 * (dT_model - GIGUERE_DIRECT_DTAD_7T) / GIGUERE_DIRECT_DTAD_7T
    err_vs_indirect = 100 * (dT_model - GIGUERE_INDIRECT_MAXWELL_DTAD_7T) / GIGUERE_INDIRECT_MAXWELL_DTAD_7T
    overestimate_factor_direct = dT_model / GIGUERE_DIRECT_DTAD_7T
    overestimate_factor_indirect = dT_model / GIGUERE_INDIRECT_MAXWELL_DTAD_7T
    papers_own_overestimate_factor = GIGUERE_INDIRECT_MAXWELL_DTAD_7T / GIGUERE_DIRECT_DTAD_7T

    lines.append(f"Model peak DeltaT_ad at 7 T: {dT_model:.2f} K (at T={T_peak:.1f} K)")
    lines.append(f"Giguere et al. DIRECT measurement at 7 T: {GIGUERE_DIRECT_DTAD_7T:.1f} K "
                 f"(Clausius-Clapeyron cross-check: {GIGUERE_CLAUSIUS_CLAPEYRON_DTAD_7T:.1f} K)")
    lines.append(f"Giguere et al. INDIRECT (Maxwell-relation, 'giant') value at 7 T: "
                 f"{GIGUERE_INDIRECT_MAXWELL_DTAD_7T:.1f} K")
    lines.append("")
    lines.append(f"Model error vs. DIRECT measurement:   {err_vs_direct:+.0f}%  "
                 f"(model overestimates by {overestimate_factor_direct:.2f}x)")
    lines.append(f"Model error vs. INDIRECT (Maxwell):    {err_vs_indirect:+.0f}%  "
                 f"(model overestimates by {overestimate_factor_indirect:.2f}x)")
    lines.append(f"For reference, Giguere et al.'s OWN indirect-vs-direct overestimate: "
                 f"{papers_own_overestimate_factor:.2f}x ({GIGUERE_INDIRECT_MAXWELL_DTAD_7T:.1f} K / "
                 f"{GIGUERE_DIRECT_DTAD_7T:.1f} K)")
    lines.append("")
    lines.append("CONCLUSION: the model, calibrated to the widely-quoted 'giant' peak DeltaS_M")
    lines.append("(~18 J/(kg K) at 5T, itself a Maxwell-relation/indirect literature value),")
    lines.append("overestimates Giguere et al.'s DIRECTLY measured DeltaT_ad by a wider margin")
    lines.append(f"({overestimate_factor_direct:.2f}x) than the paper's own indirect-vs-direct gap")
    lines.append(f"({papers_own_overestimate_factor:.2f}x). This is NOT a contradiction of the paper --")
    lines.append("it is additive: Giguere et al.'s finding is that Maxwell-relation DeltaS_M")
    lines.append("overestimates DeltaT_ad for a first-order transition; honesty flag #1 in")
    lines.append("first_order_mce.py separately notes this model's lattice-only C_lattice")
    lines.append("denominator (no latent-heat correction) is inappropriate for a first-order")
    lines.append("transition. Both effects point the same way and appear to compound here.")
    lines.append("")
    lines.append("This model should NOT be refit to match Giguere's direct DeltaT_ad instead --")
    lines.append("doing so would abandon its documented calibration to the peak DeltaS_M literature")
    lines.append("value, and this 0-D lattice-only-C_p framework cannot simultaneously match both")
    lines.append("(that would require resolving the transition's latent heat, which honesty flag #1")
    lines.append("already flags as out of scope for this model). Instead: DeltaT_ad predictions from")
    lines.append("this module should be treated as upper-bound-ish, roughly 2-2.5x optimistic vs. a")
    lines.append("direct measurement at the one field/composition point checked here, and any")
    lines.append("downstream design conclusion (e.g. cascade.py's graded-cascade capacity/COP numbers)")
    lines.append("that depends on this module's DeltaT_ad should be read with that correction in mind")
    lines.append("(see core.first_order_mce.composition_tuned_material's apply_giguere_correction flag).")

    text = "\n".join(lines) + "\n"
    with open(out_path, "w") as f:
        f.write(text)
    if verbose:
        print(text)
        print(f"Wrote {out_path}")
    return {
        "T_peak_K": T_peak,
        "model_dTad_7T_K": dT_model,
        "direct_dTad_7T_K": GIGUERE_DIRECT_DTAD_7T,
        "indirect_dTad_7T_K": GIGUERE_INDIRECT_MAXWELL_DTAD_7T,
        "overestimate_factor_vs_direct": overestimate_factor_direct,
        "overestimate_factor_vs_indirect": overestimate_factor_indirect,
        "papers_own_overestimate_factor": papers_own_overestimate_factor,
    }


# Correction factor to bring the model's DeltaT_ad (and hence DeltaS-driven
# downstream quantities) in line with Giguere et al.'s DIRECT measurement,
# derived from the single 7 T cross-check above. This is an EMPIRICAL
# correction from one field/composition point, not a physical derivation --
# treat it as a documented, honest fudge factor, not a validated model
# extension. Computed once here (not hardcoded blind) so it stays traceable
# to run_validation()'s numbers if the calibration in first_order_mce.py
# ever changes.
_T_peak, _dT_model_7T = _model_peak_dTad(7.0)
DTAD_CORRECTION_FACTOR = GIGUERE_DIRECT_DTAD_7T / _dT_model_7T  # ~0.41-0.42


# --- Paper-Mining Pass Part 3, §2: Pecharsky & Gschneidner (1997)'s own
#     Gd5Si2Ge2-vs-Gd peak DeltaT_ad ratio, a SECOND independent primary
#     source (heat-capacity-based, not pulse-field-thermometry-based like
#     Giguere et al. above), and a SECOND independent field point (2T/5T,
#     not the 7T DTAD_CORRECTION_FACTOR was fit to) ---
#
# Source: Pecharsky & Gschneidner, "Giant Magnetocaloric Effect in
# Gd5Si2Ge2", Phys. Rev. Lett. 78, 4494 (1997) -- already cited in this
# codebase for GD5SI2GE2_FIRST_ORDER's Tc/J, but this specific number
# (read directly from the PDF, not assumed) had not previously been
# pulled from its own text: "The DeltaT_ad values of Gd5Si2Ge2 are larger
# than the corresponding DeltaT_ad values for Gd by about 30%, comparing
# the peak values, regardless of the temperature." This is a
# HEAT-CAPACITY-DERIVED comparison (S(T) integrated from C(T,H) for BOTH
# materials, same paper, same Fig. 6, field changes of 0->2T and 0->5T),
# independent of both the direct-measurement method (Giguere et al. 1999)
# and the Maxwell-relation/DeltaS_M route this codebase's Landau model is
# itself calibrated to.
PECHARSKY_1997_PEAK_RATIO = 1.30  # Gd5Si2Ge2 peak DeltaT_ad / Gd peak DeltaT_ad, ~field-independent per the paper's own "regardless of the temperature" framing


def run_pecharsky_ratio_check(verbose=True):
    """Checks the model's own peak-DeltaT_ad(Gd5Si2Ge2)/peak-DeltaT_ad(Gd)
    ratio against Pecharsky & Gschneidner (1997)'s ~1.30 figure, at BOTH
    fields their own Fig. 6 comparison uses (2T and 5T) -- for BOTH the
    raw (uncorrected) model and the Giguere-et-al.-derived
    DTAD_CORRECTION_FACTOR-corrected model, since it's not obvious in
    advance which one this independent check should land closer to.

    ACTUAL RESULT (do not silently update this docstring to hide an
    unfavorable finding -- same precedent run_validation()'s own docstring
    above and core.validation.run_curie_shift_check() set): at 5T, the
    RAW (uncorrected) model's ratio comes out at ~1.24 -- close to
    Pecharsky & Gschneidner's ~1.30, an unexpected agreement given it was
    never fit to this number. But applying DTAD_CORRECTION_FACTOR (fit to
    a single 7T DIRECT measurement) drags the SAME ratio down to ~0.51 at
    5T (and ~0.87 at 2T) -- i.e. the corrected model predicts Gd5Si2Ge2
    UNDERPERFORMS pure Gd, which contradicts the basic "giant" MCE premise
    the whole first_order_mce.py module exists to represent. This is a
    genuine, informative limitation of DTAD_CORRECTION_FACTOR, not a
    contradiction dismissed as noise: a single-field (7T) empirical
    correction, when extrapolated down to 2-5T, OVERCORRECTS. The
    correction's own docstring already says "treat as an honest fudge
    factor from a single field/composition point, not a validated model
    extension" -- this check is the concrete evidence for exactly that
    caveat, from an independent source and field range. NOT used to
    change DTAD_CORRECTION_FACTOR's value or composition_tuned_material()'s
    default (see that function's own docstring in first_order_mce.py for
    why re-fitting again would just repeat the same single-point-
    calibration problem one level up).
    """
    from core.mce_material import GADOLINIUM
    rows = []
    for B in (2.0, 5.0):
        T_peak_gd5, dT_gd5_raw = _model_peak_dTad(B)
        dT_gd5_corrected = dT_gd5_raw * DTAD_CORRECTION_FACTOR
        Ts = np.linspace(270.0, 320.0, 1001)
        dT_gd_curve = GADOLINIUM.delta_T_adiabatic(Ts, B / mu0)
        i = int(np.argmax(dT_gd_curve))
        dT_gd_peak = float(dT_gd_curve[i])
        ratio_raw = dT_gd5_raw / dT_gd_peak
        ratio_corrected = dT_gd5_corrected / dT_gd_peak
        rows.append({
            "field_T": B, "gd5si2ge2_peak_raw_K": round(dT_gd5_raw, 2),
            "gd5si2ge2_peak_corrected_K": round(dT_gd5_corrected, 2),
            "gd_peak_K": round(dT_gd_peak, 2),
            "ratio_raw": round(ratio_raw, 2), "ratio_corrected": round(ratio_corrected, 2),
        })
        if verbose:
            print(f"{B:.0f}T: Gd5Si2Ge2 peak dTad raw={dT_gd5_raw:5.2f}K "
                  f"corrected={dT_gd5_corrected:5.2f}K  |  Gd peak dTad={dT_gd_peak:5.2f}K  |  "
                  f"ratio raw={ratio_raw:.2f}  corrected={ratio_corrected:.2f}  "
                  f"(Pecharsky & Gschneidner 1997: ~{PECHARSKY_1997_PEAK_RATIO:.2f})")
    if verbose:
        print("Finding: the RAW model's ratio is closer to Pecharsky & Gschneidner's "
              "~1.30 than the Giguere-corrected model's ratio is -- applying the single-"
              "field (7T) DTAD_CORRECTION_FACTOR at 2T/5T overcorrects to the point of "
              "predicting Gd5Si2Ge2 underperforms Gd (ratio < 1), contradicting the "
              "'giant' MCE premise. See this function's docstring.")
    return rows


def run_latent_heat_validation(verbose=True):
    """checks GD5SI2GE2_FIRST_ORDER_LATENT_HEAT (the
    field-tracked-latent-heat-Cp-spike variant, see that instance's own
    block comment in first_order_mce.py for full derivation/citations)
    against BOTH cross-checks already in this module -- Giguere et al.'s
    7T direct DeltaT_ad, and Pecharsky & Gschneidner's ~1.30 Gd5Si2Ge2/Gd
    peak-ratio at 2T/5T -- alongside the existing GD5SI2GE2_FIRST_ORDER
    (no latent heat) and DTAD_CORRECTION_FACTOR-corrected numbers, so all
    three treatments are visible side by side.

    ACTUAL RESULT (same "report the real finding" precedent as this
    module's other two check functions above):
      - 7T vs. Giguere's direct 10.0K: RAW model 24.17K -> LATENT-HEAT
        model 19.01K -- a genuine ~36% closure of the gap using only
        literature-grounded (L, sigma), NOT a full fix, and NOT tuned to
        hit 10.0K (see GD5SI2GE2_FIRST_ORDER_LATENT_HEAT's own comment for
        why sigma was not re-tuned further).
      - 2T/5T vs. Pecharsky & Gschneidner's ~1.30 Gd5Si2Ge2/Gd ratio: the
        latent-heat model gives 0.70 (2T) / 0.91 (5T) / 0.81 (7T) -- closer
        to 1.0 (and so less severely wrong) than the flat
        DTAD_CORRECTION_FACTOR's 0.87 (2T) / 0.51 (5T), but STILL below
        1.0 at every field checked, i.e. it ALSO predicts Gd5Si2Ge2
        underperforms plain Gd, which still contradicts the "giant" MCE
        premise Pecharsky & Gschneidner's own ~1.30 finding represents.
        This is a real, partial improvement over the flat correction
        factor (less wrong, not right), not a resolution of
        run_pecharsky_ratio_check()'s own documented finding.

    CONCLUSION: latent heat is a genuine, literature-grounded physical
    mechanism this model was missing, and adding it (correctly, tracking
    the field-dependent transition location rather than a fixed Tc) closes
    real ground on BOTH independent checks in this module without being
    tuned to either one. It does not fully resolve either gap. Do not
    read GD5SI2GE2_FIRST_ORDER_LATENT_HEAT as a validated replacement for
    GD5SI2GE2_FIRST_ORDER or for the DTAD_CORRECTION_FACTOR path --
    downstream code should keep using whichever of those two the calling
    context already used (see ROADMAP.md for the decision not to
    swap either default in this pass)."""
    from core.mce_material import GADOLINIUM
    from core.first_order_mce import GD5SI2GE2_FIRST_ORDER, GD5SI2GE2_FIRST_ORDER_LATENT_HEAT

    lines = []
    lines.append("GD5SI2GE2_FIRST_ORDER_LATENT_HEAT vs. GD5SI2GE2_FIRST_ORDER (no latent heat)")
    lines.append("=" * 78)

    Ts_wide = np.linspace(260.0, 300.0, 1601)
    dT_raw = GD5SI2GE2_FIRST_ORDER.delta_T_adiabatic(Ts_wide, 7.0 / mu0)
    dT_lh = GD5SI2GE2_FIRST_ORDER_LATENT_HEAT.delta_T_adiabatic(Ts_wide, 7.0 / mu0)
    peak_raw = float(np.max(dT_raw))
    peak_lh = float(np.max(dT_lh))
    gap_closed_pct = 100 * (peak_raw - peak_lh) / (peak_raw - GIGUERE_DIRECT_DTAD_7T)

    lines.append(f"7T peak DeltaT_ad: raw={peak_raw:.2f}K latent-heat={peak_lh:.2f}K "
                 f"Giguere direct target={GIGUERE_DIRECT_DTAD_7T:.1f}K")
    lines.append(f"Gap to direct target closed by latent heat: {gap_closed_pct:.1f}% "
                 f"(NOT a full fix)")
    lines.append("")

    rows = []
    for B in (2.0, 5.0, 7.0):
        Ts = np.linspace(260.0, 320.0, 2001)
        dT_gd5_lh = GD5SI2GE2_FIRST_ORDER_LATENT_HEAT.delta_T_adiabatic(Ts, B / mu0)
        dT_gd = GADOLINIUM.delta_T_adiabatic(Ts, B / mu0)
        peak_gd5_lh = float(np.max(dT_gd5_lh))
        peak_gd = float(np.max(dT_gd))
        ratio_lh = peak_gd5_lh / peak_gd
        rows.append({"field_T": B, "gd5si2ge2_latent_heat_peak_K": round(peak_gd5_lh, 2),
                     "gd_peak_K": round(peak_gd, 2), "ratio_latent_heat": round(ratio_lh, 2)})
        lines.append(f"{B:.0f}T: Gd5Si2Ge2(latent-heat) peak={peak_gd5_lh:5.2f}K "
                     f"Gd peak={peak_gd:5.2f}K ratio={ratio_lh:.2f}  "
                     f"(Pecharsky & Gschneidner 1997: ~{PECHARSKY_1997_PEAK_RATIO:.2f}; "
                     f"still <1.0 -- see this function's docstring)")

    text = "\n".join(lines)
    if verbose:
        print(text)
    return {
        "peak_dTad_7T_raw_K": peak_raw,
        "peak_dTad_7T_latent_heat_K": peak_lh,
        "gap_closed_pct": gap_closed_pct,
        "pecharsky_ratio_rows": rows,
    }


if __name__ == "__main__":
    run_validation()
    print("\n" + "=" * 78)
    print("Extension: Gd5Si2Ge2/Gd peak DeltaT_ad ratio vs. Pecharsky & Gschneidner "
          "(1997) (Paper-Mining Pass Part 3, S2)")
    print("=" * 78)
    run_pecharsky_ratio_check()
    print("\n" + "=" * 78)
    print("latent-heat C_p spike (GD5SI2GE2_FIRST_ORDER_LATENT_HEAT)")
    print("=" * 78)
    run_latent_heat_validation()