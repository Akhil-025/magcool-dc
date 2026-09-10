"""
passive_regenerator_analysis.py
================================
 addition: exercises core/baseline_cooling.py's
`augmented_regenerator_cop()` / `passive_regenerator_augmentation()` (the
new functions added for this phase) across this repo's existing material
library (core/mce_material.py) at the representative ASHRAE data-center
operating point already used throughout this repo (T_cold=291.15K,
span=10K -- see main.py's REPRESENTATIVE_SPAN_K / step 4).

Scope and honesty flag
-----------------------
UPDATE: corrected after direct re-check. See core/baseline_cooling.py's
own docstring block for the full book-access honesty flag. Tishin &
Spichkin (2003) has no embedded text layer (pdfplumber extraction still
returns 0 characters -- that part was correct), but OCR (tesseract) has
now actually been run across Sect. 11.1 (passive magnetic regenerators,
book pp.359-374: rare-earth intermetallics and rare-earth-metal
regenerator materials, Table 11.1 peak heat-capacity data) and Sect.
11.2's general-consideration equations (pp.374-390: Carnot/AMR-cycle
entropy balance eq. 11.2-11.20, N_tu-effectiveness relations, and
Barclay & Sarangi's (1984) regenerator-geometry design guidance -- see
regenerator_specific_area_design_point() below for what that guidance
actually was and how it's used here).

A specific earlier claim is corrected here directly: a prior pass
asserted this OCR would recover "a real penetration-depth formula
(eq. 11.1)" for the passive-regenerator section. That formula does not
exist in this book -- checked directly by OCR-ing and text-searching
the entirety of chapter 11 (pp.359-418, both the passive-regenerator
and active-refrigeration sections) for "penetration depth", "diffusion
depth" and "thermal wave" with zero matches. Equation (11.1) itself is
an unlabelled Carnot-cycle heat-absorption relation (Q_c = T_cold *
Delta S_m), not a penetration-depth formula. No penetration-depth
function has been added to this repo as a result -- adding one would
have meant fabricating it. What the OCR pass DID recover and IS usable
is the real regenerator design guidance now in
regenerator_specific_area_design_point(): the porosity/particle-
diameter/specific-area/pressure-drop trade-off Barclay & Sarangi (1984)
worked through for packed-bed, tube, and plate regenerator geometries.

This module remains a design-exploration / comparison tool, not a
validated benchmark-backed result -- the same disposition
core/thermal_diode_analysis.py (thermal diodes) and
core/fluid_mce_analysis.py (magnetocaloric fluids) gave their own new
modules, for the same underlying reason (no benchmark device for this
specific configuration exists in data/amr_experimental_benchmarks.csv,
which is solid-AMR-only).

Only core/mce_material.py's second-order (mean-field/Brillouin)
MagnetocaloricMaterial instances are used -- GADOLINIUM, GD5SI2GE2 (its
mean-field parameterization, retained per that module's own "parameter
library entry only" caveat), and LACAMNO3 -- because only that class
exposes total_heat_capacity()'s lambda-anomaly term. First-order
materials (core/first_order_mce.py's FirstOrderMCEMaterial) do not model
a Curie-point heat-capacity peak at all (their giant MCE comes from a
latent-heat-like entropy discontinuity instead), so they are not
meaningful "passive regenerator" candidates under this specific
mechanism and are deliberately excluded rather than silently coerced.
"""

from dataclasses import dataclass
from core.mce_material import GADOLINIUM, GD5SI2GE2, LACAMNO3
from core.baseline_cooling import (
    vapor_compression_cop, augmented_regenerator_cop,
    MAX_COP_GAIN_AT_FULL_EFFECTIVENESS,
)

T_COLD_K = 291.15   # 18 C, matches main.py's REPRESENTATIVE operating point
SPAN_K = 10.0
T_HOT_K = T_COLD_K + SPAN_K

CANDIDATE_MATERIALS = (GADOLINIUM, GD5SI2GE2, LACAMNO3)


def compare_candidate_materials(T_cold=T_COLD_K, T_hot=T_HOT_K, verbose=True):
    """Runs augmented_regenerator_cop() for every material in
    CANDIDATE_MATERIALS against the same base vapor-compression COP, so the
    only thing that differs row-to-row is how well each material's own
    Curie temperature aligns with [T_cold, T_hot]. Returns (base_cop,
    rows) where rows is a list of PassiveRegeneratorResult, sorted by
    descending augmented_COP."""
    base = vapor_compression_cop(T_cold, T_hot)
    results = [augmented_regenerator_cop(base.COP, m, (T_cold, T_hot))
               for m in CANDIDATE_MATERIALS]
    results.sort(key=lambda r: r.augmented_COP, reverse=True)
    if verbose:
        print(f"Base vapor-compression COP at Tc={T_cold:.2f}K, Th={T_hot:.2f}K "
              f"(span={T_hot - T_cold:.1f}K): {base.COP:.4f} "
              f"(Carnot={base.COP_carnot:.4f}, eta_2nd_law={base.second_law_eff:.3f})")
        for r in results:
            print(f"  {r.material_name:32s}  Tc_material={_material_tc(r):6.1f}K "
                  f"eps: {r.eps_baseline:.3f} -> {r.eps_augmented:.3f}  "
                  f"(delta={r.delta_eps:+.3f})   COP: {r.base_COP:.4f} -> "
                  f"{r.augmented_COP:.4f}  ({r.cop_gain_fraction:+.2%})")
    return base, results


def _material_tc(result):
    for m in CANDIDATE_MATERIALS:
        if m.name == result.material_name:
            return m.Tc
    return float("nan")


def span_sweep(spans_K=(5.0, 10.0, 15.0, 20.0), T_cold=T_COLD_K, verbose=True):
    """Repeats compare_candidate_materials() across a few spans at the same
    T_cold, to show the alignment effect directly: as span widens, the
    [T_cold, T_hot] window increasingly overlaps a fixed material's Curie
    temperature (or stops overlapping it), so delta_eps and the resulting
    COP gain should move accordingly rather than staying flat -- checked
    here, not assumed."""
    rows = []
    for span in spans_K:
        base, results = compare_candidate_materials(T_cold, T_cold + span, verbose=False)
        best = results[0]
        rows.append({
            "span_K": span, "base_COP": base.COP, "best_material": best.material_name,
            "best_augmented_COP": best.augmented_COP, "best_delta_eps": best.delta_eps,
            "best_cop_gain_fraction": best.cop_gain_fraction,
        })
        if verbose:
            print(f" span={span:5.1f}K base_COP={base.COP:.4f}  "
                  f"best={best.material_name:32s}  delta_eps={best.delta_eps:+.3f}  "
                  f"augmented_COP={best.augmented_COP:.4f} ({best.cop_gain_fraction:+.2%})")
    return rows


def run_passive_regenerator_analysis(out_path="results/passive_regenerator_analysis.txt"):
    import io, contextlib
    buf = io.StringIO()
    with contextlib.redirect_stdout(buf):
        print("=" * 90)
        print("PHASE 21: passive/hybrid magnetic regenerator augmentation of a")
        print("conventional (vapor-compression) gas cycle -- see")
        print("core/baseline_cooling.py's own docstring block for the honesty")
        print("flag (Tishin Ch.11 OCR'd -- passive-regenerator materials and design")
        print("guidance recovered; no digitized effectiveness/COP curve for THIS")
        print("module's VCC-augmentation use case, and no penetration-depth formula")
        print("exists in the source -- see this module's own docstring)")
        print("and for why the effectiveness-to-COP mapping is an illustrative,")
        print("literature-range-anchored ceiling rather than a fitted coefficient.")
        print("=" * 90)

        print(f"\n--- Candidate-material comparison at the representative ASHRAE point "
              f"(T_cold={T_COLD_K:.2f}K, span={SPAN_K}K) ---")
        base, results = compare_candidate_materials()

        print(f"\n--- Span sweep at fixed T_cold={T_COLD_K:.2f}K "
              "(does the alignment effect move as expected as the window widens?) ---")
        sweep_rows = span_sweep()

        best = results[0]
        worst = results[-1]
        print("\n--- Conclusion ---")
        print(f"At the representative operating point, {best.material_name} (its own "
              f"Curie temperature sits inside [{T_COLD_K:.1f}, {T_HOT_K:.1f}]K) gives the "
              f"largest passive-regenerator boost: eps {best.eps_baseline:.3f} -> "
              f"{best.eps_augmented:.3f}, COP {best.base_COP:.4f} -> "
              f"{best.augmented_COP:.4f} ({best.cop_gain_fraction:+.2%}), while "
              f"{worst.material_name} (Curie temperature far outside this window) gives "
              f"{worst.cop_gain_fraction:+.2%} -- confirming, in this repo's own model, "
              "the plan's own framing that the benefit is an alignment effect, not a "
              "fixed per-material bonus. Every reported gain is capped by construction at "
              f"{MAX_COP_GAIN_AT_FULL_EFFECTIVENESS:.0%} (the illustrative full-effectiveness "
              "ceiling -- see the honesty flag above), so these numbers should be read as "
              "'this mechanism could plausibly be worth up to X%, IF the underlying "
              "literature-range ceiling holds for a magnetically-augmented regenerator "
              "specifically, which has not been separately confirmed' -- not as a "
              "validated device-level COP prediction.")

    text = buf.getvalue()
    print(text, end="")
    if out_path:
        import os
        os.makedirs(os.path.dirname(out_path), exist_ok=True)
        with open(out_path, "w") as f:
            f.write(text)
    return {"base": base, "candidate_results": results, "span_sweep": sweep_rows}


# ---------------------------------------------------------------------------
# Regenerator specific-area / particle-diameter design guidance -- Barclay &
# Sarangi (1984), as reported in Tishin & Spichkin (2003) Sect. 11.2.2
# (book p.373), recovered by the OCR pass described in this module's own
# docstring above. This is real, citable design content (NOT a penetration-
# depth formula -- see the docstring's correction of that earlier claim),
# additive and independently testable, same pattern as the tariff work in
# core/economics.py.
# ---------------------------------------------------------------------------

# Typical overall porosity (void fraction) for a randomly packed bed of
# spherical regenerator particles, per Barclay & Sarangi (1984) as quoted
# in Tishin & Spichkin (2003) p.373 ("The overall porosity ... was chosen
# to be 0.4 for all cases under consideration (this value is typical for
# a randomly packed bed of spherical particles)").
TYPICAL_PACKED_BED_POROSITY = 0.4

# Barclay (1991)'s own reported regenerator design point, as quoted in
# Tishin & Spichkin (2003) p.373: N_tu=500 (a high-effectiveness target)
# at an operational frequency of 0.5 Hz requires a specific area of
# 25000 m^-1, achievable with particles of 240 um diameter. Used below
# purely as a sanity-check reference point, not as an input to any
# calculation.
_BARCLAY_1991_REFERENCE_POINT = {
    "N_tu": 500, "frequency_Hz": 0.5,
    "specific_area_per_m": 25000.0, "particle_diameter_m": 240e-6,
}


def regenerator_specific_area_per_m(particle_diameter_m):
    """Surface-to-volume ratio (specific area, m^-1) of a spherical
    regenerator particle: a = 6 / d. This is the elementary sphere
    geometry relation, not itself a Tishin-specific formula -- but it is
    exactly the relation that reproduces Tishin & Spichkin (2003)'s own
    quoted Barclay (1991) design point (240 um -> 25000 m^-1, checked in
    check_against_barclay_1991_reference_point() below), which is why it
    is the "specific area" this repo uses when citing that design
    guidance, rather than some other packed-bed correlation that would
    NOT reproduce the book's own reported number."""
    if particle_diameter_m <= 0:
        raise ValueError("particle_diameter_m must be positive")
    return 6.0 / particle_diameter_m


def particle_diameter_for_specific_area(specific_area_per_m):
    """Inverse of regenerator_specific_area_per_m(): the spherical
    particle diameter (m) needed to hit a target specific area (m^-1)."""
    if specific_area_per_m <= 0:
        raise ValueError("specific_area_per_m must be positive")
    return 6.0 / specific_area_per_m


def check_against_barclay_1991_reference_point(tol_frac=0.01):
    """Verifies regenerator_specific_area_per_m() reproduces Tishin &
    Spichkin (2003)'s own quoted Barclay (1991) design point (240 um
    particles -> 25000 m^-1 specific area, at N_tu=500 / 0.5 Hz) to
    within tol_frac. Returns (matches: bool, computed_specific_area).
    This is a literature cross-check, not a fitted coefficient -- if it
    ever fails, that means the OCR'd reference number or this function's
    formula disagree and BOTH should be re-examined, not silently
    reconciled."""
    ref = _BARCLAY_1991_REFERENCE_POINT
    computed = regenerator_specific_area_per_m(ref["particle_diameter_m"])
    matches = abs(computed - ref["specific_area_per_m"]) <= tol_frac * ref["specific_area_per_m"]
    return matches, computed


@dataclass
class RegeneratorDesignPoint:
    particle_diameter_m: float
    specific_area_per_m: float
    porosity: float
    note: str


def regenerator_specific_area_design_point(particle_diameter_m,
                                            porosity=TYPICAL_PACKED_BED_POROSITY):
    """Returns a RegeneratorDesignPoint for a packed bed of spherical
    regenerator particles of the given diameter, at the given porosity
    (defaulting to Barclay & Sarangi's (1984) typical value of 0.4, per
    Tishin & Spichkin (2003) p.373). Reports specific area via
    regenerator_specific_area_per_m(); porosity is carried through for
    reference/reporting only -- it does NOT enter that formula (the
    book's own reference point is matched by the bare 6/d sphere
    relation, not a porosity-scaled packed-bed correlation -- see that
    function's own docstring), since decreasing particle diameter
    increases specific area but also increases pressure drop across the
    regenerator (requiring additional fluid-pump work), a qualitative
    trade-off Tishin & Spichkin (2003) p.373 states explicitly but does
    not give a quantitative pressure-drop correlation for -- NOT modeled
    quantitatively here, stated as an open gap rather than silently
    ignored."""
    a = regenerator_specific_area_per_m(particle_diameter_m)
    return RegeneratorDesignPoint(
        particle_diameter_m=particle_diameter_m,
        specific_area_per_m=a,
        porosity=porosity,
        note=("Decreasing particle diameter increases specific area but "
              "also increases regenerator pressure drop (additional pump "
              "work) -- Tishin & Spichkin (2003) p.373, Barclay & Sarangi "
              "(1984) -- not quantitatively modeled here."),
    )


if __name__ == "__main__":
    run_passive_regenerator_analysis()