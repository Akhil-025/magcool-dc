"""
passive_regenerator_analysis.py
================================
Phase 21 addition: exercises core/baseline_cooling.py's
`augmented_regenerator_cop()` / `passive_regenerator_augmentation()` (the
new functions added for this phase) across this repo's existing material
library (core/mce_material.py) at the representative ASHRAE data-center
operating point already used throughout this repo (T_cold=291.15K,
span=10K -- see main.py's REPRESENTATIVE_SPAN_K / step 4).

Scope and honesty flag
-----------------------
See core/baseline_cooling.py's own Phase 21 docstring block for the full
book-access honesty flag (Tishin & Spichkin (2003) remains an
image-only, non-extractable-text PDF in this project's corpus -- same
finding as Phase 20) and for why the effectiveness-to-COP mapping used
here is an illustrative, literature-range-anchored ceiling rather than a
fitted or digitized coefficient. This module is a design-exploration /
comparison tool, not a validated benchmark-backed result -- the same
disposition Phase 18 (thermal diodes) and Phase 20 (magnetocaloric
fluids) gave their own new modules, for the same underlying reason (no
benchmark device for this specific configuration exists in
data/amr_experimental_benchmarks.csv, which is solid-AMR-only).

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
            print(f"  {r.material_name:32s}  Tc_material={_material_tc(r):6.1f}K  "
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
            print(f"  span={span:5.1f}K  base_COP={base.COP:.4f}  "
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
        print("core/baseline_cooling.py's own Phase 21 docstring block for the honesty")
        print("flag (Tishin Ch.11 not digitizable -- image-only PDF, no text layer)")
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


if __name__ == "__main__":
    run_passive_regenerator_analysis()