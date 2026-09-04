"""
nanocomposite_material.py
===========================
 (see phase_plan.md and ROADMAP.md): adds a genuinely new
CANDIDATE MATERIAL FAMILY -- an engineered multi-phase nanocomposite blend
-- alongside Gd / Gd5Si2Ge2 / La(Fe,Si)13Hy / (Mn,Fe)2(P,Si) in
core/material_family_comparison.py.

HONESTY FLAG -- book access and scope (checked directly, not assumed)
-----------------------------------------------------------------------
Tishin & Spichkin (2003) Sect. 2.9 and Ch. 10 (superparamagnetic /
nanocomposite / molecular-cluster systems) are this item's own named
sources per phase_plan.md. Re-confirmed directly for this pass: pdfplumber
extracts zero characters from every page of this project's copy sampled
(0, 1, 2, 50, 51) -- the same image-only-PDF finding already recorded for
Tishin Ch. 11 and Sect. 2.8 . Those sections'
specific content -- whatever superparamagnetic/nanocomposite model or
material data they report -- could not be read or digitized here.

phase_plan.md itself flags this item as "lower priority -- needs new
composition-tuning data your corpus doesn't currently have digitized."
This module honors that limitation rather than working around it: it does
NOT invent new composition/Tc/DeltaS_M data for a superparamagnetic
nanoparticle system (no such data exists anywhere in this repo's corpus,
and none could be read from the book). It also does NOT attempt true
superparamagnetic physics (Neel relaxation, a blocking-temperature-
dependent reduction of the effective moment, particle-size-distribution
effects on M(H,T)) -- that would need particle-size-distribution and
anisotropy-constant data this project does not have.

What IS implemented, and why it is a legitimate, narrower reading of
"nanocomposite": the broader MCE engineering literature (e.g. Smith et
al., "Materials challenges for high performance magnetocaloric
refrigeration devices", Adv. Energy Mater. 2 (2012) 1288-1318; Bahl et
al.'s reviews of layered/graded AMR beds) documents ENGINEERED composite
regenerators -- deliberately blending multiple composition-tuned phases
of the SAME base material family, each with a slightly different Tc, to
broaden/flatten the working temperature range of a single regenerator
layer -- as a real, established design strategy, distinct from
core/cascade.py's discrete-STAGE Curie grading (which places a single
sharp-Tc material per physically separate stage). This module implements
that specific, narrower reading: a `WeightedMaterialEnsemble` of 3
composition-tuned La(Fe,Si)13Hy phases (core/cascade.py's own already-
validated, already-tunable LAFESIH_FAMILY -- chosen because
material_family_comparison.py's own Track A2 ranking already found
it the best-performing tunable family at the ASHRAE point), spread
+/- NANOCOMPOSITE_SPREAD_K around a target composition, blended with a
fixed triangular weighting -- rather than any new, undigitized material
data.

This deliberately reuses core/inhomogeneous_broadening.py's core insight -- multiple Tc-shifted phases combined into one
effective response -- but for a DIFFERENT physical situation (a few
discrete, deliberately-engineered phases spread over a design-chosen
range, vs. item 1's many random grains distributed by manufacturing
inhomogeneity) and, necessarily, a DIFFERENT mixing rule (see
WeightedMaterialEnsemble's own docstring for why).
"""

import dataclasses

import numpy as np

from core.first_order_mce import (
    LAFESIH_FIRST_ORDER, LAFESIH_TC_MIN_K, LAFESIH_TC_MAX_K,
    lafesih_composition_tuned_material,
)
from core.mce_material import GADOLINIUM
from core.cascade import GradedFamily, _target_composition_for_peak

NANOCOMPOSITE_SPREAD_K = 4.0
# Illustrative, design-chosen spread (not a digitized or fitted value --
# see this module's own honesty flag above: no experimental nanocomposite
# curve exists in this repo's corpus to fit against). Kept modest (a few
# K) relative to core/inhomogeneous_broadening.py's swept 0-5K RANDOM
# grain-broadening range, since this represents a deliberately narrow,
# design-chosen spread meant to flatten the working range without
# sacrificing much peak DeltaT_ad -- the opposite intent from item 1's
# broadening (which is an unavoidable manufacturing side-effect, not a
# design choice), even though the underlying math (multiple Tc-shifted
# phases combined) rhymes with it.
NANOCOMPOSITE_WEIGHTS = (0.25, 0.5, 0.25)  # triangular, heaviest at center


class WeightedMaterialEnsemble:
    """A fixed, explicit weighted blend of N (possibly heterogeneous)
    magnetocaloric material objects sharing the common
    `delta_T_adiabatic(T, H_final, H_initial=0.0)` interface used
    throughout this repo (core.mce_material.MagnetocaloricMaterial and
    core.first_order_mce.FirstOrderMCEMaterial both implement it).

    Mixing rule -- deliberately different from
    core.inhomogeneous_broadening.BroadenedMagnetocaloricMaterial's
    entropy/heat-capacity-level averaging, and this difference is stated
    rather than papered over: FirstOrderMCEMaterial instances carry a
    per-phase `dTad_correction` (calibrated against Giguere et al.'s
    DIRECT DeltaT_ad measurement -- see core/giguere_validation.py) that
    is a correction to the WHOLE DeltaT_ad ratio -T*dS/C, not to entropy
    or heat capacity individually. Averaging entropy and heat capacity
    separately (as item 1 does, correctly, for MagnetocaloricMaterial's
    C_lattice-only mean-field model which carries no such correction)
    would silently discard that per-phase correction here. So this
    ensemble instead mixes at the DeltaT_ad level directly:

        DeltaT_ad_ensemble(T) = sum_i w_i * DeltaT_ad_phase_i(T)

    -- a coarser approximation than item 1's (it does not model a shared
    physical sample temperature relaxing across phases with different
    local heat capacities), but the one consistent with how this repo's
    giant-MCE phases are already validated and used everywhere else
    (core/cascade.py, core/material_family_comparison.py both call each
    phase's own delta_T_adiabatic() directly, never its raw entropy/heat
    capacity separately).

    delta_S_isothermal() is still provided (a simple weighted average) as
    a diagnostic/reporting quantity -- it is NOT what feeds
    delta_T_adiabatic() here, unlike in item 1's ensemble.
    """

    def __init__(self, materials, weights, name=None):
        weights = np.asarray(weights, dtype=float)
        if not np.isclose(weights.sum(), 1.0, atol=1e-6):
            raise ValueError(f"weights must sum to 1.0, got {weights.sum():.6f}")
        if len(materials) != len(weights):
            raise ValueError("materials and weights must be the same length")
        self.materials = list(materials)
        self.weights = weights
        self.name = name or (
            f"Nanocomposite blend ({len(materials)} phases, "
            f"Tc={[round(m.Tc, 1) for m in materials]}K)"
        )
        # Weighted-average hysteresis loss, so a nanocomposite built from
        # phases that DO carry hysteresis_loss_J_per_kg values
        # keeps that cost visible to AMRSystem's _hysteresis_power_W()
        # (which uses getattr(..., 0.0), so this attribute is optional but
        # honored automatically if present -- see core/amr_cycle.py).
        self.hysteresis_loss_J_per_kg = float(np.sum(
            weights * np.array([getattr(m, "hysteresis_loss_J_per_kg", 0.0) for m in self.materials])
        ))

    @property
    def Tc(self):
        """Weighted-mean Tc across phases -- informational only; nothing
        in this repo's AMRSystem/cascade machinery reads an ensemble's own
        .Tc directly (see this module's docstring for the interface
        actually required), but GradedFamily.reference_material and
        similar callers elsewhere sometimes inspect .Tc for logging."""
        return float(np.sum(self.weights * np.array([m.Tc for m in self.materials])))

    def delta_S_isothermal(self, T, H_final, H_initial=0.0):
        """Weighted-average entropy change across phases -- a reporting/
        diagnostic quantity only; NOT what delta_T_adiabatic() uses (see
        class docstring)."""
        T = np.atleast_1d(np.asarray(T, dtype=float))
        out = np.zeros_like(T)
        for w, m in zip(self.weights, self.materials):
            out = out + w * np.asarray(m.delta_S_isothermal(T, H_final, H_initial)).ravel()
        return out

    def delta_T_adiabatic(self, T, H_final, H_initial=0.0):
        T = np.atleast_1d(np.asarray(T, dtype=float))
        out = np.zeros_like(T)
        for w, m in zip(self.weights, self.materials):
            out = out + w * np.asarray(m.delta_T_adiabatic(T, H_final, H_initial)).ravel()
        return out


def nanocomposite_tuned_material(Tc_center_K, spread_K=NANOCOMPOSITE_SPREAD_K,
                                  weights=NANOCOMPOSITE_WEIGHTS,
                                  base_tuned_fn=lafesih_composition_tuned_material):
    """Builds a 3-phase engineered nanocomposite blend of LAFESIH-family
    phases at Tc_center-spread_K, Tc_center, Tc_center+spread_K, weighted
    by `weights` (default triangular). Raises ValueError (via
    base_tuned_fn) if any of the three phase Tc's falls outside
    LAFESIH_TC_MIN_K/_MAX_K -- callers should stay within
    NANOCOMPOSITE_TC_MIN_K/_MAX_K below to avoid that."""
    Tc_values = [Tc_center_K - spread_K, Tc_center_K, Tc_center_K + spread_K]
    phases = [base_tuned_fn(tc) for tc in Tc_values]
    return WeightedMaterialEnsemble(
        phases, weights,
        name=f"Nanocomposite (LAFESIH 3-phase blend, center Tc={Tc_center_K:.1f}K, "
             f"spread=+/-{spread_K:.1f}K)")


# The nanocomposite's own outer phases must themselves stay within LAFESIH's
# documented tunability window, so its effective window is tightened by
# spread_K on each side (rather than adding per-phase fallback logic).
NANOCOMPOSITE_TC_MIN_K = LAFESIH_TC_MIN_K + NANOCOMPOSITE_SPREAD_K
NANOCOMPOSITE_TC_MAX_K = LAFESIH_TC_MAX_K - NANOCOMPOSITE_SPREAD_K

NANOCOMPOSITE_FAMILY = GradedFamily(
    name="Nanocomposite (LAFESIH 3-phase blend)",
    tuned_fn=nanocomposite_tuned_material,
    tc_min=NANOCOMPOSITE_TC_MIN_K, tc_max=NANOCOMPOSITE_TC_MAX_K,
    # reference_material seeds _target_composition_for_peak()'s initial
    # peak-vs-Tc offset guess (core/cascade.py) -- a single LAFESIH phase
    # is a reasonable seed since the blend's own peak-vs-center-Tc offset
    # should be close to (though not identical to) a single phase's own
    # offset, and the root-finder self-corrects from there.
    reference_material=LAFESIH_FIRST_ORDER, fallback_material=GADOLINIUM,
)


def run_robustness_check(design_span_K=10.0, off_design_spans_K=(5.0, 10.0, 15.0, 20.0),
                          T_cold_K=291.15, mu0H_max=2.0, mass_per_stage=5.0,
                          out_path="results/nanocomposite_robustness.txt",
                          verbose=True):
    """Point-performance comparisons (material_family_comparison.py) tune
    EVERY candidate freshly at each span's own T_mid -- the fairest-shot
    comparison for "which composition should I pick for THIS known span",
    but it cannot show whether a blend's breadth is worth anything, since
    a perfectly retuned single phase always wins at its own design point
    by construction (blending trades peak height for width -- see this
    module's own docstring). The actual value proposition the broader MCE
    literature (Smith et al. 2012 and similar) attributes to composite/
    graded regenerators is ROBUSTNESS: performance when the built device's
    fixed composition ends up operating away from its design point (e.g.
    real ambient conditions drifting from the design span).

    This function checks that directly: builds ONE nanocomposite blend and
    ONE single LAFESIH phase, BOTH tuned once at design_span_K's own
    T_mid, then evaluates both (WITHOUT retuning) at every span in
    off_design_spans_K via core.cascade.run_cascade's 1-stage case."""
    from core.cascade import run_cascade

    T_mid_design = T_cold_K + design_span_K / 2.0
    Tc_design = _target_composition_for_peak(T_mid_design, mu0H_max, NANOCOMPOSITE_FAMILY)
    nanocomposite = nanocomposite_tuned_material(Tc_design)
    single_phase = lafesih_composition_tuned_material(Tc_design)

    lines = []

    def log(s=""):
        if verbose:
            print(s)
        lines.append(s)

    log("=" * 90)
    log("PHASE 22 ITEM 2 (follow-up): does the nanocomposite blend's breadth pay off")
    log("when the FIXED, already-built composition operates away from its design span?")
    log("=" * 90)
    log(f"Both candidates tuned ONCE at design_span_K={design_span_K:.1f}K's own "
        f"T_mid={T_mid_design:.2f}K (Tc_design={Tc_design:.2f}K), then evaluated "
        f"WITHOUT retuning at every span below:")
    log("")
    header = f"{'span_K':>7}{'nanocomp_COP':>14}{'nanocomp_Qc_W':>15}{'single_COP':>12}{'single_Qc_W':>13}"
    log(header)

    rows = []
    for span in off_design_spans_K:
        r_nano = run_cascade(T_cold_K, span, 1, material=nanocomposite,
                              mu0H_max=mu0H_max, mass_per_stage=mass_per_stage)
        r_single = run_cascade(T_cold_K, span, 1, material=single_phase,
                                mu0H_max=mu0H_max, mass_per_stage=mass_per_stage)
        nano_cop = r_nano["COP_cascade"] if r_nano["feasible"] else 0.0
        nano_qc = r_nano["Qc_W"] if r_nano["feasible"] else 0.0
        single_cop = r_single["COP_cascade"] if r_single["feasible"] else 0.0
        single_qc = r_single["Qc_W"] if r_single["feasible"] else 0.0
        rows.append({"span_K": span, "nanocomposite_COP": nano_cop, "nanocomposite_Qc_W": nano_qc,
                     "single_phase_COP": single_cop, "single_phase_Qc_W": single_qc,
                     "is_design_span": span == design_span_K})
        log(f"{span:>7.1f}{nano_cop:>14.3f}{nano_qc:>15.1f}{single_cop:>12.3f}{single_qc:>13.1f}"
            f"{'  <- design span' if span == design_span_K else ''}")

    off_design_rows = [r for r in rows if not r["is_design_span"]]
    both_feasible = [r for r in off_design_rows
                      if r["nanocomposite_Qc_W"] > 0 and r["single_phase_Qc_W"] > 0]
    nano_only_feasible = [r for r in off_design_rows
                           if r["nanocomposite_Qc_W"] > 0 and r["single_phase_Qc_W"] == 0]
    single_only_feasible = [r for r in off_design_rows
                             if r["single_phase_Qc_W"] > 0 and r["nanocomposite_Qc_W"] == 0]
    neither_feasible = [r for r in off_design_rows
                         if r["nanocomposite_Qc_W"] == 0 and r["single_phase_Qc_W"] == 0]
    nano_wins_qc = sum(1 for r in both_feasible
                        if r["nanocomposite_Qc_W"] > r["single_phase_Qc_W"])

    log("")
    log(f"Off-design feasibility (Qc>0): nanocomposite feasible at "
        f"{len(off_design_rows) - len(single_only_feasible) - len(neither_feasible)}/"
        f"{len(off_design_rows)} off-design span(s); single phase feasible at "
        f"{len(off_design_rows) - len(nano_only_feasible) - len(neither_feasible)}/"
        f"{len(off_design_rows)}.")
    if nano_only_feasible:
        log(f"At {len(nano_only_feasible)} off-design span(s) "
            f"({[r['span_K'] for r in nano_only_feasible]}), the single phase collapses to "
            f"Qc=0 (span_fraction feasibility cutoff, core/amr_cycle.py) while the "
            f"nanocomposite still delivers positive Qc -- its broader working range covers a "
            f"span the sharply-tuned single phase cannot.")
    if single_only_feasible:
        log(f"At {len(single_only_feasible)} off-design span(s) "
            f"({[r['span_K'] for r in single_only_feasible]}), the reverse holds: the single "
            f"phase remains feasible while the nanocomposite collapses to Qc=0.")

    if both_feasible:
        log(f"Where BOTH remain feasible ({len(both_feasible)} span(s)), the nanocomposite's "
            f"Qc exceeds the single phase's Qc in {nano_wins_qc}/{len(both_feasible)} case(s) "
            f"-- generally the single phase wins on raw Qc/COP where both work, since it isn't "
            f"paying the peak-height cost of blending.")

    if nano_only_feasible and not single_only_feasible:
        if both_feasible:
            both_feasible_note = (
                f" Where both remain feasible ({len(both_feasible)} span(s)), the single "
                f"phase still wins on raw Qc/COP ({nano_wins_qc}/{len(both_feasible)} for "
                f"the nanocomposite), consistent with blending trading peak height for width."
            )
        else:
            both_feasible_note = (
                " No off-design span in this set left both candidates feasible "
                "simultaneously, so no raw Qc/COP comparison could be made there; the design "
                "span itself (where the single phase is optimally tuned and wins outright) "
                "is the only point where both are feasible."
            )
        conclusion = (
            f"The nanocomposite's breadth DOES pay off in one specific, checkable way: at "
            f"{len(nano_only_feasible)} off-design span(s) narrower than the design span, "
            f"the sharply-tuned single phase collapses to Qc=0 (its no-load DeltaT_ad no "
            f"longer covers a mismatched span) while the nanocomposite -- whose working range "
            f"is deliberately spread -- still delivers positive Qc."
            + both_feasible_note +
            " So the genuine finding here is ROBUSTNESS TO NARROWING (avoiding a "
            "catastrophic Qc=0 failure), not raw performance, and this is a first-pass "
            "finding at ONE spread value and ONE design/off-design span set, not a general "
            "claim."
        )
    elif not off_design_rows:
        conclusion = "No off-design spans were checked (off_design_spans_K only contained the design span)."
    else:
        conclusion = (
            f"No clean robustness story emerges from this specific span set: nanocomposite-"
            f"only feasibility occurred at {len(nano_only_feasible)} span(s), single-phase-"
            f"only feasibility at {len(single_only_feasible)}, both-feasible at "
            f"{len(both_feasible)} (nanocomposite ahead on Qc in {nano_wins_qc} of those), "
            f"and neither feasible at {len(neither_feasible)}. See the per-span table above "
            f"for the specific pattern; this should be read as a first-pass finding at ONE "
            f"spread value, not a general claim about composite regeneration."
        )
    log("CONCLUSION: " + conclusion)

    import os
    os.makedirs(os.path.dirname(out_path) or ".", exist_ok=True)
    with open(out_path, "w") as f:
        f.write("\n".join(lines) + "\n")
    if verbose:
        print(f"\nWrote {out_path}")

    return {"rows": rows, "Tc_design_K": Tc_design, "conclusion": conclusion}