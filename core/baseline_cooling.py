"""
baseline_cooling.py
====================
Reference COP models for two conventional data-center cooling technologies
used as benchmarks for magnetic cooling:

1. Vapor-compression CRAC/CRAH (air-cooled, room/row-based)
2. Direct liquid cooling (cold-plate, facility-water-cooled)

Both are modeled as fractions of the reverse-Carnot COP, with the fraction
(Lorenz/second-law efficiency) taken from published data-center cooling
efficiency studies rather than assumed:

    - Vapor-compression DX/chiller plants in data centers: second-law
      efficiency ~ 0.35-0.45 of Carnot for packaged CRAC units, up to
      ~0.5-0.55 for well-optimized chilled-water plants.
      (ASHRAE Datacom Series; Shah, Bash & Patel, "Cooling and Power
      Considerations for Chips", ASME (2004); Ebrahimi, Jones & Fleischer,
      Renew. Sustain. Energy Rev. 31 (2014) 622-638 - DC cooling review)

    - Direct liquid cooling (cold plate, facility water loop): removes the
      need for a low-temperature chiller stage entirely for a large fraction
      of the load (can use free/economizer cooling at facility-water
      temperatures ~18-32 C per ASHRAE TC9.9 W-class envelopes), so its
      *effective* COP is reported as the compressor-side COP only when
      mechanical cooling is still needed, otherwise pump-only COP is very
      high (>20). We model both the "mechanical assist" and "economizer"
      regimes.
      (ASHRAE TC9.9, "Liquid Cooling Guidelines for Datacom Equipment
      Centers", 2nd ed. 2021; Ellsworth, Iyengar et al., IEEE ITherm
      proceedings, various years)
"""

import numpy as np
from dataclasses import dataclass


@dataclass
class CoolingResult:
    technology: str
    Tc: float
    Th: float
    COP: float
    COP_carnot: float
    second_law_eff: float


def carnot_cop(Tc, Th):
    return Tc / (Th - Tc) if Th > Tc else np.inf


def vapor_compression_cop(Tc, Th, eta_2nd_law=0.42, source_note="packaged DX CRAC"):
    """eta_2nd_law default 0.42 is representative of packaged DX CRAC units
    per Ebrahimi et al. (2014) data-center cooling review; use 0.50-0.55 for
    optimized chilled-water plants."""
    cc = carnot_cop(Tc, Th)
    return CoolingResult("Vapor-compression (%s)" % source_note, Tc, Th,
                          eta_2nd_law * cc, cc, eta_2nd_law)


def liquid_cooling_cop(Tc, Th, economizer_hours_fraction=0.6,
                        mechanical_eta_2nd_law=0.42, pump_equivalent_cop=25.0):
    """Direct/indirect liquid cooling: for `economizer_hours_fraction` of
    operating hours the facility water loop can reject heat directly (dry
    cooler / cooling tower) without a compressor -- effective COP is very
    high (pump + fan power only, modeled here at pump_equivalent_cop, a
    literature-informed placeholder per ASHRAE TC9.9 W-class case studies).
    For the remaining hours mechanical cooling (chiller) is engaged at the
    same second-law efficiency as vapor compression. Result is an
    hours-weighted annual average COP."""
    cc = carnot_cop(Tc, Th)
    cop_mech = mechanical_eta_2nd_law * cc
    cop_avg = (economizer_hours_fraction * pump_equivalent_cop +
               (1 - economizer_hours_fraction) * cop_mech)
    eff_2nd_law_blended = cop_avg / cc if np.isfinite(cc) and cc > 0 else 0.0
    return CoolingResult("Liquid cooling (blended, %.0f%% economizer hrs)"
                          % (economizer_hours_fraction * 100),
                          Tc, Th, cop_avg, cc, eff_2nd_law_blended)


# ---------------------------------------------------------------------------
# Phase 21 -- passive/hybrid magnetic regenerator augmentation of a
# conventional (vapor-compression) gas cycle.
# ---------------------------------------------------------------------------
#
# HONESTY FLAG (book access -- same tier as Phases 17-20's own flags). The
# ROADMAP.md Phase 21 plan's stated data sources were Tishin & Spichkin
# (2003) Sect. 11.1 (passive magnetic regenerators used inside conventional
# gas-cycle refrigerators), Sect. 11.2.3 (magnetically-augmented gas
# regenerators) and Sect. 11.2.4 (hybrid magnetic working bodies). As
# already documented in this project's Phase 20 ROADMAP.md entry, this
# project's copy of Tishin & Spichkin (2003) is a scanned, image-only PDF
# with NO extractable text layer (confirmed again for this pass) -- none of
# those sections' own equations, reported effectiveness-vs-alignment curves,
# or COP figures could be digitized. What is implemented below is instead,
# exactly as the Phase 21 plan itself anticipated ("this doesn't need new
# physics or new benchmark data -- it recombines your existing
# mce_material.py entropy/heat-capacity curves with your existing
# baseline_cooling.py gas-cycle correlations in a new way"): a reuse of data
# this repo already computes, not a reproduction of Tishin's own numbers.
#
# Physical picture. A "passive" magnetic regenerator is not actively
# magnetized/demagnetized in step with the flow (that would make it an AMR,
# already covered by core/amr_cycle.py) -- instead its lambda-anomaly heat
# capacity near its own Curie temperature is exploited passively to boost
# the *thermal mass* of a conventional gas-cycle's internal regenerator
# (e.g. a vapor-compression system's liquid-suction heat exchanger), the
# same way core/thermal.py's regenerator_effectiveness() already treats
# higher solid-side heat capacity (relative to the fluid-side heat-capacity
# flow rate) as reducing its utilization term U and therefore raising
# effectiveness eps for the SAME mass/frequency/flow geometry. This module
# reuses that existing eps(cp_solid) relationship (via thermal.py's new,
# backward-compatible cp_solid override) rather than inventing a second
# effectiveness model.
#
# How effectiveness maps to COP. There is no digitized source (Tishin or
# otherwise) in this pass's corpus for exactly how much a real gas cycle's
# COP improves as its internal-regenerator effectiveness rises. The general
# phenomenon -- internal (liquid-suction) heat exchangers raising
# vapor-compression COP by a modest amount -- is well established in the
# refrigeration literature (see e.g. ASHRAE Handbook -- Refrigeration,
# internal heat exchanger chapter; commonly-cited ranges are on the order of
# a few percent up to roughly 10%, refrigerant- and design-dependent), so
# `MAX_COP_GAIN_AT_FULL_EFFECTIVENESS` below is set from that generic,
# non-magnetic-specific literature range as an illustrative CEILING (the
# gain if eps rose from the non-magnetic baseline all the way to eps=1),
# scaled linearly by this module's own delta_eps -- explicitly NOT a fitted
# coefficient and NOT specific to a magnetically-augmented regenerator
# (same evidentiary tier as Phase 16's hysteresis_loss_J_per_kg and Phase
# 18's actuation_energy_J_per_cycle placeholders: a round, cited-range
# illustration, not a calibrated fit).

from core.thermal import regenerator_effectiveness, CP_SOLID_GD  # noqa: E402

MAX_COP_GAIN_AT_FULL_EFFECTIVENESS = 0.08  # illustrative ceiling, see docstring above


@dataclass
class PassiveRegeneratorResult:
    technology: str
    material_name: str
    T_cold: float
    T_hot: float
    base_COP: float
    augmented_COP: float
    cop_gain_fraction: float
    eps_baseline: float
    eps_augmented: float
    delta_eps: float
    cp_solid_baseline_J_kgK: float
    cp_solid_augmented_J_kgK: float


def _mean_total_heat_capacity(material, T_cold, T_hot, H_field=0.0, n_points=25):
    """Average, over [T_cold, T_hot], of `material`'s own total (lattice +
    magnetic-anomaly) specific heat from core/mce_material.py's
    MagnetocaloricMaterial.total_heat_capacity() -- the SAME quantity every
    existing validation figure in this repo already computes (e.g.
    validation.py's Dan'kov cross-check), just averaged over a temperature
    window and reused here for a new purpose. If T_cold >= T_hot (a caller
    error) falls back to the single point T_cold rather than raising, since
    linspace(T_cold, T_hot, n) with T_hot <= T_cold is still well-defined
    for n=1 but not informative for n>1."""
    if T_hot <= T_cold:
        Ts = np.array([T_cold])
    else:
        Ts = np.linspace(T_cold, T_hot, n_points)
    C = material.total_heat_capacity(Ts, H=H_field)
    return float(np.mean(C))


def _mean_lattice_heat_capacity(material, T_cold, T_hot, n_points=25):
    """Same window-averaging as _mean_total_heat_capacity(), but the
    LATTICE-ONLY term (material.lattice_heat_capacity()) -- i.e. the same
    material with its own magnetic lambda-anomaly switched off. This is
    used as the DEFAULT baseline in passive_regenerator_augmentation()
    below (see that function's own docstring for why: comparing two
    DIFFERENT materials' total_heat_capacity against one shared flat
    constant conflates each material's own bulk lattice heat capacity --
    which varies a lot with atoms-per-formula-unit and molar mass, and has
    nothing to do with magnetic alignment -- with the actual effect this
    module is trying to isolate)."""
    if T_hot <= T_cold:
        Ts = np.array([T_cold])
    else:
        Ts = np.linspace(T_cold, T_hot, n_points)
    C = material.lattice_heat_capacity(Ts)
    return float(np.mean(C))


def passive_regenerator_augmentation(passive_regenerator_material, T_cold, T_hot,
                                      mass_regenerator=2.0, frequency=1.0, mdot=0.08,
                                      H_field=0.0, cp_solid_baseline=None,
                                      n_T_points=25):
    """Compares regenerator effectiveness (core/thermal.py's own
    regenerator_effectiveness(), same geometry/frequency/flow both times)
    for a `passive_regenerator_material`-based regenerator WITHOUT its own
    magnetic lambda-anomaly (lattice-only heat capacity -- i.e. the same
    material used far from its own Curie point, or a hypothetical
    non-magnetic analog of the same bulk composition) against the SAME
    material WITH its own total_heat_capacity(T) (lattice + magnetic
    anomaly), both averaged over [T_cold, T_hot].

    IMPLEMENTATION NOTE (a real decision made during this phase, stated
    rather than left implicit): the baseline deliberately uses the SAME
    material's own lattice_heat_capacity(), not a shared flat constant
    (e.g. CP_SOLID_GD) applied across every candidate material. An early
    version of this function used the latter and produced a materially
    wrong ranking -- e.g. La0.7Ca0.3MnO3 (5 atoms per formula unit, a much
    larger Dulong-Petit lattice heat capacity per kg than Gd) ranked ABOVE
    Gd even though Gd's own Curie temperature (294K) sits inside the
    representative [291.15K, 301.15K] operating window and
    La0.7Ca0.3MnO3's (267K) does not -- because the flat-constant
    comparison was mostly measuring each material's bulk lattice
    properties, not magnetic alignment. Using each material's own
    lattice-only heat capacity as its own baseline isolates the magnetic-
    anomaly contribution specifically, matching the Phase 21 plan's own
    framing ("the augmentation factor is a function of how well the
    material's heat capacity peak aligns with the gas cycle's cold-end
    temperature") rather than rewarding heavy, multi-atom formula units.
    `cp_solid_baseline` remains available as an explicit override (e.g.
    CP_SOLID_GD) for a caller who specifically wants the flat-reference
    comparison instead; default None selects the same-material lattice-only
    baseline described above."""
    cp_augmented = _mean_total_heat_capacity(passive_regenerator_material, T_cold, T_hot,
                                              H_field=H_field, n_points=n_T_points)
    cp_baseline = (cp_solid_baseline if cp_solid_baseline is not None else
                   _mean_lattice_heat_capacity(passive_regenerator_material, T_cold, T_hot,
                                                n_points=n_T_points))
    eps_baseline = regenerator_effectiveness(mass_regenerator, frequency, mdot,
                                              cp_solid=cp_baseline)["eps"]
    eps_augmented = regenerator_effectiveness(mass_regenerator, frequency, mdot,
                                               cp_solid=cp_augmented)["eps"]
    return {
        "cp_solid_baseline_J_kgK": cp_baseline,
        "cp_solid_augmented_J_kgK": cp_augmented,
        "eps_baseline": eps_baseline,
        "eps_augmented": eps_augmented,
        "delta_eps": eps_augmented - eps_baseline,
    }


def augmented_regenerator_cop(base_cop, passive_regenerator_material, T_range,
                               mass_regenerator=2.0, frequency=1.0, mdot=0.08,
                               H_field=0.0, cp_solid_baseline=None,
                               max_cop_gain_at_full_effectiveness=MAX_COP_GAIN_AT_FULL_EFFECTIVENESS,
                               n_T_points=25):
    """Phase 21 deliverable, named and shaped exactly as the plan
    specified: `augmented_regenerator_cop(base_cop, passive_regenerator_material,
    T_range)`. `base_cop` is normally `vapor_compression_cop(Tc, Th).COP`
    (this module's own existing function -- see main.py step 15 for a
    worked call). `T_range=(T_cold, T_hot)` sets both the regenerator's own
    operating window (for _mean_total_heat_capacity's alignment check) and
    the reported Tc/Th on the returned result. Returns a
    PassiveRegeneratorResult; augmented_COP >= base_COP always, since
    delta_eps is clipped at 0 (a passive regenerator whose material is
    poorly aligned with the operating window should reduce to the
    unaugmented baseline, not degrade it -- this module does not attempt to
    model a *worse*-than-baseline passive regenerator)."""
    T_cold, T_hot = T_range
    aug = passive_regenerator_augmentation(
        passive_regenerator_material, T_cold, T_hot, mass_regenerator, frequency, mdot,
        H_field, cp_solid_baseline, n_T_points)
    cop_gain_fraction = max_cop_gain_at_full_effectiveness * max(aug["delta_eps"], 0.0)
    augmented_cop = base_cop * (1.0 + cop_gain_fraction)
    return PassiveRegeneratorResult(
        technology="Vapor-compression + passive magnetic regenerator (%s)"
                   % passive_regenerator_material.name,
        material_name=passive_regenerator_material.name,
        T_cold=T_cold, T_hot=T_hot,
        base_COP=base_cop, augmented_COP=augmented_cop,
        cop_gain_fraction=cop_gain_fraction,
        eps_baseline=aug["eps_baseline"], eps_augmented=aug["eps_augmented"],
        delta_eps=aug["delta_eps"],
        cp_solid_baseline_J_kgK=aug["cp_solid_baseline_J_kgK"],
        cp_solid_augmented_J_kgK=aug["cp_solid_augmented_J_kgK"],
    )


# ---------------------------------------------------------------------------
# Phase 23 -- elastocaloric energy conversion as a static comparison
# reference point (NOT a simulated system -- see docstring below for why).
# ---------------------------------------------------------------------------
#
# HONESTY FLAG (book access -- same tier as Phases 17-22's own flags).
# phase_plan.md's own Phase 23 data source is Kitanovski et al. (2015)
# Ch. 10 Sect. 10.3 ("Elastocaloric Energy Conversion", per that book's own
# table of contents pp. 438-446). Checked directly for this pass: this
# project's copy of Kitanovski et al. (2015) is a 30-page excerpt (cover,
# front matter, preface, full table of contents, and the opening pages of
# Chapter 1 only, confirmed by pdfplumber's own extracted page count) --
# Chapter 10 is simply not present in the file this repo has access to.
# Tishin & Spichkin (2003) does not cover elastocalorics at all (the book
# predates the field's modern development; no elastocaloric chapter appears
# in its own table of contents) and is separately an image-only scan with
# no text layer besides (already re-confirmed in Phases 20-22's own
# honesty flags). So, exactly as phase_plan.md's own Phase 23 entry
# anticipated ("using published elastocaloric COP/exergy-efficiency
# figures as a static reference point"), the values below come from
# external, independently-published, peer-reviewed literature located by
# this pass's own search -- NOT from either of this repo's two source
# books, and NOT reproducing any specific book's own numbers.
#
# Sourced anchors (device/SYSTEM-level COP -- not the narrower COP_mat
# material figure-of-merit some elastocaloric papers report separately,
# which is not directly comparable to this repo's own COP_electrical):
#   - Qian, Catalini, Muehlbauer, Liu, Mevada, Hou, Hwang, Radermacher &
#     Takeuchi, "High-performance multimode elastocaloric cooling system",
#     Science 380, 722-727 (2023): simulated steady-state SYSTEM COP = 5.8,
#     at up to a 22.5 K span -- the closest published anchor to this
#     repo's own 5-20K ASHRAE sweep range. The paper's own initial
#     hardware measurements did not yet include every system loss, so the
#     authors themselves report this simulated figure as the one still
#     requiring future experimental verification -- stated here rather
#     than smoothed over.
#   - Li, Hua & Sun, "Continuous and efficient elastocaloric air cooling by
#     coil-bending", Nat. Commun. 14, 7982 (2023), DOI: 10.1038/s41467-
#     023-43611-6: MEASURED device-level system COP = 3.7, but at a much
#     narrower ~0.9-1 K temperature drop -- a genuinely different, far
#     smaller operating span than this repo's own sweep, so this anchor
#     is weaker evidence for this repo's own span range specifically than
#     the Qian et al. figure is. Kept as the LOW end of the reported
#     range for that reason, not averaged in as an equal-weight data
#     point. Phase 30 citation-audit correction: this was previously
#     misattributed here as "Wu et al." -- the paper's actual authors are
#     Xueshi Li, Peng Hua & Qingping Sun (Hong Kong Univ. of Science and
#     Technology); no author named Wu appears on it. The journal, volume,
#     page, year, and COP=3.7/~1K-span figures were all independently
#     re-verified via web search this pass and are correct; only the
#     author name was wrong.
#
# What this deliberately is NOT: a span-dependent model, and NOT a
# simulation of an elastocaloric AMR-analog cycle the way core/amr_cycle.py
# simulates magnetic AMR. Neither literature anchor above reports a
# COP(span) curve over anything resembling this repo's own 5-20K sweep
# (Qian et al.'s own device reaches spans up to 22.5K but at a single
# reported operating COP, not a swept curve; Li, Hua & Sun's own COP=3.7
# figure is reported at a ~1K span). Fitting a COP(span) curve from two
# single-point anchors measured at different, non-overlapping spans would
# invent a slope this repo has no data for. So, exactly per phase_plan.md's
# own framing ("similar treatment to how Carnot COP is already just a
# reference line in plots.py fig08, not a simulated system"), the function
# below returns ONE flat representative COP (plus the [low, high]
# literature range it was drawn from) meant to be plotted as a single
# horizontal reference line across the whole span sweep -- not a per-span
# calculation, and not claimed to be one.

ELASTOCALORIC_COP_LOW = 3.7    # Li, Hua & Sun (2023) Nat. Commun. 14, 7982 -- MEASURED, ~1K span
ELASTOCALORIC_COP_HIGH = 5.8   # Qian et al. (2023) Science 380, 722-727 -- SIMULATED, up to 22.5K span
ELASTOCALORIC_COP_SOURCE_NOTE = (
    "External literature (NOT this repo's two source books -- see "
    "core/baseline_cooling.py's Phase 23 honesty flag): Qian et al. (2023) "
    "Science 380, 722-727, DOI 10.1126/science.adg7043 (simulated "
    "steady-state system COP=5.8, up to a 22.5K span) and Li, Hua & Sun "
    "(2023) Nat. Commun. 14, 7982, DOI 10.1038/s41467-023-43611-6 "
    "(measured system COP=3.7, at a much narrower ~1K span; corrected "
    "Phase 30 -- previously misattributed here as 'Wu et al.'). A single "
    "static reference point, not a span-dependent simulation."
)


@dataclass
class ElastocaloricReferenceResult:
    technology: str
    COP_representative: float
    COP_low: float
    COP_high: float
    source_note: str


def elastocaloric_reference_cop():
    """Phase 23 deliverable, named and shaped as the plan's own
    `elastocaloric_reference_cop()` lookup (not a simulation -- see the
    module-level honesty flag above). Returns a STATIC literature
    reference point for use as a fourth comparison entry alongside
    Carnot / vapor-compression / liquid-cooling, e.g. as a horizontal
    reference line in core/plots.py's fig08 and as a reference row/column
    in main.py step 4's own comparison table.

    Takes no arguments (T_cold/T_hot/span) on purpose: unlike
    vapor_compression_cop()/liquid_cooling_cop(), which compute an actual
    Carnot-fraction COP from the requested temperatures, this function has
    no digitized or literature-sourced elastocaloric COP(span) relation to
    evaluate (see honesty flag) -- it always returns the same static
    reference regardless of what temperatures a caller has in mind, and
    callers should read COP_low/COP_high before treating COP_representative
    as more precise than the literature it was drawn from actually is.

    COP_representative is the geometric mean of the two sourced anchors
    (chosen over an arithmetic mean since COP is a ratio-like efficiency
    figure -- the physically meaningful "middle" of two COP values is
    their geometric, not arithmetic, mean)."""
    cop_rep = float(np.sqrt(ELASTOCALORIC_COP_LOW * ELASTOCALORIC_COP_HIGH))
    return ElastocaloricReferenceResult(
        technology="Elastocaloric (literature reference)",
        COP_representative=cop_rep,
        COP_low=ELASTOCALORIC_COP_LOW,
        COP_high=ELASTOCALORIC_COP_HIGH,
        source_note=ELASTOCALORIC_COP_SOURCE_NOTE,
    )