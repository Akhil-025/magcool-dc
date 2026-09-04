"""
antiperovskite_material.py
============================
 (see ROADMAP.md): adds a fourth composition-tunable candidate
material family -- Ga1-xCMn3+x, a Mn-based antiperovskite -- alongside
GD_FAMILY / LAFESIH_FAMILY / MNFEPSI_FAMILY in core/cascade.py and
core/material_family_comparison.py.

Where this came from
---------------------
A "10 realistic materials" pitch document supplied by the user (no
citations of its own) named "GaCMn3" as a first-order antiperovskite
candidate with a transition "around 296K" and RCP~2.1 J/cm^3 at 2T.
Checked directly against the literature before writing any code (per
this repo's own no-fabricated-numbers standard -- see ROADMAP.md
entry for the full per-material verification writeup covering all 10
claims in that document, not just this one):

  - The STOICHIOMETRIC compound Mn3GaC undergoes its large-MCE
    antiferromagnetic<->ferromagnetic transition at Tt=159-165K, NOT
    296K (Lewis et al., J. Phys.: Condens. Matter (2015); Wang et al.,
    Europhys. Lett. 85, 47004 (2009); multiple independent sources agree
    on this number). The document's "around 296K" claim for this specific
    compound is wrong.
  - A carbon-deficient variant, GaC0.78Mn3, does have an FM Tc near 295K,
    but its own reported entropy change there is only -DeltaS_M=3.7 J/kg/K
    at 5T (Wang et al., arXiv:0905.1773) -- far below giant-MCE magnitude,
    and this specific off-stoichiometry composition's Tc-tunability is not
    characterized across a range in the accessible literature.
  - A DIFFERENT, well-characterized composition series, Ga(1-x)CMn(3+x)
    (x=0, 0.06, 0.07, 0.08; Wang, Tong, Sun et al., J. Appl. Phys. 105,
    083907 (2009), arXiv:0905.1777), is a genuinely good match to what the
    document was gesturing at: SECOND-ORDER (the paper explicitly reports
    "no observable hysteresis during cooling and warming... indicates the
    existence of a second-order transition"), Tc tunable from 250K (x=0)
    to 323.5K (x=0.08) by Mn-excess doping, with the x=0.07 composition
    (Tc=296.5K) giving the series' largest relative cooling power (RCP=2.1
    J/cm^3 at 45 kOe = 4.5T -- the SAME 2.1 J/cm^3 figure the document
    quoted, but at 4.5T, not 2T as the document stated). "Inexpensive and
    innoxious [non-toxic] raw materials" is the paper's own phrase.

This module implements THAT family (Ga1-xCMn3+x), not the stoichiometric
Mn3GaC the document named, since it is the literature-grounded, genuinely
room-temperature-tunable, hysteresis-free candidate -- reusing the SAME
mean-field/Brillouin machinery core.mce_material.MagnetocaloricMaterial
already implements for Gd and La0.7Ca0.3MnO3, via that class's existing
`with_Tc()` method (added for 's inhomogeneous-broadening
ensemble; reused here unchanged, no new physics machinery needed).

HONESTY FLAG -- what is and is not calibrated here
----------------------------------------------------
GADOLINIUM and LACAMNO3 (core/mce_material.py) each have (J, g, theta_D)
chosen/verified so the model's own predicted DeltaT_ad/DeltaS_M lands
close to a DIRECTLY digitized literature DeltaS_M(T,H) curve or table.
No such digit-for-digit target exists in the sources located for THIS
family: the accessible search results report Tc(x) and a single RCP
figure in J/cm^3 (not J/kg, and RCP is an integrated area, not a peak
DeltaS_M value), not a peak |DeltaS_M| in J/(kg K) at a stated field.
Converting RCP (J/cm^3) to a peak DeltaS_M (J/(kg K)) target would need
this specific alloy's density, which was not located either -- inventing
one to force a calibration would replace one unsourced number with
another. So, UNLIKE GADOLINIUM/LACAMNO3, this material's peak
DeltaS_M/DeltaT_ad MAGNITUDE is NOT calibrated against a literature
number -- only:
  (a) Tc and its tunable range (250.0-323.5K, four real measured
      compositions), and
  (b) the second-order/zero-hysteresis character
are literature-grounded. J=1.5 (an illustrative effective spin for
octahedral Mn in this antiperovskite class, of the same rough order as
LACAMNO3's own J=2.0 "effective Mn moment, approximate") and
theta_D=380K (order-of-magnitude placeholder, same range as this
repo's other 3d-transition-metal-based room-temperature MCE compounds --
GD5SI2GE2's 200K, LACAMNO3's 400K -- no direct measurement located) are
both explicitly-flagged placeholders, not fitted values. Treat this
family's Qc/COP numbers in material_family_comparison.py as a WEAKER,
uncalibrated-magnitude estimate relative to GD_FAMILY/LAFESIH_FAMILY/
MNFEPSI_FAMILY -- its Tc-window and hysteresis-free ranking are the
literature-grounded parts, not its absolute cooling capacity.

No hysteresis_loss_J_per_kg field is needed here (unlike
core.first_order_mce.FirstOrderMCEMaterial): this reuses
core.mce_material.MagnetocaloricMaterial, whose second-order/mean-field
construction is -- like GADOLINIUM -- genuinely (not approximately) free
of thermal hysteresis, which is directly confirmed for this specific
composition series by Wang et al.'s own no-observable-hysteresis finding
quoted above (not merely assumed by the model's mathematical form, as it
is for GADOLINIUM).
"""

import dataclasses

from core.mce_material import MagnetocaloricMaterial, GADOLINIUM

# Ga(1-x)C Mn(3+x), x=0.07 (the composition Wang et al. 2009 report as
# this series' RCP-maximizing, and the one closest to typical ASHRAE
# data-center supply temperatures). Molar mass computed directly from
# standard atomic weights (Ga=69.723, C=12.011, Mn=54.938 g/mol) for this
# exact stoichiometry -- not a literature figure, but not fabricated
# either (arithmetic on the composition, which itself IS literature-
# sourced), same treatment core/first_order_mce.py already gives
# GD5SI2GE2_FIRST_ORDER's M_molar.
_GA = 69.723e-3
_C = 12.011e-3
_MN = 54.938e-3
_X_REF = 0.07
M_MOLAR_GA1XCMN3X_REF = (1 - _X_REF) * _GA + 1 * _C + (3 + _X_REF) * _MN  # kg/mol

GA1XCMN3X_REF = MagnetocaloricMaterial(
    name="Ga0.93CMn3.07 (antiperovskite, x=0.07)",
    Tc=296.5,           # K, Wang et al. 2009's own measured value for x=0.07
    J=1.5,               # effective moment, ILLUSTRATIVE -- see module honesty flag
    g=2.0,
    M_molar=M_MOLAR_GA1XCMN3X_REF,
    theta_D=380.0,       # K, ILLUSTRATIVE placeholder -- see module honesty flag
    n_atoms_per_fu=5,     # (1-x)+1+(3+x) = 5 atoms/formula unit, exact for any x
    source="Wang, Tong, Sun, Zhu, Luo, Li, Song, Yang & Dai, J. Appl. Phys. 105, "
           "083907 (2009) [arXiv:0905.1777] -- Tc(x) measured for x=0, 0.06, 0.07, "
           "0.08 (250.0, 281.5, 296.5, 323.5K); explicitly reports 'no observable "
           "hysteresis during cooling and warming... second-order transition'; "
           "largest RCP=2.1 J/cm^3 at 45 kOe (4.5T) for x=0.07. J/g/theta_D are "
           "NOT from this paper -- see module docstring honesty flag: this "
           "family's Tc-tunability and hysteresis-free character are literature- "
           "grounded, its DeltaS_M/DeltaT_ad MAGNITUDE is not.",
)

# Measured endpoints of Wang et al.'s own four-composition series --
# treated as the family's documented tunability window, same convention
# core.first_order_mce.MNFEPSI_TC_MIN_K/_MAX_K already uses for a
# directly-measured (not literature-survey) composition range.
GA1XCMN3X_TC_MIN_K = 250.0   # x=0.0
GA1XCMN3X_TC_MAX_K = 323.5   # x=0.08


def ga1xcmn3x_composition_tuned_material(Tc_target_K, name=None):
    """Returns a MagnetocaloricMaterial representing a hypothetical
    Mn-excess-tuned Ga1-xCMn3+x antiperovskite with Curie temperature
    Tc_target_K, for use as a GradedFamily.tuned_fn (the antiperovskite
    analog of core.first_order_mce.mnfepsi_composition_tuned_material()).

    Reuses MagnetocaloricMaterial.with_Tc()  rather than
    reimplementing the same "hold every other parameter fixed, vary only
    Tc (and its derived Weiss constant)" pattern -- appropriate here
    specifically because this is a SECOND-ORDER family (with_Tc() was
    designed for, and is only really physically apt for, mean-field
    second-order materials -- unlike the Landau first-order families,
    whose composition-tuned_fn helpers must also carry a peak-vs-Tc
    offset and an empirical dTad_correction that with_Tc() does not
    handle).

    Raises ValueError if Tc_target_K falls outside
    GA1XCMN3X_TC_MIN_K/_MAX_K (Wang et al.'s own measured x=0 to x=0.08
    range) -- same fail-fast convention as the other three tunable
    families in this repo."""
    if not (GA1XCMN3X_TC_MIN_K <= Tc_target_K <= GA1XCMN3X_TC_MAX_K):
        raise ValueError(
            f"Tc_target_K={Tc_target_K:.1f}K is outside the directly-measured "
            f"tunability range for the Ga1-xCMn3+x family "
            f"({GA1XCMN3X_TC_MIN_K:.1f}-{GA1XCMN3X_TC_MAX_K:.1f}K -- Wang et al., "
            f"J. Appl. Phys. 105, 083907 (2009), x=0 to x=0.08)."
        )
    mat = GA1XCMN3X_REF.with_Tc(Tc_target_K)
    if name is not None:
        mat = dataclasses.replace(mat, name=name)
    return mat


if __name__ == "__main__":
    print("Ga1-xCMn3+x antiperovskite family -- Tc tunability check")
    for tc in (250.0, 281.5, 296.5, 323.5):
        m = ga1xcmn3x_composition_tuned_material(tc)
        print(f" Tc={tc:.1f}K -> {m.name}, lambda={m.lam:.4e}")
    print("\nNOTE: this family's Tc values reproduce Wang et al. (2009)'s own "
          "measured x=0/0.06/0.07/0.08 points exactly (by construction -- Tc is "
          "the only tuned parameter). Peak DeltaS_M/DeltaT_ad magnitude is NOT "
          "independently calibrated against that paper's own reported RCP -- "
          "see this module's docstring honesty flag.")