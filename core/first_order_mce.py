"""
first_order_mce.py
====================
Implements an extended Landau free-energy model for Gd5Si2Ge2, replacing
the mean-field/Brillouin treatment with a model capable of describing the
first-order magnetostructural transition responsible for the giant
magnetocaloric effect (Pecharsky & Gschneidner, Phys. Rev. Lett. 78, 4494
(1997)).

Model
-----
Free energy per magnetic ion, in units of k_B*Tc, as a function of reduced
magnetization m = M/M_sat and reduced temperature tau = T/Tc:

    f(m, tau, h) = (A/2)*(tau-1)*m^2 + (B/4)*m^4 + (C/6)*m^6 - h*m

This is the standard "extended Landau theory" (quadratic-quartic-sextic
expansion) used for first-order magnetocaloric materials when B < 0 (the
quartic term must be negative to produce a discontinuous jump in the
equilibrium m at h=0, tau=1; the positive sextic term C is required for
free-energy stability). This is the same qualitative approach as the
classical Bean & Rodbell (1962) volume-strain-coupled model, but expressed
directly as a phenomenological order-parameter expansion rather than via
self-consistent lattice strain -- a standard simplification used e.g. in
de Oliveira & von Ranke's review (Phys. Rep. 489 (2010) 89-159, Section 4)
of first-order MCE models, which is also the source cited for the
mean-field-theory limitations flagged elsewhere in this codebase.

Reduced field:
    h = g*J*mu_B*(mu0*H) / (kB*Tc)     (same natural scale as the Brillouin
                                          argument used in mce_material.py)

Equilibrium m(tau,h): real root of  A*(tau-1)*m + B*m^3 + C*m^5 = h  that
GLOBALLY MINIMIZES f (not just any stationary point) -- this correctly
selects the equilibrium/reversible branch through the first-order jump,
consistent with treating S as a state function for the ideal AMR cycle
(hysteresis/irreversibility at the transition is a real effect this
simplified treatment does not capture -- flagged below).

Entropy (envelope theorem, since df/dm=0 at equilibrium):
    S(tau,h)/  (N*kB) = -(A/2) * m(tau,h)^2   [+ const, cancels in DeltaS]
    DeltaS_M(tau,h) = -(A/2) * [m(tau,h)^2 - m(tau,0)^2] * N*kB

Calibration: (A, B, C) = (10, -4, 8) were found by grid search to reproduce
the most consistently cited literature peak entropy change for Gd5Si2Ge2,
|DeltaS_M| ~ 18 J/(kg K) at mu0*DeltaH = 5 T near Tc = 276-278 K (Gschneidner
& Pecharsky review value, cross-checked against the original 1997 PRL and
the Journal of Superconductivity and Novel Magnetism (2019) DFT+Monte Carlo
reproduction, which reports 9.97 J/(kg K) at 5 T -- literature itself spans
roughly 10-18.5 J/(kg K) depending on sample preparation/purity, per
Pecharsky & Gschneidner's own note that "optimally prepared" samples show
the largest effect. The calibration here targets the upper (widely-quoted
"~18 J/kgK") end -- treat this as ONE defensible calibration choice, not
the only one, and note the real spread when citing this in the paper.

**Honesty flags**:
  1. delta_T_adiabatic here uses DeltaT ~ -T*DeltaS_M/C_lattice, using ONLY
     the Debye lattice heat capacity (unlike GADOLINIUM's second-order
     treatment in mce_material.py, which also adds a magnetic lambda-anomaly
     C_mag). That correction is appropriate for continuous transitions; for
     a first-order transition the physically correct denominator involves
     the transition's latent heat structure, which this 0-D model does not
     resolve. Treat ΔT_ad from this module as an upper-bound-ish estimate,
     not a validated number -- there is no direct system-level benchmark
     for it in this codebase (see data/amr_experimental_benchmarks.csv,
     which has none for Gd5Si2Ge2).
    2. UPDATE -- this cross-check has now been done (see
    core/giguere_validation.py, results/giguere_validation.txt). This
    model's peak DeltaT_ad at 7 T (~24 K) overestimates Giguere et al.'s
    (Phys. Rev. Lett. 83, 2262 (1999)) DIRECTLY measured value (10.0 K) by
    ~2.4x -- worse than the ~1.49x gap that paper itself found between
    Maxwell-relation ("indirect") and direct DeltaT_ad for this same
    first-order transition. The two effects are consistent, not
    contradictory: Giguere et al.'s finding explains PART of the gap
    (this model is calibrated to the same kind of indirect DeltaS_M value
    their Maxwell-relation method produces); honesty flag #1's lattice-only
    C_lattice denominator plausibly explains the rest. This model was NOT
    refit to close the gap -- doing so would abandon its documented
    calibration to the peak DeltaS_M literature value, and a 0-D
    lattice-only-C_p framework cannot match both simultaneously (that needs
    the transition's latent heat, out of scope here). Instead, an
    EMPIRICAL correction factor derived from this one cross-check point
    (core.giguere_validation.DTAD_CORRECTION_FACTOR, ~0.41) is available
    via the `dTad_correction` field below and is applied by default in
    `composition_tuned_material()`, so downstream design conclusions (e.g.
    cascade.py's Curie-graded cascade) are not built on the ~2.4x-optimistic
    raw number. Treat that correction as an honest fudge factor from a
    single field/composition point, not a validated model extension.
"""

import numpy as np
from dataclasses import dataclass

kB = 1.380649e-23
muB = 9.2740100783e-24
NA = 6.02214076e23
mu0 = 4 * np.pi * 1e-7


@dataclass
class FirstOrderMCEMaterial:
    name: str
    Tc: float
    J: float
    g: float
    M_molar: float
    theta_D: float
    n_atoms_per_fu: int
    A: float
    B: float
    C: float
    source: str = ""
    dTad_correction: float = 1.0
    # Multiplicative correction applied to delta_T_adiabatic() only (NOT to
    # delta_S_isothermal(), which Giguere et al. (1999) show is reasonably
    # consistent with independent Clausius-Clapeyron entropy for this
    # transition -- the discrepancy this module's honesty flag #2 documents
    # is specifically in DeltaT_ad, i.e. in the C_lattice-only denominator).
    # Default 1.0 preserves the module's originally documented (uncorrected)
    # calibration for GD5SI2GE2_FIRST_ORDER below; composition_tuned_material()
    # sets this to core.giguere_validation.DTAD_CORRECTION_FACTOR by default.

    def __post_init__(self):
        self.N = NA / self.M_molar

    def _h_reduced(self, mu0H):
        return self.g * self.J * muB * mu0H / (kB * self.Tc)

    def _equilibrium_m(self, tau, h):
        coeffs = [self.C, 0, self.B, 0, self.A * (tau - 1), -h]
        roots = np.roots(coeffs)
        real_roots = roots[np.abs(roots.imag) < 1e-6].real
        real_roots = real_roots[np.abs(real_roots) <= 1.5]
        if len(real_roots) == 0:
            return 0.0

        def f(m):
            return (0.5 * self.A * (tau - 1) * m ** 2 + 0.25 * self.B * m ** 4
                    + (self.C / 6) * m ** 6 - h * m)
        vals = [f(m) for m in real_roots]
        return real_roots[int(np.argmin(vals))]

    def delta_S_isothermal(self, T, H_final, H_initial=0.0):
        T = np.atleast_1d(np.asarray(T, dtype=float))
        out = np.zeros_like(T)
        for i, Ti in enumerate(T):
            tau = Ti / self.Tc
            # H_final/H_initial are passed in A/m (same convention as
            # mce_material.py), so mu0*H converts to Tesla for _h_reduced
            h_f = self._h_reduced(mu0 * H_final)
            h_i = self._h_reduced(mu0 * H_initial)
            m_f = self._equilibrium_m(tau, h_f)
            m_i = self._equilibrium_m(tau, h_i)
            s_f = -0.5 * self.A * m_f ** 2
            s_i = -0.5 * self.A * m_i ** 2
            out[i] = (s_f - s_i) * self.N * kB
        return out

    def lattice_heat_capacity(self, T, n_debye_points=400):
        T = np.atleast_1d(np.asarray(T, dtype=float))
        R = 8.314462618
        c_molar = np.zeros_like(T)
        for i, Ti in enumerate(T):
            Ti = max(Ti, 1.0)
            xmax = self.theta_D / Ti
            xs = np.linspace(1e-4, xmax, n_debye_points)
            integrand = (xs ** 4 * np.exp(xs)) / (np.expm1(xs) ** 2)
            trapz_fn = getattr(np, "trapezoid", None) or np.trapz
            integral = trapz_fn(integrand, xs)
            c_molar[i] = 9 * self.n_atoms_per_fu * R * (Ti / self.theta_D) ** 3 * integral
        return c_molar / self.M_molar

    def delta_T_adiabatic(self, T, H_final, H_initial=0.0):
        T = np.atleast_1d(np.asarray(T, dtype=float))
        dS = self.delta_S_isothermal(T, H_final, H_initial)
        C = self.lattice_heat_capacity(T)
        return self.dTad_correction * (-T * dS / C)


GD5SI2GE2_FIRST_ORDER = FirstOrderMCEMaterial(
    name="Gd5Si2Ge2 (first-order Landau model)",
    Tc=276.0, J=3.5, g=2.0,
    M_molar=(5 * 157.25 + 2 * 28.085 + 2 * 72.63) * 1e-3,
    theta_D=200.0, n_atoms_per_fu=9,
    A=10.0, B=-4.0, C=8.0,
    source="Landau coefficients calibrated to peak |DeltaS_M|~18 J/(kg K) "
           "at 5T (Pecharsky & Gschneidner 1997; Gschneidner & Pecharsky "
           "review); NOT independently validated against a second dataset "
           "(see module docstring honesty flag #2).",
)

# --- La(Fe,Si)13Hy (itinerant-electron metamagnetic giant-MCE family) ---
#
# This is the material actually used in the Astronautics_rotary_2014 benchmark
# row (data/amr_experimental_benchmarks.csv, Jacobs et al., Int. J. Refrig. 37
# (2014) 84-91), which validation_system.py previously ran against GADOLINIUM
# as an explicitly-flagged stand-in "because that material is not yet included
# in the material library." It now is.
#
# Composition & calibration target: La(Fe0.90Si0.10)13H1.1, Tc=287K, peak
# |DeltaS_M| ~ 31 J/(kg K) at mu0*DeltaH = 0-5T (indirect, Maxwell-relation),
# with an indirectly-estimated peak DeltaT_ad ~ 15.4 K at 5T -- Fujieda,
# Fujita & Fukamichi, Appl. Phys. Lett. 81, 1276 (2002) and related
# La(Fe,Si)13Hy literature, the composition and Tc/DeltaS_M values most
# commonly cited for the hydrogenated, room-temperature-tuned La(Fe,Si)13Hy
# family this device uses (the real Astronautics beds were SIX
# Curie-temperature-graded layers spanning roughly 304-316K -- see e.g. Bahl
# et al., Int. J. Refrig. 74 (2017) 22-29 for the layer Tc range of a similar
# Astronautics-style device; this single-Tc=287K entry is therefore a
# representative single-layer material, not a reconstruction of the actual
# 6-layer graded bed, matching the level of simplification
# validation_system.py already uses for every other single-material row in
# the benchmark CSV).
#
# Molar mass: 1 La (138.905) + 11.7 Fe (0.90*13, 55.845) + 1.3 Si
# (0.10*13, 28.0855) + 1.1 H (1.008) = 829.91 g/mol.
# n_atoms_per_fu = 15 (1 La + 13 (Fe,Si) + ~1.1 H, rounded).
#
# (A, B, C) = (15.0, -6.0, 12.0) found by the same grid-search calibration
# method as GD5SI2GE2_FIRST_ORDER (same B/A=-0.4, C/A=0.8 ratio, A rescaled),
# targeting the peak |DeltaS_M|~31 J/(kg K) at 5T above: this reproduces
# -30.9 J/(kg K) at T=297.8K (field-shifted peak, ~11K above the nominal
# Tc=287K -- same field-shift effect documented for GD5SI2GE2_FIRST_ORDER
# above and in giant_mce_analysis.py).
#
# theta_D=350K is NOT a literature-measured value for this specific
# composition -- a targeted search for this addition did not turn up a
# directly reported Debye temperature for La(Fe,Si)13Hy, only low-T specific
# heat studies (e.g. Phys. Rev. B 94, 134405 (2016)) that extract a
# Sommerfeld coefficient and note theta_D shifts with hydrogenation without
# tabulating a room-temperature-relevant number. 350K is an order-of-magnitude
# placeholder from comparable Fe-intermetallics (e.g. FeGe, theta_D~348-390K);
# treat lattice_heat_capacity() and therefore delta_T_adiabatic() for this
# material as correspondingly less trustworthy than delta_S_isothermal(),
# same caveat structure as GD5SI2GE2_FIRST_ORDER's own honesty flags.
#
# dTad_correction is left at the class default (1.0, uncorrected): the
# Giguere et al. (1999) empirical correction factor used by
# composition_tuned_material() below was derived from a Gd5Si2Ge2-specific
# direct-vs-indirect DeltaT_ad comparison and has NOT been shown to transfer
# to this different first-order compound family -- applying it here would be
# fabricating a validation that doesn't exist. With that said, this model's
# raw peak DeltaT_ad at 5T (see __main__ block below) comes out to ~21.9K,
# noticeably above the ~15.4K indirect literature estimate; the same
# lattice-only-C_p honesty flag #1 above is the most likely explanation,
# unconfirmed absent a direct measurement to check against (no Giguere-style
# cross-check dataset for this material was located for this pass).
LAFESIH_FIRST_ORDER = FirstOrderMCEMaterial(
    name="La(Fe0.90Si0.10)13H1.1 (first-order Landau model)",
    Tc=287.0, J=3.5, g=2.0,
    M_molar=(138.905 + 11.7 * 55.845 + 1.3 * 28.0855 + 1.1 * 1.008) * 1e-3,
    theta_D=350.0, n_atoms_per_fu=15,
    A=15.0, B=-6.0, C=12.0,
    source="Composition/Tc/DeltaS_M target: Fujieda, Fujita & Fukamichi, Appl. "
           "Phys. Lett. 81, 1276 (2002) and related La(Fe,Si)13Hy MCE literature "
           "(peak |DeltaS_M|~31 J/(kg K), indirect DeltaT_ad~15.4K at 5T near "
           "Tc=287K); device this calibrates against: Jacobs et al., Int. J. "
           "Refrig. 37 (2014) 84-91 (Astronautics rotary AMR). Landau (A,B,C) "
           "and theta_D are this-repo calibrations/placeholders, NOT literature "
           "values -- see the block comment above for exact provenance of each "
           "parameter and its honesty flags. NOT independently validated "
           "against a second dataset (same caveat as GD5SI2GE2_FIRST_ORDER).",
)

# --- Composition tunability of the Gd5(SixGe1-x)4 family, for the
#     Curie-graded cascade (ROADMAP.md Phase 7 open item) ---
#
# Literature range for GIANT (first-order) MCE character in this family,
# read from the sources cited in giant_mce_analysis.py plus a targeted
# search performed for this roadmap item (not fabricated from the two
# Gd5Si2Ge2/Gd5Si4 endpoints alone):
#   - Pecharsky & Gschneidner, Appl. Phys. Lett. 70, 3299 (1997): giant MCE
#     ordering temperature tunable ~20-275/290 K via the Si:Ge ratio for
#     x in [0, 0.5]; alloying with Ga extends the top of that range to
#     ~290 K.
#   - Pecharsky & Gschneidner, Phys. Rev. Lett. 78, 4494 (1997): Gd5Si4
#     (x=1) itself orders at 335 K, but that composition is OUTSIDE the
#     giant/first-order regime (giant MCE is reported for x <= 0.5) -- it
#     is a normal second-order ferromagnet, not a member of the tunable
#     giant-MCE family this cascade needs.
# So the DOCUMENTED giant-MCE Tc window for this family is ~20-290 K, not
# 20-335 K. This matters directly for the DC application: ASHRAE's 18-27C
# (291.15-300.15 K) supply range sits AT OR ABOVE this documented ceiling,
# which is exactly the tension giant_mce_analysis.py already flagged as an
# open materials-research question -- see cascade.py's run_graded_cascade().
GIANT_MCE_TC_MIN_K = 20.0
GIANT_MCE_TC_MAX_K = 290.0  # upper end only reached via Ga-alloying; ~275-276K without


def composition_tuned_material(Tc_target_K, apply_giguere_correction=True, name=None):
    """Returns a FirstOrderMCEMaterial representing a hypothetical
    composition-tuned Gd5(SixGe1-x)4(-Ga) alloy with Curie temperature
    Tc_target_K, for use in a Curie-graded cascade.

    IMPORTANT SIMPLIFYING ASSUMPTION: only Tc is shifted. The Landau
    coefficients (A, B, C), Debye temperature, molar mass and atom count
    are all held fixed at the Gd5Si2Ge2 values calibrated in this module --
    composition-specific (A, B, C) for other x are not available in the
    literature reviewed for this project. This means the PEAK DeltaS_M and
    the shape/width of the transition are assumed identical across the
    graded stages, which is very unlikely to be exactly true (Pecharsky &
    Gschneidner's own data shows the peak DeltaS_M itself varies with x,
    not just Tc). Treat this as a first-order (pun intended) approximation
    for exploring the cascade CONCEPT, not a quantitatively validated
    per-stage material model.

    Raises ValueError if Tc_target_K falls outside the documented giant-MCE
    range for this family (GIANT_MCE_TC_MIN_K to GIANT_MCE_TC_MAX_K) --
    this function will not silently extrapolate a "giant" MCE material into
    a temperature range where the literature does not support one existing.
    """
    if not (GIANT_MCE_TC_MIN_K <= Tc_target_K <= GIANT_MCE_TC_MAX_K):
        raise ValueError(
            f"Tc_target_K={Tc_target_K:.1f}K is outside the documented "
            f"giant-MCE range for the Gd5(SixGe1-x)4(-Ga) family "
            f"({GIANT_MCE_TC_MIN_K:.0f}-{GIANT_MCE_TC_MAX_K:.0f}K; Pecharsky & "
            f"Gschneidner, Appl. Phys. Lett. 70, 3299 (1997)). Gd5Si4 (x=1, "
            f"Tc=335K) is NOT a giant-MCE composition (second-order transition), "
            f"so this range cannot be extended by extrapolating toward it."
        )
    dTad_correction = 1.0
    if apply_giguere_correction:
        from core.giguere_validation import DTAD_CORRECTION_FACTOR
        dTad_correction = DTAD_CORRECTION_FACTOR
    return FirstOrderMCEMaterial(
        name=name or f"Gd5(SixGe1-x)4-type, composition-tuned to Tc={Tc_target_K:.1f}K",
        Tc=Tc_target_K, J=GD5SI2GE2_FIRST_ORDER.J, g=GD5SI2GE2_FIRST_ORDER.g,
        M_molar=GD5SI2GE2_FIRST_ORDER.M_molar, theta_D=GD5SI2GE2_FIRST_ORDER.theta_D,
        n_atoms_per_fu=GD5SI2GE2_FIRST_ORDER.n_atoms_per_fu,
        A=GD5SI2GE2_FIRST_ORDER.A, B=GD5SI2GE2_FIRST_ORDER.B, C=GD5SI2GE2_FIRST_ORDER.C,
        dTad_correction=dTad_correction,
        source="Composition-tuned analog of Gd5Si2Ge2 -- (A,B,C)/theta_D/M_molar held "
               "fixed at the Gd5Si2Ge2 calibration (approximation, see docstring); Tc "
               "tunability range from Pecharsky & Gschneidner, Appl. Phys. Lett. 70, "
               "3299 (1997); dTad_correction from core.giguere_validation "
               f"({'applied, ' + str(round(dTad_correction, 3)) if apply_giguere_correction else 'NOT applied'}).",
    )


# --- Composition tunability of the La(Fe,Si)13Hy family, for a Curie-graded
#     bed matching the REAL Astronautics_rotary_2014 device (Jacobs et al.,
#     Int. J. Refrig. 37 (2014) 84-91), which used six layers graded roughly
#     303.6-316.2K (see LAFESIH_FIRST_ORDER's block comment above; Bahl et
#     al., Int. J. Refrig. 74 (2017) 22-29, cites this range for a similar
#     Astronautics-style device). This is the direct follow-up to the Phase 9
#     finding in ROADMAP.md that a single Tc=287K material can't reproduce
#     that device -- same mechanism, and same "composition_tuned_material"
#     pattern, as the Gd5(SixGe1-x)4 family above.
#
# Tc tunability range: hydrogenation and Si content (plus Co/Mn/Al
# substitution) are reported across the La(Fe,Si)13Hy literature to tune Tc
# from as low as ~190K (unhydrogenated, low-Si La(Fe,Si)13, e.g. x~0.10-0.12,
# TC~195K -- PMC 10938420, "A Short Review on the Evolution of Magnetocaloric
# La(Fe,Si)13") up to the ~330-340K range with heavier hydrogenation/Co
# co-doping used to push the transition above room temperature for practical
# devices (see e.g. the review "La(Fe,Si/Al)13-based materials with
# exceptional magnetic functionalities", oaepublish 2024, and MDPI Magnetism
# 7(1):13 (2021) "Tuning the Magnetocaloric Properties of La(Fe,Si)13
# Compounds by Chemical Substitution and Light Element Insertion" for the
# general substitution/hydrogenation tuning mechanisms). 190-340K is used
# here as the documented window; unlike the Gd5(SixGe1-x)4 case, this was
# not independently re-derived from a single pair of endpoint papers, so
# treat the exact edges as approximate (a few tens of K), not as precisely
# sourced as GIANT_MCE_TC_MIN_K/_MAX_K above.
LAFESIH_TC_MIN_K = 190.0
LAFESIH_TC_MAX_K = 340.0


def lafesih_composition_tuned_material(Tc_target_K, name=None):
    """Returns a FirstOrderMCEMaterial representing a hypothetical
    composition/hydrogenation-tuned La(Fe,Si)13Hy alloy with Curie
    temperature Tc_target_K, for use in a Curie-graded bed -- the
    La(Fe,Si)13Hy analog of composition_tuned_material() above.

    SAME simplifying assumption as composition_tuned_material(): only Tc is
    shifted. (A, B, C), theta_D, M_molar and n_atoms_per_fu are all held
    fixed at the LAFESIH_FIRST_ORDER calibration -- composition-specific
    Landau coefficients for other La(Fe,Si)13Hy Si:H ratios were not found
    in the literature reviewed for this addition. This is very likely less
    accurate here than for the Gd5(SixGe1-x)4 case above: real La(Fe,Si)13Hy
    peak |DeltaS_M| is known to vary substantially with both Si content and
    H loading (not just Tc), and the real device's layers are also reported
    to differ in more than Tc alone. Treat this as exploring whether a
    Curie-graded LAYERED bed is even the right STRUCTURE to explain the
    Astronautics device's performance, not as a quantitatively validated
    per-layer material model.

    Unlike composition_tuned_material(), no Giguere-style empirical
    dTad_correction is applied (or available) -- LAFESIH_FIRST_ORDER's own
    dTad_correction default (1.0, uncorrected) carries through unchanged.

    Raises ValueError if Tc_target_K falls outside LAFESIH_TC_MIN_K to
    LAFESIH_TC_MAX_K.
    """
    if not (LAFESIH_TC_MIN_K <= Tc_target_K <= LAFESIH_TC_MAX_K):
        raise ValueError(
            f"Tc_target_K={Tc_target_K:.1f}K is outside the documented "
            f"tunability range for the La(Fe,Si)13Hy family "
            f"({LAFESIH_TC_MIN_K:.0f}-{LAFESIH_TC_MAX_K:.0f}K -- see the "
            f"block comment above LAFESIH_TC_MIN_K for sourcing)."
        )
    return FirstOrderMCEMaterial(
        name=name or f"La(Fe,Si)13Hy-type, composition-tuned to Tc={Tc_target_K:.1f}K",
        Tc=Tc_target_K, J=LAFESIH_FIRST_ORDER.J, g=LAFESIH_FIRST_ORDER.g,
        M_molar=LAFESIH_FIRST_ORDER.M_molar, theta_D=LAFESIH_FIRST_ORDER.theta_D,
        n_atoms_per_fu=LAFESIH_FIRST_ORDER.n_atoms_per_fu,
        A=LAFESIH_FIRST_ORDER.A, B=LAFESIH_FIRST_ORDER.B, C=LAFESIH_FIRST_ORDER.C,
        dTad_correction=LAFESIH_FIRST_ORDER.dTad_correction,
        source="Composition-tuned analog of La(Fe0.90Si0.10)13H1.1 -- (A,B,C)/"
               "theta_D/M_molar held fixed at the LAFESIH_FIRST_ORDER calibration "
               "(approximation, see docstring); Tc tunability range is a general "
               "literature reading (see LAFESIH_TC_MIN_K/_MAX_K comment above), "
               "not a single-paper-sourced window like the Gd5(SixGe1-x)4 case.",
    )


if __name__ == "__main__":
    mu0_ = 4 * np.pi * 1e-7
    print("First-order Landau model calibration check, Gd5Si2Ge2")
    print("(the transition temperature shifts with field, so the peak |DeltaS_M| "
          "is found by scanning T near Tc=276K rather than evaluating at Tc alone --")
    print(" evaluating only at the fixed nominal Tc understates the peak, since the "
          "actual field-shifted peak sits a few K above Tc; see giant_mce_analysis.py)")
    Ts = np.linspace(260.0, 300.0, 401)
    for B_T in [1, 2, 5]:
        H = B_T / mu0_
        dS_scan = GD5SI2GE2_FIRST_ORDER.delta_S_isothermal(Ts, H)
        i_peak = int(np.argmin(dS_scan))
        T_peak = Ts[i_peak]
        dS_peak = dS_scan[i_peak]
        dT_peak = GD5SI2GE2_FIRST_ORDER.delta_T_adiabatic(np.array([T_peak]), H)[0]
        dS_at_Tc = GD5SI2GE2_FIRST_ORDER.delta_S_isothermal(np.array([276.0]), H)[0]
        print(f"  {B_T}T: peak dS={dS_peak:.2f} J/(kg K) at T={T_peak:.1f}K   "
              f"dTad(at peak)={dT_peak:.2f} K   "
              f"(dS evaluated at fixed Tc=276K instead: {dS_at_Tc:.2f} J/(kg K), "
              "an understatement -- see note above)")
    print("\nTarget: peak dS ~ -18 J/(kg K) at 5T (Pecharsky & Gschneidner 1997 review "
          "value) -- matched by the scanned peak above, not by evaluating at the fixed "
          "nominal Tc.")

    print("\nFirst-order Landau model calibration check, La(Fe0.90Si0.10)13H1.1")
    Ts_la = np.linspace(272.0, 312.0, 401)
    for B_T in [1, 2, 5]:
        H = B_T / mu0_
        dS_scan = LAFESIH_FIRST_ORDER.delta_S_isothermal(Ts_la, H)
        i_peak = int(np.argmin(dS_scan))
        T_peak = Ts_la[i_peak]
        dS_peak = dS_scan[i_peak]
        dT_peak = LAFESIH_FIRST_ORDER.delta_T_adiabatic(np.array([T_peak]), H)[0]
        print(f"  {B_T}T: peak dS={dS_peak:.2f} J/(kg K) at T={T_peak:.1f}K   "
              f"dTad(at peak)={dT_peak:.2f} K")
    print("\nTarget: peak dS ~ -31 J/(kg K), dTad(indirect) ~ 15.4 K at 5T near Tc=287K "
          "(Fujieda, Fujita & Fukamichi 2002 and related La(Fe,Si)13Hy literature -- "
          "see the block comment above LAFESIH_FIRST_ORDER for the full honesty-flag "
          "list, in particular that the dTad match is expected to be worse than the "
          "dS match, same as for GD5SI2GE2_FIRST_ORDER.)")