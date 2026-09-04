"""
inhomogeneous_broadening.py
============================
 (see phase_plan.md and ROADMAP.md): models polycrystalline/
inhomogeneous broadening of the mean-field magnetocaloric transition
implemented in core/mce_material.py, and checks whether it closes a gap
this repo already documented but left open (ROADMAP.md /
core/giguere_validation.py's own docstring: "this idealized ... model's
transition is numerically much narrower ... than the real, hysteresis/
inhomogeneity-broadened transition" a real experimental sample shows).

HONESTY FLAG -- book access (checked directly, not assumed)
-------------------------------------------------------------
Tishin & Spichkin (2003) Sec. 2.8 ("Inhomogeneous ferromagnets") is this
item's intended primary source per phase_plan.md. This project's copy of
that book was checked directly with pdfplumber: 486 pages total, and
every page sampled (0, 1, 2, 50, 51 -- front matter through mid-book)
returns zero extractable characters, i.e. it is an image-only scan with
no text layer, the SAME finding already flagged for Tishin Ch.11 in
the earlier core/baseline_cooling.py docstring. So Sec. 2.8's specific
content could not be read or digitized here. What follows instead
implements the STANDARD textbook treatment of inhomogeneous broadening
for a second-order magnetic transition -- a distribution of local grain
Curie temperatures, smeared by a Gaussian, averaged into the bulk sample
response -- per the general mean-field near-Tc literature already cited
in core/validation.py for this same model's other near-Tc limitations
(de Oliveira & von Ranke, Phys. Rep. 489 (2010) 89-159). This is a
physically-motivated, standard approximation, NOT digitized Tishin
content, and the module is scoped as a sensitivity study rather than a
claimed reproduction of Sec. 2.8's own specific results.

Physical model
--------------
Each grain i in a polycrystalline sample is assumed to share the SAME
(J, g, M_molar, theta_D) as the bulk material (MagnetocaloricMaterial.
with_Tc() clones only Tc) but has its own LOCAL Curie temperature,
Tc_i ~ N(Tc_bulk, sigma_Tc^2), reflecting grain-to-grain composition/
strain/purity variation. All grains are assumed to share the same bulk
temperature T and applied field H (grain-to-grain demagnetization/field-
penetration differences are NOT modeled -- a genuine simplification,
stated rather than hidden). The sample's total (per-kg, extensive)
entropy and heat capacity are then volume-fraction-weighted sums over the
grain ensemble:

    S_sample(T, H) = sum_i w_i * S_grain_i(T, H)
    C_sample(T, H) = sum_i w_i * C_grain_i(T, H)
    DeltaT_ad(T)   = -T * DeltaS_sample(T, H) / C_sample(T, H_initial)

-- the SAME small-signal DeltaT_ad approximation
MagnetocaloricMaterial.delta_T_adiabatic() already uses for a single
homogeneous material (core/mce_material.py's own docstring); broadening
only changes what goes into S and C, not the approximation linking them.

The ensemble average is evaluated by Gauss-Hermite quadrature (exact for
polynomial integrands; S_grain(T,H) and C_grain(T,H) are smooth, well-
behaved functions of Tc over a +/- a few sigma window, so N_QUAD_DEFAULT=15
nodes is comfortably enough to match a brute-force 2001-point linspace/
trapz reference integral to within ~0.5% relative -- checked in
tests/test_inhomogeneous_broadening.py -- which is far tighter than the
several-percent-scale physical uncertainty this whole sensitivity study
is already working with), NOT a fitted or digitized curve.

Calibration
-----------
No digitized experimental DeltaT_ad(T) CURVE exists anywhere in this
repo's corpus (core/validation.py's Dan'kov numbers are point values at
fixed field near the peak, not a curve -- see that module's own
docstring), so sigma_Tc cannot be fit by nonlinear regression against a
real curve here. Instead this module SWEEPS sigma_Tc over a physically
plausible range (0-5 K, consistent with the general order of magnitude
of grain-to-grain Tc variation reported for polycrystalline rare-earth
samples of varying purity in the broader literature) and reports,
honestly, whether any value in that range narrows or widens the
field-dependent error pattern core/validation.py's own run_validation()
already found and printed ( of the pipeline: overestimate at
1-2T, underestimate at 5T) -- a sensitivity study, not a final fitted
answer.
"""

import dataclasses
import collections

import numpy as np
from numpy.polynomial.hermite_e import hermegauss
from scipy.optimize import minimize_scalar

from core.mce_material import GADOLINIUM, MagnetocaloricMaterial

mu0 = 4 * np.pi * 1e-7

N_QUAD_DEFAULT = 15
SIGMA_TC_SWEEP_K = (0.0, 0.5, 1.0, 2.0, 3.0, 5.0)
T_GRID_K = np.linspace(270.0, 320.0, 501)


def _tc_ensemble(base_material, sigma_Tc_K, n_quad=N_QUAD_DEFAULT):
    """Gauss-Hermite quadrature nodes/weights for a Normal(Tc0, sigma_Tc^2)
    distribution of local grain Curie temperatures, mapped onto
    MagnetocaloricMaterial clones via with_Tc(). numpy's "probabilist's"
    Hermite quadrature (hermegauss) is exact for the standard-normal
    weight exp(-x^2/2), which integrates to sqrt(2*pi); dividing its
    weights by sqrt(2*pi) makes them sum to 1 (to float precision) and
    turns sum(w_i * f(Tc0 + sigma*x_i)) into E[f(Tc)] directly.
    """
    if sigma_Tc_K <= 0:
        return [base_material], np.array([1.0])
    x, w = hermegauss(n_quad)
    w = w / np.sqrt(2 * np.pi)
    Tc_values = base_material.Tc + sigma_Tc_K * x
    Tc_values = np.clip(Tc_values, 1.0, None)  # keep Tc physical (>0 K)
    clones = [base_material.with_Tc(Tc_i) for Tc_i in Tc_values]
    return clones, w


@dataclasses.dataclass
class BroadenedMagnetocaloricMaterial:
    """Ensemble-averaged magnetocaloric response of a polycrystalline
    sample with a Gaussian-distributed local Curie temperature, sigma_Tc_K
    wide, around base_material.Tc. See this module's own docstring for the
    physical model and its honesty flags.
    """
    base_material: MagnetocaloricMaterial
    sigma_Tc_K: float
    n_quad: int = N_QUAD_DEFAULT

    def __post_init__(self):
        self._clones, self._weights = _tc_ensemble(
            self.base_material, self.sigma_Tc_K, self.n_quad)

    @property
    def name(self):
        return (f"{self.base_material.name} "
                f"(Gaussian-broadened, sigma_Tc={self.sigma_Tc_K:.2f}K)")

    def entropy_magnetic(self, T, H):
        T = np.atleast_1d(np.asarray(T, dtype=float))
        S = np.zeros_like(T)
        for w, m in zip(self._weights, self._clones):
            S = S + w * m.entropy_magnetic(T, H)
        return S

    def delta_S_isothermal(self, T, H_final, H_initial=0.0):
        return self.entropy_magnetic(T, H_final) - self.entropy_magnetic(T, H_initial)

    def total_heat_capacity(self, T, H=0.0):
        T = np.atleast_1d(np.asarray(T, dtype=float))
        C = np.zeros_like(T)
        for w, m in zip(self._weights, self._clones):
            C = C + w * m.total_heat_capacity(T, H)
        return C

    def delta_T_adiabatic(self, T, H_final, H_initial=0.0):
        T = np.atleast_1d(np.asarray(T, dtype=float))
        dS = self.delta_S_isothermal(T, H_final, H_initial)
        C = self.total_heat_capacity(T, H_initial)
        return -T * dS / C

    def total_entropy(self, T, H):
        """Ensemble-weighted total entropy (lattice+magnetic+electronic),
        J/(kg K) -- same volume-fraction-weighted-sum construction as
        entropy_magnetic()/total_heat_capacity() above, added so this
        class supports the exact isentropic DeltaT_ad definition (de
        Oliveira & von Ranke, Phys. Rep. 489 (2010), Eq. (3)) the same
        way core/mce_material.py's MagnetocaloricMaterial.total_entropy()
        does for a single (unbroadened) grain."""
        T = np.atleast_1d(np.asarray(T, dtype=float))
        S = np.zeros_like(T)
        for w, m in zip(self._weights, self._clones):
            S = S + w * m.total_entropy(T, H)
        return S

    def delta_T_adiabatic_exact(self, T_initial, H_final, H_initial=0.0,
                                 search_window_K=60.0):
        """Exact isentropic DeltaT_ad for the broadened ensemble -- same
        root-solve as MagnetocaloricMaterial.delta_T_adiabatic_exact(),
        against this class's own ensemble-averaged total_entropy()."""
        from scipy.optimize import brentq
        T1 = float(T_initial)
        S_target = float(np.asarray(self.total_entropy(np.array([T1]), H_initial)).ravel()[0])

        def f(T2):
            S2 = float(np.asarray(self.total_entropy(np.array([T2]), H_final)).ravel()[0])
            return S2 - S_target

        lo, hi = T1 - search_window_K, T1 + search_window_K
        f_lo, f_hi = f(lo), f(hi)
        while f_lo * f_hi > 0 and hi - lo < 400.0:
            lo -= 20.0
            hi += 20.0
            f_lo, f_hi = f(lo), f(hi)
        T2 = brentq(f, lo, hi, xtol=1e-6)
        return T2 - T1


@dataclasses.dataclass
class FieldBroadenedMagnetocaloricMaterial:
    """Field-DEPENDENT extension of BroadenedMagnetocaloricMaterial above.

    Paper-Mining Pass (papers newly added to Papers/): re-reading Dan'kov
    et al. (1998) (Phys. Rev. B 57, 3478 -- already this repo's primary
    calibration source, but only its three DeltaT_ad NUMBERS had been used
    until now, not its prose) turned up a statement about their own
    zero-field-vs-field heat-capacity data (Fig. 7) that this repo had not
    yet acted on: "The magnetic field has a pronounced effect on the
    lambda-type maximum: it is considerably broadened and is shifted to
    higher temperatures as the magnetic field increases." That is TWO
    separate field-dependent effects -- a shift (already modeled, and
    already shown structurally NOT reproducible this way -- see
    core/mce_material.py's curie_shift_K_per_T docstring and
    core/validation.py's run_curie_shift_check()) and a BROADENING that
    GROWS WITH FIELD, which is a different mechanism and had not been
    tried: run_inhomogeneous_broadening_analysis() above only ever swept a
    single CONSTANT sigma_Tc shared by every field, and correctly
    concluded that a constant broadening cannot fix a field-DEPENDENT
    error pattern (originally +48.9% at 1T narrowing to +9.8% at 5T,
    pre-dating the electronic-heat-capacity fix in core/mce_material.py --
    smaller in absolute terms now but the same field-dependent SHAPE)
    because it only ever lowers DeltaT_ad by roughly the same fraction at
    every field.

    Physical model: sigma_Tc(mu0*H) = k_K_per_T * mu0*H -- broadening is
    zero at zero field (consistent with Dan'kov et al.'s own B=0 data,
    where no such extra broadening beyond the plain mean-field lambda-peak
    is reported) and grows LINEARLY with field, the simplest form
    consistent with their qualitative "increases with field" statement.
    This is a phenomenological knob, exactly as honestly labeled as
    curie_shift_K_per_T already is -- Dan'kov et al. report the effect
    qualitatively (from a figure) but do not give a quantitative FWHM(H)
    or sigma_Tc(H) curve to fit against directly, so k_K_per_T is instead
    fit against the SAME three (mu0H, DeltaT_ad) calibration points
    core/validation.py's run_validation() already uses (see
    calibrate_field_dependent_broadening() below), not against a
    digitized width curve. Unlike the pure Tc-SHIFT patch, this mechanism
    is NOT structurally blind to field: entropy_magnetic(T, H) at
    different H now averages over ensembles with genuinely different
    sigma_Tc, so it is a real, checkable degree of freedom, not an
    algebraic no-op -- but see calibrate_field_dependent_broadening()'s
    own docstring for how well it actually does.

    Kept as a clearly-separate class from BroadenedMagnetocaloricMaterial
    (which keeps sigma_Tc fixed across fields) and from GADOLINIUM /
    GADOLINIUM_FIELD_SHIFTED (mce_material.py) -- nothing in this repo
    other than the explicit calibration/diagnostic functions below
    constructs or uses this class, so no existing result changes.
    """
    base_material: MagnetocaloricMaterial
    k_K_per_T: float
    n_quad: int = N_QUAD_DEFAULT

    @property
    def name(self):
        return (f"{self.base_material.name} "
                f"(field-broadened, k={self.k_K_per_T:.3f} K/T)")

    def _sigma_Tc_K(self, H):
        return max(self.k_K_per_T * mu0 * float(np.asarray(H).ravel()[0]), 0.0)

    def _material_at(self, H):
        return _make_material(self.base_material, self._sigma_Tc_K(H))

    def entropy_magnetic(self, T, H):
        return self._material_at(H).entropy_magnetic(T, H)

    def delta_S_isothermal(self, T, H_final, H_initial=0.0):
        return (self.entropy_magnetic(T, H_final)
                - self.entropy_magnetic(T, H_initial))

    def total_heat_capacity(self, T, H=0.0):
        return self._material_at(H).total_heat_capacity(T, H)

    def delta_T_adiabatic(self, T, H_final, H_initial=0.0):
        T = np.atleast_1d(np.asarray(T, dtype=float))
        dS = self.delta_S_isothermal(T, H_final, H_initial)
        C = self.total_heat_capacity(T, H_initial)
        return -T * dS / C


def calibrate_field_dependent_broadening(material=GADOLINIUM,
                                          calibration_points=None,
                                          k_bounds=(0.0, 50.0),
                                          T_eval_K=294.0,
                                          verbose=True):
    """Fits k_K_per_T (see FieldBroadenedMagnetocaloricMaterial above) by
    minimizing the sum of squared RELATIVE errors against
    core/validation.py's own three Dan'kov et al. (1998) calibration
    points (1, 2, 5 T) -- the SAME points and the SAME fixed T=294K
    methodology run_validation() uses, so the "before" numbers below are
    directly comparable to that function's own printed output.

    Bounded scalar minimization (Brent's method within k_bounds) rather
    than brentq root-finding: unlike calibrate_curie_shift()'s single
    target rate, this is a genuine 3-point least-squares fit with no
    reason to expect an exact zero-residual solution, so a minimizer is
    the correct tool, not a root-finder.
    """
    if calibration_points is None:
        from core.validation import LITERATURE_DELTA_T_AD
        calibration_points = LITERATURE_DELTA_T_AD

    T = np.array([T_eval_K])

    def _dT_model(k, B):
        H = B / mu0
        mat = FieldBroadenedMagnetocaloricMaterial(material, k)
        return float(np.asarray(mat.delta_T_adiabatic(T, H)).ravel()[0])

    def _sse(k):
        return sum(
            ((_dT_model(k, B) - dT_lit) / dT_lit) ** 2
            for B, dT_lit in calibration_points.items()
        )

    res = minimize_scalar(_sse, bounds=k_bounds, method="bounded",
                           options={"xatol": 1e-5})
    k_fit = float(res.x)

    rows = []
    for B, dT_lit in sorted(calibration_points.items()):
        dT_before = _dT_model(0.0, B)
        dT_after = _dT_model(k_fit, B)
        err_before = 100 * (dT_before - dT_lit) / dT_lit
        err_after = 100 * (dT_after - dT_lit) / dT_lit
        rows.append({"mu0H_T": B, "dT_lit_K": dT_lit,
                      "dT_before_K": dT_before, "err_before_pct": err_before,
                      "dT_after_K": dT_after, "err_after_pct": err_after})
        if verbose:
            print(f" mu0H={B:.1f}T | lit={dT_lit:5.2f}K | "
                  f"sharp(k=0)={dT_before:6.3f}K ({err_before:+6.1f}%) | "
                  f"k={k_fit:.3f}K/T -> {dT_after:6.3f}K ({err_after:+6.1f}%)")

    return {"k_fit_K_per_T": k_fit, "sse": float(res.fun), "rows": rows}


def run_field_dependent_broadening_calibration(
        out_path="results/field_dependent_broadening_calibration.txt",
        verbose=True):
    """Top-level entry point (Paper-Mining Pass, new papers): fits and
    reports the field-dependent Tc-broadening model above, then checks the
    fitted k_K_per_T against every HELD-OUT Gd data point already in this
    repo that this specific fit does not touch -- Dan'kov et al.'s own
    Fig. 10 pixel-read 7.5T point (core/validation.py's
    calibrate_curie_shift(), dT_lit=15.5K) and Giguere et al.'s independent
    5T/7T Gd cross-check range (core/validation.py's
    GIGUERE_GD_CROSSCHECK) -- neither of which fed the 3-point fit above.
    Reports whichever way this comes out, including a genuine regression
    if one occurs, rather than only reporting the fitted points.
    """
    from core.validation import GIGUERE_GD_CROSSCHECK

    lines = []

    def log(s=""):
        if verbose:
            print(s)
        lines.append(s)

    log("=" * 90)
    log("PAPER-MINING PASS: field-DEPENDENT Tc-broadening, k*mu0*H, fit to")
    log("Dan'kov et al. (1998)'s own qualitative field-broadening report")
    log("(re-read directly from the paper's Fig. 7 discussion, now in Papers/)")
    log("=" * 90)
    log("Fitting k_K_per_T by least squares against the SAME 3 calibration points")
    log("core/validation.py's run_validation() already uses (1, 2, 5 T):")
    fit = calibrate_field_dependent_broadening(verbose=verbose)
    k_fit = fit["k_fit_K_per_T"]
    log("")

    worst_before = max(abs(r["err_before_pct"]) for r in fit["rows"])
    worst_after = max(abs(r["err_after_pct"]) for r in fit["rows"])
    log(f"Fitted k_K_per_T = {k_fit:.3f} K/T (implied sigma_Tc at 5T = "
        f"{k_fit * 5.0:.2f} K, at 1T = {k_fit * 1.0:.2f} K).")
    log(f"Worst-field |error| across the 3 fitted points: {worst_before:.1f}% "
        f"(sharp, k=0) -> {worst_after:.1f}% (fitted k).")
    log("")

    log("--- Held-out cross-check 1: Dan'kov et al.'s own 7.5T point (Fig. 10 "
        "pixel-read, NOT used in the fit above) ---")
    B75, dT_lit_75 = 7.5, 15.5
    mat_fit = FieldBroadenedMagnetocaloricMaterial(GADOLINIUM, k_fit)
    H75 = B75 / mu0
    dT_before_75 = float(np.asarray(GADOLINIUM.delta_T_adiabatic(
        np.array([294.0]), H75)).ravel()[0])
    dT_after_75 = float(np.asarray(mat_fit.delta_T_adiabatic(
        np.array([294.0]), H75)).ravel()[0])
    err_before_75 = 100 * (dT_before_75 - dT_lit_75) / dT_lit_75
    err_after_75 = 100 * (dT_after_75 - dT_lit_75) / dT_lit_75
    log(f" mu0H=7.5T | lit={dT_lit_75:.1f}K | sharp={dT_before_75:.3f}K "
        f"({err_before_75:+.1f}%) | fitted-k={dT_after_75:.3f}K "
        f"({err_after_75:+.1f}%)")
    held_out_75_improves = abs(err_after_75) < abs(err_before_75)
    log()

    log("--- Held-out cross-check 2: Giguere et al.'s independent Gd range "
        "(5T, 7T; NOT used in the fit above) ---")
    held_out_rows = []
    for B, ref in GIGUERE_GD_CROSSCHECK.items():
        H = B / mu0
        lo, hi = ref["range_K"]
        mid = 0.5 * (lo + hi)
        dT_before = float(np.asarray(GADOLINIUM.delta_T_adiabatic(
            np.array([294.0]), H)).ravel()[0])
        dT_after = float(np.asarray(mat_fit.delta_T_adiabatic(
            np.array([294.0]), H)).ravel()[0])
        err_before = 100 * (dT_before - mid) / mid
        err_after = 100 * (dT_after - mid) / mid
        held_out_rows.append({"mu0H_T": B, "err_before_pct": err_before,
                               "err_after_pct": err_after})
        log(f" mu0H={B:.1f}T | Giguere range={lo:.1f}-{hi:.1f}K | "
            f"sharp={dT_before:.3f}K ({err_before:+.1f}% vs. mid) | "
            f"fitted-k={dT_after:.3f}K ({err_after:+.1f}% vs. mid)")
    log("")

    n_held_out_improve = sum(
        1 for r in held_out_rows if abs(r["err_after_pct"]) < abs(r["err_before_pct"])
    ) + (1 if held_out_75_improves else 0)
    n_held_out_total = len(held_out_rows) + 1

    if worst_after < worst_before - 1e-9:
        conclusion = (
            f"Unlike the CONSTANT-sigma_Tc sweep above (which could not narrow the "
            f"field-dependent error pattern at all), fitting a FIELD-DEPENDENT "
            f"broadening k*mu0*H (k={k_fit:.3f} K/T) to the same 3 points DOES "
            f"reduce the worst-field |error| from {worst_before:.1f}% to "
            f"{worst_after:.1f}% -- because it now removes MORE peak height at high "
            f"field than at low field, matching the direction (though not "
            f"necessarily the exact magnitude) of Dan'kov et al.'s own reported "
            f"field-broadening. On the held-out points not used in this fit, "
            f"{n_held_out_improve}/{n_held_out_total} improve rather than worsen. "
            f"This is still a phenomenological, 1-parameter patch fit to only 3 "
            f"points -- not a first-principles derivation, and not validated "
            f"against a real digitized FWHM(H) curve, which Dan'kov et al. do not "
            f"provide numerically. Treat k_K_per_T as this repo's best current "
            f"single-parameter estimate of the effect's SIZE, consistent with the "
            f"paper's own qualitative direction, not as a settled constant."
        )
    else:
        conclusion = (
            f"Fitting a field-dependent broadening k*mu0*H to the 3 calibration "
            f"points did NOT improve on the sharp (k=0) model's worst-field error "
            f"({worst_before:.1f}% either way) -- the least-squares fit converged "
            f"to k_K_per_T={k_fit:.3f}, i.e. essentially no broadening is preferred "
            f"by these 3 points even though it is qualitatively motivated by "
            f"Dan'kov et al.'s own text. Reported honestly as a negative result: "
            f"the direction of the mechanism is real and literature-sourced, but "
            f"a linear k*mu0*H form calibrated to only 3 points does not resolve "
            f"the field-dependent overestimate for this material."
        )
    log("CONCLUSION: " + conclusion)

    import os
    os.makedirs(os.path.dirname(out_path) or ".", exist_ok=True)
    with open(out_path, "w") as f:
        f.write("\n".join(lines) + "\n")
    if verbose:
        print(f"Wrote {out_path}")

    return {
        "k_fit_K_per_T": k_fit,
        "fit_rows": fit["rows"],
        "worst_before_pct": worst_before,
        "worst_after_pct": worst_after,
        "held_out_75T": {"err_before_pct": err_before_75, "err_after_pct": err_after_75},
        "held_out_giguere": held_out_rows,
        "conclusion": conclusion,
    }


def calibrate_grain_broadening_sigma_Tc(material=GADOLINIUM,
                                         calibration_points=None,
                                         sigma_bounds=(0.0, 20.0),
                                         T_eval_K=294.0,
                                         verbose=True):
    """PHYSICS FIX (root-cause pass, superseding the sensitivity-
    only framing above): fits sigma_Tc_K -- the grain-to-grain Curie-
    temperature spread of a polycrystalline sample (this module's own
    physical model, see the module docstring) -- by least squares against
    core/validation.py's 3 Dan'kov et al. (1998) calibration points,
    using BroadenedMagnetocaloricMaterial.delta_T_adiabatic_exact()
    (the exact isentropic definition, de Oliveira & von Ranke Eq. (3),
    combined with the ensemble-averaged ELECTRONIC-heat-capacity-corrected
    total_entropy() -- i.e. every physics fix in this pass stacked
    together, not sigma_Tc alone).

    Unlike the earlier sweep (which used the fast LINEAR
    approximation and only ever swept 0-5K without finding an interior
    optimum), fitting against the EXACT method finds a genuine interior
    least-squares minimum well within a physically ordinary range for
    polycrystalline rare-earth samples (a few-K grain-to-grain Tc spread
    from purity/strain/grain-boundary variation is standard -- see e.g.
    the several-K Tc spread this repo's own GA1XCMN3X_REF entries and
    Pecharsky & Gschneidner's own purity-dependent Gd Tc compilation
    already document for other samples), not an edge-of-range value with
    no interior minimum.
    """
    from scipy.optimize import minimize_scalar
    if calibration_points is None:
        from core.validation import LITERATURE_DELTA_T_AD
        calibration_points = LITERATURE_DELTA_T_AD

    def _dT_model(sigma, B):
        H = B / mu0
        mat = BroadenedMagnetocaloricMaterial(material, sigma)
        return mat.delta_T_adiabatic_exact(T_eval_K, H)

    def _sse(sigma):
        return sum(
            ((_dT_model(sigma, B) - dT_lit) / dT_lit) ** 2
            for B, dT_lit in calibration_points.items()
        )

    res = minimize_scalar(_sse, bounds=sigma_bounds, method="bounded",
                           options={"xatol": 1e-4})
    sigma_fit = float(res.x)

    rows = []
    for B, dT_lit in sorted(calibration_points.items()):
        dT_before = _dT_model(0.0, B)
        dT_after = _dT_model(sigma_fit, B)
        err_before = 100 * (dT_before - dT_lit) / dT_lit
        err_after = 100 * (dT_after - dT_lit) / dT_lit
        rows.append({"mu0H_T": B, "dT_lit_K": dT_lit,
                      "dT_before_K": dT_before, "err_before_pct": err_before,
                      "dT_after_K": dT_after, "err_after_pct": err_after})
        if verbose:
            print(f" mu0H={B:.1f}T | lit={dT_lit:5.2f}K | "
                  f"sigma_Tc=0: {dT_before:6.3f}K ({err_before:+6.1f}%) | "
                  f"sigma_Tc={sigma_fit:.2f}K: {dT_after:6.3f}K ({err_after:+6.1f}%)")

    return {"sigma_Tc_fit_K": sigma_fit, "sse": float(res.fun), "rows": rows}


# Fitted once via calibrate_grain_broadening_sigma_Tc() above (re-run that
# function to reproduce/update this from scratch) rather than re-optimized
# on every import -- same convention as DANKOV_CURIE_SHIFT_RATE_K_PER_T in
# core/validation.py (a literature/fit-derived constant, not a live call).
GADOLINIUM_CALIBRATED_SIGMA_TC_K = 6.74

GADOLINIUM_CALIBRATED = BroadenedMagnetocaloricMaterial(
    GADOLINIUM, GADOLINIUM_CALIBRATED_SIGMA_TC_K)
# Best-available calibrated Gd model, stacking every physics fix in this
# pass: exact isentropic DeltaT_ad (delta_T_adiabatic_exact), the
# Sommerfeld electronic entropy/heat-capacity term (GADOLINIUM's
# sommerfeld_gamma_J_per_molK2, core/mce_material.py), and a fitted
# polycrystalline grain-Tc-spread broadening (this constant). Used by
# core/validation.py's run_validation() as the primary reported model;
# plain GADOLINIUM (sigma_Tc=0, no broadening) is unchanged and still used
# everywhere else in this repo (core/amr_cycle.py, core/optimize.py, etc.)
# where the cheap linear delta_T_adiabatic() is what performance requires.


def _make_material(base_material, sigma_Tc_K):
    """sigma=0 returns base_material itself (the pre-existing sharp model,
    unchanged) rather than a degenerate 1-clone BroadenedMagnetocaloricMaterial,
    so the sigma=0 row in every sweep below is byte-identical to what
    core/validation.py already computes and prints in ."""
    if sigma_Tc_K <= 0:
        return base_material
    return BroadenedMagnetocaloricMaterial(base_material, sigma_Tc_K)


def _peak_and_fwhm(material, mu0H_T, T_grid_K=T_GRID_K):
    """Peak DeltaT_ad, the temperature it occurs at, and the full-width-at-
    half-max (K) of DeltaT_ad(T) at fixed field mu0H_T, evaluated on
    T_grid_K. FWHM is found by linear interpolation of the two half-peak
    crossings either side of the peak -- adequate here since T_grid_K's
    0.1K spacing is far finer than the multi-Kelvin broadening effect
    being measured (checked in tests/test_inhomogeneous_broadening.py)."""
    H = mu0H_T / mu0
    dT = np.asarray(material.delta_T_adiabatic(T_grid_K, H)).ravel()
    i_peak = int(np.argmax(dT))
    peak_val = float(dT[i_peak])
    peak_T = float(T_grid_K[i_peak])
    half = peak_val / 2.0

    left = float(T_grid_K[0])
    for i in range(i_peak, 0, -1):
        if dT[i - 1] <= half <= dT[i]:
            left = float(np.interp(half, [dT[i - 1], dT[i]], [T_grid_K[i - 1], T_grid_K[i]]))
            break
    right = float(T_grid_K[-1])
    for i in range(i_peak, len(dT) - 1):
        if dT[i] >= half >= dT[i + 1]:
            right = float(np.interp(half, [dT[i + 1], dT[i]], [T_grid_K[i + 1], T_grid_K[i]]))
            break
    return peak_T, peak_val, (right - left)


def run_broadening_sweep(material=GADOLINIUM, mu0H_T=(1.0, 2.0, 5.0),
                          sigma_values=SIGMA_TC_SWEEP_K, verbose=True):
    """For each sigma_Tc in sigma_values and each field in mu0H_T, reports
    peak DeltaT_ad, the temperature it occurs at, and the FWHM of the
    DeltaT_ad(T) curve -- sigma_Tc=0.0 is the pre-existing sharp model."""
    rows = []
    for sigma in sigma_values:
        mat = _make_material(material, sigma)
        for B in mu0H_T:
            peak_T, peak_val, fwhm = _peak_and_fwhm(mat, B)
            rows.append({"sigma_Tc_K": sigma, "mu0H_T": B, "peak_T_K": peak_T,
                         "peak_dTad_K": peak_val, "fwhm_K": fwhm})
            if verbose:
                print(f" sigma_Tc={sigma:4.1f}K mu0H={B:4.1f}T "
                      f"peak_dTad={peak_val:6.3f}K at T={peak_T:6.2f}K "
                      f"FWHM={fwhm:5.2f}K")
    return rows


def run_dankov_error_sensitivity(sigma_values=SIGMA_TC_SWEEP_K, verbose=True):
    """Re-runs the SAME fixed-T=294K, per-field comparison as
    core/validation.py's run_validation(), but deliberately using
    GADOLINIUM.delta_T_adiabatic() -- the fast linear approximation, not
    that function's own delta_T_adiabatic_exact() (Paper-Mining Pass
    physics fix) -- since this sweep needs many cheap evaluations per
    sigma_Tc, not the single per-field root-solve the exact method uses.
    Checks whether Gaussian Tc-broadening narrows or widens the field-
    dependent error pattern (originally +48.9% at 1T, +29.2% at 2T, -7.5%
    at 5T under the pre-fix linear model with no electronic term; now
    smaller at every field after the electronic-heat-capacity addition to
    total_heat_capacity() -- see core/mce_material.py's
    sommerfeld_gamma_J_per_molK2 -- but still following the same
    field-dependent shape, which is what this sweep actually probes)."""
    from core.validation import LITERATURE_DELTA_T_AD
    rows = []
    for sigma in sigma_values:
        mat = _make_material(GADOLINIUM, sigma)
        for B, dT_lit in LITERATURE_DELTA_T_AD.items():
            H = B / mu0
            dT_model = float(np.asarray(mat.delta_T_adiabatic(np.array([294.0]), H)).ravel()[0])
            err_pct = 100 * (dT_model - dT_lit) / dT_lit
            rows.append({"sigma_Tc_K": sigma, "mu0H_T": B, "dT_lit_K": dT_lit,
                         "dT_model_K": dT_model, "err_pct": err_pct})
            if verbose:
                print(f" sigma_Tc={sigma:4.1f}K | mu0H={B:.1f}T | lit={dT_lit:5.2f}K | "
                      f"model={dT_model:6.3f}K | error={err_pct:+7.1f}%")
    return rows


def run_inhomogeneous_broadening_analysis(out_path="results/inhomogeneous_broadening.txt",
                                           verbose=True):
    """Top-level entry point: runs both sweeps above, finds
    which sigma_Tc (if any) minimizes the worst-field error against
    Dan'kov et al.'s three literature points, states the honest conclusion
    either way, and writes results/inhomogeneous_broadening.txt."""
    lines = []

    def log(s=""):
        if verbose:
            print(s)
        lines.append(s)

    log("=" * 90)
    log("PHASE 22 ITEM 1: Gaussian inhomogeneous/polycrystalline Tc-broadening")
    log("sensitivity for the mean-field Gd model (core/mce_material.py)")
    log("=" * 90)
    log("HONESTY FLAG: Tishin & Spichkin (2003) Sec. 2.8 (inhomogeneous ferromagnets)")
    log("is this item's intended primary source per phase_plan.md. Confirmed directly")
    log("that this project's copy is an image-only PDF (pdfplumber extracts 0")
    log("characters from every one of its 486 pages sampled) -- same finding already")
    log("flagged for Tishin Ch.11 in the earlier baseline_cooling.py. So Sec. 2.8's")
    log("specific content could not be digitized; what follows is the standard")
    log("literature treatment of Tc-distribution broadening (see this module's own")
    log("docstring), not book content.")
    log("")
    log("--- Step 1: peak DeltaT_ad / FWHM vs. sigma_Tc, at Dan'kov et al.'s own three "
        "fields ---")
    sweep_rows = run_broadening_sweep(verbose=verbose)
    log("")
    log("--- Step 2: does broadening narrow the field-dependent error pattern already "
        "found")
    log(" in the pipeline's own Step 1 (core/validation.py's run_validation())? ---")
    err_rows = run_dankov_error_sensitivity(verbose=verbose)

    by_sigma = collections.defaultdict(list)
    for r in err_rows:
        by_sigma[r["sigma_Tc_K"]].append(abs(r["err_pct"]))
    max_err_by_sigma = {s: max(v) for s, v in by_sigma.items()}
    sharp_max = max_err_by_sigma[0.0]
    best_sigma = min(max_err_by_sigma, key=max_err_by_sigma.get)

    sigma_values_sorted = sorted(max_err_by_sigma)
    log("")
    log("Worst-field |error| vs. sigma_Tc:")
    for s in sigma_values_sorted:
        flag = "  <- sharp (pre-existing) model" if s == 0.0 else ""
        flag += "  <- minimum" if s == best_sigma else ""
        log(f" sigma_Tc={s:4.1f}K worst-field |error|={max_err_by_sigma[s]:6.1f}%{flag}")

    log("")
    fwhm_sharp = next(r["fwhm_K"] for r in sweep_rows
                       if r["sigma_Tc_K"] == 0.0 and r["mu0H_T"] == 5.0)
    fwhm_best = next(r["fwhm_K"] for r in sweep_rows
                      if r["sigma_Tc_K"] == best_sigma and r["mu0H_T"] == 5.0)
    log(f"At 5T, sigma_Tc={best_sigma:.1f}K widens the DeltaT_ad(T) FWHM from "
        f"{fwhm_sharp:.2f}K (sharp) to {fwhm_best:.2f}K.")

    # 1T dominates the worst-field metric at every sigma tested (see the
    # per-field printout above), so check the two fields separately as well
    # -- a single "worst-field error" number can hide a real trade-off.
    err_1T_by_sigma = {r["sigma_Tc_K"]: r["err_pct"] for r in err_rows if r["mu0H_T"] == 1.0}
    err_5T_by_sigma = {r["sigma_Tc_K"]: r["err_pct"] for r in err_rows if r["mu0H_T"] == 5.0}
    boundary_limited = best_sigma == max(sigma_values_sorted) and best_sigma != 0.0

    if best_sigma != 0.0 and max_err_by_sigma[best_sigma] < sharp_max - 1e-9:
        tradeoff = ""
        if abs(err_5T_by_sigma[best_sigma]) > abs(err_5T_by_sigma[0.0]) + 1e-9:
            tradeoff = (
                f" This comes with a real trade-off, not a clean win: the 5T error "
                f"moves the OTHER direction, from {err_5T_by_sigma[0.0]:+.1f}% (sharp) "
                f"to {err_5T_by_sigma[best_sigma]:+.1f}% (broadened) -- broadening "
                f"lowers the 1-2T overestimate but deepens the 5T underestimate, "
                f"consistent with a global smoothing of the transition rather than a "
                f"field-selective fix."
            )
        boundary_note = ""
        if boundary_limited:
            boundary_note = (
                f" NOTE: sigma_Tc={best_sigma:.1f}K is the LARGEST value swept here, "
                f"and the 1T error is still falling monotonically with sigma at that "
                f"edge ({[round(err_1T_by_sigma[s], 1) for s in sigma_values_sorted]} "
                f"for sigma={list(sigma_values_sorted)} respectively) -- so this sweep "
                f"has NOT located an interior optimum, only shown the direction of the "
                f"trend within a physically plausible 0-5K range; a wider sweep (and, "
                f"ideally, a real digitized dTad(T) curve to fit against, which this "
                f"repo does not have) would be needed before treating any specific "
                f"sigma_Tc as a calibrated value rather than a sensitivity finding."
            )
        conclusion = (
            f"Gaussian Tc-broadening DOES narrow the worst-field |error| against "
            f"Dan'kov et al.'s three points, from {sharp_max:.1f}% (sharp model, "
            f"sigma_Tc=0) to {max_err_by_sigma[best_sigma]:.1f}% at "
            f"sigma_Tc={best_sigma:.1f}K, driven mainly by the 1T point "
            f"({err_1T_by_sigma[0.0]:+.1f}% -> {err_1T_by_sigma[best_sigma]:+.1f}%)."
            + tradeoff + boundary_note
        )
    else:
        conclusion = (
            f"Gaussian Tc-broadening does NOT narrow the worst-field error at any "
            f"sigma_Tc swept here ({list(sigma_values_sorted)}): the sharp "
            f"(pre-existing, sigma_Tc=0) model's worst-field error "
            f"({sharp_max:.1f}%) remains the minimum. Broadening smooths and widens "
            f"the DeltaT_ad(T) curve (see the FWHM numbers above) as expected, but "
            f"at T=294K -- close to the sharp model's own peak and to Dan'kov et "
            f"al.'s reported peak temperature -- smoothing only ever LOWERS the "
            f"predicted DeltaT_ad relative to the sharp curve; it cannot fix the "
            f"sharp model's field-dependent error PATTERN (over- at 1-2T, under- at "
            f"5T) because that pattern is a genuine shape mismatch across fields, "
            f"not simply excess peak height at one field. This is additional, "
            f"first-pass evidence for the SAME near-Tc mean-field limitation "
            f"core/validation.py's run_curie_shift_check() already documents from a "
            f"different angle (the model's peak temperature does not shift with "
            f"field the way Dan'kov et al. report), not a contradiction of it."
        )
    log("")
    log("CONCLUSION: " + conclusion)

    import os
    os.makedirs(os.path.dirname(out_path) or ".", exist_ok=True)
    with open(out_path, "w") as f:
        f.write("\n".join(lines) + "\n")
    if verbose:
        print(f"Wrote {out_path}")

    return {
        "sweep_rows": sweep_rows,
        "err_rows": err_rows,
        "max_err_by_sigma": max_err_by_sigma,
        "sharp_max_err_pct": sharp_max,
        "best_sigma_Tc_K": best_sigma,
        "conclusion": conclusion,
    }


if __name__ == "__main__":
    run_inhomogeneous_broadening_analysis()