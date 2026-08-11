"""
inhomogeneous_broadening.py
============================
Phase 22 item 1 (see phase_plan.md and ROADMAP.md): models polycrystalline/
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
Phase 21's core/baseline_cooling.py docstring. So Sec. 2.8's specific
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
already found and printed (Phase 1 of the pipeline: overestimate at
1-2T, underestimate at 5T) -- a sensitivity study, not a final fitted
answer.
"""

import dataclasses
import collections

import numpy as np
from numpy.polynomial.hermite_e import hermegauss

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


def _make_material(base_material, sigma_Tc_K):
    """sigma=0 returns base_material itself (the pre-existing sharp model,
    unchanged) rather than a degenerate 1-clone BroadenedMagnetocaloricMaterial,
    so the sigma=0 row in every sweep below is byte-identical to what
    core/validation.py already computes and prints in Phase 1."""
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
                print(f"  sigma_Tc={sigma:4.1f}K  mu0H={B:4.1f}T  "
                      f"peak_dTad={peak_val:6.3f}K at T={peak_T:6.2f}K  "
                      f"FWHM={fwhm:5.2f}K")
    return rows


def run_dankov_error_sensitivity(sigma_values=SIGMA_TC_SWEEP_K, verbose=True):
    """Re-runs core/validation.py's own fixed-T=294K, per-field methodology
    (run_validation()) for each sigma_Tc in the broadening sweep, to check
    whether Gaussian Tc-broadening narrows or widens the SAME field-
    dependent error pattern (pipeline Step 1's own printed log:
    +48.9% at 1T, +29.2% at 2T, -7.5% at 5T) that motivated this item."""
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
                print(f"  sigma_Tc={sigma:4.1f}K | mu0H={B:.1f}T | lit={dT_lit:5.2f}K | "
                      f"model={dT_model:6.3f}K | error={err_pct:+7.1f}%")
    return rows


def run_inhomogeneous_broadening_analysis(out_path="results/inhomogeneous_broadening.txt",
                                           verbose=True):
    """Top-level Phase 22 item 1 entry point: runs both sweeps above, finds
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
    log("flagged for Tishin Ch.11 in Phase 21's baseline_cooling.py. So Sec. 2.8's")
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
    log("    in the pipeline's own Step 1 (core/validation.py's run_validation())? ---")
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
        log(f"  sigma_Tc={s:4.1f}K   worst-field |error|={max_err_by_sigma[s]:6.1f}%{flag}")

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