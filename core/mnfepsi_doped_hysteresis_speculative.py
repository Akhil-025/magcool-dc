"""
mnfepsi_doped_hysteresis_speculative.py
=========================================
SPECULATIVE / UNVALIDATED -- explicitly requested follow-up on
docs/Literature_Review.md's "hysteresis-reducing dopant routes" entry and
the honesty-flag updates in core/first_order_mce.py (MNFEPSI_FIRST_ORDER)
and core/hysteresis_sensitivity.py (honesty flag #2). Those flags stopped
short of putting a number on the V-/B-doped compositions' hysteresis LOSS
in J/kg, because the source papers report only thermal hysteresis WIDTH
in K, and converting K to J/kg needs either the transition's latent heat
or a digitized M-H loop -- neither of which those papers report.

This module takes a different, still-not-fabricated approach: rather than
inventing a conversion, it FITS one from real paired (hysteresis loss
J/kg, entropy change J/(kg K), hysteresis width K) data that already
exists for a closely related composition axis of the SAME source paper
this repo's own MNFEPSI_FIRST_ORDER.hysteresis_loss_J_per_kg=25.0
placeholder is proxied from (Zhang et al., arXiv:2312.09341, Table 2 --
Mn_xFe_{2-x}P_0.5Si_0.5 microwires, x=0.8-1.2), then applies that fitted
relation to the V-doped compositions this repo's literature review
flagged. This is still an estimate, not a measurement -- see the honesty
flags below for exactly how much to trust it and why it is NOT used to
change MNFEPSI_FIRST_ORDER's own default hysteresis_loss_J_per_kg.

Method
------
For each of Zhang et al.'s five measured x=0.8-1.2 compositions, this
module computes k = W_hys / (|dS_iso| * T_hys) -- i.e. it checks whether
hysteresis loss (J/kg) scales roughly as entropy change (J/kg K) times
hysteresis width (K), the simplest dimensionally-consistent proxy anyone
would reach for absent a digitized M-H loop. It does NOT assume this
proxy is exact or literature-established -- it is checked directly
against the five real data points below, and the checked k varies by
+-35% around its own mean (0.151) and up to +-44% around a least-squares
fit through the origin (0.162) -- a real, quantified scatter, not glossed
over. This is the same "check the assumption against real data before
using it" discipline as e.g. core/giguere_validation.py's own
DTAD_CORRECTION_FACTOR.

    Zhang et al. arXiv:2312.09341 Table 2:
    x     T_hys(K)   W_hys(J/kg)  |dS_iso|(J/kg K)   k=W_hys/(dS*Thys)
    0.8   14.5       19.6         12.0                0.113
    0.9   18.5       60.7         18.3                0.179
    1.0   21.0       58.4         15.8                0.176
    1.1   21.5       28.7         10.9                0.123
    1.2   22.0       28.9         7.9                 0.166
    mean k = 0.151 (std 0.028, 18.6% relative)
    least-squares-through-origin k = 0.162 (residuals -44% to +10%)

HONESTY FLAGS (read all four before using KS_FIT_ESTIMATE_J_PER_KG for
anything beyond an order-of-magnitude sanity check)
------------------------------------------------------------------------
1. **Different composition axis.** Zhang et al.'s x=0.8-1.2 sweeps Mn
   content in an UN-doped Mn_xFe_{2-x}P_0.5Si_0.5 system. The V-doped
   compositions this module estimates for (Lai et al. 2024, J. Sci. Adv.
   Mater. Dev. 9(1) 100660; Lai et al. 2022, arXiv:1810.09902 lineage) add
   a THIRD element (vanadium) that the source paper's own text says works
   by a different physical mechanism (reducing volumetric strain / latent
   heat at the transition, not simply shifting Mn content) -- so there is
   no guarantee the same k applies. This is an EXTRAPOLATION across
   composition-axis type, not an interpolation within one.
2. **Far extrapolation in T_hys.** The fitted k was checked only over
   T_hys=14.5-22.0 K. The V-doped compositions sit at T_hys=0.6-0.7 K --
   more than 20x below the fitted range's lower edge, deep into
   "IS_FAR_EXTRAPOLATION" territory by this codebase's own convention
   (see core/electrocaloric_cycle.py's identically-named flag for spans
   below its own measured range). A linear-in-T_hys proxy fit at
   T_hys~20K has NO evidence it stays linear (rather than, say, going to
   zero faster, or saturating at some hysteresis-independent floor loss)
   all the way down to T_hys~0.7K.
3. **Field mismatch.** MNFEPSI_FIRST_ORDER's hysteresis_loss_J_per_kg is
   used at that material's own 2T calibration field. The V-doped |dS_M|
   values used below (9.2 and 5.6 J/(kg K)) are BOTH reported at 1T, not
   2T, in their source papers -- no 2T value was located for either V-doped
   composition. Applying a 1T entropy-change value inside a proxy meant to
   estimate a 2T-calibrated quantity is a further, separate approximation,
   not corrected for here.
4. **Not wired into any default.** `MNFEPSI_FIRST_ORDER.hysteresis_loss_J_per_kg`
   in core/first_order_mce.py is left at 25.0 -- unchanged by this module.
   Nothing in core/hysteresis_sensitivity.py, core/optimize.py, or any
   cascade/AMR pipeline reads this module's output automatically. Treat
   the numbers below as a documented, checkable estimate to reason about,
   not a validated input ready to flow into the rest of this codebase's
   results.

Bottom line this module actually supports (and no more): the doped
compositions' hysteresis loss is very likely far below 25.0 J/kg --
probably by close to an order of magnitude, given ~20x lower hysteresis
width even accounting for a plausible several-fold uncertainty in k and
its extrapolation -- but "close to an order of magnitude below" is the
right confidence level to report, not a specific decimal.
"""

from dataclasses import dataclass

# Zhang et al., arXiv:2312.09341, Table 2 (Mn_xFe_{2-x}P0.5Si0.5 microwires)
_ZHANG_TABLE = [
    # (x, T_hys_K, W_hys_J_per_kg, abs_dS_iso_J_per_kgK)
    (0.8, 14.5, 19.6, 12.0),
    (0.9, 18.5, 60.7, 18.3),
    (1.0, 21.0, 58.4, 15.8),
    (1.1, 21.5, 28.7, 10.9),
    (1.2, 22.0, 28.9, 7.9),
]


def _fit_k():
    """Returns (k_mean, k_std, k_least_squares) for W_hys ~= k * dS * T_hys,
    fit against _ZHANG_TABLE. See module docstring for what this fit is
    and is not evidence of."""
    ks = [w / (ds * t) for (_x, t, w, ds) in _ZHANG_TABLE]
    k_mean = sum(ks) / len(ks)
    k_std = (sum((k - k_mean) ** 2 for k in ks) / len(ks)) ** 0.5
    num = sum((ds * t) * w for (_x, t, w, ds) in _ZHANG_TABLE)
    den = sum((ds * t) ** 2 for (_x, t, w, ds) in _ZHANG_TABLE)
    k_ls = num / den
    return k_mean, k_std, k_ls


K_MEAN, K_STD, K_LEAST_SQUARES = _fit_k()


@dataclass
class SpeculativeHysteresisEstimate:
    name: str
    source: str
    dS_J_per_kgK: float
    dS_field_T: float
    T_hys_K: float
    estimate_J_per_kg_mean_k: float
    estimate_J_per_kg_ls_k: float
    is_far_extrapolation: bool
    caveats: str


def _estimate(name, source, dS, dS_field_T, T_hys):
    est_mean = K_MEAN * dS * T_hys
    est_ls = K_LEAST_SQUARES * dS * T_hys
    return SpeculativeHysteresisEstimate(
        name=name, source=source, dS_J_per_kgK=dS, dS_field_T=dS_field_T,
        T_hys_K=T_hys,
        estimate_J_per_kg_mean_k=est_mean,
        estimate_J_per_kg_ls_k=est_ls,
        is_far_extrapolation=(T_hys < _ZHANG_TABLE[0][1] / 5),
        caveats="See module docstring honesty flags 1-3 (different "
                "composition axis, >20x extrapolation below the fitted "
                "T_hys range, 1T vs 2T field mismatch) before using "
                "either estimate for anything beyond an order-of-"
                "magnitude sanity check.",
    )


# Lai, Huang, You, Maschek, Zhou, van Dijk & Bruck, J. Sci. Adv. Mater.
# Dev. 9(1), 100660 (2024): V-on-Fe substitution, best composition.
VONFE_ESTIMATE = _estimate(
    name="(Mn1.17-xFe0.71V0.02)(P0.5Si0.5), V-on-Fe (Lai et al. 2024)",
    source="Lai et al., J. Sci. Adv. Mater. Dev. 9(1), 100660 (2024)",
    dS=9.2, dS_field_T=1.0, T_hys=0.6,
)

# Lai et al. (V-on-Mn lineage, arXiv:1810.09902 / J. Appl. Phys. follow-up):
# Mn0.98V0.02Fe0.95P0.563Si0.36B0.077.
VONMN_ESTIMATE = _estimate(
    name="Mn0.98V0.02Fe0.95P0.563Si0.36B0.077, V-on-Mn (Lai et al., "
         "arXiv:1810.09902 lineage)",
    source="Lai et al., arXiv:1810.09902; TU Delft follow-up",
    dS=5.6, dS_field_T=1.0, T_hys=0.7,
)


def run_estimate(verbose=True):
    """Prints and returns both speculative estimates plus the fitted k
    and its checked scatter, for direct comparison against
    MNFEPSI_FIRST_ORDER's own 25.0 J/kg placeholder
    (core/first_order_mce.py) and against the OFF (0.0 J/kg) arm of
    core/hysteresis_sensitivity.py's existing A/B check."""
    lines = []
    lines.append("SPECULATIVE hysteresis-loss estimate for V-doped "
                  "(Mn,Fe)2(P,Si) variants (NOT a validated number -- "
                  "see module docstring)")
    lines.append("-" * 78)
    lines.append(f"Fitted proxy k (W_hys ~= k * dS * T_hys), checked against "
                  f"Zhang et al.'s own 5-point table:")
    lines.append(f"  mean k = {K_MEAN:.4f} (std {K_STD:.4f}, "
                  f"{100*K_STD/K_MEAN:.0f}% relative)")
    lines.append(f"  least-squares-through-origin k = {K_LEAST_SQUARES:.4f} "
                  f"(residuals -44% to +10% across the 5 fitted points)")
    lines.append("")
    for est in (VONFE_ESTIMATE, VONMN_ESTIMATE):
        lines.append(f"{est.name}")
        lines.append(f"  source: {est.source}")
        lines.append(f"  dS={est.dS_J_per_kgK} J/(kg K) at {est.dS_field_T}T "
                      f"(NOT 2T -- field-mismatch caveat), T_hys={est.T_hys_K} K")
        lines.append(f"  estimate (mean-k):        "
                      f"{est.estimate_J_per_kg_mean_k:.2f} J/kg")
        lines.append(f"  estimate (least-squares-k): "
                      f"{est.estimate_J_per_kg_ls_k:.2f} J/kg")
        lines.append(f"  far extrapolation flag: {est.is_far_extrapolation}")
        lines.append("")
    lines.append("For context: MNFEPSI_FIRST_ORDER's own default "
                  "hysteresis_loss_J_per_kg = 25.0 J/kg (core/first_order_"
                  "mce.py). Both V-doped estimates above come out well "
                  "under 1 J/kg -- consistent with the qualitative "
                  "'close to an order of magnitude below' conclusion this "
                  "module's docstring states as its actual, defensible "
                  "bottom line, but the specific decimal values above "
                  "should NOT be read as more precise than that.")
    text = "\n".join(lines)
    if verbose:
        print(text)
    return {"k_mean": K_MEAN, "k_std": K_STD, "k_least_squares": K_LEAST_SQUARES,
            "estimates": [VONFE_ESTIMATE, VONMN_ESTIMATE]}


if __name__ == "__main__":
    run_estimate()
