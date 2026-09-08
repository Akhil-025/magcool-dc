"""
hybrid_solid_state_regenerator.py
===================================
NEW CANDIDATE DEVICE ARCHITECTURE (deliberately NOT a material blend --
see the honesty flag below for why that lever was already ruled out).

Implements Lin, Wang, Dai, Qiao, Zhou, Zhao, Hu & Shen, "A full
solid-state conceptual magnetocaloric refrigerator based on hybrid
regeneration", The Innovation (Camb) 5(4):100645 (2024),
doi:10.1016/j.xinn.2024.100645, PMCID PMC11192848 -- found via a
literature search run specifically because the previously-tried lever
(core/nanocomposite_material.py, blending several MCM phases) was
checked and confirmed NOT to raise COP at a design point (only
robustness to off-design operation, per this repo's own README).

WHY THIS IS A DIFFERENT LEVER THAN MATERIAL BLENDING
-------------------------------------------------------
core/nanocomposite_material.py's WeightedMaterialEnsemble mixes several
magnetocaloric PHASES at the DeltaT_ad level; the resulting blend is
mathematically bounded by its best constituent's own peak (a weighted
average cannot exceed its largest term), which is exactly why the
existing tests/README document it underperforming a single tuned phase
at its own design span. That result is a property of averaging entropy
curves, not a limitation of this repo's implementation -- so blending
harder, or blending different families, cannot get around it.

The lever this module implements does not touch the magnetocaloric
material's entropy curve at all. It replaces the heat-transfer FLUID
that this repo's entire existing AMRSystem/cascade framework assumes
(core/amr_cycle.py, core/thermal.py, core/loss_model.py -- all built
around mdot, particle_diameter, pumping power, dead volume) with
alternating solid slices of a high-thermal-conductivity material (HTCM:
Cu, Ag, or an Al-diamond composite) interleaved with the magnetocaloric
material (MCM) itself. Both layers rotate past each other so BOTH
develop their own stable temperature gradient -- a "hybrid magnetic
regenerator" (HMR), as opposed to AMR (only the MCM carries a gradient)
or PMR (only the fluid does). This removes pumping power, viscous/
Darcy-Weisbach losses, and dead-volume losses from the loss budget
ENTIRELY (there is no fluid), replacing them with a different parasitic
channel (mechanical friction between the rotating layers) that the
source paper itself measures and that this module carries forward
rather than dropping.

HONESTY FLAG -- what is and is not first-principles here
-------------------------------------------------------------------
The source paper's headline numbers (its Table 1: exergy efficiency
eta=73.1% at 1 Hz, 55.3% at 5 Hz, 54.2% at 10 Hz, for Gd/Cu slices at
mu0H=1.35 T, Ns=24 sectors, its own optimized utilization factor Uf=1.0)
come from a COMSOL finite-element simulation of coupled transient heat
conduction between rotating slices, calibrated against the paper's own
SQUID-VSM-measured Gd/LaFeSiCoSi (LFS) delta_T_ad and M-H data.
Reproducing that FEA is out of scope here (no mesh, no transient solve
in this repo). What IS implemented, honestly split into two tiers:

  TIER 1 -- first-principles, from the paper's own closed-form equations
  and real material properties (nothing invented):
    - utilization_factor() = Eq. 1 of the source paper.
    - thermal_diffusion_length() = Eq. 2 of the source paper.
    - diffusion_sufficiency_ok() checks the paper's own stated
      sufficiency criterion (slice thickness <= diffusion length) that
      it uses to explain why eta/Rf collapse above ~1 Hz for Gd.

  TIER 2 -- a curve FIT to the paper's own three reported data points
  (its Table 1 "this work" rows, Gd/Cu, Ns=24, Uf=1.0: f in {1, 5,
  10} Hz), NOT a re-derivation of the FEA:
    - eta_hmr(frequency_Hz) and delta_T_span_hmr(frequency_Hz)
      piecewise-linearly interpolate those three points; evaluating
      outside [1, 10] Hz is extrapolation and is flagged in the
      returned result (extrapolated=True), the same discipline this
      repo already applies in
      core/mnfepsi_doped_hysteresis_speculative.py's own >20x-
      extrapolation flag. This is explicitly a fit to three literature
      points, not a physical derivation.

  Friction/air-gap parasitic loss is an OPTIONAL derating, using the
  paper's own two reported points at f=10 Hz (100 nm air gap: -13%
  cooling power, -26% exergy efficiency; ~0% at 10 nm) -- default OFF
  (frictionless limit, matching the paper's own headline Table 1
  numbers), so callers see both the paper's best-case number and a
  friction-derated one explicitly. This derating is ONLY validated at
  f=10 Hz in the source paper; using it at another frequency is flagged.

This module is purely additive: it imports nothing from and is imported
by nothing in core/amr_cycle.py, core/optimize.py, core/cascade.py, or
core/loss_model.py, so it changes no existing result. Run standalone:
    python -m core.hybrid_solid_state_regenerator

INTEGRATION NOTE (main.py step "19."): main.py calls compare_to_vcc()
at this repo's own representative operating point (T_cold=291.15K,
span=REPRESENTATIVE_SPAN_K=10.0K -- the exact same point step 4's own
representative_row is drawn from) as a new, purely additive pipeline
stage. It writes only results/hybrid_solid_state_regenerator.txt and
does not read from or write to any variable any other stage consumes
(no existing stage's numbers change; representative_row itself is only
READ, not modified). See main.py's own stages list for where it's
wired in.
"""

import dataclasses

import numpy as np

from core.baseline_cooling import carnot_cop, vapor_compression_cop

# ---------------------------------------------------------------------------
# Real material properties (not invented). Sources given per-value.
# ---------------------------------------------------------------------------

MATERIAL_PROPERTIES = {
    # Gd: density is this repo's own existing GD_FAMILY.density_kg_m3
    # (core/cascade.py, 7900 kg/m^3, standard literature value for
    # elemental Gd). Specific heat ~236 J/(kg K) and thermal conductivity
    # ~10.5 W/(m K) near room temperature are standard literature values
    # for elemental Gd (Dan'kov, Tishin, Pecharsky & Gschneidner, "Magnetic
    # phase transitions and the magnetothermal properties of gadolinium",
    # Phys. Rev. B 57, 3478 (1998) -- already cited elsewhere in this
    # repo, e.g. main.py's material-level validation step -- report Gd's
    # zero-field heat capacity near Tc in this range).
    "Gd": dict(density_kg_m3=7900.0, cp_J_per_kgK=236.0, kappa_W_per_mK=10.5,
               delta_T_ad_K=4.0,  # @ 1.35 T, source paper's own SQUID-VSM measurement
               source="Dan'kov et al. PRB 57, 3478 (1998) [Cp, kappa]; "
                      "Lin et al., Innovation 5(4):100645 (2024) [delta_T_ad @1.35T]"),
    # LFS = La(Fe0.92Co0.08)11.7Si1.3, the source paper's own second MCM.
    # Density/Cp/kappa taken from the source paper's own Table S1
    # (supplemental information) is NOT re-digitized here (image-only
    # figure in the PDF this pass had access to) -- this repo's own
    # LAFESIH_FAMILY.density_kg_m3 (core/cascade.py, 7300 kg/m^3, a
    # different but closely related La(Fe,Si)13Hy composition) is reused
    # instead as the best available density proxy. Cp/kappa for LFS are
    # therefore NOT populated here (would need the supplemental table);
    # only delta_T_ad (a headline number stated in the paper's own main
    # text, not a supplemental-only figure) is populated.
    "LFS": dict(density_kg_m3=7300.0, cp_J_per_kgK=None, kappa_W_per_mK=None,
                delta_T_ad_K=1.78,  # @ 1.35 T, source paper's own SQUID-VSM measurement
                source="core/cascade.py LAFESIH_FAMILY density (density proxy only); "
                       "Lin et al., Innovation 5(4):100645 (2024) [delta_T_ad @1.35T]"),
    # Cu: standard elemental physical constants (CRC Handbook of
    # Chemistry and Physics; not specific to this paper).
    "Cu": dict(density_kg_m3=8960.0, cp_J_per_kgK=385.0, kappa_W_per_mK=401.0,
               source="CRC Handbook of Chemistry and Physics (standard elemental values)"),
    # Ag: standard elemental physical constants.
    "Ag": dict(density_kg_m3=10490.0, cp_J_per_kgK=235.0, kappa_W_per_mK=429.0,
               source="CRC Handbook of Chemistry and Physics (standard elemental values)"),
}


def utilization_factor(rho_htcm, t_htcm, cp_htcm, rho_mcm, t_mcm, cp_mcm):
    """Eq. 1 of Lin et al. (2024): ratio of the HTCM slice's heat
    capacity (areal, per unit face area) to the MCM slice's own. The
    paper's own finding (checked here only algebraically, not
    re-simulated): Uf=1.0 is the FEA-found optimum, stays optimal across
    1-10 Hz and across all four MCM/HTCM combinations the paper tested
    (Gd/LFS x Cu/Ag/ALC) -- this function only computes Uf for a given
    geometry/material choice; it does not re-derive that 1.0 is optimal."""
    return (rho_htcm * t_htcm * cp_htcm) / (rho_mcm * t_mcm * cp_mcm)


def thermal_diffusion_length(kappa_W_per_mK, rho_kg_m3, cp_J_per_kgK, frequency_Hz):
    """Eq. 2 of Lin et al. (2024): L = sqrt(alpha * delta_t), with
    thermal diffusivity alpha = kappa/(rho*cp) and heat-exchange time
    delta_t = 1/(96*f) (the paper's own stated fraction of one cycle
    period available for MCM<->HTCM heat exchange, from its own
    Ns-sector rotary geometry). Returns L in meters."""
    alpha = kappa_W_per_mK / (rho_kg_m3 * cp_J_per_kgK)
    delta_t = 1.0 / (96.0 * frequency_Hz)
    return float(np.sqrt(alpha * delta_t))


def diffusion_sufficiency_ok(thickness_m, kappa_W_per_mK, rho_kg_m3,
                              cp_J_per_kgK, frequency_Hz):
    """The paper's own stated sufficiency criterion (Results, "Impact of
    working frequency"): heat exchange between MCM and HTCM sectors is
    sufficient only while the slice thickness stays <= the thermal
    diffusion length at that frequency. Returns (ok: bool, L_m: float,
    thickness_over_L: float) so callers can see the margin, not just a
    pass/fail."""
    L = thermal_diffusion_length(kappa_W_per_mK, rho_kg_m3, cp_J_per_kgK, frequency_Hz)
    return thickness_m <= L, L, thickness_m / L


# ---------------------------------------------------------------------------
# TIER 2: curve fit to the paper's own three reported (frequency, eta,
# delta_T_span) points -- Gd/Cu, Ns=24, Uf=1.0 (its Table 1's own "this
# work" rows). NOT re-derived from FEA -- see module docstring.
# ---------------------------------------------------------------------------

_HMR_GDCU_FREQ_HZ = np.array([1.0, 5.0, 10.0])
_HMR_GDCU_ETA = np.array([0.731, 0.553, 0.542])            # exergy efficiency (COP/COP_carnot)
_HMR_GDCU_SPAN_K = np.array([47.0, 38.0, 26.0])             # no-load max span, K
_HMR_GDCU_SCP_KW_PER_KG = np.array([0.9, 4.3, 8.3])         # specific cooling power, kW/kg

# Single extra point at 10 Hz for LFS/Cu (0.1 mm LFS / 0.13 mm Cu slices),
# the paper's own reported alternative-MCM case -- used only as a spot
# comparison, not fit into the frequency curve above (only one frequency
# is reported for this combination).
_HMR_LFSCU_10HZ = dict(eta=0.533, span_K=14.0, scp_kW_per_kg=9.8)


@dataclasses.dataclass
class HMRResult:
    technology: str
    frequency_Hz: float
    eta: float                  # exergy efficiency, COP/COP_carnot
    COP: float
    COP_carnot: float
    delta_T_span_max_K: float   # this frequency's own no-load max span (informational)
    extrapolated: bool          # True if frequency_Hz fell outside the fitted [1,10] Hz range
    source_note: str


def eta_hmr(frequency_Hz):
    """Interpolates/extrapolates the source paper's own three (f, eta)
    points (Gd/Cu, Ns=24, Uf=1.0). Returns (eta, extrapolated: bool).
    Linear interpolation within [1, 10] Hz; linear extrapolation from the
    nearest edge slope outside it (flagged via the returned bool -- the
    paper's own data shows eta declining roughly monotonically with f in
    this range, but nothing supports the shape of that decline beyond
    the fitted points, particularly below 1 Hz where AMR/PMR devices
    conventionally operate and this HMR was not evaluated)."""
    extrapolated = not (_HMR_GDCU_FREQ_HZ[0] <= frequency_Hz <= _HMR_GDCU_FREQ_HZ[-1])
    eta = float(np.interp(frequency_Hz, _HMR_GDCU_FREQ_HZ, _HMR_GDCU_ETA))
    return eta, extrapolated


def delta_T_span_max_hmr(frequency_Hz):
    extrapolated = not (_HMR_GDCU_FREQ_HZ[0] <= frequency_Hz <= _HMR_GDCU_FREQ_HZ[-1])
    span = float(np.interp(frequency_Hz, _HMR_GDCU_FREQ_HZ, _HMR_GDCU_SPAN_K))
    return span, extrapolated


# Paper's own two reported friction/air-gap points at f=10 Hz (Table 1's
# air-gap rows + the main text's separate measured-friction paragraph).
# ONLY validated at 10 Hz -- applying elsewhere is an assumption, flagged
# by the function below.
_AIR_GAP_DERATING_10HZ = {
    0: dict(cooling_power_frac=1.0, eta_frac=1.0),      # no air gap (thermal grease contact)
    10: dict(cooling_power_frac=1.0, eta_frac=1.0),     # 10 nm gap: "nearly unchanged"
    100: dict(cooling_power_frac=1.0 - 0.13, eta_frac=1.0 - 0.26),  # 100 nm gap
}


def apply_air_gap_derating(cop, eta, air_gap_nm, frequency_Hz):
    """Derates a headline (frictionless-limit) COP/eta using the source
    paper's own two measured/simulated air-gap points, both reported
    ONLY at f=10 Hz. air_gap_nm must be one of {0, 10, 100} (the paper's
    own three tested values) -- this function does not interpolate
    between them (no intermediate point is reported) or extrapolate
    beyond 100 nm. Returns (cop_derated, eta_derated, note)."""
    if air_gap_nm not in _AIR_GAP_DERATING_10HZ:
        raise ValueError(
            f"air_gap_nm={air_gap_nm} not one of the paper's own tested "
            f"values {sorted(_AIR_GAP_DERATING_10HZ)} -- no derating factor "
            f"available (not interpolated/extrapolated; would be invented)."
        )
    d = _AIR_GAP_DERATING_10HZ[air_gap_nm]
    note = (f"Air-gap derating ({air_gap_nm} nm) applied from the source paper's "
            f"own f=10 Hz measurement/simulation")
    if not np.isclose(frequency_Hz, 10.0):
        note += (f"; APPLIED HERE AT f={frequency_Hz:.2f} Hz, which the paper did NOT "
                  f"test -- this is an unvalidated assumption of frequency-independence, "
                  f"not a reported result.")
    return cop * d["cooling_power_frac"], eta * d["eta_frac"], note


def hmr_cop(T_cold_K, T_hot_K, frequency_Hz=1.0, air_gap_nm=None):
    """COP of the Gd/Cu HMR device at (T_cold_K, T_hot_K), via the
    source paper's own Eq. 4 (eta = COP/COP_carnot) rearranged as
    COP = eta(f) * COP_carnot(Tc, Th) -- i.e. this function does NOT
    independently recompute COP from Qc/Wm (Eq. 3); it uses the paper's
    own already-computed eta(f) (TIER 2 fit above) applied to THIS
    caller's own (Tc, Th), exactly the same "static efficiency fraction
    of Carnot" pattern this repo's own vapor_compression_cop()/
    liquid_cooling_cop() already use for their conventional baselines --
    so the comparison is apples-to-apples in structure, not just in units.

    air_gap_nm: optional, one of {0, 10, 100} -- see
    apply_air_gap_derating(); None (default) = frictionless headline
    number, matching the paper's own Table 1."""
    cc = carnot_cop(T_cold_K, T_hot_K)
    eta, extrapolated = eta_hmr(frequency_Hz)
    span_max, _ = delta_T_span_max_hmr(frequency_Hz)
    cop = eta * cc
    note = ("Lin et al., Innovation 5(4):100645 (2024), Gd/Cu HMR, Ns=24, Uf=1.0, "
            "eta(f) fit to 3 reported points (1,5,10 Hz)")
    if extrapolated:
        note += f" -- EXTRAPOLATED to f={frequency_Hz:.2f} Hz (outside fitted [1,10] Hz range)"
    if air_gap_nm is not None:
        cop, eta, air_note = apply_air_gap_derating(cop, eta, air_gap_nm, frequency_Hz)
        note += "; " + air_note
    return HMRResult(
        technology="Hybrid solid-state magnetic regenerator (HMR, Gd/Cu)",
        frequency_Hz=frequency_Hz, eta=eta, COP=cop, COP_carnot=cc,
        delta_T_span_max_K=span_max, extrapolated=extrapolated, source_note=note,
    )


def compare_to_vcc(T_cold_K, span_K, frequencies_Hz=(1.0, 5.0, 10.0),
                    eta_2nd_law_vcc=0.42, out_path="results/hybrid_solid_state_regenerator.txt",
                    verbose=True):
    """Runs the actual comparison this module exists to answer: does the
    HMR's COP, at the source paper's own reported operating points, close
    any of the gap to vapor-compression cooling at this repo's own
    ASHRAE-representative (T_cold, span)? Computed honestly -- if the
    answer is no, this function reports that, not a massaged version."""
    T_hot_K = T_cold_K + span_K
    vcc = vapor_compression_cop(T_cold_K, T_hot_K, eta_2nd_law=eta_2nd_law_vcc)

    lines = []

    def log(s=""):
        if verbose:
            print(s)
        lines.append(s)

    log("=" * 88)
    log("Hybrid solid-state magnetic regenerator (HMR) vs. vapor-compression, "
        f"T_cold={T_cold_K:.2f}K, span={span_K:.1f}K (T_hot={T_hot_K:.2f}K)")
    log("Source: Lin et al., Innovation (Camb) 5(4):100645 (2024), Gd/Cu, Ns=24, Uf=1.0")
    log("=" * 88)
    log(f"{'freq(Hz)':>9}{'eta_HMR':>10}{'COP_HMR':>10}{'COP_carnot':>12}"
        f"{'COP_VCC':>10}{'HMR/VCC':>10}{'span_max_K':>12}{'flag':>14}")

    rows = []
    for f in frequencies_Hz:
        r = hmr_cop(T_cold_K, T_hot_K, frequency_Hz=f)
        ratio = r.COP / vcc.COP if vcc.COP > 0 else float("nan")
        flag = "EXTRAPOLATED" if r.extrapolated else ""
        span_note = (" (device's own max span < requested span -- "
                      "infeasible at this span, informational only)"
                      if r.delta_T_span_max_K < span_K else "")
        log(f"{f:>9.1f}{r.eta:>10.3f}{r.COP:>10.3f}{r.COP_carnot:>12.3f}"
            f"{vcc.COP:>10.3f}{ratio:>10.3f}{r.delta_T_span_max_K:>12.1f}{flag:>14}{span_note}")
        rows.append(dict(frequency_Hz=f, eta_HMR=r.eta, COP_HMR=r.COP,
                          COP_carnot=r.COP_carnot, COP_VCC=vcc.COP,
                          HMR_over_VCC=ratio, span_max_K=r.delta_T_span_max_K,
                          extrapolated=r.extrapolated, span_feasible=r.delta_T_span_max_K >= span_K))

    best = max(rows, key=lambda r: r["COP_HMR"])
    feasible_rows = [r for r in rows if r["span_feasible"]]

    log("")
    log("*** APPLES-TO-ORANGES WARNING (read before trusting the ratio above) ***")
    log(
        "The source paper's own COP = Qc/Wm (its Eq. 3) is MAGNETIC-CYCLE-WORK ONLY -- no "
        "drive-motor power to rotate the layers, no permanent-magnet-system overhead, no "
        "controls/sensing power. That is the exact same quantity this repo's own "
        "AMRCycleResult.COP field holds (Qc/W_mag, its docstring literally calls it "
        "'ideal magnetic-cycle-only COP'), as distinct from AMRCycleResult.COP_electrical "
        "(Qc/(W_mag+W_parasitic), which is what this repo's own main.py step 4 actually "
        "compares against VCC). vapor_compression_cop()'s eta_2nd_law=0.42, by contrast, IS "
        "a real INSTALLED-SYSTEM second-law efficiency (Ebrahimi et al. 2014's data-center "
        "review) -- it already includes the compressor motor's own electrical losses. So the "
        "HMR/VCC ratio computed above compares an IDEAL, parasitic-free numerator against a "
        "REAL, parasitic-INCLUDED denominator -- structurally the same mismatch this repo's "
        "own diagnose_cop_gap.py exists to catch elsewhere (there, swapping a naive parasitic "
        "assumption for a calibrated one moved COP by ~2x with the identical W_mag). The "
        "source paper reports no drive-motor/magnet-system electrical draw to correct this "
        "with (it explicitly scopes COP to Wm only) -- so this module's own COP_HMR should be "
        "read as an UPPER BOUND on what a real HMR device's electrical COP could achieve, not "
        "an apples-to-apples electrical-COP figure. Contains one real, HMR-specific parasitic "
        "channel not present in fluid AMR at all: inter-layer mechanical friction (see "
        "apply_air_gap_derating(), validated only at 10 Hz: a 100 nm air gap alone costs 26% "
        "of exergy efficiency -- friction/contact losses in this all-solid design are not "
        "necessarily smaller than fluid AMR's pumping losses, just a DIFFERENT channel)."
    )
    log("")
    if best["HMR_over_VCC"] >= 1.0:
        conclusion = (
            f"At its best reported operating point ({best['frequency_Hz']:.1f} Hz), the HMR's "
            f"own eta={best['eta_HMR']:.3f}-implied IDEAL COP ({best['COP_HMR']:.2f}) exceeds "
            f"this repo's VCC ELECTRICAL COP baseline ({vcc.COP:.2f}) at this (T_cold, span) -- "
            f"ratio {best['HMR_over_VCC']:.2f}x -- but per the warning above, this compares an "
            f"ideal (motor/magnet-overhead-free) number against a real installed one, so it is "
            f"NOT evidence the gap is actually closed; it says only that the underlying "
            f"magnetic-cycle thermodynamics of this architecture leave enough headroom (73.1% "
            f"vs. this repo's Gd AMR ideal-COP numbers -- see README's own reported ideal-vs-"
            f"electrical COP gap) that closing the real electrical gap is not obviously "
            f"foreclosed the way it would be if even the IDEAL number already trailed VCC's "
            f"REAL number. This is the headline, frictionless-limit number from a single 2024 "
            f"conceptual/simulation paper (COMSOL FEA + measured M-H data, a sub-millimeter "
            f"lab-scale conceptual device -- not yet a built full-scale prototype)."
        )
    else:
        conclusion = (
            f"Even at its best reported operating point ({best['frequency_Hz']:.1f} Hz), and "
            f"even comparing this architecture's IDEAL (parasitic-free) COP against VCC's REAL "
            f"electrical COP -- the more favorable of the two comparisons for HMR -- the HMR's "
            f"own eta={best['eta_HMR']:.3f}-implied COP ({best['COP_HMR']:.2f}) does NOT exceed "
            f"this repo's VCC baseline ({vcc.COP:.2f}) at this (T_cold, span) -- ratio "
            f"{best['HMR_over_VCC']:.2f}x. Reported honestly, not massaged: this architecture "
            f"would not close the gap even before accounting for its own real-world losses."
        )
    log("CONCLUSION: " + conclusion)
    if not feasible_rows:
        log(f"CAVEAT: none of the tested frequencies' own no-load max span reaches the "
            f"requested span_K={span_K:.1f}K at Ns=24 -- the device as reported would need a "
            f"larger Ns (more sectors) or a lower frequency than tested to physically reach "
            f"this span at all; the COP numbers above are informational, not span-feasible.")
    log("")
    log("Compare against this repo's OWN fluid-based AMR baseline at the same point via "
        "main.py step 4 (results/comparison_table.csv) -- not recomputed here to avoid "
        "duplicating that pipeline; see this module's own docstring for how to read the two "
        "results together (frictionless-limit conceptual device vs. this repo's calibrated, "
        "loss-model-derated existing AMR numbers).")

    import os
    os.makedirs(os.path.dirname(out_path) or ".", exist_ok=True)
    with open(out_path, "w") as f_out:
        f_out.write("\n".join(lines) + "\n")
    if verbose:
        print(f"\nWrote {out_path}")

    return {"rows": rows, "vcc": vcc, "conclusion": conclusion}


if __name__ == "__main__":
    # Repo's own representative ASHRAE point (T_cold=18C=291.15K) swept
    # across its own 5-20K ASHRAE span range, same range main.py step 4 uses.
    T_COLD_K = 291.15
    for span in (5.0, 10.0, 15.0, 20.0):
        compare_to_vcc(T_COLD_K, span)
        print()
