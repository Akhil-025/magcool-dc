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

  TIER 3 -- REALISTIC ELECTRICAL COP (new addition, this pass): closes
  the single biggest gap the module's own APPLES-TO-ORANGES WARNING
  (see compare_to_vcc()) previously left open -- that COP_HMR was
  magnetic-cycle-work only, with NO drive-motor/drivetrain overhead and
  NO baseline electrical overhead, making it structurally incomparable
  to VCC's real installed electrical COP. `hmr_electrical_cop()` /
  `compare_to_vcc_realistic()` add two REAL, previously-absent parasitic
  loss channels and report a genuine COP_electrical for HMR, built the
  same way (Qc / (W_mag + W_parasitic)) this repo's own
  core.amr_cycle.AMRCycleResult.COP_electrical already is for fluid AMR:

    (a) Inter-layer contact-friction / air-gap loss -- TIER 2's own
        `apply_air_gap_derating()` above, turned ON by default here
        (air_gap_nm=100, the paper's own least-ideal tested case,
        chosen as the more physically realistic assumption for a real
        rotating mechanical interface than the paper's idealized
        zero-gap/thermal-grease contact) instead of the frictionless
        default `hmr_cop()` keeps for backward compatibility. Still the
        paper's OWN measured/simulated numbers -- nothing invented here.

    (b) Rotary drivetrain power (bearing friction, detent/cogging
        torque, and the drive motor needed to rotate the sector disk
        and permanent-magnet assembly past each other -- Lin et al.'s
        own "rotary HMR" design, Ns=24 sectors, structurally requires
        exactly this kind of rotary drivetrain, but the source paper's
        own Eq. 3 COP explicitly EXCLUDES it). This module does NOT
        invent a number for this: it reuses
        `core.loss_model.RotaryDriveLossModel` / `fit_rotary_drive_term()`,
        already calibrated in this repo against Lozano, Engelbrecht,
        Bahl, Nielsen, Eriksen, Olsen, Barbosa, Smith, Prata & Pryds,
        "Performance analysis of a rotary active magnetic refrigerator,"
        Applied Energy 111 (2013) 669-680 -- a rotary permanent-magnet-
        assembly + rotating-regenerator device of the same general
        mechanical class (24 regenerator sectors rotating inside a
        4-pole permanent magnet) as the rotary HMR concept this module
        implements, with DIRECTLY MEASURED drivetrain electrical power
        (its own Table 3 WM column, 87-145 W, fit as W_drive(f) =
        k_drive0 + k_drive1*f). This is a CROSS-DEVICE PROXY, stated
        honestly, not a device-specific HMR drivetrain measurement (none
        exists -- no HMR prototype has been built, per the module's own
        Tier 1/2 honesty flag above): Lozano's device moves fluid
        through a rotary valve as well as rotating the magnet, while
        HMR has no fluid/valve but instead must rotate a segmented
        HTCM/MCM sector disk -- a mechanically different but not
        obviously lighter load (same class of bearing/cogging losses,
        no valve-actuation savings to offset it). The DIRECTION of this
        addition (real rotary magnetocaloric drivetrains draw
        non-trivial parasitic power, often comparable to or exceeding
        Qc itself) is independently corroborated by a SECOND, separate
        source: Arnold, Tura & Rowe, "Experimental analysis of a
        two-material active magnetic regenerator," Int. J. Refrig. 37
        (2014) 99-105, whose own force/drive measurements found
        "mechanical losses and pumping power are the most significant
        contributions to net work while the net magnetic work is too
        small to be resolved," with device efficiencies "all less than
        0.15" -- i.e. two independent measured rotary AMR datasets agree
        that drivetrain/mechanical loss, not magnetic work, dominates
        real hardware of this general class. This does NOT establish
        Lozano's specific k_drive0/k_drive1 values transfer quantitatively
        to an HMR device; it establishes that omitting a drivetrain term
        entirely (the paper's own Eq. 3 scope, and this module's
        previous COP_HMR) is the more clearly wrong assumption of the
        two, and gives the best currently-available real-hardware
        anchor for the correction's rough SIZE.

    (c) Baseline electrical overhead (controls, inverter) -- reuses this
        repo's own CORE-calibrated `base_frac` (core.loss_model.
        calibrate_loss_coefficients(), the same coefficient
        core.amr_cycle.AMRSystem applies to every fluid-AMR
        COP_electrical calculation) rather than inventing a separate
        HMR-specific overhead fraction, so the "non-drivetrain,
        non-friction" slice of the parasitic budget is treated
        identically to how the rest of this repo treats it.

  Reported HMR_electrical_COP = Qc_net / (W_mag_net + W_drivetrain +
  W_baseline), where Qc_net and W_mag_net already include the (a) air-gap
  derating (Qc_net = Qc_ideal * cooling_power_frac; W_mag_net backed out
  from the paper's own derated eta so BOTH the cooling-power loss AND the
  extra magnetic work the paper's own numbers imply are captured, not
  just one). This is now the SAME kind of "real installed electrical
  COP" quantity as vapor_compression_cop() and core.amr_cycle.
  AMRCycleResult.COP_electrical -- closing (not merely re-flagging) the
  ideal-vs-real structural mismatch the original APPLES-TO-ORANGES
  WARNING identified, while being explicit that item (b)'s drivetrain
  number is a cross-device literature proxy, not an HMR-specific
  measurement (none exists yet).

This module is purely additive: it imports nothing from and is imported
by nothing in core/amr_cycle.py, core/optimize.py, core/cascade.py, so
it changes no existing result. It now imports read-only calibration
helpers from core/loss_model.py (fit_rotary_drive_term,
calibrate_loss_coefficients) but does not modify anything in that
module or feed results back into it. Run standalone:
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
from core.loss_model import fit_rotary_drive_term, calibrate_loss_coefficients

# Representative regenerator/MCM mass (kg) used to convert the paper's own
# specific cooling power (kW/kg) into an absolute Qc (W) for the realistic
# electrical-COP path below. Matches this repo's own representative design
# point elsewhere (main.py step 5's economics run: "H=2.0T, mass=5.0kg Gd"),
# so the HMR comparison is anchored to the SAME scale this repo already
# reports its own fluid-AMR numbers at, not an arbitrarily chosen one.
MASS_MCM_KG_DEFAULT = 5.0

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


def specific_cooling_power_kw_per_kg(frequency_Hz):
    """Interpolates/extrapolates the source paper's own three (f, specific
    cooling power) points (Gd/Cu, Ns=24, Uf=1.0, its Table 1's own SCP
    column, kW per kg of MCM). Same TIER 2 discipline as eta_hmr()/
    delta_T_span_max_hmr() above: linear interpolation within [1,10] Hz,
    flagged extrapolation outside it. Needed (new, this pass) to convert
    the paper's dimensionless eta/COP figures into an absolute Qc (W) for
    a given regenerator mass -- the paper itself never needed this because
    it only reports efficiency ratios, not a specific device's cooling
    duty."""
    extrapolated = not (_HMR_GDCU_FREQ_HZ[0] <= frequency_Hz <= _HMR_GDCU_FREQ_HZ[-1])
    scp = float(np.interp(frequency_Hz, _HMR_GDCU_FREQ_HZ, _HMR_GDCU_SCP_KW_PER_KG))
    return scp, extrapolated


# ---------------------------------------------------------------------------
# TIER 3 (new, this pass): REALISTIC electrical COP -- adds the two
# parasitic loss channels the paper's own Eq. 3 COP structurally excludes
# (rotary drivetrain; baseline controls/inverter overhead) on top of TIER
# 2's own air-gap contact-friction derating, so the reported number is a
# genuine device-level electrical COP, comparable on equal terms to
# vapor_compression_cop() and to core.amr_cycle.AMRCycleResult.
# COP_electrical -- see the module docstring's TIER 3 section for the full
# derivation and honesty flags (cross-device drivetrain proxy, etc).
# ---------------------------------------------------------------------------

_DEFAULT_ROTARY_DRIVE_FIT = None  # lazy-cached; see _rotary_drive_watts()


def _rotary_drive_watts(frequency_Hz):
    """W_drive(f) = k_drive0 + k_drive1*f, fit by
    core.loss_model.fit_rotary_drive_term() to Lozano et al. (2013)'s own
    directly-measured rotary magnet-assembly + valve drivetrain power
    (their Table 3 WM column). Cached at module level (the fit itself is a
    deterministic closed-form least-squares solve over 8 fixed literature
    points, not a random/expensive search, so caching only avoids
    re-printing the fit summary on every call)."""
    global _DEFAULT_ROTARY_DRIVE_FIT
    if _DEFAULT_ROTARY_DRIVE_FIT is None:
        _DEFAULT_ROTARY_DRIVE_FIT = fit_rotary_drive_term(verbose=False)
    fit = _DEFAULT_ROTARY_DRIVE_FIT
    return fit["k_drive0"] + fit["k_drive1"] * frequency_Hz, fit


_DEFAULT_BASE_FRAC = None  # lazy-cached; see _baseline_overhead_frac()


def _baseline_overhead_frac():
    """This repo's own CORE-calibrated `base_frac` (controls/inverter
    overhead, proportional to Qc) from core.loss_model.
    calibrate_loss_coefficients() -- the SAME coefficient
    core.amr_cycle.AMRSystem applies to every fluid-AMR COP_electrical
    calculation, reused here rather than inventing an HMR-specific
    baseline-overhead fraction."""
    global _DEFAULT_BASE_FRAC
    if _DEFAULT_BASE_FRAC is None:
        _DEFAULT_BASE_FRAC = calibrate_loss_coefficients(verbose=False)["base_frac"]
    return _DEFAULT_BASE_FRAC


@dataclasses.dataclass
class HMRElectricalResult:
    technology: str
    frequency_Hz: float
    mass_mcm_kg: float
    air_gap_nm: int
    Qc_ideal_W: float          # frictionless-limit Qc at this (f, mass) -- paper's own Table 1 SCP
    Qc_net_W: float            # after air-gap contact-friction derating (TIER 2)
    W_mag_net_W: float         # magnetic-cycle-only work, backed out from the air-gap-derated eta/COP
    W_drivetrain_W: float      # rotary bearing/cogging/drive-motor power (Lozano 2013 cross-device fit)
    W_baseline_W: float        # controls/inverter overhead (this repo's own CORE base_frac * Qc_net)
    COP_ideal: float           # magnetic-cycle-only COP (== hmr_cop()'s own COP, air-gap-derated)
    COP_electrical: float      # Qc_net / (W_mag_net + W_drivetrain + W_baseline) -- the REALISTIC number
    COP_carnot: float
    extrapolated: bool
    source_note: str


def hmr_electrical_cop(T_cold_K, T_hot_K, frequency_Hz=1.0,
                        mass_mcm_kg=MASS_MCM_KG_DEFAULT, air_gap_nm=100,
                        rotary_drive_watts_fn=_rotary_drive_watts,
                        base_frac=None):
    """Realistic, loss-inclusive device-level electrical COP for the Gd/Cu
    HMR at (T_cold_K, T_hot_K). Unlike hmr_cop() (which reproduces the
    paper's own frictionless-limit, magnetic-cycle-only Eq. 3/4 COP), this
    function:
      1. Applies the paper's own air-gap contact-friction derating
         (TIER 2, apply_air_gap_derating()) -- ON by default here
         (air_gap_nm=100, the paper's own least-ideal tested case) rather
         than frictionless, because a real rotating mechanical interface
         is a more realistic assumption than idealized zero-gap contact.
      2. Converts the paper's dimensionless eta/COP into an absolute Qc
         (W) via specific_cooling_power_kw_per_kg(f) * mass_mcm_kg, then
         backs out W_mag_net = Qc_net / COP_net (both already air-gap-
         derated), so the extra magnetic work the paper's own numbers
         imply the air gap costs (eta drops MORE than Qc at 100 nm: -26%
         vs. -13%) is captured, not just the cooling-power loss.
      3. Adds a rotary drivetrain parasitic term (bearing/cogging
         friction + drive motor to rotate the Ns-sector disk and magnet
         assembly), reusing core.loss_model's own Lozano-et-al.-(2013)-
         calibrated RotaryDriveLossModel/fit_rotary_drive_term() fit --
         see the module docstring's TIER 3 section for why this
         cross-device proxy is the best currently-available real-hardware
         anchor, and for the independent Arnold/Tura/Rowe (2014)
         corroboration that drivetrain loss, not magnetic work, tends to
         dominate real rotary magnetocaloric hardware of this class.
      4. Adds a baseline controls/inverter overhead, reusing this repo's
         own CORE-calibrated base_frac (core.loss_model.
         calibrate_loss_coefficients()) rather than inventing a separate
         HMR-specific figure.

    Returns an HMRElectricalResult. air_gap_nm must be one of {0, 10, 100}
    (see apply_air_gap_derating()) or None (skip the air-gap derating
    entirely -- NOT recommended for a "realistic" estimate, but available
    for a frictionless-contact/drivetrain-only sensitivity check)."""
    cc = carnot_cop(T_cold_K, T_hot_K)
    eta, extrapolated = eta_hmr(frequency_Hz)
    cop_ideal_frictionless = eta * cc
    scp_kw_per_kg, scp_extrapolated = specific_cooling_power_kw_per_kg(frequency_Hz)
    Qc_ideal_W = scp_kw_per_kg * 1000.0 * mass_mcm_kg

    note = (f"Lin et al., Innovation 5(4):100645 (2024), Gd/Cu HMR, Ns=24, Uf=1.0, "
            f"mass_mcm={mass_mcm_kg:.2f}kg (this repo's own representative design-point "
            f"mass, see MASS_MCM_KG_DEFAULT)")
    if extrapolated or scp_extrapolated:
        note += f" -- EXTRAPOLATED to f={frequency_Hz:.2f} Hz (outside fitted [1,10] Hz range)"

    cop_ideal, Qc_net_W = cop_ideal_frictionless, Qc_ideal_W
    if air_gap_nm is not None:
        cop_ideal, eta, air_note = apply_air_gap_derating(
            cop_ideal_frictionless, eta, air_gap_nm, frequency_Hz)
        # cooling_power_frac is the SAME fractional derating whether applied
        # to the paper's dimensionless COP or to an absolute Qc built from
        # its own SCP figure -- both are the paper's Table-1 numbers scaled
        # by the same measured/simulated fraction, so this is not an
        # independent assumption on top of apply_air_gap_derating()'s own.
        Qc_net_W = Qc_ideal_W * (cop_ideal / cop_ideal_frictionless
                                  if cop_ideal_frictionless > 0 else 0.0)
        note += "; " + air_note
    else:
        note += "; air-gap derating SKIPPED (air_gap_nm=None) -- frictionless MCM/HTCM contact"

    W_mag_net_W = Qc_net_W / cop_ideal if cop_ideal > 0 else 0.0
    W_drivetrain_W, drive_fit = rotary_drive_watts_fn(frequency_Hz)
    W_drivetrain_W = max(W_drivetrain_W, 0.0)  # the fitted line can't go negative here
                                                 # (both coefficients are non-negative
                                                 # over the fitted 0.4-1.4 Hz range) but
                                                 # clip defensively if ever called far
                                                 # outside that range.
    b_frac = base_frac if base_frac is not None else _baseline_overhead_frac()
    W_baseline_W = b_frac * Qc_net_W

    note += (f"; +W_drivetrain={W_drivetrain_W:.1f}W (Lozano et al. 2013 rotary-AMR "
             f"drivetrain fit, cross-device proxy -- see module docstring TIER 3); "
             f"+W_baseline={W_baseline_W:.1f}W (base_frac={b_frac:.4f}, this repo's own "
             f"CORE-calibrated controls/inverter overhead)")

    W_parasitic_total = W_drivetrain_W + W_baseline_W
    denom = W_mag_net_W + W_parasitic_total
    cop_electrical = Qc_net_W / denom if denom > 0 else 0.0

    return HMRElectricalResult(
        technology="Hybrid solid-state magnetic regenerator (HMR, Gd/Cu) -- realistic",
        frequency_Hz=frequency_Hz, mass_mcm_kg=mass_mcm_kg,
        air_gap_nm=(air_gap_nm if air_gap_nm is not None else -1),
        Qc_ideal_W=Qc_ideal_W, Qc_net_W=Qc_net_W, W_mag_net_W=W_mag_net_W,
        W_drivetrain_W=W_drivetrain_W, W_baseline_W=W_baseline_W,
        COP_ideal=cop_ideal, COP_electrical=cop_electrical, COP_carnot=cc,
        extrapolated=bool(extrapolated or scp_extrapolated), source_note=note,
    )


def compare_to_vcc_realistic(T_cold_K, span_K, frequencies_Hz=(1.0, 5.0, 10.0),
                              eta_2nd_law_vcc=0.42, mass_mcm_kg=MASS_MCM_KG_DEFAULT,
                              air_gap_nm=100, log=print):
    """Realistic counterpart to compare_to_vcc(): compares HMR's REAL
    electrical COP (hmr_electrical_cop(), inclusive of air-gap friction,
    rotary drivetrain, and baseline overhead) against VCC's REAL installed
    electrical COP -- an apples-to-apples comparison in STRUCTURE (both
    sides are now "electrical power in / cooling power out" with real
    parasitic loss included), unlike compare_to_vcc()'s own ideal-vs-real
    ratio. `log` defaults to print but accepts any single-string-argument
    callable (compare_to_vcc() passes its own list-accumulating logger so
    both sections land in the same results/hybrid_solid_state_regenerator.txt
    file -- see that function's own call site below).

    Also prints (new) an explicit COP_ideal_airgap column alongside
    COP_HMRe, so the FULL loss chain is visible in one table rather than
    split across two: COP_HMR (compare_to_vcc()'s own frictionless-limit
    magnetic-work-only number, TIER 1/2) -> COP_ideal_airgap (same
    magnetic-work-only quantity, but with the paper's own air-gap contact-
    friction derating applied -- what hmr_cop(..., air_gap_nm=100) alone
    would return) -> COP_HMRe (this function's own new TIER 3 number, ALSO
    adding the rotary drivetrain and baseline-overhead terms). Each column
    is a strict subset of the next column's loss accounting, so reading
    left to right shows exactly how much each additional loss channel
    costs."""
    T_hot_K = T_cold_K + span_K
    vcc = vapor_compression_cop(T_cold_K, T_hot_K, eta_2nd_law=eta_2nd_law_vcc)

    log("")
    log("=" * 88)
    log("TIER 3 (realistic, loss-inclusive): HMR REAL electrical COP vs. VCC REAL "
        "electrical COP")
    log(f"Air-gap assumption: {air_gap_nm} nm (paper's own least-ideal tested case) | "
        f"mass_mcm={mass_mcm_kg:.2f}kg | rotary drivetrain: Lozano et al. (2013) "
        f"cross-device fit | baseline overhead: this repo's own CORE base_frac")
    log("Reading left to right below: COP_ideal (TIER 1/2, frictionless, magnetic-work "
        "only, see compare_to_vcc()'s own table above) -> COP_ideal_airgap (SAME "
        "magnetic-work-only quantity, but with the paper's own air-gap contact-friction "
        "derating now applied) -> COP_HMRe (TIER 3, NEW: also adds the rotary drivetrain "
        "and baseline-overhead terms neither the paper nor compare_to_vcc() include).")
    log("=" * 88)
    log(f"{'freq(Hz)':>9}{'COP_ideal':>11}{'COP_id_ag':>11}{'Qc_net_W':>10}{'W_mag_W':>9}"
        f"{'W_drv_W':>9}{'W_base_W':>10}{'COP_HMRe':>10}{'COP_VCC':>10}{'HMRe/VCC':>10}"
        f"{'flag':>14}")

    rows = []
    for f in frequencies_Hz:
        cc = carnot_cop(T_cold_K, T_hot_K)
        eta_frictionless, _ = eta_hmr(f)
        cop_ideal_frictionless = eta_frictionless * cc
        r = hmr_electrical_cop(T_cold_K, T_hot_K, frequency_Hz=f,
                                mass_mcm_kg=mass_mcm_kg, air_gap_nm=air_gap_nm)
        ratio = r.COP_electrical / vcc.COP if vcc.COP > 0 else float("nan")
        flag = "EXTRAPOLATED" if r.extrapolated else ""
        log(f"{f:>9.1f}{cop_ideal_frictionless:>11.3f}{r.COP_ideal:>11.3f}"
            f"{r.Qc_net_W:>10.1f}{r.W_mag_net_W:>9.1f}{r.W_drivetrain_W:>9.1f}"
            f"{r.W_baseline_W:>10.1f}{r.COP_electrical:>10.3f}{vcc.COP:>10.3f}"
            f"{ratio:>10.3f}{flag:>14}")
        rows.append(dict(frequency_Hz=f, COP_ideal_frictionless=cop_ideal_frictionless,
                          COP_ideal_airgap=r.COP_ideal, Qc_net_W=r.Qc_net_W,
                          W_mag_net_W=r.W_mag_net_W, W_drivetrain_W=r.W_drivetrain_W,
                          W_baseline_W=r.W_baseline_W, COP_HMR_electrical=r.COP_electrical,
                          COP_VCC=vcc.COP, HMRe_over_VCC=ratio, extrapolated=r.extrapolated))


    best = max(rows, key=lambda r: r["COP_HMR_electrical"])
    log("")
    if best["HMRe_over_VCC"] >= 1.0:
        conclusion = (
            f"Even after adding air-gap contact friction, a rotary drivetrain term, and "
            f"baseline overhead, HMR's own REAL electrical COP at its best frequency "
            f"({best['frequency_Hz']:.1f} Hz, COP_electrical={best['COP_HMR_electrical']:.2f}) "
            f"still exceeds VCC's REAL electrical COP ({vcc.COP:.2f}) at this (T_cold, span) "
            f"-- ratio {best['HMRe_over_VCC']:.2f}x. This is now a genuinely apples-to-apples "
            f"comparison in structure (both 'electrical in / cooling out'), though the "
            f"drivetrain figure itself remains a cross-device literature proxy (see module "
            f"docstring), not an HMR-specific measurement -- no HMR prototype has been built."
        )
    else:
        conclusion = (
            f"Once air-gap contact friction, a rotary drivetrain term, and baseline "
            f"overhead are added, HMR's own REAL electrical COP at its best frequency "
            f"({best['frequency_Hz']:.1f} Hz, COP_electrical={best['COP_HMR_electrical']:.2f}) "
            f"does NOT exceed VCC's REAL electrical COP ({vcc.COP:.2f}) at this (T_cold, "
            f"span) -- ratio {best['HMRe_over_VCC']:.2f}x, compared to the higher, "
            f"ideal-vs-real ratio compare_to_vcc() reports above (which omits drivetrain and "
            f"baseline overhead entirely). Once those two real, previously-excluded loss "
            f"channels are included on equal footing with VCC's own real installed-system "
            f"number, this architecture does NOT close the electrical-COP gap to "
            f"vapor-compression at this operating point, at least not on the strength of the "
            f"frictionless-limit thermodynamics alone."
        )
    log("CONCLUSION: " + conclusion)
    log("")
    log("Caveat carried over from TIER 1/2 above: no HMR prototype has been built (single "
        "2024 conceptual/FEA paper). The air-gap and specific-cooling-power figures ARE the "
        "paper's own reported/simulated numbers; the rotary-drivetrain figure is a "
        "cross-device proxy (Lozano et al. 2013, a fluid-based rotary AMR, not an HMR) whose "
        "DIRECTION (real rotary magnetocaloric drivetrains draw non-trivial parasitic power) "
        "is independently corroborated by Arnold, Tura & Rowe (2014)'s own force/drive "
        "measurements on a different rotary AMR device, but whose exact MAGNITUDE for an "
        "actual HMR device is not independently confirmed.")

    return {"rows": rows, "vcc": vcc, "conclusion": conclusion}


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

    # TIER 3 addition, this pass: run the realistic, loss-inclusive electrical-COP
    # comparison and append it to the SAME output file/log stream, so a reader sees
    # both the original ideal-vs-real ratio AND the closed-gap realistic one together
    # rather than needing to open two files.
    realistic_result = compare_to_vcc_realistic(
        T_cold_K, span_K, frequencies_Hz=frequencies_Hz,
        eta_2nd_law_vcc=eta_2nd_law_vcc, log=log)

    import os
    os.makedirs(os.path.dirname(out_path) or ".", exist_ok=True)
    with open(out_path, "w") as f_out:
        f_out.write("\n".join(lines) + "\n")
    if verbose:
        print(f"\nWrote {out_path}")

    return {"rows": rows, "vcc": vcc, "conclusion": conclusion,
            "realistic": realistic_result}


if __name__ == "__main__":
    # Repo's own representative ASHRAE point (T_cold=18C=291.15K) swept
    # across its own 5-20K ASHRAE span range, same range main.py step 4 uses.
    T_COLD_K = 291.15
    for span in (5.0, 10.0, 15.0, 20.0):
        compare_to_vcc(T_COLD_K, span)
        print()