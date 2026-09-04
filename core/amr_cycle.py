"""
amr_cycle.py
============
0-D (lumped) Active Magnetic Regenerator (AMR) cycle model.

The AMR cycle (Barclay, 1982; standard reference cycle for room-temperature
magnetic refrigeration) has four processes analogous to a regenerative
Brayton cycle, but with a solid magnetocaloric regenerator bed instead of a
gas:

    1. Adiabatic magnetization   : bed heats by DeltaT_ad (H: 0 -> H_max)
    2. Cold-to-hot fluid flow    : fluid rejects heat to hot reservoir (Qh)
    3. Adiabatic demagnetization : bed cools by DeltaT_ad (H: H_max -> 0)
    4. Hot-to-cold fluid flow    : fluid absorbs heat from cold reservoir (Qc)

This module wraps an MagnetocaloricMaterial (see mce_material.py) with a
regenerator effectiveness / utilization treatment (Engelbrecht & Bahl, 2010;
Kitanovski et al. 2015, Ch. 4) to give cooling capacity, COP and required
magnetic work per cycle as functions of:
    - temperature span (Th - Tc)
    - peak field mu0*H
    - operating frequency f
    - fluid utilization factor U (heat capacity ratio, fluid/regenerator)
    - regenerator effectiveness eps (NTU-based, from thermal.py)
    - blow fraction (fraction of the cycle period spent in cold-to-hot
      flow vs. hot-to-cold flow) -- see BLOW_FRACTION_MASCHE and
      AMRSystem's blow_fraction parameter below (Paper-Mining Pass
      recommendation #1). Default 0.5 (symmetric blow) reproduces this
      model's original behavior exactly.

This is a first-order performance model (matches the level of the AMR
"characteristic curve" approach used in Tusek et al., Int. J. Refrig. 33
(2010) and Nielsen et al., Int. J. Refrig. 34 (2011) 603-616) — good enough
for system-level COP comparison against vapor-compression / liquid cooling,
NOT a replacement for a full 2-D/3-D COMSOL regenerator-bed solve (see
COMSOL_setup.md in the roadmap for that follow-on step).

 (ROADMAP.md) added an AMR cycle-topology switch (`cycle_type`,
default "brayton" = previous behavior unchanged): "ericsson" and
"carnot" variants apply small, documented multipliers to specific cooling
power and the second-law efficiency ceiling, intended to reproduce the
QUALITATIVE ranking Carnot-like >= Ericsson-like >= Brayton-like described
in the AMR-cycle-comparison literature -- see CYCLE_TYPE_FACTORS below for
the honesty flag on why these are illustrative multipliers rather than a
digitization of Kitanovski et al.'s own closed-form relations.

 (ROADMAP.md) added an optional `thermal_diode` parameter
(default None = previous behavior unchanged): when a
`core.thermal_diode.MechanicalContactDiode` instance is supplied,
its (illustrative -- see that module's honesty flag)
`actuation_energy_J_per_cycle * frequency` switching power is added to
W_parasitic, the same additive-parasitic-load accounting used
for hysteresis loss. See `_diode_switching_power_W()` and
`core/thermal_diode_analysis.py`'s docstring for what this module
deliberately does NOT claim (no frequency-ceiling relaxation is
modeled -- see that finding's own writeup).
"""

import numpy as np
from dataclasses import dataclass
from core.mce_material import MagnetocaloricMaterial


# --- AMR cycle topology (, ROADMAP.md) ---
#
# Source intent: Kitanovski et al. (2015) Sect. 4.1.1-4.1.4 ("Characteristics
# of an Ericsson-like AMR Cycle" / "...Hybrid Brayton-Ericsson-like AMR
# Cycle" / "...Carnot-like AMR Cycle" / "Maximum Specific Cooling Power in
# the AMR Cycle", book pp. 104-109) derive closed-form relations for how
# each cycle topology's magnetization/demagnetization-vs-flow phasing
# changes the achievable specific cooling power and second-law efficiency.
#
# HONESTY FLAG (read before trusting the numbers below): this project's own
# copy of that book (see Kitanovski_et_al...pdf in this repo's source
# corpus) is a 30-page front-matter/Chapter-1/table-of-contents excerpt --
# it does NOT include pp. 104-109, so Sect. 4.1.1-4.1.4's actual closed-form
# equations were never available to digitize into this module, unlike
# (e.g.) Ch. 1's thermodynamic relations, which this project's copy does
# contain. The three multipliers below are therefore NOT a reproduction of
# Kitanovski's own formulas. They encode only the QUALITATIVE, well-
# established ranking of the three cycle types (Brayton-like: adiabatic
# magnetization/demagnetization with the fluid static, isofield blows only,
# already this model's previous behavior; Ericsson-like: field change
# happens under continuous fluid contact, so the magnetization/
# demagnetization legs also exchange heat rather than being "wasted"
# adiabatic excursions, which the general AMR-cycle-comparison literature
# -- and this project's own planning note -- describes as
# improving both specific cooling power and second-law efficiency relative
# to Brayton; Carnot-like: the idealized reversible reference bound, a
# theoretical upper limit rather than a claim that a real regenerator bed
# can achieve it) as small, monotonic, illustrative multipliers on this
# model's existing Qc and eta_2nd_law formulas -- NOT a benchmark-fitted or
# book-digitized result. Treat any downstream conclusion that depends on
# the exact SIZE of the Ericsson/Carnot uplift (as opposed to their
# ORDERING relative to Brayton) as provisional until pp. 104-109 are
# actually available to check against.
CYCLE_TYPES = ("brayton", "ericsson", "carnot")

CYCLE_TYPE_FACTORS = {
    # qc_multiplier   : relative change in specific cooling power (applied
    #                   in cooling_capacity()) vs. this model's Brayton-like
    #                   baseline, at fixed span/eps/field/frequency.
    # eta_uplift       : relative change in the eta_2nd_law ceiling used by
    #                   magnetic_work() (0.35 + 0.20*eps at brayton) --
    #                   still passes through the existing np.clip(...,
    #                   0.02, 0.95), so "carnot" approaches but does not
    #                   exceed that ceiling.
    "brayton":  {"qc_multiplier": 1.00, "eta_uplift": 1.00},
    "ericsson": {"qc_multiplier": 1.12, "eta_uplift": 1.15},
    "carnot":   {"qc_multiplier": 1.30, "eta_uplift": 1.35},
}


# --- Flow-waveform asymmetry (blow fraction), Paper-Mining Pass
#     recommendation #1 ---
#
# Source: Masche, Liang, Engelbrecht & Bahl, "Improving magnetic cooling
# efficiency and pulldown by varying flow profiles," Applied Thermal
# Engineering 215 (2022) 118945. DTU rotary AMR device, 13 trapezoidal beds,
# 295g Gd spheres/bed, solenoid-valve-controlled blow fraction (the fraction
# of the cycle period during which fluid flows cold-to-hot vs hot-to-cold).
#
# This model previously had no notion of blow fraction at all -- Qc and the
# magnetic-cycle second-law efficiency implicitly assumed a symmetric 50/50
# split between the two flow directions. The two data points below are the
# ONLY reported operating condition (T_span=16K, U=0.32, f=1.4Hz) with
# reported values at two different blow fractions, so BLOW_FRACTION_MASCHE
# is a two-point calibration, not a full characteristic curve:
BLOW_FRACTION_MASCHE = {
    "source": "Masche, Liang, Engelbrecht & Bahl, Appl. Thermal Eng. 215 (2022) "
              "118945 -- DTU rotary AMR, 13 trapezoidal beds, 295g Gd spheres/bed",
    "operating_point": {"T_span_K": 16.0, "U": 0.32, "frequency_Hz": 1.4},
    "low": {"blow_fraction": 0.250, "Qc_W": 70.0, "exergy_eff": 0.026},
    "best_found": {"blow_fraction": 0.416, "Qc_W": 330.0, "exergy_eff": 0.174},
    # Best blow fraction found across BOTH the 6K and 16K spans the paper
    # tested (i.e. not span-specific, per the source paper's own framing).
    # Lower blow fractions favor faster pulldown (~30% faster to reach ~14K
    # span at a lower blow fraction) -- a transient-dynamics effect this
    # steady-state 0-D cycle model does not represent and is NOT calibrated
    # against here.
}


def _blow_fraction_multiplier(blow_fraction, value_at_low, value_at_peak,
                                bf_low=BLOW_FRACTION_MASCHE["low"]["blow_fraction"],
                                bf_peak=BLOW_FRACTION_MASCHE["best_found"]["blow_fraction"],
                                bf_symmetric=0.5):
    """Relative correction factor for a blow-fraction-sensitive quantity
    (Qc or the second-law efficiency), calibrated to a parabola through the
    two Masche et al. (2022) data points at bf_low and bf_peak, and
    NORMALIZED so that bf_symmetric=0.5 (this model's pre-existing implicit
    symmetric-blow assumption) returns a multiplier of 1.0 -- i.e. this is
    a RELATIVE correction layered on top of the model's existing behavior,
    not an absolute reproduction of the DTU device.

    Honesty flags:
      1. Only TWO (blow_fraction, value) points are available, at ONE
         operating condition (T_span=16K, U=0.32, f=1.4Hz). A parabola
         needs three points to be fit outright, so bf_peak is additionally
         TREATED as the location of the true continuous-curve maximum --
         the source paper only reports it as "best blow fraction found
         across both spans tested," not as a confirmed peak. This is an
         explicit modeling choice, not a literature-derived shape.
      2. This calibration is applied identically regardless of T_span, U or
         frequency -- the source paper only reports the two-point comparison
         at ONE (T_span, U, f) combination, so extrapolating the SHAPE of
         this curve to other operating points is unvalidated.
      3. Values are clipped to stay within [0.05, 3.0]x, since the raw
         parabola diverges (goes negative, then blows up) far from bf_peak
         -- physically nonsensical outside roughly the tested 0.1-0.6 window.
    """
    k = (value_at_peak - value_at_low) / (bf_low - bf_peak) ** 2

    def value(bf):
        return max(value_at_peak - k * (bf - bf_peak) ** 2, 0.0)

    baseline = value(bf_symmetric)
    if baseline <= 0:
        return 1.0
    return float(np.clip(value(blow_fraction) / baseline, 0.05, 3.0))


@dataclass
class AMRCycleResult:
    T_span: float
    Qc: float          # W, cooling capacity
    Qh: float           # W, heat rejected
    W_mag: float         # W, net magnetic (thermodynamic-cycle) work input
    W_parasitic: float    # W, pump + motor-drive overhead (see note below) --
                            # also includes material thermal-
                            # hysteresis loss (_hysteresis_power_W()), 0.0
                            # for GADOLINIUM/previous materials.
                            # also includes thermal-diode
                            # actuation switching power
                            # (_diode_switching_power_W()), 0.0 unless a
                            # thermal_diode is supplied
    COP: float            # ideal magnetic-cycle-only COP (Qc / W_mag)
    COP_electrical: float  # device-level electrical COP (Qc / (W_mag + W_parasitic))
                             # -- this is the number comparable to published
                             # "COPe" / device COP figures, and to the
                             # vapor-compression/liquid-cooling baselines in
                             # baseline_cooling.py, which are also electrical.
    exergy_eff: float    # second-law efficiency vs Carnot (magnetic-cycle-only)


class AMRSystem:
    PUMP_MOTOR_EFFICIENCY_LITERATURE = 0.6
    #  addition. Two independent industry sources agree small
    # (fractional-scale, lab/prototype-appropriate) centrifugal pumps run
    # 50-70% wire-to-water efficiency (Pumps & Systems, "How to Define &
    # Measure Centrifugal Pump Efficiency": "typical efficiencies are 55
    # percent for small pumps"; Linquip, "Calculation of Pump Efficiency":
    # "smaller pumps typically fall into the range of 50 to 70 percent")
    # -- 0.6 is their midpoint, and independently corroborates (rather than
    # being derived from) the SAME 0.5-0.7 range a user-supplied document
    # separately proposed for this exact parameter. See
    # `_geometry_pumping_power_W()`'s own docstring for exactly where this
    # is applied and, just as importantly, where it is deliberately NOT
    # applied.
    #
    # NOT the class default: `pump_motor_efficiency` defaults to 1.0 (see
    # __init__ below) -- i.e. idealized, no pump/motor loss, EXACTLY the
    # previous behavior -- matching this repo's own established
    # discipline of adding a new, real, literature-
    # grounded capability as an OPT-IN parameter rather than silently
    # changing every existing caller's numeric output. Every existing
    # optimize.py NSGA-III run always sets particle_diameter (it is a
    # design variable, bounds [0.05,2.0]mm), so changing this default
    # would silently change every production Pareto front this repo has
    # ever generated. Callers that want the literature-grounded 0.6 value
    # (e.g. `optimize.py`, if/when it opts in) must pass
    # `pump_motor_efficiency=AMRSystem.PUMP_MOTOR_EFFICIENCY_LITERATURE`
    # explicitly -- see ROADMAP.md for the concrete
    # follow-up this leaves open.

    def __init__(self, material: MagnetocaloricMaterial, mu0H_max: float,
                 mass_regenerator: float, frequency: float,
                 fluid_cp: float = 4186.0, fluid_mdot: float = 0.05,
                 regenerator_effectiveness: float = 0.85,
                 parasitic_fraction: float = 0.15,
                 loss_model=None, use_ntu_thermal_model: bool = False,
                 blow_fraction: float = 0.5,
                 particle_diameter: float = None,
                 bed_cross_section_area: float = 0.002,
                 hypereg_n_parallel: int = None,
                 cycle_type: str = "brayton",
                 thermal_diode=None,
                 pump_motor_efficiency: float = 1.0,
                 no_load_span_override: float = None):
        """
        material               : MagnetocaloricMaterial instance
        mu0H_max                : peak applied field, Tesla
        mass_regenerator        : kg of magnetocaloric material in the bed
        frequency                : AMR cycle frequency, Hz
        fluid_cp                 : heat transfer fluid specific heat, J/(kg K)
        fluid_mdot                : fluid mass flow rate, kg/s
        regenerator_effectiveness : NTU-based regenerator effectiveness (0-1),
                                     from thermal.py NTU correlation
        parasitic_fraction        : pump + magnet-motor-drive electrical
                                     overhead, as a fraction of Qc, ADDED ON
                                     TOP of the ideal magnetic-cycle work to
                                     get device-level electrical COP. The default 0.15 is a literature-calibrated
                                     value based on two comparably sized lab
                                     devices (DTU rotary Gd: 0.171, Tusek
                                     single-bed Gd: 0.118 - see
                                     core/validation_system.py). The large
                                     Astronautics naval-cooler prototype
                                     implied 0.453, which Jacobs et al. (2014)
                                     attribute explicitly to "electrical
                                     components with mediocre efficiency" at
                                     that scale/vintage. Treat 0.15 as an
                                     optimistic lab-scale figure rather than a
                                     production-hardware guarantee, and widen
                                     it in any economics sensitivity study.
                                     IGNORED if loss_model is provided.
        loss_model                : optional core.loss_model.StateDependentLossModel.
                                    If provided, W_parasitic is computed as a
                                    function of (frequency, mu0H_max,
                                    fluid_mdot, Qc) instead of the constant
                                    parasitic_fraction*Qc. This restores
                                    field-, frequency-, and flow-dependent
                                    electrical losses to COP_electrical.
        blow_fraction              : fraction of the AMR cycle period during
                                     which fluid flows cold-to-hot (vs.
                                     hot-to-cold) -- a real degree of freedom
                                     this model previously had no notion of
                                     (Masche et al. 2022; see
                                     BLOW_FRACTION_MASCHE above). Default
                                     0.5 (symmetric blow) exactly reproduces
                                     this model's pre-existing behavior --
                                     passing a different value applies a
                                     RELATIVE correction to Qc and the
                                     second-law efficiency via
                                     _blow_fraction_multiplier(), calibrated
                                     to the one reported (Qc, exergy_eff)
                                     pair at each of two tested blow
                                     fractions. See
                                     _blow_fraction_multiplier's docstring
                                     for the honesty flags on this
                                     calibration (two points, one operating
                                     condition, extrapolated shape).
        particle_diameter          :  addition. Packed-sphere-bed
                                     particle diameter, m. Default None
                                     preserves ALL previous behavior
                                     exactly (regenerator_effectiveness is
                                     used as before if
                                     use_ntu_thermal_model=False, or the
                                     NTU model's own default particle
                                     diameter is used if True). When set
                                     AND use_ntu_thermal_model=True, this
                                     (a) feeds core.thermal.
                                     regenerator_effectiveness()'s
                                     particle_diameter, coupling
                                     regenerator geometry to thermal
                                     effectiveness the way
                                     geometry_analysis.py already
                                     demonstrates, and (b) computes a
                                     geometry-explicit hydraulic pumping
                                     power (core.thermal.
                                     pumping_power_packed_bed(), Tusek et
                                     al. 2013) that REPLACES loss_model's
                                     generic k_pump*mdot**2 term (via
                                     StateDependentLossModel.
                                     parasitic_power's
                                     pumping_power_override) rather than
                                     adding to it -- see that function's
                                     docstring for why: k_pump is
                                     CORE-calibrated against real devices'
                                     TOTAL parasitic power at their own
                                     (unknown-to-this-model) geometries, so
                                     adding a geometry-explicit hydraulic
                                     term on top would double-count the
                                     pumping-loss channel. If particle_
                                     diameter is set but no loss_model is
                                     provided, the geometry-explicit
                                     pumping term is NOT applied (the
                                     constant parasitic_fraction*Qc model
                                     has no separate pumping component to
                                     substitute) -- a documented, not
                                     silent, limitation.
        bed_cross_section_area      : m^2, passed through to the NTU/
                                     pumping-power geometry calculations
                                     when particle_diameter is set. Default
                                     0.002 m^2 matches thermal.py's and
                                     geometry_analysis.py's own default
                                     (a representative ~5x4cm bed face).
        hypereg_n_parallel          :  addition. If set (and
                                     particle_diameter and loss_model are
                                     also set), uses core.thermal.
                                     pumping_power_packed_bed_hypereg()
                                     instead of the conventional-series
                                     pumping_power_packed_bed() -- i.e.
                                     models a Hypereg-style parallel-
                                     hydraulic regenerator split into this
                                     many sub-regenerators (Klinar et al.
                                     2024; see results/hypereg_findings.md
                                     and core/hypereg_analysis.py). Default
                                     None reproduces the conventional
                                     (non-Hypereg) geometry-explicit
                                     pumping term.
        cycle_type                  :  addition. One of "brayton"
                                     (default), "ericsson", or "carnot" --
                                     see CYCLE_TYPE_FACTORS above and its
                                     honesty flag. Default "brayton"
                                     preserves ALL previous behavior
                                     exactly (qc_multiplier=eta_uplift=1.0).
                                     Raises ValueError for any other value.
        thermal_diode                :  addition. Optional
                                     core.thermal_diode.MechanicalContactDiode
                                     instance. Default None preserves ALL
                                     previous behavior exactly (no
                                     switching-power term is added to
                                     W_parasitic). When supplied, its
                                     switching_power_W(frequency) --
                                     an illustrative, unbenchmarked
                                     actuation cost, see that module's
                                     honesty flag -- is added to
                                     W_parasitic unconditionally, the same
                                     accounting pattern used for
                                     hysteresis loss. This parameter does
                                     NOT change cooling_capacity() or
                                     magnetic_work() in any way: no
                                     frequency-ceiling relaxation or
                                     rectification-ratio heat-transfer
                                     benefit is modeled here (see
                                     core/thermal_diode_analysis.py's
                                     docstring for why, and for what a
                                     rectification-ratio-driven heat
                                     -transfer benefit would require that
                                     this repo does not yet have the data
                                     to support).
        no_load_span_override        : opt-in, additive fix for the
                                     "no regenerative amplification" gap
                                     documented in
                                     results/regenerative_amplification_diagnostic.txt
                                     and core/regenerator_1d.py. Default
                                     None preserves ALL pre-existing
                                     behavior exactly: cooling_capacity()
                                     clamps Qc to zero at
                                     T_span == 2*dTad_noload (the
                                     material's own single-blow adiabatic
                                     dT at mid-bed temperature), which
                                     cannot represent the temperature-
                                     profile buildup a real packed bed
                                     achieves over many cycles. When set
                                     (a span in K), THIS value replaces
                                     2*dTad_noload as the span at which Qc
                                     reaches zero -- letting the system
                                     reach spans a single-blow model
                                     structurally cannot. Qc's MAGNITUDE
                                     at zero span is unaffected (still set
                                     by dTad_noload -- amplification
                                     extends how far the bed's gradient
                                     can be pushed, not the material's own
                                     per-blow capacity). Intended to be
                                     populated from
                                     core.regenerator_1d.regenerative_span_cap()
                                     -- a genuine multi-cycle transient
                                     simulation, NOT a formula -- computed
                                     ONCE per design point and passed in
                                     here, since that simulation takes
                                     tens of seconds per call and cannot
                                     be run inside cooling_capacity()
                                     itself (called many times per
                                     optimization/sweep run). HONESTY
                                     FLAG, read before using: that
                                     simulation's own validation
                                     (results/regenerator_1d_validation.txt)
                                     shows directionally-inconsistent
                                     error against the three directly-
                                     measured benchmarks (+112%, -92%,
                                     -61%) -- it is NOT yet an
                                     independently-calibrated quantitative
                                     model. This parameter exists so the
                                     capability is available and testable
                                     (see
                                     validation_system.run_regenerative_amplification_override_check(),
                                     wired into main.py as an additive
                                     diagnostic step, NOT used by any
                                     default caller -- optimize.py,
                                     cascade.py, and every existing test
                                     are unaffected), not because the
                                     override is ready to replace the old
                                     cap as this project's default.
        """
        if cycle_type not in CYCLE_TYPE_FACTORS:
            raise ValueError(
                f"cycle_type must be one of {CYCLE_TYPES}, got {cycle_type!r}")
        self.mat = material
        self.mu0H_max = mu0H_max
        self.m_reg = mass_regenerator
        self.f = frequency
        self.cp_f = fluid_cp
        self.mdot_f = fluid_mdot
        self.eps = regenerator_effectiveness
        self.parasitic_fraction = parasitic_fraction
        self.loss_model = loss_model
        self.use_ntu_thermal_model = use_ntu_thermal_model
        self.blow_fraction = blow_fraction
        self.particle_diameter = particle_diameter
        self.bed_cross_section_area = bed_cross_section_area
        self.hypereg_n_parallel = hypereg_n_parallel
        self.cycle_type = cycle_type
        self.thermal_diode = thermal_diode
        self.pump_motor_efficiency = pump_motor_efficiency
        self.no_load_span_override = no_load_span_override
        self._last_ntu_info = None

    def _cycle_type_factor(self):
        """ addition. Returns this system's {"qc_multiplier",
        "eta_uplift"} dict from CYCLE_TYPE_FACTORS -- see that constant's
        docstring for the honesty flags on where these numbers do (and do
        not) come from."""
        return CYCLE_TYPE_FACTORS[self.cycle_type]

    def _hysteresis_power_W(self):
        """ addition. Returns the parasitic electrical power (W)
        attributable to irreversible thermal-hysteresis loss in the
        regenerator material itself, computed as:

            W_hys = material.hysteresis_loss_J_per_kg * mass_regenerator * frequency

        i.e. the per-kg energy dissipated over one FULL field-up/field-down
        hysteresis loop (see core.first_order_mce.FirstOrderMCEMaterial's
        hysteresis_loss_J_per_kg field docstring for where these numbers
        come from and their honesty flags) times the kg of material in the
        bed times the number of loops per second. getattr(..., 0.0) means
        this returns exactly 0.0 for GADOLINIUM (mce_material.py's
        MagnetocaloricMaterial has no such field at all -- a genuinely,
        not just approximately, hysteresis-free second-order transition)
        and for any FirstOrderMCEMaterial instance that predates ,
        so every previous caller and test gets IDENTICAL numbers to
        before this addition.

        Modeling simplification, stated plainly: this treats hysteresis
        loss as scaling with frequency the same way eddy-current loss
        does (once per cycle, material-intrinsic, independent of Qc/span)
        rather than folding it into magnetic_work()'s Carnot-referenced
        ideal-work calculation or eta_2nd_law -- i.e. it is accounted for
        as an ADDITIONAL parasitic electrical load, the same accounting
        choice already used for eddy/pump/base losses in
        core.loss_model.StateDependentLossModel. This is a reasonable
        first-order treatment, not a claim that hysteresis loss is
        electrically identical in character to eddy-current or pumping
        loss -- see ROADMAP.md for the full discussion
        of why this accounting choice was made over folding it into
        eta_2nd_law instead.
        """
        hysteresis_loss_J_per_kg = getattr(self.mat, "hysteresis_loss_J_per_kg", 0.0)
        return hysteresis_loss_J_per_kg * self.m_reg * self.f

    def _diode_switching_power_W(self):
        """ addition. Returns the parasitic electrical power (W)
        to actuate this system's `thermal_diode` (a
        core.thermal_diode.MechanicalContactDiode) once per AMR cycle, or
        exactly 0.0 if thermal_diode is None (the default) -- so every
        previous caller and test gets IDENTICAL numbers to before this
        addition, the same backward-compatibility guarantee
        used for particle_diameter/cycle_type/etc. Unlike
        _hysteresis_power_W(), this does NOT scale with mass_regenerator:
        see core.thermal_diode.MechanicalContactDiode's docstring for why
        actuation energy is modeled as a per-diode, not per-kg, quantity."""
        if self.thermal_diode is None:
            return 0.0
        return self.thermal_diode.switching_power_W(self.f)

    def _blow_fraction_qc_multiplier(self):
        ref = BLOW_FRACTION_MASCHE
        return _blow_fraction_multiplier(
            self.blow_fraction,
            value_at_low=ref["low"]["Qc_W"], value_at_peak=ref["best_found"]["Qc_W"])

    def _blow_fraction_eta_multiplier(self):
        ref = BLOW_FRACTION_MASCHE
        return _blow_fraction_multiplier(
            self.blow_fraction,
            value_at_low=ref["low"]["exergy_eff"], value_at_peak=ref["best_found"]["exergy_eff"])

    def _effective_eps(self):
        """If use_ntu_thermal_model is enabled, compute regenerator
        effectiveness from the NTU model (core/thermal.py) instead of using
        the prescribed constant value. This allows regenerator mass to
        influence cooling capacity. when particle_diameter is
        also set, it is passed through so geometry (not just mass/
        frequency/mdot) affects the NTU calculation too -- see __init__'s
        docstring."""
        if not self.use_ntu_thermal_model:
            return self.eps
        from core.thermal import regenerator_effectiveness as ntu_eps
        kwargs = dict(bed_cross_section_area=self.bed_cross_section_area)
        if self.particle_diameter is not None:
            kwargs["particle_diameter"] = self.particle_diameter
        info = ntu_eps(self.m_reg, self.f, self.mdot_f, **kwargs)
        self._last_ntu_info = info
        return info["eps"]

    def _geometry_pumping_power_W(self):
        """ addition. Returns the geometry-explicit hydraulic
        pumping power (W) if particle_diameter is set, else None (meaning
        "no override -- use loss_model's generic k_pump*mdot**2 term
        unchanged", the previous behavior). Uses the Hypereg parallel-
        hydraulic variant (core.thermal.pumping_power_packed_bed_hypereg)
        instead of the conventional series-flow one if hypereg_n_parallel
        is also set -- see results/hypereg_findings.md and
        core/hypereg_analysis.py.

         addition: the value returned here is divided by
        `self.pump_motor_efficiency` (default 1.0, i.e. no change, unless
        the caller explicitly opts into a lower value -- see
        PUMP_MOTOR_EFFICIENCY_LITERATURE's own comment above for the
        literature-grounded 0.6 figure available for that opt-in) before being returned,
        converting core.thermal.pumping_power_packed_bed()'s own IDEALIZED
        hydraulic power (dP * volumetric flow, no pump/motor losses -- see
        that function's own docstring) into an estimate of the ELECTRICAL
        power actually drawn to achieve it -- the physically correct
        quantity for a `pumping_power_override` that ultimately feeds
        COP_electrical (electrical power in / cooling power out).

        WHY THIS ONLY APPLIES HERE, NOT TO loss_model's generic
        k_pump*mdot**2 term: that term is CORE-calibrated (NNLS fit) to
        real AMR devices' own directly-reported ELECTRICAL parasitic power
        -- i.e. any real device's own pump/motor inefficiency is ALREADY
        baked into the fitted k_pump coefficient by construction, since the
        fit target IS electrical power, not idealized hydraulic power.
        Dividing that term by an efficiency AGAIN would double-count the
        same physical loss. This only applies to the geometry-explicit
        path (this function), which starts from a purely idealized
        hydraulic-power calculation with no efficiency loss represented at
        all -- confirmed by inspection of pumping_power_packed_bed()'s own
        docstring ("no pump/motor efficiency is applied").

        geometry_analysis.py's own direct calls to
        core.thermal.pumping_power_packed_bed()/pumping_power_parallel_plate()
        (its Table-3-comparison diagnostics) are UNCHANGED by this --
        those calls do not go through AMRSystem/this method at all, and
        that module's own docstring already explicitly states its
        "augmented COP" is idealized-hydraulic and NOT meant to be read as
        a production electrical-COP estimate; this phase does not touch
        that intentional idealization."""
        if self.particle_diameter is None:
            return None
        from core.thermal import pumping_power_packed_bed, pumping_power_packed_bed_hypereg
        if self.hypereg_n_parallel is not None:
            info = pumping_power_packed_bed_hypereg(
                self.mdot_f, particle_diameter=self.particle_diameter,
                bed_cross_section_area=self.bed_cross_section_area,
                mass_regenerator=self.m_reg,
                n_parallel_subregenerators=self.hypereg_n_parallel)
        else:
            info = pumping_power_packed_bed(
                self.mdot_f, particle_diameter=self.particle_diameter,
                bed_cross_section_area=self.bed_cross_section_area,
                mass_regenerator=self.m_reg)
        return info["P_pump_W"] / self.pump_motor_efficiency

    def _geometry_eddy_power_W(self):
        """ addition. Returns the geometry-explicit intragranular
        eddy-current power (W) if particle_diameter is set, else 0.0
        (meaning "no additional term -- use loss_model's CORE-calibrated
        k_eddy*f**2*mu0H**2 support-structure term alone, unchanged", the
        previous behavior). Unlike `_geometry_pumping_power_W()`
        (which OVERRIDES the generic k_pump term), this is ADDED on top of
        k_eddy -- see loss_model.StateDependentLossModel.parasitic_power()'s
        `intragranular_eddy_power_W` parameter docstring for why the two
        eddy channels are physically additive rather than alternatives."""
        if self.particle_diameter is None:
            return 0.0
        from core.thermal import intragranular_eddy_power
        return intragranular_eddy_power(
            self.f, self.mu0H_max, particle_diameter=self.particle_diameter,
            mass_regenerator=self.m_reg)

    def cooling_capacity(self, T_cold, T_span):
        """Cooling capacity Qc (W) at a given no-load DeltaT_ad and imposed
        span, using the regenerator-effectiveness degradation model:
            Qc = eps * mdot*cp * (DeltaT_ad_local - T_span/2) ... averaged
        which reduces to the standard 'characteristic curve' shape: Qc is
        maximum at zero span and falls roughly linearly to zero at the
        no-load span (Nielsen et al. 2011).

        MODEL LIMITATION (Track A3, documented rather than "fixed" -- see
        ROADMAP.md): `span_fraction = max(0, 1 - T_span/(2*dTad_noload))`
        is a LINEAR approximation of that characteristic curve. It is
        exact at its two anchor points (span_fraction=1 at zero span,
        =0 at the no-load span) but produces a sharper, straight-line
        cutoff near the no-load span limit than a real AMR device would
        show -- published Qc(span) curves (e.g. Nielsen et al. 2011;
        Tusek et al. 2013) typically round off gradually near their
        no-load limit rather than hitting a hard corner. No literature
        source for the exact fall-off shape near that limit was found in
        this project's corpus, so this module intentionally keeps the
        documented linear clamp rather than inventing an undocumented
        smoothing function -- a soft cutoff with no citation would be a
        downgrade (unfounded precision), not a fix. Treat Qc/COP values
        within roughly the last ~10-15% of a material's no-load span as a
        conservative lower bound rather than a precise prediction.

        SPAN CAP: `self.no_load_span_override`, if set (see __init__), is
        used in place of 2*dTad_noload for the point where span_fraction
        reaches zero -- see that parameter's own docstring for what this
        does and does not fix, and its honesty flag before relying on it.
        Default None (unset) leaves this method's behavior completely
        unchanged from before that parameter existed."""
        T_hot = T_cold + T_span
        T_mid = 0.5 * (T_cold + T_hot)
        H = self.mu0H_max / (4 * np.pi * 1e-7)
        dTad_noload = float(self.mat.delta_T_adiabatic(np.array([T_mid]), H)[0])
        if dTad_noload <= 0:
            return 0.0, dTad_noload
        span_cap = (self.no_load_span_override if self.no_load_span_override is not None
                    else 2 * dTad_noload)
        span_fraction = max(0.0, 1.0 - T_span / span_cap) if span_cap > 0 else 0.0
        eps = self._effective_eps()
        Qc = eps * self.mdot_f * self.cp_f * dTad_noload * span_fraction
        Qc *= self._blow_fraction_qc_multiplier()
        Qc *= self._cycle_type_factor()["qc_multiplier"]
        return max(Qc, 0.0), dTad_noload

    def cooling_capacity_span_sweep(self, T_cold, spans, n_dense_points=400):
        """Evaluates cooling_capacity() over the requested `spans` PLUS an
        internal dense scan grid (0 to max(spans), n_dense_points points),
        then applies a running-minimum monotonicity clamp: Qc at a given
        span can never exceed Qc at any smaller span already evaluated in
        the same sweep.

        WHY THIS EXISTS (see core/validation_system.py's
        diagnose_qc_feasibility_reopening() for the full mechanism, and
        core/mce_material.py's magnetic_heat_capacity() docstring for the
        root cause): the mean-field magnetic heat capacity has a genuine
        finite-jump discontinuity approaching Tc from below. Because
        cooling_capacity() evaluates dTad_noload at a single T_mid per
        call with no knowledge of neighboring spans, a raw single-point Qc
        evaluation can REPORT A LARGER value at a bigger span than at a
        smaller one, right where T_mid crosses that discontinuity --
        physically backwards, since a real device's achievable cooling
        capacity cannot increase as the demanded span widens.

        WHY A DENSE GRID, NOT JUST THE CALLER'S OWN SPANS: if the reopening
        artifact's positive-to-negative-to-positive excursion falls
        entirely BETWEEN two widely-spaced requested spans, a running
        minimum computed only over those sparse points would never see the
        intervening zero and would under-clamp (report the reopened value
        as if it were legitimate). The internal dense grid guarantees the
        clamp sees the discontinuity regardless of how sparse the caller's
        own `spans` are.

        This is a POST-HOC MONOTONICITY CLAMP, not a new physical
        mechanism and not a fit to any literature value -- it can only
        reduce a reported Qc relative to a smaller span already evaluated
        in the sweep, never invent or increase a number.
        `cooling_capacity()` itself is completely unchanged (same
        signature, same return value, same default behavior) -- this is
        an additive, opt-in method; no existing caller anywhere in this
        repo uses it, so no existing result changes because of it.

        Returns a list of dicts, one per requested span (in the order
        given), each with span_K, Qc_raw_W (the unclamped
        cooling_capacity() value), Qc_W (the clamped value actually safe
        to report/plot), and dTad_noload_K (from cooling_capacity() at
        that exact span, unaffected by the clamp)."""
        spans = [float(s) for s in spans]
        max_span = max(spans) if spans else 0.0
        dense_grid = np.linspace(0.0, max_span, n_dense_points) if max_span > 0 else np.array([0.0])
        combined = np.unique(np.concatenate([dense_grid, np.array(spans)]))
        combined.sort()

        raw_Qc = np.empty_like(combined)
        dTads = np.empty_like(combined)
        for i, s in enumerate(combined):
            Qc_i, dTad_i = self.cooling_capacity(T_cold, float(s))
            raw_Qc[i] = Qc_i
            dTads[i] = dTad_i
        running_min_Qc = np.minimum.accumulate(raw_Qc)

        # combined is sorted+unique, so searchsorted gives an exact index
        # for every requested span (each is itself a member of combined).
        idx_map = np.searchsorted(combined, spans)
        out = []
        for s, i in zip(spans, idx_map):
            out.append({"span_K": s, "Qc_raw_W": float(raw_Qc[i]),
                        "Qc_W": float(running_min_Qc[i]),
                        "dTad_noload_K": float(dTads[i])})
        return out

    def magnetic_work(self, T_cold, T_span, Qc):
        """Net magnetic work input per unit time (W). Approximated from the
        entropy generated by finite-effectiveness regeneration plus the
        ideal (Carnot-referenced) work for the delivered Qc, following the
        second-law decomposition in Kitanovski et al. (2015), Ch. 6:
            W = Qc * (Th/Tc - 1) / eta_2nd_law
        where eta_2nd_law captures AMR irreversibilities (regenerator
        mismatch, viscous dissipation, demagnetization losses) and is taken
        as a literature-informed 0.35-0.55 for well-designed lab-scale AMRs
        (Tusek et al. 2010; Eriksen et al. 2015, Int. J. Refrig. 58)."""
        T_hot = T_cold + T_span
        carnot_work = Qc * (T_hot / T_cold - 1.0) if T_cold > 0 else np.inf
        eta_2nd_law = 0.35 + 0.20 * self._effective_eps()  # 0.35 at eps=0 .. 0.52 at eps=0.85
        eta_2nd_law *= self._cycle_type_factor()["eta_uplift"]
        eta_2nd_law *= self._blow_fraction_eta_multiplier()
        eta_2nd_law = float(np.clip(eta_2nd_law, 0.02, 0.95))
        W = carnot_work / max(eta_2nd_law, 1e-3)
        return W, eta_2nd_law

    def run(self, T_cold, T_span) -> AMRCycleResult:
        Qc, dTad = self.cooling_capacity(T_cold, T_span)
        W, eta2 = self.magnetic_work(T_cold, T_span, Qc)
        if self.loss_model is not None:
            pump_override = self._geometry_pumping_power_W()
            eddy_intragranular = self._geometry_eddy_power_W()
            W_parasitic = self.loss_model.parasitic_power(
                self.f, self.mu0H_max, self.mdot_f, Qc,
                pumping_power_override=pump_override,
                intragranular_eddy_power_W=eddy_intragranular)
        else:
            W_parasitic = self.parasitic_fraction * Qc
        # hysteresis loss is added HERE, unconditionally, rather
        # than threaded through loss_model.parasitic_power() -- this
        # deliberately catches BOTH the loss_model and the constant-
        # parasitic_fraction branches above with one code path (e.g.
        # cascade.py's _single_stage() baseline helper does not pass a
        # loss_model at all; it would otherwise silently miss this term).
        # Returns 0.0 for GADOLINIUM and for any previous
        # FirstOrderMCEMaterial -- see _hysteresis_power_W()'s docstring.
        W_parasitic += self._hysteresis_power_W()
        # thermal-diode actuation switching power, added HERE
        # unconditionally for the same reason (catches both the
        # loss_model and constant-parasitic_fraction branches above with
        # one code path). Returns 0.0 unless thermal_diode is supplied --
        # see _diode_switching_power_W()'s docstring. Deliberately does
        # NOT change Qc or W_mag: no rectification-ratio heat-transfer
        # benefit or frequency-ceiling relaxation is modeled here (see
        # core/thermal_diode_analysis.py's docstring for why).
        W_parasitic += self._diode_switching_power_W()
        Qh = Qc + W
        COP = Qc / W if W > 0 else 0.0
        COP_electrical = Qc / (W + W_parasitic) if (W + W_parasitic) > 0 else 0.0
        T_hot = T_cold + T_span
        COP_carnot = T_cold / (T_hot - T_cold) if T_hot > T_cold else np.inf
        exergy_eff = COP / COP_carnot if np.isfinite(COP_carnot) and COP_carnot > 0 else 0.0
        return AMRCycleResult(T_span=T_span, Qc=Qc, Qh=Qh, W_mag=W,
                               W_parasitic=W_parasitic, COP=COP,
                               COP_electrical=COP_electrical, exergy_eff=exergy_eff)

    def characteristic_curve(self, T_cold, spans):
        return [self.run(T_cold, s) for s in spans]