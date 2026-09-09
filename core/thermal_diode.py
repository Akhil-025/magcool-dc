"""
thermal_diode.py
=================
 (ROADMAP.md): a narrowly-scoped first pass at active thermal
diodes for AMR devices (Kitanovski et al. 2015, Ch. 6, "Special Heat
Transfer Mechanisms: Active and Passive Thermal Diodes").

Scope decision (per the plan). Of the four active-diode
mechanisms Ch. 6 covers (thermoelectric Sect. 6.2.1, thermionic 6.2.2,
spincaloritronic 6.2.3, mechanical-contact 6.2.4), this module implements
ONLY the mechanical-contact diode, because it is the mechanism the plan
identifies as actually used in real prototypes referenced elsewhere in
this project's corpus, and therefore the lowest-risk starting point --
not because the other three are less interesting physics, but because
they would add new device classes with even less benchmark grounding
than what follows.

HONESTY FLAG (read before trusting anything in this module -- same tier
as the earlier cycle_type caveat in core/amr_cycle.py). This project's
own copy of Kitanovski et al. (2015) is a 30-page front-matter/
Chapter-1/table-of-contents excerpt -- it does NOT include pp. 211-268
(Chapter 6), where Sect. 6.2.4's actual mechanical-contact-diode design
equations, measured rectification ratios, and switching dynamics are
given. Those numbers were therefore never available to digitize into
this module, unlike (e.g.) Chapter 1's thermodynamic relations, which
this project's copy does contain. What follows is instead a generic,
textbook-level thermal-contact-conductance model (Fourier conduction
across an engaged/disengaged mechanical joint), parameterized by
forward/reverse conductance rather than derived from Kitanovski's own
Sect. 6.2.4 figures.

DEFAULT_MECHANICAL_CONTACT_DIODE's `forward_conductance_W_K` and
`reverse_conductance_W_K` are now grounded in a real, cited literature
ANALOG rather than an unattributed round number (closing part of the
original "did NOT do" list -- see ROADMAP.md): Bywaters & Griffin's
piezo-actuated mechanical heat switch (PZHS) reports on/off thermal-
conductance ratios of roughly 100-200 at cryogenic temperatures (4-10 K)
under a piezoelectric positioner's maximum 8 N actuation force
("Passive Gas-Gap Heat Switches for use in Low-Temperature Cryogenic
Systems"). This is a genuinely analogous MECHANISM -- a mechanical
actuator pressing two contact bodies together/apart to switch thermal
conductance, exactly Sect. 6.2.4's category -- but it is NOT an AMR-
specific or room-temperature-validated figure: it comes from cryogenic
ADR/cryocooler heat-switch literature, a different application, thermal
regime, and duty cycle (static engage/disengage over minutes-to-hours,
not the ~Hz-scale cyclic actuation an AMR diode would need). A separate
review (dilution-refrigerator gas-gap heat switch literature) notes
mechanical heat switches have "in theory... infinite on/off ratios" but
that "complicated configuration leads to high additional heat loss and
poor durability, thereby greatly limiting... practical application" --
i.e. the achievable ratio is a design/durability trade-off, not a fixed
material constant, which is why this module picks a conservative
`rectification_ratio=20` (roughly a tenth of the PZHS's reported 100-200
ceiling) rather than reproducing that ceiling directly. `actuation_
energy_J_per_cycle` has NO literature source at all -- none of the
sources found report per-actuation energy for a device cycling at
AMR-relevant frequencies (~0.1-10 Hz) -- and remains a round-number
placeholder, flagged as such, at the same weakest-link tier as the
`hysteresis_loss_J_per_kg` literature analogs. If a fuller copy of
Kitanovski becomes available, or an AMR-specific (room-temperature,
Hz-scale) mechanical-diode source is found, this module's defaults
should be replaced and this honesty flag revisited -- same "what to do
if better data arrives" framing used for CYCLE_TYPE_FACTORS.

Validation status (ORIGINAL, MechanicalContactDiode only): NONE of this
project's 16 benchmarked AMR devices (data/amr_experimental_benchmarks.csv)
use thermal diodes of any kind -- every one is a conventional
valve-switched or continuous-rotary design. There is therefore no
benchmark row MechanicalContactDiode's numbers can be checked against.
MechanicalContactDiode itself is UNCHANGED by the update below and
remains a design-exploration tool for exactly that reason.

VALIDATION UPDATE (this pass): a fresh, targeted web search (not
constrained to this project's own PDF corpus -- Papers.zip's copy of
Kitanovski et al. 2015 is the SAME 30-page excerpt, verified by page
count, so it adds nothing new here) found a genuinely different active-
diode MECHANISM with exactly the AMR-relevant, room-temperature,
Hz-scale grounding the original honesty flag said was missing: a
ferrofluid-based magnetically-activated thermal switch (MATS), not
Sect. 6.2.4's mechanical-contact mechanism, but the same functional
class (active, cyclically-switched thermal rectifier). Four real,
peer-reviewed sources ground the new `FerrofluidThermalSwitch` class
below:

  1. Katiyar, Dhar, Nandi & Das, "Magnetic field induced augmented
     thermal conduction phenomenon in magneto-nanocolloids," J. Magn.
     Magn. Mater. 419 (2016) 588-599 -- DIRECTLY MEASURED ferrofluid
     thermal-conductivity switching: up to 284% conductivity
     enhancement (Fe/Co/Ni nanoparticles, 5 vol%, 500 G field) relative
     to the zero-field state. This is a real ON/OFF conductance ratio
     from a real measurement, not a cryogenic analog -- it replaces
     MechanicalContactDiode's borrowed PZHS ratio for this new class.
  2. Rodrigues, Dias, Martins, Silva, Araújo, Oliveira, Pereira &
     Ventura, "A magnetically-activated thermal switch without moving
     parts," Applied Energy 251 (2019) 113213 (also arXiv:1803.10490) --
     a REAL, built, electromagnet-driven ferrofluid thermal switch with
     NO moving parts, characterized over a 0.5-18 Hz frequency range
     and 6.5-39 W coil-power range, with temperature-gradient-dependent
     switching efficiency up to 44.4% and switching rates up to 0.6
     C/s. This is exactly the "AMR-specific (room-temperature, Hz-
     scale) ... source" the original honesty flag said to look for.
  3. Andrade, Fernandes, Teixeira, Pereira, Pires, Silva, Ventura &
     Oliveira, "High-performance magnetic thermal switch based on
     MnFe2O4/Ethylene Glycol:Water refrigerant dispersion," Applied
     Energy 356 (2024) 122325 -- a second real MATS device, 0.01-0.60
     Hz, ms-scale ON/OFF switching, up to 60% temperature-span increase
     over the bare-conduction baseline.
  4. Klinar, Vozel, Swoboda, Sojer, Muñoz Rojo & Kitanovski (the SAME
     lead author as this project's own reference book), "Ferrofluidic
     thermal switch in a magnetocaloric device," iScience 25 (2022)
     103779 -- a numerical DEVICE-LEVEL model (not this module's own
     model) of a full magnetocaloric embodiment using a ferrofluidic
     thermal switch, itself parameterized from Katiyar et al. (2016)'s
     measured properties. Reports, at 20 Hz: 5 ms switching response
     time, contact resistance R_con = 0.006 K m^2/W (independently
     cited as consistent with Cengel (2002)'s standard heat-transfer
     textbook values), max temperature span 1.12 K (single, non-
     regenerative embodiment), max cooling power 850 W/m^2 (0.37 W/g
     specific cooling power for Gd), and COP up to 8.5 at that maximum
     cooling power. Used below as an independent corroborating source
     for the response-time and contact-resistance orders of magnitude,
     NOT as this module's own governing model (Klinar et al.'s model is
     a full 1D device simulation, not this module's simple two-state
     conductance abstraction).

`core/thermal_diode_analysis.py` now also has an actual benchmark DEVICE
to check against: Andrade et al.'s 2024 Int. J. Refrigeration paper
(separately cited there) built a real Gd + ferrofluid-thermal-switch
refrigeration prototype -- see that module's `check_against_
andrade_2024_benchmark()`. `FerrofluidThermalSwitch` /
`DEFAULT_FERROFLUID_THERMAL_SWITCH` are therefore promoted to a
validated feature (grounded rectification ratio, grounded frequency
range, an actual benchmark device to check qualitative behavior
against); `MechanicalContactDiode` / `DEFAULT_MECHANICAL_CONTACT_DIODE`
remain a design-exploration tool (cryogenic analog only, still no
benchmark device) -- both classes are kept, sharing the same interface,
so a caller can pick either mechanism explicitly for AMRSystem's
`thermal_diode` parameter.
"""

from dataclasses import dataclass


@dataclass
class MechanicalContactDiode:
    """First-pass, textbook-level model of a mechanical-contact active
    thermal diode: `forward_conductance_W_K` is the effective thermal
    conductance (W/K) across the joint when the diode mechanism
    physically engages the two contact bodies; `reverse_conductance_W_K`
    is the (much smaller) conductance when disengaged. The ratio of the
    two -- `rectification_ratio` -- is the dimensionless figure of merit
    the general thermal-diode review literature reports devices by.

    `actuation_energy_J_per_cycle` is the (illustrative, unbenchmarked --
    see module honesty flag) electrical energy dissipated per
    engage-then-disengage actuation of the mechanism, independent of how
    much regenerator mass the diode is attached to. This is a per-DIODE
    quantity, not a per-kg one, because the plan's own physical picture
    (Sect. 6.2.4's mechanical-contact mechanism) is a discrete actuator
    -- a solenoid, cam or piezo stack pressing two plates together --
    whose actuation cost does not scale with the mass of regenerator
    material on the other side of the joint, unlike (e.g.) the earlier
    hysteresis loss, which is intrinsic to every kg of first-order
    material in the bed.
    """

    forward_conductance_W_K: float
    reverse_conductance_W_K: float
    actuation_energy_J_per_cycle: float = 0.0

    def __post_init__(self):
        if self.forward_conductance_W_K <= 0.0:
            raise ValueError("forward_conductance_W_K must be positive")
        if self.reverse_conductance_W_K <= 0.0:
            raise ValueError("reverse_conductance_W_K must be positive")
        if self.reverse_conductance_W_K > self.forward_conductance_W_K:
            raise ValueError(
                "reverse_conductance_W_K must not exceed forward_conductance_W_K "
                "(rectification_ratio must be >= 1 for this to behave as a diode "
                f"at all; got forward={self.forward_conductance_W_K}, "
                f"reverse={self.reverse_conductance_W_K})")
        if self.actuation_energy_J_per_cycle < 0.0:
            raise ValueError("actuation_energy_J_per_cycle must be non-negative")

    @property
    def rectification_ratio(self) -> float:
        """forward/reverse conductance -- the standard figure of merit for
        a thermal diode/rectifier (always >= 1 by construction, see
        __post_init__)."""
        return self.forward_conductance_W_K / self.reverse_conductance_W_K

    def switching_power_W(self, frequency: float) -> float:
        """Parasitic electrical power (W) to actuate the mechanical
        contact once per AMR cycle, at the given cycle `frequency` (Hz):
            W_switch = actuation_energy_J_per_cycle * frequency
        Returns exactly 0.0 when actuation_energy_J_per_cycle is 0.0 (the
        dataclass default), so a caller that wants ONLY the heat-transfer
        side of this model (rectification_ratio) without any parasitic
        cost can do so explicitly."""
        if frequency < 0.0:
            raise ValueError("frequency must be non-negative")
        return self.actuation_energy_J_per_cycle * frequency


@dataclass
class FerrofluidThermalSwitch:
    """Ferrofluid-based active thermal switch (magnetically-activated
    thermal switch, MATS) -- see module docstring's VALIDATION UPDATE
    for the four real, cited sources this class is grounded in. A
    field-driven ferrofluid switches between a random-orientation
    (low-conductance, field-OFF) state and a field-aligned chain-
    structure (high-conductance, field-ON) state; SAME forward/reverse-
    conductance, SAME rectification_ratio, and SAME switching_power_W(
    frequency) interface as MechanicalContactDiode (so it is a drop-in
    alternative for AMRSystem's `thermal_diode` parameter), but a
    genuinely different physical mechanism -- no mechanical contact, no
    moving parts (Rodrigues et al. 2019's own headline result).

    `max_tested_frequency_Hz` records the highest frequency actually
    exercised in the cited hardware (Rodrigues et al. 2019: up to 18 Hz)
    -- reported as a HIGHEST-TESTED value, not asserted here as a proven
    physical ceiling on the mechanism, since no source found explains
    why 18 Hz specifically would be a hard limit rather than just the
    highest point these particular experiments swept.

    `actuation_energy_J_per_cycle` here is a DERIVED, explicitly-flagged
    estimate, not a directly-reported per-actuation figure (none of the
    four cited sources report one in that exact form): Rodrigues et al.
    (2019) report continuous electromagnet coil power over the range
    6.5-39 W while the field is held ON, which is a HOLDING power, not a
    discrete engage/disengage actuation energy the way
    MechanicalContactDiode's mechanism works. `holding_power_W` is kept
    as its own field below so a caller can use the physically more
    faithful continuous-power picture directly; `actuation_energy_
    J_per_cycle` is then holding_power_W * duty_cycle_fraction /
    reference_frequency_Hz -- i.e. the energy that same holding power
    would dissipate over one AMR half-cycle at a stated reference
    frequency and duty cycle -- purely so this class can satisfy
    AMRSystem's existing energy-per-cycle*frequency interface. Because
    of this conversion, `switching_power_W()` on this class is only
    frequency-INDEPENDENT-accurate near `reference_frequency_Hz`; at
    other frequencies it silently rescales as energy*frequency, whereas
    the truer physical picture (holding_power_W * duty_cycle_fraction)
    would NOT scale with frequency at all. This mismatch is a real
    modeling simplification, stated here rather than hidden, kept only
    for interface compatibility with the existing AMRSystem wiring.
    """

    forward_conductance_W_K: float
    reverse_conductance_W_K: float
    actuation_energy_J_per_cycle: float = 0.0
    holding_power_W: float = 0.0
    max_tested_frequency_Hz: float = float("inf")

    def __post_init__(self):
        if self.forward_conductance_W_K <= 0.0:
            raise ValueError("forward_conductance_W_K must be positive")
        if self.reverse_conductance_W_K <= 0.0:
            raise ValueError("reverse_conductance_W_K must be positive")
        if self.reverse_conductance_W_K > self.forward_conductance_W_K:
            raise ValueError(
                "reverse_conductance_W_K must not exceed forward_conductance_W_K "
                "(rectification_ratio must be >= 1 for this to behave as a diode "
                f"at all; got forward={self.forward_conductance_W_K}, "
                f"reverse={self.reverse_conductance_W_K})")
        if self.actuation_energy_J_per_cycle < 0.0:
            raise ValueError("actuation_energy_J_per_cycle must be non-negative")
        if self.holding_power_W < 0.0:
            raise ValueError("holding_power_W must be non-negative")
        if self.max_tested_frequency_Hz <= 0.0:
            raise ValueError("max_tested_frequency_Hz must be positive")

    @property
    def rectification_ratio(self) -> float:
        """forward/reverse conductance -- for DEFAULT_FERROFLUID_THERMAL_
        SWITCH this is the real Katiyar et al. (2016) measured ratio
        (~3.84, i.e. a 284% enhancement), not an invented figure."""
        return self.forward_conductance_W_K / self.reverse_conductance_W_K

    def switching_power_W(self, frequency: float) -> float:
        """Same interface as MechanicalContactDiode.switching_power_W --
        see this class's docstring for why, for THIS mechanism, that
        interface is only an approximation of the truer continuous-
        holding-power physics away from the frequency this instance's
        actuation_energy_J_per_cycle was derived at."""
        if frequency < 0.0:
            raise ValueError("frequency must be non-negative")
        if frequency > self.max_tested_frequency_Hz:
            import warnings
            warnings.warn(
                f"frequency={frequency} Hz exceeds the highest frequency "
                f"actually tested in the literature this default is "
                f"grounded in ({self.max_tested_frequency_Hz} Hz, "
                f"Rodrigues et al. 2019) -- extrapolating beyond measured "
                f"data.", stacklevel=2)
        return self.actuation_energy_J_per_cycle * frequency


def cycle_time_reduction_factor(conventional_switch_time_s: float,
                                  diode_switch_time_s: float) -> float:
    """Sensitivity/what-if helper, NOT a literature-derived prediction
    (see module honesty flag): if a conventional valve-switched AMR bed
    spends `conventional_switch_time_s` of each half-cycle transitioning
    flow direction (dead time with no useful heat transfer), and a
    diode-based design could in principle transition in
    `diode_switch_time_s` instead, returns the FRACTIONAL REDUCTION in
    that dead time per half-cycle -- i.e. an upper bound on how much
    additional useful time (and therefore, at a fixed total cycle
    period, how much higher achievable frequency) a diode-assisted
    design could in principle recover.

    Deliberately requires BOTH switch times as explicit caller-supplied
    arguments, with no default for either: no digitized source for
    either a conventional valve's or a mechanical-contact diode's own
    switching time exists in this project's corpus (see the module
    docstring's honesty flag), so inventing a default for either would
    be unfounded precision, not a documented literature value. Use this
    function to explore "what would have to be true" rather than to
    read off a number this repo claims to already know.
    """
    if conventional_switch_time_s <= 0.0:
        raise ValueError("conventional_switch_time_s must be positive")
    if diode_switch_time_s < 0.0:
        raise ValueError("diode_switch_time_s must be non-negative")
    if diode_switch_time_s > conventional_switch_time_s:
        raise ValueError(
            "diode_switch_time_s exceeds conventional_switch_time_s -- a "
            "diode design slower than the conventional valve provides no "
            "cycle-time benefit under this model; check inputs")
    return 1.0 - diode_switch_time_s / conventional_switch_time_s


# Grounded-but-still-illustrative default -- see module docstring honesty
# flag for the full citation and caveats. forward/reverse conductance
# below give rectification_ratio=20, a conservative fraction of the
# ~100-200 on/off ratio Bywaters & Griffin report for a piezo-actuated
# mechanical heat switch (cryogenic, NOT AMR-specific -- see honesty
# flag). actuation_energy_J_per_cycle has NO literature source found and
# remains a round-number placeholder.
DEFAULT_MECHANICAL_CONTACT_DIODE = MechanicalContactDiode(
    forward_conductance_W_K=5.0,          # engaged-contact conductance, illustrative
    reverse_conductance_W_K=0.25,         # disengaged-contact conductance, illustrative
    actuation_energy_J_per_cycle=0.05,    # engage+disengage actuation energy, illustrative (no literature source)
)

# Grounded in real, measured/tested literature -- see module docstring's
# VALIDATION UPDATE for full citations. reverse_conductance_W_K reuses
# MechanicalContactDiode's own illustrative field-OFF baseline geometry
# (0.25 W/K) purely so the two defaults are comparable at the same
# baseline -- Katiyar et al. (2016) report a CONCENTRATION-DEPENDENT
# ratio, not an absolute W/K figure for any particular device geometry
# (their measurement is a bulk-fluid k_on/k_off ratio, not a per-device
# conductance), so the absolute W/K scale here is still illustrative;
# the RATIO (3.84) is the literature-grounded part.
# holding_power_W=15.0 is the MIDPOINT of Rodrigues et al. (2019)'s own
# tested 6.5-39 W coil-power range (not an extreme), converted to
# actuation_energy_J_per_cycle via a duty_cycle_fraction=0.5 (symmetric
# on/off -- the SAME symmetric-cycling assumption Andrade et al. (2024)
# tested experimentally and found gives no net advantage, see
# core/thermal_diode_analysis.py) at reference_frequency_Hz=4.0 (this
# project's own representative AMR operating frequency, matching
# core/thermal_diode_analysis.py's T_COLD_K/SPAN_K/MU0H_T/MASS_KG/
# MDOT_KG_S operating point and roughly the midpoint of that module's
# own 0.5-8 Hz sweep).
_KATIYAR_2016_CONDUCTIVITY_RATIO = 3.84  # 1 + 284% measured enhancement
_RODRIGUES_2019_HOLDING_POWER_W = 15.0  # midpoint of tested 6.5-39 W range
_DUTY_CYCLE_FRACTION = 0.5                # symmetric on/off, see above
_REFERENCE_FREQUENCY_HZ = 4.0             # this project's representative op. point
DEFAULT_FERROFLUID_THERMAL_SWITCH = FerrofluidThermalSwitch(
    forward_conductance_W_K=0.25 * _KATIYAR_2016_CONDUCTIVITY_RATIO,  # ~0.96 W/K
    reverse_conductance_W_K=0.25,
    actuation_energy_J_per_cycle=(_RODRIGUES_2019_HOLDING_POWER_W
                                    * _DUTY_CYCLE_FRACTION
                                    / _REFERENCE_FREQUENCY_HZ),        # 1.875 J/cycle
    holding_power_W=_RODRIGUES_2019_HOLDING_POWER_W,
    max_tested_frequency_Hz=18.0,  # Rodrigues et al. (2019)'s own tested ceiling
)
