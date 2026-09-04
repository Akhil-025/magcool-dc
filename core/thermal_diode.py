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

Validation status: NONE of this project's 16 benchmarked AMR devices
(data/amr_experimental_benchmarks.csv) use thermal diodes of any kind --
every one is a conventional valve-switched or continuous-rotary design.
There is therefore no benchmark row this module's numbers can be checked
against, and `core/thermal_diode_analysis.py` (the validation
deliverable) says so explicitly rather than forcing a fit. Treat this
module as a design-exploration tool, not a validated feature -- exactly
the disposition the plan itself recommended for this item.
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