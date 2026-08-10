"""
fluid_mce_cycle.py
====================
Phase 20: magnetocaloric fluids (ferrofluid / magnetorheological particle
suspension) as an alternative working-body class, motivated by
ROADMAP.md's Phase 20 plan item ("Magnetocaloric fluids as an
alternative working-body class").

STRUCTURAL SCOPE (per the plan's own scoping decision, repeated here for
anyone reading only this file): this is a NEW SIBLING to
core.amr_cycle.AMRSystem, `FerrofluidMCESystem`, not a parameter bolted
onto AMRSystem. A magnetocaloric fluid flows continuously and is
magnetized/demagnetized IN PLACE as it passes through hot- and cold-side
heat exchangers -- there is no separate solid regenerator bed and no
separate heat-transfer fluid; the working substance and the heat-
transfer medium are the same continuously-flowing suspension. Forcing
that through AMRSystem's packed-bed/parallel-plate geometry assumptions
(core/thermal.py) would misrepresent the architecture, per the plan's
own recommendation.

HONESTY FLAG #1 (book access). This project's copy of Kitanovski et al.
(2015) is a 30-page front-matter/Chapter-1-only excerpt -- it does NOT
include Chapter 5 (Magnetocaloric Fluids, pp. 167-209), where Sect. 5.2
(ferrofluid/magnetorheological rheology), Sect. 5.4/5.5 (device design
notes) live, the sections the Phase 20 plan itself named as the source
for "rheology-adjusted relations." None of that chapter's own equations
or coefficients could be digitized here. This project's copy of Tishin &
Spichkin (2003) is a scanned, image-only 486-page PDF with NO extractable
text layer (`pypdf` returns zero characters from every page checked in
this pass) -- its actual chapter contents could not be searched or
digitized either, so the Phase 20 plan's own claim that "Tishin doesn't
cover this (its Ch. 11 passive-regenerator focus is solid-state)" is
taken at face value here, UNVERIFIED by this pass, not independently
confirmed.

What follows instead is a first-principles model built from three
STANDARD, textbook physical relations, none of which require either
book:
  1. A mixture (lumped) heat-capacity dilution argument for the
     suspension's effective adiabatic temperature change (see
     `suspension_delta_T_adiabatic()` below) -- ordinary calorimetry,
     not a magnetocaloric-specific result.
  2. The Krieger-Dougherty equation for suspension viscosity (Krieger &
     Dougherty, Trans. Soc. Rheol. 3 (1959) 137-152) -- a standard,
     widely-used rheology correlation, generically attributable rather
     than digitized from either book.
  3. The Darcy-Weisbach equation for pipe pressure drop (standard
     textbook fluid mechanics), used instead of core/thermal.py's
     packed-bed correlations because a continuous-flow fluid loop has no
     packed bed -- see `pumping_power_pipe_flow()`.

HONESTY FLAG #2 (benchmark availability, found during this pass's
literature search -- web search only, no full-text fetch or
digitization of any paper found). No clean, directly-usable Qc/COP
benchmark for a magnetocaloric-FLUID-as-WORKING-BODY refrigeration
device was found in this project's corpus or in this pass's search
(data/amr_experimental_benchmarks.csv is solid-AMR only). Two ADJACENT,
but structurally DIFFERENT, results surfaced and are deliberately NOT
conflated with this module's own topic:
  * Andrade, Fernandes, Silva, Teixeira & Pereira, "Magnetic
    refrigeration enhanced by magnetically-activated thermal switch: an
    experimental proof-of-concept" (2024) couples a SOLID Gd regenerator
    with a ferrofluid used as a THERMAL SWITCH (a heat-transfer-
    enhancement/rectification role, closer in spirit to Phase 18's
    thermal-diode item than to this module's working-body topic) -- its
    own abstract reports no COP advantage from the ferrofluid over plain
    Gd under symmetric cycling.
  * A gallium-based magnetocaloric liquid-metal ferrofluid study (2017)
    tests cooling performance on an electric transformer coil, not a
    refrigeration cycle with a reported Qc/COP figure comparable to this
    repo's benchmark convention.
  Neither is used as a calibration point here. This module is therefore
  explicitly a DESIGN-EXPLORATION / comparison tool, not a validated
  feature -- the same disposition Phase 18 gave
  core/thermal_diode_analysis.py.

Physics used here
-------------------
Suspension = ferromagnetic particles (volume fraction phi) dispersed in
a non-magnetic carrier liquid (default: water, reusing
core.thermal.water_properties()).

1. Adiabatic temperature change (dilution model). When the particle
   phase's own field-dependent entropy change releases/absorbs heat,
   that heat is shared across the ENTIRE suspension's heat capacity
   (particles + carrier), not just the particles' own mass -- ordinary
   mixture calorimetry, not a magnetocaloric-specific assumption:

       dTad_suspension = dTad_pure(T, H) * (phi*rho_p*cp_p) /
                          (phi*rho_p*cp_p + (1-phi)*rho_c*cp_c)

   This is a MONOTONICALLY INCREASING function of phi (dTad_suspension
   -> dTad_pure as phi -> 1, -> 0 as phi -> 0) -- the "MCE-intensity-vs-
   phi" side of the tradeoff the Phase 20 plan named.
2. Suspension viscosity (Krieger-Dougherty): mu_suspension = mu_carrier *
   (1 - phi/phi_max)^(-[eta]*phi_max), with phi_max=0.63 (random close
   packing) and intrinsic viscosity [eta]=2.5 (Einstein value for rigid
   spheres) as defaults -- the "viscosity-vs-phi" side of the tradeoff,
   diverging as phi -> phi_max.
3. No regeneration. Unlike AMRSystem (where a solid regenerator bed lets
   the achievable span exceed a single stage's own dTad via regenerative
   amplification -- see amr_cycle.py's own characteristic-curve
   discussion), a single-pass continuously-flowing fluid loop has no
   regenerator: the full imposed span must fit within ONE adiabatic
   swing. `cooling_capacity()` below therefore uses
   `span_fraction = max(0, 1 - T_span/dTad_suspension)` (a full-swing
   denominator), NOT AMRSystem's regenerated `1 - T_span/(2*dTad)`
   half-swing denominator -- see that function's own docstring for why
   this is a genuine structural difference, not a copy-paste of
   AMRSystem's formula.

Limitations, stated rather than hidden
-----------------------------------------
* `eta_2nd_law_fluid` (see `FerrofluidMCESystem.__init__`) is a fixed,
  illustrative, UNCALIBRATED second-law efficiency (no benchmark device
  exists to calibrate it against -- see HONESTY FLAG #2), unlike
  AMRSystem's own `eta_2nd_law` term, which cites two specific literature
  devices.
* No demagnetizing-field correction is applied to the particles' own
  field-dependent MCE inside the suspension (a real effect for magnetic
  nanoparticle suspensions, not modeled here).
* The representative pipe-loop geometry (`pipe_diameter_m`,
  `pipe_length_m`) is illustrative, not device-specific -- see
  `pumping_power_pipe_flow()`'s own docstring.
* Brownian/thermal relaxation of particle magnetic moments (superpara-
  magnetic behavior at small particle sizes) is not modeled; particles
  are treated as if they carry the SAME field-dependent entropy change
  as the corresponding bulk MagnetocaloricMaterial, scaled only by mass
  fraction -- a simplification, not a nanoparticle-physics model.
"""
import numpy as np
from dataclasses import dataclass
from core.mce_material import MagnetocaloricMaterial

DEFAULT_PHI_MAX = 0.63       # random close packing, standard value
DEFAULT_INTRINSIC_VISCOSITY = 2.5  # Einstein value for rigid spheres


def krieger_dougherty_viscosity(mu_carrier, phi, phi_max=DEFAULT_PHI_MAX,
                                  intrinsic_viscosity=DEFAULT_INTRINSIC_VISCOSITY):
    """Krieger & Dougherty (1959) suspension viscosity relation. Diverges
    (mu -> infinity) as phi -> phi_max, a genuine physical ceiling on
    achievable particle loading this function does not clip -- callers
    sweeping phi should expect rapidly growing viscosity/pumping power
    near phi_max, not a wall."""
    if not (0 <= phi < phi_max):
        raise ValueError(f"phi must satisfy 0 <= phi < phi_max ({phi_max})")
    return mu_carrier * (1 - phi / phi_max) ** (-intrinsic_viscosity * phi_max)


def suspension_effective_properties(phi, rho_particle=5180.0, cp_particle=670.0,
                                      carrier="water", T_K=300.0):
    """Volume-fraction-weighted suspension density/heat-capacity, plus the
    dilution factor used by `suspension_delta_T_adiabatic()`.
    `rho_particle`/`cp_particle` default to magnetite (Fe3O4)'s standard
    reference values (5180 kg/m^3, ~670 J/(kg K)) -- the most common
    ferrofluid particle material -- not a claim about any specific
    magnetocaloric particle composition; a caller modeling a different
    particle material should pass its own values."""
    from core.thermal import water_properties
    if carrier == "water":
        fluid = water_properties(T_K)
    else:
        raise ValueError(f"unknown carrier {carrier!r}; only 'water' is implemented")
    rho_c, cp_c, mu_c = fluid["rho"], fluid["cp"], fluid["mu"]

    rho_susp = phi * rho_particle + (1 - phi) * rho_c
    particle_heat_capacity = phi * rho_particle * cp_particle
    carrier_heat_capacity = (1 - phi) * rho_c * cp_c
    total_heat_capacity = particle_heat_capacity + carrier_heat_capacity
    cp_susp = total_heat_capacity / rho_susp if rho_susp > 0 else 0.0
    dilution_factor = (particle_heat_capacity / total_heat_capacity
                        if total_heat_capacity > 0 else 0.0)
    return {
        "rho_susp_kg_m3": rho_susp,
        "cp_susp_J_kgK": cp_susp,
        "mu_carrier_Pa_s": mu_c,
        "dilution_factor": dilution_factor,
    }


def suspension_delta_T_adiabatic(material: MagnetocaloricMaterial, T, mu0H_final,
                                   phi, rho_particle=5180.0, cp_particle=670.0,
                                   mu0H_initial=0.0, carrier="water", T_K=300.0):
    """Eq. in module docstring item 1: the suspension's own effective
    adiabatic temperature change, diluted from `material`'s pure-material
    dTad by the mixture heat-capacity dilution_factor from
    `suspension_effective_properties()`. Monotonically increasing in
    phi -- see module docstring."""
    mu0 = 4 * np.pi * 1e-7
    H_final = mu0H_final / mu0
    H_initial = mu0H_initial / mu0
    dTad_pure = float(material.delta_T_adiabatic(np.array([T]), H_final, H_initial)[0])
    props = suspension_effective_properties(phi, rho_particle, cp_particle, carrier, T_K)
    return dTad_pure * props["dilution_factor"]


def pumping_power_pipe_flow(mdot, mu_susp, rho_susp, pipe_diameter_m=0.01,
                              pipe_length_m=1.0):
    """Darcy-Weisbach pipe pressure drop -> hydraulic pumping power (W),
    standard textbook fluid mechanics (module docstring item 3), used
    instead of core/thermal.py's packed-bed correlations since a
    continuous-flow fluid loop has no packed bed. Friction factor: 64/Re
    for laminar flow (Re<2300), Blasius correlation
    f=0.316*Re^-0.25 for turbulent flow (Re>=2300) -- both standard,
    textbook, not requiring either project book. `pipe_diameter_m`/
    `pipe_length_m` are a representative, illustrative loop geometry, NOT
    device-specific (see module docstring limitations)."""
    if mdot <= 0:
        return 0.0
    area = np.pi * (pipe_diameter_m / 2) ** 2
    velocity = mdot / (rho_susp * area)
    Re = rho_susp * velocity * pipe_diameter_m / mu_susp if mu_susp > 0 else np.inf
    if Re < 2300:
        f = 64 / max(Re, 1e-6)
    else:
        f = 0.316 * Re ** -0.25
    delta_P = f * (pipe_length_m / pipe_diameter_m) * (rho_susp * velocity ** 2 / 2)
    volumetric_flow = mdot / rho_susp
    return delta_P * volumetric_flow  # W = Pa * m^3/s


@dataclass
class FluidMCEResult:
    T_span: float
    Qc: float
    Qh: float
    W_mag: float
    W_parasitic: float
    COP: float
    COP_electrical: float
    exergy_eff: float
    dTad_suspension_K: float
    dilution_factor: float
    viscosity_Pa_s: float


class FerrofluidMCESystem:
    """New sibling to core.amr_cycle.AMRSystem -- see module docstring
    for why this is a separate class rather than a parameter on
    AMRSystem."""

    def __init__(self, material: MagnetocaloricMaterial, mu0H_max: float,
                 particle_volume_fraction: float, fluid_mdot: float,
                 rho_particle: float = 5180.0, cp_particle: float = 670.0,
                 pipe_diameter_m: float = 0.01, pipe_length_m: float = 1.0,
                 eta_2nd_law_fluid: float = 0.20,
                 phi_max: float = DEFAULT_PHI_MAX,
                 intrinsic_viscosity: float = DEFAULT_INTRINSIC_VISCOSITY,
                 carrier: str = "water", T_K: float = 300.0):
        """
        material                   : MagnetocaloricMaterial (bulk particle
                                      material properties -- see module
                                      docstring limitation on Brownian/
                                      superparamagnetic effects)
        mu0H_max                    : peak applied field, Tesla
        particle_volume_fraction    : phi, 0 <= phi < phi_max
        fluid_mdot                  : suspension mass flow rate, kg/s
        rho_particle, cp_particle   : particle material density/specific
                                      heat; default to magnetite (Fe3O4)
                                      reference values -- see
                                      `suspension_effective_properties()`
        pipe_diameter_m, pipe_length_m : representative, illustrative
                                      pipe-loop geometry for the Darcy-
                                      Weisbach pumping-power term (see
                                      module docstring limitation -- NOT
                                      device-specific)
        eta_2nd_law_fluid           : fixed, UNCALIBRATED second-law
                                      efficiency for the non-regenerated
                                      Brayton-like fluid cycle (see
                                      module docstring HONESTY FLAG #2 --
                                      no benchmark device exists to
                                      calibrate this against, unlike
                                      AMRSystem's own eta_2nd_law, which
                                      cites two literature devices).
                                      Default 0.20 is deliberately LOWER
                                      than AMRSystem's regenerated
                                      0.35-0.52 range, reflecting the
                                      lack of regeneration (module
                                      docstring physics item 3) -- an
                                      illustrative choice, not a fitted
                                      one.
        phi_max, intrinsic_viscosity : Krieger-Dougherty parameters, see
                                      `krieger_dougherty_viscosity()`
        carrier                     : carrier fluid, only "water"
                                      implemented (reuses
                                      core.thermal.water_properties())
        T_K                         : carrier fluid property evaluation
                                      temperature
        """
        if not (0 <= particle_volume_fraction < phi_max):
            raise ValueError(
                f"particle_volume_fraction must satisfy 0 <= phi < phi_max ({phi_max})")
        self.mat = material
        self.mu0H_max = mu0H_max
        self.phi = particle_volume_fraction
        self.mdot = fluid_mdot
        self.rho_particle = rho_particle
        self.cp_particle = cp_particle
        self.pipe_diameter_m = pipe_diameter_m
        self.pipe_length_m = pipe_length_m
        self.eta_2nd_law_fluid = eta_2nd_law_fluid
        self.phi_max = phi_max
        self.intrinsic_viscosity = intrinsic_viscosity
        self.carrier = carrier
        self.T_K = T_K

    def _suspension_props(self):
        return suspension_effective_properties(
            self.phi, self.rho_particle, self.cp_particle, self.carrier, self.T_K)

    def cooling_capacity(self, T_cold, T_span):
        """See module docstring physics item 3 for why this uses a
        FULL-swing denominator (1 - T_span/dTad), not AMRSystem's
        regenerated half-swing (1 - T_span/(2*dTad)) -- a genuine
        structural consequence of having no regenerator, not an
        arbitrary choice."""
        T_hot = T_cold + T_span
        T_mid = 0.5 * (T_cold + T_hot)
        dTad = suspension_delta_T_adiabatic(
            self.mat, T_mid, self.mu0H_max, self.phi, self.rho_particle,
            self.cp_particle, carrier=self.carrier, T_K=self.T_K)
        if dTad <= 0:
            return 0.0, dTad
        span_fraction = max(0.0, 1.0 - T_span / dTad)
        props = self._suspension_props()
        Qc = self.mdot * props["cp_susp_J_kgK"] * dTad * span_fraction
        return max(Qc, 0.0), dTad

    def magnetic_work(self, T_cold, T_span, Qc):
        """Non-regenerated Brayton-like reference work -- same second-law
        decomposition structure as AMRSystem.magnetic_work(), but with
        the fixed, illustrative `eta_2nd_law_fluid` instead of
        AMRSystem's effectiveness-dependent formula (there is no
        regenerator effectiveness in this architecture -- see module
        docstring)."""
        T_hot = T_cold + T_span
        carnot_work = Qc * (T_hot / T_cold - 1.0) if T_cold > 0 else np.inf
        eta = float(np.clip(self.eta_2nd_law_fluid, 0.02, 0.95))
        W = carnot_work / max(eta, 1e-3)
        return W, eta

    def run(self, T_cold, T_span) -> FluidMCEResult:
        Qc, dTad = self.cooling_capacity(T_cold, T_span)
        W, eta2 = self.magnetic_work(T_cold, T_span, Qc)
        props = self._suspension_props()
        mu_susp = krieger_dougherty_viscosity(
            props["mu_carrier_Pa_s"], self.phi, self.phi_max, self.intrinsic_viscosity)
        W_parasitic = pumping_power_pipe_flow(
            self.mdot, mu_susp, props["rho_susp_kg_m3"],
            self.pipe_diameter_m, self.pipe_length_m)
        Qh = Qc + W
        COP = Qc / W if W > 0 else 0.0
        COP_electrical = Qc / (W + W_parasitic) if (W + W_parasitic) > 0 else 0.0
        T_hot = T_cold + T_span
        COP_carnot = T_cold / (T_hot - T_cold) if T_hot > T_cold else np.inf
        exergy_eff = COP / COP_carnot if np.isfinite(COP_carnot) and COP_carnot > 0 else 0.0
        return FluidMCEResult(
            T_span=T_span, Qc=Qc, Qh=Qh, W_mag=W, W_parasitic=W_parasitic,
            COP=COP, COP_electrical=COP_electrical, exergy_eff=exergy_eff,
            dTad_suspension_K=dTad, dilution_factor=props["dilution_factor"],
            viscosity_Pa_s=mu_susp)

    def characteristic_curve(self, T_cold, spans):
        return [self.run(T_cold, s) for s in spans]