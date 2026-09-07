"""
barocaloric_material.py
========================
Clausius-Clapeyron model of the (colossal) barocaloric effect in the
plastic crystal neopentylglycol (NPG), structured the same way this repo
models other first-order caloric effects (see elastocaloric_material.py's
docstring for the shared pattern, and mce_material.py for the magnetic
analogue).

Physics
-------
NPG undergoes a pressure-driven order-disorder solid-solid phase
transition. As for the elastocaloric case, the transformed fraction under
an applied pressure is approximated with this repo's own documented
linear-clamp approximation (see amr_cycle.py's span_fraction and
elastocaloric_material.py's transformed_fraction -- same functional form,
different driving variable):
    x(P) = clip(P / P_sat, 0, 1)
    Delta_T_ad(P) = T * Delta_S * x(P) / c_p

Parameter sources
------------------
  - Isothermal entropy change |Delta_S| ~ 389 J/(kg K) at an applied
    pressure of 45.0 MPa -- the original "colossal barocaloric effect"
    measurement on NPG (Li, B. et al., Nature 567 (2019) 506-510, as
    reported via NSFC's own English-language research highlight, which
    this repo's search process could access directly; the original Nature
    paper itself was not available to this project). This is an order of
    magnitude larger than typical magnetocaloric/electrocaloric entropy
    changes -- the actual reason plastic crystals attracted interest as
    VCC replacements.
  - Transition temperature T0 ~ 308 K (NPG's solid I/solid II transition
    near 40 C at ambient pressure; commonly cited across the NPG
    barocaloric literature, e.g. the PMC-hosted "Colossal barocaloric
    effects near room temperature in plastic crystals of neopentylglycol"
    paper).
  - Saturation pressure P_sat: this repo uses 45 MPa, i.e. the SAME
    pressure at which the 389 J/(kg K) figure was measured -- meaning
    x(P_sat)=1 by construction at that specific point. This is a modeling
    choice flagged explicitly: it is NOT an independently-sourced
    "transformation complete" pressure, unlike sigma_sat in
    elastocaloric_material.py. Treat Delta_T_ad predictions here as
    anchored correctly AT 45 MPa and only qualitatively extrapolated
    elsewhere.
  - Specific heat c_p ~ 2200 J/(kg K) for NPG in its plastic-crystal phase
    (representative organic-solid value; NPG-specific calorimetry citing
    this exact figure was not directly available to this project -- FLAGGED
    as an approximation, following this repo's own convention of flagging
    where a source could not be directly verified, e.g. LIMITATIONS.md's
    treatment of image-only/inaccessible references).
  - Pressure hysteresis width Delta_P_hys ~ 15-30 MPa is commonly reported
    for pure NPG (colossal-BC papers report hysteresis narrowing
    dramatically -- to near zero -- at higher applied pressure in some
    preparations; this repo uses 20 MPa as a representative, NOT
    zero-hysteresis, midpoint appropriate for a realistic device rather
    than the best-case reversible limit some papers highlight).

HONESTY FLAG: unlike elastocaloric_material.py, this module rests on a
single directly-verified entropy-change data point (389 J/kg K @ 45 MPa)
plus several representative-but-not-directly-verified constants (c_p,
hysteresis width, T0). It should be read as substantially less certain
than the elastocaloric module, and is treated that way in
alternative_caloric_comparison.py's reported uncertainty.
"""

import numpy as np
from dataclasses import dataclass


@dataclass
class BarocaloricMaterial:
    name: str
    T0: float
    delta_S_J_per_kgK: float   # measured AT P_sat, see docstring
    P_sat_MPa: float
    hysteresis_pressure_MPa: float
    c_p_J_per_kgK: float = 2200.0
    source: str = ""

    def transformed_fraction(self, P_MPa):
        return float(np.clip(P_MPa / self.P_sat_MPa, 0.0, 1.0))

    def delta_T_ad(self, P_MPa, T=None):
        T = self.T0 if T is None else T
        x = self.transformed_fraction(P_MPa)
        return T * self.delta_S_J_per_kgK * x / self.c_p_J_per_kgK

    def hysteresis_loss_J_per_kg(self, P_MPa):
        """Dissipated work per unit mass per half-cycle from pressure
        hysteresis: hysteresis_pressure_MPa * fractional volume change.
        NPG's transition volume change is ~4-8% (commonly cited across the
        colossal-BC literature); this repo uses 5% as a representative
        midpoint, same flagging convention as elastocaloric_material.py's
        transformation_strain."""
        delta_V_over_V = 0.05
        x = self.transformed_fraction(P_MPa)
        rho = 1060.0  # NPG density, kg/m^3, representative organic-solid value
        work_density_J_per_m3 = self.hysteresis_pressure_MPa * 1e6 * delta_V_over_V * x
        return work_density_J_per_m3 / rho


NPG_plastic_crystal = BarocaloricMaterial(
    name="Neopentylglycol (NPG) plastic crystal",
    T0=308.0,
    delta_S_J_per_kgK=389.0,
    P_sat_MPa=45.0,
    hysteresis_pressure_MPa=20.0,
    source="Entropy change: Li, B. et al., Nature 567 (2019) 506-510, as "
           "reported via NSFC English-language research highlight "
           "(nsfc.gov.cn); T0 and general NPG barocaloric behavior: "
           "'Colossal barocaloric effects near room temperature in plastic "
           "crystals of neopentylglycol', PMC6472423. c_p and hysteresis "
           "width are representative approximations, NOT directly "
           "NPG-verified -- see module docstring honesty flag.",
)


def validate_against_literature(verbose=True):
    """At P_sat=45 MPa this model returns Delta_T_ad by construction from
    the measured Delta_S -- this is a consistency check on the arithmetic,
    NOT an independent validation (see docstring: P_sat was CHOSEN to be
    the measurement pressure). A genuine validation would need a second,
    independently-measured Delta_T_ad at a DIFFERENT pressure, which was
    not located during this project's search -- flagged as an open item,
    not silently skipped."""
    P = 45.0
    dT = NPG_plastic_crystal.delta_T_ad(P)
    if verbose:
        print(f"  P={P:.0f} MPa -> Delta_T_ad={dT:.2f} K "
              f"(consistency check only -- see docstring: this is NOT an "
              f"independent validation, P_sat was set equal to the "
              f"measurement pressure)")
    return dT


if __name__ == "__main__":
    print("NPG barocaloric Delta_T_ad (consistency check, see docstring caveat):")
    validate_against_literature()
