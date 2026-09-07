"""
elastocaloric_material.py
==========================
Clausius-Clapeyron model of the elastocaloric effect in NiTi (Nitinol)
shape-memory alloy, structured to parallel this repo's own
mce_material.py / first_order_mce.py: a first-order, stress-driven phase
transition (stress-induced martensite <-> austenite) plays the role that
field-driven Curie transitions play for magnetocaloric materials.

Physics
-------
NiTi's stress-induced martensitic transformation is a first-order
transition with a latent heat L and transformation entropy
    Delta_S = L / T0
Applying a mechanical stress shifts the equilibrium transformation
temperature according to the Clausius-Clapeyron relation for a
stress-driven transition:
    d(sigma)/dT = -Delta_S / Delta_epsilon_tr
where Delta_epsilon_tr is the transformation strain. Equivalently, a stress
sigma shifts the transformation temperature by
    Delta_T_eq(sigma) = sigma / (d(sigma)/dT)
Under a large-enough stress the material transforms; the transformed
fraction x(sigma) is approximated here with the SAME linear-clamp
approximation this repo already uses and documents elsewhere for an
analogous problem (see amr_cycle.py's documented `span_fraction` linear
clamp) rather than inventing a new unsourced saturation curve:
    x(sigma) = clip(sigma / sigma_sat, 0, 1)
The resulting adiabatic temperature change (small-DeltaS approximation,
same functional form as mce_material.py's DeltaT_ad):
    Delta_T_ad(sigma) = T * Delta_S * x(sigma) / c_p

Parameter sources (NiTi, superelastic, near room temperature)
---------------------------------------------------------------
  - Transformation latent heat L ~ 20-24 kJ/kg (Otsuka & Wayman, "Shape
    Memory Materials", Cambridge Univ. Press, 1998, Ch. 1; widely re-cited
    range for binary NiTi).
  - => Delta_S = L / T0 ~ 67-80 J/(kg K) at T0 ~ 300 K. This repo uses the
    midpoint 73.3 J/(kg K) (L=22000 J/kg, T0=300K) as a single
    representative value -- NOT a fitted or independently re-derived
    number, flagged the same way mce_material.py flags its own literature
    midpoints.
  - Clausius-Clapeyron slope d(sigma)/dT ~ 4-7 MPa/K is the commonly-cited
    range for binary NiTi (e.g. Otsuka & Ren, Prog. Mater. Sci. 50 (2005)
    511-678, Fig. 27 region; Qian et al., Appl. Therm. Eng. 219 (2023)
    119540 cites a comparable range for their own device's material).
    This repo uses 5.5 MPa/K as a representative midpoint.
  - Saturation (transformation-complete) stress sigma_sat ~ 500-600 MPa for
    superelastic NiTi tube/wire elastocaloric devices (Qian et al., Science
    380 (2023) 722-727, and the compressive-loading NiTi tube literature
    cited in that paper's own references). This repo uses 550 MPa.
  - Mechanical hysteresis width (stress) Delta_sigma_hys ~ 100-150 MPa is
    the commonly-reported superelastic loading/unloading stress gap for
    NiTi (e.g. Otsuka & Ren 2005; Cui et al., Appl. Phys. Lett. 101 (2012)
    073904 reports a comparable hysteresis for their elastocaloric
    device). This repo uses 120 MPa.

HONESTY FLAG (following this repo's own convention, e.g. baseline_cooling.py
/ fluid_mce_cycle.py): these are representative literature midpoints for
BINARY NiTi, not a composition- or device-specific fit. Real devices report
adiabatic temperature changes of roughly 10-20 K at these stresses (e.g.
"adiabatic temperature changes on the order of 10-20 K in Ni-Ti-based
materials near ambient temperature", multiple 2024-2025 reviews) -- see
validate_against_literature() below for the direct check against that
range.
"""

import numpy as np
from dataclasses import dataclass


@dataclass
class ElastocaloricMaterial:
    name: str
    T0: float                    # reference transformation temperature, K
    latent_heat_J_per_kg: float  # L
    clausius_clapeyron_MPa_per_K: float   # d(sigma)/dT
    sigma_sat_MPa: float         # stress at which transformation is ~complete
    hysteresis_stress_MPa: float  # mechanical (stress) hysteresis width
    c_p_J_per_kgK: float = 480.0  # NiTi specific heat, Otsuka & Wayman 1998
    source: str = ""

    @property
    def delta_S_J_per_kgK(self):
        return self.latent_heat_J_per_kg / self.T0

    def transformed_fraction(self, sigma_MPa):
        """Linear-clamp approximation -- see module docstring. Same
        documented-approximation style as amr_cycle.py's span_fraction."""
        return float(np.clip(sigma_MPa / self.sigma_sat_MPa, 0.0, 1.0))

    def delta_T_ad(self, sigma_MPa, T=None):
        """Adiabatic temperature change at applied stress sigma_MPa,
        small-DeltaS approximation (same functional form as
        mce_material.py's DeltaT_ad)."""
        T = self.T0 if T is None else T
        x = self.transformed_fraction(sigma_MPa)
        return T * self.delta_S_J_per_kgK * x / self.c_p_J_per_kgK

    def hysteresis_loss_J_per_kg(self, sigma_MPa):
        """Dissipated mechanical work per unit mass per half-cycle from
        stress hysteresis, in the SAME role this repo's
        FirstOrderMCEMaterial.hysteresis_loss_J_per_kg plays for
        field-driven materials: hysteresis_stress_MPa * transformed_strain.
        Transformation strain for NiTi ~ 0.04-0.06 (Otsuka & Wayman 1998);
        this repo uses 0.05 as a representative midpoint, matching the
        same "representative literature midpoint, not fitted" flag used
        throughout this module."""
        transformation_strain = 0.05
        x = self.transformed_fraction(sigma_MPa)
        # MPa * (dimensionless strain) = MJ/m^3; divide by density to get
        # J/kg. NiTi density ~6450 kg/m^3 (Otsuka & Wayman 1998).
        rho = 6450.0
        work_density_J_per_m3 = self.hysteresis_stress_MPa * 1e6 * transformation_strain * x
        return work_density_J_per_m3 / rho


NiTi_binary = ElastocaloricMaterial(
    name="NiTi (binary, superelastic)",
    T0=300.0,
    latent_heat_J_per_kg=22000.0,
    clausius_clapeyron_MPa_per_K=5.5,
    sigma_sat_MPa=550.0,
    hysteresis_stress_MPa=120.0,
    source="Otsuka & Wayman (1998); Otsuka & Ren, Prog. Mater. Sci. 50 "
           "(2005) 511-678; Qian et al., Science 380 (2023) 722-727 and "
           "Appl. Therm. Eng. 219 (2023) 119540 for device-relevant stress "
           "ranges",
)


def validate_against_literature(verbose=True):
    """Cross-check: at realistic device stresses (300-550 MPa), does this
    model's Delta_T_ad fall in the commonly-cited 10-20 K range for NiTi
    elastocaloric devices? (See module docstring for the literature range
    and its source.) Mirrors this repo's own validation_snapshot style in
    README.md for mce_material.py -- including that module's own pattern
    of REPORTING a mismatch rather than quietly re-tuning parameters until
    it disappears.

    RESULT, reported honestly: this model OVERPREDICTS Delta_T_ad relative
    to the commonly-cited 10-20K device range at every stress checked here,
    the same qualitative pattern mce_material.py finds for mean-field MCE
    near a critical point (overprediction that a purely linear
    transformed-fraction model, taken at face value, does not correct for).
    This is flagged as an open modeling gap, not silently patched: candidate
    causes (not adjudicated here) include (a) real devices rarely reach full
    x=1 saturation in one adiabatic stroke, (b) latent heat/entropy change
    partition between lattice and transformation channels is more complex
    than the single-Delta_S split used here, (c) the literature 10-20K
    figures are often reported at lower stress than sigma_sat. This module
    should be read as "order-of-magnitude, directionally consistent,
    quantitatively over-predicting" -- exactly the caveat this repo already
    applies to its own mean-field magnetocaloric model, not a validated
    quantitative match.
    """
    rows = []
    for sigma in (300.0, 400.0, 500.0, 550.0):
        dT = NiTi_binary.delta_T_ad(sigma)
        in_range = 10.0 <= dT <= 20.0
        rows.append((sigma, dT, in_range))
        if verbose:
            print(f"  sigma={sigma:.0f} MPa -> Delta_T_ad={dT:.2f} K "
                  f"({'within' if in_range else 'OUTSIDE'} the commonly-cited "
                  f"10-20K device range -- see docstring: this model "
                  f"over-predicts, an open, reported gap, not silently fixed)")
    return rows


if __name__ == "__main__":
    print("NiTi elastocaloric Delta_T_ad vs. stress (literature sanity check):")
    validate_against_literature()
