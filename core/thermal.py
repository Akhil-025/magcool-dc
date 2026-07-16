"""
thermal.py
===========
Phase 4: NTU-based packed-bed regenerator effectiveness model, replacing the
fixed `regenerator_effectiveness=0.85` placeholder used through Phases 1-3.
This directly targets the gap optimize.py's Phase 3 run exposed: with a
fixed eps, `mass_regenerator` had no effect on cooling capacity at all, so
the optimizer trivially minimized it. Here, mass (via bed volume and
heat-transfer area) and frequency/flow (via NTU and utilization) both feed
into eps, so bigger/smaller regenerators have a real, physically motivated
consequence.

Model chain (packed sphere bed, the geometry used in most published AMR
prototypes -- Tusek et al. 2011; Trevizoli & Barbosa 2017 review):

  1. Bed geometry from mass_regenerator, particle diameter d_p and packing
     porosity phi (default 0.365, mid-range of experimentally reported
     0.36-0.37 for packed Gd sphere beds -- Trevizoli et al., Appl. Therm.
     Eng. (2016) porous-matrix geometry comparison; Tusek et al. (2011)):
         V_bed = mass_regenerator / (rho_Gd * (1 - phi))
         specific surface area a = 6*(1-phi)/d_p        [m^2 interstitial
                                                            area per m^3 bed]
         A_total = a * V_bed

  2. Interstitial heat transfer coefficient h from the Wakao & Kaguei (1982)
     packed-bed correlation (widely used in the AMR modeling literature,
     e.g. the packed-bed AMR models cited in Trevizoli & Barbosa's review):
         Nu = 2 + 1.1 * Re^0.6 * Pr^(1/3)
         Re = rho_f * u_s * d_p / mu_f      (u_s = superficial velocity)
         h = Nu * k_f / d_p

  3. NTU = h * A_total / (mdot * cp_f)

  4. Utilization factor (fraction of solid thermal mass "swept" by fluid
     per half-cycle blow), standard AMR definition (Engelbrecht 2010;
     Nielsen et al. 2011):
         U = (mdot * cp_f) / (2 * f * mass_regenerator * cp_solid)

  5. Regenerator effectiveness, balanced periodic-flow regenerator
     approximation (Kays & London, "Compact Heat Exchangers", 3rd ed.,
     1984 -- standard handbook formula for a rotary/matrix regenerator with
     roughly equal blow/counter-blow thermal capacity, adjusted for
     utilization following the same qualitative form used in the AMR
     literature, e.g. Engelbrecht's thesis / Nielsen et al. 2011 Fig. 3-4
     eps-vs-NTU-and-U curves):
         eps = NTU / (NTU + 2)  * (1 - 0.3*U)      [U-degradation term is a
                                                       simple, literature-
                                                       motivated but NOT
                                                       independently fit
                                                       correction -- flagged
                                                       below]
     clipped to [0, 0.97].

**Honesty flag**: step 5's U-degradation factor (1 - 0.3*U) is a
qualitatively-motivated placeholder (higher utilization does reduce
effectiveness in the published eps-NTU-U curves), not a coefficient fit to
data the way loss_model.py's terms are. Phase 5/6 should replace it with a
digitized fit to Nielsen et al. (2011) or Trevizoli et al. (2016)'s actual
eps-NTU-U curves.
"""

import numpy as np

RHO_GD = 7900.0          # kg/m^3, gadolinium density (standard literature value)
CP_SOLID_GD = 236.0        # J/(kg K), approx Gd specific heat near room temp
                             # (Dan'kov et al. 1998 report C_p peaking near
                             # 300 J/kg/K at Tc; 236 J/kg/K is representative
                             # of the broader near-room-temperature range)


def water_properties(T_K=300.0):
    """Simplified constant water properties near room temperature (adequate
    for this 0-D estimate; a full model would use IAPWS correlations)."""
    return {"rho": 997.0, "cp": 4186.0, "mu": 8.9e-4, "k": 0.606}


def regenerator_effectiveness(mass_regenerator, frequency, mdot,
                                particle_diameter=0.0005, porosity=0.365,
                                bed_cross_section_area=0.002, T_K=300.0):
    """Returns (eps, NTU, utilization, h, Re) for a packed-sphere-bed AMR
    regenerator. bed_cross_section_area (m^2) sets superficial velocity from
    mdot; default 0.002 m^2 (~ a 5x4 cm bed face) is representative of the
    lab-scale devices in data/amr_experimental_benchmarks.csv."""
    fluid = water_properties(T_K)
    V_bed = mass_regenerator / (RHO_GD * (1 - porosity))
    a_specific = 6 * (1 - porosity) / particle_diameter   # m^2/m^3
    A_total = a_specific * V_bed

    u_s = mdot / (fluid["rho"] * bed_cross_section_area)   # superficial velocity, m/s
    Re = fluid["rho"] * u_s * particle_diameter / fluid["mu"]
    Pr = fluid["mu"] * fluid["cp"] / fluid["k"]
    Nu = 2 + 1.1 * (max(Re, 1e-6) ** 0.6) * (Pr ** (1 / 3))
    h = Nu * fluid["k"] / particle_diameter

    NTU = h * A_total / (mdot * fluid["cp"]) if mdot > 0 else 0.0
    U = (mdot * fluid["cp"]) / (2 * frequency * mass_regenerator * CP_SOLID_GD) \
        if (frequency > 0 and mass_regenerator > 0) else np.inf

    eps_base = NTU / (NTU + 2)
    eps = eps_base * max(0.0, 1 - 0.3 * min(U, 1.0))
    eps = float(np.clip(eps, 0.0, 0.97))
    return {"eps": eps, "NTU": NTU, "U": U, "h_W_m2K": h, "Re": Re, "A_total_m2": A_total}


if __name__ == "__main__":
    print("Regenerator effectiveness sweep vs. mass_regenerator (f=1Hz, mdot=0.08kg/s)")
    for m in [0.5, 1, 2, 5, 10, 15]:
        r = regenerator_effectiveness(m, frequency=1.0, mdot=0.08)
        print(f"  mass={m:5.1f}kg  NTU={r['NTU']:6.2f}  U={r['U']:6.3f}  eps={r['eps']:.3f}")
    print("\nRegenerator effectiveness sweep vs. frequency (mass=2kg, mdot=0.08kg/s)")
    for f in [0.25, 0.5, 1, 2, 4]:
        r = regenerator_effectiveness(2.0, frequency=f, mdot=0.08)
        print(f"  f={f:5.2f}Hz  NTU={r['NTU']:6.2f}  U={r['U']:6.3f}  eps={r['eps']:.3f}")
