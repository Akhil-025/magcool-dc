"""
thermal.py
==========
NTU-based packed-bed regenerator effectiveness model for active magnetic
regenerator (AMR) systems.

This model estimates regenerator effectiveness from packed-bed geometry,
heat-transfer correlations, and thermal utilization, allowing regenerator
mass, operating frequency, and fluid flow rate to influence thermal
performance.

Model chain
-----------
Packed-sphere-bed geometry (representative of many experimental AMR
regenerators):

1. Bed geometry

       V_bed = mass_regenerator / (rho_Gd * (1 - phi))

       specific surface area:
           a = 6*(1-phi)/d_p

       total heat-transfer area:
           A_total = a * V_bed

where phi is the packing porosity and d_p is the particle diameter.

2. Convective heat transfer

The interstitial heat-transfer coefficient is computed using the
Wakao–Kaguei packed-bed correlation:

       Nu = 2 + 1.1 Re^0.6 Pr^(1/3)

with

       Re = rho_f u_s d_p / mu_f

and

       h = Nu k_f / d_p

3. Number of transfer units

       NTU = h A_total / (m_dot c_p,f)

4. Utilization factor

The utilization ratio compares the fluid thermal capacity moved during
each half-cycle with the thermal capacity of the regenerator:

       U = (m_dot c_p,f) /
           (2 f m_reg c_p,solid)

5. Regenerator effectiveness

A balanced periodic-flow regenerator approximation is used:

       eps = NTU/(NTU + 2) x (1 - 0.3 U)

with the result clipped to the interval [0, 0.97].

References
----------
Geometry and packed-bed concepts:
    • Tusek et al.
    • Trevizoli & Barbosa (2017)

Heat-transfer correlation:
    • Wakao & Kaguei (1982)

Regenerator theory:
    • Kays & London, Compact Heat Exchangers, 3rd ed. (1984)
    • Engelbrecht (2010)
    • Nielsen et al. (2011)

Limitations
-----------
The utilization correction

       (1 - 0.3 U)

is a phenomenological approximation intended to reproduce the qualitative
reduction in effectiveness observed at higher utilization. The coefficient
0.3 is literature-motivated but has not been calibrated against digitized
experimental effectiveness curves and should therefore be regarded as an
engineering approximation rather than a validated empirical fit.
"""

import numpy as np

RHO_GD = 7900.0              # kg/m^3, gadolinium density (standard literature value)
CP_SOLID_GD = 236.0          # J/(kg K), approx Gd specific heat near room temp
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
