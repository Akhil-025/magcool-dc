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
    • Tusek, Kitanovski, Poredos (2013), Int. J. Refrig. 36, 1456-1464
      (geometry-dependent effectiveness and pumping power, Phase 7)
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

Phase 7 addition -- geometry-dependent pumping power (packed-bed AND
parallel-plate)
-----------------------------------------------------------------------
Tusek, Kitanovski, Poredos, "Geometrical optimization of packed-bed and
parallel-plate active magnetic regenerators", Int. J. Refrig. 36 (2013)
1456-1464 (`Papers/Optimization/Geometrical optimization of packed-bed
and parallel-plate active magnetic regenerators.pdf`) reports optimum
packed-bed sphere diameters of 0.07-0.17 mm and optimum parallel-plate
spacings of 0.035-0.075 mm -- both far smaller than this module's
pre-Phase-7 default `particle_diameter=0.0005 m` (0.5 mm) -- because
smaller particles/gaps raise the heat-transfer coefficient (more
recoverable effectiveness) but ALSO raise the viscous pressure drop for a
given mass-flow rate, so the paper's Figs. 3-5 show a genuine trade-off
optimum in specific cooling load and COP vs. geometry.

Sweeping the pre-Phase-7 `regenerator_effectiveness()`'s
`particle_diameter` alone (see `core/geometry_analysis.py`) shows `eps`
rising monotonically as particle_diameter shrinks, with no optimum --
because this module only ever computed the heat-transfer (NTU) side of
the trade-off; the pumping-power side was handled separately in
`loss_model.py` as `W_pump = k_pump * mdot**2`, a device-calibrated term
with NO particle-diameter or plate-geometry dependence at all. So the
pre-Phase-7 model structure could not, even in principle, reproduce the
paper's Fig. 3-5 trade-off shape. `pressure_drop_packed_bed()`,
`pumping_power_packed_bed()`, and the parallel-plate equivalents below
close that gap using the paper's own friction-factor correlations
(Eqs. 5-6) and hydraulic-diameter definition (Eq. 7, itself sourced by
the paper to Kays & London 1984), combined with the standard
Darcy-Weisbach pressure-drop relation dP = f (L/d_h) (rho u_s^2 / 2)
(the paper cites Kays & London 1984 for the friction factors but does
not itself reproduce a dP formula in the text extracted from the PDF;
this is the standard form associated with those correlations, not a
value taken verbatim from the paper). This new, geometry-explicit
hydraulic pumping power is an IDEALIZED (no pump/motor efficiency
losses) supplementary estimate for exploring the geometry trade-off in
`geometry_analysis.py` -- it does NOT replace the calibrated
`loss_model.StateDependentLossModel`, which remains the production
parasitic-power estimate because it is fitted to real device data.
"""

import numpy as np

RHO_GD = 7900.0              # kg/m^3, gadolinium density (standard literature value)
CP_SOLID_GD = 236.0          # J/(kg K), approx Gd specific heat near room temp
K_SOLID_GD = 10.5            # W/(m K), gadolinium thermal conductivity near room temp
                              # (standard literature value; used by
                              # core/regenerator_1d.py for axial solid-bed
                              # conduction -- see that module's docstring).
                             # (Dan'kov et al. 1998 report C_p peaking near
                             # 300 J/kg/K at Tc; 236 J/kg/K is representative
                             # of the broader near-room-temperature range)
GD_SIGMA_E_S_PER_M = 7.6e5   # S/m, gadolinium electrical conductivity near room
                             # temp. Two independent sources agree within ~4%:
                             # periodictable.org's tabulated resistivity of
                             # 1310 nOhm*m (-> sigma=1/rho=7.63e5 S/m), and a
                             # magnetic-refrigeration-device patent's own stated
                             # value of 0.736e6 S/m used for an eddy-loss
                             # calculation on Gd (US Patent, "Magnetic structure
                             # and magnetic air-conditioning and heating device
                             # using same") -- 7.6e5 S/m used here as their
                             # midpoint. See intragranular_eddy_power()'s own
                             # docstring for the eddy-loss formula this feeds.


def water_properties(T_K=300.0):
    """Simplified constant water properties near room temperature (adequate
    for this 0-D estimate; a full model would use IAPWS correlations)."""
    return {"rho": 997.0, "cp": 4186.0, "mu": 8.9e-4, "k": 0.606}


def regenerator_effectiveness(mass_regenerator, frequency, mdot,
                                particle_diameter=0.0005, porosity=0.365,
                                bed_cross_section_area=0.002, T_K=300.0,
                                cp_solid=None):
    """Returns (eps, NTU, utilization, h, Re) for a packed-sphere-bed AMR
    regenerator. bed_cross_section_area (m^2) sets superficial velocity from
    mdot; default 0.002 m^2 (~ a 5x4 cm bed face) is representative of the
    lab-scale devices in data/amr_experimental_benchmarks.csv.

    cp_solid (Phase 21 addition, core/baseline_cooling.py's
    passive_regenerator_augmentation()): optional override for the solid
    regenerator specific heat used in the utilization term U, J/(kg K).
    Default None reproduces the exact pre-Phase-21 behavior (module-level
    CP_SOLID_GD, a fixed Gd-near-room-temperature value) for every existing
    caller. Passing a temperature-averaged *total* (lattice + magnetic-
    anomaly) heat capacity from core/mce_material.py's own
    MagnetocaloricMaterial.total_heat_capacity() lets a passive regenerator's
    Curie-point heat-capacity peak reduce U (raise buffering capacity per
    cycle) relative to a conventional non-magnetic regenerator material at
    the same mass/frequency/flow -- same additive-override discipline as
    Phase 15's pumping_power_override and Phase 18's thermal_diode=None."""
    fluid = water_properties(T_K)
    cp_solid_eff = CP_SOLID_GD if cp_solid is None else cp_solid
    V_bed = mass_regenerator / (RHO_GD * (1 - porosity))
    a_specific = 6 * (1 - porosity) / particle_diameter   # m^2/m^3
    A_total = a_specific * V_bed

    u_s = mdot / (fluid["rho"] * bed_cross_section_area)   # superficial velocity, m/s
    Re = fluid["rho"] * u_s * particle_diameter / fluid["mu"]
    Pr = fluid["mu"] * fluid["cp"] / fluid["k"]
    Nu = 2 + 1.1 * (max(Re, 1e-6) ** 0.6) * (Pr ** (1 / 3))
    h = Nu * fluid["k"] / particle_diameter

    NTU = h * A_total / (mdot * fluid["cp"]) if mdot > 0 else 0.0
    U = (mdot * fluid["cp"]) / (2 * frequency * mass_regenerator * cp_solid_eff) \
        if (frequency > 0 and mass_regenerator > 0) else np.inf

    eps_base = NTU / (NTU + 2)
    eps = eps_base * max(0.0, 1 - 0.3 * min(U, 1.0))
    eps = float(np.clip(eps, 0.0, 0.97))
    return {"eps": eps, "NTU": NTU, "U": U, "h_W_m2K": h, "Re": Re, "A_total_m2": A_total}


def pressure_drop_packed_bed(mdot, particle_diameter=0.0005, porosity=0.365,
                              bed_cross_section_area=0.002, mass_regenerator=2.0,
                              T_K=300.0):
    """Viscous pressure drop (Pa) across a packed-sphere-bed AMR, using the
    friction-factor correlation Eq. (5) of Tusek, Kitanovski, Poredos
    (2013), Int. J. Refrig. 36, 1456-1464:

        f = 23.462 * Re^-0.6716,   10 < Re < 5e5

    and the standard Darcy-Weisbach relation dP = f (L/d_h) (rho u_s^2/2)
    with the hydraulic diameter from the paper's Eq. (7),
    d_h = 4 V_bed eps / A_total (equivalent to the common packed-bed form
    d_h = (2/3) d_p eps/(1-eps) once A_total = 6(1-eps)/d_p * V_bed is
    substituted). Returns a dict with dP, Re, f, u_s, d_h, L for inspection.
    """
    fluid = water_properties(T_K)
    V_bed = mass_regenerator / (RHO_GD * (1 - porosity))
    a_specific = 6 * (1 - porosity) / particle_diameter
    A_total = a_specific * V_bed
    L = V_bed / bed_cross_section_area
    d_h = 4 * V_bed * porosity / A_total  # Eq. (7)

    u_s = mdot / (fluid["rho"] * bed_cross_section_area)
    Re = fluid["rho"] * u_s * particle_diameter / fluid["mu"]
    Re_c = max(Re, 1e-6)
    f = 23.462 * Re_c ** -0.6716  # Eq. (5), valid 10 < Re < 5e5
    dP = f * (L / d_h) * (fluid["rho"] * u_s ** 2 / 2)
    return {"dP_Pa": dP, "Re": Re, "f": f, "u_s_m_s": u_s, "d_h_m": d_h, "L_m": L}


def pumping_power_packed_bed(mdot, particle_diameter=0.0005, porosity=0.365,
                              bed_cross_section_area=0.002, mass_regenerator=2.0,
                              T_K=300.0):
    """Idealized hydraulic pumping power (W) = dP * volumetric flow rate,
    for a packed-sphere-bed AMR. No pump/motor efficiency is applied (see
    module docstring) -- this is a supplementary, geometry-explicit
    estimate for `geometry_analysis.py`, not the production
    `loss_model.StateDependentLossModel` parasitic-power term."""
    fluid = water_properties(T_K)
    info = pressure_drop_packed_bed(mdot, particle_diameter, porosity,
                                     bed_cross_section_area, mass_regenerator, T_K)
    Q_vol = mdot / fluid["rho"]
    P_pump = info["dP_Pa"] * Q_vol
    info["P_pump_W"] = P_pump
    return info


def intragranular_eddy_power(frequency, mu0H, particle_diameter=0.0005,
                              mass_regenerator=2.0, sigma_e=GD_SIGMA_E_S_PER_M):
    """Phase 27: geometry-explicit eddy-current power dissipated WITHIN the
    MCM particles/plates themselves (W), as a function of particle_diameter
    (or, for a parallel-plate bed, plate_thickness passed through this same
    argument -- see honesty flag below).

    THIS IS A DIFFERENT LOSS CHANNEL from loss_model.StateDependentLossModel's
    own k_eddy*f**2*mu0H**2 term. That term is CORE-calibrated from real AMR
    devices' aggregate parasitic power and is explicitly documented (see
    loss_model.py's module docstring, citing Kitanovski et al. 2015 Ch. 6) as
    representing "eddy-current losses in the magnet/regenerator SUPPORT
    STRUCTURE" -- i.e. conducting metal parts of the magnet assembly and
    housing, not the MCM itself. It has no way to represent the physically
    distinct mechanism a user-supplied document asked this repo to add:
    eddy-current self-heating WITHIN the magnetocaloric material's own
    particles/plates as they cut through dB/dt, which DOES depend on
    particle/plate size (smaller, more electrically-segmented MCM pieces
    dissipate less per unit volume) in a way the support-structure term
    cannot capture. This function is that second, additive channel -- see
    `core.loss_model.StateDependentLossModel.parasitic_power()`'s
    `eddy_power_override` parameter for how the two combine (by ADDITION,
    not replacement -- unlike `pumping_power_override`, since eddy-current
    loss in the support structure and eddy-current loss in the MCM itself
    are physically simultaneous, not alternative estimates of the same
    thing).

    Formula: Pe_volume [W/m^3] = (pi**2/6) * sigma_e * L**2 * f**2 * Bmax**2,
    the standard classical eddy-current loss density for a thin conducting
    slab/lamination of thickness L in a sinusoidal field of peak flux
    density Bmax and frequency f (the SAME pi**2/6 coefficient appears,
    with a real Gd worked example, in a magnetic-refrigeration-device
    patent -- "Magnetic structure and magnetic air-conditioning and heating
    device using same" -- and is standard textbook eddy-current-in-
    laminations physics, not something specific to that patent). Total
    power = Pe_volume * V_MCM, where V_MCM = mass_regenerator / RHO_GD (the
    SOLID material volume -- no porosity factor, since eddy dissipation
    only occurs in the solid MCM, not the fluid-filled pores).

    HONESTY FLAG: this formula is derived for a PLATE/lamination geometry.
    Applied here to `particle_diameter` (a packed-SPHERE-bed parameter) as
    well as to plate thickness, as a deliberate order-of-magnitude
    APPROXIMATION -- classical eddy-current theory gives a sphere a
    DIFFERENT numeric prefactor than a thin slab (the two geometries are
    not the same problem), but no independently-citable sphere-specific
    coefficient was located to use instead of guessing one. Treat this
    function's packed-bed output as approximate in magnitude, not as a
    validated sphere-specific result -- the frequency-squared/field-squared/
    length-squared SCALING (which is what actually matters for the
    optimizer's geometry trade-off, per the motivating document) is robust
    across both geometries; the absolute prefactor is not.

    mu0H here is mu0*H (Tesla), matching this repo's convention elsewhere
    (e.g. StateDependentLossModel.parasitic_power()'s own mu0H argument)."""
    V_MCM = mass_regenerator / RHO_GD
    Pe_volume = (np.pi ** 2 / 6.0) * sigma_e * particle_diameter ** 2 * frequency ** 2 * mu0H ** 2
    return Pe_volume * V_MCM


def regenerator_effectiveness_parallel_plate(mass_regenerator, frequency, mdot,
                                              plate_thickness=0.0005, plate_spacing=0.0002,
                                              bed_cross_section_area=0.002, T_K=300.0):
    """Parallel-plate analogue of `regenerator_effectiveness()`, using the
    Nickolay & Martin (2002) laminar-entry Nusselt correlation (Eq. 4 of
    Tusek et al. 2013), read directly off the rasterized page image
    because pdftotext garbled the exponents in this equation:

        Nu = [ (7.541)^3.592 + (1.841 * Gz^(1/3))^3.592 ]^(1/3.592),
        Gz < 1e5

    a Churchill-Usagi-style blend between the fully-developed
    constant-wall-temperature limit (Nu_inf = 7.541) and the Leveque-type
    thermal-entry-length asymptote (1.841 Gz^(1/3)), with Gz = (d_h/L) Re
    Pr. Same balanced-periodic-flow effectiveness approximation used for
    the packed bed. Geometry follows a parallel-plate unit cell of fluid
    gap `plate_spacing` (r) and solid plate `plate_thickness` (d):
    porosity eps = r/(r+d), specific surface area a = 2/(r+d) (two wetted
    faces per unit-cell period), hydraulic diameter d_h = 2*r
    (parallel-plate slot limit of Eq. 7). Returns the same dict shape as
    `regenerator_effectiveness()`."""
    fluid = water_properties(T_K)
    porosity = plate_spacing / (plate_spacing + plate_thickness)
    V_bed = mass_regenerator / (RHO_GD * (1 - porosity))
    a_specific = 2.0 / (plate_spacing + plate_thickness)
    A_total = a_specific * V_bed
    L = V_bed / bed_cross_section_area
    d_h = 2 * plate_spacing

    u_s = mdot / (fluid["rho"] * bed_cross_section_area)
    Re = fluid["rho"] * u_s * d_h / fluid["mu"]
    Pr = fluid["mu"] * fluid["cp"] / fluid["k"]
    Re_c = max(Re, 1e-6)
    Gz = max((d_h / L) * Re_c * Pr, 1e-6)
    n = 3.592
    Nu = (7.541 ** n + (1.841 * Gz ** (1 / 3)) ** n) ** (1 / n)  # Eq. (4)
    h = Nu * fluid["k"] / d_h

    NTU = h * A_total / (mdot * fluid["cp"]) if mdot > 0 else 0.0
    U = (mdot * fluid["cp"]) / (2 * frequency * mass_regenerator * CP_SOLID_GD) \
        if (frequency > 0 and mass_regenerator > 0) else np.inf

    eps_base = NTU / (NTU + 2)
    eps = eps_base * max(0.0, 1 - 0.3 * min(U, 1.0))
    eps = float(np.clip(eps, 0.0, 0.97))
    return {"eps": eps, "NTU": NTU, "U": U, "h_W_m2K": h, "Re": Re,
            "A_total_m2": A_total, "d_h_m": d_h, "L_m": L, "porosity": porosity}


def pumping_power_parallel_plate(mdot, plate_thickness=0.0005, plate_spacing=0.0002,
                                  bed_cross_section_area=0.002, mass_regenerator=2.0,
                                  T_K=300.0):
    """Idealized hydraulic pumping power (W) for a parallel-plate AMR,
    using the laminar friction factor f = 24/Re (Eq. 6 of Tusek et al.
    2013, valid Re < 2300) with the same Darcy-Weisbach relation and
    hydraulic diameter as `regenerator_effectiveness_parallel_plate()`.
    Idealized (no pump/motor efficiency) -- see module docstring."""
    fluid = water_properties(T_K)
    porosity = plate_spacing / (plate_spacing + plate_thickness)
    V_bed = mass_regenerator / (RHO_GD * (1 - porosity))
    L = V_bed / bed_cross_section_area
    d_h = 2 * plate_spacing

    u_s = mdot / (fluid["rho"] * bed_cross_section_area)
    Re = fluid["rho"] * u_s * d_h / fluid["mu"]
    Re_c = max(Re, 1e-6)
    f = 24.0 / Re_c  # Eq. (6), Re < 2300
    dP = f * (L / d_h) * (fluid["rho"] * u_s ** 2 / 2)
    Q_vol = mdot / fluid["rho"]
    P_pump = dP * Q_vol
    return {"dP_Pa": dP, "P_pump_W": P_pump, "Re": Re, "f": f,
            "u_s_m_s": u_s, "d_h_m": d_h, "L_m": L, "porosity": porosity}


# =============================================================================
# Phase 15 addition: Hypereg parallel-hydraulic pressure-drop reduction
# =============================================================================
#
# Source: Klinar, K. et al., "Perspectives and Energy Applications of
# Magnetocaloric, Pyromagnetic, Electrocaloric and Pyroelectric Materials",
# Adv. Energy Mater. 14, 2401739 (2024), Section "Future Heat Transfer and
# Regenerator Principles" (Figs. 18-20). Read directly from the PDF for
# this pass (Papers/Reviews/... Klinar ... .pdf) rather than relying on
# search-snippet-level knowledge, per the Phase 15 plan for this item.
#
# What the paper actually claims (paraphrased, not quoted -- see
# results/hypereg_findings.md for the full findings note): "Hypereg" is a
# newly-patented (Klinar et al., patent [275] in the review), NOT-YET-
# published-elsewhere regenerator concept whose entire mechanism is
# HYDRAULIC, not electromagnetic: instead of fluid flowing in SERIES
# through one long regenerator bed (pressure-drop length L_dp), it flows
# in PARALLEL through several shorter sub-regenerator beds fed by a
# shared propulsion system, so the pressure-drop-relevant length drops to
# L_dp/n for n parallel sub-beds (the paper's own illustrative Figure 19
# example uses n=4, i.e. L2_dp = (1/4) L1_dp). This is explicitly a
# PRESSURE-DROP / PUMPING-POWER effect (this module's domain: Darcy-flow
# pressure drop, pumping_power_packed_bed()), NOT an eddy-current effect
# (core.loss_model.StateDependentLossModel.k_eddy's domain) -- answering
# the Phase 15 plan's research question directly: it is a thermal.py-style
# hydraulic-geometry change, not a loss_model.py electromagnetic-
# calibration change.
#
# HONESTY FLAGS (read before trusting any number from this function):
#   1. The review states this is the AUTHORS' OWN newly-unveiled, patented
#      concept, described as "our initial assessment" -- there is no
#      published experimental or simulated pressure-drop reduction factor
#      in the paper, and no operating AMR prototype using it exists yet.
#      The n=4 sub-regenerator count is ONLY an illustrative example in
#      Figure 19, not a validated or recommended design value.
#   2. This function does not change the regenerator's heat-transfer
#      effectiveness (regenerator_effectiveness()) at all -- splitting one
#      bed into n parallel, individually-shorter sub-beds at the SAME
#      total mass and total cross-sectional flow area leaves NTU
#      (proportional to heat-transfer area / mdot*cp) structurally
#      unchanged in this 0-D model; the benefit modeled here is PURELY
#      the reduced pumping power (and, downstream, the higher frequency
#      that reduced pumping power can afford), matching the paper's own
#      framing that Hypereg targets pressure drop specifically, not heat
#      transfer.
#   3. No literature source in this corpus gives the "cost" of achieving
#      n-way parallelization (e.g. n separate flow headers, n sets of
#      fluidic diodes/valves -- Figure 20's own b1-b3/c1-c3/d1-d3
#      alternatives) in either $ or additional dead volume/parasitic
#      leakage -- this function is a PUMPING-POWER-ONLY estimate, not a
#      full engineering trade study of the n sub-regenerator designs the
#      paper sketches.
def pumping_power_packed_bed_hypereg(mdot, particle_diameter=0.0005, porosity=0.365,
                                       bed_cross_section_area=0.002, mass_regenerator=2.0,
                                       T_K=300.0, n_parallel_subregenerators=4):
    """Hypereg-style parallel-hydraulic pumping power (W): identical to
    `pumping_power_packed_bed()` except the pressure-drop length L is
    divided by `n_parallel_subregenerators` (Klinar et al. 2024 Fig. 19;
    default n=4 matches the paper's own illustrative example -- see the
    honesty flags in the section docstring above, especially flag #1).
    Total regenerator mass, cross-section area, and particle diameter are
    unchanged from the conventional-series case; only the flow path is
    reconfigured, so `regenerator_effectiveness()` (heat-transfer side) is
    unaffected -- callers should pair this with the SAME eps as the
    conventional case, not a re-derived one (see honesty flag #2)."""
    if n_parallel_subregenerators < 1:
        raise ValueError("n_parallel_subregenerators must be >= 1")
    fluid = water_properties(T_K)
    info = pressure_drop_packed_bed(mdot, particle_diameter, porosity,
                                     bed_cross_section_area, mass_regenerator, T_K)
    L_hypereg = info["L_m"] / n_parallel_subregenerators
    dP_hypereg = info["f"] * (L_hypereg / info["d_h_m"]) * (fluid["rho"] * info["u_s_m_s"] ** 2 / 2)
    Q_vol = mdot / fluid["rho"]
    P_pump_hypereg = dP_hypereg * Q_vol
    return {
        "dP_Pa": dP_hypereg, "P_pump_W": P_pump_hypereg,
        "Re": info["Re"], "f": info["f"], "u_s_m_s": info["u_s_m_s"],
        "d_h_m": info["d_h_m"], "L_m": L_hypereg,
        "n_parallel_subregenerators": n_parallel_subregenerators,
        "P_pump_W_conventional_series": pumping_power_packed_bed(
            mdot, particle_diameter, porosity, bed_cross_section_area,
            mass_regenerator, T_K)["P_pump_W"],
    }


if __name__ == "__main__":
    print("Regenerator effectiveness sweep vs. mass_regenerator (f=1Hz, mdot=0.08kg/s)")
    for m in [0.5, 1, 2, 5, 10, 15]:
        r = regenerator_effectiveness(m, frequency=1.0, mdot=0.08)
        print(f"  mass={m:5.1f}kg  NTU={r['NTU']:6.2f}  U={r['U']:6.3f}  eps={r['eps']:.3f}")
    print("\nRegenerator effectiveness sweep vs. frequency (mass=2kg, mdot=0.08kg/s)")
    for f in [0.25, 0.5, 1, 2, 4]:
        r = regenerator_effectiveness(2.0, frequency=f, mdot=0.08)
        print(f"  f={f:5.2f}Hz  NTU={r['NTU']:6.2f}  U={r['U']:6.3f}  eps={r['eps']:.3f}")