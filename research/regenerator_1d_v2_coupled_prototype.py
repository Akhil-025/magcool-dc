"""
Standalone prototype of a properly-coupled 1D AMR regenerator solver,
following Nielsen et al. 2011 (Int. J. Refrig. 34, 603-616), Eqs. 5-8:
  solid:  rho_s cp_s dT_s/dt = d/dx(k_s dT_s/dx) + Q_MCE + Q_HT
  fluid:  rho_f cp_f (dT_f/dt + u dT_f/dx) = -Q_HT      (axial fluid conduction neglected, standard)
  Q_HT = h*A_HT*(T_s - T_f)

Design: at each timestep dt,
  1. Freeze T_s, compute the QUASI-STEADY fluid temperature profile through
     the bed for this dt (exact exponential solution per node given fixed
     T_s -- same trick the original code used, valid since fluid response
     is much faster than solid).
  2. Update T_s via an IMPLICIT (backward-Euler) tridiagonal solve that
     handles axial conduction AND the convective heat source SIMULTANEOUSLY
     in one linear solve, not as two sequential/split operators.
This directly avoids the operator-splitting-order issue diagnosed earlier,
and implicit conduction is unconditionally stable, so grid refinement is
not expected to require re-tuning stability substeps.

======================================================================
STATUS (Paper-Mining Pass, cross-session): built and run, but does NOT
resolve issue #8 (regenerator span undershoot). Root cause NOT found.

Hypotheses tested and ELIMINATED by this prototype and the session that
built it (kept here so they aren't re-tested from scratch):
  1. Operator-splitting order (conduction vs. convection applied
     sequentially, in the wrong order, in core/regenerator_1d.py) --
     this prototype couples them in one implicit solve instead; span is
     still ~0, so splitting order is not the cause.
  2. Axial-conduction magnitude -- separately checked, not the driver.
  3. NTU magnitude -- the pre-existing calibrated flow rate implied
     NTU~487-1949, far outside the realistic NTU~10-30 range reported for
     working AMR devices (Nielsen et al. 2010, 2D AMR parameter study).
     Re-tested this SAME coupled solver at mdot values giving a realistic
     NTU (~0.03-0.07 kg/s) -- span is STILL ~0 across that whole range.
     So the NTU magnitude was a real, separate bug worth knowing about,
     but is not, by itself, the explanation for the span undershoot.
  4. Grid resolution -- not the primary driver either, given (1)-(3)
     above already fail at multiple resolutions.

Not imported or called by anything in core/ or main.py; not covered by
tests/. Kept as a reference implementation only, so future work on issue
#8 starts from "four hypotheses already ruled out" instead of zero.

Independently re-run against this repo's current core/mce_material.py
(Paper-Mining Pass verification pass): the loop as written here is not
even numerically stable over enough cycles -- it diverges to NaN via
Brillouin-function overflow (uncontrolled temperature drift with no
clamping/damping in this quick-prototype harness), reinforcing that this
is exploratory, not a drop-in replacement for anything.
See LIMITATIONS.md Item 1.3.
======================================================================
"""
import numpy as np

RHO_GD = 7900.0

def thomas_solve(a, b, c, d):
    """Solve a tridiagonal system. a=sub-diag, b=diag, c=super-diag, d=RHS.
    a[0] and c[-1] are unused. Standard Thomas algorithm."""
    n = len(b)
    cp = np.zeros(n); dp = np.zeros(n)
    cp[0] = c[0]/b[0]
    dp[0] = d[0]/b[0]
    for i in range(1, n):
        m = b[i] - a[i]*cp[i-1]
        cp[i] = c[i]/m if i < n-1 else 0.0
        dp[i] = (d[i]-a[i]*dp[i-1])/m
    x = np.zeros(n)
    x[-1] = dp[-1]
    for i in range(n-2, -1, -1):
        x[i] = dp[i]-cp[i]*x[i+1]
    return x


def simulate_v2(material, mu0H_max, mass_total, frequency, mdot,
                 n_nodes=20, porosity=0.2857, bed_cross_section_area=3.9e-4,
                 plate_thickness=0.00025, plate_spacing=0.0001,
                 T_K_for_ntu=290.0, max_cycles=2000, tol=1e-5,
                 n_substeps_per_blow=30):
    import core.regenerator_1d as r1d
    import core.thermal as th

    fluid = th.water_properties(T_K_for_ntu)
    cp_f = fluid["cp"]
    k_f = fluid["k"]
    cp_solid = r1d.CP_SOLID_GD
    k_solid = r1d.K_SOLID_GD

    V_bed = mass_total/(RHO_GD*(1-porosity))
    L_bed = V_bed/bed_cross_section_area
    dx = L_bed/n_nodes
    m_node = mass_total/n_nodes

    k_eff_axial = r1d._parallel_plate_effective_axial_conductivity(porosity, k_f, k_solid)

    ntu_result = th.regenerator_effectiveness_parallel_plate(
        mass_total, frequency, mdot, plate_thickness=plate_thickness,
        plate_spacing=plate_spacing, bed_cross_section_area=bed_cross_section_area,
        T_K=T_K_for_ntu)
    NTU_total = ntu_result["NTU"]
    NTU_node = NTU_total/n_nodes
    decay = 1 - np.exp(-NTU_node)

    cycle_period = 1.0/frequency
    tau_blow = cycle_period/2.0
    dt = tau_blow/n_substeps_per_blow
    H_Am = mu0H_max/(4*np.pi*1e-7)

    G = k_eff_axial*bed_cross_section_area/dx
    alpha_dt = dt/(m_node*cp_solid)

    def implicit_conduction_solve(T_s_old, Q_source):
        n = n_nodes
        a = np.zeros(n); b = np.zeros(n); c = np.zeros(n); d = np.zeros(n)
        for i in range(n):
            neighbors = 0
            if i > 0: neighbors += 1
            if i < n-1: neighbors += 1
            b[i] = 1 + alpha_dt*G*neighbors
            if i > 0: a[i] = -alpha_dt*G
            if i < n-1: c[i] = -alpha_dt*G
            d[i] = T_s_old[i] + alpha_dt*Q_source[i]
        return thomas_solve(a, b, c, d)

    T = np.full(n_nodes, T_K_for_ntu)
    T_fluid_cold_end = T_K_for_ntu
    T_fluid_hot_end = T_K_for_ntu
    span_history = []

    for cycle in range(max_cycles):
        T = T + material.delta_T_adiabatic(T, H_Am)

        T_in = T_fluid_cold_end
        for _ in range(n_substeps_per_blow):
            T_f = T_in
            Q = np.zeros(n_nodes)
            for i in range(n_nodes):
                T_f_out = T[i] - (T[i]-T_f)*decay
                Q[i] = mdot*cp_f*(T_f - T_f_out)
                T_f = T_f_out
            T_f_exit = T_f
            T = implicit_conduction_solve(T, Q)
        T_fluid_hot_end = T_f_exit

        T = T - material.delta_T_adiabatic(T, H_Am)

        T_in = T_fluid_hot_end
        for _ in range(n_substeps_per_blow):
            T_f = T_in
            Q = np.zeros(n_nodes)
            for i in reversed(range(n_nodes)):
                T_f_out = T[i] - (T[i]-T_f)*decay
                Q[i] = mdot*cp_f*(T_f - T_f_out)
                T_f = T_f_out
            T_f_exit = T_f
            T = implicit_conduction_solve(T, Q)
        T_fluid_cold_end = T_f_exit

        span_history.append(float(T[-1]-T[0]))
        if cycle > 20:
            recent = span_history[-10:]
            if (max(recent)-min(recent)) < tol:
                return {"span_K": float(np.mean(recent)), "converged": True, "n_cycles": cycle+1}
    return {"span_K": float(np.mean(span_history[-10:])), "converged": False,
            "n_cycles": max_cycles, "last10": span_history[-10:]}
