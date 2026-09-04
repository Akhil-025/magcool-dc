"""
Development/validation script for an Oguchi (Bethe-Peierls) two-spin
cluster correction to the mean-field Brillouin model already in
core/mce_material.py.

Ising-type coupling (Sz-Sz only) is used deliberately, matching the level
of approximation the EXISTING mean-field model already operates at (its
own Brillouin function is the exact self-consistent solution of a
single-site Sz-only Weiss-field Hamiltonian) -- this is a genuine
short-range-correlation correction to that SAME starting point, not a
switch to a different (e.g. full Heisenberg) model.

Cavity/Bethe-Peierls derivation (site with z nearest neighbors):
  H_pair = -J_ex * m1*m2 - (h + h_cav) * (m1 + m2)
  h_cav = (z-1) * J_ex * <m>          [self-consistent cavity field]
  Z_pair(T) = sum_{m1,m2} exp[ (J_ex*m1*m2 + (h+h_cav)*(m1+m2)) / kT ]
  <m> = (1/2) * sum_{m1,m2} (m1+m2) * exp[...] / Z_pair     [self-consistency]

Free energy per site (cavity method, avoids double-counting the z-1
cavity bonds used to build h_cav):
  f/N = -kT * [ (z/2) ln(Z_pair) - (z-1) ln(Z_cav) ]
  Z_cav = sum_m exp[ (h+h_cav)*m / kT ]        <- SAME local field as the pair sees

VALIDATION PLAN (in order -- do not trust this on Gd until step 1 passes):
  1. z -> infinity, J_ex -> 0 holding z*J_ex = 3*kB*Tc_MF/(m_max*(m_max+1))
     fixed (the standard mean-field normalization) MUST recover the
     standard mean-field Brillouin Tc and M(T) EXACTLY -- if this fails,
     there is a derivation/sign/prefactor bug, full stop.
  2. For finite z, Oguchi Tc must come out LOWER than mean-field Tc
     (a structural, model-independent property of every correctly-built
     cluster correction -- short-range correlations always suppress
     ordering relative to mean field, never enhance it).
  3. Specific heat near Tc must be a rounded, finite peak, not a
     discontinuous jump.
Only after all three pass does this get applied to real Gd parameters.

======================================================================
STATUS (Paper-Mining Pass, cross-session): built and internally
validated (passes all 3 checks above) but ULTIMATELY NOT ADOPTED --
when checked against real Gd data it does not resolve the target
near-Tc discrepancy any better than the mean-field model already in
core/mce_material.py, and the session that built this went on to find a
DIRECT, DATA-LEVEL reason no cluster/scaling correction of this kind can
work here: fitting Dan'kov et al.'s own three published DeltaT_ad values
(1, 2, 5T) for the field-scaling exponent n implies n=0.90-1.01 (best
fit ~0.946) -- nowhere near mean-field's n=2/3 (0.667) OR 3D Heisenberg
critical scaling's n=0.637. Both theoretical exponents are further from
the data than each other, meaning Gd at 1-5T is simply not in the
asymptotic critical-scaling regime this class of correction (Oguchi
cluster OR Franco-style universal-curve exponent scaling) is built for.
Kept here as a reference implementation -- genuinely useful for anyone
revisiting issues #1/#2/#3 in the future, so the same two dead ends
aren't re-walked from scratch -- NOT imported or called by anything in
core/ or main.py, and not covered by tests/. See LIMITATIONS.md Item 1.1c.
======================================================================
"""
import numpy as np
from scipy.optimize import brentq

kB = 1.380649e-23


def m_values(J):
    n = int(round(2*J+1))
    return np.linspace(-J, J, n)


def solve_pair(J_ex, h, z, J, T, m_guess=0.0, max_iter=200, tol=1e-12):
    """Self-consistently solve <m> for the Oguchi pair cluster at temperature T
    (K), external reduced field h (energy units, i.e. already g*muB*mu0*H),
    exchange J_ex (energy units), coordination z, spin quantum number J."""
    ms = m_values(J)
    m = m_guess
    kT = kB*T
    for _ in range(max_iter):
        h_cav = (z-1)*J_ex*m
        htot = h + h_cav
        M1, M2 = np.meshgrid(ms, ms, indexing='ij')
        E = -(J_ex*M1*M2 + htot*(M1+M2))
        w = np.exp(-(E - E.min())/kT)
        Z = w.sum()
        m_new = 0.5*np.sum((M1+M2)*w)/Z
        if abs(m_new - m) < tol:
            m = m_new
            break
        m = 0.5*m + 0.5*m_new
    return m, Z, htot


def free_energy_per_site(J_ex, h, z, J, T, m):
    """CORRECTED (2nd derivation) via the Bogoliubov/variational route:
    partition the lattice into N/2 non-overlapping pairs, each pair's own
    bond treated exactly, every OTHER (z-1 per site) bond treated via mean
    field m; F <= F_0 + <H-H_0>_0 gives, after the double-counting
    correction is worked through exactly:
        f(T,m) = -(kT/2)*ln(Z_pair) + (z-1)/2 * J_ex * m^2
    VERIFIED (not assumed): d f/dm = (z-1)*J_ex*(m - <m>_pair(m)), so this
    f is stationary in m EXACTLY where m equals the pair-average <m>_pair(m)
    -- i.e. exactly the self-consistency equation solve_m_robust() already
    solves."""
    kT = kB*T
    ms = m_values(J)
    htot = h + (z-1)*J_ex*m
    M1, M2 = np.meshgrid(ms, ms, indexing='ij')
    E = -(J_ex*M1*M2 + htot*(M1+M2))
    Emin = E.min()
    logZpair = np.log(np.exp(-(E-Emin)/kT).sum()) - Emin/kT
    f = -(kT/2.0)*logZpair + 0.5*(z-1)*J_ex*m**2
    return f


def find_Tc_oguchi(J_ex, z, J, T_lo, T_hi, h=0.0):
    def m_at_T(T):
        m, Z, ht = solve_pair(J_ex, h, z, J, T, m_guess=0.3*J)
        return m
    lo, hi = T_lo, T_hi
    for _ in range(60):
        mid = 0.5*(lo+hi)
        if m_at_T(mid) > 1e-4*J:
            lo = mid
        else:
            hi = mid
    return 0.5*(lo+hi)


def _m_new_from_m(m, J_ex, h, z, J, T):
    ms = m_values(J)
    kT = kB*T
    h_cav = (z-1)*J_ex*m
    htot = h + h_cav
    M1, M2 = np.meshgrid(ms, ms, indexing='ij')
    E = -(J_ex*M1*M2 + htot*(M1+M2))
    w = np.exp(-(E - E.min())/kT)
    Z = w.sum()
    return 0.5*np.sum((M1+M2)*w)/Z


def find_Tc_oguchi_robust(J_ex, z, J, T_lo, T_hi, h=0.0, eps=1e-6):
    """Robust Tc via linear-stability (marginal susceptibility) criterion."""
    def slope_at_zero(T):
        m_plus = _m_new_from_m(eps, J_ex, h, z, J, T)
        return m_plus/eps
    def f(T):
        return slope_at_zero(T) - 1.0
    return brentq(f, T_lo, T_hi, xtol=1e-8)


def solve_m_robust(J_ex, h, z, J, T, m_prev_hint=None):
    """Robust self-consistent <m> solver via root-finding on the residual
    m_new(m) - m over the physical range, NOT naive fixed-point iteration."""
    def resid(m):
        return _m_new_from_m(m, J_ex, h, z, J, T) - m
    grid = np.linspace(J, 1e-6, 400)
    r_prev = resid(grid[0])
    for i in range(1, len(grid)):
        r = resid(grid[i])
        if np.sign(r) != np.sign(r_prev) and r_prev != 0:
            try:
                root = brentq(resid, grid[i], grid[i-1], xtol=1e-10)
                if root > 1e-4:
                    return root
            except Exception:
                pass
        r_prev = r
    return 0.0


def compute_F_S_C(J_ex, h, z, J, T_array):
    """Compute free energy, entropy (via central finite difference on F),
    and specific heat (via central finite difference on S) across a
    temperature array, using the ROBUST self-consistent m solver above."""
    F = np.zeros_like(T_array)
    m_arr = np.zeros_like(T_array)
    for i, T in enumerate(T_array):
        m = solve_m_robust(J_ex, h, z, J, T)
        m_arr[i] = m
        F[i] = free_energy_per_site(J_ex, h, z, J, T, m)
    dT = T_array[1]-T_array[0]
    S = -np.gradient(F, dT)
    C = T_array * np.gradient(S, dT)
    return F, S, C, m_arr
