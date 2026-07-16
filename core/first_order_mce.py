"""
first_order_mce.py
====================
Implements an extended Landau free-energy model for Gd5Si2Ge2, replacing
the mean-field/Brillouin treatment with a model capable of describing the
first-order magnetostructural transition responsible for the giant
magnetocaloric effect (Pecharsky & Gschneidner, Phys. Rev. Lett. 78, 4494
(1997)).

Model
-----
Free energy per magnetic ion, in units of k_B*Tc, as a function of reduced
magnetization m = M/M_sat and reduced temperature tau = T/Tc:

    f(m, tau, h) = (A/2)*(tau-1)*m^2 + (B/4)*m^4 + (C/6)*m^6 - h*m

This is the standard "extended Landau theory" (quadratic-quartic-sextic
expansion) used for first-order magnetocaloric materials when B < 0 (the
quartic term must be negative to produce a discontinuous jump in the
equilibrium m at h=0, tau=1; the positive sextic term C is required for
free-energy stability). This is the same qualitative approach as the
classical Bean & Rodbell (1962) volume-strain-coupled model, but expressed
directly as a phenomenological order-parameter expansion rather than via
self-consistent lattice strain -- a standard simplification used e.g. in
de Oliveira & von Ranke's review (Phys. Rep. 489 (2010) 89-159, Section 4)
of first-order MCE models, which is also the source cited for the
mean-field-theory limitations flagged elsewhere in this codebase.

Reduced field:
    h = g*J*mu_B*(mu0*H) / (kB*Tc)     (same natural scale as the Brillouin
                                          argument used in mce_material.py)

Equilibrium m(tau,h): real root of  A*(tau-1)*m + B*m^3 + C*m^5 = h  that
GLOBALLY MINIMIZES f (not just any stationary point) -- this correctly
selects the equilibrium/reversible branch through the first-order jump,
consistent with treating S as a state function for the ideal AMR cycle
(hysteresis/irreversibility at the transition is a real effect this
simplified treatment does not capture -- flagged below).

Entropy (envelope theorem, since df/dm=0 at equilibrium):
    S(tau,h)/  (N*kB) = -(A/2) * m(tau,h)^2   [+ const, cancels in DeltaS]
    DeltaS_M(tau,h) = -(A/2) * [m(tau,h)^2 - m(tau,0)^2] * N*kB

Calibration: (A, B, C) = (10, -4, 8) were found by grid search to reproduce
the most consistently cited literature peak entropy change for Gd5Si2Ge2,
|DeltaS_M| ~ 18 J/(kg K) at mu0*DeltaH = 5 T near Tc = 276-278 K (Gschneidner
& Pecharsky review value, cross-checked against the original 1997 PRL and
the Journal of Superconductivity and Novel Magnetism (2019) DFT+Monte Carlo
reproduction, which reports 9.97 J/(kg K) at 5 T -- literature itself spans
roughly 10-18.5 J/(kg K) depending on sample preparation/purity, per
Pecharsky & Gschneidner's own note that "optimally prepared" samples show
the largest effect. The calibration here targets the upper (widely-quoted
"~18 J/kgK") end -- treat this as ONE defensible calibration choice, not
the only one, and note the real spread when citing this in the paper.

**Honesty flags**:
  1. delta_T_adiabatic here uses DeltaT ~ -T*DeltaS_M/C_lattice, using ONLY
     the Debye lattice heat capacity (unlike GADOLINIUM's second-order
     treatment in mce_material.py, which also adds a magnetic lambda-anomaly
     C_mag). That correction is appropriate for continuous transitions; for
     a first-order transition the physically correct denominator involves
     the transition's latent heat structure, which this 0-D model does not
     resolve. Treat ΔT_ad from this module as an upper-bound-ish estimate,
     not a validated number -- there is no direct system-level benchmark
     for it in this codebase (see data/amr_experimental_benchmarks.csv,
     which has none for Gd5Si2Ge2).
    2. No independent validation dataset exists in this repository for this
    model (unlike GADOLINIUM's Dan'kov et al. (1998) validation in
    validation.py). It is calibrated to a single literature value
    (peak DeltaS_M) and has not been cross-checked against an independent
    experimental dataset. Future validation should include the direct
    DeltaT_ad measurements of Giguère et al. (Phys. Rev. Lett. 83, 2262
    (1999)) before using this model for quantitative design conclusions.
"""

import numpy as np
from dataclasses import dataclass

kB = 1.380649e-23
muB = 9.2740100783e-24
NA = 6.02214076e23
mu0 = 4 * np.pi * 1e-7


@dataclass
class FirstOrderMCEMaterial:
    name: str
    Tc: float
    J: float
    g: float
    M_molar: float
    theta_D: float
    n_atoms_per_fu: int
    A: float
    B: float
    C: float
    source: str = ""

    def __post_init__(self):
        self.N = NA / self.M_molar

    def _h_reduced(self, mu0H):
        return self.g * self.J * muB * mu0H / (kB * self.Tc)

    def _equilibrium_m(self, tau, h):
        coeffs = [self.C, 0, self.B, 0, self.A * (tau - 1), -h]
        roots = np.roots(coeffs)
        real_roots = roots[np.abs(roots.imag) < 1e-6].real
        real_roots = real_roots[np.abs(real_roots) <= 1.5]
        if len(real_roots) == 0:
            return 0.0

        def f(m):
            return (0.5 * self.A * (tau - 1) * m ** 2 + 0.25 * self.B * m ** 4
                    + (self.C / 6) * m ** 6 - h * m)
        vals = [f(m) for m in real_roots]
        return real_roots[int(np.argmin(vals))]

    def delta_S_isothermal(self, T, H_final, H_initial=0.0):
        T = np.atleast_1d(np.asarray(T, dtype=float))
        out = np.zeros_like(T)
        for i, Ti in enumerate(T):
            tau = Ti / self.Tc
            # H_final/H_initial are passed in A/m (same convention as
            # mce_material.py), so mu0*H converts to Tesla for _h_reduced
            h_f = self._h_reduced(mu0 * H_final)
            h_i = self._h_reduced(mu0 * H_initial)
            m_f = self._equilibrium_m(tau, h_f)
            m_i = self._equilibrium_m(tau, h_i)
            s_f = -0.5 * self.A * m_f ** 2
            s_i = -0.5 * self.A * m_i ** 2
            out[i] = (s_f - s_i) * self.N * kB
        return out

    def lattice_heat_capacity(self, T, n_debye_points=400):
        T = np.atleast_1d(np.asarray(T, dtype=float))
        R = 8.314462618
        c_molar = np.zeros_like(T)
        for i, Ti in enumerate(T):
            Ti = max(Ti, 1.0)
            xmax = self.theta_D / Ti
            xs = np.linspace(1e-4, xmax, n_debye_points)
            integrand = (xs ** 4 * np.exp(xs)) / (np.expm1(xs) ** 2)
            trapz_fn = getattr(np, "trapezoid", None) or np.trapz
            integral = trapz_fn(integrand, xs)
            c_molar[i] = 9 * self.n_atoms_per_fu * R * (Ti / self.theta_D) ** 3 * integral
        return c_molar / self.M_molar

    def delta_T_adiabatic(self, T, H_final, H_initial=0.0):
        T = np.atleast_1d(np.asarray(T, dtype=float))
        dS = self.delta_S_isothermal(T, H_final, H_initial)
        C = self.lattice_heat_capacity(T)
        return -T * dS / C


GD5SI2GE2_FIRST_ORDER = FirstOrderMCEMaterial(
    name="Gd5Si2Ge2 (first-order Landau model)",
    Tc=276.0, J=3.5, g=2.0,
    M_molar=(5 * 157.25 + 2 * 28.085 + 2 * 72.63) * 1e-3,
    theta_D=200.0, n_atoms_per_fu=9,
    A=10.0, B=-4.0, C=8.0,
    source="Landau coefficients calibrated to peak |DeltaS_M|~18 J/(kg K) "
           "at 5T (Pecharsky & Gschneidner 1997; Gschneidner & Pecharsky "
           "review); NOT independently validated against a second dataset "
           "(see module docstring honesty flag #2).",
)

if __name__ == "__main__":
    mu0_ = 4 * np.pi * 1e-7
    print("First-order Landau model calibration check, Gd5Si2Ge2, T=Tc=276K")
    for B_T in [1, 2, 5]:
        H = B_T / mu0_
        dS = GD5SI2GE2_FIRST_ORDER.delta_S_isothermal(np.array([276.0]), H)
        dT = GD5SI2GE2_FIRST_ORDER.delta_T_adiabatic(np.array([276.0]), H)
        print(f"  {B_T}T: dS={dS[0]:.2f} J/(kg K)   dTad={dT[0]:.2f} K")
    print("\nTarget: dS ~ -18 J/(kg K) at 5T (Pecharsky & Gschneidner 1997 review value)")
