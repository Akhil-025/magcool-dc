"""
mce_material.py
================
Mean-field (molecular-field / Brillouin) model of the magnetocaloric effect (MCE)
for room-temperature magnetic refrigeration materials.

Physics
-------
For a localized-moment ferromagnet (Gd and Gd-based alloys behave close to this
limit), the magnetization follows the Brillouin function:

    M(T, H) = N g mu_B J * B_J(x),      x = g mu_B J mu0(H + lambda*M) / (kB T)

where lambda is the Weiss molecular-field constant, fixed by the Curie
temperature: lambda = 3 kB Tc / (N g^2 mu_B^2 J(J+1) mu0).

The magnetic entropy is obtained from the Brillouin free energy:

    S_M(T, H) = N kB [ln(sinh((2J+1)x/2J) / sinh(x/2J)) - x*B_J(x)]

Isothermal entropy change:      DeltaS_M(T, H) = S_M(T, H) - S_M(T, 0)
Adiabatic temperature change:   DeltaT_ad(T, H) = -T/C_lattice(T) * DeltaS_M(T,H)
                                 (small-DeltaS approximation, standard in AMR
                                 literature, e.g. Tishin & Spichkin 2003;
                                 Kitanovski et al., "Magnetocaloric Energy
                                 Conversion", Springer 2015)

Calibration targets (literature, gadolinium, polycrystalline, ~294 K):
    mu0*DeltaH = 1 T  -> DeltaT_ad ~ 3.0-3.3 K   (Pecharsky & Gschneidner, 1999)
    mu0*DeltaH = 2 T  -> DeltaT_ad ~ 6.1-6.6 K
    mu0*DeltaH = 5 T  -> DeltaT_ad ~ 14-15 K
    Peak isothermal DeltaS_M at 5 T ~ -9.5 J/kg/K (Pecharsky & Gschneidner,
        Phys. Rev. Lett. 78, 4494 (1997) for the giant-MCE Gd5Si2Ge2 family;
        pure Gd baseline ~ -4.2 to -4.8 J/kg/K at 2 T near Tc)

These are used only to validate the mean-field parameters below (see
validation.py) — this module does not hard-code the answer, it computes it
from J, g, Tc and a Debye lattice heat capacity.
"""

import numpy as np
from dataclasses import dataclass

kB = 1.380649e-23      # J/K
muB = 9.2740100783e-24  # J/T
NA = 6.02214076e23
mu0 = 4 * np.pi * 1e-7


@dataclass
class MagnetocaloricMaterial:
    name: str
    Tc: float           # Curie temperature, K
    J: float             # total angular momentum quantum number
    g: float              # Lande g-factor
    M_molar: float        # molar mass, kg/mol
    theta_D: float       # Debye temperature, K (lattice heat capacity)
    n_atoms_per_fu: int = 1   # atoms per formula unit contributing lattice modes
    source: str = ""

    def __post_init__(self):
        self.N = NA / self.M_molar          # spins per kg
        # Weiss molecular field constant from Tc (mean-field relation)
        self.lam = (3 * kB * self.Tc) / (
            self.N * (self.g ** 2) * (muB ** 2) * self.J * (self.J + 1) * mu0
        )

    # ---- Brillouin function and its use in a self-consistent M(T,H) solve ----
    @staticmethod
    def _brillouin(x, J):
        x = np.asarray(x, dtype=float)
        out = np.zeros_like(x)
        small = np.abs(x) < 1e-8
        a = (2 * J + 1) / (2 * J)
        b = 1 / (2 * J)
        out[~small] = a / np.tanh(a * x[~small]) - b / np.tanh(b * x[~small])
        out[small] = ((J + 1) / (3 * J)) * x[small]  # series limit near x=0
        return out

    def magnetization(self, T, H, tol=1e-10, max_iter=500):
        """Self-consistent solve of M(T,H) via fixed-point iteration on the
        Brillouin function (mean-field molecular field theory)."""
        T = np.atleast_1d(np.asarray(T, dtype=float))
        Msat = self.N * self.g * muB * self.J
        M = np.full_like(T, Msat * 0.5)
        for _ in range(max_iter):
            x = (self.g * muB * self.J * mu0 * (H + self.lam * M)) / (kB * np.maximum(T, 1e-6))
            M_new = Msat * self._brillouin(x, self.J)
            if np.max(np.abs(M_new - M)) < tol * Msat:
                M = M_new
                break
            M = 0.5 * M + 0.5 * M_new  # damped update for stability near Tc
        return M

    def entropy_magnetic(self, T, H):
        """Magnetic entropy per kg, J/(kg K), from Brillouin free energy."""
        T = np.atleast_1d(np.asarray(T, dtype=float))
        M = self.magnetization(T, H)
        x = (self.g * muB * self.J * mu0 * (H + self.lam * M)) / (kB * np.maximum(T, 1e-6))
        a = (2 * self.J + 1) / (2 * self.J)
        b = 1 / (2 * self.J)
        with np.errstate(divide="ignore", invalid="ignore"):
            term = np.log(np.sinh(np.maximum(a * x, 1e-12)) / np.sinh(np.maximum(b * x, 1e-12)))
        S = self.N * kB * (term - x * self._brillouin(x, self.J))
        return S

    def delta_S_isothermal(self, T, H_final, H_initial=0.0):
        """Isothermal magnetic entropy change, J/(kg K), applying field H_initial->H_final."""
        return self.entropy_magnetic(T, H_final) - self.entropy_magnetic(T, H_initial)

    def lattice_heat_capacity(self, T, n_debye_points=400):
        """Debye lattice specific heat, J/(kg K)."""
        T = np.atleast_1d(np.asarray(T, dtype=float))
        R = 8.314462618  # J/mol/K
        c_molar = np.zeros_like(T)
        for i, Ti in enumerate(T):
            Ti = max(Ti, 1.0)
            xmax = self.theta_D / Ti
            xs = np.linspace(1e-4, xmax, n_debye_points)
            integrand = (xs ** 4 * np.exp(xs)) / (np.expm1(xs) ** 2)
            trapz_fn = getattr(np, "trapezoid", None) or np.trapz
            integral = trapz_fn(integrand, xs)
            c_molar[i] = 9 * self.n_atoms_per_fu * R * (Ti / self.theta_D) ** 3 * integral
        c_kg = c_molar / self.M_molar
        return c_kg

    def magnetic_heat_capacity(self, T, H=0.0, dT=0.5):
        """Magnetic (lambda-anomaly) contribution to heat capacity, J/(kg K),
        via C_mag = T * dS_M/dT at fixed field. This term peaks sharply at Tc
        and is what a pure-Debye-lattice estimate omits; including it is
        required to reproduce measured total C(T) near Tc (see e.g. Dan'kov
        et al., Phys. Rev. B 57, 3478 (1998), experimental Gd heat capacity)."""
        T = np.atleast_1d(np.asarray(T, dtype=float))
        S_plus = self.entropy_magnetic(T + dT, H)
        S_minus = self.entropy_magnetic(T - dT, H)
        return T * (S_plus - S_minus) / (2 * dT)

    def total_heat_capacity(self, T, H=0.0):
        """C_total = C_lattice + C_magnetic (+ small electronic term, ~ a few
        J/kg/K for Gd, neglected here). This is the physically appropriate
        denominator for DeltaT_ad, not the lattice term alone."""
        return self.lattice_heat_capacity(T) + self.magnetic_heat_capacity(T, H)

    def delta_T_adiabatic(self, T, H_final, H_initial=0.0):
        """Adiabatic temperature change, K, using DeltaT ~ -T*DeltaS_M / C_total(T,H_initial)
        (standard small-signal AMR approximation, e.g. Kitanovski et al. 2015,
        'Magnetocaloric Energy Conversion', Ch. 2)."""
        T = np.atleast_1d(np.asarray(T, dtype=float))
        dS = self.delta_S_isothermal(T, H_final, H_initial)
        C = self.total_heat_capacity(T, H_initial)
        return -T * dS / C


# --- Materials library (parameters from published crystallographic /
#     magnetic characterization data) ---
GADOLINIUM = MagnetocaloricMaterial(
    name="Gd (polycrystalline)",
    Tc=294.0,          # K, standard literature value
    J=7.0 / 2.0,        # 4f^7, S=7/2, L=0 -> J=7/2
    g=2.0,
    M_molar=157.25e-3,  # kg/mol
    theta_D=169.0,       # K, Debye temperature of Gd
    n_atoms_per_fu=1,
    source="Pecharsky & Gschneidner, J. Magn. Magn. Mater. 200 (1999) 44-56; "
           "Tishin & Spichkin, 'The Magnetocaloric Effect and its Applications', IOP (2003)",
)

GD5SI2GE2 = MagnetocaloricMaterial(
    name="Gd5Si2Ge2 (giant MCE)",
    Tc=276.0,
    J=7.0 / 2.0,
    g=2.0,
    M_molar=(5 * 157.25 + 2 * 28.085 + 2 * 72.63) * 1e-3,
    theta_D=200.0,
    n_atoms_per_fu=9,
    source="Pecharsky & Gschneidner, Phys. Rev. Lett. 78, 4494 (1997)",
)
# HONESTY FLAG (found in Phase 4, core/cascade.py): the mean-field/Brillouin
# framework used throughout mce_material.py is built for second-order
# (continuous) magnetic transitions, which is valid for pure Gd. Gd5Si2Ge2's
# "giant" MCE comes from a first-order, coupled magnetostructural phase
# transition (Pecharsky & Gschneidner 1997) -- physics a continuous
# mean-field/Brillouin model cannot capture. Running GD5SI2GE2 through
# delta_T_adiabatic() here UNDERPREDICTS its real effect by roughly an order
# of magnitude (model gives ~1 K at 2 T near its own Tc vs. the several-K-
# to-double-digit effects reported experimentally near a first-order
# transition). Treat GD5SI2GE2 in this codebase as a materials-library
# placeholder, not yet a validated giant-MCE model — a proper treatment
# needs a Bean-Rodbell or Landau free-energy model with magnetoelastic
# coupling, deferred to Phase 5 (see ROADMAP.md).

LACAMNO3 = MagnetocaloricMaterial(
    name="La0.7Ca0.3MnO3 (perovskite manganite)",
    Tc=267.0,
    J=2.0,   # effective Mn moment, approximate
    g=2.0,
    M_molar=208.9e-3,
    theta_D=400.0,
    n_atoms_per_fu=5,
    source="Guo et al., Appl. Phys. Lett. 78, 1142 (1997); Phan & Yu, "
           "J. Magn. Magn. Mater. 308 (2007) 325-340 (review)",
)
