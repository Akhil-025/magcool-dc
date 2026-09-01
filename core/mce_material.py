"""
mce_material.py
================
Mean-field (molecular-field / Brillouin) model of the magnetocaloric effect
for materials exhibiting second-order magnetic phase transitions. This model
is appropriate for elemental gadolinium and similar compounds but is not
intended for first-order giant-MCE materials such as Gd5Si2Ge2, which are
modeled separately using a Landau free-energy formulation.

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
    Peak isothermal DeltaS_M at 5 T ~ -18 J/kg/K (Pecharsky & Gschneidner,
        Phys. Rev. Lett. 78, 4494 (1997) for the giant-MCE Gd5Si2Ge2 family;
        pure Gd baseline ~ -4.2 to -4.8 J/kg/K at 2 T near Tc). Phase 35
        correction: this line previously read "~-9.5 J/kg/K", inconsistent
        with giguere_validation.py's own correctly-cited ~18 J/(kg K) for
        the SAME quantity/paper -- re-checked directly against Ref. [Phys.
        Rev. Lett. 78, 4494 (1997)]'s Fig. 4 (now in this repo's
        Papers/Magnetocaloric effect and materials physics/), which shows
        the 0-5 T Gd5Si2Ge2 curve peaking at very close to 18 J/(kg K) near
        276 K -- 9.5 J/kg/K is not supported by the figure at 5 T (it is
        closer to the paper's own 0-2 T curve's peak instead). This
        constant is documentation/context only -- not read by any function
        in this module or first_order_mce.py -- so the error had no
        effect on any computed result, only on this docstring's own
        accuracy.

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
    curie_shift_K_per_T: float = 0.0
    # Phenomenological field-shift of the Curie point, d(Tc_eff)/d(mu0*H),
    # K per Tesla. Added directly in response to run_curie_shift_check()'s
    # (core/validation.py) documented null result: the plain Brillouin
    # mean-field construction below has NO mechanism for the ordering
    # temperature itself to move with applied field (DeltaS_M(T,H) is built
    # from the same free-energy form at every field, so its peak-location
    # symmetry cannot shift under that construction alone -- see that
    # function's docstring for the full argument). This parameter does NOT
    # emerge from the Brillouin physics; it is a deliberate, explicit patch:
    # the molecular-field constant used when evaluating a NONZERO applied
    # field is derived from an effective Tc_eff(H) = Tc + curie_shift_K_per_T
    # * mu0*H rather than the fixed zero-field Tc (see _lambda_for_field()
    # below). Default 0.0 exactly reproduces this class's original behavior
    # for every material that doesn't opt in -- GADOLINIUM, GD5SI2GE2 and
    # LACAMNO3 below are all UNCHANGED by this addition. A separate,
    # explicitly-shifted instance (GADOLINIUM_FIELD_SHIFTED, see below) is
    # what core/validation.py's calibration/sensitivity check actually
    # exercises, precisely so this addition cannot silently change any
    # existing result that depends on the plain GADOLINIUM instance.
    #
    # Whether Tc_eff(H) is the physically "correct" mechanism behind
    # Dan'kov et al.'s reported ~6 K/T shift is NOT established here --
    # real explanations in the literature (de Oliveira & von Ranke, Phys.
    # Rep. 489 (2010) 89-159, already cited by run_curie_shift_check())
    # invoke short-range-correlation/critical-fluctuation physics well
    # outside a Weiss mean-field treatment. This is a phenomenological knob
    # fit to reproduce the reported SHIFT RATE, not a first-principles
    # derivation of it -- see core/validation.py's calibrate_curie_shift()
    # for the fit and its own honest accounting of what this trades away.

    def __post_init__(self):
        self.N = NA / self.M_molar          # spins per kg
        # Weiss molecular field constant from Tc (mean-field relation)
        self.lam = (3 * kB * self.Tc) / (
            self.N * (self.g ** 2) * (muB ** 2) * self.J * (self.J + 1) * mu0
        )

    def _lambda_for_field(self, H):
        """Molecular-field constant to use when evaluating field H.

        Returns the fixed self.lam unchanged (bit-for-bit) whenever
        curie_shift_K_per_T is 0.0 -- the default, and the value every
        pre-existing material in this module's library uses. Only
        materials that explicitly opt into a nonzero curie_shift_K_per_T
        get a field-dependent, Tc_eff(H)-derived molecular-field constant
        instead. H is expected to be a scalar (A/m) -- every call site in
        this repo passes a scalar applied field, never a per-point array
        (confirmed by inspection of every magnetization()/entropy_magnetic()
        call site in core/), so no array-broadcast handling is implemented.
        """
        if self.curie_shift_K_per_T == 0.0:
            return self.lam
        Tc_eff = max(self.Tc + self.curie_shift_K_per_T * (mu0 * float(H)), 1.0)
        return (3 * kB * Tc_eff) / (
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

    @staticmethod
    def _brillouin_deriv(x, J):
        """dB_J/dx, needed for Newton's method in magnetization(). Analytic
        derivative of the Brillouin function: B_J(x) = a*coth(ax) - b*coth(bx)
        => B_J'(x) = -a^2/sinh^2(ax) + b^2/sinh^2(bx), with the x->0 series
        limit B_J'(0) = (J+1)/(3J) (consistent with _brillouin's small-x
        branch, whose slope this must match)."""
        x = np.asarray(x, dtype=float)
        out = np.zeros_like(x)
        small = np.abs(x) < 1e-8
        a = (2 * J + 1) / (2 * J)
        b = 1 / (2 * J)
        # sinh(z) overflows float64 for z > ~710; beyond z > ~20 the term is
        # already <1e-17 and can be treated as 0 without affecting tol=1e-10
        # convergence in magnetization()'s Newton iteration.
        az, bz = np.clip(a * x[~small], -20, 20), np.clip(b * x[~small], -20, 20)
        sinh_a, sinh_b = np.sinh(az), np.sinh(bz)
        term_a = np.where(np.abs(a * x[~small]) > 20, 0.0, (a ** 2) / sinh_a ** 2)
        term_b = np.where(np.abs(b * x[~small]) > 20, 0.0, (b ** 2) / sinh_b ** 2)
        out[~small] = -term_a + term_b
        out[small] = (J + 1) / (3 * J)
        return out

    def magnetization(self, T, H, tol=1e-10, max_iter=100):
        """Self-consistent solve of M(T,H) via Newton's method on
        F(M) = M - Msat*B_J(x(M)) = 0 (mean-field molecular field theory).

        NOTE (performance): this used to be a damped fixed-point iteration
        (M <- 0.5*M + 0.5*Msat*B_J(x(M))), which is only linearly convergent
        and, near the Curie temperature where the iteration map's effective
        slope approaches 1, took ~300-500 iterations per call to reach
        tol=1e-10 (profiled: ~90k Brillouin evaluations for 200 calls to
        this function, i.e. ~450 iters/call). That single hot loop was the
        reason core/optimize.py's NSGA-III run (2,400 AMR evaluations) took
        minutes instead of seconds. Newton's method is quadratically
        convergent and reaches the same tolerance in ~5-15 iterations
        (verified against the old fixed-point result to <1e-8 relative
        difference across the T range used by validation.py)."""
        T = np.atleast_1d(np.asarray(T, dtype=float))
        T_safe = np.maximum(T, 1e-6)
        Msat = self.N * self.g * muB * self.J
        c = (self.g * muB * self.J * mu0) / (kB * T_safe)
        lam_H = self._lambda_for_field(H)
        M = np.full_like(T_safe, Msat * 0.5)
        for _ in range(max_iter):
            x = c * (H + lam_H * M)
            F = M - Msat * self._brillouin(x, self.J)
            Fp = 1.0 - Msat * self._brillouin_deriv(x, self.J) * c * lam_H
            Fp_safe = np.where(np.abs(Fp) > 1e-12, Fp, 1e-12)
            M_new = np.clip(M - F / Fp_safe, 0.0, Msat)
            if np.max(np.abs(M_new - M)) < tol * Msat:
                M = M_new
                break
            M = M_new
        return M

    @staticmethod
    def _log_sinh_ratio(a, b, x):
        """log(sinh(a*x)/sinh(b*x)) for x >= 0, correct in both limits:

        BUG THIS REPLACES: the previous implementation computed this as
        log(sinh(max(a*x,1e-12))/sinh(max(b*x,1e-12))) -- flooring the two
        sinh *arguments* independently. As x -> 0 (e.g. H=0 above/at Tc,
        where the self-consistent M is exactly 0), both a*x and b*x are
        floored to the *same* 1e-12, so the ratio collapses to sinh(1e-12)/
        sinh(1e-12) = 1 -> log(1) = 0. The correct limit is log(a/b) =
        log(2J+1) (the full spin-degeneracy entropy, N*kB*ln(2J+1)), not 0.
        This silently zeroed out the H=0 entropy reference point for any T
        at or above Tc, which then fed into delta_S_isothermal /
        delta_T_adiabatic as an ~N*kB*ln(2J+1) (~110 J/kg/K for Gd) error --
        large enough to flip C_total negative and produce unphysical
        results (see LITERATURE_REVIEW.md / this fix's changelog entry for
        the T=Tc trace that exposed it).

        For large x this also avoids overflow via log(sinh(z)) ~= z - ln(2).
        """
        x = np.asarray(x, dtype=float)
        out = np.empty_like(x)
        small = x < 1e-4
        xs = x[small]
        # sinh(z) = z + z^3/6 + O(z^5)  =>  ratio -> (a/b)*(1 + (a^2-b^2)x^2/6)
        out[small] = np.log(a / b) + (a ** 2 - b ** 2) * xs ** 2 / 6.0
        xl = x[~small]
        az, bz = a * xl, b * xl
        log_sinh_a = np.where(az > 20, az - np.log(2.0), np.log(np.sinh(az)))
        log_sinh_b = np.where(bz > 20, bz - np.log(2.0), np.log(np.sinh(bz)))
        out[~small] = log_sinh_a - log_sinh_b
        return out

    def entropy_magnetic(self, T, H):
        """Magnetic entropy per kg, J/(kg K), from Brillouin free energy."""
        T = np.atleast_1d(np.asarray(T, dtype=float))
        M = self.magnetization(T, H)
        lam_H = self._lambda_for_field(H)
        x = (self.g * muB * self.J * mu0 * (H + lam_H * M)) / (kB * np.maximum(T, 1e-6))
        a = (2 * self.J + 1) / (2 * self.J)
        b = 1 / (2 * self.J)
        term = self._log_sinh_ratio(a, b, x)
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
        et al., Phys. Rev. B 57, 3478 (1998), experimental Gd heat capacity).

        MODEL LIMITATION: at H=0, this mean-field/Brillouin construction has
        C_mag(T,0) fall to exactly zero for T >= Tc (the self-consistent
        magnetization is exactly 0 in the paramagnetic phase at zero field,
        so entropy_magnetic(T,0) is a T-independent constant there) while
        approaching a finite nonzero value as T -> Tc from below. That is a
        genuine finite-jump discontinuity in idealized mean-field theory
        (the textbook Weiss-theory specific-heat jump at a second-order
        transition, not a numerical artifact of the dT finite difference --
        verified by re-checking at dT=0.02K and finding the same jump, only
        narrower). Real materials round this off via short-range
        correlations/critical fluctuations this model does not include (same
        root cause as validation.py's documented near-Tc DeltaT_ad
        overprediction and run_curie_shift_check()'s null Curie-shift
        result).

        CONSEQUENCE for delta_T_adiabatic(): because C_total's denominator
        drops sharply right above Tc while delta_S_isothermal(T,H) is still
        substantial there, delta_T_adiabatic(T,H) rises steeply in a narrow
        band just above Tc, on top of its normal single-peaked approach to
        that same region. core/amr_cycle.py's AMRSystem.cooling_capacity()
        evaluates delta_T_adiabatic at only one point (T_mid = T_cold +
        span/2) per span; if T_mid sweeps through that steep band as span
        grows, the feasibility margin (2*dTad_noload(T_mid) - span) can go
        negative (Qc clipped to 0) and then swing positive again at a
        LARGER span before finally going negative for good -- which is
        unphysical (a real device's achievable cooling capacity cannot
        increase as the demanded span widens). See
        core/validation_system.py's diagnose_qc_feasibility_reopening() for
        a reusable check of this specific failure mode and
        run_tusek_multipoint_curve_validation()'s use of it to attribute a
        large curve-validation miss (Tusek AMR(A), V*=0.95, span=12.23K:
        +787% vs. literature) to this cause rather than leaving it as an
        unexplained number.
        """
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

    def with_Tc(self, new_Tc):
        """Returns a new MagnetocaloricMaterial identical to this one except
        for Tc (and the derived Weiss constant lambda, recomputed by
        __post_init__ for the new Tc). Every other parameter (J, g, M_molar,
        theta_D, n_atoms_per_fu) is shared unchanged.

        Added for Phase 22 item 1 (core/inhomogeneous_broadening.py): a
        polycrystalline/inhomogeneous sample is modeled as an ensemble of
        grains whose LOCAL Curie temperature is distributed around the
        bulk-reported Tc (grain-to-grain composition/strain variation), each
        grain otherwise behaving as the same mean-field material -- so the
        ensemble only ever needs to vary Tc, never J/g/M_molar/theta_D.
        """
        import dataclasses
        return dataclasses.replace(self, Tc=float(new_Tc))


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
# MODEL LIMITATION
#
# The mean-field/Brillouin framework implemented in this module describes
# second-order (continuous) magnetic phase transitions and is therefore
# appropriate for materials such as elemental gadolinium.
#
# Gd5Si2Ge2 undergoes a first-order magnetostructural transition, which is
# responsible for its giant magnetocaloric effect (Pecharsky &
# Gschneidner, Phys. Rev. Lett. 78, 4494 (1997)). The continuous
# mean-field model cannot reproduce the discontinuous entropy and
# magnetization changes associated with that transition and therefore
# substantially underestimates the material's ΔS_M and ΔT_ad.
#
# The GD5SI2GE2 definition below is retained only as a parameter library
# entry. Quantitative simulations of this material should instead use the
# first-order Landau model implemented in first_order_mce.py, which is
# specifically designed for giant-MCE materials.

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

# GADOLINIUM_FIELD_SHIFTED: identical to GADOLINIUM in every parameter
# except curie_shift_K_per_T, which core.validation.calibrate_curie_shift()
# fits to reproduce Dan'kov et al. (1998)'s reported ~6 K/T Curie-point
# field-shift rate over their own stated 2-7.5 T range (see that module's
# docstring). Deliberately kept as a SEPARATE instance rather than mutating
# GADOLINIUM in place -- every other module in this repo (amr_cycle.py,
# cascade.py, plots.py, validation_system.py, ...) imports GADOLINIUM
# directly and its existing calibration/validation numbers (run_validation(),
# the whole system-level benchmark suite) all assume the ORIGINAL,
# curie_shift_K_per_T=0.0 behavior. curie_shift_K_per_T=0.0 here as a
# placeholder value only; core.validation.calibrate_curie_shift() overwrites
# it via dataclasses.replace() when the fit is actually run, and
# core/validation.py's own module-level constant
# GADOLINIUM_FIELD_SHIFTED_RATE_K_PER_T records whatever rate the last fit
# converged to. See core/validation.py's run_curie_shift_check_v2() for the
# resulting peak-shift comparison AND the honest accounting of what this
# phenomenological patch trades away elsewhere (it is a Tc_eff(H) knob, not
# new physics -- see the curie_shift_K_per_T field docstring above).
GADOLINIUM_FIELD_SHIFTED = MagnetocaloricMaterial(
    name="Gd (polycrystalline, phenomenological Curie-shift patch)",
    Tc=GADOLINIUM.Tc, J=GADOLINIUM.J, g=GADOLINIUM.g, M_molar=GADOLINIUM.M_molar,
    theta_D=GADOLINIUM.theta_D, n_atoms_per_fu=GADOLINIUM.n_atoms_per_fu,
    source=GADOLINIUM.source + " -- with an added phenomenological "
           "curie_shift_K_per_T fit to Dan'kov et al. (1998)'s reported Curie-"
           "point field-shift rate (see core/validation.py's "
           "calibrate_curie_shift()); NOT part of the original calibration "
           "and NOT used by any other module in this repo.",
    curie_shift_K_per_T=0.0,  # overwritten by calibrate_curie_shift()
)