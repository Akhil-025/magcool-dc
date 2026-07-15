"""
amr_cycle.py
============
0-D (lumped) Active Magnetic Regenerator (AMR) cycle model.

The AMR cycle (Barclay, 1982; standard reference cycle for room-temperature
magnetic refrigeration) has four processes analogous to a regenerative
Brayton cycle, but with a solid magnetocaloric regenerator bed instead of a
gas:

    1. Adiabatic magnetization   : bed heats by DeltaT_ad (H: 0 -> H_max)
    2. Cold-to-hot fluid flow    : fluid rejects heat to hot reservoir (Qh)
    3. Adiabatic demagnetization : bed cools by DeltaT_ad (H: H_max -> 0)
    4. Hot-to-cold fluid flow    : fluid absorbs heat from cold reservoir (Qc)

This module wraps an MagnetocaloricMaterial (see mce_material.py) with a
regenerator effectiveness / utilization treatment (Engelbrecht & Bahl, 2010;
Kitanovski et al. 2015, Ch. 4) to give cooling capacity, COP and required
magnetic work per cycle as functions of:
    - temperature span (Th - Tc)
    - peak field mu0*H
    - operating frequency f
    - fluid utilization factor U (heat capacity ratio, fluid/regenerator)
    - regenerator effectiveness eps (NTU-based, from thermal.py)

This is a first-order performance model (matches the level of the AMR
"characteristic curve" approach used in Tusek et al., Int. J. Refrig. 33
(2010) and Nielsen et al., Int. J. Refrig. 34 (2011) 603-616) — good enough
for system-level COP comparison against vapor-compression / liquid cooling,
NOT a replacement for a full 2-D/3-D COMSOL regenerator-bed solve (see
COMSOL_setup.md in the roadmap for that follow-on step).
"""

import numpy as np
from dataclasses import dataclass
from core.mce_material import MagnetocaloricMaterial


@dataclass
class AMRCycleResult:
    T_span: float
    Qc: float          # W, cooling capacity
    Qh: float           # W, heat rejected
    W_mag: float         # W, net magnetic (thermodynamic-cycle) work input
    W_parasitic: float    # W, pump + motor-drive overhead (see note below)
    COP: float            # ideal magnetic-cycle-only COP (Qc / W_mag)
    COP_electrical: float  # device-level electrical COP (Qc / (W_mag + W_parasitic))
                             # -- this is the number comparable to published
                             # "COPe" / device COP figures, and to the
                             # vapor-compression/liquid-cooling baselines in
                             # baseline_cooling.py, which are also electrical.
    exergy_eff: float    # second-law efficiency vs Carnot (magnetic-cycle-only)


class AMRSystem:
    def __init__(self, material: MagnetocaloricMaterial, mu0H_max: float,
                 mass_regenerator: float, frequency: float,
                 fluid_cp: float = 4186.0, fluid_mdot: float = 0.05,
                 regenerator_effectiveness: float = 0.85,
                 parasitic_fraction: float = 0.15):
        """
        material               : MagnetocaloricMaterial instance
        mu0H_max                : peak applied field, Tesla
        mass_regenerator        : kg of magnetocaloric material in the bed
        frequency                : AMR cycle frequency, Hz
        fluid_cp                 : heat transfer fluid specific heat, J/(kg K)
        fluid_mdot                : fluid mass flow rate, kg/s
        regenerator_effectiveness : NTU-based regenerator effectiveness (0-1),
                                     from thermal.py NTU correlation
        parasitic_fraction        : pump + magnet-motor-drive electrical
                                     overhead, as a fraction of Qc, ADDED ON
                                     TOP of the ideal magnetic-cycle work to
                                     get device-level electrical COP. The
                                     default 0.15 is the Phase 2 calibrated
                                     value against two comparably-sized lab
                                     devices (DTU rotary Gd: 0.171, Tusek
                                     single-bed Gd: 0.118 - see
                                     core/validation_system.py). The large
                                     Astronautics naval-cooler prototype
                                     implied 0.453, which Jacobs et al. (2014)
                                     attribute explicitly to "electrical
                                     components with mediocre efficiency" at
                                     that scale/vintage - treat 0.15 as an
                                     optimistic lab-scale figure, not a
                                     production-hardware guarantee, and widen
                                     it in any economics sensitivity study.
        """
        self.mat = material
        self.mu0H_max = mu0H_max
        self.m_reg = mass_regenerator
        self.f = frequency
        self.cp_f = fluid_cp
        self.mdot_f = fluid_mdot
        self.eps = regenerator_effectiveness
        self.parasitic_fraction = parasitic_fraction

    def cooling_capacity(self, T_cold, T_span):
        """Cooling capacity Qc (W) at a given no-load DeltaT_ad and imposed
        span, using the regenerator-effectiveness degradation model:
            Qc = eps * mdot*cp * (DeltaT_ad_local - T_span/2) ... averaged
        which reduces to the standard 'characteristic curve' shape: Qc is
        maximum at zero span and falls roughly linearly to zero at the
        no-load span (Nielsen et al. 2011)."""
        T_hot = T_cold + T_span
        T_mid = 0.5 * (T_cold + T_hot)
        H = self.mu0H_max / (4 * np.pi * 1e-7)
        dTad_noload = float(self.mat.delta_T_adiabatic(np.array([T_mid]), H)[0])
        if dTad_noload <= 0:
            return 0.0, dTad_noload
        span_fraction = max(0.0, 1.0 - T_span / (2 * dTad_noload))
        Qc = self.eps * self.mdot_f * self.cp_f * dTad_noload * span_fraction
        return max(Qc, 0.0), dTad_noload

    def magnetic_work(self, T_cold, T_span, Qc):
        """Net magnetic work input per unit time (W). Approximated from the
        entropy generated by finite-effectiveness regeneration plus the
        ideal (Carnot-referenced) work for the delivered Qc, following the
        second-law decomposition in Kitanovski et al. (2015), Ch. 6:
            W = Qc * (Th/Tc - 1) / eta_2nd_law
        where eta_2nd_law captures AMR irreversibilities (regenerator
        mismatch, viscous dissipation, demagnetization losses) and is taken
        as a literature-informed 0.35-0.55 for well-designed lab-scale AMRs
        (Tusek et al. 2010; Eriksen et al. 2015, Int. J. Refrig. 58)."""
        T_hot = T_cold + T_span
        carnot_work = Qc * (T_hot / T_cold - 1.0) if T_cold > 0 else np.inf
        eta_2nd_law = 0.35 + 0.20 * self.eps  # 0.35 at eps=0 .. 0.52 at eps=0.85
        W = carnot_work / max(eta_2nd_law, 1e-3)
        return W, eta_2nd_law

    def run(self, T_cold, T_span) -> AMRCycleResult:
        Qc, dTad = self.cooling_capacity(T_cold, T_span)
        W, eta2 = self.magnetic_work(T_cold, T_span, Qc)
        W_parasitic = self.parasitic_fraction * Qc
        Qh = Qc + W
        COP = Qc / W if W > 0 else 0.0
        COP_electrical = Qc / (W + W_parasitic) if (W + W_parasitic) > 0 else 0.0
        T_hot = T_cold + T_span
        COP_carnot = T_cold / (T_hot - T_cold) if T_hot > T_cold else np.inf
        exergy_eff = COP / COP_carnot if np.isfinite(COP_carnot) and COP_carnot > 0 else 0.0
        return AMRCycleResult(T_span=T_span, Qc=Qc, Qh=Qh, W_mag=W,
                               W_parasitic=W_parasitic, COP=COP,
                               COP_electrical=COP_electrical, exergy_eff=exergy_eff)

    def characteristic_curve(self, T_cold, spans):
        return [self.run(T_cold, s) for s in spans]
