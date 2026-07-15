"""
baseline_cooling.py
====================
Reference COP models for the two incumbent data-center cooling technologies
that magnetic cooling is being benchmarked against:

1. Vapor-compression CRAC/CRAH (air-cooled, room/row-based)
2. Direct liquid cooling (cold-plate, facility-water-cooled)

Both are modeled as fractions of the reverse-Carnot COP, with the fraction
(Lorenz/second-law efficiency) taken from published data-center cooling
efficiency studies rather than assumed:

    - Vapor-compression DX/chiller plants in data centers: second-law
      efficiency ~ 0.35-0.45 of Carnot for packaged CRAC units, up to
      ~0.5-0.55 for well-optimized chilled-water plants.
      (ASHRAE Datacom Series; Shah, Bash & Patel, "Cooling and Power
      Considerations for Chips", ASME (2004); Ebrahimi, Jones & Fleischer,
      Renew. Sustain. Energy Rev. 31 (2014) 622-638 - DC cooling review)

    - Direct liquid cooling (cold plate, facility water loop): removes the
      need for a low-temperature chiller stage entirely for a large fraction
      of the load (can use free/economizer cooling at facility-water
      temperatures ~18-32 C per ASHRAE TC9.9 W-class envelopes), so its
      *effective* COP is reported as the compressor-side COP only when
      mechanical cooling is still needed, otherwise pump-only COP is very
      high (>20). We model both the "mechanical assist" and "economizer"
      regimes.
      (ASHRAE TC9.9, "Liquid Cooling Guidelines for Datacom Equipment
      Centers", 2nd ed. 2021; Ellsworth, Iyengar et al., IEEE ITherm
      proceedings, various years)
"""

import numpy as np
from dataclasses import dataclass


@dataclass
class CoolingResult:
    technology: str
    Tc: float
    Th: float
    COP: float
    COP_carnot: float
    second_law_eff: float


def carnot_cop(Tc, Th):
    return Tc / (Th - Tc) if Th > Tc else np.inf


def vapor_compression_cop(Tc, Th, eta_2nd_law=0.42, source_note="packaged DX CRAC"):
    """eta_2nd_law default 0.42 is representative of packaged DX CRAC units
    per Ebrahimi et al. (2014) data-center cooling review; use 0.50-0.55 for
    optimized chilled-water plants."""
    cc = carnot_cop(Tc, Th)
    return CoolingResult("Vapor-compression (%s)" % source_note, Tc, Th,
                          eta_2nd_law * cc, cc, eta_2nd_law)


def liquid_cooling_cop(Tc, Th, economizer_hours_fraction=0.6,
                        mechanical_eta_2nd_law=0.42, pump_equivalent_cop=25.0):
    """Direct/indirect liquid cooling: for `economizer_hours_fraction` of
    operating hours the facility water loop can reject heat directly (dry
    cooler / cooling tower) without a compressor -- effective COP is very
    high (pump + fan power only, modeled here at pump_equivalent_cop, a
    literature-informed placeholder per ASHRAE TC9.9 W-class case studies).
    For the remaining hours mechanical cooling (chiller) is engaged at the
    same second-law efficiency as vapor compression. Result is an
    hours-weighted annual average COP."""
    cc = carnot_cop(Tc, Th)
    cop_mech = mechanical_eta_2nd_law * cc
    cop_avg = (economizer_hours_fraction * pump_equivalent_cop +
               (1 - economizer_hours_fraction) * cop_mech)
    eff_2nd_law_blended = cop_avg / cc if np.isfinite(cc) and cc > 0 else 0.0
    return CoolingResult("Liquid cooling (blended, %.0f%% economizer hrs)"
                          % (economizer_hours_fraction * 100),
                          Tc, Th, cop_avg, cc, eff_2nd_law_blended)


def compare_all(amr_result, Tc, Th, **kwargs):
    """Convenience: returns dict of technology -> CoolingResult/AMRCycleResult
    for the same (Tc, Th) operating point."""
    vcc = vapor_compression_cop(Tc, Th)
    liq = liquid_cooling_cop(Tc, Th)
    return {
        "Magnetic (AMR)": amr_result,
        "Vapor-compression": vcc,
        "Liquid cooling": liq,
    }
