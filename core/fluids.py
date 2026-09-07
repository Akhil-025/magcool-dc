"""
fluids.py
=========
Heat-transfer-fluid (HTF) property library for AMR system models.

WHY THIS EXISTS
---------------
Every other module in this repo (thermal.py, regenerator_1d.py,
fluid_mce_cycle.py, geometry_analysis.py, amr_cycle.py) previously called
a single hardcoded function, `thermal.water_properties()`, which always
returned pure-water properties. Real magnetocaloric-regenerator hardware
essentially never runs on pure water, because Gd (and most other MCM
candidates) corrodes in prolonged contact with it. Every prototype in
this project's own Papers/ corpus that reports its actual heat-transfer
fluid uses water plus a corrosion-inhibiting glycol fraction instead:

  * Nielsen et al.'s DTU rotary AMR ("Active magnetic regenerator
    refrigeration with rotary multi-bed technology") and "Improving
    magnetic cooling efficiency and pulldown by varying flow profiles":
    90 vol% deionized water + 10 vol% mono-ethylene-glycol.
  * "Development of a novel rotary magnetic refrigerator" and the
    high-frequency rotating AMR paper: 80 vol% distilled/deionised water +
    20 vol% automotive antifreeze (ethylene glycol).
  * "The performance of a large-scale rotary magnetic refrigerator" and
    the nine-layer-regenerator paper: water "with addition of an
    anti-corrosion agent" (unspecified ratio, same family).

This module exposes those literature-grounded options (plus a
low-toxicity propylene-glycol alternative and pure ethanol, included for
comparison only -- see each entry's description) as named, swappable
fluids, instead of silently hardcoding one choice.

PROPERTY VALUES
---------------
Constant, evaluated near 300 K / 27 degC -- the SAME "adequate for a 0-D
estimate" idealization the original `water_properties()` used (T_K is
accepted for interface compatibility and future temperature-dependent
correlations, but does not currently vary the returned values). Sources:
standard water/glycol engineering property compilations (ASHRAE
Fundamentals-style glycol tables; Sun & Teja, J. Chem. Eng. Data 2003,
48(1) 198-202 for ethylene-glycol/water density-viscosity-conductivity;
Bohne, Fischer & Obermeier, Ber. Bunsenges. Phys. Chem. 1984, 88(8)
739-742 for EG-water transport properties), cross-checked against the
specific volume fractions this project's own prototype papers report.
These are representative literature values, not a per-mixture
first-principles derivation -- treat them at the same precision as the
original single-fluid water_properties() (order-of-percent, not exact).

SEE ALSO: core/fluid_selection_optimization.py, which sweeps these
fluids through AMRSystem at this project's standard baseline operating
point to compare COP/Qc/pumping-power trade-offs and justifies
DEFAULT_FLUID below; results/fluid_selection_optimization.txt has the
numeric comparison.
"""

FLUID_LIBRARY = {
    "water": {
        "rho": 997.0, "cp": 4186.0, "mu": 8.9e-4, "k": 0.606,
        "freeze_C": 0.0,
        "description": (
            "Pure water (original idealized baseline). Not used by any "
            "real AMR prototype in this project's corpus -- real hardware "
            "always adds a corrosion inhibitor, see other entries."
        ),
    },
    "water_eg10": {
        "rho": 1013.0, "cp": 4020.0, "mu": 1.05e-3, "k": 0.58,
        "freeze_C": -4.0,
        "description": (
            "90 vol% deionized water / 10 vol% mono-ethylene-glycol -- "
            "DTU rotary AMR heat-transfer fluid (Nielsen et al.; "
            "'Improving magnetic cooling efficiency...' paper), added to "
            "prevent Gd corrosion."
        ),
    },
    "water_eg20": {
        "rho": 1028.0, "cp": 3860.0, "mu": 1.4e-3, "k": 0.53,
        "freeze_C": -8.0,
        "description": (
            "80 vol% distilled water / 20 vol% automotive antifreeze "
            "(ethylene glycol) -- the most common mixture in this "
            "project's AMR literature ('Development of a novel rotary "
            "magnetic refrigerator'; high-frequency rotating AMR paper)."
        ),
    },
    "water_pg30": {
        "rho": 1013.0, "cp": 3800.0, "mu": 2.3e-3, "k": 0.47,
        "freeze_C": -13.0,
        "description": (
            "70 vol% water / 30 vol% propylene glycol -- non-toxic "
            "corrosion-inhibited alternative to ethylene glycol. Not "
            "reported in this project's prototype corpus; included as a "
            "low-toxicity design option."
        ),
    },
    "ethanol": {
        "rho": 789.0, "cp": 2440.0, "mu": 1.2e-3, "k": 0.167,
        "freeze_C": -114.0,
        "description": (
            "Pure ethanol -- only relevant for sub-zero-span AMR "
            "variants; poor cp/k make it a weak choice at room "
            "temperature. Included for completeness/comparison only, not "
            "used by any prototype in this project's corpus."
        ),
    },
}

FLUID_NAMES = tuple(FLUID_LIBRARY)

# Set from core/fluid_selection_optimization.py's comparison at this
# project's standard baseline operating point (T_COLD_K=291, SPAN_K=10,
# see core/optimize.py) -- see results/fluid_selection_optimization.txt
# for the numeric COP/Qc/pumping-power trade-off this choice is based on.
#
# RESULT: pure water actually comes out ~1.4% ahead of water_eg10 on
# COP_electrical at this operating point (real, not negligible, but
# small) -- see the results file. Pure water is NOT a viable choice for
# real hardware, though: it is not used by any actual AMR prototype in
# this project's Papers/ corpus, all of which run a glycol mixture for
# corrosion protection (see FLUID_LIBRARY entries above). Of the
# corrosion-protected options, water_eg10 (10% ethylene glycol) gives up
# the least COP_electrical relative to the (unrealistic) pure-water
# ceiling -- less than water_eg20, water_pg30, or ethanol -- so it is
# the recommended REALISTIC default here.
#
# NOTE: every existing call site in this repo pins fluid="water" or
# fluid_name="water" EXPLICITLY as its own function-level default (see
# thermal.py/regenerator_1d.py/fluid_mce_cycle.py/geometry_analysis.py/
# amr_cycle.py) rather than reading this constant. That is intentional:
# fluid_cp feeds directly into Qc (Watts) for every AMRSystem caller,
# not just NTU-thermal-model ones, so flipping every function's default
# to this constant would shift every existing calibrated Qc/COP number
# in the repo (validation_system.py's Watt-level fits against real
# devices included) by several percent -- a real behavior change, not a
# no-op, and inconsistent with this repo's own "opt-in, don't silently
# change existing numbers" convention (see e.g. amr_cycle.py's
# pump_motor_efficiency docstring for the same pattern elsewhere).
# Callers doing NEW, realistic system-level work should pass
# fluid="water_eg10" (or reference this constant) explicitly.
DEFAULT_FLUID = "water_eg10"


def fluid_properties(fluid="water", T_K=300.0):
    """Returns {rho, cp, mu, k} (kg/m^3, J/(kg K), Pa s, W/(m K)) for the
    named fluid. `fluid` must be a key of FLUID_LIBRARY -- raises
    ValueError with the valid options otherwise, rather than silently
    falling back to water. T_K is accepted for interface compatibility
    with the original water_properties(T_K) and for future
    temperature-dependent correlations; it does not currently vary the
    returned values (see module docstring)."""
    if fluid not in FLUID_LIBRARY:
        raise ValueError(
            f"unknown fluid {fluid!r}; options: {sorted(FLUID_LIBRARY)}")
    entry = FLUID_LIBRARY[fluid]
    return {"rho": entry["rho"], "cp": entry["cp"],
            "mu": entry["mu"], "k": entry["k"]}
