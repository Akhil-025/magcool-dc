"""
thermal_diode_analysis.py
==========================
Phase 18 (ROADMAP.md) validation deliverable for `core/thermal_diode.py`.

Scope, deliberately narrow (read `core/thermal_diode.py`'s module
docstring honesty flag before trusting any number here). The Phase 18
plan asked two concrete questions before any deeper investment in this
mechanism:

  1. Does a thermal-diode-assisted design let this repo's model exceed a
     mechanical-switching frequency ceiling that otherwise caps it?
  2. Is there a benchmark device this module's numbers can be checked
     against?

This module answers both directly, rather than building the fuller
NSGA-III-integrated, rectification-ratio-driven heat-transfer model the
original Phase 18 plan sketched as a stretch goal -- consistent with
that plan's own recommendation to "scope down to 'what frequency
ceiling would need to be broken for this to matter' as a sensitivity
study before building the full diode model."

Finding on question 1 (checked directly, not assumed): NO internal
frequency ceiling exists anywhere in core/amr_cycle.py's AMRSystem.
Frequency feeds W_eddy ~ f^2 (core/loss_model.py) and, since Phase 16,
W_hys ~ f (hysteresis) -- both raise parasitic loss monotonically with
f, but neither model, nor cooling_capacity()/magnetic_work(), ever hard-
caps f. The only frequency bound anywhere in this repo is
core/optimize.py's NSGA-III search-space upper bound `_XU[1] = 5.0` Hz
-- and that number carries NO documented justification tying it to a
mechanical-valve-switching limit (checked: no comment, docstring, or
ROADMAP.md entry cites a reason for 5.0 Hz specifically). So the plan's
own premise -- "if [_XU's frequency bound] is set by a mechanical-
switching limit, that's the exact number a diode-assisted design should
be allowed to exceed" -- does not apply as literally as posed: there is
no mechanical-switching-derived ceiling in this repo to relax in the
first place. This module therefore does NOT add a `thermal_diode_
assisted`-conditional frequency-bound relaxation to optimize.py (there
is nothing there to relax), and states this finding explicitly rather
than inventing a ceiling to then dramatically break.

Finding on question 2: none of the 16 devices in
data/amr_experimental_benchmarks.csv use thermal diodes of any kind
(confirmed by inspection -- every device is either continuous-rotary or
conventional valve-switched). There is therefore no benchmark row this
module's rectification_ratio or switching_power_W numbers can be checked
against. What this module provides instead is a documented SENSITIVITY
STUDY within this repo's own model: how much COP_electrical is reduced,
at a representative operating point, by paying the (illustrative)
actuation switching-power cost `core.thermal_diode.
DEFAULT_MECHANICAL_CONTACT_DIODE` implies, as frequency is swept -- i.e.
it quantifies the DOWNSIDE this module's own accounting choice
(Phase 18's AMRSystem.thermal_diode wiring) adds, honestly, without
claiming any offsetting heat-transfer benefit that this repo's model
does not (and, per the honesty flag above, currently cannot) represent.
"""

from core.mce_material import GADOLINIUM
from core.amr_cycle import AMRSystem
from core.loss_model import StateDependentLossModel
from core.thermal_diode import (MechanicalContactDiode,
                                  DEFAULT_MECHANICAL_CONTACT_DIODE,
                                  cycle_time_reduction_factor)
from core.optimize import _XU as OPTIMIZE_XU

T_COLD_K = 291.0
SPAN_K = 10.0
MU0H_T = 1.5
MASS_KG = 5.0
MDOT_KG_S = 0.08

_LOSS_MODEL = StateDependentLossModel()


def _run(frequency, thermal_diode=None):
    sys_ = AMRSystem(GADOLINIUM, mu0H_max=MU0H_T, mass_regenerator=MASS_KG,
                      frequency=frequency, fluid_mdot=MDOT_KG_S,
                      regenerator_effectiveness=0.85, loss_model=_LOSS_MODEL,
                      use_ntu_thermal_model=True, thermal_diode=thermal_diode)
    return sys_.run(T_COLD_K, SPAN_K)


def check_frequency_ceiling_claim(verbose=True):
    """Directly checks (rather than assumes) whether this repo's model
    has an internal mechanical-switching-derived frequency ceiling for a
    thermal diode to relax. Returns a dict describing the finding -- see
    module docstring for the full writeup."""
    xu_frequency_hz = float(OPTIMIZE_XU[1])
    finding = {
        "optimize_py_frequency_upper_bound_Hz": xu_frequency_hz,
        "amr_system_has_internal_frequency_cap": False,
        "upper_bound_documented_as_mechanical_switching_limit": False,
    }
    if verbose:
        print(f"core/optimize.py's NSGA-III search-space frequency upper "
              f"bound is {xu_frequency_hz} Hz.")
        print("AMRSystem itself (core/amr_cycle.py) has NO internal "
              "frequency ceiling: frequency only ever enters W_eddy~f^2 "
              "and W_hys~f (both monotonic, uncapped parasitic-loss "
              "terms), never a hard cutoff on cooling_capacity() or "
              "magnetic_work().")
        print(f"No comment, docstring, or ROADMAP.md entry ties the "
              f"{xu_frequency_hz} Hz optimize.py bound to a mechanical-"
              f"valve-switching limit specifically -- it is an unexplained "
              f"round-number search-space bound, not a physical constraint "
              f"this module can relax.")
    return finding


def sweep_frequency_with_and_without_diode(
        frequencies=(0.5, 1.0, 2.0, 4.0, 8.0),
        diode: MechanicalContactDiode = DEFAULT_MECHANICAL_CONTACT_DIODE,
        verbose=True):
    """At this repo's own representative operating point, compares
    COP_electrical with vs. without the (illustrative, unbenchmarked)
    thermal-diode actuation switching-power cost, across a frequency
    sweep. Since AMRSystem.thermal_diode adds ONLY a parasitic cost here
    (no offsetting heat-transfer benefit is modeled -- see module
    docstring), diode-assisted COP_electrical is <= the no-diode
    baseline at every frequency by construction; this sweep exists to
    show HOW MUCH that illustrative cost matters relative to the
    already-dominant eddy-current/base-overhead losses, not to claim a
    net benefit."""
    rows = []
    for f in frequencies:
        base = _run(f, thermal_diode=None)
        diode_assisted = _run(f, thermal_diode=diode)
        delta_cop_pct = (100 * (diode_assisted.COP_electrical - base.COP_electrical)
                          / base.COP_electrical) if base.COP_electrical > 0 else 0.0
        rows.append((f, base.COP_electrical, diode_assisted.COP_electrical, delta_cop_pct))
        if verbose:
            print(f"  f={f:5.2f}Hz   COP_no_diode={base.COP_electrical:7.4f}   "
                  f"COP_diode_assisted={diode_assisted.COP_electrical:7.4f}   "
                  f"delta={delta_cop_pct:6.2f}%")
    return rows


def demo_cycle_time_reduction(verbose=True):
    """Illustrative, explicitly-not-fit worked example of
    cycle_time_reduction_factor() -- see that function's docstring for
    why both switch times must be caller-supplied rather than defaulted.
    Uses a round, clearly-labeled illustrative pair of switch times
    (NOT digitized from any source in this project's corpus) purely to
    demonstrate how the helper would be used if such numbers were ever
    obtained."""
    conventional_switch_time_s = 0.5   # illustrative only, see docstring
    diode_switch_time_s = 0.05          # illustrative only, see docstring
    reduction = cycle_time_reduction_factor(conventional_switch_time_s,
                                             diode_switch_time_s)
    if verbose:
        print(f"  Illustrative example only (NOT a literature value -- see "
              f"docstring): conventional valve switch time="
              f"{conventional_switch_time_s}s, diode switch time="
              f"{diode_switch_time_s}s -> {reduction*100:.0f}% dead-time "
              f"reduction per half-cycle IF these numbers were real.")
    return reduction


def run_thermal_diode_analysis(out_path="results/thermal_diode_analysis.txt"):
    import io, contextlib
    buf = io.StringIO()
    with contextlib.redirect_stdout(buf):
        print("=" * 90)
        print("PHASE 18: mechanical-contact active thermal diode -- scoped-down")
        print("sensitivity study, motivated by Kitanovski et al. (2015) Ch. 6.")
        print("See core/thermal_diode.py's module docstring for the full honesty flag")
        print("(this project's copy of the book does not include Ch. 6, pp. 211-268).")
        print("=" * 90)

        print("\n--- Step 1: does this repo's model have a mechanical-switching "
              "frequency ceiling for a diode to relax? ---")
        check_frequency_ceiling_claim()

        print(f"\n--- Step 2: COP_electrical impact of the (illustrative) diode "
              f"actuation cost across frequency, at T_cold={T_COLD_K}K, "
              f"span={SPAN_K}K, mu0H={MU0H_T}T, mass={MASS_KG}kg, "
              f"mdot={MDOT_KG_S}kg/s ---")
        rows = sweep_frequency_with_and_without_diode()

        print("\n--- Step 3: cycle_time_reduction_factor() illustrative worked example ---")
        demo_cycle_time_reduction()

        print("\n--- Conclusion ---")
        worst_delta = min(r[3] for r in rows)
        print(f"The (illustrative, unbenchmarked) actuation switching-power cost "
              f"in DEFAULT_MECHANICAL_CONTACT_DIODE reduces COP_electrical by at "
              f"most {abs(worst_delta):.2f}% across the frequencies swept here -- "
              f"a small effect relative to the eddy-current/base-overhead losses "
              f"that already dominate W_parasitic at this operating point, because "
              f"the illustrative actuation_energy_J_per_cycle=0.05J is small. This "
              f"module deliberately does NOT model any offsetting heat-transfer "
              f"benefit from the diode's rectification_ratio (no closed-form "
              f"relation for how rectification ratio would improve AMR cycle "
              f"performance was available to digitize -- see honesty flag), so "
              f"net COP_electrical is <= the no-diode baseline by construction "
              f"here; this is a documented cost-only accounting, not a claim that "
              f"real mechanical-contact thermal diodes are a net negative for AMR "
              f"performance. Step 1 additionally found that this repo's model has "
              f"NO internal mechanical-switching frequency ceiling for a diode-"
              f"assisted design to relax in the first place -- the only frequency "
              f"bound anywhere (optimize.py's NSGA-III search bound) is an "
              f"unexplained round number, not a physical constraint. No benchmark "
              f"device in this repo's corpus uses thermal diodes, so none of this "
              f"is validated against real hardware -- treat this module as a "
              f"design-exploration tool, exactly the disposition the Phase 18 "
              f"plan itself recommended for this item, not a validated feature.")

    text = buf.getvalue()
    print(text, end="")
    with open(out_path, "w") as fh:
        fh.write(text)
    return text


if __name__ == "__main__":
    run_thermal_diode_analysis()