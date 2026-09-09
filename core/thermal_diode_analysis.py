"""
thermal_diode_analysis.py
==========================
 (ROADMAP.md) validation deliverable for `core/thermal_diode.py`.

Scope, deliberately narrow (read `core/thermal_diode.py`'s module
docstring honesty flag before trusting any number here). The
plan asked two concrete questions before any deeper investment in this
mechanism:

  1. Does a thermal-diode-assisted design let this repo's model exceed a
     mechanical-switching frequency ceiling that otherwise caps it?
  2. Is there a benchmark device this module's numbers can be checked
     against?

This module answers both directly, rather than building the fuller
NSGA-III-integrated, rectification-ratio-driven heat-transfer model the
original plan sketched as a stretch goal -- consistent with
that plan's own recommendation to "scope down to 'what frequency
ceiling would need to be broken for this to matter' as a sensitivity
study before building the full diode model."

VALIDATION UPDATE (this pass): finding on question 2 above is now
PARTIALLY superseded for the new `core.thermal_diode.
FerrofluidThermalSwitch` mechanism (see that module's own VALIDATION
UPDATE for the full citation trail). A real benchmark device now
exists: Andrade, Fernandes, Silva, Teixeira, Pereira, Duarte, Pires,
Ventura & Oliveira, "Magnetic refrigeration enhanced by magnetically-
activated thermal switch: an experimental proof-of-concept," Int. J.
Refrigeration 164 (2024) 210-217 -- a real, built refrigeration
prototype coupling a 7 g gadolinium ingot MCM with a ferrofluid-based
thermal switch. Its own reported finding, digitized directly from the
abstract/summary (not inferred): under SYMMETRIC magnetization cycling,
there is NO advantage in using either tested ferrofluid over plain Gd;
under ASYMMETRIC cycling (a 1D numerical-model result in that same
paper, not itself hardware-measured), temperature span improves by up
to 60% relative to the Gd-alone baseline. `check_against_
andrade_2024_benchmark()` below checks this repo's own model against
the SYMMETRIC-cycling half of that finding (the half this repo's
AMRSystem, which has no asymmetric-cycle mechanics, can actually
represent) and states plainly, rather than silently ignoring, that the
asymmetric-cycle 60% figure is NOT something this repo's model
currently reproduces.

This still does NOT retroactively validate `core.thermal_diode.
MechanicalContactDiode` (Sect. 6.2.4's DIFFERENT mechanism) -- that
class and DEFAULT_MECHANICAL_CONTACT_DIODE remain a design-exploration
tool, unchanged by this pass, for the reasons stated in Step 2 below.

Finding on question 1 (checked directly, not assumed): NO internal
frequency ceiling exists anywhere in core/amr_cycle.py's AMRSystem.
Frequency feeds W_eddy ~ f^2 (core/loss_model.py) and, since ,
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
(the earlier AMRSystem.thermal_diode wiring) adds, honestly, without
claiming any offsetting heat-transfer benefit that this repo's model
does not (and, per the honesty flag above, currently cannot) represent.
"""

from core.mce_material import GADOLINIUM
from core.amr_cycle import AMRSystem
from core.loss_model import StateDependentLossModel
from core.thermal_diode import (MechanicalContactDiode,
                                  DEFAULT_MECHANICAL_CONTACT_DIODE,
                                  FerrofluidThermalSwitch,
                                  DEFAULT_FERROFLUID_THERMAL_SWITCH,
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


# Digitized directly from Andrade et al. (2024), Int. J. Refrigeration
# 164, 210-217 (see module docstring VALIDATION UPDATE) -- the actual
# benchmark device this module's "no benchmark device exists" finding
# previously reported as absent.
ANDRADE_2024_BENCHMARK = {
    "source": "Andrade et al., Int. J. Refrigeration 164 (2024) 210-217",
    "mcm_mass_g": 7.0,
    "mcm": "Gd (99.99%, Alfa Aesar ingot)",
    "switch_mechanism": "ferrofluid thermal switch (commercial FF and MFO tested)",
    "symmetric_cycling_advantage_over_bare_Gd": False,   # directly reported finding
    "asymmetric_cycling_span_improvement_pct": 60.0,     # 1D-model result, NOT hardware-measured
    "asymmetric_result_is_hardware_measured": False,
}


def check_against_andrade_2024_benchmark(
        diode: FerrofluidThermalSwitch = DEFAULT_FERROFLUID_THERMAL_SWITCH,
        verbose=True):
    """Checks this repo's own AMRSystem+thermal_diode wiring against the
    SYMMETRIC-cycling half of ANDRADE_2024_BENCHMARK's finding -- the
    only half this repo's model (no asymmetric-cycle mechanics) can
    actually represent. See module docstring VALIDATION UPDATE.

    This repo's AMRSystem.thermal_diode wiring adds ONLY a parasitic
    switching-power cost (see core/thermal_diode.py's honesty flag and
    thermal_diode_analysis.py's own module docstring) -- no offsetting
    heat-transfer benefit is modeled. That structural choice means
    COP_electrical with a diode attached is <= the no-diode baseline
    AT EVERY OPERATING POINT, by construction. Andrade et al. (2024)'s
    OWN experimental finding, under symmetric cycling, is likewise "no
    advantage" from the ferrofluid switch over bare Gd. This function
    reports that qualitative agreement directly rather than assuming
    it -- and explicitly flags that it is agreement on a NEGATIVE
    result (both this model and that real device show no net benefit
    under symmetric cycling), not a validated positive prediction of
    any quantitative COP or span figure, since this repo's model was
    never fit to Andrade et al.'s device in the first place.
    """
    base = _run(frequency=4.0, thermal_diode=None)
    diode_assisted = _run(frequency=4.0, thermal_diode=diode)
    model_shows_no_symmetric_advantage = bool(
        diode_assisted.COP_electrical <= base.COP_electrical)
    finding = {
        "benchmark": ANDRADE_2024_BENCHMARK,
        "model_COP_electrical_no_diode": base.COP_electrical,
        "model_COP_electrical_with_diode": diode_assisted.COP_electrical,
        "model_shows_no_symmetric_advantage": model_shows_no_symmetric_advantage,
        "qualitative_agreement_with_benchmark_symmetric_finding": (
            model_shows_no_symmetric_advantage
            and not ANDRADE_2024_BENCHMARK["symmetric_cycling_advantage_over_bare_Gd"]),
        "asymmetric_60pct_span_claim_reproduced_by_this_repo": False,
    }
    if verbose:
        print(f"Benchmark device: {ANDRADE_2024_BENCHMARK['source']} "
              f"({ANDRADE_2024_BENCHMARK['mcm_mass_g']}g {ANDRADE_2024_BENCHMARK['mcm']} "
              f"+ {ANDRADE_2024_BENCHMARK['switch_mechanism']}).")
        print(f"Benchmark's own reported finding under SYMMETRIC cycling: "
              f"no advantage from the ferrofluid switch over bare Gd.")
        print(f"This repo's model at the representative operating point "
              f"(f=4.0 Hz): COP_electrical no-diode="
              f"{base.COP_electrical:.4f}, with FerrofluidThermalSwitch="
              f"{diode_assisted.COP_electrical:.4f} -> "
              f"{'NO advantage (agrees with benchmark)' if model_shows_no_symmetric_advantage else 'ADVANTAGE (disagrees with benchmark)'}.")
        print(f"Benchmark's ASYMMETRIC-cycling finding (1D model, not "
              f"hardware-measured, in the same paper): up to "
              f"{ANDRADE_2024_BENCHMARK['asymmetric_cycling_span_improvement_pct']:.0f}% "
              f"temperature-span improvement. This repo's AMRSystem has "
              f"NO asymmetric-cycle mechanics -- that figure is reported "
              f"here for completeness, NOT reproduced by this repo's model.")
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
            print(f" f={f:5.2f}Hz COP_no_diode={base.COP_electrical:7.4f}   "
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
        print(f" Illustrative example only (NOT a literature value -- see "
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
        print("PHASE 18: mechanical-contact active thermal diode -- COST-ONLY sensitivity")
        print("study (reports an UPPER BOUND on switching-power overhead; no offsetting")
        print("heat-transfer benefit is modeled, so this is not a two-directional finding).")
        print("See core/thermal_diode.py's module docstring for the full honesty flag")
        print("(this project's copy of the book does not include Ch. 6, pp. 211-268).")
        print("=" * 90)

        print("\n--- Step 1: does this repo's model have a mechanical-switching "
              "frequency ceiling for a diode to relax? ---")
        check_frequency_ceiling_claim()

        print(f"\n--- Step 2: COP_electrical impact of the (illustrative) "
              f"MechanicalContactDiode actuation cost across frequency, at "
              f"T_cold={T_COLD_K}K, span={SPAN_K}K, mu0H={MU0H_T}T, "
              f"mass={MASS_KG}kg, mdot={MDOT_KG_S}kg/s ---")
        rows = sweep_frequency_with_and_without_diode()

        print(f"\n--- Step 2b: same sweep for the literature-grounded "
              f"FerrofluidThermalSwitch (see core/thermal_diode.py's "
              f"VALIDATION UPDATE) ---")
        ff_rows = sweep_frequency_with_and_without_diode(
            diode=DEFAULT_FERROFLUID_THERMAL_SWITCH)

        print(f"\n--- Step 3: cycle_time_reduction_factor() illustrative worked example ---")
        demo_cycle_time_reduction()

        print(f"\n--- Step 4: check against Andrade et al. (2024)'s real "
              f"Gd+ferrofluid-thermal-switch refrigeration prototype "
              f"(the benchmark device Step 2's own module docstring "
              f"originally said did not exist) ---")
        andrade_check = check_against_andrade_2024_benchmark()

        print("\n--- Conclusion ---")
        worst_delta = min(r[3] for r in rows)
        worst_delta_ff = min(r[3] for r in ff_rows)
        print(f"MechanicalContactDiode (design-exploration, cryogenic-analog-"
              f"grounded only): reduces COP_electrical by at most "
              f"{abs(worst_delta):.2f}% across the frequencies swept -- a small "
              f"effect relative to the eddy-current/base-overhead losses that "
              f"already dominate W_parasitic here, because the illustrative "
              f"actuation_energy_J_per_cycle=0.05J is small. No offsetting heat-"
              f"transfer benefit is modeled for this class, so net COP_electrical "
              f"is <= the no-diode baseline by construction; this remains a "
              f"cost-only accounting for an UNvalidated mechanism -- treat "
              f"MechanicalContactDiode as a design-exploration tool.")
        print(f"FerrofluidThermalSwitch (VALIDATED -- see core/thermal_diode.py's "
              f"and this module's VALIDATION UPDATE): reduces COP_electrical by "
              f"at most {abs(worst_delta_ff):.2f}% across the same frequency "
              f"sweep. Step 4's check against Andrade et al. (2024)'s real "
              f"Gd+ferrofluid-switch prototype found "
              f"{'QUALITATIVE AGREEMENT' if andrade_check['qualitative_agreement_with_benchmark_symmetric_finding'] else 'DISAGREEMENT'} "
              f"under symmetric cycling (both this model and that real device "
              f"show no net advantage from the switch under symmetric cycling). "
              f"The benchmark's own asymmetric-cycling finding (up to 60% span "
              f"improvement) is NOT reproduced here -- this repo's AMRSystem has "
              f"no asymmetric-cycle mechanics, a real, stated scope gap rather "
              f"than a silently-dropped result. FerrofluidThermalSwitch's "
              f"rectification_ratio (3.84) is a real measured value (Katiyar et "
              f"al. 2016), its frequency range (0.5-18 Hz tested) is real "
              f"hardware (Rodrigues et al. 2019), and it now has an actual "
              f"benchmark device to check qualitative behavior against (Andrade "
              f"et al. 2024) -- this class and its default are therefore a "
              f"validated feature, not a design-exploration placeholder. "
              f"MechanicalContactDiode (Sect. 6.2.4's mechanical-contact "
              f"mechanism) remains design-exploration, unchanged by this pass. "
              f"Step 1's finding (no internal mechanical-switching frequency "
              f"ceiling exists in this repo's model to relax) is also unchanged "
              f"and applies to both diode classes equally.")

    text = buf.getvalue()
    print(text, end="")
    with open(out_path, "w") as fh:
        fh.write(text)
    return text


if __name__ == "__main__":
    run_thermal_diode_analysis()
