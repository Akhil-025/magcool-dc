"""
electrocaloric_cycle.py
=========================
Device-level electrical-COP model for electrocaloric (PST relaxor-
ferroelectric MLCC) cooling, added to sit alongside barocaloric_cycle.py
and elastocaloric_cycle.py in alternative_caloric_comparison.py.

WHY THIS MODULE LOOKS DIFFERENT FROM THE OTHER TWO
----------------------------------------------------
barocaloric_cycle.py and elastocaloric_cycle.py both build UP from a
per-unit-mass material model (Delta_S, Delta_T_ad(driving field), a
hysteresis-loss term) to a device COP. Doing the same thing here would
require an independently-verified Delta_S-vs-field curve for the exact
PST MLCC batches used in the papers below -- this project's search did
NOT locate that (the original papers' entropy/field data was not
directly accessible, only their headline device-level results). Rather
than back-fill a plausible-looking Delta_S and c_p to force a "correct"
per-stage number (which would be fabricating a number precise enough to
look verified when it isn't -- exactly the trap this repo's other
material modules explicitly flag against, e.g. barocaloric_material.py's
c_p flag), this module works top-down instead: it takes the papers' own
reported DEVICE-level span/power/COP figures as the source of truth and
fits a single scaling law to them, following this repo's own documented
observation about its OWN vapor-compression model (see
alternative_caloric_comparison.py's docstring: "Your repo's own VCC
model runs consistently at about 42% of Carnot across the whole 5-20K
span") -- i.e. constant fraction-of-Carnot efficiency is already this
repo's chosen way to extrapolate a single-technology COP across a span
range, not a new assumption invented for electrocaloric specifically.

Source (measured, not simulated, this is the ONE calibration point)
----------------------------------------------------------------------
Li, W. et al., "A high-performance electrocaloric refrigerator", Science
382, 6669 (2023). Double-loop electrocaloric heat pump, PST MLCC stack.
  - Span: 20.9 K (covers the full 5-20 K ASHRAE range studied elsewhere
    in this repo in a single measured device).
  - Cooling power: 2.1 W. The paper's originally reported figure was
    4.2 W; a 2025 erratum halved this after a measurement correction --
    this repo uses the CORRECTED (post-erratum) 2.1 W figure and flags
    the correction explicitly, since it shows even this paper's own
    metrology needed a public fix.
  - COP: reaches 54% of Carnot efficiency (the paper's own abstract
    figure, after accounting for fluid-pumping energy, PROVIDED the
    electrical work is properly recovered). A EurekAlert press release
    quoted 64% -- this module uses the paper's own 54% figure as the
    more trustworthy of the two, and flags the discrepancy.

Other electrocaloric device papers located (same PST MLCC family,
recorded as claims in alternative_caloric_comparison.py rather than
used for calibration here, since neither reports a clean end-to-end COP
number this module can fit to):
  - Torelló, A. et al., Science 370, 6512 (2020): active electrocaloric
    regenerator, 13.0 K span, regeneration factor 6, but measured at only
    0.08 Hz -- specific cooling power just 0.012 W/g, i.e. COP achieved
    by starving the device of power (the same pattern this repo's
    barocaloric_cycle.py flags for its own 1 mHz simulated target).
  - Meng, Y. et al., Nature Energy 5 (2020) 996: cascade electrocaloric
    cooling device (several material stages, each covering a slice of
    the total span) -- the electrocaloric analogue of this repo's own
    Curie-graded cascade.py, but this project's search did not locate a
    clean device-level COP figure to compare here.

MODEL AND ITS LIMITS (read before trusting a number away from 20.9 K)
------------------------------------------------------------------------
    COP_electrical(span) = measured_fraction_of_carnot * T_cold / span

This reproduces the ONE measured point (span=20.9K) by construction --
it is NOT validated at any other span (unlike elastocaloric_cycle.py,
which is calibrated against TWO independently measured benchmarks). Any
span other than 20.9 K in run_span_sweep() below is an EXTRAPOLATION
under the assumption that this device's fraction-of-Carnot efficiency
holds outside the demonstrated range; it is flagged as such in every
result row via `is_extrapolated`.

A SECOND, STRONGER FLAG for narrow spans: the other electrocaloric
device this project located (Torello et al. 2020) demonstrated a 13.0 K
span, so [13.0, 20.9] K is at least bracketed by two real electrocaloric
devices even though only one reports a COP. Below 13.0 K there is no
electrocaloric device evidence at all -- the fraction-of-Carnot
assumption is doing much more work there, and the resulting numbers
(e.g. COP > 30 at a 5 K span) should be read as considerably LESS
trustworthy than the already-uncertain 13-20.9 K range, not as a
straightforward "electrocaloric wins at every span" result. Rows with
span < 13.0 K are additionally flagged via `is_far_extrapolation`.
Cooling POWER (2.1 W) is not scaled by
this model at all -- it stays a fixed, tiny, device-scale number
reported alongside COP, and is NOT claimed to hold at any other span or
to represent anything close to data-center kW-scale cooling. Treat this
module the way barocaloric_cycle.py asks its own single-point
calibration to be treated: as the best available honest estimate from
one real data point, not as a validated span-dependent physics model.
"""

from dataclasses import dataclass
from typing import Optional


@dataclass
class ElectrocaloricCycleResult:
    span_K: float
    T_cold: float
    carnot_cop: float
    fraction_of_carnot: float
    COP_electrical: float
    is_measured_span: bool          # True only at the single calibration span
    is_extrapolated: bool           # True everywhere else
    is_far_extrapolation: bool      # True below 13.0K -- outside ANY device evidence, see docstring
    measured_cooling_power_W: Optional[float]  # only populated at the measured span
    source: str


# Lower bound of the span range actually bracketed by real electrocaloric
# device evidence (Torello et al. 2020's 13.0K span, even though that paper
# doesn't report a clean COP). Below this, the constant-fraction-of-Carnot
# assumption has no device evidence backing it at all -- see docstring.
DEMONSTRATED_SPAN_FLOOR_K = 13.0


# The ONE measured, corrected (post-erratum) calibration point.
MEASURED_DEVICE = {
    "span_K": 20.9,
    "cooling_power_W": 2.1,            # corrected value, see docstring (was 4.2 W pre-erratum)
    "cooling_power_W_pre_erratum": 4.2,
    "fraction_of_carnot": 0.54,        # paper's own abstract figure
    "fraction_of_carnot_press_release": 0.64,  # EurekAlert figure -- flagged discrepancy, not used
    "source": "Li, W. et al., Science 382, 6669 (2023); cooling-power "
              "erratum (2025) halved the original 4.2 W figure to 2.1 W.",
}


def run_cycle(span_K, T_cold=291.15, fraction_of_carnot=None):
    """COP_electrical(span) = fraction_of_carnot * T_cold / span, per the
    module docstring. Defaults to the paper's own 54% figure, not the
    64% press-release figure -- see docstring."""
    if fraction_of_carnot is None:
        fraction_of_carnot = MEASURED_DEVICE["fraction_of_carnot"]

    carnot_cop = T_cold / span_K
    cop_electrical = fraction_of_carnot * carnot_cop

    is_measured_span = abs(span_K - MEASURED_DEVICE["span_K"]) < 1e-9
    measured_power = MEASURED_DEVICE["cooling_power_W"] if is_measured_span else None

    return ElectrocaloricCycleResult(
        span_K=span_K,
        T_cold=T_cold,
        carnot_cop=carnot_cop,
        fraction_of_carnot=fraction_of_carnot,
        COP_electrical=cop_electrical,
        is_measured_span=is_measured_span,
        is_extrapolated=not is_measured_span,
        is_far_extrapolation=(span_K < DEMONSTRATED_SPAN_FLOOR_K),
        measured_cooling_power_W=measured_power,
        source=MEASURED_DEVICE["source"],
    )


def run_span_sweep(spans_K=(5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20),
                    T_cold=291.15, fraction_of_carnot=None, verbose=True):
    rows = [run_cycle(span, T_cold=T_cold, fraction_of_carnot=fraction_of_carnot)
            for span in spans_K]

    if verbose:
        print("Electrocaloric (PST MLCC) COP_electrical vs. span "
              "(calibrated against ONE MEASURED device at 20.9 K -- "
              "every other row is an EXTRAPOLATION, see module docstring):")
        for r in rows:
            if r.is_measured_span:
                tag = "MEASURED"
            elif r.is_far_extrapolation:
                tag = "FAR EXTRAPOLATION -- below any device evidence"
            else:
                tag = "extrapolated"
            print(f"  span={r.span_K:5.1f}K  Carnot={r.carnot_cop:6.2f}  "
                  f"COP_electrical={r.COP_electrical:6.2f}  [{tag}]")

    return rows


def sensitivity_to_press_vs_paper_figure(span_K=20.9, T_cold=291.15):
    """Reports the COP under both the paper's own 54% figure and the
    64% figure quoted in secondary press coverage, so downstream code
    doesn't have to silently pick one -- see module docstring on why
    54% (the paper's own number) is the default used everywhere else."""
    paper = run_cycle(span_K, T_cold=T_cold,
                       fraction_of_carnot=MEASURED_DEVICE["fraction_of_carnot"])
    press = run_cycle(span_K, T_cold=T_cold,
                       fraction_of_carnot=MEASURED_DEVICE["fraction_of_carnot_press_release"])
    return {
        "paper_54pct": paper.COP_electrical,
        "press_release_64pct": press.COP_electrical,
    }


if __name__ == "__main__":
    run_span_sweep()
    s = sensitivity_to_press_vs_paper_figure()
    print(f"\nAt the measured 20.9K span: paper's own 54% figure -> COP="
          f"{s['paper_54pct']:.2f}; EurekAlert's 64% figure -> COP="
          f"{s['press_release_64pct']:.2f} (discrepancy flagged, paper figure used by default).")
