# Findings note: Hypereg high-frequency regenerator

**Source read directly for this pass**: Klinar, K. et al., "Perspectives
and Energy Applications of Magnetocaloric, Pyromagnetic, Electrocaloric
and Pyroelectric Materials", *Advanced Energy Materials* 14, 2401739
(2024) — `Papers/Reviews/Advanced Energy Materials - 2024 - Klinar -
Perspectives and Energy Applications of Magnetocaloric  Pyromagnetic
.pdf`, Section "Future Heat Transfer and Regenerator Principles"
(Figs. 18–21, main text pp. 26–29 of the PDF).

Only the Hypereg-relevant pages were read in depth for this pass, per the
plan's own scope ("read the Hypereg-relevant section directly"), not the
full 36-page review.

## The question this note answers

Does a Hypereg-style design reduce the eddy-current loss coefficient
`k_eddy` (`core/loss_model.py`'s domain), or does it tolerate higher
frequency because of a different **heat-exchange / hydraulic geometry**
(`core/thermal.py`'s domain)? These are different code changes, and the
plan explicitly asked this be resolved before touching either module.

## What the paper actually says

Hypereg is presented as the review authors' **own, newly-patented**
concept (patent citation [275] in the review), described in their own
words as unveiled "for the first time" in this review — there is no
separate, independently peer-reviewed Hypereg paper to cross-check
against, and no built prototype is described. The mechanism, in the
authors' own framing, is purely **hydraulic**: conventional regenerators
pass fluid through one bed in a long **series** path (pressure-drop
length `L_dp`); Hypereg instead splits the bed into several shorter
**sub-regenerators** connected **in parallel**, fed by a single shared
propulsion system, with fluid oscillating through each. Because pressure
drop scales with flow-path length, splitting into *n* parallel sub-beds
of equal total mass reduces the pressure-drop length to roughly `L_dp/n`
— the paper's own illustrative Figure 19 example uses `n=4`.
Figures 20–21 sketch several possible ways to realize the parallel split
(Tesla-valve-style fluidic diodes, resonant cavities, sweeping fluidic
oscillators, or a cross-flow layered arrangement with a separate
unidirectional carrier fluid) but do not analyze or compare them
quantitatively.

**Nowhere does the paper attribute the frequency benefit to a change in
eddy-current losses, magnet-circuit design, or any electromagnetic
mechanism.** The stated goal is explicitly to enable higher-frequency
operation (the review's own heading: "high-frequency (<25 Hz for
liquids) active regeneration principle") by cutting the **pumping-power**
penalty that would otherwise make high-frequency (and therefore
high-flow) operation impractical — not by suppressing eddy currents in
the magnet/regenerator support structure, which is a magnetic-circuit
effect this hydraulic-only concept does not touch.

## Answer

This is a **`core/thermal.py`-domain change** (pressure drop / pumping
power), **not** a `core/loss_model.py` change. No eddy-current
recalibration is justified by this source. `k_eddy` in
`core.loss_model.StateDependentLossModel` is left exactly as calibrated.

## What was (and wasn't) implemented as a result

- **Implemented**: `core.thermal.pumping_power_packed_bed_hypereg()` —
  identical to the existing `pumping_power_packed_bed()` but with the
  pressure-drop length divided by a configurable
  `n_parallel_subregenerators` (default 4, matching the paper's own
  Figure 19 example). See that function's docstring for the full set of
  honesty flags — most importantly, this is a **pumping-power-only**
  estimate; heat-transfer effectiveness is unchanged (splitting one bed
  into equal-mass parallel sub-beds does not, in this 0-D model, change
  NTU), and no cost/complexity penalty for building the n-way
  parallelization hardware itself (Fig. 20's fluidic diodes/oscillators)
  is modeled.
- **Deliberately NOT implemented**: any change to `regenerator_
  effectiveness()`, `StateDependentLossModel`, or a "Hypereg material" —
  none of those would be supported by what this source actually says.
  Also not implemented: a quantitative frequency-vs-`n` optimum, since
  the paper gives no validated pressure-drop-reduction data to calibrate
  against (only one illustrative `n=4` example) — see
  `core/hypereg_analysis.py` for what IS shown (a qualitative, clearly
  labeled sensitivity sweep, not a claimed optimum).

## Status

Left partially open by design: the pumping-power mechanism is now
represented (closing the "is this a real, implementable idea in this
model's terms" question), but a genuine design recommendation (what
`n_parallel_subregenerators`, at what frequency, is actually worthwhile)
would need real pressure-drop or prototype data this early-stage, only
just-published concept does not yet have. This mirrors how several
existing ROADMAP items (e.g. items B5, B7) are left open rather than
resolved with invented precision.