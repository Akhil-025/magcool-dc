## What Magnetocaloric Cooling Can Honestly Claim for Data Centers

*Every number in this document is reproduced from this repo's own already-computed,
already-cited results (`results/comparison_table.csv`, `regime_crossover_analysis.py`,
`emissions.py`, `beverage_cooler_validation.py`, `commercial_landscape.py`) or from a
peer-reviewed/press source those modules already cite. No new modeling was done to
write this — it is a reframing of existing findings, not a new result.*

### 1. The COP case does not exist. Say so plainly.

Two independent searches of this repo's own validated AMR model, run against
this repo's own baseline vapor-compression (VCC) model, both come back null:

- **`run_cop_crossover_search()`** swept span 3–30 K and VCC second-law
  efficiency from 0.25 (worse than any real data-center compressor) to 0.55
  (well-optimized chilled water), plus a broad AMR design grid at every
  point. AMR never exceeds VCC — not even at VCC's own worst-case setting.
- **`run_material_and_field_crossover_search()`** repeated the search using
  this repo's own best-ranked giant-MCE material family
  (La(Fe,Si)₁₃Hᵧ), composition-tuned Curie-graded cascades, and fields up to
  7 T (superconducting-magnet territory, well past the 1.5–2 T range any
  compact data-center unit could realistically use). Still null, at every
  span tested:

  | Span | Best AMR COP found | VCC COP (η=0.42) |
  |---|---|---|
  | 5 K  | 9.70 | 24.36 |
  | 10 K | 7.05 | 11.97 |
  | 15 K | 5.43 | 7.84  |
  | 20 K | 4.32 | 5.77  |

  The interesting part isn't just that AMR loses — the search never even
  *wants* the high fields it was given access to. At every span the winning
  design sits at the lowest field tried (2 T), because in this repo's own
  calibrated loss model, field-scaling parasitic losses (eddy currents,
  hysteresis) grow at least as fast as the extra adiabatic-temperature-change
  headroom a bigger magnet buys. That's a structural property of the physics
  here, not a search that gave up too early.

- **The wider caloric-technology survey** (`alternative_caloric_comparison.py`)
  finds the same pattern outside magnetics too: elastocaloric, barocaloric,
  and electrocaloric literature COP claims are almost all simulations or
  projections, not measured devices. The one genuine exception — Li et al.
  2023's electrocaloric heat pump, a real measured 54%-of-Carnot device COP
  — only delivers 2.1 W of cooling power, trivial next to kilowatt-scale
  data-center needs, and has not been scaled up. **No caloric technology
  surveyed in this repo, measured or simulated, has a credible path to
  beating VCC's COP at data-center scale and span.** A paper or pitch that
  leads with a COP-superiority claim is making a claim this repo's own
  model — and the wider literature it draws from — does not support.

### 2. Refrigerant elimination does not rescue the comparison on emissions either

It's tempting to reach for "no refrigerant, so lower total emissions" as a
fallback argument once COP superiority is off the table. This repo checked
that directly, using the same COPs its own real-device comparison already
computed (AMR=1.76, VCC=6.66, at the Magnotherm Eclipse cabinet's 0.4 kW
operating point):

| Technology | Refrigerant (tCO2e/yr) | Operational (tCO2e/yr) | Total (tCO2e/yr) |
|---|---|---|---|
| Magnetic (AMR), no refrigerant | 0.000 | 0.990 | **0.990** |
| Vapor-compression (R-410A, GWP 2088) | 0.013 | 0.262 | **0.275** |
| Liquid cooling | 0.009 | 0.116 | **0.125** |

AMR's total emissions come out **3.6× higher** than vapor-compression's at
this COP gap. Operational (energy-driven) emissions dominate the
refrigerant-leak term by roughly two orders of magnitude here — exactly as
`emissions.py`'s own docstring warns before you even run the numbers. This
holds even using a high-GWP refrigerant (R-410A, GWP 2088); the gap would
only widen against a lower-GWP alternative like R-32 (GWP 675) or an HFO
(R-1234ze, GWP 7). **If AMR's COP trails VCC's by this much, eliminating
refrigerant does not make the deployment lower-carbon — it makes it
higher-carbon**, unless the grid supplying it is unusually clean or COP
parity is closed some other way.

### 3. What the real, deployed magnetocaloric systems actually show

The one segment where magnetocaloric cooling is commercially real, not
aspirational, is beverage refrigeration — and it is instructive precisely
because neither of its two independently-published data points claims COP
superiority:

- **Magnotherm Eclipse** (real deployment: 11-week in-store pilot at a
  German REWE supermarket, May–Sept 2025; 10–20 more units rolling out in
  2026): a 0.4 kW cabinet held 4–5°C using **15% less energy than the
  incumbent R290 (propane) unit at the same duty**. That's a genuine,
  real-world efficiency win — but it's a comparison against a small,
  already-low-GWP refrigerant system in a niche (small-cabinet
  refrigeration) where vapor-compression's own advantages are weakest, not
  a data-center-scale result, and Magnotherm has not published enough of
  the cabinet's internals (field, frequency, mass, flow rate) to check it
  against this repo's own physics model directly.
- **Polaris** (peer-reviewed, the first CE-certified commercial magnetic
  beverage cooler): reports, directly and quantitatively, a **plug-in
  (whole-device electrical) COP of 1.0** and a second-law efficiency of
  5.4% at 0.8 T / 15 K span — nowhere near vapor-compression parity. The
  paper's own headline finding is that real, market-available low-flow
  pumps (~20% efficient) impose such a large parasitic penalty that
  *plug-in* COP, not thermodynamic COP, has to drive the design — the same
  "parasitic losses matter as much as the underlying caloric physics" theme
  this repo's own loss model was built around.
- **Vendor claims at data-center scale** (Magnotherm Stellar, ~125 kW;
  Cooltech Applications, 10–15 kW class, claimed COP 5–6) exist, but are
  trade-press announcements with no public datasheet — not on the same
  evidentiary footing as Polaris, and not independently verified anywhere
  in this repo's corpus. Cooltech's own earlier device (`Cooltech_France_
  2016` in `data/amr_experimental_benchmarks.csv`) does not calibrate
  against this repo's model at any flow rate, an open, honestly-reported
  gap — a reason for caution before citing the newer claim at face value,
  not a reason to dismiss it outright.

None of these is a COP-superiority result. What they are is evidence that
magnetocaloric cooling is buildable, certifiable, and commercially real
*today*, in a segment where the efficiency bar (small refrigerant charge,
already-low-GWP incumbent) is achievable — which is a fundamentally
different, and more defensible, claim than "it will out-COP a chiller."

### 4. The actual, defensible pitch

Put together, this repo's own numbers support a specific, narrower claim —
and it's a real one, not a consolation prize:

> **Magnetocaloric cooling's case for data centers is refrigerant
> elimination and regulatory/environmental parity at matched COP — not COP
> superiority.**

What that means concretely, and what it doesn't:

- **It doesn't mean "lower emissions regardless of COP."** Section 2 shows
  that's false if the COP gap stays open — refrigerant elimination only
  helps once AMR's COP is close enough to VCC's that operational emissions
  aren't swamping the comparison. The honest sequencing is: close the COP
  gap first (or find an application where the COP bar is lower, per
  Section 3), *then* refrigerant elimination becomes a genuine net
  environmental win rather than a wash.
- **It does mean a real, live regulatory tailwind.** High-GWP refrigerants
  like R-410A (GWP 2088, still the incumbent this repo's own baseline model
  uses) are already being phased down under the Kigali Amendment to the
  Montreal Protocol and the EU F-Gas Regulation. A refrigerant-free
  technology that reaches COP *parity* — not superiority — with
  vapor-compression sidesteps that entire regulatory trajectory (leak
  reporting, phase-down quotas, eventual bans on specific blends) rather
  than needing to out-compress it. This repo's model doesn't need to show
  AMR winning on COP for that argument to hold — it needs to show AMR
  closing the gap enough that Section 2's operational-emissions penalty
  stops dominating, which is a materially lower bar than outright COP
  superiority.
- **It does mean leading with the segment where this is already true.**
  Section 3's beverage-refrigeration deployments are the existence proof:
  a real, small-scale, low-span, low-GWP-incumbent niche where AMR's COP
  penalty is small enough (or the baseline is already efficient enough)
  that "no refrigerant" is a genuine, not just theoretical, selling point.
  Data-center cooling, per Section 1, is not yet that niche in this
  repo's own model, at any span or material tried.

**What would change this document's conclusion**, stated plainly rather
than left implicit: either (a) a real, non-trivial-scale AMR device COP
that closes most of the gap in Section 1's table — Section 1 already
checked every material, field, and staging combination in this repo's own
model and found none that does it, so this would need either a physics
lever genuinely outside what was searched (see
`docs/Literature_Review.md`'s new "hysteresis-exploiting multicaloric
cycle" and "hysteresis-reducing dopant" entries for the two most credible
untested candidates, neither yet quantified enough to model), or (b) a
grid carbon-intensity or refrigerant-GWP scenario extreme enough to flip
Section 2's ratio even at the current COP gap — checked directly here to
be false at the R-410A/current-COP-gap baseline, not assumed.
