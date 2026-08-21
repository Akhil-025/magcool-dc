"""
commercial_landscape.py
========================
Phase 24 addition.

Magnetocaloric cooling for data centers stopped being a purely theoretical
question during this project's own lifetime: Magnotherm (DE) announced
"Stellar", a ~125 kW refrigerant-free AMR system explicitly marketed at
data-center cooling, and Cooltech Applications (FR) -- already one of this
repo's own experimental benchmark devices in
`data/amr_experimental_benchmarks.csv` -- has been reported at 10-15 kW
class output with claimed COP ~5-6 under data-center-like conditions.

This module does three things, each deliberately modest in scope:

1. Records these public claims as a structured, source-flagged dataset
   (`COMMERCIAL_SYSTEMS`) -- NOT re-validated against a datasheet, since
   neither vendor publishes one; every number here is a *press-release or
   trade-press claim*, flagged as such, not a measured benchmark point on
   the same footing as `data/amr_experimental_benchmarks.csv`.
2. Explicitly disambiguates "magnetocaloric cooling" from "magnetic-bearing
   chiller" (Johnson Controls YDAM, Munters Circlemiser, etc.) -- these use
   frictionless magnetic-bearing compressors but are still ordinary
   vapor-compression thermodynamically, and get conflated with
   magnetocaloric systems in casual search results.
3. Confronts this repo's own model against the one commercial claim it CAN
   compare against on a like-for-like basis (Cooltech), and reports the gap
   honestly rather than leaving it as a footnote -- this repo's own
   `data/amr_experimental_benchmarks.csv` already flags "Cooltech_France_
   2016" as not calibrating at any flow rate under the CORE loss model, so
   a reviewer aware of Cooltech's data-center-specific claims would ask
   exactly this question.

Sources (trade press / vendor announcements, 2024-2026 -- NOT peer-reviewed;
treat every number here as a claim to be engaged with, not a validated
fact):
  - Magnotherm, "Stellar" product announcement and trade-press coverage
    (2025-2026), magnotherm.com and industry press (e.g. Data Center
    Dynamics, Cooling Post).
  - Cooltech Applications trade-press coverage of data-center-oriented
    demonstrations (2024-2026), various industry outlets; underlying device
    physics already partially represented in this repo via
    `Cooltech_France_2016` in `data/amr_experimental_benchmarks.csv`
    (Vasile & Muller, 2006 lineage device family).
  - Johnson Controls YDAM / York magnetic-bearing centrifugal chiller and
    Munters Circlemiser: vendor literature; both are vapor-compression /
    desiccant-wheel systems respectively, NOT magnetocaloric, included here
    only to document the naming collision.
"""

from dataclasses import dataclass
from typing import Optional

from core.amr_cycle import AMRSystem
from core.mce_material import GADOLINIUM
from core.loss_model import StateDependentLossModel


@dataclass
class CommercialClaim:
    name: str
    technology_class: str            # "magnetocaloric-AMR" or "not-magnetocaloric (name collision)"
    claimed_capacity_kW: Optional[float]
    claimed_COP: Optional[float]
    claimed_COP_range: Optional[str]
    target_application: str
    status_2026: str
    source_note: str
    is_magnetocaloric: bool


COMMERCIAL_SYSTEMS = [
    CommercialClaim(
        name="Magnotherm Stellar",
        technology_class="magnetocaloric-AMR",
        claimed_capacity_kW=125.0,
        claimed_COP=None,
        claimed_COP_range=None,
        target_application="data center / industrial cooling, refrigerant-free",
        status_2026="announced 2025, positioned for 2025-2026 deployment; "
                     "no independently-audited datasheet found in this "
                     "project's corpus -- vendor/trade-press claim only",
        source_note="Magnotherm product announcement + trade press "
                     "(Data Center Dynamics / Cooling Post-class coverage), "
                     "2025-2026. First commercial AMR system explicitly "
                     "marketed at the data-center use case this repo models.",
        is_magnetocaloric=True,
    ),
    CommercialClaim(
        name="Cooltech Applications (data-center-oriented unit)",
        technology_class="magnetocaloric-AMR",
        claimed_capacity_kW=12.5,   # midpoint of reported 10-15 kW range
        claimed_COP=5.5,            # midpoint of reported 5-6 range
        claimed_COP_range="5-6",
        target_application="data center cooling (claimed)",
        status_2026="reported demonstrations 2024-2026; same vendor lineage "
                     "as this repo's own 'Cooltech_France_2016' experimental "
                     "benchmark device (data/amr_experimental_benchmarks.csv), "
                     "but this specific data-center-oriented unit is NOT the "
                     "same physical device already in this repo's calibration "
                     "set -- treat capacity/COP numbers as an independent, "
                     "unverified claim, not a re-measurement of the existing "
                     "benchmark row.",
        source_note="Trade-press coverage, 2024-2026, vendor claim, no public "
                     "datasheet located in this project's corpus.",
        is_magnetocaloric=True,
    ),
    CommercialClaim(
        name="Johnson Controls YDAM (York magnetic-bearing chiller)",
        technology_class="not-magnetocaloric (name collision)",
        claimed_capacity_kW=None,
        claimed_COP=None,
        claimed_COP_range=None,
        target_application="data center / commercial HVAC chilled water",
        status_2026="commercially mature, widely deployed",
        source_note="Vendor literature. Uses frictionless magnetic-bearing "
                     "centrifugal compressors -- ordinary vapor-compression "
                     "refrigerant cycle, magnetic bearings only reduce "
                     "mechanical friction loss. Included ONLY to document "
                     "that 'magnetic chiller' search results are dominated "
                     "by this unrelated technology, not magnetocaloric AMR.",
        is_magnetocaloric=False,
    ),
    CommercialClaim(
        name="Munters Circlemiser",
        technology_class="not-magnetocaloric (name collision)",
        claimed_capacity_kW=None,
        claimed_COP=None,
        claimed_COP_range=None,
        target_application="desiccant dehumidification / data center air handling",
        status_2026="commercially mature",
        source_note="Vendor literature. Desiccant-wheel dehumidification "
                     "technology; does not use magnetic fields or the "
                     "magnetocaloric effect at all. Included only as a "
                     "second documented naming-collision case.",
        is_magnetocaloric=False,
    ),
]


def disambiguation_note():
    """One-paragraph, citable disambiguation for a paper's introduction."""
    return (
        "IMPORTANT TERMINOLOGY NOTE: 'magnetic cooling' in commercial and "
        "search-engine contexts frequently refers to magnetic-bearing "
        "vapor-compression chillers (e.g. Johnson Controls YDAM, York "
        "magnetic-bearing centrifugal chillers) or desiccant systems "
        "(e.g. Munters Circlemiser) -- none of which use the magnetocaloric "
        "effect. These are conventional refrigerant-cycle or desiccant "
        "systems with magnetic bearings reducing mechanical friction only. "
        "This work concerns magnetocaloric cooling specifically: solid-state, "
        "refrigerant-free heat pumping via the field-induced entropy change "
        "of a magnetic working material (active magnetic regenerator, AMR). "
        "Commercial systems in this narrower, correct sense include "
        "Magnotherm's 'Stellar' and Cooltech Applications' data-center-"
        "oriented units (both engaged with quantitatively below)."
    )


def model_prediction_at_cooltech_point(mu0H_max=1.5, frequency=2.0,
                                        mass_regenerator=8.0,
                                        fluid_mdot=0.10,
                                        T_cold_K=291.15, span_K=10.0,
                                        regenerator_effectiveness=0.85,
                                        verbose=True):
    """Run this repo's own AMRSystem at a scaled-up, Cooltech-class
    operating point (higher mass/flow than the small 2016 lab prototype,
    reasoned to approach a ~10 kW-class unit) and report COP_electrical
    against Cooltech's own claimed COP=5-6, quantifying the gap honestly.

    This uses the CORE (3-point) loss-model calibration -- the SAME
    calibration `data/amr_experimental_benchmarks.csv` already flags as
    NOT calibrating against Cooltech_France_2016 at any tested flow rate.
    So a mismatch here is an *expected*, already-documented model
    limitation being surfaced explicitly for the paper, not a new finding.
    """
    sys_ = AMRSystem(
        material=GADOLINIUM, mu0H_max=mu0H_max, mass_regenerator=mass_regenerator,
        frequency=frequency, fluid_cp=4186.0, fluid_mdot=fluid_mdot,
        regenerator_effectiveness=regenerator_effectiveness,
        loss_model=StateDependentLossModel(), use_ntu_thermal_model=True,
    )
    res = sys_.run(T_cold_K, span_K)  # run(T_cold, T_span) -- span, not T_hot
    cooltech = next(c for c in COMMERCIAL_SYSTEMS if c.name.startswith("Cooltech"))
    gap_pct = (100 * (res.COP_electrical - cooltech.claimed_COP) / cooltech.claimed_COP
               if cooltech.claimed_COP else float("nan"))
    if verbose:
        print(f"Model prediction at Cooltech-class point "
              f"(span={span_K}K, mu0H={mu0H_max}T, f={frequency}Hz, "
              f"mdot={fluid_mdot}kg/s, mass={mass_regenerator}kg):")
        print(f"  COP_electrical (this model, CORE calibration) = {res.COP_electrical:.2f}")
        print(f"  Qc = {res.Qc:.1f} W")
        print(f"  Cooltech claimed COP (trade press, unverified) = "
              f"{cooltech.claimed_COP_range}")
        print(f"  Gap vs. claimed midpoint: {gap_pct:+.1f}%")
        print("  HONEST FRAMING FOR THE PAPER: this repo's own "
              "amr_experimental_benchmarks.csv already documents that the "
              "CORE loss-model calibration does not reproduce "
              "Cooltech_France_2016 at any tested flow rate. A large gap "
              "here is a continuation of that already-flagged limitation, "
              "not a new discrepancy -- and it should be stated as such, "
              "not silently omitted. Plausible explanations to discuss: "
              "(a) vendor claims may reflect best-case/marketing operating "
              "points rather than continuous duty; (b) the CORE loss model "
              "is calibrated on 2010s lab prototypes and may not extrapolate "
              "to a more recent, differently-engineered commercial unit; "
              "(c) this repo's flow/mass scaling to '10 kW class' is a "
              "reasoned estimate, not a datasheet-matched configuration.")
    return {
        "model_COP_electrical": res.COP_electrical,
        "model_Qc_W": res.Qc,
        "cooltech_claimed_COP_range": cooltech.claimed_COP_range,
        "gap_pct_vs_claimed_midpoint": gap_pct,
    }


def write_commercial_landscape_report(path="results/commercial_landscape.txt"):
    lines = ["Commercial / state-of-the-art magnetocaloric cooling landscape",
             "=" * 70, "", disambiguation_note(), "",
             "Magnetocaloric-AMR systems targeting this repo's use case:", "-" * 60]
    for c in COMMERCIAL_SYSTEMS:
        if not c.is_magnetocaloric:
            continue
        lines.append(f"* {c.name}")
        lines.append(f"    capacity: {c.claimed_capacity_kW} kW  "
                      f"COP: {c.claimed_COP_range or c.claimed_COP}")
        lines.append(f"    status: {c.status_2026}")
        lines.append(f"    source: {c.source_note}")
        lines.append("")
    lines.append("Naming-collision systems (NOT magnetocaloric, documented for disambiguation):")
    lines.append("-" * 60)
    for c in COMMERCIAL_SYSTEMS:
        if c.is_magnetocaloric:
            continue
        lines.append(f"* {c.name} -- {c.source_note}")
    lines.append("")
    gap = model_prediction_at_cooltech_point(verbose=False)
    lines.append("This model vs. Cooltech's own data-center-oriented claim:")
    lines.append("-" * 60)
    lines.append(f"  model COP_electrical = {gap['model_COP_electrical']:.2f}")
    lines.append(f"  Cooltech claimed COP range = {gap['cooltech_claimed_COP_range']}")
    lines.append(f"  gap vs. claimed midpoint = {gap['gap_pct_vs_claimed_midpoint']:+.1f}%")
    lines.append("  (see docstring / run this module directly for full honest framing)")
    with open(path, "w") as f:
        f.write("\n".join(lines) + "\n")
    print(f"Wrote {path}")


if __name__ == "__main__":
    print(disambiguation_note())
    print()
    model_prediction_at_cooltech_point()
    print()
    write_commercial_landscape_report()
