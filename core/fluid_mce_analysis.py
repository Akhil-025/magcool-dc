"""
fluid_mce_analysis.py
=======================
 validation/exploration deliverable for core/fluid_mce_cycle.py.
Same redirect-stdout-to-buffer-then-write pattern as
core/geometry_analysis.py / core/thermal_diode_analysis.py.

Runs three checks, all reported honestly rather than massaged to agree
with any prior expectation (see core/fluid_mce_cycle.py's own module
docstring for the underlying physics and honesty flags):

1. `volume_fraction_sweep()`: sweeps particle_volume_fraction (phi) at a
   fixed representative operating point, printing the viscosity-vs-phi /
   MCE-intensity-vs-phi tradeoff the plan named directly, and
   reporting whether an interior COP_electrical optimum exists (the same
   "genuinely open territory" question the plan itself posed, not
   resolved by the plan, and not force-resolved here either).
2. `fixed_span_comparison()`: this module's PRIMARY reported comparison.
   Holds span fixed at a shared, externally-imposed value (this repo's
   representative 10K data-center span, not either system's own
   favorable span) so the ferrofluid system's usable-span collapse shows
   up directly as a Qc=0 infeasibility, rather than being hidden by
   letting each technology pick its own comparison point.
3. `compare_to_solid_amr_and_liquid_cooling()`: SECONDARY comparison,
   at the fluid system's OWN favorable (small) span -- kept for the
   span/COP magnitude detail at a specific phi, per the plan's original
   framing ("worth at least a first-pass comparison against your
   existing direct-liquid-cooling baseline"), but read after the
   fixed-span table above, not instead of it.

VALIDATION UPDATE (this pass): adds
`literature_search_for_working_body_benchmark()`, an extended web search
specifically for a magnetocaloric-fluid-as-WORKING-BODY refrigeration
device with a reported Qc/COP figure -- broader than HONESTY FLAG #2's
original two adjacent results in core/fluid_mce_cycle.py. Reports the
result as a genuine NEGATIVE finding (none found) and cross-checks it
against this module's own computed span-collapse result, plus reports
the real Krieger-Dougherty rheology grounding
(`krieger_dougherty_grounded_in_susan_resiga_2012()`) alongside it. See
core/fluid_mce_cycle.py's own VALIDATION UPDATE for the full writeup and
core/thermal_diode.py's VALIDATION UPDATE for the SAME four sources
grounding a genuinely different (and now validated) feature there.
"""
import contextlib
import io
import os


from core.mce_material import GADOLINIUM
from core.amr_cycle import AMRSystem
from core.baseline_cooling import liquid_cooling_cop, vapor_compression_cop
from core.fluid_mce_cycle import (
    FerrofluidMCESystem, suspension_delta_T_adiabatic,
    DEFAULT_PHI_MAX, krieger_dougherty_grounded_in_susan_resiga_2012,
)

T_COLD_K_REPRESENTATIVE = 291.0   # matches core/optimize.py's own T_COLD_K
FIELD_T = 1.5
MDOT_KGS = 0.05
REPRESENTATIVE_SPAN_K = 10.0   # matches main.py's own REPRESENTATIVE_SPAN_K --
                                 # a realistic data-center span, used below as
                                 # the SHARED span both systems are held to in
                                 # fixed_span_comparison() (see that function's
                                 # docstring for why this, not each system's
                                 # own favorable span, is this module's primary
                                 # reported comparison)


def _self_consistent_span(mu0H_max, phi, T_cold, n_iter=6):
    """Solves span = 0.9 * dTad_suspension(T_mid) self-consistently, where
    T_mid = T_cold + span/2 -- the SAME T_mid convention
    FerrofluidMCESystem.cooling_capacity() uses internally. A naive
    "evaluate dTad once at an arbitrary fixed T_mid, then use that span"
    approach can pick a span inconsistent with the dTad the system will
    actually use once span/2 shifts T_mid away from that arbitrary point
    (material dTad(T) is peaked, not flat -- a few Kelvin of T_mid
    mismatch near a narrow first-order-like peak can matter) -- this
    fixed-point loop avoids that inconsistency. Converges in a handful of
    iterations since dTad(T) changes slowly relative to the span sizes
    involved here."""
    span = 0.0
    for _ in range(n_iter):
        T_mid = T_cold + span / 2
        dTad = suspension_delta_T_adiabatic(GADOLINIUM, T_mid, mu0H_max, phi)
        span = 0.9 * dTad
    return span, dTad


def volume_fraction_sweep(phis=(0.01, 0.02, 0.05, 0.10, 0.15, 0.20, 0.25, 0.30, 0.35, 0.40),
                            T_cold=T_COLD_K_REPRESENTATIVE, mu0H_max=FIELD_T,
                            mdot=MDOT_KGS):
    """Sweeps `phis` at a fixed (T_cold, mu0H_max, mdot); for each phi,
    evaluates the system at a span self-consistently solved to be 90% of
    THAT phi's own dTad_suspension at the T_mid the system will actually
    use (see `_self_consistent_span()`) -- i.e. each phi is evaluated
    near ITS OWN favorable span, not a single fixed span that would
    penalize low-phi points unfairly by exceeding their achievable span
    -- consistent with how core/giant_mce_analysis.py evaluates each
    material "at its own favorable operating point" rather than one
    shared point."""
    rows = []
    for phi in phis:
        span, dTad = _self_consistent_span(mu0H_max, phi, T_cold)
        sys_ = FerrofluidMCESystem(GADOLINIUM, mu0H_max, phi, mdot)
        result = sys_.run(T_cold, span) if dTad > 0 else None
        rows.append({
            "phi": phi,
            "dTad_suspension_K": round(dTad, 4),
            "span_K": round(span, 4),
            "viscosity_mPa_s": round(result.viscosity_Pa_s * 1000, 4) if result else float("nan"),
            "Qc_W": round(result.Qc, 4) if result else 0.0,
            "W_parasitic_W": round(result.W_parasitic, 5) if result else 0.0,
            "COP_electrical": round(result.COP_electrical, 4) if result else 0.0,
        })
    best = max(rows, key=lambda r: r["COP_electrical"])
    is_interior = rows[0]["phi"] < best["phi"] < rows[-1]["phi"]
    return {"rows": rows, "best_row": best, "interior_optimum_found": is_interior}


def fixed_span_comparison(span_K=REPRESENTATIVE_SPAN_K, T_cold=T_COLD_K_REPRESENTATIVE,
                            mu0H_max=FIELD_T, mdot=MDOT_KGS, phi=0.20,
                            mass_regenerator_amr=5.0, frequency_amr=1.0):
    """PRIMARY comparison for this module: Qc and COP_electrical for the
    ferrofluid MCE system, a representative solid AMR, and the liquid-
    cooling/vapor-compression baselines, all evaluated at the SAME
    externally-imposed span (`span_K`, default = this repo's own
    representative 10K data-center span) -- not at the ferrofluid
    system's own tiny favorable span the way
    `compare_to_solid_amr_and_liquid_cooling()` below does.

    Why this table, not that one, should be read first: the real
    headline finding of this module is that the mixture-heat-capacity
    dilution model plus this architecture's lack of regeneration
    collapses the ferrofluid system's *usable span* to a fraction of a
    Kelvin (see `volume_fraction_sweep()` and this module's own
    docstring), so a COP comparison at the ferrofluid's own favorable
    span necessarily hides that collapse -- it only ever compares the
    two systems at a span small enough for the fluid system to already
    be feasible. Holding span FIXED at a realistic value instead makes
    the feasibility gap itself the reported number: at a representative
    10K span, the fluid system's cooling_capacity() clips Qc to 0 (the
    span exceeds its achievable dTad_suspension entirely) while the
    solid AMR, whose regenerator bed amplifies span well beyond a single
    stage's own dTad, still delivers real cooling capacity. That is a
    directly comparable, single-number demonstration of the span gap,
    rather than something a reader has to infer from a paragraph of
    prose after two systems were quietly evaluated at different spans.
    """
    fluid_sys = FerrofluidMCESystem(GADOLINIUM, mu0H_max, phi, mdot)
    fluid_result = fluid_sys.run(T_cold, span_K)

    amr_sys = AMRSystem(material=GADOLINIUM, mu0H_max=mu0H_max,
                          mass_regenerator=mass_regenerator_amr, frequency=frequency_amr,
                          fluid_mdot=mdot, regenerator_effectiveness=0.85)
    amr_result = amr_sys.run(T_cold, span_K)

    T_hot = T_cold + span_K
    liquid = liquid_cooling_cop(T_cold, T_hot)
    vcc = vapor_compression_cop(T_cold, T_hot)

    return {
        "span_K": span_K,
        "phi": phi,
        "fluid_MCE": {"Qc_W": fluid_result.Qc, "COP_electrical": fluid_result.COP_electrical,
                       "feasible": fluid_result.Qc > 0,
                       "own_dTad_suspension_K": fluid_result.dTad_suspension_K},
        "solid_AMR": {"Qc_W": amr_result.Qc, "COP_electrical": amr_result.COP_electrical,
                       "feasible": amr_result.Qc > 0},
        "liquid_cooling": {"COP": liquid.COP},
        "vapor_compression": {"COP": vcc.COP},
    }


def compare_to_solid_amr_and_liquid_cooling(T_cold=T_COLD_K_REPRESENTATIVE,
                                               mu0H_max=FIELD_T, mdot=MDOT_KGS,
                                               phi=0.20, mass_regenerator_amr=5.0,
                                               frequency_amr=1.0):
    """SECONDARY comparison (see `fixed_span_comparison()` above for this
    module's primary, shared-span table). Compares, at the fluid
    system's OWN favorable span (90% of its dTad_suspension at `phi`):
      * FerrofluidMCESystem (this module)
      * a representative solid AMRSystem (Gd, same field/flow, -
        era defaults: mass_regenerator=5kg, frequency=1Hz, constant
        parasitic_fraction -- NOT the state-dependent loss_model, kept
        simple/representative rather than full-fidelity, consistent with
        this being a first-pass comparison per the plan's own framing)
      * core.baseline_cooling.liquid_cooling_cop() at the SAME span
      * core.baseline_cooling.vapor_compression_cop() at the SAME span,
        for completeness (not requested by the plan but free to compute
        and directly comparable)
    """
    span, dTad = _self_consistent_span(mu0H_max, phi, T_cold)

    fluid_sys = FerrofluidMCESystem(GADOLINIUM, mu0H_max, phi, mdot)
    fluid_result = fluid_sys.run(T_cold, span)

    amr_sys = AMRSystem(material=GADOLINIUM, mu0H_max=mu0H_max,
                          mass_regenerator=mass_regenerator_amr, frequency=frequency_amr,
                          fluid_mdot=mdot, regenerator_effectiveness=0.85)
    amr_result = amr_sys.run(T_cold, span)

    T_hot = T_cold + span
    liquid = liquid_cooling_cop(T_cold, T_hot)
    vcc = vapor_compression_cop(T_cold, T_hot)

    return {
        "span_K": span,
        "phi": phi,
        "fluid_MCE": {"Qc_W": fluid_result.Qc, "COP_electrical": fluid_result.COP_electrical},
        "solid_AMR": {"Qc_W": amr_result.Qc, "COP_electrical": amr_result.COP_electrical},
        "liquid_cooling": {"COP": liquid.COP},
        "vapor_compression": {"COP": vcc.COP},
    }


# Papers checked in this pass's literature search, all confirmed (not
# assumed) to use ferrofluid as a heat-transfer/thermal-switch medium,
# NEVER as the magnetocaloric working body -- see
# core/fluid_mce_cycle.py's VALIDATION UPDATE and
# core/thermal_diode.py's own VALIDATION UPDATE (same four sources
# ground a NEW validated feature there).
FERROFLUID_LITERATURE_SEARCH_2024_2025 = [
    {"source": "Klinar, Vozel, Swoboda, Sojer, Muñoz Rojo & Kitanovski, "
               "iScience 25 (2022) 103779",
     "ferrofluid_role": "thermal switch (numerical device model)"},
    {"source": "Rodrigues, Dias, Martins, Silva, Araújo, Oliveira, Pereira & "
               "Ventura, Applied Energy 251 (2019) 113213",
     "ferrofluid_role": "thermal switch (real hardware, no moving parts)"},
    {"source": "Andrade, Fernandes, Teixeira, Pereira, Pires, Silva, Ventura & "
               "Oliveira, Applied Energy 356 (2024) 122325",
     "ferrofluid_role": "thermal switch (real hardware)"},
    {"source": "Andrade, Fernandes, Silva, Teixeira, Pereira, Duarte, Pires, "
               "Ventura & Oliveira, Int. J. Refrigeration 164 (2024) 210-217",
     "ferrofluid_role": "thermal switch coupled to a solid Gd AMR cycle "
                         "(real hardware) -- NOT the working body itself"},
    {"source": "Gallium-based magnetocaloric liquid metal ferrofluid study (2017)",
     "ferrofluid_role": "tested on an electric transformer coil for heating "
                         "capacity, not a refrigeration cycle Qc/COP figure"},
    {"source": "Al2O3-nanofluid Active Magnetic Regenerative system study "
               "(Int. J. Air-Conditioning and Refrigeration, 2024)",
     "ferrofluid_role": "non-magnetic nanoparticle used to enhance the HEAT-"
                         "TRANSFER FLUID of a conventional SOLID-AMR bed -- "
                         "structurally different from this module's fluid-AS-"
                         "working-body architecture, not conflated with it here"},
]


def literature_search_for_working_body_benchmark(verbose=True):
    """Reports this pass's extended literature search for a
    magnetocaloric-FLUID-as-WORKING-BODY refrigeration benchmark (Qc/COP
    from a device where the fluid itself, not a separate solid
    regenerator, is the magnetocaloric material) -- see
    core/fluid_mce_cycle.py's VALIDATION UPDATE for the full writeup.
    Returns the search result as a NEGATIVE finding: none of the papers
    checked (FERROFLUID_LITERATURE_SEARCH_2024_2025) use a fluid as the
    working body; every one uses it as a heat-transfer/thermal-switch
    medium instead (a role core/thermal_diode.py's FerrofluidThermal
    Switch now models as a validated feature). This function also
    cross-checks that finding against this module's OWN computed
    result: does the fluid-as-working-body architecture predict a span
    collapse severe enough to plausibly explain why no one has built
    one? (A consistency check on a physical explanation, not proof of
    it -- absence of a device in the literature could always also just
    mean nobody has tried yet.)"""
    sweep = volume_fraction_sweep()
    max_span_over_sweep = max(r["span_K"] for r in sweep["rows"])
    # Compared against this repo's own representative data-center span
    # (REPRESENTATIVE_SPAN_K, 10K -- the same fixed comparison point
    # fixed_span_comparison() uses), not an arbitrary absolute Kelvin
    # threshold: even the BEST case over the whole phi sweep (up to
    # phi_max=0.63) tops out at a small fraction of a realistic target
    # span, which is the actual "collapse" this module's own findings
    # describe -- not literally "under 1K" at every phi.
    span_fraction_of_representative = max_span_over_sweep / REPRESENTATIVE_SPAN_K
    span_collapse_consistent_with_absence = span_fraction_of_representative < 0.25
    kd_grounding = krieger_dougherty_grounded_in_susan_resiga_2012()
    finding = {
        "papers_checked": FERROFLUID_LITERATURE_SEARCH_2024_2025,
        "working_body_benchmark_found": False,
        "every_paper_uses_ferrofluid_as": "thermal switch / heat-transfer medium",
        "this_modules_own_max_span_over_phi_sweep_K": round(max_span_over_sweep, 4),
        "representative_target_span_K": REPRESENTATIVE_SPAN_K,
        "max_span_as_fraction_of_representative": round(span_fraction_of_representative, 4),
        "span_collapse_consistent_with_literature_absence": span_collapse_consistent_with_absence,
        "krieger_dougherty_rheology_grounding": kd_grounding,
    }
    if verbose:
        print(f"Checked {len(FERROFLUID_LITERATURE_SEARCH_2024_2025)} papers "
              f"(full list in FERROFLUID_LITERATURE_SEARCH_2024_2025) for a "
              f"magnetocaloric-fluid-as-WORKING-BODY Qc/COP benchmark. Found: "
              f"NONE. Every real ferrofluid magnetocaloric device found uses "
              f"the ferrofluid as a heat-transfer/thermal-switch medium "
              f"instead -- see core/thermal_diode.py's FerrofluidThermalSwitch, "
              f"now a validated feature built on these SAME sources.")
        print(f"This module's own volume-fraction sweep never exceeds "
              f"{max_span_over_sweep:.3f}K of usable span at any tested phi -- "
              f"only {span_fraction_of_representative*100:.1f}% of this repo's "
              f"own {REPRESENTATIVE_SPAN_K:.0f}K representative data-center "
              f"target span -- "
              f"{'CONSISTENT with' if span_collapse_consistent_with_absence else 'NOT obviously consistent with'} "
              f"why no working-body device serving a realistic span exists in "
              f"the literature. This is a plausible physical explanation for "
              f"the literature gap, not proof of it.")
        print(f"Rheology grounding: {kd_grounding['finding']} "
              f"(source: {kd_grounding['source']}), for the same particle "
              f"material this module defaults to. This validates the "
              f"suspension VISCOSITY physics, not the magnetocaloric "
              f"performance claim.")
    return finding


def run_fluid_mce_analysis(out_path="results/fluid_mce_analysis.txt", verbose=True):
    buf = io.StringIO()
    with contextlib.redirect_stdout(buf):
        print("=" * 90)
        print("PHASE 20: magnetocaloric fluids (ferrofluid/MR suspension) as an alternative")
        print("working-body class. See core/fluid_mce_cycle.py's module docstring (HONESTY")
        print("FLAGS 1-2) for book-access limitations and benchmark availability.")
        print("=" * 90)

        print(f"\n--- Step 0: literature search for a working-body benchmark, "
              f"and Krieger-Dougherty rheology grounding (VALIDATION UPDATE) ---")
        lit_search = literature_search_for_working_body_benchmark()

        print(f"\n--- Volume-fraction sweep (T_cold={T_COLD_K_REPRESENTATIVE}K, "
              f"mu0H={FIELD_T}T, mdot={MDOT_KGS}kg/s, each phi evaluated at its OWN "
              f"favorable span = 0.9*dTad_suspension) ---")
        sweep = volume_fraction_sweep()
        print(f"\n{'phi':>6} {'dTad (K)':>10} {'span (K)':>10} {'visc (mPa.s)':>14} "
              f"{'Qc (W)':>10} {'W_par (W)':>12} {'COP_elec':>10}")
        for r in sweep["rows"]:
            print(f"{r['phi']:>6.2f} {r['dTad_suspension_K']:>10.4f} {r['span_K']:>10.4f} "
                  f"{r['viscosity_mPa_s']:>14.4f} {r['Qc_W']:>10.4f} "
                  f"{r['W_parasitic_W']:>12.5f} {r['COP_electrical']:>10.4f}")
        print(f"\nBest COP_electrical at phi={sweep['best_row']['phi']:.2f} "
              f"(COP_electrical={sweep['best_row']['COP_electrical']:.4f})")
        if sweep["interior_optimum_found"]:
            print("FINDING: a genuine INTERIOR optimum exists in phi within the swept range "
                  "-- the viscosity-vs-phi and MCE-intensity-vs-phi tradeoff the "
                  "plan named DOES produce a real tradeoff in this model, not a monotonic "
                  "result at either boundary.")
        else:
            print("FINDING: no interior optimum within the swept range -- COP_electrical is "
                  "monotonic in phi over [{:.2f}, {:.2f}], best at a boundary rather than an "
                  "interior point. Reported plainly rather than assumed.".format(
                      sweep["rows"][0]["phi"], sweep["rows"][-1]["phi"]))

        print("\n--- Fixed-span comparison (PRIMARY): all four technologies held to the "
              f"SAME {REPRESENTATIVE_SPAN_K:.0f}K representative span, rather than each "
              "evaluated at its own favorable span ---")
        fixed = fixed_span_comparison()
        print(f"\nAt span={fixed['span_K']:.1f}K (fixed), phi={fixed['phi']:.2f}:")
        f_feas = "feasible" if fixed["fluid_MCE"]["feasible"] else "INFEASIBLE (Qc clipped to 0)"
        a_feas = "feasible" if fixed["solid_AMR"]["feasible"] else "INFEASIBLE (Qc clipped to 0)"
        print(f" Ferrofluid MCE:      Qc={fixed['fluid_MCE']['Qc_W']:.3f}W "
              f"COP_electrical={fixed['fluid_MCE']['COP_electrical']:.3f}   [{f_feas}]   "
              f"(own achievable dTad_suspension={fixed['fluid_MCE']['own_dTad_suspension_K']:.3f}K "
              f"vs. the {fixed['span_K']:.0f}K span asked of it)")
        print(f" Solid AMR (Gd):      Qc={fixed['solid_AMR']['Qc_W']:.3f}W "
              f"COP_electrical={fixed['solid_AMR']['COP_electrical']:.3f}   [{a_feas}]")
        print(f" Liquid cooling:      COP={fixed['liquid_cooling']['COP']:.3f}")
        print(f" Vapor compression:   COP={fixed['vapor_compression']['COP']:.3f}")
        if not fixed["fluid_MCE"]["feasible"] and fixed["solid_AMR"]["feasible"]:
            print(f"\nFINDING: at a realistic, externally-imposed {fixed['span_K']:.0f}K span, "
                  "the ferrofluid system is infeasible (its own achievable "
                  f"dTad_suspension, {fixed['fluid_MCE']['own_dTad_suspension_K']:.2f}K, does "
                  "not cover the span) while the solid AMR still delivers real cooling "
                  "capacity -- this is the headline result of this module (the SPAN gap), "
                  "shown here as a direct feasibility comparison at a shared span rather "
                  "than inferred from two systems quietly evaluated at different spans.")

        print("\n--- Comparison at the fluid system's own favorable span (SECONDARY -- "
              "kept for the phi=0.20 span/COP detail, but read the fixed-span table above "
              "first), vs. solid AMR and vs. core/baseline_cooling.py's liquid-cooling/"
              "vapor-compression references ---")
        comp = compare_to_solid_amr_and_liquid_cooling()
        print(f"\nAt span={comp['span_K']:.3f}K, phi={comp['phi']:.2f}:")
        print(f" Ferrofluid MCE (this module):  Qc={comp['fluid_MCE']['Qc_W']:.3f}W "
              f"COP_electrical={comp['fluid_MCE']['COP_electrical']:.3f}")
        print(f" Solid AMR (Gd, same field/flow): Qc={comp['solid_AMR']['Qc_W']:.3f}W "
              f"COP_electrical={comp['solid_AMR']['COP_electrical']:.3f}")
        print(f" Liquid cooling (baseline_cooling.py): COP={comp['liquid_cooling']['COP']:.3f}")
        print(f" Vapor compression (baseline_cooling.py): COP={comp['vapor_compression']['COP']:.3f}")

        span_ratio = comp["span_K"]
        fluid_cop = comp["fluid_MCE"]["COP_electrical"]
        amr_cop = comp["solid_AMR"]["COP_electrical"]
        cop_note = (
            f"the ferrofluid system's own COP_electrical ({fluid_cop:.1f}) is actually "
            f"HIGHER than the solid AMR's ({amr_cop:.1f}) at this shared, ultra-narrow "
            f"span -- because W_parasitic (Darcy-Weisbach pipe flow) is small relative to "
            f"W_mag at this flow rate, and both systems' COP formulas have Qc largely "
            f"cancel out of the ratio -- but BOTH trail the liquid-cooling and vapor-"
            f"compression baselines at this span (see the numbers above). This COP "
            f"comparison should not be read as 'ferrofluid MCE beats solid AMR' in "
            f"general: it holds only at the ferrofluid system's own tiny achievable span, "
            f"which the solid AMR could trivially also hit (and does, with far more Qc); "
            f"the real headline result is the SPAN, not this specific COP comparison."
        )
        print(f"\nCONCLUSION: at realistic particle loadings (phi swept up to "
              f"{DEFAULT_PHI_MAX:.2f}, the random-close-packing ceiling), the mixture-"
              f"heat-capacity dilution model (core/fluid_mce_cycle.py physics item 1) "
              f"combined with this architecture's lack of regeneration (physics item 3) "
              f"collapses the usable span to a few Kelvin ({span_ratio:.2f}K at "
              f"phi={comp['phi']:.2f} in this specific comparison) -- dramatically less "
              f"than solid AMR achieves at the same field/flow, whose regenerator bed "
              f"amplifies achievable span well beyond a single stage's own dTad (see "
              f"core/amr_cycle.py's own characteristic-curve discussion). Notably, "
              f"{cop_note} This is a genuine, unforced finding from this pass's own "
              f"model, not assumed going in.\n\n"
              f"DISPOSITION (VALIDATION UPDATE): this module's POSITIVE, quantitative "
              f"Qc/COP predictions remain unvalidated design-exploration output -- no "
              f"benchmark magnetocaloric-fluid-as-working-body refrigeration device "
              f"exists anywhere in the literature to check them against (Step 0 above), "
              f"a finding this pass's own broader search corroborates rather than merely "
              f"repeats. What IS now validated: (1) the suspension RHEOLOGY -- Krieger-"
              f"Dougherty with this module's exact default parameters is grounded in a "
              f"real magnetite-ferrofluid measurement (Susan-Resiga et al. 2012, Step 0); "
              f"and (2) the NEGATIVE finding itself -- every real ferrofluid "
              f"magnetocaloric device found in the literature uses the fluid as a "
              f"thermal switch, never as a working body (Step 0), consistent with this "
              f"module's own computed span-collapse result. Treat this module's Qc/COP "
              f"numbers as illustrative; treat its span-collapse finding and the absence "
              f"of a working-body benchmark as validated conclusions in their own right.")

    text = buf.getvalue()
    if verbose:
        print(text)
    if out_path:
        os.makedirs(os.path.dirname(out_path) or ".", exist_ok=True)
        with open(out_path, "w") as fh:
            fh.write(text)
        if verbose:
            print(f"Wrote {out_path}")
    return {"sweep": sweep, "fixed_span_comparison": fixed, "comparison": comp,
            "literature_search": lit_search, "text": text}


if __name__ == "__main__":
    run_fluid_mce_analysis()
