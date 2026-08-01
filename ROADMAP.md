# Roadmap to madcool-suite parity

## Phase 1 — done
- [x] Mean-field MCE material model (`mce_material.py`), validated against
      Dan'kov et al. (1998) Gd data (`validation.py`)
- [x] 0-D AMR cycle model (`amr_cycle.py`)
- [x] Vapor-compression & liquid-cooling baseline COP models (`baseline_cooling.py`)
- [x] Order-of-magnitude economics comparison (`economics.py`)
- [x] Comparison driver across the ASHRAE thermal envelope (`main.py`)
- [x] Literature review with citations (`LITERATURE_REVIEW.md`)

## Phase 2 — system-level validation & sensitivity — done
- [x] Digitized published AMR prototype data (`data/amr_experimental_benchmarks.csv`:
      Astronautics/Jacobs 2014, DTU rotary Gd, Tušek single-bed Gd) as
      system-level COP/span validation targets
- [x] Calibrate-then-validate methodology (`validation_system.py`) —
      revealed that published COP figures are *electrical* (pump + motor
      overhead included), not thermodynamic-cycle-only; added a calibrated
      `parasitic_fraction` to `amr_cycle.py` and got the two comparable
      lab devices to ~10% agreement
- [x] Sobol/variance-based sensitivity analysis (`sensitivity.py` →
      `results/sobol_results.txt`) — found electrical COP is ~99.9%
      sensitive to `parasitic_fraction` alone in the current model
      structure, a genuine finding that motivates Phase 3
- [x] Response-surface (RSM) surrogate for cooling capacity Qc
      (`rsm.py` → `results/rsm_coefficients.txt`, R²=0.86 held-out)

## Phase 3 — state-dependent losses, optimization & multi-stage design — done
- [x] State-dependent loss model (`loss_model.py`): eddy ~ f²H², pumping ~ ṁ²,
      base overhead ~ Qc, calibrated against the 3 Phase 2 benchmark devices
      (exactly-determined 3-point fit — flagged as needing more data)
- [x] Re-ran Sobol with the new model: frequency (ST=0.68), flow (ST=0.22)
      and field (ST=0.09) now carry real sensitivity, resolving the Phase 2
      diagnostic finding
- [x] NSGA-III multi-objective optimization (`optimize.py`, COP vs. Qc vs.
      cost) → `results/pareto_front.csv` — exposed that `mass_regenerator`
      has no effect on Qc in the current model, so every Pareto design sits
      at the mass floor (a real finding, not hidden)
- [x] Multi-stage cascade design (`cascade.py`) → `results/cascade_comparison.csv`
      — confirms staging recovers the span range lost above 16 K, but ALSO
      shows both single- and multi-stage Gd/2T AMR trail vapor-compression
      and liquid cooling on electrical COP across the whole ASHRAE range at
      this design point. This is the project's current honest conclusion.

## Phase 4 — close the model gaps Phase 3 exposed — done
- [x] **NTU-based `thermal.py`**: packed-sphere-bed regenerator effectiveness
      from Wakao & Kaguei (1982) Nusselt correlation + utilization-factor
      degradation, wired into `amr_cycle.py` via `use_ntu_thermal_model`.
      Confirmed the fix: NSGA-III Pareto front mass values now spread
      1.0-14.5 kg instead of pinning at the floor
- [x] **Attempted the higher-field/giant-MCE re-run** — and found the real
      blocker: `mce_material.py`'s mean-field/Brillouin framework can't
      capture Gd5Si2Ge2's first-order-transition giant MCE (underpredicts
      ΔT_ad by ~10x), so `cascade.py` correctly returns zero capacity for it
      at the DC operating point rather than a misleadingly small nonzero
      number. Flagged explicitly in `mce_material.py`. Gd's own T_c=294K
      (21°C) already sits inside the ASHRAE 18-27°C range, which is worth
      noting in the paper as a reason Gd, not an exotic alloy, is the
      practical near-term candidate.

## Phase 5 — first-order-transition materials, economics/emissions — done
- [x] **Bean-Rodbell-family model**: `first_order_mce.py` implements an
      extended (6th-order) Landau free-energy model, calibrated to Gd5Si2Ge2's
      literature peak entropy change (~-18 J/(kg K) at 5T) — resolves the
      Phase 4 blocker with a physically appropriate framework, not a patch
- [x] **Giant-MCE vs. Gd formal comparison** (`giant_mce_analysis.py` →
      `results/giant_mce_analysis.txt`): confirmed the Curie-matching
      principle cleanly — Gd5Si2Ge2's peak-effect window (286.4K) sits ~5K
      below the ASHRAE range and collapses to zero there, but performs
      strongly (Qc=5319W, COP_elec=7.76) at its own matched point. Does NOT
      overturn the Phase 1-4 conclusion; identifies the concrete next step
      (composition-tuned Gd5(SixGe1-x)4, literature-documented as tunable)
- [x] **Refrigerant-free GWP/emissions comparison** (`emissions.py`):
      quantified honestly — real categorical benefit (zero leak/phase-out
      risk) but does not overturn the emissions comparison until AMR's COP
      gap closes; report both numbers plainly in the paper

## Phase 6 — extended validation stress-test + grounded economics — done
- [x] **Expanded `data/amr_experimental_benchmarks.csv`**: found a fully
      usable 4th device (Okamura & Hirano 2013, Qc=200W/COP=2.5/5K span/
      1.1T/1kg Gd) via targeted search. Added to `loss_model.py` as an
      EXTENDED calibration set.
- [x] **Stress-tested the loss model — and it failed the stress test.**
      The 4-point least-squares fit gives negative (unphysical) k_pump and
      base_frac, and leave-one-out CV shows errors up to +1639% — four
      orders of magnitude of device scale (6.5W to 2502W) apparently can't
      be pooled into one linear loss model. Reported via
      `run_extended_diagnostic()`, NOT silently absorbed into the
      production default, which stays on the stable 3-point CORE fit.
- [x] Attempted to add the Risoe/DTU 2011 (30K span) device too — could not
      calibrate; consistent with, and further evidence for, the Phase 1
      single-stage span-collapse finding.
- [x] **Grounded `economics.py`** in an actual costing study: Bjørk, Bahl &
      Smith (Int. J. Refrig. 34 (2011) 1805-1816), $40/kg magnet + $20/kg
      MCM, replacing the earlier loosely-sourced placeholder. Wired into
      `optimize.py`'s cost objective via `economics.material_cost()`.

## Phase 7 — remaining open items before the paper is fully evidence-complete
- [x] **Loss-model solver fix + scale-term hypothesis tested (not confirmed)**:
      switched `loss_model.py`'s calibration from unconstrained least
      squares (with negative coefficients clipped to zero post-hoc) to
      non-negative least squares (NNLS, Lawson & Hanson 1974,
      `scipy.optimize.nnls`), which solves the physically-constrained
      problem directly. This removes Phase 6's negative k_eddy/base_frac
      on the EXTENDED 4-point set and improves the worst leave-one-out
      error from +1639% to +682% — real progress, but still an
      order-of-magnitude miss, so the CORE 3-point fit stays the
      production default. Then tested the specific hypothesis Phase 6
      proposed — a device-*size* scale term — via
      `analyze_parasitic_fraction_scaling()`: sorting the four devices by
      Qc shows the parasitic fraction is **not monotonic in scale** (the
      smallest device, Tušek 6.5W, has one of the *lowest* fractions at
      11.8%; the largest, Astronautics 2502W, has the *highest* at 45.3%,
      which its own source paper attributes to that specific device's
      "mediocre" electrical-component efficiency, not a generic size
      effect). **Correction to the roadmap**: a size/scale term is not
      supported by the data in hand. What the EXTENDED set's instability
      more plausibly reflects is heterogeneous drivetrain
      topology/component-efficiency class across devices — the concrete
      next step is more benchmark devices with independently reported
      component efficiencies, not a size-dependent loss term.
- [x] **Curie-graded cascade with independent literature validation.**
      Added `core/giguere_validation.py`: read Giguere et al. (Phys. Rev.
      Lett. 83, 2262 (1999)) directly from the PDF already in this repo
      and cross-checked `first_order_mce.py`'s model against their DIRECT
      DeltaT_ad measurement (10.0 K at 7T for Gd5Si2Ge2, Clausius-Clapeyron
      cross-check 9.9K) rather than the Maxwell-relation "indirect" value
      the model was calibrated to. Result: the model overestimates the
      DIRECT measurement by ~2.4x — worse than the ~1.49x gap the paper
      itself found between its own indirect and direct methods, consistent
      with (and additive to) honesty flag #1's lattice-only-C_p concern.
      The model was NOT refit to close this gap (a 0-D lattice-only-C_p
      framework can't match both the peak DeltaS_M calibration and the
      direct DeltaT_ad simultaneously — that needs the transition's latent
      heat, out of scope). Instead, added an empirical correction factor
      (`dTad_correction`, ~0.41) derived from this one cross-check point,
      applied by default whenever a composition-tuned material is built.
      Added `composition_tuned_material(Tc_target)`, restricted to the
      LITERATURE-DOCUMENTED giant-MCE composition window (~20-290K,
      Pecharsky & Gschneidner APL 70, 3299 (1997); note Gd5Si4 at x=1,
      Tc=335K, is explicitly NOT a giant/first-order composition, so the
      window cannot be extended toward it) — requesting a Tc outside that
      range raises rather than silently extrapolating. Implemented
      `run_graded_cascade()`/`compare_graded_cascade()` in `cascade.py`:
      each stage's composition is matched to its own local operating
      temperature via an iterative peak search (needed because the
      model's peak DeltaT_ad sits off the nominal Tc, and that offset
      turned out NOT to be translation-invariant as Tc shifts — a real
      bug caught and fixed while building this, not assumed away).
      **Result at the ASHRAE range**: for 35/64 span x stage-count
      combinations swept (5-20K span, 1-4 stages), every stage's needed
      composition stays within the documented 20-290K window; the other
      29 have one or more stages exceed 290K and fall back to plain Gd for
      that stage only (larger spans/more stages push the hottest stage
      above the ceiling). At 10K span / 3 stages, the graded cascade
      delivers Qc=2381W vs. plain Gd's 1258W at essentially the same COP
      (3.83 vs 3.75) — same "bigger MCE buys Qc, not COP" pattern
      `giant_mce_analysis.py` already found for the fixed-composition
      case, now confirmed for the graded case in the application-relevant
      range. La(Fe,Si)13Hy was NOT added (no composition-tunability or
      independent-validation dataset for it was located in this repo or
      review) — that half of the original item remains open. Also flagged
      (not hidden): this idealized 6th-order Landau fit's transition is
      numerically much narrower (sub-1K) than the real, hysteresis/
      inhomogeneity-broadened transition Giguere et al.'s Fig. 3 shows
      (~10-15K wide) — a small number of individual span/stage-count
      cells show a stage's Qc collapsing to ~0 from residual
      peak-alignment error despite a nominally in-range composition. See
      `results/giguere_validation.txt` and
      `results/graded_cascade_comparison.csv`.
- [~] **Curve-level validation — partial, honest step, not the full item.**
      *Correction (see below):* an earlier pass reported Tušek 2010 and
      Nielsen 2011 as unavailable in this repository — that was wrong; both
      papers are present (`Papers/AMR Theory and Modeling/Development of a
      rotary magnetic refrigerator.pdf`; `Papers/AMR systems and
      prototypes/Review on numerical modeling...pdf`). What remains true is
      that neither contains the multi-point experimental characteristic
      curves this item needs (Tušek 2010 is a device-construction paper;
      Nielsen 2011 is a numerical-modeling review), so no digitization of
      characteristic curves was fabricated from them; that part of this
      item is still genuinely open. What *was* done: 3 of the 5 benchmark
      devices (Astronautics, DTU,
      Risø/DTU) turn out to already have a second (span, Qc) data point
      for the same physical device in `amr_experimental_benchmarks.csv`
      (a zero-span max-capacity or max-span zero-capacity reading). Added
      a `device_group` column to link these, and
      `run_curve_validation()` in `validation_system.py` now calibrates
      mdot at the normal operating point as before and predicts Qc at the
      companion span with that same fit — a real, independent check of
      predicted curve *shape* (companion point unused in calibration).
      Astronautics: +17.6% at zero span. DTU: predicts 0.0W at its
      reported no-load span (exact match). Risø/DTU: still fails to
      calibrate at all (consistent with the existing finding). Tušek and
      Okamura remain single-point. This is 2-point shape validation, not
      the full multi-point published curve — that part stays open below.
- [~] **Two more papers located in `Papers/` and worked through — a real
      device added, a real curve-source found but NOT digitized.**
      Searched the `Papers/AMR Theory and Modeling/` folder directly
      rather than assuming its contents from citations already in
      `LITERATURE_REVIEW.md`, and found two papers relevant to this open
      item that were not yet in `amr_experimental_benchmarks.csv`:
      - **Lozano, Capovilla, Trevizoli, Engelbrecht, Bahl, Barbosa Jr.
        (2016), "Development of a novel rotary magnetic refrigerator"**
        (POLO/UFSC), Int. J. Refrig., doi:10.1016/j.ijrefrig.2016.04.005.
        Its Table 3 gives 8 real (frequency, flow rate, span, Qc, COP)
        operating points as clean numbers — no digitization needed. Added
        all 8 plus the abstract's zero-span/no-load-span endpoints as a
        new `Lozano_POLO_UFSC_2016` device_group (mu0H=0.88T from the
        paper's measured H_high-H_low, mass=1.7kg Gd spheres, Sec.2/
        Table1). While adding these rows, found and fixed two pre-existing
        CSV bugs that had gone unnoticed because the file happened to end
        right after them: an unterminated quote on the Okamura row (was
        silently swallowing whatever row came after it into one giant
        field) and a blank trailing line that crashed `float()` once rows
        followed it. Extended `run_curve_validation()` with a guard so
        device_groups where every row independently varies its own
        frequency (a real multi-point sweep, not a fixed-condition span
        sweep) are routed to `run_system_validation()` instead of the
        2-point anchor/companion logic, which would otherwise silently
        compare across mismatched operating conditions. **Result**: only
        4 of the 8 rows (r4, r6, r7, r8) calibrate within mdot in
        [1e-6,5] kg/s; the other 4 report "no calibration found", now
        printed rather than silently dropped (fixed a real bug in
        `calibrate_and_check()` where the failure branch returned without
        printing in verbose mode). Of the 4 that do calibrate, COP errors
        are far larger than any other benchmark device (+590% to +805%)
        — Lozano's own paper reports COP=0.37-0.83, calling the results
        "modest ... in comparison with established cooling technologies",
        because its motor power (87-145W) is comparable to or exceeds its
        own Qc (61-120W). Fed these 4 points into `loss_model.py` as a
        new `CALIBRATION_POINTS_FURTHER_EXTENDED` set
        (`run_further_extended_diagnostic()`): the naive worst leave-one-
        out error actually *improves* (+167% vs. EXTENDED's +682%), but
        this is not genuine progress — it happens because Tušek's small
        W_parasitic becomes easier to hit with more points anchoring the
        fit. The more informative result is that all 4 Lozano points,
        held out individually, are underpredicted by a consistent ~89-91%
        (tight cluster, not noise) — real evidence that this device sits
        in a distinct, worse motor/inverter efficiency *class* a
        state-variable-only loss model (frequency, field, mdot, Qc) can't
        represent, strengthening rather than overturning the existing
        Phase 7 conclusion. CORE 3-point fit remains the production
        default. New tests in `tests/test_validation_system.py` and
        `tests/test_loss_model.py` cover all of this; full run
        saved to `results/loss_model_diagnostics.txt`.
      - **Tušek, Kitanovski, Zupan, Prebil, Poredoš (2013), "A
        comprehensive experimental analysis of gadolinium active magnetic
        regenerators"**, Appl. Therm. Eng. 53, 57-66. This turns out to be
        the paper the existing `Tusek_singlebed_Gd_2010` row's citation
        ("Tušek et al., Appl. Therm. Eng. 2011 dataset, 2010 device") was
        actually pointing at — but its own Sec. 3.1/Table 1 reports a
        DIFFERENT field (1.15T, not 1.69T) and none of its 6 AMR masses
        (0.093-0.176kg) matches 0.196kg exactly, so it's additional
        independent data, not a source/confirmation of that existing row;
        the true source of the 1.69T/0.196kg point is still unidentified.
        Its Figs. 10-11 (Qc-vs-span and COP-vs-span, 3 AMR geometries x 3
        flow ratios, 9 lines) ARE genuine multi-point published
        characteristic curves in the Tušek-2010/Nielsen-2011 sense this
        item originally asked for — **but digitizing them requires
        pixel-calibrated marker extraction from the figure images, which
        was NOT done here** (a rough, explicitly-non-authoritative visual
        read is in `results/tusek_ate2013_figs_notes.md` for whoever picks
        this up next; fabricating a precise-looking table from an
        uncalibrated eyeball read would be worse than leaving it open).
        So: this half of the item is now unblocked on *obtaining* the
        source (the paper is in `Papers/`) but still open on actually
        *digitizing* it.
      *Correction:* Nielsen et al. (2011) IS present in this repository
      (`Papers/AMR systems and prototypes/Review on numerical modeling of
      active magnetic regenerators for room temperature applications.pdf`
      — confirmed from its own title page: K.K. Nielsen et al., Int. J.
      Refrig. 34 (2011) 603-616). An earlier pass reported it as missing;
      that was a search error (only page-0 text was checked, and the
      author line falls past the truncation point used), not a genuine gap
      in `Papers/`. Having now read it: this paper is a review of numerical
      *modeling methods* (governing equations, MCE implementation, thermal
      losses, viscous dissipation, field-change schemes), not a compilation
      of experimental characteristic curves — so it still doesn't supply
      the Qc/COP-vs-span data this item originally asked for. That part of
      the original item remains open, but for a different, now-correct
      reason: the paper is here but isn't the right kind of source, rather
      than absent altogether.
- [~] Pixel-calibrated digitization of Tušek et al. (2013) Figs. 10-11 —
      **partial progress, not complete.** Extracted the two embedded
      images at native resolution (476x1093/1095 px; confirmed via
      pdfplumber that the page has 0 vector paths, so this has to be pixel
      work, not a vector-graphics shortcut). OCR (tesseract) confirmed,
      with 90-96% confidence, that each figure is 3 stacked subplots (one
      per AMR geometry A/B/F), each labeled "cooling capacity [W]" with a
      shared x-axis "temperature span [K]" ticked at 5/10/15 K — a real,
      checkable x-axis calibration. Row-wise dark-pixel density locates
      candidate y-gridlines at ~70px spacing. **Where it stopped:** y-axis
      numeric tick labels didn't OCR reliably at this resolution, so the
      y-scale rests on gridline spacing alone, not a confirmed absolute
      value; and series separation (9 overlapping lines per subplot, 3
      geometries x 3 flow ratios, grayscale with no color cues) was not
      attempted, since that requires sub-pixel marker-shape matching or a
      human with a point-and-click digitizer — auto-extracting a
      confident-looking 9-series table from blind pixel statistics here
      would be exactly the false-precision trap already flagged above, not
      genuine progress. Full notes in `results/tusek_ate2013_figs_notes.md`.
      Nielsen (2011) is in hand too, per the correction above, but is not a
      source of this kind of curve data, so it drops off this specific
      to-do.
- [~] **Full-system cost — partial, honest step, not the full item.** Searched
      for a published $ breakdown of AMR heat-exchanger/pump/motor-drive/
      controls capital cost and found none beyond the two Bjørk-group
      materials-cost studies already cited; that specific gap remains
      genuinely open, not resolved here. Confirmed this directly from the
      primary source itself: Bjørk et al., "Determining the minimum mass
      and cost of a magnetic refrigerator" (already in `Papers/Economics/`)
      states in its own Discussion section that motor and pump costs for
      the refrigeration system were not included in its analysis — the gap
      is an explicit, acknowledged
      scope limit of the source paper, not something missed in an earlier
      search. Also found a related-but-not-equivalent paper not yet in
      `Papers/`: Teyber, Trevizoli, Christiaanse, Govindappa, Niknia, Rowe,
      "Permanent Magnet Design for Magnetic Heat Pumps using Total Cost
      Minimization," J. Magn. Magn. Mater. 442 (2017) 87-96 — optimizes
      magnet topology against a combined capital+operating thermoeconomic
      cost-rate balance, so it does model pump-driven operating cost as
      part of its objective function, but it's still a magnet-topology
      cost study, not a bottom-up HX/pump/motor/controls hardware BOM.
      What *was* added previously: a second Bjørk-group paper, "The lifetime cost of a magnetic refrigerator"
      (Bjørk, Bahl & Nielsen, Int. J. Refrig. 63 (2016) 48-62), quantifies
      device *operating* cost (electricity, $0.10/kWh) over a stated device
      lifetime -- a real cost component previously entirely absent from
      `economics.py`. New `lifetime_cost()` combines this with the existing
      `material_cost()` floor and returns each piece separately. The 2016
      paper explicitly states, like its 2011 predecessor, that "actual
      manufacturing, transportation, maintenance and auxiliary systems are
      ignored" -- so HX/pump/motor/controls/enclosure hardware CAPEX stays
      unquantified; a genuine bottom-up BOM study for that hardware was not
      found and is not fabricated here.
- [x] **Optional COMSOL 2-D/3-D regenerator-bed setup guide — done as a
      guide, not as a run model.** No COMSOL license was available in this
      environment, so no 2-D/3-D solve was actually performed or
      validated; what's provided is a setup specification
      (`COMSOL_SETUP_GUIDE.md`) grounded directly in this repo's existing
      physics — same porosity/particle-diameter/Wakao-Kaguei correlation
      as `thermal.py`, the same post-Phase-8 entropy-solver code path for
      the MCE source term, and an explicit "degeneracy check" step
      (reproduce one calibrated benchmark device's Qc/COP against
      `core.validation_system`'s printed numbers before trusting anything
      new the 2-D model shows). Whoever picks this up next should treat
      every number in that guide as unvalidated until that check passes.
- [x] **Geometry-dependent pumping power (packed-bed AND parallel-plate),
      closing a real gap exposed by a newly-supplied paper.** Tušek,
      Kitanovski, Poredoš, "Geometrical optimization of packed-bed and
      parallel-plate active magnetic regenerators," Int. J. Refrig. 36
      (2013) 1456-1464 (`Papers/Optimization/Geometrical optimization of
      packed-bed and parallel-plate active magnetic regenerators.pdf`) —
      a DIFFERENT Tušek 2013 paper from the "comprehensive experimental
      analysis" one used elsewhere in Phase 7 above, already present in
      this repo. Confirmed directly (not assumed) that `thermal.py`'s
      pre-existing `regenerator_effectiveness(particle_diameter=...)` has
      no optimum: `eps` rises monotonically as particle diameter shrinks,
      because the heat-transfer (NTU) side and the pumping-loss side
      (`loss_model.py`'s `W_pump = k_pump*mdot**2`) were entirely
      decoupled — the latter has no particle-diameter or plate-geometry
      dependence at all — so the pre-existing model structure could not,
      even in principle, reproduce this paper's reported trade-off optimum
      (Table 3: sphere diameter 0.07mm Qc-optimal/0.17mm COP-optimal;
      parallel-plate spacing 0.035mm/0.075mm), despite those literature
      values sitting 3-7x below this repo's pre-existing hardcoded default
      (`particle_diameter=0.0005m`, i.e. 0.5mm). Added
      `pressure_drop_packed_bed()`/`pumping_power_packed_bed()` (Tušek et
      al.'s friction factor f=23.462·Re^-0.6716 and hydraulic diameter
      d_h=4·V_bed·eps/A_total, Eqs. 5&7) and, since no parallel-plate model
      of any kind existed in this repo before,
      `regenerator_effectiveness_parallel_plate()` (Nickolay & Martin 2002
      Churchill-Usagi-blended Nusselt correlation, Eq. 4 — read directly
      off the rasterized page image because `pdftotext` garbled this
      equation's exponents) and `pumping_power_parallel_plate()` (laminar
      f=24/Re, Eq. 6). New `core/geometry_analysis.py` couples these to
      `amr_cycle.py` and confirms a genuine interior COP optimum now
      exists vs. both packed-bed sphere diameter and parallel-plate
      spacing, at a fixed representative mdot. **A real complication found
      and reported rather than hidden:** initially tried re-optimizing
      mdot per geometry, mirroring the paper's own stated methodology, but
      found this repo's 2nd-law magnetic-work model
      (`W_mag=Qc·(Th/Tc-1)/eta_2nd_law`) makes single-objective COP
      maximization degenerate — COP rises monotonically as mdot (and Qc
      with it) falls toward zero, since `eta_2nd_law` saturates once `eps`
      hits its 0.97 clip. `check_free_mdot_cop_is_degenerate()` documents
      this; it's a separate, genuine reason (beyond the geometry gap this
      item closes) `optimize.py` is right to treat COP and Qc as competing
      multi-objective targets rather than a single COP search. The
      resulting optimum diameter/spacing found here should NOT be expected
      to numerically match the paper's own Table 3 — this repo's operating
      point (291K/10K span, 2kg, 1Hz, 1.5T) differs from the paper's fixed
      conditions (278-293K/15K span, 40x10xL mm outer dims, 0.5/3Hz), and
      the pumping power used is idealized hydraulic power with no pump/
      motor efficiency — the value is the qualitative confirmation that a
      genuine geometry optimum now exists in the model where none could
      before. New tests in `tests/test_thermal.py` and
      `tests/test_geometry_analysis.py`; full run saved to
      `results/geometry_optimization_analysis.txt`.

## Phase 8 — code-quality pass: solver performance, an entropy bug, and a test suite — done
- [x] **Fixed a real performance bug**: `mce_material.magnetization()`'s
      damped fixed-point solver took ~300-500 iterations per call near Tc
      (profiled: ~90k Brillouin evaluations for 200 calls). This was the
      reason `optimize.py`'s NSGA-III run (2,400 AMR evaluations) never
      finished in reasonable time. Replaced with Newton's method on the
      same self-consistency equation — ~5-15 iterations, ~19x faster per
      AMR evaluation, bit-identical M(T,H) values everywhere checked.
- [x] **Fixed a real physics bug this exposed**: `entropy_magnetic()`
      computed `log(sinh(a·x)/sinh(b·x))` by flooring the two sinh
      *arguments* independently at 1e-12. As x→0 (H=0 at/above Tc, where
      the self-consistent M is exactly 0) this silently collapsed the
      zero-field entropy to 0 instead of the correct max-entropy value
      `N·kB·ln(2J+1)` (~110 J/kg/K for Gd). The bug was invisible before
      because the old, imprecise solver never actually reached the exact
      M=0 solution — it stalled at a small nonzero residual that happened
      to avoid the floor. With the faster/more accurate Newton solver, the
      bug surfaced immediately: at the `optimize.py` operating point
      (T_cold=291K, span=10K → T_mid≈296K, just above Tc=294K) it produced
      ΔS=+105 J/(kg·K) (wrong sign, ~20x too large) and a *negative* total
      heat capacity, and the optimizer converged to a degenerate all-zero
      Pareto front. Fixed with a proper small-x series expansion of the
      log-sinh ratio. Net effect on the Dan'kov et al. (1998) validation
      (`validation.py`): 1T error moved from +13.6% to +48.9% (worse-looking
      but now correct — the old "better" number was two bugs partially
      cancelling), 5T improved from -9.7% to -7.5%. `optimize.py` now
      finishes in ~9s and returns ~35 non-degenerate Pareto-optimal designs
      instead of hanging.
- [x] Added a pytest suite (`tests/`, 28 tests, run in <2s) covering the
      core physics (saturation bounds, sign of ΔS, C_total>0, Dulong-Petit
      high-T limit, the two regressions above), the AMR cycle (COP
      ordering, NTU mass sensitivity, a solver-speed regression guard),
      and the supporting modules (economics, emissions, baseline COPs,
      loss-model calibration). Previously there was no test suite at all —
      only the ad hoc `validation.py` / `validation_system.py` scripts,
      which print numbers for a human to read rather than asserting on
      them. Run with `pytest -q` after `pip install -r requirements-dev.txt`.

## Phase 9 — added La(Fe,Si)13Hy to the materials library — done
- [x] **Added `LAFESIH_FIRST_ORDER`** (`core/first_order_mce.py`), the
      material that was previously missing entirely: `validation_system.py`
      had been running the `Astronautics_rotary_2014` benchmark row (Jacobs
      et al., Int. J. Refrig. 37 (2014) 84-91) against `GADOLINIUM` as an
      explicitly-flagged stand-in "because that material is not yet included
      in the material library." It now uses the real material -- same
      first-order Landau treatment as `GD5SI2GE2_FIRST_ORDER`, since
      La(Fe,Si)13Hy's giant MCE comes from an itinerant-electron metamagnetic
      (first-order-like) transition, not a second-order one, so it does NOT
      belong in the Brillouin/mean-field model (`mce_material.py`) alongside
      GADOLINIUM despite the surface-level "do it like Gd" framing -- see
      that module's own MODEL LIMITATION note for why. Composition/Tc/peak
      ΔS_M calibrated to La(Fe0.90Si0.10)13H1.1 (Tc=287K, |ΔS_M|~31 J/(kg·K)
      at 5T; Fujieda, Fujita & Fukamichi, Appl. Phys. Lett. 81, 1276 (2002)
      and related literature), by the same grid-search method used for
      `GD5SI2GE2_FIRST_ORDER`. Full parameter-by-parameter provenance and
      honesty flags (Debye temperature is an unsourced Fe-intermetallic
      placeholder; no Giguere-style ΔT_ad cross-check exists for this
      material, so no correction factor is applied) are in the block
      comment above the definition -- deliberately as detailed as the
      Gd/Gd5Si2Ge2 entries, per the request that motivated this phase.
- [x] **Wired it into `validation_system.py`**: `_material_for_row()` now
      selects `LAFESIH_FIRST_ORDER` for any benchmark row whose `material`
      column contains "La" (currently just `Astronautics_rotary_2014`),
      GADOLINIUM otherwise. Also added `T_COLD_LAFESIH_K=305K` (32°C
      measured cold-side inlet, Jacobs et al. 2014 Sec. 2) so this device
      stops silently reusing the Gd-centered `T_COLD_ASSUMED_K=289K` default,
      which was never appropriate for a naval-electronics cooler designed to
      run well above room temperature.
- [x] **Genuine finding, not hidden**: with the real material and the
      corrected T_cold, `Astronautics_rotary_2014` now reports "no
      calibration found" instead of the (misleading) numeric fit it got
      under the Gd stand-in. Root cause: the calibrated single-Tc=287K
      material's first-order transition is only a couple of K wide even
      after the field-induced shift (peak ΔT_ad ~298-299K at 1.44T), while
      this device's real bed is SIX layers graded 304-316K (Bahl et al.,
      Int. J. Refrig. 74 (2017) 22-29, for a similar Astronautics-style
      device) -- at the device's actual T_mid≈310.5K the model predicts
      ΔT_ad≈0.02K, nowhere near enough to cover an 11K span. This is the
      same "single material can't stand in for a graded multi-layer bed"
      limitation the codebase already treats as a first-class result
      elsewhere (see Phase 3's cascade.py and `first_order_mce.py`'s own
      `composition_tuned_material()`/`GIANT_MCE_TC_MIN_K`/`_MAX_K` machinery
      for exactly this kind of Curie-graded staging) -- it was NOT
      "fixed" by retuning Tc to force a fit, which would have fabricated a
      validation that doesn't exist. **Follow-up, same phase**: this WAS
      then actually done (see the addendum below) rather than left as a
      future item -- properly validating this device needed a 6-layer
      graded AMR model, and `cascade.py` already had the staging machinery
      for exactly this (previously hardcoded to the Gd5(SixGe1-x)4 family).
- [x] Updated `tests/test_validation_system.py`'s
      `test_run_system_validation_still_returns_four_point_results` for the
      above (7 calibrated original-set devices instead of 4 successes +
      Astronautics; total row count unchanged at 13) and reran the full
      suite (`pytest -q`, 55/55 passing) plus `core/first_order_mce.py`'s
      own `__main__` calibration self-check (now prints both Gd5Si2Ge2 and
      La(Fe,Si)13Hy calibration targets side by side).

**Addendum — generalized `cascade.py`'s Curie-graded bed to a pluggable
material family, and used it to actually test the single-layer finding
above:**
- [x] **Generalized `cascade.py`**: `_target_composition_for_peak()`,
      `run_graded_cascade()` and `compare_graded_cascade()` were hardcoded
      to the Gd5(SixGe1-x)4(-Ga) family (`composition_tuned_material()`,
      `GIANT_MCE_TC_MIN_K`/`_MAX_K`). Extracted a `GradedFamily` dataclass
      (tuned_fn, tc_min/tc_max, reference material, fallback material) so
      the same Curie-grading mechanism works for any composition-tunable
      first-order family. `GD_FAMILY` reproduces the original behavior
      exactly (default argument -- verified with
      `test_run_graded_cascade_gd_family_default_matches_explicit_family`);
      `LAFESIH_FAMILY` is new, backed by
      `lafesih_composition_tuned_material()` and
      `LAFESIH_TC_MIN_K`/`_MAX_K` (190-340K, a general literature-tunability
      reading, less precisely sourced than the Gd family's single-paper
      endpoints -- see `first_order_mce.py`'s comment for the caveat).
- [x] **`validate_astronautics_graded_bed()`** (`cascade.py`): builds the
      real 6-layer bed (layer Tc's spread evenly across the device's
      reported 303.6-316.2K range, since Jacobs et al. (2014) doesn't
      tabulate individual layer compositions) at the device's actual
      operating point, calibrates `fluid_mdot` to the reported Qc=2502W
      (same calibrate-then-validate methodology as
      `validation_system.py`), and compares the resulting COP against the
      reported 1.9. Result: **COP=1.69, -11.1% error** -- comparable in
      magnitude to the Gd-device errors `validation_system.py` already
      reports elsewhere, and a genuinely different outcome from the flat
      "no calibration found" the single-layer material gave this device:
      the graded-BED STRUCTURE, not just the material swap, was the
      missing piece. Still an approximation (evenly-spread layer Tc's, and
      COP is compared at calibrated mdot rather than predicted from
      scratch) -- stated as such in the function's own docstring, not
      oversold.
- [x] **Found and fixed two real bugs surfaced while building this**,
      both now regression-tested in the new `tests/test_cascade.py`:
      1. Replacing the original fixed-point Tc search with an
         `scipy.optimize.brentq` root-find (needed because the fixed-point
         update oscillates and fails to converge for LAFESIH_FAMILY's much
         narrower transition, silently leaving some stages at ~0 cooling
         capacity) initially searched down to `Tc_min-40K`, which for
         `GD_FAMILY` (Tc_min=20K) reached negative-Kelvin territory and hit
         a genuine numerical artifact in this Landau model: `DeltaT_ad`
         diverges as T->0K because the lattice heat capacity ->0
         (checked directly: T=1K gives `DeltaT_ad`~30K for
         `GD5SI2GE2_FIRST_ORDER`, a pre-existing model quirk that was
         simply never probed before since prior searches stayed near room
         temperature). Fixed by flooring the search range at 100K, comfortably
         clear of that artifact for both families.
      2. The `_peak_temperature()` search (a single 1401-point grid pass)
         became the actual runtime bottleneck once it was called from
         inside a root-finder instead of a fixed 6-iteration loop --
         `compare_graded_cascade()`'s full sweep went from tractable to
         several-minutes-plus. Replaced with a coarse-then-fine two-pass
         search (coarse full-range pass, then a zoomed fine pass around
         the coarse peak) -- same resolution near the peak, ~6x fewer
         `delta_T_adiabatic` evaluations. Full `python -m core.cascade`
         demo run: ~110s (was untested/unbounded before this pass; the
         module had no prior benchmark to regress against).
- [x] **Added `tests/test_cascade.py`** (9 tests) -- `cascade.py` had NO
      test coverage before this addendum. Covers: `GD_FAMILY` backward
      compatibility, both families' root-finder convergence (including a
      regression test for bug #1 above), both families' feasibility at
      realistic operating points, `composition_tuned_material()`/
      `lafesih_composition_tuned_material()` out-of-range errors, and the
      headline Astronautics graded-bed COP result. Full suite: `pytest -q`,
      **64/64 passing** (was 55/55 before this addendum).
- [x] Added the graded-bed validation to `cascade.py`'s own `__main__` demo
      output, alongside the pre-existing Gd5(SixGe1-x)4 sweep, AND to
      `main.py`'s end-to-end pipeline as a new **Step 7c** (a user running
      `python main.py` was otherwise seeing step 7b's Gd5(SixGe1-x)4 sweep
      but not this addendum's result at all -- `main.py` calls the standalone
      module functions directly, so a module-level addition doesn't
      automatically show up there without being wired in).
- [x] **Runtime note**: Step 7b (`compare_graded_cascade`, the 5-20K x
      1-4-stage Gd5(SixGe1-x)4 sweep) takes on the order of 1-2 minutes
      depending on hardware (observed: ~43s in this sandbox, ~126s on a
      user's Windows machine) -- each of its ~160 stage-instances needs a
      `scipy.optimize.brentq` root-find over a two-pass peak search (see
      the `_peak_temperature`/`_target_composition_for_peak` performance
      work above). Step 7c (the single Astronautics 6-layer validation,
      not a full sweep) takes ~15-25s on its own. Neither step was timed
      or bounded before this addendum (`cascade.py` had no prior benchmark
      to regress against); this is stated here as a known characteristic,
      not chased further given the accuracy fixes it was needed for.