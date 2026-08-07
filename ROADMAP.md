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
        read is in `data/tusek_ate2013_figs/notes.md` for whoever picks
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
- [x] Pixel-calibrated digitization of Tušek et al. (2013) Figs. 10-11 —
      **COMPLETE as of the Group A completion pass, see below.** Extracted the two embedded
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
      genuine progress. Full notes in `data/tusek_ate2013_figs/notes.md`.
      Nielsen (2011) is in hand too, per the correction above, but is not a
      source of this kind of curve data, so it drops off this specific
      to-do.
      **UPDATE (Group A completion pass):** the above "where it stopped"
      blockers are both resolved. y-axis gridlines calibrate cleanly (8/6/
      4/2/0 W for Fig.10, 40/30/20/10/0 COP for Fig.11, ~35px/gridline-step,
      consistent across all 6 panels) once verified against the actual
      extracted bitmap rather than OCR of the tick labels. Series separation
      was done by a human, one panel and one marker at a time, exactly as
      flagged as unavoidable here — using automated blob-centroid detection
      (thresholding + erosion + connected components) to get precise pixel
      coordinates, then visually confirming/correcting each point's series
      assignment against zoomed crops (catching several line-crossing
      artifacts that were not real markers). All 9 series in both figures
      are now digitized: `data/tusek_ate2013_figs/fig10_data.csv`,
      `fig11_data.csv`, full methodology and a documented residual
      uncertainty (including one flagged, unresolved point-count mismatch
      between the two figures at V*=0.95/AMR F) in `notes.md`.
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

## Phase 10 — Paper-mining pass: blow-fraction asymmetry & (Mn,Fe)2(P,Si) family — done
Cross-referenced this repo's current state against `Papers/` for content the
roadmap had not yet mined; see `Paper_Mining_Recommendations.md` for the full
pass. Two of its four findings were concrete, numerically-anchored model
additions; the other two are documentation-only flags (below).

- [x] **Flow-waveform asymmetry (blow fraction)** added to `amr_cycle.py`
      (`AMRSystem.blow_fraction`, `BLOW_FRACTION_MASCHE`,
      `_blow_fraction_multiplier()`) — a real degree of freedom the model
      previously had no notion of at all (Qc/second-law efficiency
      implicitly assumed a symmetric 50/50 cold-to-hot/hot-to-cold split).
      Calibrated to the ONLY two-point comparison in the source paper
      (Masche, Liang, Engelbrecht & Bahl, Appl. Thermal Eng. 215 (2022)
      118945 — DTU rotary AMR, 13 trapezoidal beds, 295g Gd spheres/bed, at
      T_span=16K/U=0.32/f=1.4Hz): blow fraction 25.0%→41.6% raised Qc
      70W→330W (4.7x) and second-law efficiency 2.6%→17.4%. Default
      `blow_fraction=0.5` exactly reproduces every pre-existing result
      (multiplier=1.0 at the symmetric baseline); verified numerically that
      the model reproduces the reported 4.71x relative Qc swing between the
      two tested blow fractions. Honesty flags (stated in
      `_blow_fraction_multiplier`'s docstring, not smoothed over): only two
      points at one operating condition are available, so the parabola's
      shape and its extrapolation to other T_span/U/frequency combinations
      are unvalidated; the "best found" point is treated as the curve's
      true peak, which the source paper does not itself claim.
      Also added as a **6th NSGA-III decision variable** in `optimize.py`
      (bounds [0.1, 0.6], bracketing/widening the paper's tested
      0.25-0.416 window) — the optimizer independently converges toward
      blow_fraction≈0.37-0.43 across the Pareto front, landing close to the
      paper's own 0.416 "best found" value, a reasonable sanity check on
      the calibration (not itself independent validation, since it's the
      same calibration data driving both).
      `tests/test_amr_cycle.py`: 3 new tests (default-is-symmetric,
      reproduces the 4.71x Qc swing, best-found beats symmetric-default).
- [x] **Added `(Mn,Fe)2(P,Si)` as a third pluggable `GradedFamily`**
      (`MNFEPSI_FAMILY` in `cascade.py`, `MNFEPSI_FIRST_ORDER` +
      `mnfepsi_composition_tuned_material()` in `first_order_mce.py`),
      alongside `GD_FAMILY`/`LAFESIH_FAMILY` — confirmed the Phase 9
      `GradedFamily` generalization is genuinely pluggable: no changes to
      `run_graded_cascade()` or `_target_composition_for_peak()` were
      needed. Source: Hanggai, Yibole, Guillou, Kwakernaak, van Dijk &
      Brück, Acta Materialia 302 (2026) 121677 — Fe-rich melt-spun
      Mn0.60+xFe1.3-xP0.66-ySi0.34+y ribbons (0<=x<=0.08, x=2y). Tc window
      295.3-331.2K is DIRECTLY MEASURED across five real compositions
      (Table 1 of the source paper) — unlike GD_FAMILY's Ga-alloying
      endpoint or LAFESIH_FAMILY's general literature-survey window — and
      sits almost entirely AT OR ABOVE the ASHRAE 291.15-300.15K
      data-center range (the parent composition's 295.3K already falls
      inside it), the opposite tension from GD_FAMILY's ceiling sitting
      just below that range. Calibrated (A,B,C)=(1.16,-0.464,0.928) —
      same B/A=-0.4, C/A=0.8 ratio as `GD5SI2GE2_FIRST_ORDER`/
      `LAFESIH_FIRST_ORDER`, only A rescaled — by grid search to
      reproduce the source paper's cross-validated (calorimetric AND
      magnetization) peak |DeltaS_M|~17.6 J/(kg K) at 2T (NOT the 5T used
      for the other two families' calibration targets; this is the field
      the source paper's own measurements were made at). No Giguere-style
      dTad_correction is applied — same honesty flag `LAFESIH_FIRST_ORDER`
      carries: the source paper reports DeltaS_M (indirect, Maxwell
      relation/calorimetric), not a direct DeltaT_ad measurement, and no
      independent cross-check paper for this family is in the corpus.
      `tests/test_first_order_mce.py` + `tests/test_cascade.py`: 6 new
      tests (peak-DeltaS calibration match, dTad_correction defaults,
      out-of-range rejection, Tc-only shift assumption, root-finder
      convergence, in-range graded-cascade feasibility). Full suite:
      **137/137 passing** (was 125/125 before this pass).
- [x] **Data now available; geometry-explicit calibration still NOT
      wired up — partial.** `Papers/AMR Theory and Modeling/...Tušek...
      2013 comprehensive experimental analysis...` reports a 20K span,
      1.15T, ~25%-porosity parallel-plate AMR (the largest reported
      parallel-plate span at that field at time of publication) — this is
      the ONLY candidate in the corpus that could validate
      `geometry_analysis.py`'s `regenerator_effectiveness_parallel_plate()`
      model against a real device (every row in
      `data/amr_experimental_benchmarks.csv` today is packed-bed or
      layered packed-bed). The exact (Qc, COP) pair needed for a full
      calibrate-then-validate CSV row is only in Figs. 10-11, the SAME
      digitization already flagged as open in Phase 7 — this entry exists
      only to name the parallel-plate validation gap as the SPECIFIC
      target for whoever picks that digitization up, not to re-flag the
      digitization itself.
      **UPDATE (Group A completion pass):** the digitization is done (see
      Phase 7's entry above) and the paper's directly-stated 19.8K/0W
      span-ceiling point (Section 3.2/Fig.6, V*=0.365, freq=0.3Hz) was
      added as `Tusek_singlebed_Gd_2010_spanceiling` in
      `amr_experimental_benchmarks.csv`. **What is still NOT done, and
      remains genuinely open:** actually calling
      `regenerator_effectiveness_parallel_plate()` with this device's own
      geometry (Table 1: 0.1mm plate spacing, 0.25mm plate thickness,
      dh=0.2mm, porosity=0.2564, heat-transfer area=0.1395 m²) and
      comparing its predicted effectiveness against an effectiveness
      backed out of the measured (span, Qc) data — that needs
      `core/geometry_analysis.py` wiring plus a fluid-property assumption
      for the water/ethylene-glycol mixture (density/viscosity/specific
      heat), neither of which was fabricated here. This row only supplies
      the previously-missing data point; the geometry-explicit validation
      itself is a good next step for whoever picks this up.
- [x] **Confirmed duplicate, no action**: `Papers/AMR Theory and
      Modeling/Performance evaluation of a nine-layer active
      regenerator.pdf` and `Papers/AMR systems and prototypes/The
      performance of a large-scale rotary magnetic refrigerator.pdf` are
      the same accepted-manuscript PDF (Jacobs, Auringer, Boeder, Chell,
      Komorowski, Leonard, Russek & Zimm, Int. J. Refrig., DOI
      10.1016/j.ijrefrig.2013.09.025 — identical title/DOI/abstract,
      already cited in `data/amr_experimental_benchmarks.csv` as
      `Astronautics_rotary_2014` with its actual 3042W/2502W@11K-span
      numbers). Noted here only so neither filename is searched again for
      "nine-layer" content that isn't actually in it.


## Phase 11 — Paper-mining pass, Part 2: deeper search — done
Went back through the remaining unmined papers a second time (tables,
in-text numeric callouts, secondary reviews' own citation tables) rather
than abstract-level skims; see `Paper_Mining_Recommendations_Part2.md` for
the full pass. Three findings became concrete additions; one is a
documentation-only note (below).

- [x] **Extended `validation.py`'s Gd checks to 7T** using Giguere et al.
      (1999)'s own pure-Gd methods-section cross-check paragraph (a paper
      already in this repo, previously mined only for its Gd5Si2Ge2 content
      by `giguere_validation.py`) — `GIGUERE_GD_CROSSCHECK`,
      `run_giguere_gd_extension()`. Confirmed the paper's own numbers
      directly from the PDF (not just the recommendations doc): "on high
      purity Gd agrees with the value of AMES laboratory within 1K (10.5
      and 11.5K, respectively, both for 5T)... For 7T, our value (12 and
      13K for industrial- and high-purity Gd, respectively) agrees well
      with that of Brown (14K)." Honest result: the model (calibrated to
      Dan'kov et al.'s HIGHER 5T value, 14.6K) overestimates relative to
      Giguere et al.'s Gd range at both 5T (13.5K model vs. 10.5-11.5K,
      +22.8% vs. midpoint) and 7T (16.7K model vs. 12-13K, +33.7% vs.
      midpoint) — a real, reported disagreement between two published Gd
      measurements, not a bug in this repo; both numbers are reported, not
      reconciled.
- [x] **Added a Curie-point field-shift check** against Dan'kov et al.
      (1998)'s own reported "~6 K/T above 2T, up to 7.5T" rate —
      `DANKOV_CURIE_SHIFT_RATE_K_PER_T`, `run_curie_shift_check()`. This is
      a genuine held-out prediction check (the peak-of-ΔT_ad(T) location is
      an EMERGENT output of mce_material.py's self-consistent M(T,H) Newton
      solve, not a hardcoded input) — and it does NOT pass: using a
      bounded-Brent sub-Kelvin-precision peak locator (not a coarse grid,
      to rule out a resolution artifact) across 12 points from 2-7.5T, the
      model's peak-ΔT_ad temperature comes out pinned at 294.5K regardless
      of field — a fitted shift rate of ~0 K/T, not ~6 K/T. Documented as a
      genuine limitation of this specific mean-field/Weiss-molecular-field
      Brillouin formulation (see run_curie_shift_check()'s docstring for
      why), not smoothed over or hidden.
- [x] **Added the Chubu Electric/Toshiba two-field-point row** to
      `data/amr_experimental_benchmarks.csv` (`ChubuToshiba_Gd_2016_4T`/
      `_2T`) — a genuine field-sensitivity pair (4T→100W/26K,
      2T→40W/24K at fixed geometry), the only such pair in the benchmark
      set (every other device here is single-field). SECONDARY SOURCE,
      same caveat as the existing `Okamura_Hirano_2013` row: cited via
      Kamran, Ahmad & Wang, Renew. Sustain. Energy Rev. 133 (2020) 110247,
      Table 2 (original ref [69] not in this repo's `Papers/`). The
      source table gives volume (V_reg=484cm³, n_reg=2) and no mass or
      COP — `mass_MCM_kg` is therefore an ESTIMATE (V_reg×n_reg×RHO_GD×
      (1-porosity), reusing `core/thermal.py`'s own existing RHO_GD=7900
      kg/m³ and porosity=0.365 constants, not an invented number), and
      COP is left blank (capacity/span-only row, same treatment as this
      CSV's existing zerospan/maxspan rows).
      Added `run_field_sensitivity_check()` to `validation_system.py` —
      the field-axis analog of the existing `run_curve_validation()` —
      to actually EXERCISE this pair rather than leave it decorative:
      calibrates mdot to the 4T anchor point, reuses that calibrated
      system at 2T, and compares the predicted Qc against the 2T
      companion row. Honest result: the 4T anchor (26K span, 4.856kg Gd,
      0.167Hz) itself does NOT calibrate — `dTad_noload` at 4T for this
      mass/frequency is too small relative to a 26K span for ANY mdot in
      [1e-6, 5] kg/s to reach the reported 100W, so the check correctly
      reports "no calibration found," the SAME honest outcome already
      documented for `Risoe_DTU_Gd_2011` and (pre-Phase-9)
      `Astronautics_rotary_2014`. Kept in the codebase anyway: this is a
      real, informative finding about the model's achievable-span ceiling
      at this device's scale, not a wasted addition.
      `tests/test_validation.py` (new file, 8 tests) and
      `tests/test_validation_system.py` (+4 tests, +1 exemption-logic
      generalization for the pre-existing device-group identity test).
      Also wired both `validation.py` extensions and
      `run_field_sensitivity_check()` into `main.py`'s pipeline (stages 1
      and 2) so they run as part of the normal full pipeline, not just
      standalone. Full suite: **148/148 passing** (was 137/137 before this
      pass).
- [ ] **Flagged, not built** (everything else in Part 1's §1 comparative-
      prototype table): the Institute of Tech./Chubu near-zero-span
      extreme (540W at 0.2K span, 844cm³×4 regenerators — a genuinely
      different high-mass operating regime than anything currently
      calibrated) and the second, independent Riso Lab data point (23cm³,
      24 regenerators, 2.25Hz, 1.24T — does NOT match the existing
      `Risoe_DTU_Gd_2011` row's parameters, so likely a different Risø/DTU
      device generation, not a duplicate) were both left OUT of the CSV.
      Neither was on Part 2's own "updated priority list" (only the Chubu
      Electric/Toshiba pair was), and the Riso point's primary source
      (refs [18,50,62,72-74] in the Kamran/Ahmad/Wang review) isn't in
      this repo's `Papers/` — adding it would need identifying which of
      five possible citations is the real source first. Noted here so
      neither gets silently forgotten, not built without a clearer
      go-ahead on scope.


## Phase 12 — Paper-mining pass, Part 3: remaining papers + reference books — done
Went through everything not yet touched in Parts 1-2: the economics paper in
full, both "Development of a rotary..." papers, the 1997 discovery paper,
the solid-state caloric cooling review, and both reference books; see
`Paper_Mining_Recommendations_Part3.md`. Two findings became concrete
additions; three are documentation-only notes (below).

- [x] **Added Cooltech 2013 (42K span stress test) and DTU MagQueen
      (LAFESIH cross-check) rows** to `data/amr_experimental_benchmarks.csv`.
      Both confirmed directly from the source PDF (Greco, Aprea, Maiorino
      & Masselli, Int. J. Refrigeration (2019), Table 2 AND body text, not
      table-only) before adding, not just taken from the recommendations
      doc's summary.
      `Cooltech_2013_rotary`: 42K span is the largest in this benchmark
      set (next is Risoe_DTU_Gd_2011 at 30K). No mass reported in the
      source (unlike Part 2's Kamran table, this one has no volume/mass
      column at all) -- left blank, falls back to
      `calibrate_and_check()`'s existing mass=1.0kg default, flagged as
      illustrative-only. No COP reported either (capacity/span-only row).
      Its own stress test does NOT calibrate at any mdot in [1e-6,5]kg/s
      -- consistent with the model's existing struggles at large spans
      (Risoe 30K also fails to calibrate).
      `DTU_MagQueen_2018`: the source paper reports this is a HEAT PUMP
      (Qh=1500W heating power, COP_h=5), not a Qc/COP_c pair -- Qc_W=1200
      and COP=4.0 in the CSV are DERIVED via the standard Qh=Qc+W identity
      (Qc=Qh*(1-1/COP_h), COP_c=COP_h-1), clearly flagged as derived, not
      measured, in the row's own source note. Material field
      ("La(Fe,Mn,Si)13Hz spheres") correctly routes to LAFESIH_FIRST_ORDER
      via the existing `_material_for_row()` "La" substring match -- the
      first LAFESIH-material benchmark point independent of
      Astronautics_rotary_2014. Also does NOT calibrate (same class of
      finding as Astronautics itself).
      Added `run_capacity_only_calibration_check()` to
      `validation_system.py` -- `run_system_validation()` silently skips
      any row without a reported COP (correct for COP comparison, but it
      meant Cooltech's 42K stress test would otherwise produce NO reported
      result at all). Reuses the existing `_calibrate_mdot()` (span/Qc/
      field/mass/frequency only, no COP needed) to report calibration
      reachability for every capacity-only row in the CSV, not just the
      new ones. Wired into `main.py`'s pipeline stage 2 alongside the
      other two validation_system.py checks.
      `tests/test_validation_system.py`: +6 tests, +1 hardcoded-count-test
      update (13->14 total rows in run_system_validation(), since
      DTU_MagQueen_2018 now has a derived COP and enters that count as a
      "no calibration found" result).
- [x] **Cross-checked the Gd5Si2Ge2 ΔT_ad correction factor** against
      Pecharsky & Gschneidner (1997)'s own text (confirmed directly from
      the PDF, garbled font-encoding notwithstanding): "the DeltaT_ad
      values of Gd5Si2Ge2 are larger than the corresponding DeltaT_ad
      values for Gd by about 30%, comparing the peak values, regardless
      of the temperature" -- a SECOND independent primary source
      (heat-capacity-based, not pulse-field-thermometry-based like
      Giguere et al. 1999) AND a second, independent field range (2T/5T,
      confirmed from the same paper's Fig. 6, vs. the single 7T point
      `DTAD_CORRECTION_FACTOR` was fit to).
      Added `run_pecharsky_ratio_check()` to `core/giguere_validation.py`.
      Genuinely interesting, non-obvious result: the RAW (uncorrected)
      model's peak-DeltaT_ad(Gd5Si2Ge2)/peak-DeltaT_ad(Gd) ratio at 5T
      (~1.24) lands close to Pecharsky & Gschneidner's ~1.30 -- an
      unexpected agreement, since the model was never fit to this number.
      But applying `DTAD_CORRECTION_FACTOR` (fit at 7T) drags that SAME
      ratio down to ~0.51 at 5T (and ~0.87 at 2T) -- i.e. the CORRECTED
      model predicts Gd5Si2Ge2 underperforms plain Gd, contradicting the
      entire "giant" MCE premise. This is not treated as grounds to
      change `DTAD_CORRECTION_FACTOR` (re-fitting a two-point correction
      would just repeat the same single-point-calibration problem one
      level up) -- it's documented as concrete evidence that the
      correction should not be read as field-independent, added to
      `first_order_mce.py`'s existing honesty-flag block as item 3.
      `tests/test_giguere_validation.py`: +5 tests.
- [x] **Checked against the primary source, "18K" not confirmed.**
      (Paper-Mining Pass Part 6 follow-up.) `txt`-extracted the actual
      Jacobs et al. paper ("The performance of a large-scale rotary
      magnetic refrigerator," now in this repo's Papers/) and searched it
      directly for any 18K span mention: none found. Its own reported
      spans are 11.0 K (headline, 2502W), 12.0 K (design target/measured
      operating point, mentioned four times), and 16.0 K (a single mention,
      in a measurement-vs-model-prediction comparison context, not
      presented as a device capability figure). "18K" does not appear
      anywhere in the primary text. Confirms the earlier caution was
      warranted -- the Greco et al. review table's "18K" entry is not
      independently verifiable from the source it's presumably drawn from,
      and should NOT be used as a third Astronautics data point. No CSV or
      code change (there wasn't one to make -- this closes the "should I
      trust this number" question with "no," not "yes, here are the
      values to add").
- [x] **Confirmed identity of the two "rotary refrigerator development"
      papers, no code impact**: `Development of a novel rotary magnetic
      refrigerator.pdf` = Lozano et al. (2016), already the primary
      source behind the `Lozano_POLO_UFSC_2016` CSV rows -- no new
      content. `Development of a rotary magnetic refrigerator.pdf` =
      Tušek, Zupan, Šarlah, Prebil & Poredoš (2010), an EARLIER Ljubljana
      prototype's mechanical/magnet-design paper with zero Qc/span
      numbers -- design and "pros and cons" (shaft-seal leakage,
      magnet-structure weight, assembly complexity) only. Noted here so
      nobody expects performance data from this specific file; added to
      `Literature_Review.md` as qualitative engineering-realism context
      (real mechanical parasitic losses `amr_cycle.py`'s idealized cycle
      doesn't capture), not a numeric validation source.
- [x] **Economics paper (Bjørk et al. 2011) confirmed fully mined, no
      action**: checked the full text, not just the abstract already
      behind `economics.py`. The only content not already reflected in
      `COST_MCM_PER_KG`/`COST_MAGNET_PER_KG` is Fig. 9 (cost vs. operating
      frequency, qualitative trend only, not a digitizable value) --
      wouldn't move the already-documented open BOM-cost gap even if
      digitized. No code change.
- [ ] **Reference books, flagged not acted on**: Kitanovski et al.,
      *Magnetocaloric Energy Conversion* (2015) -- the corpus copy is
      front-matter + ~9 pages of Ch.1 only (not the full book); Chapters
      4/7/9 (AMR performance, prototypes-by-country, costs) are listed in
      the TOC but not present. Tishin & Spichkin, *The Magnetocaloric
      Effect and its Applications* (2003) -- 486 pages, scanned images
      with NO OCR text layer (confirmed via `pdfplumber`, zero pages
      return extractable text); likely the richest materials-property
      compendium in the corpus, but OCR'ing the whole book wasn't
      attempted given the effort/yield tradeoff vs. primary sources
      already mined. Per the recommendations doc: if specific data tables
      from this book are wanted, OCR can be run on a targeted page
      range/topic rather than the whole 486 pages -- not done here absent
      a specific target.

## Phase 13 — Paper-mining pass, Part 6: DTU_rotary_Gd_2016 citation traced and corrected — done
Parts 4 and 5 (referenced throughout `core/loss_model.py` and
`data/amr_experimental_benchmarks.csv`'s row notes, but never logged here)
fixed stale `mdot` values in `CALIBRATION_POINTS_CORE` and added the
`DTU_Eriksen_rotary_Gd_2015` row from a genuine primary source,
respectively -- noted here for the record since this phase builds
directly on both. Part 6 tracked down and fixed the one citation flagged
throughout the repo as "unverified/unlocated": the `DTU_rotary_Gd_2016`
row (818 W, 10.1 K span, COP=4.2), cited only as "Bahl/Eriksen/
Engelbrecht, rotary AMR - ScienceDirect (2016)".

- [x] **Traced the real paper and found the CSV numbers don't match it.**
      The user supplied Eriksen's 2016 DTU PhD thesis ("Active magnetic
      regenerator refrigeration with rotary multi-bed technology"), whose
      Chapter 6 is exactly the paper behind the citation -- D. Eriksen,
      K. Engelbrecht, C.R.H. Bahl, R. Bjørk, "Exploring the efficiency
      potential for an active magnetic regenerator," Sci. Technol. Built
      Environ. 22(5) (2016) 527-533 (ref [20] in Masche et al. 2021,
      already in this repo's Papers/). Confirmed directly from the thesis
      text, not a secondary source: "a maximum second-law efficiency of
      18% was obtained at a cooling load of 81.5 W, resulting in a
      temperature span of 15.5 K and a COP of 3.6" at fAMR=0.61 Hz,
      1.13 T, 1.7 kg Gd -- not 818 W/10.1 K/COP=4.2 at 1.4 Hz/1.44 T. The
      real device ("MAGGIE") is the SAME physical prototype as the
      already-present `DTU_Eriksen_rotary_Gd_2015` row (Eriksen et al.,
      Int. J. Refrigeration 2015), just a later paper reporting a
      different (lower-frequency, higher-span, higher-COP) operating
      point -- confirmed by matching Chapter 3's device description
      (12x NdFeB magnet blocks, 11-compartment Curie-graded Gd +
      Gd(100-x)Yx regenerator, 1.13 T, 1.7 kg) against the 2015 paper's.
- [x] **Extracted a directly-measured loss breakdown as a bonus**: the
      thesis's Table 6.2 (Sec. 6.5.3) gives shaft power (14.0 W), total
      pumping power (8.9 W, split 4.1 W regenerator + 4.8 W external
      components), bearing/gear friction (1.6 W), valve friction (1.1 W),
      and Carnot (ideal) work (4.0 W) for this exact operating point --
      all directly reported, not back-calculated from Qc/COP like every
      other row in this benchmark set. Reported flow rate (V=2.5 L/min)
      also gives a directly-measured mdot, in principle avoiding the
      brentq back-calculation entirely -- flagged in the CSV row's note
      for future use, but not adopted as the CORE calibration input this
      pass (see next item for why).
- [x] **Checked the corrected point against this repo's own cycle model
      and found it doesn't calibrate.** At 1.13 T / 1.7 kg / 0.61 Hz,
      `amr_cycle.py`'s `cooling_capacity()` predicts Qc ≈ 0 at a 15.5 K
      span for ANY mdot in [1e-6, 5] kg/s -- the model's own zero-flow
      no-load span at this field/frequency already sits below 15.5 K.
      Same failure mode already documented for `Risoe_DTU_Gd_2011`
      (attributed there, and here, to the real device's Curie-graded
      11-layer bed reaching spans a single-uniform-Tc Gd approximation
      structurally cannot). Renamed the row `DTU_Eriksen_MAGGIE_2016` and
      kept it in the CSV with "no calibration found" status via the
      existing `calibrate_and_check()` path -- documented, not dropped.
- [x] **Replaced the CORE calibration point rather than deleting it.**
      Since the corrected DTU number can't fill `CALIBRATION_POINTS_CORE`'s
      3rd slot, `DTU_Eriksen_rotary_Gd_2015` (already in the CSV, already
      cited to a real primary source, and confirmed to calibrate cleanly:
      mdot=0.084666 kg/s reproduces its own literature Qc=102.8W exactly)
      was promoted into that slot instead of leaving CORE underdetermined
      or fabricating a new point. Re-ran the same
      `brentq(qc_residual, 1e-6, 5.0)` / `Wp=Qc*(1/COP_lit-1/COP_ideal)`
      procedure every other CORE point uses: mdot=0.084666 kg/s,
      Wp_required=26.18 W. NNLS refit stays non-negative
      (k_eddy=30.52, k_pump=0, base_frac=0.048).
- [x] **Found and corrected a downstream consequence**: with the DTU
      point's parasitic fraction changed from a fabricated 17.1% to a
      verified 25.5%, the 4-device EXTENDED-set parasitic-fraction-vs-Qc
      ranking (Tušek 11.7% < DTU 25.5% < Okamura 36.7% < Astronautics
      45.3%) flips from "non-monotonic" to a clean *monotonic increase*
      with device scale -- the opposite direction from the fixed-overhead/
      economies-of-scale hypothesis this analysis was built to test, so
      the qualitative conclusion (a size/scale term isn't supported) is
      unchanged, but the specific "non-monotonic" claim in
      `loss_model.py`'s docstring, `README.md`, and
      `magcool-dc_technical_walkthrough.md` was wrong post-correction and
      has been rewritten in all three places.
      `analyze_parasitic_fraction_scaling()` gained an explicit
      informational branch for the now-monotonic case instead of only
      ever printing the old "not monotonic" narrative.
- [x] **Test suite updated for the correction**: `tests/test_loss_model.py`
      (`_SELF_CONSISTENCY_SPANS`/`_MASS` dicts, the Lozano leave-one-out
      clustering threshold widened from 10.0 to 12.0 points to reflect the
      slightly shifted pooled fit, and the monotonicity test flipped from
      `test_parasitic_fraction_scaling_is_not_monotonic` to
      `test_parasitic_fraction_scaling_is_monotonically_increasing_with_qc`)
      and `tests/test_validation_system.py` (`with_cop` count drops from
      7 to 6 since the fabricated point used to "calibrate" and the real
      one doesn't, `total` stays at 15; the DTU curve-validation
      companion test is retired since no verified same-frequency
      companion span exists for the corrected point, and
      `test_curve_validation_covers_multi_point_groups`'s expected group
      set drops `DTU_rotary_Gd_2016`) -- full suite re-run and passing
      (148 collected, 146 pass; the 2 failures are pre-existing
      SALib/pymoo-optional-dependency sandbox issues unrelated to this
      change, not new).
- [ ] **Not done this pass**: the thesis's Ch.5 addendum reports a later
      no-load span of 29.2 K at fAMR=1.4 Hz for the SAME prototype --
      flagged in the CSV row's note as a possible 3rd MAGGIE data point
      for anyone picking this up, but not added as its own row since it's
      at yet another frequency than either existing MAGGIE row and can't
      serve as a same-condition curve companion for either. The directly-
      measured mdot (from the 2.5 L/min flow rate) and Wp (from Table
      6.2's power breakdown) were also not adopted as an alternative,
      non-back-calculated CORE calibration input -- doing so would break
      with every other CORE point's "hardcoded mdot reproduces hardcoded
      Qc under the current model" convention (since this point's real
      span doesn't reproduce under the model at all) and was judged too
      large a methodological change to make silently; flagged for a
      future pass if it's judged worth special-casing.
## Phase 14 — Bug fixes + Track A2 four-way material comparison + open-item decisions (this pass)

Closes ROADMAP.md's suggested "Track A" bug-fix list and makes explicit,
documented decisions on the remaining "Track B" open items (except B5,
left open — being digitized manually in WebPlotDigitizer) rather than
leaving them silently unaddressed.

- [x] **Fixed a real material-model bug**: `run_cascade_comparison()`
      (main.py step 7), `core/plots.py`'s fig20, and `core/cascade.py`'s
      `__main__` demo block all imported the mean-field
      `core.mce_material.GD5SI2GE2` (explicitly documented in that module
      as "retained as a parameter entry only," not for quantitative use)
      instead of the physically-appropriate first-order Landau model,
      `core.first_order_mce.GD5SI2GE2_FIRST_ORDER`, that
      `giant_mce_analysis.py` correctly uses. Fixed all three call sites.
      Re-ran the affected stage: `results/cascade_comparison_giant_mce.csv`
      and `fig20` are numerically UNCHANGED after the fix, because both
      materials share the same Tc=276K and the whole point of that
      comparison is that Gd5Si2Ge2's fixed Tc sits below the ASHRAE range
      regardless of which model computes it — already correctly documented
      in fig20's own title. So: a real bug, correctly fixed, but not one
      that silently changed a headline number this time.
- [x] **Built the four-way material family comparison**
      (`core/material_family_comparison.py`, wired into `main.py` as new
      step 8d, `plots.py` fig26): runs Gd, Gd5Si2Ge2 (fixed composition),
      and all three composition-TUNABLE giant-MCE families this repo
      already has (`cascade.py`'s `GD_FAMILY`/`LAFESIH_FAMILY`/
      `MNFEPSI_FAMILY`) through the same ASHRAE operating point(s)
      (T_cold=18C, 5-20K span sweep, 1-4 stage cascade), each family
      re-tuned per span to its own best composition via the same
      `_target_composition_for_peak` root-finder the graded-cascade code
      already relies on. Reports whether each family's documented Tc
      window actually covers the point needed, falling back to plain Gd
      where it doesn't (MNFEPSI_FAMILY's 295.3-331.2K window sits mostly
      AT/ABOVE the ASHRAE range and fails to cover the representative
      10K-span point; GD_FAMILY's 20-290K window is right at its ceiling
      there too). **Result at the representative 10K-span point**:
      La(Fe,Si)13Hy (tuned) ranks best (COP=7.33, Qc=4989W), ahead of
      plain Gd (COP=5.09, Qc=1443W) — this is the "which is best" ranking
      the original item asked for, previously left implicit across
      several separate analyses. See `results/material_family_comparison.csv`
      / `.txt` / fig26.
- [x] **A3 — span_fraction hard clamp**: kept the documented linear clamp
      (option (a) from this pass's own plan) rather than inventing an
      unsourced smoothing function; no literature source for the exact
      near-span-limit fall-off shape was found in this project's corpus.
      Documented the limitation directly in `AMRSystem.cooling_capacity()`'s
      docstring instead of leaving it implicit.
- [x] **A4 — repo-wide grep + full suite re-run**: confirmed no other
      quantitative use of the mean-field GD5SI2GE2 remained
      (`grep -rn "GD5SI2GE2\b" --include=*.py . | grep -v _FIRST_ORDER`
      now only matches the metadata test and the definition itself). Full
      suite: 170 passed, 0 failed (up from 168 pre-pass + 2 new test files
      for `material_family_comparison.py`; `test_plots.py`'s
      `FIGURE_FUNCTIONS`-vs-`run_all()` sync check caught fig26 needing
      registration there too, which it now has).
- [x] **B6 — full-system BOM cost: already closed, re-confirmed, no change
      needed.** `core/economics.py` already documents its CAPEX figure as
      a materials-only floor (Bjørk et al. 2011/2016), explicitly NOT a
      bottom-up HX/pump/motor/controls BOM, and states that gap in both
      the module docstring and `main.py`'s logged output. This matches
      option (b) from this pass's own plan (state CAPEX as a lower bound,
      not full TCO) — it was already the repo's position; getting real
      vendor quotes to close it fully is a data-gathering task outside
      what this pass can do.
- [x] **B7 — reference books (Tishin & Spichkin, Kitanovski et al.): no
      action, correctly left flagged.** Both are present in `Papers/
      Reference Books/`; per the original item's own criterion ("only
      worth OCR'ing if you have a specific data table you need from
      them"), no specific table was requested this pass, so they remain
      flagged rather than spec­ulatively OCR'd.
- [x] **B8 — two flagged-not-built CSV rows (Chubu near-zero-span extreme,
      second Risø/DTU point): no action, correctly left unbuilt.** No
      primary source for the Risø point was obtained this pass (still only
      known via a review paper's reference list per the existing
      ROADMAP.md entry) — per the original item's own instruction, not
      adding it from a secondary source alone.
- [x] **B5 — Tušek et al. (2013) Figs. 10-11 digitization: CLOSED (Group A
      completion pass).** Digitized manually (pixel-calibrated gridline
      detection + automated marker-centroid extraction, each point
      verified against a zoomed crop by eye) rather than in an external
      WebPlotDigitizer session — same end result, no external tool
      dependency. All 9 series across both figures are in
      `data/tusek_ate2013_figs/{fig10_data.csv,fig11_data.csv,notes.md}`.
      Wired into `validation_system.py`: the `Tusek_singlebed_Gd_2010` CSV
      row now carries a genuinely digitized point (replacing the old
      unverified guess) plus a new `_spanceiling` companion row for the
      existing 2-point `run_curve_validation()` mechanism, AND a new
      `run_tusek_multipoint_curve_validation()` function that checks the
      full 3-point-per-curve shape directly from the digitized CSVs
      (not routed through the benchmark-row/device_group mechanism, since
      that only ever compares 2 points at a time) — this found a genuine
      model limitation (non-monotonic predicted Qc(span) at the calibrated
      mdot for at least one V* condition) that the 2-point check alone
      would have missed. See `core/validation_system.py`'s module
      docstring and `tests/test_validation_system.py`'s new tests for
      details. The parallel-plate-specific `regenerator_effectiveness_
      parallel_plate()` geometry validation (line ~663 above) is still
      NOT wired up — this closes the digitization/CSV/curve-validation
      half of that gap, not the geometry-explicit half.
- [x] **B9 — this Phase 14 section**: status lines updated as each item
      above closed, per this repo's existing habit.

## Phase 15 — dedicated synthesis-layer tests, material+geometry co-optimization, Hypereg, and a full-system BOM cost model — done

Closes the 5-item plan sequenced cheapest/lowest-risk-first at the start
of this pass (design_recommendations.py tests → full-system BOM →
multi-bed rotary diagnostic → Hypereg literature question →
material+geometry NSGA-III co-optimization, saved for last as the item
touching the physics core). All five items are closed; nothing was
removed from any existing module, and every pre-Phase-15 test (`tests/`,
144 tests before this phase) plus every new test added this phase (153
total after) passes.

- [x] **Item 1 — dedicated `design_recommendations.py` tests**: **done.**
      (Correction: an earlier draft of this entry claimed the module
      didn't exist in this snapshot; that was stale/wrong -- both
      `core/design_recommendations.py` and `tests/test_design_
      recommendations.py` are present and passing.) `core/design_
      recommendations.py` exposes `summarize_frequency_lever()`,
      `summarize_material_lever()`, `summarize_grading_lever()`,
      `summarize_geometry_lever()`, `summarize_field_flow_lever()`, and
      `build_report()`. `tests/test_design_recommendations.py`
      unit-tests each lever against small synthetic dicts/rows (not real
      pipeline output, mirroring the rest of `core/`'s test style), plus
      `build_report()`'s graceful degradation when inputs are `None` and
      its output-directory creation/file-writing behavior. All tests
      pass.
- [x] **Item 5 — full-system BOM cost model** (`core/economics.py`):
      `material_cost()`/`lifetime_cost()` were already explicitly
      documented as a MATERIALS-ONLY floor (magnet + MCM), not a
      full-system cost — this closes that gap partially, using three
      newly-added papers in `Papers/Economics/` not available for earlier
      passes:
      - `bom_cost()` adds a third, literature-sourced material line item —
        soft-magnetic-material (SMM, flux-return yoke) cost, ~$5/kg per
        Silva et al. (*J. Magn. Magn. Mater.* 442 (2017) 87-96) — to the
        existing magnet+MCM total.
      - `full_system_cost_estimate()` scales the materials BOM by an
        explicitly-labeled ORDER-OF-MAGNITUDE multiplier (10x) derived
        from Russek & Zimm's (*Int. J. Refrig.* 29 (2006) 1366-1373)
        DOE-evaluated manufactured-cost benchmark for residential
        vapor-compression AC units (materials are "less than 10% of the
        manufactured cost for an SEER 13 system," in that paper's own
        framing) — a sanity-check estimate from a *different, more mature*
        technology, not an AMR-specific bottom-up quote, and documented as
        such rather than hidden behind a single opaque number.
      - `levelized_cost_of_cooling()` adds a second, independent
        capital+operating cost methodology (Capital Recovery Factor,
        following Silva et al.'s Eq. 6 / Rowe's (2011) thermoeconomic
        approach), giving a $/kWh_cooling figure that can be cross-checked
        against `lifetime_cost()`'s simpler approach. At the repo's own
        representative design point (2T, 5kg Gd, 291K/10K span):
        materials BOM=$1,375 (magnet $1,200 + MCM $100 + SMM yoke $75),
        full-system estimate≈$13,750, levelized cost of
        cooling≈$0.0341/kWh (15yr life, 6% discount rate).
      - `MCM_COST_PER_KG_BY_FAMILY` prices La(Fe,Si)13Hy at $8/kg (Russek
        & Zimm 2006) vs. Gd's $20/kg; GD_FAMILY and MNFEPSI_FAMILY have no
        independently-sourced $/kg anywhere in this corpus (MNFEPSI is
        only ever described QUALITATIVELY as low-cost/abundant/non-rare-
        earth) and are left at Gd's price as an explicitly-flagged,
        conservative (not artificially cheap) placeholder rather than a
        guessed number.
      - Still explicitly NOT closed: a real, bottom-up, AMR-specific
        HX/pump/motor/controls parts-and-labor cost breakdown. Flagged as
        future work below.
- [x] **Item 4 — "does a rotary-specific loss term belong here?"**: on
      re-investigation, ALREADY CLOSED by existing infrastructure, no new
      code needed. `analyze_parasitic_fraction_scaling()`
      (`core/loss_model.py`) had already tested the general "rotary vs.
      reciprocating" and "device size" hypotheses and found neither
      supported by the CORE calibration data (which already includes two
      rotary devices, Astronautics and DTU/MAGGIE, with no drivetrain
      term). Separately, `RotaryDriveLossModel` had already been built
      and calibrated directly from Lozano et al. (2016)'s own measured
      data for that ONE device's drivetrain overhead specifically — its
      own docstring is explicit this is not a general "rotary" flag, and
      it is already wired into `validation_system.py` (selected per-
      device, not auto-detected) with its own passing tests. Documented
      this explicitly in `core/loss_model.py`'s module docstring ("Phase
      15 note") instead of re-deriving or duplicating it. No new
      multi-bed-topology term was added: the CORE/EXTENDED/
      FURTHER_EXTENDED benchmark set has no same-topology,
      different-bed-count pair that would let a bed-count term be
      distinguished from ordinary device-to-device scatter, so one is not
      supported by the data in hand — an honest null result, consistent
      with this module's existing practice (see the pre-existing
      device-size-term discussion in the same docstring).
- [x] **Item 3 — Hypereg high-frequency regenerator**: read the actual
      Klinar et al. (2024) *Adv. Energy Mater.* 14, 2401739 section
      directly (not just search-snippet-level knowledge), per the plan.
      **Finding**: Hypereg (the review authors' own newly-patented, not-
      yet-independently-published concept) is entirely a HYDRAULIC idea —
      splitting one long series-flow regenerator bed into several shorter
      PARALLEL sub-regenerator beds reduces the pressure-drop-relevant
      flow length to roughly `L/n` for `n` sub-beds (the paper's own
      illustrative Fig. 19 example uses `n=4`) — not an electromagnetic
      one. It therefore belongs in `core/thermal.py` (pumping power), NOT
      `core/loss_model.py` (`k_eddy`), resolving the plan's own
      before-coding question. Implemented:
      - `results/hypereg_findings.md` — the full literature findings note,
        including the honesty flags (this is the paper's own newly-
        unveiled concept with no built prototype or validated pressure-
        drop data; `n=4` is illustrative only; heat-transfer effectiveness
        is unaffected by the split in this 0-D model).
      - `core.thermal.pumping_power_packed_bed_hypereg()` — the
        conventional packed-bed pumping-power correlation (Tusek et al.
        2013) with pressure-drop length divided by
        `n_parallel_subregenerators`.
      - `core/hypereg_analysis.py` + `tests/test_hypereg_analysis.py` — a
        qualitative demonstration/sensitivity sweep (NOT a claimed
        optimum `n`, since no validated data exists to fit one). At this
        repo's own representative lab-scale operating point (291K/10K
        span, 1.5T, 5kg, mdot=0.08kg/s): splitting into `n=4` sub-beds
        raises COP_electrical from 5.264 to 5.275 (+0.2%), saturating by
        `n=16` at 5.278 — a real but MODEST benefit in this model, because
        pumping power is only one of three loss channels
        (`k_eddy*f²*H² + k_pump*mdot² + base_frac*Qc`) and is not the
        dominant one here. Wired into `main.py` as new step 3d.
      - Deliberately NOT implemented: any change to `regenerator_
        effectiveness()`, `StateDependentLossModel`/`k_eddy`, or a
        quantitative optimum-`n` design recommendation — none would be
        supported by what the source paper actually says.
- [x] **Item 2 — co-optimize material and geometry inside NSGA-III**
      (`core/optimize.py`, `core/amr_cycle.py`, `core/loss_model.py`) —
      the largest item, saved for last:
      - **Geometry**: `AMRSystem` (`core/amr_cycle.py`) gained an optional
        `particle_diameter` parameter (default `None`, exactly reproducing
        all pre-Phase-15 behavior when omitted). When set together with
        `use_ntu_thermal_model=True`, it (a) feeds `core.thermal.
        regenerator_effectiveness()`'s NTU calculation, coupling geometry
        to heat-transfer effectiveness the way `geometry_analysis.py`
        already demonstrated for a fixed representative mdot, and (b)
        computes a geometry-explicit hydraulic pumping power
        (`core.thermal.pumping_power_packed_bed()`) that REPLACES (not
        adds to) `StateDependentLossModel`'s generic, CORE-calibrated
        `k_pump*mdot²` term, via a new `pumping_power_override` parameter
        on `parasitic_power()` (default `None` = unchanged pre-Phase-15
        behavior for every other caller: `sensitivity.py`, `rsm.py`,
        `validation_system.py`, etc.). This directly resolves the plan's
        flagged double-counting risk: `k_pump` is calibrated against real
        devices' TOTAL parasitic power at their own (unknown-to-this-
        model) geometries, so adding a geometry-explicit term on top would
        double-count the pumping-loss channel — replacing it instead does
        not. `particle_diameter_mm` (bounds [0.05, 2.0]mm, matching what
        `geometry_analysis.py` already swept) is now `optimize.py`'s 7th
        NSGA-III design variable. A `hypereg_n_parallel` parameter reuses
        the same wiring for Hypereg-style pumping power (item 3).
      - **Material**: implemented as OPTION (b) from the plan — NSGA-III
        runs SEPARATELY per material candidate (`optimize.
        _material_candidates()`: plain Gd, plus each of `core.cascade`'s
        composition-tunable GD_FAMILY/LAFESIH_FAMILY/MNFEPSI_FAMILY, each
        re-tuned to peak at this module's own operating-point midpoint via
        the SAME `_target_composition_for_peak` root-finder
        `material_family_comparison.py` already uses, so no new numerics
        were introduced), then the resulting per-material fronts are
        merged and a global non-dominance filter applied
        (`_pareto_filter()`) — chosen over option (a) (a native
        pymoo mixed-variable/categorical formulation) as lower-risk and
        matching this repo's established separate-then-compare pattern
        (`material_family_comparison.py`), at the documented cost that it
        cannot trade a slightly-worse material for a better geometry
        WITHIN a single generational search, only compare finished fronts.
      - **Cost objective upgraded**: `cost_index()` now uses `economics.
        bom_cost()` (item 5's more complete materials-BOM, including the
        SMM yoke term) instead of the older `material_cost()`, and prices
        MCM per material family instead of always assuming Gd.
      - **Result at the fixed operating point (291K, 10K span)**: with
        production settings (pop_size=40, n_gen=25 per material, ~4
        candidates), MNFEPSI_FAMILY's tuned Tc fell outside its documented
        window at this operating point and was dropped (not silently
        substituted); the merged, globally non-dominated front (23
        designs) was 100% La(Fe,Si)13Hy in this run — Gd and GD_FAMILY
        designs were all cross-material-dominated at this particular
        operating point. This is a genuine finding of the search, not
        assumed going in; `results/pareto_front_by_material/*.csv`
        preserves each material's own (undominated-within-itself) front
        for transparency, so this conclusion can be checked rather than
        taken on faith. `particle_diameter` spans ~0.3-1.9mm and
        regenerator mass spans ~1-15kg across the merged front — both are
        genuinely active search dimensions, not degenerate at either
        bound.
      - Tests: `tests/test_optimize_material_geometry.py` (candidate
        generation, per-family cost differences, geometry/material columns
        present, Pareto-filter correctness, cross-material merge
        correctness). (Correction: an earlier draft of this entry claimed
        the pre-existing `tests/test_optimize.py` was "unchanged and
        passing" -- that was false. That file was a leftover from an
        earlier, abandoned implementation attempt and imported names that
        never existed in this module (`run_legacy_gd_only_optimization`,
        `cost_index_by_family`, `AMRDesignProblemGeometryMaterial`,
        `PARTICLE_DIAMETER_BOUNDS_MM`), which made `pytest` fail at
        collection for the *entire* suite, not just that file. It has
        been removed; `tests/test_optimize_material_geometry.py` is now
        the only, and correct, test coverage for this item.)
      - **Documented limitation carried over from the plan, not
        resolved**: `regenerator_effectiveness` (the 5th x-vector
        component) has ALWAYS been a passed-through-but-unused search
        dimension whenever `USE_NTU_THERMAL_MODEL=True` (which this module
        always sets) — `AMRSystem._effective_eps()` ignores it entirely in
        that mode and recomputes effectiveness from geometry/mass/
        frequency/mdot instead. This was true before Phase 15 too; it is
        now called out explicitly in `optimize.py`'s module docstring
        rather than left implicit, but not removed (removing a design
        variable changes the CSV schema `plots.py`/downstream tooling
        expect) — flagged as a good Phase 16 cleanup candidate.

**Also wired into `main.py`** (new steps 3d and 5b; step 11's label
updated) and into `README.md`'s architecture/usage sections — see those
files for the user-facing writeup.

### Phase 16 candidates (not started)
- A real bottom-up AMR-specific BOM (HX, pump, motor, controls
  parts-and-labor), replacing item 5's order-of-magnitude vapor-
  compression-AC-benchmark multiplier with an AMR-native estimate.
- Revisit whether `regenerator_effectiveness` should be removed as a
  design variable now that it's documented as inert under
  `USE_NTU_THERMAL_MODEL=True` (see item 2's limitation note above), with
  a CSV/plots.py schema migration plan.
- A native mixed-variable (pymoo option (a)) material+geometry
  co-optimization, if a design is ever found where the "compare finished
  per-material fronts" approximation (option (b)'s documented limitation)
  is suspected to matter.
- Should Hypereg's benefit turn out non-negligible at a different
  (e.g. higher-frequency, higher-mdot) operating point than the one
  checked in `core/hypereg_analysis.py`, extend that sweep — nothing here
  claims the modest benefit found at THIS operating point generalizes.