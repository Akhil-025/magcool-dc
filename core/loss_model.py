"""
loss_model.py
=============
Implements a state-dependent parasitic loss model for the AMR system,
replacing the constant parasitic-power assumption with losses that depend
on operating conditions. The model was motivated by sensitivity analysis,
which showed that a constant parasitic term made electrical COP almost
entirely insensitive to magnetic field, frequency, and fluid flow rate.

Loss mechanisms modeled (functional forms are standard scaling laws; the
three coefficients are fitted, not derived from first principles):

    W_eddy  = k_eddy * f² * (mu0H)²
        Eddy-current losses in the magnet/regenerator support structure
        scale with (dB/dt)² and therefore approximately with f²H²
        (Kitanovski et al., 2015, Ch. 6).

    W_pump  = k_pump * mdot²
        Pumping power for porous-media flow through the regenerator.
        Assuming Darcy-flow behaviour, ΔP ∝ mdot and therefore
        pump power ∝ mdot².

    W_base  = base_frac * Qc
        Baseline electrical overhead (controls, inverter, bearings),
        assumed proportional to cooling duty.

    W_parasitic = W_eddy + W_pump + W_base

Calibration
-----------
The coefficients (k_eddy, k_pump, base_frac) are obtained by fitting a
3-parameter linear model to experimentally reported AMR devices with
published operating conditions and calibrated flow rates. The CORE set
(3 devices, 3 unknowns) is exactly determined.

Phase 7 update: the fit is solved with non-negative least squares (NNLS,
Lawson & Hanson 1974, `scipy.optimize.nnls`) rather than unconstrained
least squares with post-hoc clipping of negative coefficients to zero.
NNLS is the correct tool once the physical constraint is "loss
coefficients cannot be negative": it solves the *constrained* problem
directly (the best non-negative fit), whereas clipping an unconstrained
solution afterwards does not minimize any well-defined objective under
that constraint and can silently mask an unstable fit. For the CORE
3-point set the unconstrained optimum is already non-negative, so NNLS
and the old lstsq-then-clip approach agree exactly (checked by
`tests/test_loss_model.py`). For the EXTENDED 4-point set, NNLS
removes the negative k_eddy/base_frac Phase 6 reported (`k_eddy=0` is
now the constrained optimum instead of a clipped negative value), but
leave-one-out error predicting the smallest device (Tušek, 5.3 W --
Phase 31: Qc/Wp corrected from an old, unverified 6.5W/0.76W guessed
point to a genuinely digitized one, see CALIBRATION_POINTS_CORE below)
from the other three is now ~333% (down from ~680% under the old,
pre-Phase-31 Tušek point, and from +1639% with plain lstsq before that,
but still an order-of-magnitude miss) -- confirming Phase 6's conclusion
that pooling four orders of magnitude of device scale needs a structural
change to the model, not just a better-behaved solver or better-sourced
calibration data.

`analyze_parasitic_fraction_scaling()` checks the natural next
hypothesis -- that this is a simple device-size effect, i.e. that
smaller devices carry proportionally more fixed (non-Qc-scaling)
overhead. Sorting the four devices by Qc shows this does NOT hold in
the fixed-overhead/economies-of-scale direction: the smallest device
(Tušek, 5.3 W -- Phase 31: updated from the old 6.5W point, see above)
has the *lowest* parasitic fraction (13.8%), and the
largest (Astronautics, 2502 W) has the *highest* (45.3%) -- the
opposite of what a fixed-overhead story predicts (which would put
small devices at the top, not the bottom). Paper-Mining Pass Part 6
correction: with the DTU point corrected from the old fabricated
818 W/17.1% figure to the verified 102.8 W/25.5% figure (see this
module's correction note above), the four-device ranking by Qc
(Tušek 5.3 W/13.8%, DTU 102.8 W/25.5%, Okamura 200 W/36.7%,
Astronautics 2502 W/45.3%) is now cleanly MONOTONICALLY INCREASING
with device scale -- the opposite trend from a fixed-overhead story,
and no longer "non-monotonic in between" as the old fabricated DTU
figure had made it appear. This strengthens rather than reverses the
original conclusion: a fixed-overhead/economies-of-scale size term is
still not what's happening (that would predict the fraction falling
with scale, not rising), it just turns out the actual relationship
between the four verified devices is a clean trend in the opposite
direction rather than a scattered non-pattern. The Astronautics figure
is independently flagged by its own source paper as reflecting
"mediocre" electrical-component efficiency at that scale (a
device-specific engineering choice), not a generic size law, and with
only 4 points spanning orders of magnitude in scale and design era,
a genuinely monotonic trend by itself is not strong evidence of a real
physical size law either. So a size/scale *term* -- as the Phase 6
write-up speculated -- is still not adopted here; what varies
device-to-device looks more like motor/inverter efficiency class and
drivetrain topology (rotary vs. reciprocating single-bed) than raw
cooling capacity. Flagged as a
correction to the Phase 7 roadmap: more benchmark devices with
independently reported component efficiencies, not a size term, is the
concrete next step.

A fifth candidate dataset (Risø/DTU 2011, 30 K span) could not be
calibrated because the corresponding operating point does not yield
positive cooling capacity under the present AMR model.

Phase 15 note -- "does a rotary-specific loss term belong here?" (this
question was re-raised going into Phase 15; re-checking the existing
analysis below shows it was already answered, no new code needed): YES,
in the specific, data-supported sense already implemented as
`RotaryDriveLossModel` below, and NO in the broader "flag every rotary
device" sense the question could also be read as. `analyze_parasitic_
fraction_scaling()` already tested the more general hypothesis (device
SIZE, not topology) and found no size term supported by the data (see
above). Separately, and specifically, `fit_rotary_drive_term()` /
`RotaryDriveLossModel` already found and calibrated a genuine additional
loss channel (mechanical drivetrain power to move a rotary permanent-
magnet assembly + rotary valve) directly from Lozano et al. (2016)'s own
measured WM data -- but its own docstring is explicit that this is NOT a
general "rotary AMR" flag: Astronautics and DTU are ALSO rotary devices
and are already well-represented by the CORE calibration without any
drivetrain term, because what distinguishes Lozano is that device's own
drivetrain overhead relative to its own small Qc (an early-generation,
unoptimized prototype the source paper itself calls "modest ... in
comparison with established cooling technologies"), not rotary operation
per se. So: the CORE 3-point set (which already includes two rotary
devices, Astronautics and DTU/MAGGIE) remains the production default;
`RotaryDriveLossModel` remains an available, data-justified PER-DEVICE
option for devices in Lozano's drivetrain-dominated class, selected
explicitly by the caller (as `validation_system.py` already does for the
Lozano rows) rather than auto-detected from a "rotary" flag. No new
multi-bed-topology term is added this pass: the CORE/EXTENDED/
FURTHER_EXTENDED benchmark set does not contain a same-topology,
different-bed-count pair that would let a bed-count term be distinguished
from ordinary device-to-device scatter, so inventing one would not be
supported by data in hand -- consistent with this module's existing
practice of stating an honest null result (see the size-term discussion
above) rather than adding unsupported structure.

Paper-Mining Pass Part 6 correction: the CORE point previously labeled
"DTU_rotary_Gd_2016" (818 W, 10.1 K span, COP=4.2, cited only as
"Bahl/Eriksen/Engelbrecht, rotary AMR - ScienceDirect (2016)") was never
actually located in this repo's Papers/ and has now been checked directly
against the real paper behind that citation -- D. Eriksen, K. Engelbrecht,
C.R.H. Bahl, R. Bjørk, "Exploring the efficiency potential for an active
magnetic regenerator," Sci. Technol. Built Environ. 22(5) (2016) 527-533
(reproduced as Chapter 6 of Eriksen's 2016 DTU PhD thesis, now in this
repo's Papers/). The real paper's own headline result is 81.5 W at a
15.5 K span with COP=3.6 (18% second-law efficiency) at fAMR=0.61 Hz,
1.13 T, 1.7 kg Gd -- not 818 W/10.1 K/COP=4.2. See
data/amr_experimental_benchmarks.csv's DTU_Eriksen_MAGGIE_2016 row for
the full correction note, including a directly-measured loss breakdown
(shaft power, pump power split, Carnot work) taken straight from the
thesis's Table 6.2 -- no back-calculation needed for those numbers.

That corrected point does NOT calibrate under this module's cycle model,
however: at 1.13 T / 1.7 kg / 0.61 Hz, amr_cycle.py's cooling_capacity()
predicts Qc ≈ 0 at a 15.5 K span for any mdot in [1e-6, 5] kg/s -- the
model's own zero-flow no-load span at this field/frequency already sits
below 15.5 K, the same failure mode already documented above for
Risø/DTU 2011. This is attributed to the real device's Curie-graded
11-layer regenerator (Gd + three Gd(100-x)Yx alloys) reaching spans a
single-uniform-Tc Gd approximation structurally cannot reach, not a
data error -- the row is kept in the CSV as a documented
non-calibrating point, not silently dropped.

The 3rd CORE slot vacated by this correction is now filled by
DTU_Eriksen_rotary_Gd_2015 (10.2 K span, 102.8 W, COP=3.1 at 0.75 Hz --
the SAME physical prototype, "MAGGIE", at an earlier paper's
lower-span operating point), which was already in this repo's benchmark
set with a verified primary citation and DOES calibrate cleanly under
the model. The parasitic-fraction comparison two paragraphs above (Astronautics
45.3%, Okamura 36.7%, old-DTU 17.1%, Tušek 11.8%) is retained as
historical context but the "old-DTU 17.1%" figure describes the now-
retracted 818 W point and should not be treated as current.

Follow-up (this session, after amr_cycle.py gained no_load_span_override):
DTU_Eriksen_MAGGIE_2016's own 81.5W/15.5K/COP=3.6 point, the one just
described as non-calibrating, DOES calibrate once the override supplies a
real (simulated, not assumed) span cap above 15.5K -- see
CALIBRATION_POINTS_CORE_PLUS_MAGGIE_HIGHSPAN below and
run_core_plus_maggie_highspan_diagnostic(). Honest leave-one-out result:
this 4th point -- unlike EXTENDED/FURTHER_EXTENDED's larger-scale-jump
points, which do NOT generalize -- is the SAME hardware as an existing
CORE point at a different condition, and its own held-out fold predicts
within +30% (vs. the ~250-700% every other fold, old or new, still
shows), while modestly IMPROVING (not degrading) the two folds it shares
with CORE. Real, partial progress on the "which kind of additional data
would actually help" question this module's docstring poses above --
still not enough to promote past CORE as the production default, since
Astronautics and Tušek's own folds remain order-of-magnitude misses
either way.
"""

import numpy as np
from scipy.optimize import nnls

# (device, f_Hz, mu0H_T, mdot_kg_s, Qc_W, W_parasitic_required_W)
# CORE: the stable, well-behaved 3-point set used as the production default.
#
# FIX (Paper-Mining Pass Part 4): all three mdot values below were stale.
# Verified by directly checking whether they reproduce their own literature
# Qc under amr_cycle.py's CURRENT cooling_capacity() -- none did:
#   Astronautics: old mdot=1.0854 -> predicted 10734W, not 2502W (4.3x off)
#   DTU:          old mdot=0.3251 -> predicted 1343W,  not  818W (1.64x off)
#   Tusek:        old mdot=0.0045 -> predicted   12W,  not  6.5W (1.86x off)
# i.e. the cooling-capacity model was changed at some point after these
# were hardcoded (most likely the NTU thermal-model work) and this
# calibration set was never re-synced -- meaning every downstream user of
# the default StateDependentLossModel() (Sobol 9b, the RSM surrogate,
# NSGA-III) had been running against a loss model fit to inputs the cycle
# model itself could no longer reproduce. Re-calibrated below via the same
# brentq(qc_residual, 1e-6, 5.0) procedure validation_system.py uses, with
# GADOLINIUM material and T_cold=289K for all three (matching the original
# CORE set's convention -- this is a deliberately-separate, Gd-only
# approximation from validation_system.py's newer LAFESIH_FIRST_ORDER-based
# system-validation row for the same physical Astronautics device; see that
# module's docstring for why the two are intentionally different).
# W_parasitic_required recomputed at the same time (Wp = Qc*(1/COP_lit -
# 1/COP_ideal)); interestingly these barely moved even though mdot did
# (COP_ideal turns out to be fairly stable across the corrected mdot
# range), so the NNLS fit's target vector was mostly fine -- only the
# mdot^2 pumping-term *inputs* (the A matrix) were wrong.
# `tests/test_loss_model.py::test_core_calibration_points_are_self_consistent`
# guards against this drifting again silently.
CALIBRATION_POINTS_CORE = [
    # CITATION AUDIT (Phase 30, CITATION_AUDIT_PHASE30.md items 1-2):
    # Astronautics_rotary_2014 and DTU_Eriksen_rotary_Gd_2015 below were
    # independently re-verified this pass directly against their PDF
    # source text (Jacobs et al. and Eriksen et al. respectively, both
    # present in Papers/) -- field, frequency, Qc, and COP all match the
    # primary source word-for-word for both points. Tusek_singlebed_Gd_2010
    # (below) was ALSO independently checked and its existing field/
    # frequency discrepancy flag (see that point's own long comment) was
    # CONFIRMED, not newly discovered: the source paper's own text reports
    # 1.15T/0.3Hz throughout, "1.69" does not appear anywhere in it.
    ("Astronautics_rotary_2014", 4.0, 1.44, 0.252999, 2502.0, 1133.70),
    # FIX (Paper-Mining Pass Part 6): replaces the fabricated/unlocated
    # "DTU_rotary_Gd_2016" point (818W, mdot=0.198062, Wp=139.79 -- these
    # numbers do not correspond to any real published operating point, see
    # this module's docstring and data/amr_experimental_benchmarks.csv's
    # DTU_Eriksen_MAGGIE_2016 row for the correction). DTU_Eriksen_rotary_Gd_2015
    # is the same physical prototype ("MAGGIE") at a genuinely verified,
    # primary-sourced operating point (Eriksen et al., Int. J. Refrigeration
    # 2015) that calibrates cleanly under the current cooling_capacity()
    # model. mdot and Wp_required recomputed with the same brentq(
    # qc_residual, 1e-6, 5.0)/Wp=Qc*(1/COP_lit-1/COP_ideal) procedure used
    # for every other CORE point, GADOLINIUM material, T_cold=289K:
    #   mdot_cal = 0.084666 kg/s  ->  Qc_model = 102.8W (exact match)
    #   COP_ideal = 14.73  ->  Wp_required = 102.8*(1/3.1 - 1/14.73) = 26.18W
    ("DTU_Eriksen_rotary_Gd_2015", 0.75, 1.13, 0.084666, 102.8, 26.18),
    # Tusek: FIXED (Phase 31). Previously used the pre-correction field/mass/
    # frequency (1.69T, 0.196kg, 0.25Hz) as a deliberate stopgap, because the
    # corrected, paper-verified operating point (1.15T, 0.1763kg, 0.3Hz) did
    # not calibrate at all against the old guessed (span=15K, Qc=6.5W) pair
    # -- see this module's git history / the long comment this replaces for
    # the full original diagnosis. That guessed pair has since been replaced
    # in data/amr_experimental_benchmarks.csv with a genuinely pixel-
    # digitized point from the source paper's own Figs. 10-11 (AMR-A,
    # V*=0.95: span=7.26K, Qc=5.27W, COP=5.38 -- see
    # data/tusek_ate2013_figs/fig10_data.csv, fig11_data.csv, and
    # data/tusek_ate2013_figs/tusek_ate2013_figs_notes.md for the full
    # pixel-calibration/uncertainty methodology, +-0.15K span / +-0.15W Qc).
    # That digitized point DOES calibrate cleanly at the correct 1.15T field
    # -- recomputed here via the same brentq(qc_residual, 1e-6, 5.0)/
    # Wp=Qc*(1/COP_lit-1/COP_ideal) procedure used for every other CORE
    # point, GADOLINIUM material, T_cold=289K:
    #   mdot_cal = 0.007351 kg/s  ->  Qc_model = 5.270W (exact match)
    #   COP_ideal = 20.70  ->  Wp_required = 5.27*(1/5.38 - 1/20.70) = 0.725W
    # (mdot_cal=0.0074 independently cross-checked against the CSV row's own
    # note, which reports the same value to 2 sig figs from a separate
    # computation -- good agreement.) This CORE point now uses the SAME
    # verified field/mass/frequency as data/amr_experimental_benchmarks.csv's
    # Tusek_singlebed_Gd_2010 row, closing the discrepancy
    # CITATION_AUDIT_PHASE30.md flagged between this module and that CSV.
    ("Tusek_singlebed_Gd_2010", 0.3, 1.15, 0.007351, 5.27, 0.7250),
]
# EXTENDED: CORE + Okamura & Hirano (2013). Retained only as a
# diagnostic comparison; not used as the default calibration.
CALIBRATION_POINTS_EXTENDED = CALIBRATION_POINTS_CORE + [
    # frequency not reported in the secondary source for this device --
    # 1.0 Hz placeholder (see data/amr_experimental_benchmarks.csv note)
    ("Okamura_Hirano_2013", 1.0, 1.1, 0.0502, 200.0, 0.367 * 200.0),
]
# FURTHER_EXTENDED: EXTENDED + 4 points from Lozano et al. (2016), the
# POLO/UFSC rotary device (see data/amr_experimental_benchmarks.csv and
# validation_system.py). Of the 8 reported (frequency, flow, span, Qc, COP)
# operating points, only these 4 calibrate at all within mdot in [1e-6,5]
# kg/s -- the other 4 (r1, r2, r3, r5) reported a Qc the model cannot reach
# at that span/field/mass under any flow rate, and are reported as
# calibration failures, not silently dropped (see
# validation_system.run_system_validation() output). Each Wp_required here
# comes from that device's OWN calibrated mdot: Wp = Qc * (1/COP_lit -
# 1/COP_ideal), same formula used to build CORE/EXTENDED above.
#
# CITATION AUDIT FLAG (Phase 30, CITATION_AUDIT_PHASE30.md item 4 --
# UPDATED/RESOLVED by a follow-up web search pass, see item 4a below):
# the "2016" citation itself is CORRECT. Lozano, Capovilla, Trevizoli,
# Bahl, Engelbrecht, Nielsen, Barbosa, "Development of a novel rotary
# magnetic refrigerator", Int. J. Refrig. 68 (2016) 187-197, is a real,
# genuinely 2016-dated paper (confirmed via DTU's own publication
# database and ScienceDirect) -- a DIFFERENT paper from the Lozano et al.
# Int. J. Refrig. 37 (2014) 92-98 paper that happened to be the only
# Lozano PDF in this project's local Papers/ corpus. The earlier version
# of this flag incorrectly checked the 2016 numbers against the wrong
# (2014) paper and found no match, which is expected -- they are
# different devices. The 2016 paper's own reported figures (max zero-load
# span 12K at 1Hz/150L/h; max zero-span cooling power 150W at 0.8Hz/
# 200L/h; 80.4W thermal load producing a 7.1K span; ~1T field; ~1.7kg Gd
# spheres in 8 regenerator-bed pairs) are consistent in device class,
# field, and order of magnitude with the four (Qc, COP) pairs used below,
# though the exact per-point r4/r6/r7/r8 table values were not
# individually re-verified against the paper's own data tables in this
# pass (the full PDF was not accessible to fetch directly). Net effect:
# the citation is real and correctly dated; only the exact per-row digits
# remain formally unverified, a much smaller gap than "wrong year."
CALIBRATION_POINTS_FURTHER_EXTENDED = CALIBRATION_POINTS_EXTENDED + [
    ("Lozano_POLO_UFSC_2016_r4", 0.4, 0.88, 0.3244, 62.5, 1.684 * 62.5),
    ("Lozano_POLO_UFSC_2016_r6", 0.8, 0.88, 0.0430, 81.2, 1.505 * 81.2),
    ("Lozano_POLO_UFSC_2016_r7", 0.4, 0.88, 0.0207, 80.8, 1.291 * 80.8),
    ("Lozano_POLO_UFSC_2016_r8", 0.8, 0.88, 0.0308, 120.4, 1.180 * 120.4),
]
# CORE_PLUS_MAGGIE_HIGHSPAN: CORE + DTU_Eriksen_MAGGIE_2016 (81.5W, 15.5K span,
# COP=3.6 -- the SAME physical prototype as the CORE set's own
# DTU_Eriksen_rotary_Gd_2015 point, 102.8W, at a different, later-paper
# operating point; see that CSV row's own note). Unlike EXTENDED/
# FURTHER_EXTENDED (which add devices spanning orders of magnitude in scale
# from CORE -- exactly what this module's docstring already found doesn't
# generalize, needing "a structural change to the model, not just more
# points"), this point is the SAME hardware at a COMPARABLE scale to an
# existing CORE point -- the cleanest possible test of whether a genuine
# 4th, non-exactly-determined point (4 points, 3 unknowns, an actual
# regression for once) helps or not.
#
# THIS POINT COULD NOT BE CALIBRATED until this session's amr_cycle.py
# no_load_span_override addition: at 1.13T/1.7kg/0.61Hz, cooling_capacity()'s
# default 2*dTad_noload cap sits below the reported 15.5K span for ANY mdot
# (documented in this CSV row's own note and, before that fix, in this
# module's docstring above) -- a single-uniform-Tc Gd approximation
# structurally cannot reach a span the real Curie-graded 11-layer bed
# reaches. With no_load_span_override=21.04K (core.regenerator_1d.
# regenerative_span_cap() for this device's field/mass/frequency -- see
# results/regenerative_amplification_override_check.txt), it calibrates
# cleanly: mdot_cal=0.014650 kg/s reproduces Qc=81.5W exactly, giving
# COP_ideal=9.695 -> Wp_required=81.5*(1/3.6 - 1/9.695)=14.23W. This
# dependency is a property of HOW this point's own (mdot, Wp) numbers were
# derived, not of the loss-model fit itself: once fixed as (f, H, mdot, Qc,
# Wp) numbers below, calibrate_loss_coefficients()/leave_one_out_cv() use
# them exactly like any other point, with no override needed at fit time.
# A directly-measured loss breakdown is ALSO available for this device
# (Eriksen 2016 PhD thesis Table 6.2: shaft power 14.0W including bearing/
# valve friction, pumping power 8.9W, total input 22.9W -- see the CSV
# row's note) but is NOT used here: it isn't obvious how to map "shaft"
# vs. "pump" boundaries in the real drivetrain onto this model's three
# specific terms (eddy/pump/base) without guessing, so the same
# back-calculated Wp=Qc*(1/COP_lit-1/COP_ideal) formula every other CORE
# point already uses is applied instead, for an apples-to-apples fit.
#
# HONESTY FLAG, checked directly: the device's own DIRECTLY-MEASURED flow
# rate (2.5 L/min, from the same thesis Table 6.2 cited above) is
# 0.04167 kg/s -- 2.84x LARGER than this point's back-calculated
# mdot=0.014650 kg/s. Run at the real measured flow rate instead of the
# back-calculated one, cooling_capacity() predicts Qc=231.8W against the
# actual measured Qc=81.5W -- a ~2.8x OVER-prediction. This means this
# calibration point's mdot is doing exactly what every other CORE point's
# mdot also does (solved to reproduce the reported Qc under this model,
# not taken from the device's real flow rate) and should not be read as
# "the model reproduces this device's real operating point" -- only as
# "there exists SOME flow rate under which this model's Qc(span) formula
# passes through the reported (span, Qc)." The same caveat quietly applies
# to every other point in CALIBRATION_POINTS_CORE (none of their mdot
# values are independently verified against a reported flow rate either,
# except Lozano's WM/frequency data used by RotaryDriveLossModel below);
# this is simply the first point where a real measured flow rate happened
# to be available to check against. Not treated as disqualifying -- the
# alternative (using the real 0.04167 kg/s and accepting a ~2.8x Qc miss
# instead of an exact match) would be a strictly worse calibration input,
# not a more honest one -- but it is a genuine, quantified gap between
# "this model's own effective flow dependence" and "this device's real
# flow dependence" worth having on record for whoever next revisits how
# mdot is calibrated in this module.
CALIBRATION_POINTS_CORE_PLUS_MAGGIE_HIGHSPAN = CALIBRATION_POINTS_CORE + [
    ("DTU_Eriksen_MAGGIE_2016", 0.61, 1.13, 0.014650, 81.5, 14.23),
]
CALIBRATION_POINTS = CALIBRATION_POINTS_CORE  # backward-compat alias


def _build_system(points):
    A = np.zeros((len(points), 3))
    b = np.zeros(len(points))
    for i, (_name, f, H, mdot, Qc, Wp) in enumerate(points):
        A[i, 0] = f ** 2 * H ** 2
        A[i, 1] = mdot ** 2
        A[i, 2] = Qc
        b[i] = Wp
    return A, b


def leave_one_out_cv(points=None, verbose=True):
    """Fit on N-1 points, predict the held-out point's W_parasitic, repeat
    for each point. Reports absolute and percent error per held-out device."""
    points = points if points is not None else CALIBRATION_POINTS_CORE
    results = []
    for i in range(len(points)):
        train = [p for j, p in enumerate(points) if j != i]
        test = points[i]
        A, b = _build_system(train)
        coeffs, _resid_norm = nnls(A, b)
        k_eddy, k_pump, base_frac = coeffs
        name, f, H, mdot, Qc, Wp_true = test
        Wp_pred = k_eddy * f ** 2 * H ** 2 + k_pump * mdot ** 2 + base_frac * Qc
        err_pct = 100 * (Wp_pred - Wp_true) / Wp_true if Wp_true != 0 else float("nan")
        results.append((name, Wp_true, Wp_pred, err_pct))
        if verbose:
            print(f"  held out {name:<28} W_parasitic true={Wp_true:8.2f}W  "
                  f"predicted={Wp_pred:8.2f}W  error={err_pct:+7.1f}%")
    return results


def calibrate_loss_coefficients(points=None, verbose=True, label="CORE (production default)"):
    """Fit (k_eddy, k_pump, base_frac) via non-negative least squares (NNLS,
    Lawson & Hanson 1974). NNLS solves the constrained problem directly --
    the best-fitting non-negative coefficients -- rather than solving the
    unconstrained problem and clipping negative results afterwards, which
    does not minimize any well-defined objective under the non-negativity
    constraint. When the unconstrained optimum is already non-negative (as
    for the exactly-determined CORE 3-point set) NNLS returns the same
    solution as unconstrained least squares."""
    points = points if points is not None else CALIBRATION_POINTS_CORE
    A, b = _build_system(points)
    coeffs, resid_norm = nnls(A, b)
    k_eddy, k_pump, base_frac = coeffs
    pred = A @ coeffs
    if verbose:
        print(f"Calibrated loss-model coefficients [{label}], "
              f"{len(points)} points, 3 unknowns "
              f"({'exactly-determined' if len(points) == 3 else 'over-determined'}, "
              f"NNLS, residual norm={resid_norm:.4f}):")
        print(f"  k_eddy    = {k_eddy: .6f}  W / (Hz^2 * T^2)")
        print(f"  k_pump    = {k_pump: .6f}  W / (kg/s)^2")
        print(f"  base_frac = {base_frac: .6f}  (dimensionless, x Qc)")
        print("  Fit residuals (predicted - required W_parasitic):")
        for (name, *_, Wp_true), Wp_pred in zip(points, pred):
            print(f"    {name:<28} true={Wp_true:8.2f}W  fit={Wp_pred:8.2f}W  "
                  f"resid={Wp_pred - Wp_true:+7.2f}W")
        for c, name in zip(coeffs, ["k_eddy", "k_pump", "base_frac"]):
            if c == 0.0:
                print(f"  NOTE: {name} pinned to 0 by the non-negativity "
                      "constraint (the unconstrained optimum was negative here).")
    return {"k_eddy": k_eddy, "k_pump": k_pump,
            "base_frac": base_frac, "raw": coeffs}


def analyze_parasitic_fraction_scaling(points=None, verbose=True):
    """Tests the hypothesis (raised in the Phase 6 write-up) that the loss
    model's cross-device instability is a simple device-*size* effect --
    e.g. small devices carrying proportionally more fixed overhead. Sorts
    the benchmark devices by Qc and reports the parasitic fraction
    (W_parasitic / Qc) for each. Returns the sorted (name, Qc, fraction)
    list so the caller/tests can check the (non-)monotonicity claim."""
    points = points if points is not None else CALIBRATION_POINTS_EXTENDED
    rows = [(name, Qc, Wp / Qc) for (name, f, H, mdot, Qc, Wp) in points]
    rows.sort(key=lambda r: r[1])
    if verbose:
        print("Parasitic fraction (W_parasitic / Qc) sorted by device scale:")
        for name, Qc, frac in rows:
            print(f"    {name:<28} Qc={Qc:8.1f} W   fraction={frac:.3f}")
        fracs = [r[2] for r in rows]
        monotonic = all(fracs[i] <= fracs[i + 1] for i in range(len(fracs) - 1)) or \
                    all(fracs[i] >= fracs[i + 1] for i in range(len(fracs) - 1))
        print(f"  Monotonic in device scale? {monotonic}")
        if not monotonic:
            print(f"  CONCLUSION: no monotonic size trend in this {len(rows)}-point "
                  "set -- the smallest device (Tusek, 6.5W) does NOT have the "
                  "highest overhead fraction, and the largest (Astronautics, "
                  "2502W) does NOT have the lowest. A simple size/scale term is "
                  "not supported by the data in hand; the Astronautics outlier is "
                  "independently attributed by its own source paper to 'mediocre' "
                  "electrical-component efficiency at that scale, i.e. a "
                  "device-specific engineering choice, not a generic size law. "
                  "With the Lozano points included, the picture sharpens further: "
                  "the four highest fractions in the whole set (1.18-1.68) are "
                  "ALL Lozano points clustered in the low-to-mid Qc range, sitting "
                  "well above Okamura (200W, 0.367) and Astronautics (2502W, "
                  "0.453) -- i.e. grouped by device/paper, not ordered by scale.")
        else:
            print("  NOTE (Paper-Mining Pass Part 6): with the DTU point corrected "
                  "from its old fabricated 818W/0.171 figure to the verified "
                  "102.8W/0.255 figure, this 4-point EXTENDED set is now "
                  "monotonically INCREASING with Qc -- the opposite direction "
                  "from a fixed-overhead/economies-of-scale story (which "
                  "predicts fraction falling as Qc grows), not confirmation of "
                  "one. See this module's docstring for the full discussion; "
                  "a size/scale term is still not adopted as a result.")
    return rows


def run_extended_diagnostic():
    """Demonstrates the instability of the four-point extended fit.
    The diagnostic is provided for transparency and is not used as the
    production calibration."""
    print("=" * 90)
    print("DIAGNOSTIC: adding Okamura & Hirano (2013) as a fourth calibration point")
    print("=" * 90)
    calibrate_loss_coefficients(CALIBRATION_POINTS_EXTENDED, verbose=True,
                                  label="EXTENDED (4pt, diagnostic only, NNLS)")
    print("\n  Leave-one-out cross-validation on the EXTENDED set (NNLS per fold):")
    loo = leave_one_out_cv(CALIBRATION_POINTS_EXTENDED, verbose=True)
    worst = max(loo, key=lambda r: abs(r[3]))
    print(f"\n  CONCLUSION: switching from unconstrained lstsq to NNLS removes the "
          f"negative (unphysical) coefficients Phase 6 found, and improves the "
          f"worst leave-one-out error from +1639% to {worst[3]:+.0f}% "
          f"(held-out device: {worst[0]}) -- but that is still an order-of-"
          f"magnitude miss. A better-behaved solver alone does not make a "
          f"single linear model generalize across devices spanning 6.5W to "
          f"2502W of cooling capacity. The CORE 3-point fit remains the "
          f"production default.")
    print("\n  Testing the natural next hypothesis -- that this is a simple "
          "device-size effect:")
    analyze_parasitic_fraction_scaling(CALIBRATION_POINTS_EXTENDED, verbose=True)


def run_further_extended_diagnostic():
    """Adds the 4 calibratable Lozano et al. (2016) POLO/UFSC points to the
    EXTENDED set (8 points total spanning 5.3-2502W) and re-checks NNLS
    fit stability and the scale-monotonicity question. This device class
    is qualitatively different from CORE/EXTENDED: Lozano's own paper
    reports COP of 0.37-0.83 (vs. 1.9-4.6 for every other benchmark
    device) because its motor power (87-145W) is comparable to or exceeds
    its own Qc (61-120W) -- an early-generation, unoptimized rotary
    prototype the paper itself calls "modest ... in comparison with
    established cooling technologies". The implied parasitic fractions
    here (1.18-1.68) exceed 1.0, i.e. parasitic power alone exceeds Qc --
    outside the range of any CORE/EXTENDED device (max 0.453)."""
    print("=" * 90)
    print("DIAGNOSTIC: adding 4 calibratable Lozano POLO/UFSC (2016) points "
          "(8pt FURTHER_EXTENDED)")
    print("=" * 90)
    calibrate_loss_coefficients(CALIBRATION_POINTS_FURTHER_EXTENDED, verbose=True,
                                  label="FURTHER_EXTENDED (8pt, diagnostic only, NNLS)")
    print("\n  Leave-one-out cross-validation on the FURTHER_EXTENDED set (NNLS per fold):")
    loo = leave_one_out_cv(CALIBRATION_POINTS_FURTHER_EXTENDED, verbose=True)
    worst = max(loo, key=lambda r: abs(r[3]))
    lozano_loo = [r for r in loo if r[0].startswith("Lozano")]
    lozano_mean_err = np.mean([r[3] for r in lozano_loo])
    print(f"\n  CONCLUSION: the naive 'worst leave-one-out error' metric actually "
          f"IMPROVES with 8 points ({worst[3]:+.0f}%, held-out device: {worst[0]}) "
          f"versus the 4-point EXTENDED set's +682% -- but this is misleading, "
          f"not genuine progress: it happens because Tusek's small W_parasitic "
          f"becomes easier to hit once more low/mid-Qc points anchor the fit, "
          f"not because the model now generalizes to the Lozano device class. "
          f"The more informative number is that all 4 Lozano points, when held "
          f"out, are underpredicted by a similar amount (mean "
          f"{lozano_mean_err:+.0f}%, individually {[f'{r[3]:+.0f}%' for r in lozano_loo]}) "
          f"-- a CONSISTENT, one-directional miss, not noise. That consistency is "
          f"itself evidence for the standing hypothesis: Lozano's device sits in a "
          f"different, definably-worse motor/inverter efficiency class (its own "
          f"paper: COP 'modest ... in comparison with established cooling "
          f"technologies', motor power comparable to or exceeding Qc) that a "
          f"state-variable-only loss model (frequency, field, mdot, Qc) cannot "
          f"represent, regardless of solver. The CORE 3-point fit remains the "
          f"production default.")
    print("\n  Re-testing scale monotonicity with the enlarged 8-device set:")
    analyze_parasitic_fraction_scaling(CALIBRATION_POINTS_FURTHER_EXTENDED, verbose=True)


def run_core_plus_maggie_highspan_diagnostic():
    """Tests CALIBRATION_POINTS_CORE_PLUS_MAGGIE_HIGHSPAN -- CORE + a 4th
    point that is, unlike EXTENDED/FURTHER_EXTENDED, the SAME physical
    device as an existing CORE point (DTU_Eriksen_rotary_Gd_2015) at a
    different operating point, only reachable using this session's
    amr_cycle.py no_load_span_override addition (see that calibration set's
    own docstring above for the full derivation). This is the first
    genuinely over-determined (4 points, 3 unknowns) fit in this module's
    history where the added point is comparable in scale to CORE rather
    than spanning additional orders of magnitude -- the honest test of
    whether more data of the RIGHT kind (same device, different point) does
    better than more data of the kind already tried (different, bigger/
    smaller devices)."""
    print("=" * 90)
    print("DIAGNOSTIC: adding DTU_Eriksen_MAGGIE_2016 (same hardware as an existing")
    print("CORE point, different operating point) as a 4th calibration point")
    print("=" * 90)
    calibrate_loss_coefficients(CALIBRATION_POINTS_CORE_PLUS_MAGGIE_HIGHSPAN, verbose=True,
                                  label="CORE_PLUS_MAGGIE_HIGHSPAN (4pt, NNLS)")
    print("\n  Leave-one-out cross-validation (NNLS per fold):")
    loo_4pt = leave_one_out_cv(CALIBRATION_POINTS_CORE_PLUS_MAGGIE_HIGHSPAN, verbose=True)
    print("\n  Same three folds under the 3-point CORE set alone, for direct comparison:")
    loo_3pt = leave_one_out_cv(CALIBRATION_POINTS_CORE, verbose=True)
    err_3pt = {r[0]: r[3] for r in loo_3pt}
    lines = []
    for name, _true, _pred, err4 in loo_4pt:
        if name in err_3pt:
            lines.append(f"    {name:<28} CORE(3pt)={err_3pt[name]:+7.1f}%  "
                          f"CORE_PLUS_MAGGIE_HIGHSPAN(4pt)={err4:+7.1f}%")
        else:
            lines.append(f"    {name:<28} (new point, no 3pt baseline) "
                          f"CORE_PLUS_MAGGIE_HIGHSPAN(4pt)={err4:+7.1f}%")
    print("\n  Per-device comparison:")
    for line in lines:
        print(line)
    print(f"\n  CONCLUSION: unlike EXTENDED/FURTHER_EXTENDED (which add devices spanning "
          f"additional orders of magnitude and do NOT generalize -- see "
          f"run_extended_diagnostic()/run_further_extended_diagnostic() above), this point "
          f"is the SAME physical prototype as an existing CORE point at a different "
          f"operating condition, and genuinely helps: the new point's own held-out "
          f"leave-one-out error is {[r[3] for r in loo_4pt if r[0]=='DTU_Eriksen_MAGGIE_2016'][0]:+.0f}%, "
          f"far inside the ~250-700% range every other fold (old or new) still shows, and "
          f"the two pre-existing folds shared with CORE (Astronautics, "
          f"DTU_Eriksen_rotary_Gd_2015) both improve modestly rather than degrading. This is "
          f"real, partial progress -- NOT a fix for the underlying zero-degrees-of-freedom "
          f"problem (Astronautics and Tusek's own folds are still order-of-magnitude misses "
          f"with 4 points, same as with 3), consistent with this module's standing diagnosis: "
          f"what generalizes is more data of the SAME device class, not more devices per se. "
          f"CORE (3pt) remains the production default -- promoting a still-mostly-"
          f"order-of-magnitude-off fit to default status on the strength of one good fold "
          f"would overstate what this shows. Available as an explicit opt-in "
          f"(calibrate_loss_coefficients(CALIBRATION_POINTS_CORE_PLUS_MAGGIE_HIGHSPAN)) for "
          f"anyone specifically modeling this device or its operating regime.")


# Lozano et al. (2016), Table 3, WM column: electrical power to drive the
# rotary permanent-magnet assembly + rotary valve system (measured
# SEPARATELY from WP, the pump power -- see Table 3 and Eq. 2, COP =
# Qc/(WP+WM)). Digitized directly from the paper (Table 3 row 1-8),
# NOT derived/backed-out from Qc and COP, unlike CALIBRATION_POINTS_CORE's
# Wp_required column.
LOZANO_DRIVE_CALIBRATION_HZ_W = [
    # (f_Hz, WM_W)
    (1.4, 145.0), (0.8, 108.1), (0.8, 106.5), (0.4, 89.5),
    (0.8, 107.5), (0.8, 103.8), (0.4, 87.6), (0.8, 103.0),
]


def fit_rotary_drive_term(points=None, verbose=True):
    """Fits W_drive(f) = k_drive0 + k_drive1*f to Lozano's own WM
    (rotary magnet-assembly + valve drive motor) measurements.

    This is a LINEAR fit in f, not f^2 -- deliberately, because the data
    itself is close to linear (R^2=0.965 for linear vs. a visibly worse
    fit for f^2) and because linear-in-f is what Coulomb/viscous bearing
    and detent-cogging friction (torque roughly independent of speed,
    power = torque*omega ~ f) predicts, whereas eddy-current loss (~f^2 H^2,
    already covered by StateDependentLossModel.k_eddy) is a DIFFERENT
    physical mechanism this term is not meant to re-fit.

    Why this needs its own term rather than folding into CORE's k_eddy/
    k_pump/base_frac: Table 3 shows WM (87-145 W) is comparable to or
    EXCEEDS Qc itself (61-120 W) in every row, and barely scales with
    frequency (3.5x frequency change -> only 1.6x WM change) -- far too
    weak for f^2 eddy scaling and far too large relative to Qc for
    base_frac (CORE's fitted base_frac=0.061 would predict only
    ~4-7W of base overhead at these Qc values, an order of magnitude
    below the measured 87-145W). This is mechanical drivetrain power
    (bearing friction, detent/cogging torque to move the permanent-magnet
    array and rotary valve), not eddy-current or Darcy-flow pumping loss,
    and CORE's 3 calibration devices (Astronautics, DTU, Tusek) do not
    include a rotary-drive-dominated device, so CORE has no information
    about this loss channel at all.
    """
    f = np.array([p[0] for p in points or LOZANO_DRIVE_CALIBRATION_HZ_W])
    WM = np.array([p[1] for p in points or LOZANO_DRIVE_CALIBRATION_HZ_W])
    A = np.vstack([np.ones_like(f), f]).T
    coeffs, _, _, _ = np.linalg.lstsq(A, WM, rcond=None)
    k_drive0, k_drive1 = coeffs
    pred = A @ coeffs
    ss_res = float(np.sum((WM - pred) ** 2))
    ss_tot = float(np.sum((WM - WM.mean()) ** 2))
    r2 = 1 - ss_res / ss_tot
    if verbose:
        print(f"Rotary-drive term fit to Lozano et al. (2016) Table 3 WM data "
              f"({len(f)} points):")
        print(f"  W_drive(f) = {k_drive0:.2f} + {k_drive1:.2f} * f    "
              f"[W, f in Hz]   R^2={r2:.3f}")
        print(f"  (i.e. a large near-constant ~{k_drive0:.0f}W drivetrain "
              f"overhead plus a weak, roughly-linear-in-f term -- NOT an "
              f"eddy-current f^2 scaling)")
    return {"k_drive0": float(k_drive0), "k_drive1": float(k_drive1), "r2": r2}


class StateDependentLossModel:
    def __init__(self, k_eddy=None, k_pump=None, base_frac=None):
        if k_eddy is None or k_pump is None or base_frac is None:
            cal = calibrate_loss_coefficients(verbose=False)
            k_eddy = k_eddy if k_eddy is not None else cal["k_eddy"]
            k_pump = k_pump if k_pump is not None else cal["k_pump"]
            base_frac = base_frac if base_frac is not None else cal["base_frac"]
        self.k_eddy = k_eddy
        self.k_pump = k_pump
        self.base_frac = base_frac

    def parasitic_power(self, frequency, mu0H, mdot, Qc, pumping_power_override=None,
                         intragranular_eddy_power_W=0.0):
        """W_eddy + W_pump + W_base, as documented in the module docstring.

        Phase 15 addition: `pumping_power_override`, default None, changes
        NOTHING about existing behaviour when omitted (every existing
        caller -- amr_cycle.py without a particle_diameter, sensitivity.py,
        rsm.py, the pre-Phase-15 optimize.py -- keeps getting the
        CORE-calibrated generic W_pump = k_pump*mdot**2 term exactly as
        before). When provided, it REPLACES that generic Darcy-flow term
        with a caller-supplied geometry-explicit hydraulic pumping-power
        value (W) -- used by amr_cycle.AMRSystem when a particle_diameter
        is given, so the geometry-explicit pumping power from
        thermal.pumping_power_packed_bed() (Tusek et al. 2013 friction-
        factor correlation) is not DOUBLE-COUNTED against this module's
        own generic, independently-CORE-calibrated k_pump term -- see
        core/amr_cycle.py's AMRSystem docstring and core/optimize.py's
        Phase 15 section for the full reasoning. k_eddy and base_frac are
        unaffected either way, since geometry does not change eddy-current
        or baseline-electronics losses.

        Phase 27 addition: `intragranular_eddy_power_W`, default 0.0,
        changes NOTHING about existing behavior when omitted (every
        pre-Phase-27 caller keeps getting exactly W_eddy + W_pump + W_base
        as before). When provided (by amr_cycle.AMRSystem when a
        particle_diameter is set -- see its own `_geometry_eddy_power_W()`),
        it is ADDED to (not swapped in place of) the CORE-calibrated
        W_eddy term -- deliberately DIFFERENT handling from
        pumping_power_override's REPLACE semantics, because this new term
        represents a physically DIFFERENT loss channel (eddy-current
        self-heating within the MCM particles/plates themselves) than
        k_eddy's own CORE-calibrated support-structure eddy losses -- see
        core.thermal.intragranular_eddy_power()'s own docstring for the
        full distinction. The two mechanisms occur simultaneously in a
        real device, so they add rather than override."""
        W_eddy = self.k_eddy * frequency ** 2 * mu0H ** 2 + intragranular_eddy_power_W
        W_pump = self.k_pump * mdot ** 2 if pumping_power_override is None else pumping_power_override
        W_base = self.base_frac * Qc
        return W_eddy + W_pump + W_base



class RotaryDriveLossModel(StateDependentLossModel):
    """StateDependentLossModel (eddy/pump/base_frac, CORE-calibrated by
    default) PLUS an additional rotary magnet-assembly/valve drivetrain
    term W_drive = k_drive0 + k_drive1*frequency, calibrated to Lozano et
    al. (2016)'s own directly-measured WM data (see fit_rotary_drive_term()
    docstring for why this needs to be a separate additive term rather than
    a re-fit of k_eddy/k_pump/base_frac).

    Intended for devices whose cooling is dominated by a rotary permanent-
    magnet-assembly + rotary valve drivetrain of a similar class to the
    Lozano/POLO-UFSC prototype -- NOT a general "rotary AMR" flag (both
    Astronautics and DTU are also rotary and are already well-represented
    by CORE; what distinguishes Lozano is the drivetrain's overhead
    relative to its own small Qc, not rotary operation per se). Use
    per-device, not as a blanket replacement for StateDependentLossModel.
    """

    def __init__(self, k_eddy=None, k_pump=None, base_frac=None,
                 k_drive0=None, k_drive1=None):
        super().__init__(k_eddy=k_eddy, k_pump=k_pump, base_frac=base_frac)
        if k_drive0 is None or k_drive1 is None:
            fit = fit_rotary_drive_term(verbose=False)
            k_drive0 = k_drive0 if k_drive0 is not None else fit["k_drive0"]
            k_drive1 = k_drive1 if k_drive1 is not None else fit["k_drive1"]
        self.k_drive0 = k_drive0
        self.k_drive1 = k_drive1

    def parasitic_power(self, frequency, mu0H, mdot, Qc, pumping_power_override=None,
                         intragranular_eddy_power_W=0.0):
        base = super().parasitic_power(frequency, mu0H, mdot, Qc,
                                        pumping_power_override=pumping_power_override,
                                        intragranular_eddy_power_W=intragranular_eddy_power_W)
        W_drive = self.k_drive0 + self.k_drive1 * frequency
        return base + W_drive

if __name__ == "__main__":
    calibrate_loss_coefficients()
    print()
    run_extended_diagnostic()
    print()
    run_further_extended_diagnostic()