import pytest

from core.validation_system import (
    load_benchmarks, run_system_validation, run_curve_validation,
    run_field_sensitivity_check, run_capacity_only_calibration_check,
    run_tusek_multipoint_curve_validation, _load_tusek_curve,
    calibrate_and_check, infer_cycle_type_for_device, run_cycle_type_validation,
    diagnose_qc_feasibility_reopening, run_calibrated_gd_system_level_comparison,
)
from core.mce_material import GADOLINIUM


def test_load_benchmarks_has_device_group_column():
    rows = load_benchmarks()
    assert len(rows) > 0
    for r in rows:
        assert "device_group" in r
        assert r["device_group"]  # non-empty


def test_device_groups_with_multiple_rows_share_static_params():
    """Companion rows in the same device_group must describe the same
    physical device (only span_K/Qc_W/COP may legitimately differ) --
    UNLESS the group is a documented multi-point sweep across a specific
    operating condition (frequency for Lozano, magnetic field for the new
    Chubu Electric/Toshiba pair -- Paper-Mining Pass Part 2, §1), in which
    case that ONE condition is exempted while everything else must still
    match."""
    rows = load_benchmarks()
    groups = {}
    for r in rows:
        groups.setdefault(r["device_group"], []).append(r)
    device_identity_fields = ["material", "mu0H_T", "mass_MCM_kg", "frequency_Hz"]
    # groups where every row independently varies ONE specific operating
    # condition (a real multi-point sweep, not a fixed-condition span
    # sweep) are exempt from matching on that one field -- see the guards
    # in validation_system.run_curve_validation() (frequency) and the
    # Chubu Electric/Toshiba rows' own source notes (field).
    field_varying_groups = {
        "Lozano_POLO_UFSC_2016": "frequency_Hz",
        "ChubuToshiba_Gd_2016": "mu0H_T",
    }
    for group_name, group_rows in groups.items():
        if len(group_rows) < 2:
            continue
        first = group_rows[0]
        exempt_field = field_varying_groups.get(group_name)
        fields = [f for f in device_identity_fields if f != exempt_field]
        for r in group_rows[1:]:
            for field in fields:
                assert r[field] == first[field]


def test_run_system_validation_still_returns_four_point_results():
    """The original point-wise validation must be unaffected by
    adding the device_group column: of the original 5 devices, DTU and
    Okamura calibrate successfully; Risoe reports a "no calibration found"
    status dict rather than a numeric result, as before.

    UPDATED (Paper-Mining Pass Part 4): Tusek_singlebed_Gd_2010 now ALSO
    reports "no calibration found", not a regression -- its field was
    corrected from an unverified 1.69T to the paper's own stated 1.15T
    (Tusek et al. 2013, Abstract: "the magnetic flux density is 1.15 T"),
    and at the weaker, correct field the old (span=15K, Qc=6.5W) point is
    no longer reachable at all. This is a genuine open item (a real
    (span, Qc) pair for this device at 1.15T still needs proper
    digitization from Figs. 10-11 -- see
    results/tusek_ate2013_figs_notes.md), not a bug in this validation
    code. with_cop count therefore drops from 7 to 6.

    Astronautics_rotary_2014 is now a "no calibration found" row
    too, NOT a regression: it used to "calibrate" only because it was run
    against GADOLINIUM as an explicitly-flagged stand-in. Now that it uses
    the real LAFESIH_FIRST_ORDER material (core/first_order_mce.py,
    Tc=287K, calibrated to the commonly-cited Fujieda et al. 2002 La(Fe,Si)13Hy
    composition) and a device-appropriate T_cold=305K instead of the Gd
    default, the model correctly shows this material's first-order
    transition is far too narrow and centered too low (peak dTad near
    ~298K even accounting for the field shift) to produce useful cooling at
    this device's actual ~305-316K operating window -- exactly what you'd
    expect from representing a SIX-layer, Curie-graded bed (real layer Tc
    304-316K, see first_order_mce.py's LAFESIH_FIRST_ORDER comment) with a
    single Tc=287K material. This is a genuine finding about the
    single-layer approximation's limits, not something to paper over by
    retuning Tc to force a fit -- see ROADMAP.md .
     added 8 real (span>0, Qc, COP) Lozano POLO/UFSC (2016) rows
    (r1-r8; the zerospan/maxspan endpoint rows are filtered exactly like
    every other device's endpoint rows), of which 4 (r4, r6, r7, r8)
    successfully calibrate and 4 (r1, r2, r3, r5) report "no calibration
    found" -- so 7 total with COP results, 13 total (5 original + 8 Lozano).

    Paper-Mining Pass Part 3, §1 added DTU_MagQueen_2018 (a DERIVED
    Qc/COP pair, from the source paper's own reported heating power/COP --
    see that row's own CSV source note) with COP=4.0 -- it reports "no
    calibration found" (same as Astronautics_rotary_2014, another
    LAFESIH_FIRST_ORDER-material row), so with_cop stays at 6 but the
    total grows to 14. Cooltech_2013_rotary (added in the same pass) has
    NO reported COP, so it's filtered by calibrate_and_check() before
    reaching either count here -- see
    run_capacity_only_calibration_check() instead, which does cover it.

    ROADMAP.md Group A completion pass: Tusek_singlebed_Gd_2010's span_K/
    Qc_W/COP were replaced with a genuinely pixel-digitized point from
    Figs. 10-11 (AMR (A.), V*=0.95, span=7.26K/Qc=5.27W/COP=5.38 -- see
    data/tusek_ate2013_figs/notes.md and the row's own CSV source note).
    Unlike the placeholder it replaces, this point DOES calibrate (mdot=
    0.0074 kg/s), so with_cop RISES from 6 to 7 here (total stays 15: the
    row was already present as a "no calibration found" entry before).
    A second row, Tusek_singlebed_Gd_2010_spanceiling (the paper's directly
    -stated 19.8K/0W zero-capacity point), was also added in this pass but
    has no COP, so calibrate_and_check() filters it out entirely (same as
    every other zero-span/zero-Qc endpoint row) -- it does not appear in
    run_system_validation()'s output at all, only in run_curve_validation()
    (see test_curve_validation_covers_multi_point_groups below).

    Paper-Mining Pass Part 5 added DTU_Eriksen_rotary_Gd_2015, a genuinely
    new/independent DTU rotary Gd device from the primary paper (Eriksen
    et al., Int. J. Refrigeration 2015) -- NOT the same device as the
    existing DTU_rotary_Gd_2016 row (that citation remains unverified/
    unlocated; the numbers don't match this 2015 paper at all). It
    calibrates successfully (COP error +4.9%), so with_cop rises to 7 and
    the total to 15.

    Paper-Mining Pass Part 6 CORRECTED DTU_rotary_Gd_2016 (818W, 10.1K,
    COP=4.2), the row referenced immediately above as "unverified/
    unlocated": it has now been located (D. Eriksen, K. Engelbrecht,
    C.R.H. Bahl, R. Bjork, "Exploring the efficiency potential for an
    active magnetic regenerator," Sci. Technol. Built Environ. 22(5)
    (2016) 527-533, reproduced as Ch.6 of Eriksen's 2016 DTU PhD thesis)
    and renamed DTU_Eriksen_MAGGIE_2016 with its real reported numbers:
    81.5W at a 15.5K span, COP=3.6. It turns out that this genuine,
    verified operating point does NOT calibrate under the current
    single-Tc-Gd cycle model (span exceeds the model's reachable range at
    this device's field/mass/frequency, same failure mode as
    Risoe_DTU_Gd_2011 -- see loss_model.py's docstring for why: the real
    device's Curie-graded 11-layer bed reaches spans a single-Tc
    approximation structurally cannot). So the OLD fabricated number
    happened to calibrate (which is presumably why it was hardcoded to
    begin with) but the REAL number does not -- with_cop therefore DROPS
    from 7 to 6 (losing this device's contribution), while total stays at
    15 (the row is still present with a "no calibration found" status,
    like Astronautics_rotary_2014, Risoe_DTU_Gd_2011, DTU_MagQueen_2018,
    and 4 of the 8 Lozano rows -- present, not silently dropped).

    Gd physics-fix pass (exact isentropic DeltaT_ad + Sommerfeld
    electronic entropy/heat-capacity term, core/mce_material.py's
    GADOLINIUM.sommerfeld_gamma_J_per_molK2): adding the electronic term
    to total_heat_capacity() changes total_heat_capacity() everywhere
    GADOLINIUM is used, including the reachable-mdot search this
    calibration loop runs for every Gd-material row. That shift pushes
    Lozano_POLO_UFSC_2016_r4 (span=6.1K) just outside its previously-
    reachable mdot range -- a genuine, expected side effect of a real
    physics correction, not a bug in the calibration search itself.
    with_cop therefore drops again, from 6 to 5... but empirically only
    r4 is lost (Tusek_singlebed_Gd_2010, Okamura_Hirano_2013,
    DTU_Eriksen_rotary_Gd_2015, and Lozano r6/r7/r8 all still calibrate),
    so with_cop lands at 6, not 5 -- see this test's own assertion below
    for the current count; total stays at 15 either way (r4 is still
    present with a "no calibration found" status)."""
    results = run_system_validation()
    with_cop = [r for r in results if "COP_error_pct" in r]
    assert len(with_cop) == 6
    assert len(results) == 15


def test_curve_validation_covers_multi_point_groups():
    """Paper-Mining Pass Part 6: DTU_rotary_Gd_2016 (which had a same-
    frequency zero-Qc companion row, DTU_rotary_Gd_2016_maxspan) was
    corrected and renamed DTU_Eriksen_MAGGIE_2016 (see
    data/amr_experimental_benchmarks.csv and loss_model.py's docstring).
    The real primary source (Eriksen 2016 PhD thesis, Ch.6) does not
    report a same-frequency no-load-span companion point for its 0.61Hz/
    15.5K/81.5W operating point -- the thesis's only other no-load-span
    number (29.2K) is from a different chapter's later test campaign at
    1.4Hz, not a valid same-condition companion -- so no maxspan row was
    fabricated to replace the old one, and DTU_Eriksen_MAGGIE_2016 is a
    single-row device_group, correctly absent from curve validation
    (which requires >=2 rows per group).

    ROADMAP.md Group A completion pass: Tusek_singlebed_Gd_2010 gained a
    same-field/mass/frequency zero-Qc companion row (span-ceiling point,
    19.8K/0W, stated directly in the paper's text -- see that row's CSV
    source note), so it now joins the 2-point-companion groups here too."""
    results = run_curve_validation(verbose=False)
    group_names = {r["device_group"] for r in results}
    # the 3 groups with a fixed-condition companion span point, plus
    # Lozano_POLO_UFSC_2016, whose 8 independent (own frequency/flow)
    # rows are reported here with a "multi-point set" status (see the
    # guard in run_curve_validation()) rather than run through 2-point
    # pairing
    assert group_names == {"Astronautics_rotary_2014", "Tusek_singlebed_Gd_2010",
                            "Risoe_DTU_Gd_2011", "Lozano_POLO_UFSC_2016"}


def test_curve_validation_tusek_spanceiling_companion_matches_model():
    """Genuine finding (ROADMAP.md Group A): at the calibrated mdot fitted
    to Tusek AMR(A)'s V*=0.95 anchor point (span=7.26K, Qc=5.27W), the
    model's own predicted no-load span happens to land almost exactly on
    the paper's directly-stated 19.8K zero-capacity point -- both model
    and literature give Qc=0W there. Locked down so a future, unrelated
    change to amr_cycle.py that silently changes this doesn't go
    unnoticed."""
    results = run_curve_validation(verbose=False)
    tusek = next(r for r in results if r["device_group"] == "Tusek_singlebed_Gd_2010")
    assert tusek["companion_span_K"] == pytest.approx(19.8)
    assert tusek["companion_Qc_lit_W"] == pytest.approx(0.0)
    assert tusek["companion_Qc_model_W"] == pytest.approx(0.0, abs=0.5)


def test_curve_validation_companion_not_used_in_calibration():
    """The calibrated mdot must match the point-wise calibration at the
    anchor span -- i.e. curve validation reuses, not re-fits, mdot."""
    results = run_curve_validation(verbose=False)
    point_results = {r["device"]: r for r in run_system_validation()}
    for r in results:
        if "mdot_calibrated_kg_s" not in r or r["mdot_calibrated_kg_s"] is None:
            continue
        anchor_point = point_results.get(r["device_group"])
        if anchor_point is None:
            continue
        assert r["mdot_calibrated_kg_s"] == pytest.approx(
            anchor_point["mdot_calibrated_kg_s"], rel=1e-6)


def test_curve_validation_risoe_reports_no_calibration():
    results = run_curve_validation(verbose=False)
    risoe = next(r for r in results if r["device_group"] == "Risoe_DTU_Gd_2011")
    assert risoe.get("status") == "no calibration found at anchor point"


def test_curve_validation_no_longer_covers_dtu():
    """Paper-Mining Pass Part 6: DTU_rotary_Gd_2016's companion-span test
    (the model predicting ~0W at the device's no-load span) is retired
    along with the fabricated row it depended on -- see
    test_curve_validation_covers_multi_point_groups's docstring and
    data/amr_experimental_benchmarks.csv's DTU_Eriksen_MAGGIE_2016 row.
    This test just confirms no group calling itself DTU_rotary_Gd_2016
    lingers in curve validation output."""
    results = run_curve_validation(verbose=False)
    group_names = {r["device_group"] for r in results}
    assert "DTU_rotary_Gd_2016" not in group_names


def test_field_sensitivity_check_finds_chubu_rows():
    """The new ChubuToshiba_Gd_2016 group (Paper-Mining Pass Part 2, §1)
    should be found and produce a result dict, even though the outcome
    itself is 'no calibration found' (see the next test) -- this pins
    down that the rows/grouping are wired correctly, independent of
    whether the underlying AMR model happens to calibrate them."""
    result = run_field_sensitivity_check(verbose=False)
    assert result is not None
    assert result["device_group"] == "ChubuToshiba_Gd_2016"


def test_field_sensitivity_check_chubu_reports_honest_no_calibration():
    """Actual, documented finding: the Chubu Electric/Toshiba 4T anchor
    point (26K span at 4.856kg Gd, 0.167Hz) exceeds this single-stage 0-D
    model's max achievable span_fraction (dTad_noload at 4T is too small
    relative to the 26K span for ANY mdot to reach the reported 100W) --
    same 'no calibration found' outcome already documented for Risoe_DTU
    and Astronautics_rotary_2014, not a bug introduced by this row. Locked
    down so a future, unrelated change to amr_cycle.py that silently makes
    this calibrate doesn't go unnoticed without someone updating this test
    (and ROADMAP.md) accordingly."""
    result = run_field_sensitivity_check(verbose=False)
    assert result.get("status") == "no calibration found at anchor field"


def test_field_sensitivity_check_missing_group_returns_none():
    result = run_field_sensitivity_check(verbose=False, device_group="NoSuchDevice_2099")
    assert result is None

def test_capacity_only_check_includes_cooltech_stress_test():
    results = run_capacity_only_calibration_check(verbose=False)
    devices = {r["device"] for r in results}
    assert "Cooltech_2013_rotary" in devices


def test_capacity_only_check_excludes_rows_with_reported_cop():
    """Rows with a reported COP (e.g. DTU_MagQueen_2018) are already
    covered by run_system_validation() -- this function should not
    double-report them."""
    results = run_capacity_only_calibration_check(verbose=False)
    devices = {r["device"] for r in results}
    assert "DTU_MagQueen_2018" not in devices
    assert "Astronautics_rotary_2014" not in devices


def test_capacity_only_check_cooltech_stress_test_does_not_calibrate():
    """Actual, documented finding (Paper-Mining Pass Part 3, §1): the
    42K-span Cooltech_2013_rotary row -- the largest span in this
    benchmark set -- does NOT calibrate at any fluid flow rate in
    [1e-6, 5] kg/s, at the mass=1.0kg fallback this row uses (no mass was
    reported in the source). Locked down so a future, unrelated change to
    amr_cycle.py that silently makes this calibrate doesn't go unnoticed."""
    results = run_capacity_only_calibration_check(verbose=False)
    cooltech = next(r for r in results if r["device"] == "Cooltech_2013_rotary")
    assert cooltech["status"] == "no calibration found"


def test_dtu_magqueen_row_uses_lafesih_material_and_derived_qc():
    """DTU_MagQueen_2018's material field contains 'La', so it must be
    routed to LAFESIH_FIRST_ORDER (Astronautics_rotary_2014's material),
    per the CSV row's own source note. Qc_W=1200 is a DERIVED value
    (Qc = Qh*(1-1/COP_h) = 1500*(1-1/5)), not directly reported -- pin the
    exact derivation so a future edit can't silently substitute a
    different (undocumented) number."""
    rows = load_benchmarks()
    row = next(r for r in rows if r["device"] == "DTU_MagQueen_2018")
    assert "La" in row["material"]
    assert float(row["Qc_W"]) == pytest.approx(1500.0 * (1 - 1 / 5.0))
    assert float(row["COP"]) == pytest.approx(5.0 - 1.0)


def test_cooltech_row_has_no_reported_cop():
    """Cooltech_2013_rotary is a capacity/span-only row (no COP in the
    source) -- confirms it's correctly excluded from COP-based validation
    (run_system_validation()) and only reachable via
    run_capacity_only_calibration_check()."""
    rows = load_benchmarks()
    row = next(r for r in rows if r["device"] == "Cooltech_2013_rotary")
    assert row["COP"] == ""

def test_load_tusek_curve_returns_sorted_points():
    """AMR (A.) at V*=0.95 (fig10_data.csv) has 3 digitized points; must
    come back sorted by ascending span."""
    pts = _load_tusek_curve(amr="A", v_star=0.95)
    assert len(pts) == 3
    spans = [p[0] for p in pts]
    assert spans == sorted(spans)


def test_load_tusek_curve_missing_combo_returns_empty():
    pts = _load_tusek_curve(amr="A", v_star=99.0)
    assert pts == []


def test_tusek_multipoint_curve_validation_calibrates_at_anchor():
    """AMR (A.) V*=0.95's first (lowest-span) point is the same anchor
    point used for the Tusek_singlebed_Gd_2010 CSV row, so the two must
    calibrate to the identical mdot."""
    out = run_tusek_multipoint_curve_validation(verbose=False, amr="A", v_star=0.95)
    assert out["status"] == "ok"
    assert out["anchor_span_K"] == pytest.approx(7.26)
    point_results = {r["device"]: r for r in run_system_validation()}
    tusek_point = point_results["Tusek_singlebed_Gd_2010"]
    assert out["mdot_calibrated_kg_s"] == pytest.approx(
        tusek_point["mdot_calibrated_kg_s"], rel=1e-6)


def test_tusek_multipoint_curve_validation_predicts_both_other_points():
    out = run_tusek_multipoint_curve_validation(verbose=False, amr="A", v_star=0.95)
    spans_predicted = {p["span_K"] for p in out["predictions"]}
    assert spans_predicted == {12.23, 14.75}


def test_tusek_multipoint_curve_validation_span_reopening_now_clamped():
    """Historical note (superseded finding, kept for context): before the
    span-reopening clamp was wired into this function
    (`cooling_capacity_span_sweep()`, core/amr_cycle.py), this repo's
    single-Tc 0-D Qc(span) model was non-monotonic at this calibrated
    mdot -- it dropped to ~0W around span=8K, then rose again to a large,
    spurious local maximum (~20.5W, +910% vs. the literature 2.03W)
    before finally falling back to 0 past span~14K. A 2-point
    companion check alone would not have caught this; only a genuine
    3+-point curve-shape check exposed it.

    Current, correct behavior: the SAME underlying near-Tc discontinuity
    is still there (qc_feasibility_reopening still reports `reopens`),
    but the reported Qc at this span is now the CLAMPED value (0W, not
    the spurious 20.5W reopening peak) -- physically consistent with the
    model's own behavior at every other span in this reopened window.
    The remaining -100% error at this point is a genuine, honestly-
    reported model limitation (the model just can't produce ANY nonzero
    Qc here, vs. the literature's small positive 2.03W), not the
    reopening artifact itself. Locked down so a future, unrelated change
    to amr_cycle.py that silently un-clamps this doesn't go unnoticed."""
    out = run_tusek_multipoint_curve_validation(verbose=False, amr="A", v_star=0.95)
    mid_span_pred = next(p for p in out["predictions"] if p["span_K"] == 12.23)
    assert mid_span_pred["Qc_lit_W"] == pytest.approx(2.03)
    # Clamped: no longer a spurious overshoot far above the anchor's 5.27W.
    assert mid_span_pred["Qc_model_W"] == pytest.approx(0.0, abs=1e-6)
    assert mid_span_pred["Qc_error_pct"] == pytest.approx(-100.0, abs=0.1)
    # The clamp DID fire here (raw model would have been the spurious
    # reopening value) -- confirms this test is actually exercising the
    # clamp, not a case where clamping was a no-op.
    assert mid_span_pred["span_reopening_clamped"] is True
    assert mid_span_pred["Qc_model_raw_W"] > 10.0

    endpoint_pred = next(p for p in out["predictions"] if p["span_K"] == 14.75)
    assert endpoint_pred["Qc_lit_W"] == pytest.approx(0.0)
    assert endpoint_pred["Qc_model_W"] == pytest.approx(0.0, abs=0.5)

    # Root-cause attribution (added after tracing the above symptom to its
    # actual origin): the 12.23K miss is diagnosed as the known near-Tc
    # feasibility-margin reopening, not left as an unexplained number, and
    # is flagged as such on the prediction row itself.
    assert out["qc_feasibility_reopening"]["reopens"] is True
    assert mid_span_pred.get("flagged_near_Tc_nonmonotonicity") is True


def test_diagnose_qc_feasibility_reopening_detects_the_near_Tc_reopening():
    """Gd's own feasibility margin (2*dTad_noload(T_mid) - span), at the
    field/T_cold this Tusek AMR(A) curve uses, is known (see
    core/mce_material.py's magnetic_heat_capacity() docstring) to go
    negative, swing back positive just above Tc=294K, then go negative
    again -- confirm the diagnostic actually catches this rather than
    reporting false confidence."""
    result = diagnose_qc_feasibility_reopening(
        GADOLINIUM, T_cold=289.0, mu0H_max=1.15, span_lo=6.0, span_hi=15.0)
    assert result["reopens"] is True
    assert result["n_feasible_to_infeasible_crossings"] == 2
    assert result["n_infeasible_to_feasible_crossings"] == 1
    # the reopening should sit near where T_mid crosses Tc=294K, i.e.
    # span ~= 2*(294 - 289) = 10K
    assert 9.0 < result["first_reopen_span_K"] < 11.0
    assert result["final_reclose_span_K"] is not None


def test_diagnose_qc_feasibility_reopening_false_when_range_avoids_Tc():
    """A span range whose T_mid never approaches Tc should show a single,
    ordinary feasibility crossing (or none) -- the diagnostic shouldn't
    over-flag unrelated ranges."""
    result = diagnose_qc_feasibility_reopening(
        GADOLINIUM, T_cold=260.0, mu0H_max=1.15, span_lo=0.5, span_hi=4.0)
    assert result["reopens"] is False
    assert result["n_infeasible_to_feasible_crossings"] == 0


def test_tusek_multipoint_curve_validation_missing_curve_reports_status():
    out = run_tusek_multipoint_curve_validation(verbose=False, amr="A", v_star=99.0)
    assert out["status"] == "fewer than 2 points"


def test_tusek_fig10_and_fig11_csvs_have_matching_series():
    """fig10_data.csv (Qc) and fig11_data.csv (COP) digitize the same 9
    series (3 AMR geometries x 3 V* ratios) from the same source figures;
    every (amr, V_star) combination in one must also appear in the other."""
    import csv as _csv
    with open("data/tusek_ate2013_figs/fig10_data.csv") as f:
        combos_10 = {(r["amr"], r["V_star"]) for r in _csv.DictReader(f)}
    with open("data/tusek_ate2013_figs/fig11_data.csv") as f:
        combos_11 = {(r["amr"], r["V_star"]) for r in _csv.DictReader(f)}
    assert combos_10 == combos_11
    assert len(combos_10) == 9

# --- cycle-type validation sensitivity ---

def test_infer_cycle_type_for_device_rotary_vs_other():
    rotary_row = {"device": "Astronautics_rotary_2014",
                  "device_group": "Astronautics_rotary_2014"}
    other_row = {"device": "Tusek_singlebed_Gd_2010",
                 "device_group": "Tusek_singlebed_Gd_2010"}
    assert infer_cycle_type_for_device(rotary_row) == "ericsson"
    assert infer_cycle_type_for_device(other_row) == "brayton"


def test_infer_cycle_type_is_case_insensitive():
    row = {"device": "SOME_ROTARY_DEVICE", "device_group": "SOME_ROTARY_DEVICE"}
    assert infer_cycle_type_for_device(row) == "ericsson"


def test_calibrate_and_check_accepts_cycle_type_kwarg():
    """Regression guard for the threading of cycle_type through
    calibrate_and_check() -- a row that calibrates under brayton should
    also (independently) attempt calibration under ericsson without
    raising, and the two should not be forced to return identical numbers
    (ericsson's qc_multiplier != 1.0 changes the calibrated mdot)."""
    rows = load_benchmarks()
    row = next(r for r in rows if r["device"] == "DTU_Eriksen_rotary_Gd_2015")
    baseline = calibrate_and_check(row, verbose=False, cycle_type="brayton")
    ericsson = calibrate_and_check(row, verbose=False, cycle_type="ericsson")
    assert baseline is not None and ericsson is not None
    assert baseline["mdot_calibrated_kg_s"] != ericsson["mdot_calibrated_kg_s"]


def test_run_cycle_type_validation_returns_one_row_per_cop_target(tmp_path):
    """Every row calibrate_and_check() would treat as a COP validation
    target (span>0, reported Qc and COP) must appear exactly once in
    run_cycle_type_validation()'s output, whether or not it is inferred
    as rotary."""
    out_file = tmp_path / "cycle_type_validation.txt"
    results = run_cycle_type_validation(verbose=False, out_path=str(out_file))
    rows = load_benchmarks()
    expected_devices = {r["device"] for r in rows
                         if calibrate_and_check(r, verbose=False) is not None}
    got_devices = {r["device"] for r in results}
    assert got_devices == expected_devices
    assert out_file.exists()
    assert "PHASE 17" in out_file.read_text()


def test_run_cycle_type_validation_non_rotary_rows_unchanged():
    """A device not inferred as rotary must show cycle_type_inferred ==
    'brayton' and its COP error must exactly match the plain
    run_system_validation() baseline (no re-solve should have happened)."""
    results = run_cycle_type_validation(verbose=False, out_path=None)
    baseline_rows = run_system_validation()
    baseline_by_device = {r["device"]: r for r in baseline_rows if "COP_error_pct" in r}
    for r in results:
        if r.get("cycle_type_inferred") == "brayton" and r["device"] in baseline_by_device:
            assert r["COP_error_pct_cycle_inferred"] == pytest.approx(
                baseline_by_device[r["device"]]["COP_error_pct"])


def test_run_cycle_type_validation_out_path_none_skips_file_write(tmp_path):
    """out_path=None must not attempt any file write (used by the tests
    above and any quick interactive call)."""
    before = set(tmp_path.iterdir())
    results = run_cycle_type_validation(verbose=False, out_path=None)
    assert len(results) > 0
    assert set(tmp_path.iterdir()) == before

def test_calibrated_gd_system_level_comparison_returns_rows(tmp_path):
    """Basic structural check: runs against real, already-calibrating Gd
    benchmark rows, returns a comparison for each, and writes a report."""
    out_path = tmp_path / "calibrated_gd_system_level_comparison.txt"
    result = run_calibrated_gd_system_level_comparison(verbose=False, out_path=str(out_path))
    assert out_path.exists()
    assert len(result["rows"]) > 0
    for r in result["rows"]:
        assert r["COP_plain"] >= 0.0
        assert r["COP_calibrated"] >= 0.0
        assert r["Qc_plain_W"] >= 0.0
        assert r["Qc_calibrated_W"] >= 0.0


def test_calibrated_gd_system_level_comparison_same_mdot_both_materials():
    """The whole comparison is only meaningful if both AMRSystem instances
    for a given row share the SAME calibrated mdot -- otherwise a COP/Qc
    difference could just be a flow-rate artifact, not the material
    physics difference this function is meant to isolate. Checked
    indirectly: re-running calibrate_and_check() for the same device
    must reproduce the exact mdot this function used internally."""
    from core.validation_system import load_benchmarks, calibrate_and_check
    rows = load_benchmarks()
    checked_devices = set()
    for row in rows:
        result = calibrate_and_check(row, verbose=False)
        if result is not None and "status" not in result:
            checked_devices.add(row["device"])
    result = run_calibrated_gd_system_level_comparison(verbose=False)
    compared_devices = {r["device"] for r in result["rows"]}
    # Every device this function actually compared must be one that
    # independently calibrates via the same mechanism -- confirms this
    # function isn't inventing its own separate (and possibly
    # inconsistent) calibration path.
    assert compared_devices.issubset(checked_devices)


def test_calibrated_gd_system_level_comparison_documents_feasibility_zeroing():
    """Real, load-bearing finding from this comparison (see the
    function's own docstring): at least one benchmark device goes from a
    real, positive Qc under plain GADOLINIUM to a hard zero under
    GADOLINIUM_CALIBRATED at the identical calibrated mdot/span/field --
    this is not just a performance argument against wiring
    GADOLINIUM_CALIBRATED in system-wide, it's a correctness one. Locked
    down so a future, unrelated change to either material doesn't
    silently make this finding disappear (or appear) unnoticed."""
    result = run_calibrated_gd_system_level_comparison(verbose=False)
    zeroed = [r for r in result["rows"]
              if r["Qc_calibrated_W"] == 0.0 and r["Qc_plain_W"] > 0.0]
    assert len(zeroed) >= 1
