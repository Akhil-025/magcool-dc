import pytest

from core.validation_system import (
    load_benchmarks, run_system_validation, run_curve_validation,
    run_field_sensitivity_check, run_capacity_only_calibration_check,
)


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
    """The original point-wise validation (Phase 2/6) must be unaffected by
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

    Astronautics_rotary_2014 (Phase 9) is now a "no calibration found" row
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
    retuning Tc to force a fit -- see ROADMAP.md Phase 9.
    Phase 7 added 8 real (span>0, Qc, COP) Lozano POLO/UFSC (2016) rows
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

    Paper-Mining Pass Part 5 added DTU_Eriksen_rotary_Gd_2015, a genuinely
    new/independent DTU rotary Gd device from the primary paper (Eriksen
    et al., Int. J. Refrigeration 2015) -- NOT the same device as the
    existing DTU_rotary_Gd_2016 row (that citation remains unverified/
    unlocated; the numbers don't match this 2015 paper at all). It
    calibrates successfully (COP error +4.9%), so with_cop rises to 7 and
    the total to 15."""
    results = run_system_validation()
    with_cop = [r for r in results if "COP_error_pct" in r]
    assert len(with_cop) == 7
    assert len(results) == 15


def test_curve_validation_covers_multi_point_groups():
    results = run_curve_validation(verbose=False)
    group_names = {r["device_group"] for r in results}
    # the 3 groups with a fixed-condition companion span point, plus
    # Lozano_POLO_UFSC_2016, whose 8 independent (own frequency/flow) rows
    # are reported here with a "multi-point set" status (see the guard in
    # run_curve_validation()) rather than run through 2-point pairing
    assert group_names == {"Astronautics_rotary_2014", "DTU_rotary_Gd_2016",
                            "Risoe_DTU_Gd_2011", "Lozano_POLO_UFSC_2016"}


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


def test_curve_validation_dtu_predicts_near_zero_at_noload_span():
    """The DTU device's companion point is its reported no-load (zero
    capacity) span -- the model should predict close to zero there too."""
    results = run_curve_validation(verbose=False)
    dtu = next(r for r in results if r["device_group"] == "DTU_rotary_Gd_2016")
    assert dtu["companion_Qc_model_W"] == pytest.approx(0.0, abs=1.0)


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
    (and ROADMAP.md's Phase 11 note) accordingly."""
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