import pytest

from core.validation_system import (
    load_benchmarks, run_system_validation, run_curve_validation,
)


def test_load_benchmarks_has_device_group_column():
    rows = load_benchmarks()
    assert len(rows) > 0
    for r in rows:
        assert "device_group" in r
        assert r["device_group"]  # non-empty


def test_device_groups_with_multiple_rows_share_static_params():
    """Companion rows in the same device_group must describe the same
    physical device (only span_K/Qc_W/COP may legitimately differ)."""
    rows = load_benchmarks()
    groups = {}
    for r in rows:
        groups.setdefault(r["device_group"], []).append(r)
    device_identity_fields = ["material", "mu0H_T", "mass_MCM_kg"]
    # groups where every row independently varies its own frequency (a real
    # multi-point sweep across operating conditions, not a fixed-condition
    # span sweep) are exempt from the frequency_Hz match -- see the guard in
    # validation_system.run_curve_validation()
    multi_point_groups = {"Lozano_POLO_UFSC_2016"}
    for group_name, group_rows in groups.items():
        if len(group_rows) < 2:
            continue
        first = group_rows[0]
        fields = device_identity_fields if group_name in multi_point_groups \
            else device_identity_fields + ["frequency_Hz"]
        for r in group_rows[1:]:
            for field in fields:
                assert r[field] == first[field]


def test_run_system_validation_still_returns_four_point_results():
    """The original point-wise validation (Phase 2/6) must be unaffected by
    adding the device_group column: of the original 5 devices, DTU, Tusek
    and Okamura calibrate successfully; Risoe reports a "no calibration
    found" status dict rather than a numeric result, as before.

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
    found" -- so 7 total with COP results, 13 total (5 original + 8 Lozano)."""
    results = run_system_validation()
    with_cop = [r for r in results if "COP_error_pct" in r]
    assert len(with_cop) == 7
    assert len(results) == 13


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