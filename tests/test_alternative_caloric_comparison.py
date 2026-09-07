import math

from core.alternative_caloric_comparison import (
    CALORIC_CLAIMS,
    compare_all_at_matched_spans,
    compare_physics_models_across_spans,
)


def _fake_vcc_table():
    """Mirrors this repo's own comparison_table.csv span/VCC_COP pairs
    closely enough to test matching logic without depending on the file."""
    return {5.0: 24.46, 10.0: 12.23, 15.0: 8.15, 20.0: 6.11}


def test_every_claim_has_a_source_and_caveat():
    for claim in CALORIC_CLAIMS:
        assert claim.source, f"{claim.system_or_material} missing a source"
        assert claim.caveat, f"{claim.system_or_material} missing a caveat"


def test_nan_claims_excluded_from_comparison():
    """The Barocal vendor claim has no numeric COP and must not silently
    become a 0.0 or otherwise-fabricated comparison point."""
    vcc = _fake_vcc_table()
    results = compare_all_at_matched_spans(vcc, verbose=False)
    systems = [r["system"] for r in results]
    assert not any("Barocal Ltd" in s for s in systems)


def test_exactly_one_measured_device_cop_result():
    """As of the sources checked here, every alternative-caloric COP figure
    is a simulation/projection EXCEPT the electrocaloric Li et al. 2023
    device -- this test locks in that specific, narrow exception so a
    future edit can't silently add (or lose) a "measured win" claim
    without the test noticing."""
    vcc = _fake_vcc_table()
    results = compare_all_at_matched_spans(vcc, verbose=False)
    assert len(results) >= 1
    measured = [r for r in results if r["is_measured_device_cop"] is True]
    assert len(measured) == 1
    assert "Li et al. 2023" in measured[0]["system"]


def test_electrocaloric_li_2023_beats_vcc_at_measured_span():
    """Locks in this module's strongest honest finding: a genuinely
    MEASURED device COP (not simulated/projected) that beats VCC."""
    vcc = _fake_vcc_table()
    results = compare_all_at_matched_spans(vcc, verbose=False)
    li = [r for r in results if "Li et al. 2023" in r["system"]][0]
    assert li["is_measured_device_cop"] is True
    assert li["beats_vcc"] is True


def test_widest_span_used_when_source_gives_no_span():
    vcc = _fake_vcc_table()
    results = compare_all_at_matched_spans(vcc, verbose=False)
    unpinned = [r for r in results if not r["span_was_pinned_by_source"]]
    assert unpinned, "expected at least one claim with no source-pinned span"
    for r in unpinned:
        assert r["matched_span_K"] == max(vcc.keys())


def test_elastocaloric_qian_2023_beats_vcc_at_matched_span():
    """Locks in this module's one 'beats VCC' finding so it's visible if a
    future data update changes it either way."""
    vcc = _fake_vcc_table()
    results = compare_all_at_matched_spans(vcc, verbose=False)
    qian = [r for r in results if "Qian et al. 2023" in r["system"]][0]
    assert qian["beats_vcc"] is True


def test_barocaloric_npg_does_not_beat_vcc_at_5K():
    vcc = _fake_vcc_table()
    results = compare_all_at_matched_spans(vcc, verbose=False)
    npg = [r for r in results if "Neopentylglycol" in r["system"]][0]
    assert math.isclose(npg["matched_span_K"], 5.0)
    assert npg["beats_vcc"] is False


def test_physics_model_sweep_covers_every_span():
    vcc = _fake_vcc_table()
    results, eq_frac, bc_frac = compare_physics_models_across_spans(vcc, verbose=False)
    assert {r["span_K"] for r in results} == set(vcc.keys())
    assert eq_frac >= 0.0
    assert bc_frac >= 0.0


def test_physics_model_neither_elasto_nor_baro_beats_vcc_on_this_grid():
    """Locks in the SAME headline honest finding this module previously had
    for elastocaloric/barocaloric specifically: neither beats VCC anywhere
    on this repo's own 5-20K span grid. Electrocaloric is the one
    exception (see next test) -- but only because its single calibration
    point genuinely is a measured device win, not because this repo's
    models became more favorable to caloric cooling in general."""
    vcc = _fake_vcc_table()
    results, _eq_frac, _bc_frac = compare_physics_models_across_spans(vcc, verbose=False)
    assert not any(r["elastocaloric_beats_vcc"] for r in results)
    assert not any(r["barocaloric_beats_vcc"] for r in results)


def test_physics_model_electrocaloric_beats_vcc_at_measured_span_only():
    """Electrocaloric's calibration model beats VCC at (and near) its one
    genuinely measured span (20.9K, nearest grid point 20K) -- but this is
    an extrapolated scaling law everywhere else in the grid, not a second
    independent measurement, per electrocaloric_cycle.py's own docstring."""
    vcc = _fake_vcc_table()
    results, _eq_frac, _bc_frac = compare_physics_models_across_spans(vcc, verbose=False)
    at_20 = [r for r in results if r["span_K"] == 20.0][0]
    assert at_20["electrocaloric_beats_vcc"] is True
    assert at_20["electrocaloric_is_extrapolated"] is True  # 20K != measured 20.9K
    assert at_20["electrocaloric_is_far_extrapolation"] is False  # within [13.0, 20.9]
