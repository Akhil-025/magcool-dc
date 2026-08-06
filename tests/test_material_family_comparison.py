"""
Unit tests for core/material_family_comparison.py (Track A2 item).
"""
from core.material_family_comparison import (
    build_comparison_table, run_analysis, _tuned_candidate, REPRESENTATIVE_SPAN_K,
)
from core.cascade import GD_FAMILY, LAFESIH_FAMILY, MNFEPSI_FAMILY


def test_build_comparison_table_covers_all_five_candidates_per_span():
    rows = build_comparison_table(spans_K=(10.0,))
    candidates = {r["candidate"] for r in rows}
    assert len(candidates) == 5
    assert len(rows) == 5


def test_fixed_candidates_use_the_same_material_across_spans():
    rows = build_comparison_table(spans_K=(5.0, 20.0))
    gd_rows = [r for r in rows if r["candidate"] == "Gd (fixed)"]
    assert len({r["Tc_used_K"] for r in gd_rows}) == 1


def test_tunable_family_composition_tracks_the_operating_point():
    """A family re-tuned for a hotter span should need a higher Tc than the
    same family re-tuned for a colder span (peak_T(Tc) is monotonic)."""
    material_cold, tc_cold, _ = _tuned_candidate(GD_FAMILY, T_mid_K=293.0)
    material_hot, tc_hot, _ = _tuned_candidate(GD_FAMILY, T_mid_K=305.0)
    assert tc_hot > tc_cold


def test_family_outside_its_tc_window_falls_back_to_gd():
    """MNFEPSI_FAMILY's documented window (295.3-331.2K) sits mostly AT or
    ABOVE the ASHRAE range; targeting a point well below its window should
    report in_range=False and fall back to the family's fallback_material."""
    material, tc, in_range = _tuned_candidate(MNFEPSI_FAMILY, T_mid_K=250.0)
    assert in_range is False
    assert material is MNFEPSI_FAMILY.fallback_material


def test_family_inside_its_tc_window_does_not_fall_back():
    material, tc, in_range = _tuned_candidate(LAFESIH_FAMILY, T_mid_K=296.0)
    assert in_range is True
    assert LAFESIH_FAMILY.tc_min <= tc <= LAFESIH_FAMILY.tc_max
    assert material is not LAFESIH_FAMILY.fallback_material


def test_run_analysis_writes_output_files(tmp_path):
    out_csv = tmp_path / "material_family_comparison.csv"
    out_txt = tmp_path / "material_family_comparison.txt"
    rows = run_analysis(out_csv=str(out_csv), out_txt=str(out_txt))
    assert out_csv.exists()
    assert out_txt.exists()
    assert len(rows) > 0
    text = out_txt.read_text()
    assert "RANKED" in text
    assert "La(Fe,Si)13Hy" in text


def test_fallback_candidates_are_not_given_a_fake_independent_rank(tmp_path):
    """MNFEPSI_FAMILY falls back to plain Gd at the representative span (its
    Tc window doesn't cover that point) -- it must not appear as its own
    numbered rank line, since that would read as a distinct result when it's
    actually just Gd under a different label."""
    out_csv = tmp_path / "material_family_comparison.csv"
    out_txt = tmp_path / "material_family_comparison.txt"
    run_analysis(out_csv=str(out_csv), out_txt=str(out_txt))
    text = out_txt.read_text()
    assert "Not independently ranked" in text
    # no numbered rank line ("  N. (Mn,Fe)2(P,Si)...") for the fallback candidate
    for line in text.splitlines():
        if "(Mn,Fe)2(P,Si)" in line and line.strip()[:1].isdigit():
            assert False, f"fallback candidate given a fake numbered rank: {line!r}"


def test_representative_span_is_one_of_the_swept_spans():
    from core.material_family_comparison import SPANS_K
    assert REPRESENTATIVE_SPAN_K in SPANS_K