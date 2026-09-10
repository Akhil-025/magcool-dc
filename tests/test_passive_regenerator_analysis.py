"""
Unit tests for core/passive_regenerator_analysis.py .
"""
import os
import tempfile

import pytest

from core.passive_regenerator_analysis import (
    compare_candidate_materials,
    span_sweep,
    run_passive_regenerator_analysis,
    CANDIDATE_MATERIALS,
    T_COLD_K,
    SPAN_K,
    regenerator_specific_area_per_m,
    particle_diameter_for_specific_area,
    check_against_barclay_1991_reference_point,
    regenerator_specific_area_design_point,
    TYPICAL_PACKED_BED_POROSITY,
)
from core.mce_material import GADOLINIUM


def test_compare_candidate_materials_returns_all_candidates():
    base, results = compare_candidate_materials(verbose=False)
    assert len(results) == len(CANDIDATE_MATERIALS)
    assert base.COP > 0


def test_results_sorted_descending_by_augmented_cop():
    _, results = compare_candidate_materials(verbose=False)
    cops = [r.augmented_COP for r in results]
    assert cops == sorted(cops, reverse=True)


def test_gd_ranks_first_at_representative_point():
    """Gd's Curie temperature (294K) sits inside the representative ASHRAE
    window (291.15-301.15K); this is the concrete, checked instance of the
    'alignment' claim this module's own docstring makes."""
    _, results = compare_candidate_materials(T_cold=T_COLD_K, T_hot=T_COLD_K + SPAN_K,
                                               verbose=False)
    assert results[0].material_name == GADOLINIUM.name


def test_span_sweep_returns_one_row_per_span():
    spans = (5.0, 10.0, 20.0)
    rows = span_sweep(spans_K=spans, verbose=False)
    assert len(rows) == len(spans)
    assert [r["span_K"] for r in rows] == list(spans)


def test_run_passive_regenerator_analysis_writes_file_and_returns_dict():
    with tempfile.TemporaryDirectory() as d:
        out_path = os.path.join(d, "sub", "passive_regenerator_analysis.txt")
        result = run_passive_regenerator_analysis(out_path=out_path)
        assert os.path.exists(out_path)
        with open(out_path) as f:
            content = f.read()
        assert "PHASE 21" in content
        assert "candidate_results" in result
        assert "span_sweep" in result


def test_run_passive_regenerator_analysis_no_file_write():
    result = run_passive_regenerator_analysis(out_path=None)
    assert result["candidate_results"]


# --- regenerator specific-area design guidance (Barclay & Sarangi 1984,
#     via Tishin & Spichkin (2003) p.373, OCR'd from the book's own scan) ---

def test_specific_area_matches_barclay_1991_reference_point():
    """Tishin & Spichkin (2003) p.373 quotes Barclay (1991): N_tu=500 at
    0.5 Hz requires 25000 m^-1 specific area, achieved by 240 um
    particles. This is the concrete number the OCR pass actually
    recovered from the book -- checked here, not assumed."""
    matches, computed = check_against_barclay_1991_reference_point()
    assert matches
    assert computed == pytest.approx(25000.0, rel=1e-6)


def test_specific_area_sphere_relation():
    # a = 6/d for a sphere: halving diameter doubles specific area.
    a1 = regenerator_specific_area_per_m(200e-6)
    a2 = regenerator_specific_area_per_m(100e-6)
    assert a2 == pytest.approx(2 * a1)


def test_specific_area_diameter_are_inverses():
    d = 150e-6
    a = regenerator_specific_area_per_m(d)
    assert particle_diameter_for_specific_area(a) == pytest.approx(d)


def test_specific_area_rejects_nonpositive_diameter():
    with pytest.raises(ValueError):
        regenerator_specific_area_per_m(0)
    with pytest.raises(ValueError):
        regenerator_specific_area_per_m(-1e-6)


def test_design_point_reports_default_porosity_and_note():
    dp = regenerator_specific_area_design_point(240e-6)
    assert dp.porosity == TYPICAL_PACKED_BED_POROSITY
    assert dp.specific_area_per_m == pytest.approx(25000.0, rel=1e-6)
    assert "pressure drop" in dp.note


def test_design_point_smaller_particles_increase_specific_area():
    coarse = regenerator_specific_area_design_point(300e-6)
    fine = regenerator_specific_area_design_point(150e-6)
    assert fine.specific_area_per_m > coarse.specific_area_per_m