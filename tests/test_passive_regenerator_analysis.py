"""
Unit tests for core/passive_regenerator_analysis.py (Phase 21).
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