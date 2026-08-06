"""
Unit tests for core/design_recommendations.py (Phase 16 item).

Uses small synthetic input dicts/rows rather than real pipeline output
(the same way the rest of core/ is unit-tested) so these tests run in
milliseconds and don't depend on NSGA-III, Sobol sampling, or the
Curie-graded cascade sweep actually executing. See test_main.py's own
docstring for why the full pipeline is never run inside the test suite.
"""
import pytest

from core.design_recommendations import (
    summarize_frequency_lever,
    summarize_material_lever,
    summarize_grading_lever,
    summarize_geometry_lever,
    summarize_field_flow_lever,
    build_report,
)


# ---------------------------------------------------------------------
# Synthetic fixtures, shaped like the real objects each lever consumes.
# ---------------------------------------------------------------------

def _pareto_rows():
    """Shaped like optimize.run_optimization()'s return value."""
    return [
        {"mu0H_max_T": 2.993, "frequency_Hz": 0.322, "fluid_mdot_kgs": 0.4999,
         "mass_regenerator_kg": 11.55, "regen_effectiveness": 0.83, "blow_fraction": 0.414,
         "COP_electrical": 9.494, "Qc_W": 14186.04, "cost_index_USD": 4380.2},
        {"mu0H_max_T": 2.988, "frequency_Hz": 4.231, "fluid_mdot_kgs": 0.499,
         "mass_regenerator_kg": 11.69, "regen_effectiveness": 0.924, "blow_fraction": 0.414,
         "COP_electrical": 2.89, "Qc_W": 19649.48, "cost_index_USD": 4423.3},
        {"mu0H_max_T": 2.994, "frequency_Hz": 0.321, "fluid_mdot_kgs": 0.4984,
         "mass_regenerator_kg": 1.63, "regen_effectiveness": 0.73, "blow_fraction": 0.413,
         "COP_electrical": 9.449, "Qc_W": 13750.19, "cost_index_USD": 619.4},
    ]


def _sobol_Si(frequency_ST=0.8564):
    """Shaped like a SALib Si dict, indexed in sensitivity.PROBLEM["names"]
    order: [mu0H_max_T, frequency_Hz, fluid_mdot_kgs, regen_effectiveness,
    parasitic_fraction]."""
    from core.sensitivity import PROBLEM
    idx = PROBLEM["names"].index("frequency_Hz")
    st = [0.0] * len(PROBLEM["names"])
    st[idx] = frequency_ST
    return {"ST": st, "S1": st}


def _material_rows():
    """Shaped like material_family_comparison.build_comparison_table()'s
    return value, restricted to the representative 10K span."""
    return [
        {"candidate": "Gd (fixed)", "span_K": 10.0, "in_range": True,
         "1stage_COP": 5.09, "1stage_Qc_W": 1443.4},
        {"candidate": "Gd5(SixGe1-x)4(-Ga) (tuned)", "span_K": 10.0, "in_range": True,
         "1stage_COP": 4.95, "1stage_Qc_W": 1352.3},
        {"candidate": "La(Fe,Si)13Hy (tuned)", "span_K": 10.0, "in_range": True,
         "1stage_COP": 7.33, "1stage_Qc_W": 4989.1},
        {"candidate": "(Mn,Fe)2(P,Si) (tuned)", "span_K": 10.0, "in_range": False,
         "1stage_COP": 5.09, "1stage_Qc_W": 1443.4},
    ]


def _graded_row():
    return {"span_K": 10, "Graded_3stage_COP": 2.781, "Graded_3stage_Qc_W": 2388.1,
            "Graded_3stage_n_fallback_to_Gd": 0}


def _gd_cascade_row():
    return {"span_K": 10, "AMR_3stage_COP": 2.415, "AMR_3stage_Qc_W": 1258.0}


def _pb_best_cop_row():
    return (0.5, 821.37, 15.2681)


def _pp_best_cop_row():
    return (0.1, 819.44, 15.2680)


# ---------------------------------------------------------------------
# Lever 1: frequency
# ---------------------------------------------------------------------

def test_frequency_lever_extracts_correct_sobol_ST():
    _, data = summarize_frequency_lever(_sobol_Si(0.8564), _pareto_rows())
    assert data["frequency_ST"] == pytest.approx(0.8564)


def test_frequency_lever_picks_best_cop_and_best_qc_designs_correctly():
    text, data = summarize_frequency_lever(_sobol_Si(), _pareto_rows())
    assert data["best_cop_design"]["COP_electrical"] == 9.494
    assert data["best_qc_design"]["Qc_W"] == 19649.48
    # best-COP design should have a lower frequency than best-Qc design,
    # matching this repo's own documented frequency/loss trade-off
    assert data["best_cop_design"]["frequency_Hz"] < data["best_qc_design"]["frequency_Hz"]
    assert "OPERATING FREQUENCY" in text


def test_frequency_lever_handles_missing_sobol_data_gracefully():
    text, data = summarize_frequency_lever(None, _pareto_rows())
    assert data["frequency_ST"] is None
    assert "OPERATING FREQUENCY" in text  # header still present


def test_frequency_lever_handles_empty_pareto_rows():
    text, data = summarize_frequency_lever(_sobol_Si(), [])
    assert data["best_cop_design"] is None
    assert data["best_qc_design"] is None


# ---------------------------------------------------------------------
# Lever 2: material
# ---------------------------------------------------------------------

def test_material_lever_ranks_lafesih_first():
    text, data = summarize_material_lever(_material_rows(), representative_span_K=10.0)
    assert data["ranked"][0]["candidate"] == "La(Fe,Si)13Hy (tuned)"
    assert "La(Fe,Si)13Hy" in text


def test_material_lever_excludes_out_of_range_fallback_candidates():
    text, data = summarize_material_lever(_material_rows(), representative_span_K=10.0)
    ranked_names = [r["candidate"] for r in data["ranked"]]
    assert "(Mn,Fe)2(P,Si) (tuned)" not in ranked_names


def test_material_lever_reports_gain_relative_to_gd():
    text, _ = summarize_material_lever(_material_rows(), representative_span_K=10.0)
    assert "vs. plain Gd" in text
    assert "+44%" in text  # 7.33/5.09 - 1 = +44%


def test_material_lever_handles_no_candidates_at_span():
    text, data = summarize_material_lever(_material_rows(), representative_span_K=999.0)
    assert data["ranked"] == []
    assert "No in-range tunable candidate found" in text


# ---------------------------------------------------------------------
# Lever 3: Curie grading
# ---------------------------------------------------------------------

def test_grading_lever_reports_graded_vs_plain_gd():
    text, data = summarize_grading_lever(_graded_row(), _gd_cascade_row(), n_stages=3)
    assert "2.781" in text
    assert "2.415" in text
    assert data["graded_row"] is not None


def test_grading_lever_handles_missing_data():
    text, data = summarize_grading_lever(None, None)
    assert "not available" in text
    assert data["graded_row"] is None


# ---------------------------------------------------------------------
# Lever 4: geometry
# ---------------------------------------------------------------------

def test_geometry_lever_reports_both_interior_optima():
    text, data = summarize_geometry_lever(_pb_best_cop_row(), _pp_best_cop_row())
    assert "0.5 mm" in text
    assert "0.1 mm" in text
    assert data["pb_best_cop"][0] == 0.5
    assert data["pp_best_cop"][0] == 0.1


def test_geometry_lever_handles_partial_data():
    text, data = summarize_geometry_lever(_pb_best_cop_row(), None)
    assert data["pp_best_cop"] is None
    assert "0.5 mm" in text


# ---------------------------------------------------------------------
# Lever 5: field/flow balance
# ---------------------------------------------------------------------

def test_field_flow_lever_finds_a_knee_point_within_the_front():
    text, data = summarize_field_flow_lever(_pareto_rows())
    assert data["knee_point"] in _pareto_rows()
    assert "Knee-point" in text


def test_field_flow_lever_handles_empty_front():
    text, data = summarize_field_flow_lever([])
    assert data["knee_point"] is None
    assert "not available" in text


# ---------------------------------------------------------------------
# build_report: end-to-end assembly + graceful degradation + file writing
# ---------------------------------------------------------------------

def test_build_report_assembles_all_five_sections(tmp_path):
    out_path = tmp_path / "design_recommendations.txt"
    result = build_report(
        sobol_state_dependent_Si=_sobol_Si(),
        pareto_rows=_pareto_rows(),
        material_rows=_material_rows(),
        graded_row=_graded_row(),
        gd_cascade_row=_gd_cascade_row(),
        pb_best_cop_row=_pb_best_cop_row(),
        pp_best_cop_row=_pp_best_cop_row(),
        out_path=str(out_path),
    )
    for key in ("frequency", "material", "grading", "geometry", "field_flow"):
        assert key in result
    assert "RECOMMENDED STARTING DESIGN POINT" in result["text"]
    assert out_path.exists()
    assert out_path.read_text() == result["text"] + "\n"


def test_build_report_degrades_gracefully_with_no_inputs(tmp_path):
    """Running the synthesis with nothing computed yet (e.g. a partial or
    failed pipeline run) must not raise -- every lever should report
    itself as unavailable instead."""
    out_path = tmp_path / "design_recommendations.txt"
    result = build_report(out_path=str(out_path))
    assert result["frequency"]["frequency_ST"] is None
    assert result["frequency"]["best_cop_design"] is None
    assert result["field_flow"]["knee_point"] is None
    assert result["grading"]["graded_row"] is None
    assert "OPERATING FREQUENCY" in result["text"]
    assert out_path.exists()


def test_build_report_creates_output_directory_if_missing(tmp_path):
    nested = tmp_path / "nested" / "dir" / "design_recommendations.txt"
    build_report(pareto_rows=_pareto_rows(), out_path=str(nested))
    assert nested.exists()


def test_build_report_recommended_design_point_matches_knee_point(tmp_path):
    out_path = tmp_path / "design_recommendations.txt"
    result = build_report(pareto_rows=_pareto_rows(), out_path=str(out_path))
    knee = result["field_flow"]["knee_point"]
    assert str(knee["mu0H_max_T"]) in result["text"]
    assert str(knee["frequency_Hz"]) in result["text"]