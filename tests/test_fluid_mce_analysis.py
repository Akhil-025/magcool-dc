import os

from core.fluid_mce_analysis import (
    volume_fraction_sweep,
    fixed_span_comparison,
    compare_to_solid_amr_and_liquid_cooling,
    run_fluid_mce_analysis,
    literature_search_for_working_body_benchmark,
    FERROFLUID_LITERATURE_SEARCH_2024_2025,
)


def test_volume_fraction_sweep_returns_all_rows():
    result = volume_fraction_sweep(phis=(0.01, 0.05, 0.1, 0.2, 0.3))
    assert len(result["rows"]) == 5
    for row in result["rows"]:
        assert row["dTad_suspension_K"] >= 0.0
        assert row["span_K"] >= 0.0


def test_volume_fraction_sweep_span_grows_with_phi():
    """More particle loading -> less dilution -> larger achievable
    self-consistent span (module docstring physics item 1)."""
    result = volume_fraction_sweep(phis=(0.02, 0.1, 0.3))
    spans = [row["span_K"] for row in result["rows"]]
    assert spans[0] < spans[1] < spans[2]


def test_fixed_span_comparison_holds_span_equal_for_all_technologies():
    """The PRIMARY comparison: span is an input, not derived from either
    system's own favorable point, so it must come back unchanged."""
    comp = fixed_span_comparison(span_K=10.0)
    assert comp["span_K"] == 10.0
    assert comp["liquid_cooling"]["COP"] > 0.0
    assert comp["vapor_compression"]["COP"] > 0.0


def test_fixed_span_comparison_shows_ferrofluid_span_collapse():
    """Headline finding this module documents: at a realistic, externally
    -imposed 10K span, the ferrofluid system's own achievable
    dTad_suspension can't cover it (infeasible / Qc clipped to 0) while
    the solid AMR, whose regenerator amplifies span beyond a single
    stage's own dTad, remains feasible."""
    comp = fixed_span_comparison(span_K=10.0)
    assert comp["fluid_MCE"]["feasible"] is False
    assert comp["fluid_MCE"]["Qc_W"] == 0.0
    assert comp["solid_AMR"]["feasible"] is True
    assert comp["solid_AMR"]["Qc_W"] > 0.0


def test_compare_returns_all_technologies():
    comp = compare_to_solid_amr_and_liquid_cooling()
    assert comp["fluid_MCE"]["Qc_W"] >= 0.0
    assert comp["solid_AMR"]["Qc_W"] >= 0.0
    assert comp["liquid_cooling"]["COP"] > 0.0
    assert comp["vapor_compression"]["COP"] > 0.0
    assert comp["span_K"] > 0.0


def test_compare_solid_amr_delivers_more_cooling_capacity():
    """At the fluid system's own tiny favorable span, a solid AMR bed
    (with much more thermal mass and NTU-driven heat exchange) should
    still deliver far more absolute cooling power -- the headline
    'span collapses' finding this module's own writeup states."""
    comp = compare_to_solid_amr_and_liquid_cooling()
    assert comp["solid_AMR"]["Qc_W"] > comp["fluid_MCE"]["Qc_W"]


def test_run_fluid_mce_analysis_writes_file(tmp_path):
    out_path = str(tmp_path / "fluid_mce_analysis.txt")
    result = run_fluid_mce_analysis(out_path=out_path, verbose=False)
    assert os.path.exists(out_path)
    assert "PHASE 20" in result["text"]
    assert "sweep" in result and "comparison" in result and "fixed_span_comparison" in result
    assert "PRIMARY" in result["text"]


def test_run_fluid_mce_analysis_no_file_write_when_out_path_none():
    result = run_fluid_mce_analysis(out_path=None, verbose=False)
    assert "PHASE 20" in result["text"]


# ---- Literature search for a working-body benchmark (VALIDATION UPDATE) ----

def test_literature_search_reports_no_working_body_benchmark():
    finding = literature_search_for_working_body_benchmark(verbose=False)
    assert finding["working_body_benchmark_found"] is False
    assert finding["every_paper_uses_ferrofluid_as"] == "thermal switch / heat-transfer medium"
    assert len(finding["papers_checked"]) == len(FERROFLUID_LITERATURE_SEARCH_2024_2025)


def test_literature_search_reports_span_collapse_cross_check():
    finding = literature_search_for_working_body_benchmark(verbose=False)
    assert finding["this_modules_own_max_span_over_phi_sweep_K"] >= 0.0
    assert finding["representative_target_span_K"] > 0.0
    assert 0.0 <= finding["max_span_as_fraction_of_representative"]
    # not asserting the specific boolean value of
    # span_collapse_consistent_with_literature_absence here -- that is a
    # genuine finding, not a fixed property of the interface -- only that
    # the field is present and derived from the two span numbers above.
    assert finding["span_collapse_consistent_with_literature_absence"] == (
        finding["max_span_as_fraction_of_representative"] < 0.25)


def test_literature_search_includes_krieger_dougherty_grounding():
    finding = literature_search_for_working_body_benchmark(verbose=False)
    grounding = finding["krieger_dougherty_rheology_grounding"]
    assert "Susan-Resiga" in grounding["source"]
    assert grounding["this_module_uses_matching_defaults"] is True


def test_literature_search_papers_all_flagged_non_working_body():
    for paper in FERROFLUID_LITERATURE_SEARCH_2024_2025:
        assert "source" in paper and "ferrofluid_role" in paper
        assert "working body" not in paper["ferrofluid_role"].lower() or \
               "not" in paper["ferrofluid_role"].lower()
