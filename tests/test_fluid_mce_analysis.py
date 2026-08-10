import os

from core.fluid_mce_analysis import (
    volume_fraction_sweep,
    compare_to_solid_amr_and_liquid_cooling,
    run_fluid_mce_analysis,
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
    assert "sweep" in result and "comparison" in result


def test_run_fluid_mce_analysis_no_file_write_when_out_path_none():
    result = run_fluid_mce_analysis(out_path=None, verbose=False)
    assert "PHASE 20" in result["text"]
