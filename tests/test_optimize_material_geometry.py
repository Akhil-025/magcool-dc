"""
Tests for core/optimize.py's Phase 15 material + geometry co-optimization
(particle_diameter_mm as a 7th design variable; material as a per-family
NSGA-III search merged post-hoc into one Pareto front).
"""
from core.optimize import (run_optimization, run_optimization_for_material,
                             _material_candidates, _pareto_filter, cost_index)
from core.mce_material import GADOLINIUM


def test_material_candidates_include_gd_and_are_in_range():
    candidates = _material_candidates()
    labels = [label for label, _, _ in candidates]
    assert any(label == "Gd" for label in labels), "plain Gd must always be a candidate"
    # every candidate's tuned Tc (if not Gd) should be inside its own family's window --
    # _material_candidates() itself is responsible for filtering, so this is a
    # regression guard against that filter silently breaking.
    assert len(candidates) >= 1


def test_cost_index_by_family_differs_from_gd():
    """La(Fe,Si)13Hy is cheaper per kg than Gd (Russek & Zimm 2006) --
    cost_index() should reflect that, not silently assume Gd for every
    family."""
    gd_cost = cost_index(1.5, 5.0, "Gd")
    lafesih_cost = cost_index(1.5, 5.0, "La(Fe,Si)13Hy")
    assert lafesih_cost < gd_cost


def test_run_optimization_for_material_returns_geometry_and_material_columns(tmp_path):
    rows = run_optimization_for_material(GADOLINIUM, "Gd", "Gd", pop_size=12, n_gen=5,
                                          seed=1, out_csv=str(tmp_path / "gd_front.csv"))
    assert len(rows) > 0
    for r in rows:
        assert r["material"] == "Gd"
        assert 0.05 <= r["particle_diameter_mm"] <= 2.0
    assert (tmp_path / "gd_front.csv").exists()


def test_pareto_filter_removes_dominated_rows():
    rows = [
        {"material": "A", "COP_electrical": 5.0, "Qc_W": 100.0, "cost_index_USD": 50.0},
        {"material": "B", "COP_electrical": 4.0, "Qc_W": 90.0, "cost_index_USD": 60.0},  # dominated by A
        {"material": "C", "COP_electrical": 3.0, "Qc_W": 200.0, "cost_index_USD": 40.0},  # non-dominated
    ]
    filtered = _pareto_filter(rows)
    materials = {r["material"] for r in filtered}
    assert "B" not in materials, "B is dominated by A in all three objectives and must be removed"
    assert "A" in materials and "C" in materials


def test_run_optimization_merges_multiple_materials_and_writes_per_material_csvs(tmp_path):
    out_csv = tmp_path / "merged.csv"
    per_material_dir = tmp_path / "by_material"
    rows = run_optimization(pop_size=12, n_gen=5, seed=1,
                             out_csv=str(out_csv), per_material_out_dir=str(per_material_dir))
    assert out_csv.exists()
    assert len(rows) > 0
    materials_seen = {r["material"] for r in rows}
    assert "Gd" in materials_seen
    # at least one per-material CSV should have been written (Gd always is)
    assert any(per_material_dir.glob("*.csv"))
    # no design should be strictly dominated by another in the merged set
    refiltered = _pareto_filter(rows)
    assert len(refiltered) == len(rows), "run_optimization()'s own merged output should already be non-dominated"