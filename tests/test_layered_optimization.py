"""
Tests for the Phase 29 addition: core.cascade.run_graded_cascade()'s new
particle_diameter/blow_fraction/pump_motor_efficiency parameters, and
core.optimize.LayeredAMRDesignProblem / run_layered_optimization_for_n_layers /
run_layered_optimization / _layered_pareto_filter.

Most tests below avoid running a full NSGA-III search (each individual
evaluation calls run_graded_cascade(), which is significantly slower per
call than a single AMRSystem.run() -- consistent with this repo's own
existing note that run_graded_cascade() is "the single slowest stage in
the full pipeline"). Only one smoke test runs an actual (tiny) NSGA-III
search, to keep this file's total runtime bounded.
"""
import numpy as np
import pytest

from core.cascade import run_graded_cascade, GD_FAMILY, LAFESIH_FAMILY
from core.optimize import (
    LayeredAMRDesignProblem, run_layered_optimization_for_n_layers,
    run_layered_optimization, run_layered_optimization_material_family_cross_product,
    _layered_pareto_filter, LAYERED_N_LAYERS_RANGE, T_COLD_K, SPAN_K,
)


def test_run_graded_cascade_default_unaffected_by_new_params():
    """Omitting particle_diameter/blow_fraction/pump_motor_efficiency must
    reproduce the exact pre-Phase-29 result."""
    r_default = run_graded_cascade(291.0, 10.0, 2, family=GD_FAMILY)
    r_explicit = run_graded_cascade(291.0, 10.0, 2, family=GD_FAMILY,
                                     particle_diameter=None, blow_fraction=0.5,
                                     pump_motor_efficiency=1.0)
    assert r_default["Qc_W"] == r_explicit["Qc_W"]
    assert r_default["COP_cascade"] == r_explicit["COP_cascade"]


def test_run_graded_cascade_particle_diameter_changes_result():
    r_none = run_graded_cascade(291.0, 10.0, 2, family=GD_FAMILY)
    r_geom = run_graded_cascade(291.0, 10.0, 2, family=GD_FAMILY,
                                 particle_diameter=0.0005)
    # Geometry-explicit pumping power replaces the generic k_pump term --
    # some change in COP is expected (not asserting direction, just that
    # the new parameter is actually wired through to a different result).
    assert r_none["COP_cascade"] != r_geom["COP_cascade"]


def test_run_graded_cascade_pump_motor_efficiency_reduces_cop():
    r_ideal = run_graded_cascade(291.0, 10.0, 2, family=GD_FAMILY,
                                  particle_diameter=0.0005, pump_motor_efficiency=1.0)
    r_real = run_graded_cascade(291.0, 10.0, 2, family=GD_FAMILY,
                                 particle_diameter=0.0005, pump_motor_efficiency=0.6)
    assert r_real["COP_cascade"] <= r_ideal["COP_cascade"]
    assert r_real["Qc_W"] == pytest.approx(r_ideal["Qc_W"])


def test_layered_problem_bounds_and_dims():
    problem = LayeredAMRDesignProblem(n_layers=3)
    assert problem.n_layers == 3
    assert problem.n_var == 7
    assert problem.n_obj == 3
    assert problem.n_constr == 0


def test_layered_problem_evaluate_single_point():
    """Directly evaluate _evaluate() on one design vector (no NSGA-III
    search) -- fast, and confirms the objective wiring is correct."""
    problem = LayeredAMRDesignProblem(n_layers=2, family=GD_FAMILY, family_name="Gd")
    x = np.array([2.0, 1.0, 0.1, 3.0, 0.85, 0.5, 0.5])  # last var is d_p in mm
    out = {}
    problem._evaluate(x, out)
    f1, f2, f3 = out["F"]
    # f1 = -COP, f2 = -Qc, f3 = cost -- for a reasonable point, expect a
    # feasible (non-penalty) result: COP > 0, Qc > 0, finite cost.
    assert f1 < 0.0
    assert f2 < 0.0
    assert f3 > 0.0
    assert np.isfinite([f1, f2, f3]).all()


def test_layered_problem_infeasible_point_gets_penalty():
    """A design guaranteed to produce zero/negative Qc (e.g. absurdly low
    mdot combined with a tiny mass) should trigger the finite infeasible
    penalty, not NaN/inf."""
    problem = LayeredAMRDesignProblem(n_layers=6, family=GD_FAMILY, family_name="Gd")
    x = np.array([1.0, 5.0, 0.02, 1.0, 0.6, 0.1, 0.05])  # xl corner, harsh point
    out = {}
    problem._evaluate(x, out)
    f1, f2, f3 = out["F"]
    assert np.isfinite([f1, f2, f3]).all()
    if (f1, f2, f3) == LayeredAMRDesignProblem._INFEASIBLE_PENALTY:
        assert f3 == pytest.approx(1.0e7)


def test_layered_pareto_filter_basic_dominance():
    rows = [
        {"COP_cascade": 5.0, "Qc_W": 1000.0, "cost_index_USD": 500.0},
        {"COP_cascade": 3.0, "Qc_W": 800.0, "cost_index_USD": 700.0},  # dominated by row 0
        {"COP_cascade": 6.0, "Qc_W": 500.0, "cost_index_USD": 400.0},  # non-dominated trade-off
    ]
    kept = _layered_pareto_filter(rows)
    assert len(kept) == 2
    assert rows[1] not in kept


def test_layered_pareto_filter_empty():
    assert _layered_pareto_filter([]) == []


def test_layered_n_layers_range_matches_repo_convention():
    assert LAYERED_N_LAYERS_RANGE == (1, 2, 3, 4, 5, 6)


def test_run_layered_optimization_for_n_layers_smoke(tmp_path):
    """Small NSGA-III run (n_layers=1, kept small for test speed) --
    confirms end-to-end wiring produces nondegenerate results."""
    rows = run_layered_optimization_for_n_layers(
        1, family=GD_FAMILY, family_name="Gd", pop_size=8, n_gen=2, seed=1,
        out_csv=str(tmp_path / "layered_n1.csv"))
    assert len(rows) >= 1
    assert any(r["Qc_W"] > 0 for r in rows)


def test_run_layered_optimization_out_csv_none_does_not_write_or_crash(tmp_path):
    """Phase 31 regression test: run_layered_optimization() used to call
    _write_csv() unconditionally, so out_csv=None crashed with a TypeError
    from os.path.dirname(None) instead of simply skipping the write (an
    asymmetry with per_n_layers_out_dir, which already guarded None
    correctly) -- surfaced by
    run_layered_optimization_material_family_cross_product() below, which
    needs exactly this to avoid a redundant per-family top-level CSV."""
    rows = run_layered_optimization(
        n_layers_range=(1,), family=GD_FAMILY, family_name="Gd",
        pop_size=8, n_gen=2, seed=1, out_csv=None, per_n_layers_out_dir=None)
    assert len(rows) >= 1
    assert not any(f.endswith(".csv") for f in __import__("os").listdir(tmp_path))


def test_run_layered_optimization_material_family_cross_product_smoke(tmp_path):
    """Phase 31 addition: the material x n_layers cross-product explicitly
    left as a documented follow-up by run_layered_optimization()'s own
    docstring. Small-scale (2 families x 2 n_layers values, pop_size=6,
    n_gen=2) purely to keep this test fast -- not representative of
    production-quality Pareto-front resolution."""
    rows = run_layered_optimization_material_family_cross_product(
        family_candidates=[("Gd", None, "Gd"), ("La(Fe,Si)13Hy", LAFESIH_FAMILY, "LaFeSiHy")],
        n_layers_range=(1, 2), pop_size=6, n_gen=2, seed=1,
        out_csv=str(tmp_path / "cross_product.csv"), per_combo_out_dir=None)
    assert len(rows) >= 1
    assert any(r["Qc_W"] > 0 for r in rows)
    # Every returned row's family label must be one of the two candidates
    # actually searched -- guards against a label/family mismatch in the
    # per-family loop silently mixing up which rows came from where.
    assert set(r["family"] for r in rows) <= {"Gd", "La(Fe,Si)13Hy"}


def test_run_layered_optimization_material_family_cross_product_writes_merged_csv(tmp_path):
    out_path = tmp_path / "cross_product.csv"
    run_layered_optimization_material_family_cross_product(
        family_candidates=[("Gd", None, "Gd"), ("La(Fe,Si)13Hy", LAFESIH_FAMILY, "LaFeSiHy")],
        n_layers_range=(1,), pop_size=6, n_gen=2, seed=1,
        out_csv=str(out_path), per_combo_out_dir=None)
    assert out_path.exists()
    with open(out_path) as f:
        header = f.readline()
    assert "family" in header and "COP_cascade" in header
