"""
Smoke test for core/optimize.py's NSGA-III design search.

Uses a small population/generation count so the test runs in ~1s (this is
only possible at all because of the Newton-solver perf fix -- with the
original damped fixed-point solver even this reduced run was too slow to
use in a test suite). Regression-guards against the degenerate outcome
found while investigating the entropy bug, where every "Pareto-optimal"
design collapsed to a single point with COP_electrical=0 and Qc=0 because
the T=291K/span=10K operating point sat in an unphysical region of the
model caused by the entropy bug.
"""
from core.optimize import run_optimization


def test_optimizer_finds_nondegenerate_designs(tmp_path):
    rows = run_optimization(pop_size=20, n_gen=10, seed=1,
                             out_csv=str(tmp_path / "pareto_front.csv"))
    assert len(rows) > 1, "optimizer collapsed to a single design"
    assert any(r["Qc_W"] > 0 for r in rows), "no design has nonzero cooling capacity"
    assert any(r["COP_electrical"] > 1.0 for r in rows), "no design has a sane COP"
