from core.geometry_analysis import (
    demonstrate_earlier_model_had_no_optimum,
    check_free_mdot_cop_is_degenerate,
    sweep_packed_bed_diameter,
    sweep_parallel_plate_spacing,
)


def test_pre_phase7_gap_is_confirmed_monotonic():
    """The whole module exists to close this gap -- verify it's real."""
    _, eps_vals, monotonic = demonstrate_earlier_model_had_no_optimum(verbose=False)
    assert monotonic


def test_free_mdot_cop_maximization_is_degenerate():
    """Confirms the documented reason this module fixes mdot at a
    representative value rather than re-optimizing it per geometry."""
    _, _, monotonic = check_free_mdot_cop_is_degenerate(verbose=False)
    assert monotonic


def test_packed_bed_sweep_has_interior_cop_optimum():
    """With the geometry-coupled pumping-power term, COP_aug should show
    a genuine interior maximum (not pinned at either swept boundary)."""
    rows, best_qc_row, best_cop_row = sweep_packed_bed_diameter(verbose=False)
    diam_values = [r[0] for r in rows]
    cop_optimal_diam = best_cop_row[0]
    assert cop_optimal_diam != max(diam_values)
    assert cop_optimal_diam != min(diam_values)


def test_parallel_plate_sweep_has_interior_cop_optimum():
    rows, best_qc_row, best_cop_row = sweep_parallel_plate_spacing(verbose=False)
    spacing_values = [r[0] for r in rows]
    cop_optimal_spacing = best_cop_row[0]
    assert cop_optimal_spacing != max(spacing_values)
    assert cop_optimal_spacing != min(spacing_values)