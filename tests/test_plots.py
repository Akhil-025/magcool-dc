"""
Smoke tests for core/plots.py.

core/plots.py has no dedicated tests today -- each figure function is only
exercised by running `python plots.py` (or main.py) manually. These tests
check that every figure function in the public figure set runs without
raising and writes a non-empty PNG + PDF pair, which would catch e.g. an
import-time typo, a signature mismatch with an upstream module, or a
function that silently stopped calling save().

Figure output is redirected to a temporary directory (monkeypatching
plots.FIG_DIR / plots.RESULTS_DIR) so running this suite does not
overwrite the checked-in results/figures/ artifacts.
"""
import pytest

from core import plots

# (figure function, output basename passed to plots.save()) for every
# figure in plots.run_all()'s figure_fns list, kept in sync with that list.
FIGURE_FUNCTIONS = [
    (plots.plot_gd_validation, "fig01_gd_mce_validation"),
    (plots.plot_gd_entropy_dTad, "fig02_gd_entropy_and_dTad_vs_T"),
    (plots.plot_landau_giant_mce, "fig03_landau_giant_mce_calibration"),
    (plots.plot_giguere_validation, "fig04_giguere_direct_vs_indirect_validation"),
    (plots.plot_material_comparison, "fig05_material_comparison_dTad"),
    (plots.plot_amr_characteristic_curve, "fig06_amr_characteristic_curve"),
    (plots.plot_amr_energy_balance, "fig07_amr_energy_balance_vs_span"),
    (plots.plot_amr_vs_baselines, "fig08_amr_vs_baselines_cop"),
    (plots.plot_regenerator_effectiveness, "fig09_regenerator_effectiveness_ntu"),
    (plots.plot_geometry_packed_bed, "fig10_geometry_optimum_packed_bed"),
    (plots.plot_geometry_parallel_plate, "fig11_geometry_optimum_parallel_plate"),
    (plots.plot_loss_model_calibration, "fig12_loss_model_calibration_fit"),
    (plots.plot_parasitic_fraction_scaling, "fig13_parasitic_fraction_scaling"),
    (plots.plot_system_validation, "fig14_system_validation_scatter"),
    (plots.plot_curve_validation, "fig15_curve_validation_companion"),
    (plots.plot_sobol_sensitivity, "fig16_sobol_sensitivity_comparison"),
    (plots.plot_rsm_surrogate, "fig17_rsm_surrogate_parity"),
    (plots.plot_nsga3_pareto, "fig18_nsga3_pareto_front"),
    (plots.plot_cascade_staging_gd, "fig19_cascade_staging_gd"),
    (plots.plot_cascade_giant_vs_gd, "fig20_cascade_giant_mce_vs_gd"),
    (plots.plot_graded_cascade, "fig21_graded_cascade_performance"),
    (plots.plot_economics, "fig22_economics_tco_comparison"),
    (plots.plot_emissions, "fig23_emissions_comparison"),
    (plots.plot_giant_mce_targeting, "fig24_giant_mce_targeting_comparison"),
    (plots.plot_astronautics_validation, "fig25_astronautics_graded_bed_validation"),
    (plots.plot_material_family_comparison, "fig26_material_family_comparison"),
    (plots.plot_inhomogeneous_broadening, "fig27_inhomogeneous_tc_broadening"),
    (plots.plot_nanocomposite_robustness, "fig28_nanocomposite_offdesign_robustness"),
    (plots.plot_thermal_diode_sensitivity, "fig29_thermal_diode_sensitivity"),
    (plots.plot_fluid_mce_sweep, "fig30_fluid_mce_volume_fraction"),
    (plots.plot_passive_regenerator_alignment, "fig31_passive_regenerator_alignment"),
    (plots.plot_cycle_type_validation, "fig32_cycle_type_validation"),
    (plots.plot_hysteresis_sensitivity, "fig33_hysteresis_pareto_sensitivity"),
    (plots.plot_magnet_geometry_pareto_sensitivity, "fig34_magnet_geometry_pareto_sensitivity"),
]


@pytest.fixture(autouse=True)
def _redirect_plot_output(tmp_path, monkeypatch):
    """Point plots.RESULTS_DIR / plots.FIG_DIR at a scratch directory for
    the duration of each test, so figure/data generation here never
    touches the real results/ tree."""
    results_dir = tmp_path / "results"
    fig_dir = results_dir / "figures"
    fig_dir.mkdir(parents=True, exist_ok=True)
    monkeypatch.setattr(plots, "RESULTS_DIR", results_dir)
    monkeypatch.setattr(plots, "FIG_DIR", fig_dir)
    yield fig_dir


@pytest.mark.parametrize("fn,basename", FIGURE_FUNCTIONS,
                          ids=[basename for _, basename in FIGURE_FUNCTIONS])
def test_figure_function_runs_and_saves_png_and_pdf(fn, basename, _redirect_plot_output):
    fig_dir = _redirect_plot_output
    fn()
    png_path = fig_dir / f"{basename}.png"
    pdf_path = fig_dir / f"{basename}.pdf"
    assert png_path.exists(), f"{fn.__name__} did not write {png_path.name}"
    assert pdf_path.exists(), f"{fn.__name__} did not write {pdf_path.name}"
    assert png_path.stat().st_size > 0
    assert pdf_path.stat().st_size > 0


def test_figure_set_covers_every_function_in_run_all():
    """Keeps FIGURE_FUNCTIONS in sync with run_all()'s own figure list, so
    a newly added figure doesn't silently go untested (and a removed one
    doesn't leave a stale entry here)."""
    import inspect
    import re

    source = inspect.getsource(plots.run_all)
    run_all_names = set(re.findall(r"\bplot_\w+", source))
    tested_names = {fn.__name__ for fn, _ in FIGURE_FUNCTIONS}

    missing_from_tests = run_all_names - tested_names
    stale_in_tests = tested_names - run_all_names
    assert not missing_from_tests, (
        f"run_all() lists figure function(s) not covered here: {missing_from_tests}"
    )
    assert not stale_in_tests, (
        f"FIGURE_FUNCTIONS lists function(s) no longer in run_all(): {stale_in_tests}"
    )


def test_plot_material_family_comparison_color_count_matches_candidates(_redirect_plot_output):
    """Phase 22 item 2 added a 6th material-family candidate (the
    nanocomposite blend). plot_material_family_comparison()'s colors6
    list must cover at least as many entries as build_comparison_table()
    returns for the representative span, or matplotlib's bar(color=...)
    raises a length-mismatch error -- this is the regression guard
    core/plots.py's own comment above colors6 points to."""
    from core import material_family_comparison

    rows = material_family_comparison.build_comparison_table()
    rep = [r for r in rows
           if r['span_K'] == material_family_comparison.REPRESENTATIVE_SPAN_K]
    assert len(rep) == 6

    # Must not raise (this is the actual regression check: an earlier
    # length mismatch here raises inside matplotlib's bar()).
    plots.plot_material_family_comparison()