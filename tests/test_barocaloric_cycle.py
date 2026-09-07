from core.barocaloric_cycle import (
    run_cycle,
    run_span_sweep,
    calibrate_parasitic_fraction_to_simulation,
    SIMULATED_CALIBRATION_TARGET,
)


def test_feasible_across_full_ashrae_span_range():
    rows, _frac = run_span_sweep(spans_K=(5, 10, 15, 20))
    assert all(r.feasible for r in rows)


def test_cop_electrical_never_exceeds_cop_ideal():
    rows, _frac = run_span_sweep(spans_K=(5, 10, 15, 20))
    for r in rows:
        assert r.COP_electrical <= r.COP_ideal + 1e-9


def test_calibration_matches_simulated_target_closely():
    """Only one calibration point exists here (see module docstring: no
    measured barocaloric device COP was found), so this should match
    almost exactly rather than in a least-squares sense across multiple
    points."""
    span, cop_target, _src = SIMULATED_CALIBRATION_TARGET
    frac = calibrate_parasitic_fraction_to_simulation(verbose=False)
    r = run_cycle(span, parasitic_fraction_of_Qc=frac)
    assert abs(r.COP_electrical - cop_target) < 0.05


def test_calibrated_model_stays_below_typical_vcc_cop_range():
    """Locks in this module's headline honest finding, with the weaker
    (simulation-calibrated, not measurement-calibrated) confidence its own
    docstring describes."""
    rows, _frac = run_span_sweep(spans_K=(5, 10, 15, 20))
    for r in rows:
        assert r.COP_electrical < 6.5
