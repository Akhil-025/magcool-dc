from core.elastocaloric_cycle import (
    run_cycle,
    run_span_sweep,
    calibrate_parasitic_fraction_to_literature,
    MEASURED_BENCHMARKS,
)


def test_feasible_across_full_ashrae_span_range():
    """Unlike this repo's own AMR model (structurally capped well before
    20K, see regenerator_1d.py / the 0-D span cap), this model should stay
    feasible across the full 5-20K range given NiTi's much larger
    achievable Delta_T_ad -- this is a real, notable difference worth
    locking in a test for."""
    rows, _frac = run_span_sweep(spans_K=(5, 10, 15, 20))
    assert all(r.feasible for r in rows)


def test_cop_electrical_never_exceeds_cop_ideal():
    """Adding actuator parasitic work can only ever hurt COP, matching
    amr_cycle.py's own Qc/(W+W_parasitic) <= Qc/W invariant."""
    rows, _frac = run_span_sweep(spans_K=(5, 10, 15, 20))
    for r in rows:
        assert r.COP_electrical <= r.COP_ideal + 1e-9


def test_calibration_reproduces_measured_benchmarks_reasonably():
    frac = calibrate_parasitic_fraction_to_literature(verbose=False)
    for span, cop_lit, _src in MEASURED_BENCHMARKS:
        r = run_cycle(span, parasitic_fraction_of_Qc=frac)
        # Loose tolerance -- this is a 0-D calibration against two points,
        # not a tight per-device fit (see module docstring).
        assert abs(r.COP_electrical - cop_lit) < 1.0


def test_calibrated_model_stays_below_typical_vcc_cop_range():
    """Locks in this module's headline honest finding: calibrated against
    MEASURED literature benchmarks, elastocaloric COP_electrical stays
    well below the 6-24 VCC COP range this repo's own comparison_table.csv
    reports across the same 5-20K span."""
    rows, _frac = run_span_sweep(spans_K=(5, 10, 15, 20))
    for r in rows:
        assert r.COP_electrical < 6.0
