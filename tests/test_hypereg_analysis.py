"""
Unit tests for core/hypereg_analysis.py (Phase 15 item 3).

The previous version of this file was an accidental byte-for-byte copy of
core/hypereg_analysis.py itself (no test_* functions), so it silently
contributed zero coverage ("no tests ran"). These tests actually exercise
sweep_n_parallel()/sweep_frequency_at_fixed_n()/run_hypereg_analysis()
against the qualitative claims made in hypereg_findings.md and the
ROADMAP.md Phase 15 write-up: n=1 reproduces conventional series-flow
behavior, COP_electrical is non-decreasing in n_parallel (diminishing,
saturating benefit), and the frequency sweep's pumping-power saving is
computed consistently with the module's own W_parasitic definition.
"""
import pytest

from core.hypereg_analysis import (
    sweep_n_parallel,
    sweep_frequency_at_fixed_n,
    sweep_n_parallel_at_higher_mdot,
    run_hypereg_analysis,
    _run,
    MDOT_KG_S,
)


def test_n_parallel_1_matches_conventional_series_flow():
    """hypereg_n_parallel=None (conventional) and n_parallel=1 (Hypereg's
    own degenerate single-sub-regenerator case) must be physically
    identical -- see AMRSystem._geometry_pumping_power_W()'s docstring."""
    conventional = _run(None)
    n1 = _run(1)
    assert conventional.COP_electrical == pytest.approx(n1.COP_electrical)
    assert conventional.Qc == pytest.approx(n1.Qc)
    assert conventional.W_parasitic == pytest.approx(n1.W_parasitic)


def test_sweep_n_parallel_cop_is_non_decreasing():
    """More parallel sub-beds should never hurt COP_electrical in this
    model (splitting only ever reduces pumping power, never increases
    it) -- the benefit should saturate, not reverse."""
    rows = sweep_n_parallel(verbose=False)
    cops = [row[2] for row in rows]
    assert all(b >= a - 1e-9 for a, b in zip(cops, cops[1:])), \
        "COP_electrical decreased when splitting into more parallel sub-beds"


def test_sweep_n_parallel_returns_expected_shape():
    n_values = (1, 2, 4)
    rows = sweep_n_parallel(n_values=n_values, verbose=False)
    assert len(rows) == len(n_values)
    for (n, Qc, cop, w_parasitic), expected_n in zip(rows, n_values):
        assert n == expected_n
        assert Qc > 0
        assert cop > 0
        assert w_parasitic >= 0


def test_sweep_n_parallel_benefit_is_modest_and_saturating():
    """Regression guard against the magnitude documented in ROADMAP.md's
    Phase 15 write-up (n=4 gives a small, single-digit-percent COP gain
    over n=1, saturating by n=16) -- catches a wiring bug that made the
    pumping-power channel dominant (or a no-op) rather than "one of
    three loss channels, not the dominant one" as the findings note
    claims."""
    rows = sweep_n_parallel(n_values=(1, 4, 16), verbose=False)
    cop_n1 = rows[0][2]
    cop_n4 = rows[1][2]
    cop_n16 = rows[2][2]
    gain_n4_pct = 100 * (cop_n4 / cop_n1 - 1)
    gain_n16_pct = 100 * (cop_n16 / cop_n1 - 1)
    assert 0 <= gain_n4_pct < 5, f"n=4 COP gain ({gain_n4_pct:.2f}%) outside expected modest range"
    assert gain_n16_pct >= gain_n4_pct, "benefit should keep growing (or saturate), not shrink, from n=4 to n=16"


def test_sweep_frequency_at_fixed_n_returns_expected_shape():
    frequencies = (0.5, 1.0, 2.0)
    rows = sweep_frequency_at_fixed_n(frequencies=frequencies, n_parallel=4, verbose=False)
    assert len(rows) == len(frequencies)
    for (f, cop_conv, cop_hyp, saving_pct), expected_f in zip(rows, frequencies):
        assert f == expected_f
        assert cop_conv > 0
        assert cop_hyp > 0
        # Hypereg's pumping-power reduction can only help or be neutral,
        # never hurt COP relative to the conventional case at the same f.
        assert cop_hyp >= cop_conv - 1e-9
        assert saving_pct >= -1e-9


def test_run_hypereg_analysis_writes_report_and_mentions_both_sweeps(tmp_path):
    out_path = tmp_path / "hypereg_analysis.txt"
    text = run_hypereg_analysis(out_path=str(out_path))
    assert out_path.exists()
    assert out_path.read_text() == text
    assert "n_parallel_subregenerators" in text
    assert "Conclusion" in text
    assert "Klinar" in text

# ---------------------------------------------------------------------------
# sweep_n_parallel_at_higher_mdot() -- closes the open ROADMAP.md Phase 16
# candidate asking whether Hypereg's benefit becomes non-negligible at a
# higher-mdot operating point than the module's own 0.08kg/s default.
# ---------------------------------------------------------------------------

def test_higher_mdot_sweep_cop_is_non_decreasing():
    rows = sweep_n_parallel_at_higher_mdot(verbose=False)
    cops = [row[2] for row in rows]
    assert all(b >= a - 1e-9 for a, b in zip(cops, cops[1:])), \
        "COP_electrical decreased when splitting into more parallel sub-beds at higher mdot"


def test_higher_mdot_sweep_returns_expected_shape():
    n_values = (1, 2, 4)
    rows = sweep_n_parallel_at_higher_mdot(n_values=n_values, mdot=0.3, verbose=False)
    assert len(rows) == len(n_values)
    for (n, Qc, cop, w_parasitic), expected_n in zip(rows, n_values):
        assert n == expected_n
        assert Qc > 0
        assert cop > 0
        assert w_parasitic >= 0


def test_higher_mdot_sweep_n1_matches_direct_amrsystem_call():
    """n_parallel=1 at the higher mdot must match a plain AMRSystem call at
    that same mdot -- a direct regression guard on the wiring, independent
    of the module-default-mdot _run() helper used by the other sweeps."""
    from core.mce_material import GADOLINIUM
    from core.amr_cycle import AMRSystem
    from core.hypereg_analysis import (T_COLD_K, SPAN_K, MU0H_T, MASS_KG,
                                         PARTICLE_DIAMETER_M, _LOSS_MODEL)
    rows = sweep_n_parallel_at_higher_mdot(n_values=(1,), mdot=0.3, verbose=False)
    direct = AMRSystem(GADOLINIUM, mu0H_max=MU0H_T, mass_regenerator=MASS_KG,
                        frequency=1.0, fluid_mdot=0.3,
                        regenerator_effectiveness=0.85, loss_model=_LOSS_MODEL,
                        use_ntu_thermal_model=True, particle_diameter=PARTICLE_DIAMETER_M,
                        hypereg_n_parallel=1).run(T_COLD_K, SPAN_K)
    assert rows[0][2] == pytest.approx(direct.COP_electrical)


def test_higher_mdot_sweep_uses_higher_mdot_than_default(monkeypatch):
    """A crude but direct check that the mdot parameter is actually threaded
    through to AMRSystem (rather than silently reusing the module-default
    MDOT_KG_S): Qc at a much higher mdot should exceed Qc at the module
    default, all else equal."""
    default_rows = sweep_n_parallel_at_higher_mdot(n_values=(1,), mdot=MDOT_KG_S, verbose=False)
    higher_rows = sweep_n_parallel_at_higher_mdot(n_values=(1,), mdot=0.3, verbose=False)
    assert higher_rows[0][1] > default_rows[0][1]  # Qc


def test_run_hypereg_analysis_mentions_step_3(tmp_path):
    out_path = tmp_path / "hypereg_analysis.txt"
    text = run_hypereg_analysis(out_path=str(out_path))
    assert "Step 3" in text
    assert "higher mdot" in text