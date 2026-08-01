"""
Unit tests for core/giant_mce_analysis.py.

Previously only exercised via main.py's integration run. These tests check
the two claims the module's docstring and run_analysis() output rest on:
(1) Gd5Si2Ge2's own peak-effect temperature sits below the ASHRAE range
while Gd's sits inside it, and (2) each material performs far better when
operated at its own favorable point than when forced onto the other
material's operating point.
"""
from core.giant_mce_analysis import find_peak_temperature, run_analysis
from core.mce_material import GADOLINIUM
from core.first_order_mce import GD5SI2GE2_FIRST_ORDER
from core.amr_cycle import AMRSystem
from core.loss_model import StateDependentLossModel

ASHRAE_LOW_K = 291.0
ASHRAE_HIGH_K = 300.0


def test_gd_peak_temperature_is_inside_ashrae_range():
    peak_T_gd = find_peak_temperature(GADOLINIUM, mu0H=2.0)
    assert ASHRAE_LOW_K <= peak_T_gd <= ASHRAE_HIGH_K


def test_gd5si2ge2_peak_temperature_is_below_ashrae_range():
    """This is the central finding the module exists to establish: the
    giant-MCE material's own peak sits below the data-center supply range,
    which is why it isn't directly usable as-is."""
    peak_T_giant = find_peak_temperature(GD5SI2GE2_FIRST_ORDER, mu0H=2.0)
    assert peak_T_giant < ASHRAE_LOW_K


def test_gd5si2ge2_collapses_when_forced_onto_ashrae_point():
    """At the ASHRAE operating point, Gd5Si2Ge2's narrow first-order
    transition window is far from where it's centered, so cooling
    capacity should be much smaller than at its own favorable point."""
    loss_model = StateDependentLossModel()
    span = 10.0
    sys_giant = AMRSystem(material=GD5SI2GE2_FIRST_ORDER, mu0H_max=2.0,
                           mass_regenerator=5.0, frequency=1.0, fluid_mdot=0.08,
                           loss_model=loss_model, use_ntu_thermal_model=True)
    r_ashrae = sys_giant.run(ASHRAE_LOW_K, span)

    peak_T_giant = find_peak_temperature(GD5SI2GE2_FIRST_ORDER, mu0H=2.0)
    r_own_point = sys_giant.run(peak_T_giant - span / 2, span)

    assert r_ashrae.Qc < r_own_point.Qc


def test_gd_underperforms_at_giant_mce_favorable_point():
    """Symmetric claim: Gd, whose own peak is near the ASHRAE range,
    should do poorly if forced onto Gd5Si2Ge2's (much colder) favorable
    point."""
    loss_model = StateDependentLossModel()
    span = 10.0
    peak_T_giant = find_peak_temperature(GD5SI2GE2_FIRST_ORDER, mu0H=2.0)

    sys_gd = AMRSystem(material=GADOLINIUM, mu0H_max=2.0, mass_regenerator=5.0,
                        frequency=1.0, fluid_mdot=0.08, loss_model=loss_model,
                        use_ntu_thermal_model=True)
    r_gd_at_ashrae = sys_gd.run(ASHRAE_LOW_K, span)
    r_gd_at_giant_point = sys_gd.run(peak_T_giant - span / 2, span)

    assert r_gd_at_ashrae.Qc > r_gd_at_giant_point.Qc


def test_run_analysis_writes_output_file(tmp_path):
    out_path = tmp_path / "giant_mce_analysis.txt"
    run_analysis(out_path=str(out_path))
    assert out_path.exists()
    text = out_path.read_text()
    assert "Gd5Si2Ge2" in text
    assert "ASHRAE" in text
