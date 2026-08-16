"""
Tests for the Phase 28 addition: AMRSystem's pump_motor_efficiency
parameter and its effect on _geometry_pumping_power_W().
"""
import pytest

from core.amr_cycle import AMRSystem
from core.mce_material import GADOLINIUM
from core.loss_model import StateDependentLossModel


def _make_system(**kwargs):
    defaults = dict(material=GADOLINIUM, mu0H_max=2.0, mass_regenerator=2.0,
                     frequency=2.0, fluid_mdot=0.1, particle_diameter=0.0005,
                     use_ntu_thermal_model=True, loss_model=StateDependentLossModel())
    defaults.update(kwargs)
    return AMRSystem(**defaults)


def test_default_pump_motor_efficiency_is_one():
    """Default must be 1.0 (idealized, no loss) -- preserves every
    pre-Phase-28 caller's behavior exactly, since particle_diameter is
    always set by optimize.py's NSGA-III search."""
    sys_ = _make_system()
    assert sys_.pump_motor_efficiency == 1.0


def test_default_behavior_unchanged_from_pre_phase28():
    """A system built with the default (no pump_motor_efficiency passed)
    must give byte-for-byte the same geometry pumping power as one
    explicitly passing pump_motor_efficiency=1.0."""
    sys_default = _make_system()
    sys_explicit = _make_system(pump_motor_efficiency=1.0)
    assert (sys_default._geometry_pumping_power_W()
            == sys_explicit._geometry_pumping_power_W())


def test_literature_efficiency_constant_value():
    assert AMRSystem.PUMP_MOTOR_EFFICIENCY_LITERATURE == pytest.approx(0.6)


def test_opting_into_literature_efficiency_increases_pumping_power():
    sys_ideal = _make_system(pump_motor_efficiency=1.0)
    sys_real = _make_system(pump_motor_efficiency=AMRSystem.PUMP_MOTOR_EFFICIENCY_LITERATURE)
    P_ideal = sys_ideal._geometry_pumping_power_W()
    P_real = sys_real._geometry_pumping_power_W()
    assert P_real == pytest.approx(P_ideal / 0.6)
    assert P_real > P_ideal


def test_opting_into_literature_efficiency_decreases_cop():
    sys_ideal = _make_system(pump_motor_efficiency=1.0)
    sys_real = _make_system(pump_motor_efficiency=AMRSystem.PUMP_MOTOR_EFFICIENCY_LITERATURE)
    r_ideal = sys_ideal.run(290.0, 10.0)
    r_real = sys_real.run(290.0, 10.0)
    assert r_real.COP_electrical < r_ideal.COP_electrical
    # Qc and magnetic work are unaffected -- only the parasitic pumping
    # term changes.
    assert r_real.Qc == pytest.approx(r_ideal.Qc)
    assert r_real.W_mag == pytest.approx(r_ideal.W_mag)


def test_no_particle_diameter_unaffected_by_pump_motor_efficiency():
    """When particle_diameter is None, _geometry_pumping_power_W() returns
    None regardless of pump_motor_efficiency -- the CORE-calibrated
    generic k_pump*mdot**2 term is used unchanged (see this function's own
    docstring: applying efficiency there would double-count losses already
    baked into that fitted coefficient)."""
    sys_ = AMRSystem(material=GADOLINIUM, mu0H_max=2.0, mass_regenerator=2.0,
                      frequency=2.0, fluid_mdot=0.1, use_ntu_thermal_model=True,
                      loss_model=StateDependentLossModel(),
                      pump_motor_efficiency=AMRSystem.PUMP_MOTOR_EFFICIENCY_LITERATURE)
    assert sys_._geometry_pumping_power_W() is None
    r_default_eff = sys_.run(290.0, 10.0)
    sys_ideal = AMRSystem(material=GADOLINIUM, mu0H_max=2.0, mass_regenerator=2.0,
                           frequency=2.0, fluid_mdot=0.1, use_ntu_thermal_model=True,
                           loss_model=StateDependentLossModel(), pump_motor_efficiency=1.0)
    r_ideal = sys_ideal.run(290.0, 10.0)
    assert r_default_eff.COP_electrical == pytest.approx(r_ideal.COP_electrical)
