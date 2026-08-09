import pytest

from core.thermal_diode import (MechanicalContactDiode,
                                  DEFAULT_MECHANICAL_CONTACT_DIODE,
                                  cycle_time_reduction_factor)


def test_default_diode_has_rectification_ratio_above_one():
    assert DEFAULT_MECHANICAL_CONTACT_DIODE.rectification_ratio > 1.0


def test_rectification_ratio_is_forward_over_reverse():
    diode = MechanicalContactDiode(forward_conductance_W_K=10.0,
                                    reverse_conductance_W_K=2.0)
    assert diode.rectification_ratio == pytest.approx(5.0)


def test_reverse_conductance_may_not_exceed_forward():
    with pytest.raises(ValueError):
        MechanicalContactDiode(forward_conductance_W_K=1.0,
                                reverse_conductance_W_K=2.0)


def test_conductances_must_be_positive():
    with pytest.raises(ValueError):
        MechanicalContactDiode(forward_conductance_W_K=0.0,
                                reverse_conductance_W_K=0.0)
    with pytest.raises(ValueError):
        MechanicalContactDiode(forward_conductance_W_K=-1.0,
                                reverse_conductance_W_K=0.5)


def test_negative_actuation_energy_rejected():
    with pytest.raises(ValueError):
        MechanicalContactDiode(forward_conductance_W_K=5.0,
                                reverse_conductance_W_K=0.5,
                                actuation_energy_J_per_cycle=-1.0)


def test_switching_power_is_zero_when_actuation_energy_is_zero():
    diode = MechanicalContactDiode(forward_conductance_W_K=5.0,
                                    reverse_conductance_W_K=0.5)
    assert diode.switching_power_W(frequency=3.0) == 0.0


def test_switching_power_scales_linearly_with_frequency():
    diode = MechanicalContactDiode(forward_conductance_W_K=5.0,
                                    reverse_conductance_W_K=0.5,
                                    actuation_energy_J_per_cycle=0.1)
    assert diode.switching_power_W(2.0) == pytest.approx(0.2)
    assert diode.switching_power_W(4.0) == pytest.approx(2 * diode.switching_power_W(2.0))


def test_switching_power_rejects_negative_frequency():
    diode = MechanicalContactDiode(forward_conductance_W_K=5.0,
                                    reverse_conductance_W_K=0.5,
                                    actuation_energy_J_per_cycle=0.1)
    with pytest.raises(ValueError):
        diode.switching_power_W(-1.0)


def test_cycle_time_reduction_factor_basic():
    assert cycle_time_reduction_factor(0.5, 0.05) == pytest.approx(0.9)
    assert cycle_time_reduction_factor(1.0, 1.0) == pytest.approx(0.0)
    assert cycle_time_reduction_factor(1.0, 0.0) == pytest.approx(1.0)


def test_cycle_time_reduction_factor_rejects_invalid_inputs():
    with pytest.raises(ValueError):
        cycle_time_reduction_factor(0.0, 0.0)
    with pytest.raises(ValueError):
        cycle_time_reduction_factor(1.0, -0.1)
    with pytest.raises(ValueError):
        cycle_time_reduction_factor(0.5, 0.6)  # diode slower than conventional