import numpy as np
import pytest

from core.mce_material import GADOLINIUM
from core.amr_cycle import AMRSystem
from core.loss_model import StateDependentLossModel


def make_system(**overrides):
    kwargs = dict(material=GADOLINIUM, mu0H_max=2.0, mass_regenerator=5.0,
                  frequency=1.0, fluid_mdot=0.1, regenerator_effectiveness=0.8)
    kwargs.update(overrides)
    return AMRSystem(**kwargs)


def test_cooling_capacity_nonnegative():
    sys_ = make_system()
    Qc, dTad = sys_.cooling_capacity(T_cold=291.0, T_span=10.0)
    assert Qc >= 0.0


def test_cooling_capacity_zero_for_very_large_span():
    """Qc must fall to 0 once the imposed span swamps the no-load dTad."""
    sys_ = make_system()
    Qc_at_limit, _ = sys_.cooling_capacity(T_cold=291.0, T_span=200.0)
    assert Qc_at_limit == pytest.approx(0.0, abs=1e-6)


def test_electrical_cop_never_exceeds_ideal_cop():
    """Adding parasitic losses can only reduce COP, never increase it."""
    sys_ = make_system()
    result = sys_.run(T_cold=291.0, T_span=10.0)
    assert result.COP_electrical <= result.COP + 1e-9


def test_exergy_efficiency_between_zero_and_one():
    sys_ = make_system()
    result = sys_.run(T_cold=291.0, T_span=10.0)
    assert 0.0 <= result.exergy_eff <= 1.0 + 1e-9


def test_state_dependent_loss_model_increases_with_frequency_and_field():
    lm = StateDependentLossModel()
    low = lm.parasitic_power(frequency=0.5, mu0H=1.0, mdot=0.1, Qc=100.0)
    high = lm.parasitic_power(frequency=4.0, mu0H=3.0, mdot=0.1, Qc=100.0)
    assert high > low


def test_ntu_thermal_model_lets_mass_affect_cooling_capacity():
    """With use_ntu_thermal_model=True, more regenerator mass should improve
    (or at least not worsen) effectiveness / cooling capacity, unlike the
    constant-effectiveness model where mass has zero effect."""
    small = make_system(mass_regenerator=1.0, use_ntu_thermal_model=True)
    large = make_system(mass_regenerator=10.0, use_ntu_thermal_model=True)
    Qc_small, _ = small.cooling_capacity(291.0, 10.0)
    Qc_large, _ = large.cooling_capacity(291.0, 10.0)
    assert Qc_large >= Qc_small


def test_amr_run_is_reasonably_fast():
    """Regression guard for the Newton-solver performance fix: 100
    sequential AMR evaluations (used heavily by optimize.py / sensitivity.py
    / rsm.py) should take well under a second, not the ~6s the old damped
    fixed-point solver needed."""
    import time
    lm = StateDependentLossModel()
    t0 = time.time()
    for _ in range(100):
        sys_ = make_system(loss_model=lm, use_ntu_thermal_model=True)
        sys_.run(291.0, 10.0)
    elapsed = time.time() - t0
    assert elapsed < 2.0, f"AMR evaluations took {elapsed:.2f}s for 100 calls (expected <2s)"
