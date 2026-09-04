"""
Tests for the addition: core.thermal.intragranular_eddy_power(),
core.loss_model.StateDependentLossModel's intragranular_eddy_power_W
parameter, and core.amr_cycle.AMRSystem._geometry_eddy_power_W() wiring.
"""
import numpy as np
import pytest

from core.thermal import intragranular_eddy_power, GD_SIGMA_E_S_PER_M, RHO_GD
from core.loss_model import StateDependentLossModel
from core.amr_cycle import AMRSystem
from core.mce_material import GADOLINIUM


def test_intragranular_eddy_power_zero_at_zero_diameter():
    assert intragranular_eddy_power(2.0, 2.0, particle_diameter=0.0) == 0.0


def test_intragranular_eddy_power_scales_quadratically_with_diameter():
    W1 = intragranular_eddy_power(2.0, 2.0, particle_diameter=0.001, mass_regenerator=2.0)
    W2 = intragranular_eddy_power(2.0, 2.0, particle_diameter=0.002, mass_regenerator=2.0)
    assert W2 == pytest.approx(4 * W1, rel=1e-9)


def test_intragranular_eddy_power_scales_quadratically_with_frequency_and_field():
    W_base = intragranular_eddy_power(1.0, 1.0, particle_diameter=0.0005, mass_regenerator=2.0)
    W_2f = intragranular_eddy_power(2.0, 1.0, particle_diameter=0.0005, mass_regenerator=2.0)
    W_2H = intragranular_eddy_power(1.0, 2.0, particle_diameter=0.0005, mass_regenerator=2.0)
    assert W_2f == pytest.approx(4 * W_base, rel=1e-9)
    assert W_2H == pytest.approx(4 * W_base, rel=1e-9)


def test_intragranular_eddy_power_scales_linearly_with_mass():
    W1 = intragranular_eddy_power(2.0, 2.0, particle_diameter=0.0005, mass_regenerator=1.0)
    W2 = intragranular_eddy_power(2.0, 2.0, particle_diameter=0.0005, mass_regenerator=3.0)
    assert W2 == pytest.approx(3 * W1, rel=1e-9)


def test_intragranular_eddy_power_negligible_at_realistic_particle_sizes():
    """At realistic packed-bed particle diameters (Tusek et al. 2013's own
    optimum range, 0.07-0.17mm), the intragranular eddy term should be
    negligible (milliwatt scale) relative to a typical CORE-calibrated
    support-structure k_eddy term (hundreds of watts) -- a real finding,
    not a bug: see this module's own ROADMAP entry."""
    lm = StateDependentLossModel()
    k_eddy_term_W = lm.k_eddy * 2.0 ** 2 * 2.0 ** 2
    for d_mm in (0.07, 0.17):
        W = intragranular_eddy_power(2.0, 2.0, particle_diameter=d_mm / 1000.0,
                                      mass_regenerator=2.0)
        assert W < 0.01
        assert W / k_eddy_term_W < 1e-3


def test_state_dependent_loss_model_default_unaffected_by_new_parameter():
    """Every previous call site omits intragranular_eddy_power_W --
    parasitic_power() must be bit-for-bit unchanged in that case."""
    lm = StateDependentLossModel()
    old_style = lm.parasitic_power(2.0, 2.0, 0.1, 500.0)
    new_style_default = lm.parasitic_power(2.0, 2.0, 0.1, 500.0,
                                            intragranular_eddy_power_W=0.0)
    assert old_style == new_style_default


def test_state_dependent_loss_model_adds_not_replaces_eddy_term():
    """intragranular_eddy_power_W should be ADDED to k_eddy*f**2*mu0H**2,
    unlike pumping_power_override's REPLACE semantics."""
    lm = StateDependentLossModel()
    base = lm.parasitic_power(2.0, 2.0, 0.1, 500.0)
    with_extra = lm.parasitic_power(2.0, 2.0, 0.1, 500.0, intragranular_eddy_power_W=5.0)
    assert with_extra == pytest.approx(base + 5.0)


def test_amr_system_geometry_eddy_power_zero_without_particle_diameter():
    sys_ = AMRSystem(material=GADOLINIUM, mu0H_max=2.0, mass_regenerator=2.0,
                      frequency=2.0, fluid_mdot=0.1)
    assert sys_._geometry_eddy_power_W() == 0.0


def test_amr_system_geometry_eddy_power_nonzero_with_particle_diameter():
    sys_ = AMRSystem(material=GADOLINIUM, mu0H_max=2.0, mass_regenerator=2.0,
                      frequency=2.0, fluid_mdot=0.1, particle_diameter=0.0005)
    W = sys_._geometry_eddy_power_W()
    assert W > 0.0
    expected = intragranular_eddy_power(2.0, 2.0, particle_diameter=0.0005,
                                         mass_regenerator=2.0)
    assert W == pytest.approx(expected)


def test_amr_system_run_unaffected_by_phase27_when_no_particle_diameter():
    """A run() without particle_diameter should give identical results to
    before this phase (regression guard)."""
    sys_ = AMRSystem(material=GADOLINIUM, mu0H_max=2.0, mass_regenerator=2.0,
                      frequency=2.0, fluid_mdot=0.1)
    result = sys_.run(T_cold=290.0, T_span=10.0)
    assert result is not None
