import math

from core.barocaloric_material import NPG_plastic_crystal


def test_delta_T_ad_at_P_sat_matches_measured_entropy_change():
    """Consistency check only (see module docstring: P_sat was CHOSEN to
    equal the measurement pressure) -- verifies the arithmetic, not an
    independent validation."""
    dT = NPG_plastic_crystal.delta_T_ad(NPG_plastic_crystal.P_sat_MPa)
    expected = (NPG_plastic_crystal.T0 * NPG_plastic_crystal.delta_S_J_per_kgK
                / NPG_plastic_crystal.c_p_J_per_kgK)
    assert math.isclose(dT, expected, rel_tol=1e-9)


def test_transformed_fraction_clamped_to_unit_interval():
    assert NPG_plastic_crystal.transformed_fraction(-5.0) == 0.0
    assert NPG_plastic_crystal.transformed_fraction(1e6) == 1.0


def test_hysteresis_loss_nonnegative():
    assert NPG_plastic_crystal.hysteresis_loss_J_per_kg(20.0) >= 0.0
    assert NPG_plastic_crystal.hysteresis_loss_J_per_kg(0.0) == 0.0
