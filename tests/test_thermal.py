import pytest

from core.thermal import (
    regenerator_effectiveness,
    pressure_drop_packed_bed,
    pumping_power_packed_bed,
    regenerator_effectiveness_parallel_plate,
    pumping_power_parallel_plate,
)


def test_packed_bed_effectiveness_between_zero_and_clip():
    r = regenerator_effectiveness(2.0, frequency=1.0, mdot=0.08)
    assert 0.0 <= r["eps"] <= 0.97 + 1e-9


def test_cp_solid_none_reproduces_default_behavior():
    """ addition: cp_solid=None must reproduce the exact previous
    result (module-level CP_SOLID_GD), for every existing caller."""
    from core.thermal import CP_SOLID_GD
    r_default = regenerator_effectiveness(2.0, 1.0, 0.08)
    r_explicit_none = regenerator_effectiveness(2.0, 1.0, 0.08, cp_solid=None)
    r_explicit_same = regenerator_effectiveness(2.0, 1.0, 0.08, cp_solid=CP_SOLID_GD)
    assert r_default["eps"] == r_explicit_none["eps"]
    assert r_default["eps"] == pytest.approx(r_explicit_same["eps"])


def test_cp_solid_override_changes_utilization_and_effectiveness():
    """A higher solid heat capacity should lower the utilization term U and
    therefore raise (or leave unchanged, if already NTU-capped) eps."""
    r_low_cp = regenerator_effectiveness(2.0, 1.0, 0.08, cp_solid=50.0)
    r_high_cp = regenerator_effectiveness(2.0, 1.0, 0.08, cp_solid=5000.0)
    assert r_high_cp["U"] < r_low_cp["U"]
    assert r_high_cp["eps"] >= r_low_cp["eps"]


def test_packed_bed_pumping_power_increases_as_particle_shrinks():
    """Smaller particles -> smaller hydraulic diameter -> more viscous
    pressure drop at a fixed mdot. This is the physical relationship
    `pumping_power_packed_bed` exists to capture (Tusek et al. 2013,
    Int. J. Refrig. 36, Eqs. 5 & 7)."""
    p_large = pumping_power_packed_bed(0.08, particle_diameter=0.001, mass_regenerator=2.0)
    p_small = pumping_power_packed_bed(0.08, particle_diameter=0.0001, mass_regenerator=2.0)
    assert p_small["P_pump_W"] > p_large["P_pump_W"]


def test_packed_bed_pressure_drop_nonnegative():
    info = pressure_drop_packed_bed(0.08, particle_diameter=0.0005, mass_regenerator=2.0)
    assert info["dP_Pa"] >= 0.0
    assert info["Re"] > 0.0
    assert info["d_h_m"] > 0.0


def test_particle_diameter_alone_has_no_optimum_in_ntu_effectiveness():
    """Documents the previous gap this module's pumping-power addition
    closes: regenerator_effectiveness()'s eps rises monotonically as
    particle_diameter shrinks (no trade-off is representable without a
    coupled pumping-power cost)."""
    diam_m = [0.002, 0.001, 0.0005, 0.00025, 0.0001, 0.00005]
    eps_vals = [regenerator_effectiveness(2.0, 1.0, 0.08, particle_diameter=d)["eps"]
                for d in diam_m]
    assert all(eps_vals[i] <= eps_vals[i + 1] + 1e-9 for i in range(len(eps_vals) - 1))


def test_parallel_plate_effectiveness_between_zero_and_clip():
    r = regenerator_effectiveness_parallel_plate(2.0, frequency=1.0, mdot=0.08,
                                                  plate_thickness=0.00025,
                                                  plate_spacing=0.0001)
    assert 0.0 <= r["eps"] <= 0.97 + 1e-9


def test_parallel_plate_porosity_matches_spacing_over_period():
    r = regenerator_effectiveness_parallel_plate(2.0, frequency=1.0, mdot=0.08,
                                                  plate_thickness=0.0005,
                                                  plate_spacing=0.0002)
    assert r["porosity"] == pytest.approx(0.0002 / (0.0002 + 0.0005))


def test_parallel_plate_pumping_power_increases_as_spacing_shrinks():
    p_wide = pumping_power_parallel_plate(0.08, plate_spacing=0.0005, mass_regenerator=2.0)
    p_narrow = pumping_power_parallel_plate(0.08, plate_spacing=0.00005, mass_regenerator=2.0)
    assert p_narrow["P_pump_W"] > p_wide["P_pump_W"]


def test_parallel_plate_hydraulic_diameter_is_twice_spacing():
    r = pumping_power_parallel_plate(0.08, plate_spacing=0.0002, mass_regenerator=2.0)
    assert r["d_h_m"] == pytest.approx(2 * 0.0002)