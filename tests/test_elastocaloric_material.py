import math

from core.elastocaloric_material import NiTi_binary, validate_against_literature


def test_delta_T_ad_increases_monotonically_with_stress():
    prev = -1.0
    for sigma in (50, 150, 250, 350, 450, 550):
        dT = NiTi_binary.delta_T_ad(sigma)
        assert dT > prev
        prev = dT


def test_transformed_fraction_clamped_to_unit_interval():
    assert NiTi_binary.transformed_fraction(-10.0) == 0.0
    assert NiTi_binary.transformed_fraction(1e6) == 1.0
    assert 0.0 <= NiTi_binary.transformed_fraction(300.0) <= 1.0


def test_hysteresis_loss_nonnegative_and_saturates():
    loss_low = NiTi_binary.hysteresis_loss_J_per_kg(100.0)
    loss_high = NiTi_binary.hysteresis_loss_J_per_kg(NiTi_binary.sigma_sat_MPa)
    loss_beyond = NiTi_binary.hysteresis_loss_J_per_kg(NiTi_binary.sigma_sat_MPa * 2)
    assert loss_low >= 0.0
    assert loss_high >= loss_low
    assert math.isclose(loss_high, loss_beyond)


def test_validate_against_literature_reports_overprediction_honestly():
    """Locks in the module's own documented finding: this model
    over-predicts Delta_T_ad relative to the commonly-cited 10-20K device
    range at realistic device stresses. If a future parameter change makes
    this pass 'in range', that's a real change worth noticing -- this test
    exists so it doesn't happen silently."""
    rows = validate_against_literature(verbose=False)
    in_range_flags = [in_range for _sigma, _dT, in_range in rows]
    assert not all(in_range_flags), (
        "model no longer over-predicts at every checked stress -- "
        "re-check module docstring's honesty flag, it may need updating"
    )
