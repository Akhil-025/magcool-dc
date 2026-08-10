import numpy as np
import pytest

from core.fluid_mce_cycle import (
    krieger_dougherty_viscosity,
    suspension_effective_properties,
    suspension_delta_T_adiabatic,
    pumping_power_pipe_flow,
    FerrofluidMCESystem,
    DEFAULT_PHI_MAX,
)
from core.mce_material import GADOLINIUM


def test_krieger_dougherty_increases_with_phi():
    mu0 = krieger_dougherty_viscosity(1e-3, 0.0)
    mu1 = krieger_dougherty_viscosity(1e-3, 0.3)
    mu2 = krieger_dougherty_viscosity(1e-3, 0.5)
    assert mu0 == pytest.approx(1e-3)
    assert mu1 < mu2


def test_krieger_dougherty_diverges_near_phi_max():
    mu_far = krieger_dougherty_viscosity(1e-3, 0.3, phi_max=0.63)
    mu_near = krieger_dougherty_viscosity(1e-3, 0.62, phi_max=0.63)
    assert mu_near > 10 * mu_far


def test_krieger_dougherty_rejects_out_of_range_phi():
    with pytest.raises(ValueError):
        krieger_dougherty_viscosity(1e-3, -0.1)
    with pytest.raises(ValueError):
        krieger_dougherty_viscosity(1e-3, 0.63, phi_max=0.63)
    with pytest.raises(ValueError):
        krieger_dougherty_viscosity(1e-3, 0.70, phi_max=0.63)


def test_suspension_properties_bounds():
    props_low = suspension_effective_properties(0.01)
    props_high = suspension_effective_properties(0.4)
    # more particle material -> higher density (magnetite >> water)
    assert props_high["rho_susp_kg_m3"] > props_low["rho_susp_kg_m3"]
    # dilution factor (particle share of total heat capacity) grows with phi
    assert 0.0 <= props_low["dilution_factor"] <= 1.0
    assert props_high["dilution_factor"] > props_low["dilution_factor"]


def test_suspension_dTad_zero_at_zero_phi():
    dTad = suspension_delta_T_adiabatic(GADOLINIUM, 294.5, 1.5, 0.0)
    assert dTad == pytest.approx(0.0, abs=1e-9)


def test_suspension_dTad_monotonically_increases_with_phi():
    phis = [0.01, 0.05, 0.1, 0.2, 0.3]
    dTads = [suspension_delta_T_adiabatic(GADOLINIUM, 294.5, 1.5, p) for p in phis]
    assert all(b > a for a, b in zip(dTads, dTads[1:]))


def test_suspension_dTad_below_pure_material_dTad():
    """Dilution can only ever reduce dTad relative to the pure bulk
    material -- this is the whole point of the mixture heat-capacity
    argument (module docstring physics item 1)."""
    mu0 = 4 * np.pi * 1e-7
    dTad_pure = float(GADOLINIUM.delta_T_adiabatic(np.array([294.5]), 1.5 / mu0)[0])
    dTad_susp = suspension_delta_T_adiabatic(GADOLINIUM, 294.5, 1.5, 0.3)
    assert dTad_susp < dTad_pure


def test_pumping_power_zero_at_zero_flow():
    assert pumping_power_pipe_flow(0.0, 1e-3, 1000.0) == 0.0


def test_pumping_power_increases_with_viscosity():
    p_low = pumping_power_pipe_flow(0.05, 1e-3, 1000.0)
    p_high = pumping_power_pipe_flow(0.05, 5e-3, 1000.0)
    assert p_high > p_low


def test_ferrofluid_system_rejects_invalid_phi():
    with pytest.raises(ValueError):
        FerrofluidMCESystem(GADOLINIUM, 1.5, -0.1, 0.05)
    with pytest.raises(ValueError):
        FerrofluidMCESystem(GADOLINIUM, 1.5, DEFAULT_PHI_MAX, 0.05)


def test_ferrofluid_system_run_returns_sane_result():
    sys_ = FerrofluidMCESystem(GADOLINIUM, mu0H_max=1.5, particle_volume_fraction=0.2,
                                 fluid_mdot=0.05)
    result = sys_.run(291.0, 0.5)
    assert result.Qc >= 0.0
    assert result.W_mag >= 0.0
    assert result.W_parasitic >= 0.0
    assert 0.0 <= result.dilution_factor <= 1.0
    assert result.viscosity_Pa_s > 0.0


def test_ferrofluid_qc_falls_to_zero_beyond_dTad_span():
    """No regeneration (module docstring physics item 3): once the
    imposed span exceeds the suspension's own dTad, Qc must be zero --
    there is no regenerative amplification to fall back on."""
    sys_ = FerrofluidMCESystem(GADOLINIUM, mu0H_max=1.5, particle_volume_fraction=0.2,
                                 fluid_mdot=0.05)
    _, dTad = sys_.cooling_capacity(291.0, 0.1)
    Qc_beyond, _ = sys_.cooling_capacity(291.0, dTad * 2.0)
    assert Qc_beyond == 0.0


def test_ferrofluid_characteristic_curve_length_matches_spans():
    sys_ = FerrofluidMCESystem(GADOLINIUM, mu0H_max=1.5, particle_volume_fraction=0.2,
                                 fluid_mdot=0.05)
    spans = [0.1, 0.3, 0.5]
    curve = sys_.characteristic_curve(291.0, spans)
    assert len(curve) == len(spans)
    assert all(r.T_span == s for r, s in zip(curve, spans))


def test_carrier_water_only():
    with pytest.raises(ValueError):
        suspension_effective_properties(0.1, carrier="oil")
