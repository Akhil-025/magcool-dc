"""Phase 31 addition: core/regenerator_1d.py had no dedicated test file
even though it is imported by main.py's pipeline (see README's Tier-1
test-coverage gap). simulate_amr_1d() itself is a genuine multi-cycle
transient simulation (tens of seconds per call at production settings --
see the module's own docstring) so every test below uses drastically
reduced n_nodes/max_cycles purely to exercise the real code path fast
(confirmed ~0.03s/call at these settings); none of these check for
convergence or benchmark accuracy -- that is validate_against_benchmarks()'s
job, run separately (see results/regenerator_1d_validation.txt)."""

import numpy as np
import pytest

from core.mce_material import GADOLINIUM
from core.thermal import K_SOLID_GD
from core.regenerator_1d import (
    simulate_amr_1d,
    _apply_axial_conduction,
    _packed_bed_effective_axial_conductivity,
    _cache_key,
    _MODEL_VERSION,
)


# --- _packed_bed_effective_axial_conductivity (Phase 31's new correlation) ---

def test_packed_bed_conductivity_approaches_pure_fluid_as_porosity_to_one():
    k = _packed_bed_effective_axial_conductivity(0.999999, k_fluid=0.61, k_solid=10.5)
    assert k == pytest.approx(0.61, rel=1e-4)


def test_packed_bed_conductivity_approaches_pure_solid_as_porosity_to_zero():
    k = _packed_bed_effective_axial_conductivity(1e-9, k_fluid=0.61, k_solid=10.5)
    assert k == pytest.approx(10.5, rel=1e-4)


def test_packed_bed_conductivity_at_typical_gd_bed_porosity_is_between_bounds():
    # porosity=0.365 is this module's own default bed_porosity for Gd runs.
    k = _packed_bed_effective_axial_conductivity(0.365, k_fluid=0.61, k_solid=10.5)
    series_bound = 1.0 / ((1 - 0.365) / 10.5 + 0.365 / 0.61)  # harmonic mean (lower bound)
    parallel_bound = 0.365 * 0.61 + (1 - 0.365) * 10.5         # arithmetic mean (upper bound)
    assert series_bound < k < parallel_bound


def test_packed_bed_conductivity_monotonic_in_porosity():
    # Less porosity (more solid) should give higher effective conductivity,
    # since k_solid > k_fluid for Gd-in-water.
    k_high_porosity = _packed_bed_effective_axial_conductivity(0.6, 0.61, 10.5)
    k_low_porosity = _packed_bed_effective_axial_conductivity(0.2, 0.61, 10.5)
    assert k_low_porosity > k_high_porosity


def test_packed_bed_conductivity_degenerate_inputs_fall_back_to_fluid():
    assert _packed_bed_effective_axial_conductivity(0.0, 0.61, 10.5) == pytest.approx(0.61)
    assert _packed_bed_effective_axial_conductivity(1.0, 0.61, 10.5) == pytest.approx(0.61)
    assert _packed_bed_effective_axial_conductivity(0.365, 0.0, 10.5) > 0
    assert _packed_bed_effective_axial_conductivity(0.365, 0.61, 0.0) > 0


# --- _apply_axial_conduction ---

def test_apply_axial_conduction_conserves_total_energy():
    T = np.array([310.0, 300.0, 290.0, 280.0])
    m_node = 0.05
    cp = 450.0
    T_before_sum = T.sum()
    T_after = _apply_axial_conduction(T.copy(), dt_total=1.0, dx=0.01, bed_area=0.002,
                                       k_eff_axial=2.0, m_node=m_node, cp_solid_eff=cp)
    # Pure node-to-node conduction with no source/sink must conserve the
    # sum of temperatures across equal-mass nodes (energy conservation).
    assert T_after.sum() == pytest.approx(T_before_sum, abs=1e-6)


def test_apply_axial_conduction_smooths_a_gradient_towards_uniform():
    T = np.array([320.0, 300.0, 280.0])
    T_after = _apply_axial_conduction(T.copy(), dt_total=2.0, dx=0.01, bed_area=0.002,
                                       k_eff_axial=5.0, m_node=0.05, cp_solid_eff=450.0)
    spread_before = T.max() - T.min()
    spread_after = T_after.max() - T_after.min()
    assert spread_after < spread_before


def test_apply_axial_conduction_no_op_for_single_node_or_zero_conductivity():
    T = np.array([305.0])
    T_after = _apply_axial_conduction(T.copy(), dt_total=1.0, dx=0.01, bed_area=0.002,
                                       k_eff_axial=5.0, m_node=0.05, cp_solid_eff=450.0)
    assert np.array_equal(T_after, T)

    T2 = np.array([320.0, 280.0])
    T2_after = _apply_axial_conduction(T2.copy(), dt_total=1.0, dx=0.01, bed_area=0.002,
                                        k_eff_axial=0.0, m_node=0.05, cp_solid_eff=450.0)
    assert np.array_equal(T2_after, T2)


# --- simulate_amr_1d (reduced-scale smoke tests only, see module docstring above) ---

def test_simulate_amr_1d_returns_expected_keys_at_reduced_scale():
    r = simulate_amr_1d(GADOLINIUM, mu0H_max=1.5, mass_total=1.0, frequency=1.0,
                         mdot=0.002, n_nodes=4, max_cycles=8, tol=1.0)
    for key in ("converged", "n_cycles", "NTU_total", "span_K", "span_history_last10"):
        assert key in r


def test_simulate_amr_1d_span_is_non_negative_and_finite():
    r = simulate_amr_1d(GADOLINIUM, mu0H_max=1.5, mass_total=1.0, frequency=1.0,
                         mdot=0.002, n_nodes=4, max_cycles=8, tol=1.0)
    assert np.isfinite(r["span_K"])
    assert r["span_K"] >= 0


def test_simulate_amr_1d_zero_field_gives_near_zero_span():
    # No applied field -> no magnetocaloric effect -> no driving temperature
    # difference for the regenerator to amplify. Not exactly 0.0 in
    # practice at this test's deliberately tiny max_cycles=8 (a numerical
    # floor from the initial-condition transient not having fully settled,
    # not a real physics effect) -- tolerance is set from the observed
    # magnitude (~0.009K) at these reduced settings, not tightened to 1e-6.
    r = simulate_amr_1d(GADOLINIUM, mu0H_max=0.0, mass_total=1.0, frequency=1.0,
                         mdot=0.002, n_nodes=4, max_cycles=8, tol=1.0)
    assert r["span_K"] == pytest.approx(0.0, abs=0.05)


def test_simulate_amr_1d_k_solid_override_changes_result():
    # Phase 31: k_solid is a new parameter (mirrors the existing cp_solid
    # override pattern) -- check the override actually reaches the axial
    # conductivity calculation, not just that it's accepted and ignored.
    r_default = simulate_amr_1d(GADOLINIUM, mu0H_max=1.5, mass_total=1.0, frequency=1.0,
                                 mdot=0.002, n_nodes=4, max_cycles=8, tol=1.0)
    r_override = simulate_amr_1d(GADOLINIUM, mu0H_max=1.5, mass_total=1.0, frequency=1.0,
                                  mdot=0.002, n_nodes=4, max_cycles=8, tol=1.0,
                                  k_solid=K_SOLID_GD * 100)  # implausibly high, for contrast
    assert r_default["span_K"] != r_override["span_K"]


def test_simulate_amr_1d_higher_field_gives_higher_or_equal_span():
    r_low = simulate_amr_1d(GADOLINIUM, mu0H_max=0.5, mass_total=1.0, frequency=1.0,
                             mdot=0.002, n_nodes=4, max_cycles=8, tol=1.0)
    r_high = simulate_amr_1d(GADOLINIUM, mu0H_max=2.0, mass_total=1.0, frequency=1.0,
                              mdot=0.002, n_nodes=4, max_cycles=8, tol=1.0)
    assert r_high["span_K"] >= r_low["span_K"]


# --- cache key / model versioning ---

def test_cache_key_changes_when_model_version_changes(monkeypatch):
    # Phase 31: guards against the exact silent-stale-cache risk this
    # change introduced -- if _MODEL_VERSION is ever bumped again without
    # being folded into the cache key, this test will fail loudly instead
    # of old-physics results being silently served from disk forever.
    # _cache_key() returns an opaque SHA-256 hash (not the raw payload), so
    # the only way to check model_version is actually IN the hashed
    # payload is to vary it and confirm the resulting key changes.
    import core.regenerator_1d as regen1d
    key_args = dict(mu0H_max=1.5, mass_total=1.0, frequency=1.0, n_nodes=20,
                     mdot_search=(0.002,), extra_kwargs={"mdot": 0.002})
    monkeypatch.setattr(regen1d, "_MODEL_VERSION", 999)
    key_a = _cache_key(GADOLINIUM, **key_args)
    monkeypatch.setattr(regen1d, "_MODEL_VERSION", 1000)
    key_b = _cache_key(GADOLINIUM, **key_args)
    assert key_a != key_b


def test_cache_key_differs_for_different_field():
    key_a = _cache_key(GADOLINIUM, 1.0, 1.0, 1.0, n_nodes=20, mdot_search=(0.002,), extra_kwargs={})
    key_b = _cache_key(GADOLINIUM, 2.0, 1.0, 1.0, n_nodes=20, mdot_search=(0.002,), extra_kwargs={})
    assert key_a != key_b
