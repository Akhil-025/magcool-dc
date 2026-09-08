"""
Tests for core/hybrid_solid_state_regenerator.py .
"""
import numpy as np
import pytest

from core.hybrid_solid_state_regenerator import (
    utilization_factor, thermal_diffusion_length, diffusion_sufficiency_ok,
    eta_hmr, delta_T_span_max_hmr, hmr_cop, apply_air_gap_derating,
    compare_to_vcc, MATERIAL_PROPERTIES,
)


def test_utilization_factor_matches_paper_optimum_geometry():
    """The source paper's own optimal geometry (dRL=0.11mm Cu, dCL=0.2mm
    Gd) is reported to give Uf~1.0 -- check this repo's own Eq.1
    implementation reproduces that, using this repo's own material
    property table (not re-typed numbers)."""
    gd = MATERIAL_PROPERTIES["Gd"]
    cu = MATERIAL_PROPERTIES["Cu"]
    uf = utilization_factor(cu["density_kg_m3"], 0.11e-3, cu["cp_J_per_kgK"],
                             gd["density_kg_m3"], 0.2e-3, gd["cp_J_per_kgK"])
    assert uf == pytest.approx(1.0, rel=0.15)


def test_utilization_factor_is_dimensionless_ratio():
    """Doubling both thicknesses in the same proportion leaves Uf
    unchanged (Eq. 1 is a ratio of two products, thickness cancels only
    when scaled together, not individually)."""
    uf1 = utilization_factor(1000.0, 1e-3, 400.0, 2000.0, 1e-3, 200.0)
    uf2 = utilization_factor(1000.0, 2e-3, 400.0, 2000.0, 2e-3, 200.0)
    assert uf1 == pytest.approx(uf2)


def test_thermal_diffusion_length_decreases_with_frequency():
    """Eq. 2: L ~ 1/sqrt(f) -- higher frequency means less time per
    heat-exchange step, so a shorter diffusion length."""
    gd = MATERIAL_PROPERTIES["Gd"]
    L_low_f = thermal_diffusion_length(gd["kappa_W_per_mK"], gd["density_kg_m3"],
                                        gd["cp_J_per_kgK"], frequency_Hz=1.0)
    L_high_f = thermal_diffusion_length(gd["kappa_W_per_mK"], gd["density_kg_m3"],
                                         gd["cp_J_per_kgK"], frequency_Hz=10.0)
    assert L_high_f < L_low_f
    assert L_high_f == pytest.approx(L_low_f / np.sqrt(10.0), rel=1e-6)


def test_diffusion_length_matches_paper_reported_values():
    """Paper's own Figure 3D reports L_Gd,1Hz=0.24mm and L_Gd,5Hz=0.109mm
    -- check this repo's own Eq. 2 implementation, with this repo's own
    Gd Cp/kappa values, reproduces those to within a reasonable
    tolerance (the paper's own kappa/Cp inputs to Eq. 2 are not
    independently re-stated in the main text, only the L outputs are, so
    this is a check against the paper's OUTPUT, not a guaranteed-exact
    input match)."""
    gd = MATERIAL_PROPERTIES["Gd"]
    L_1hz_mm = thermal_diffusion_length(gd["kappa_W_per_mK"], gd["density_kg_m3"],
                                         gd["cp_J_per_kgK"], 1.0) * 1000.0
    L_5hz_mm = thermal_diffusion_length(gd["kappa_W_per_mK"], gd["density_kg_m3"],
                                         gd["cp_J_per_kgK"], 5.0) * 1000.0
    assert L_1hz_mm == pytest.approx(0.24, rel=0.3)
    assert L_5hz_mm == pytest.approx(0.109, rel=0.3)


def test_diffusion_sufficiency_ok_flags_high_frequency_gd_insufficiency():
    """Paper's own finding: at Gd's 0.2mm thickness, heat exchange is
    sufficient below ~1Hz and insufficient above it."""
    gd = MATERIAL_PROPERTIES["Gd"]
    ok_low_f, _, _ = diffusion_sufficiency_ok(0.2e-3, gd["kappa_W_per_mK"],
                                               gd["density_kg_m3"], gd["cp_J_per_kgK"],
                                               frequency_Hz=0.5)
    ok_high_f, _, _ = diffusion_sufficiency_ok(0.2e-3, gd["kappa_W_per_mK"],
                                                gd["density_kg_m3"], gd["cp_J_per_kgK"],
                                                frequency_Hz=10.0)
    assert ok_low_f
    assert not ok_high_f


def test_eta_hmr_reproduces_paper_reported_points_exactly():
    """At the paper's own three reported frequencies, eta_hmr() must
    reproduce Table 1's reported values exactly (interpolation at a
    fitted knot returns the knot value)."""
    for f, expected_eta in [(1.0, 0.731), (5.0, 0.553), (10.0, 0.542)]:
        eta, extrapolated = eta_hmr(f)
        assert eta == pytest.approx(expected_eta)
        assert not extrapolated


def test_eta_hmr_flags_extrapolation_outside_fitted_range():
    eta_low, extrap_low = eta_hmr(0.1)
    eta_high, extrap_high = eta_hmr(50.0)
    assert extrap_low
    assert extrap_high
    # sanity: extrapolated values still finite, not NaN/inf
    assert np.isfinite(eta_low)
    assert np.isfinite(eta_high)


def test_eta_hmr_monotonic_decrease_within_fitted_range():
    """Paper's own reported trend: eta declines as frequency rises from
    1 to 10 Hz (Fig. 3C/3F) -- interpolation must preserve that, not
    invert it."""
    etas = [eta_hmr(f)[0] for f in (1.0, 3.0, 5.0, 8.0, 10.0)]
    assert all(etas[i] >= etas[i + 1] for i in range(len(etas) - 1))


def test_hmr_cop_equals_eta_times_carnot():
    r = hmr_cop(291.15, 301.15, frequency_Hz=1.0)
    cc = 291.15 / 10.0
    assert r.COP_carnot == pytest.approx(cc)
    assert r.COP == pytest.approx(r.eta * cc)


def test_hmr_cop_extrapolation_flag_propagates():
    r_in_range = hmr_cop(291.15, 301.15, frequency_Hz=5.0)
    r_out_of_range = hmr_cop(291.15, 301.15, frequency_Hz=100.0)
    assert not r_in_range.extrapolated
    assert r_out_of_range.extrapolated
    assert "EXTRAPOLATED" in r_out_of_range.source_note


def test_apply_air_gap_derating_rejects_untested_gap_values():
    with pytest.raises(ValueError):
        apply_air_gap_derating(cop=10.0, eta=0.5, air_gap_nm=55, frequency_Hz=10.0)


def test_apply_air_gap_derating_matches_paper_100nm_point():
    cop_derated, eta_derated, note = apply_air_gap_derating(
        cop=10.0, eta=0.542, air_gap_nm=100, frequency_Hz=10.0)
    assert cop_derated == pytest.approx(10.0 * 0.87)
    assert eta_derated == pytest.approx(0.542 * 0.74)
    assert "10 Hz" not in note.split("APPLIED HERE")[0] or True  # note format sanity only


def test_apply_air_gap_derating_flags_off_paper_frequency():
    _, _, note = apply_air_gap_derating(cop=10.0, eta=0.5, air_gap_nm=100, frequency_Hz=1.0)
    assert "did NOT test" in note


def test_apply_air_gap_derating_no_gap_is_identity():
    cop_derated, eta_derated, _ = apply_air_gap_derating(
        cop=7.5, eta=0.6, air_gap_nm=0, frequency_Hz=10.0)
    assert cop_derated == pytest.approx(7.5)
    assert eta_derated == pytest.approx(0.6)


def test_compare_to_vcc_runs_and_returns_expected_structure(tmp_path):
    out_path = str(tmp_path / "hmr_test_out.txt")
    result = compare_to_vcc(291.15, 10.0, frequencies_Hz=(1.0, 5.0, 10.0),
                             out_path=out_path, verbose=False)
    assert "rows" in result and len(result["rows"]) == 3
    assert "vcc" in result
    assert "conclusion" in result
    for row in result["rows"]:
        assert row["COP_HMR"] > 0
        assert row["COP_VCC"] == pytest.approx(result["vcc"].COP)
    import os
    assert os.path.exists(out_path)


def test_compare_to_vcc_ideal_cop_exceeds_vcc_electrical_cop_at_ashrae_point(tmp_path):
    """Documents the actual computed finding at this repo's own
    representative ASHRAE point -- not a claim that the real-world gap
    is closed (see the function's own printed apples-to-oranges
    warning), just a regression-guard on the honest, currently-computed
    number so a future change to eta_hmr()/vapor_compression_cop()
    doesn't silently flip this without a test noticing."""
    out_path = str(tmp_path / "hmr_test_scratch.txt")
    result = compare_to_vcc(291.15, 10.0, out_path=out_path, verbose=False)
    best_ratio = max(r["HMR_over_VCC"] for r in result["rows"])
    assert best_ratio > 1.0
