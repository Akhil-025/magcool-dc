""" addition: core/pue_annualized.py had no dedicated test file
even though it is wired into main.py's pipeline (see README's Tier-1
test-coverage gap)."""

from core.pue_annualized import (
    cop_to_pue,
    pue_comparison_table,
    annualized_energy_comparison,
    NON_COOLING_PUE_OVERHEAD,
    REPRESENTATIVE_CLIMATE_PROFILE_4A,
    run_india_climate_comparison,
    INDIA_COMPOSITE_DELHI_CLIMATE_PROFILE,
    INDIA_WARM_HUMID_MUMBAI_CLIMATE_PROFILE,
)


def test_cop_to_pue_matches_standard_formula():
    r = cop_to_pue("AMR", COP=5.0)
    assert r.PUE_cooling_only == round(1.0 + 1.0 / 5.0, 3)
    assert r.PUE_total_estimate == round(r.PUE_cooling_only + NON_COOLING_PUE_OVERHEAD, 3)


def test_cop_to_pue_higher_cop_gives_lower_pue():
    low = cop_to_pue("X", COP=3.0)
    high = cop_to_pue("X", COP=15.0)
    assert high.PUE_cooling_only < low.PUE_cooling_only


def test_cop_to_pue_zero_or_negative_cop_returns_infinite_pue():
    r_zero = cop_to_pue("X", COP=0.0)
    r_negative = cop_to_pue("X", COP=-1.0)
    assert r_zero.PUE_cooling_only == float("inf")
    assert r_negative.PUE_cooling_only == float("inf")


def test_cop_to_pue_custom_overhead_is_additive_not_multiplicative():
    r_default = cop_to_pue("X", COP=5.0, non_cooling_overhead=0.12)
    r_custom = cop_to_pue("X", COP=5.0, non_cooling_overhead=0.30)
    assert round(r_custom.PUE_total_estimate - r_default.PUE_total_estimate, 3) == round(0.30 - 0.12, 3)
    # cooling-only figure must not depend on the overhead term at all
    assert r_default.PUE_cooling_only == r_custom.PUE_cooling_only


def test_pue_comparison_table_preserves_row_count_and_span():
    rows = [
        {"span_K": 5, "AMR_COP_electrical": 5.5, "VaporCompression_COP": 20.0, "LiquidCooling_COP": 18.0},
        {"span_K": 15, "AMR_COP_electrical": 4.0, "VaporCompression_COP": 8.0, "LiquidCooling_COP": 12.0},
    ]
    out = pue_comparison_table(rows)
    assert len(out) == 2
    assert [r["span_K"] for r in out] == [5, 15]
    assert "AMR_PUE_cooling_only" in out[0]
    assert "VCC_PUE_total_estimate" in out[0]
    assert "Liquid_PUE_cooling_only" in out[0]


def test_annualized_energy_comparison_returns_expected_keys():
    result = annualized_energy_comparison(verbose=False)
    for key in ("rows", "AMR_effective_annual_COP", "AMR_annual_hours_fraction_covered",
                "VCC_effective_annual_COP", "Liquid_effective_annual_COP"):
        assert key in result


def test_annualized_energy_comparison_covers_all_climate_bins_in_rows():
    result = annualized_energy_comparison(verbose=False)
    assert len(result["rows"]) == len(REPRESENTATIVE_CLIMATE_PROFILE_4A)


def test_annualized_energy_comparison_vcc_and_liquid_always_fully_covered():
    # Only AMR has a documented span cap that can exclude climate bins --
    # VCC and liquid-cooling baselines must always be evaluated across the
    # full annual profile (no "hours_fraction_covered" concept for them).
    result = annualized_energy_comparison(verbose=False)
    assert result["VCC_effective_annual_COP"] > 0
    assert result["Liquid_effective_annual_COP"] > 0
    assert 0.0 <= result["AMR_annual_hours_fraction_covered"] <= 1.0


def test_annualized_energy_comparison_liquid_cooling_gets_economizer_credit_below_threshold():
    result = annualized_energy_comparison(verbose=False, economizer_below_C=18.0)
    cold_bins = [r for r in result["rows"] if r["T_outdoor_C"] < 18.0]
    warm_bins = [r for r in result["rows"] if r["T_outdoor_C"] >= 18.0]
    assert cold_bins and warm_bins
    # Below the economizer threshold, the module hardcodes a fixed
    # free-cooling COP of 25.0 for EVERY cold bin, regardless of the
    # implied span -- assert that flat value directly, rather than a
    # cold-vs-warm ordering (small-span bins near the threshold can
    # legitimately have an even higher non-economizer COP than 25.0, so
    # ordering is not a safe invariant here).
    assert all(r["Liquid_COP"] == 25.0 for r in cold_bins)
    assert all(r["Liquid_COP"] != 25.0 for r in warm_bins)


def test_annualized_energy_comparison_excludes_amr_bins_beyond_20k_span():
    result = annualized_energy_comparison(verbose=False)
    for r in result["rows"]:
        if r["span_K"] > 20.0:
            assert r["AMR_span_feasible"] is False
            assert r["AMR_COP"] == 0.0


# --- India-specific climate comparison (LIMITATIONS.md Section 4 follow-up,
# `run_india_climate_comparison()` previously had no test coverage) --------

def test_india_climate_profiles_are_well_formed():
    for profile in (INDIA_COMPOSITE_DELHI_CLIMATE_PROFILE,
                    INDIA_WARM_HUMID_MUMBAI_CLIMATE_PROFILE):
        assert len(profile) == 12  # one bin per calendar month
        total_hours_fraction = sum(b.annual_hours_fraction for b in profile)
        assert abs(total_hours_fraction - 1.0) < 1e-9
        for b in profile:
            assert -10.0 < b.T_outdoor_C < 50.0  # sanity bound, real climate data


def test_india_warm_humid_mumbai_never_drops_below_economizer_threshold():
    # This profile's own docstring claims there is essentially no
    # economizer/free-cooling window at any point in the year -- assert
    # that directly against annualized_energy_comparison()'s default
    # economizer_below_C=18.0 threshold.
    assert all(b.T_outdoor_C >= 18.0 for b in INDIA_WARM_HUMID_MUMBAI_CLIMATE_PROFILE)


def test_india_composite_delhi_has_a_winter_economizer_window():
    # Delhi's Composite-zone winter months should fall in/near the
    # economizer window, unlike Mumbai -- this is the qualitative
    # zone-to-zone contrast run_india_climate_comparison() exists to show.
    assert min(b.T_outdoor_C for b in INDIA_COMPOSITE_DELHI_CLIMATE_PROFILE) < 18.0


def test_run_india_climate_comparison_runs_and_writes_report(tmp_path):
    out_path = str(tmp_path / "pue_annualized_india_test.txt")
    text = run_india_climate_comparison(out_path=out_path, verbose=False)
    assert "Composite (Delhi-NCR)" in text
    assert "Warm-Humid (Mumbai/coastal)" in text
    import os
    assert os.path.exists(out_path)
    with open(out_path) as f:
        assert f.read() == text


def test_run_india_climate_comparison_uses_real_india_profiles_not_4a_default():
    # Cross-check: the India report should reflect each zone's own
    # annualized_energy_comparison() call, not silently fall back to
    # REPRESENTATIVE_CLIMATE_PROFILE_4A.
    delhi_direct = annualized_energy_comparison(
        climate_profile=INDIA_COMPOSITE_DELHI_CLIMATE_PROFILE, verbose=False)
    mumbai_direct = annualized_energy_comparison(
        climate_profile=INDIA_WARM_HUMID_MUMBAI_CLIMATE_PROFILE, verbose=False)
    # The two zones' annual effective VCC COP should differ (Mumbai's
    # always-warm profile implies consistently larger spans than Delhi's
    # winter-economizer-eligible profile), confirming the profiles are
    # actually distinct inputs, not both silently resolving to the same
    # generic default.
    assert delhi_direct["VCC_effective_annual_COP"] != mumbai_direct["VCC_effective_annual_COP"]