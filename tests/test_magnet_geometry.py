"""
Unit tests for core/magnet_geometry.py .

Covers: the closed-form Halbach-cylinder relation's basic algebraic
properties (round-trip consistency, monotonicity, the genuine
super-linear-in-field mass growth ROADMAP.md's plan asked for),
input validation, and that `bjork_qualitative_check()`/
`run_magnet_geometry_analysis()` run end-to-end and report their own
honest (non-forced) finding.
"""
import numpy as np
import pytest

from core.magnet_geometry import (
    halbach_bore_field_T,
    halbach_outer_radius_for_field_m,
    halbach_magnet_mass_kg,
    bore_geometry_from_air_gap_volume,
    halbach_field_vs_mass,
    bjork_qualitative_check,
    run_magnet_geometry_analysis,
    DEFAULT_REMANENCE_T,
)


def test_bore_field_and_outer_radius_are_inverses():
    """halbach_outer_radius_for_field_m() inverts halbach_bore_field_T():
    solving for Ro at a target field and then recomputing the field from
    that Ro must reproduce the target field."""
    Ri = 0.02
    for B_target in (0.5, 1.0, 2.0, 3.0):
        Ro = halbach_outer_radius_for_field_m(B_target, Ri)
        B_check = halbach_bore_field_T(Ro, Ri)
        assert B_check == pytest.approx(B_target, rel=1e-9)


def test_bore_field_requires_outer_greater_than_inner():
    with pytest.raises(ValueError):
        halbach_bore_field_T(0.01, 0.02)


def test_outer_radius_monotonically_increases_with_field():
    """Higher target field at fixed bore radius must require a larger
    outer radius (and hence, at fixed length, more magnet mass)."""
    Ri = 0.02
    radii = [halbach_outer_radius_for_field_m(B, Ri) for B in (0.5, 1.0, 1.5, 2.0, 2.5, 3.0)]
    assert radii == sorted(radii)
    assert len(set(radii)) == len(radii)  # strictly increasing


def test_magnet_mass_is_super_linear_in_field():
    """The genuine nonlinearity ROADMAP.md's plan asked for:
    doubling the field should more than double the required magnet mass
    at a fixed air-gap geometry (mass grows like exp(2B/Br), not
    linearly)."""
    air_gap_volume = 0.001
    m1 = halbach_field_vs_mass(1.0, air_gap_volume)["magnet_mass_kg"]
    m2 = halbach_field_vs_mass(2.0, air_gap_volume)["magnet_mass_kg"]
    m3 = halbach_field_vs_mass(3.0, air_gap_volume)["magnet_mass_kg"]
    assert m2 > 2 * m1
    assert m3 > 2 * m2  # accelerating growth, not just non-linear once


def test_magnet_mass_scales_with_bed_length_at_fixed_field():
    """At a fixed target field (and hence fixed Ro/Ri ratio for a fixed
    bore radius), doubling the air-gap volume (i.e. doubling length at
    fixed cross-section) should double the magnet mass -- mass is linear
    in LENGTH even though it is super-linear in FIELD."""
    m1 = halbach_field_vs_mass(1.5, 0.001)["magnet_mass_kg"]
    m2 = halbach_field_vs_mass(1.5, 0.002)["magnet_mass_kg"]
    assert m2 == pytest.approx(2 * m1, rel=1e-9)


def test_bore_geometry_from_air_gap_volume_consistency():
    """Ri = sqrt(A/pi), L = V/A -- and Ri^2*pi*L must reproduce V."""
    A = 0.002
    V = 0.0015
    Ri, L = bore_geometry_from_air_gap_volume(V, A)
    assert Ri == pytest.approx(np.sqrt(A / np.pi))
    assert np.pi * Ri ** 2 * L == pytest.approx(A * L)
    assert A * L == pytest.approx(V, rel=1e-9)


def test_halbach_field_vs_mass_return_shape():
    result = halbach_field_vs_mass(2.0, 0.001)
    for key in ("inner_radius_m", "outer_radius_m", "length_m",
                "outer_to_inner_ratio", "magnet_mass_kg", "mu0H_target_T"):
        assert key in result
    assert result["outer_radius_m"] > result["inner_radius_m"] > 0
    assert result["magnet_mass_kg"] > 0
    assert result["mu0H_target_T"] == 2.0


@pytest.mark.parametrize("fn,args", [
    (halbach_bore_field_T, (0.01, 0.02)),           # Ro < Ri
    (halbach_outer_radius_for_field_m, (-1.0, 0.02)),  # negative field
    (halbach_magnet_mass_kg, (0.02, 0.01, 0.1)),     # Ro < Ri
    (bore_geometry_from_air_gap_volume, (-0.001, 0.002)),  # negative volume
])
def test_invalid_inputs_raise(fn, args):
    with pytest.raises(ValueError):
        fn(*args)


def test_bjork_qualitative_check_runs_and_reports_honestly():
    """bjork_qualitative_check() must run end-to-end, return a best field
    within the swept range, and its `matches_2T_claim` flag must be
    consistent with whether best_field_T is actually within 0.5 T of
    2.0 -- i.e. the function reports what it found rather than always
    claiming a match."""
    result = bjork_qualitative_check()
    assert "rows" in result and len(result["rows"]) > 0
    assert "best_field_T" in result
    fields = [r["mu0H_T"] for r in result["rows"]]
    assert result["best_field_T"] in fields
    expected_match = abs(result["best_field_T"] - 2.0) <= 0.5
    assert result["matches_2T_claim"] == expected_match
    # every row must have a finite, non-negative cost-per-K value
    for row in result["rows"]:
        assert row["cost_per_K_$"] >= 0


def test_bjork_qualitative_check_cost_per_K_rows_are_all_positive_cost():
    result = bjork_qualitative_check()
    for row in result["rows"]:
        assert row["total_cost_$"] > 0
        assert row["magnet_cost_$"] >= 0
        assert row["mcm_cost_$"] > 0


def test_run_magnet_geometry_analysis_writes_file(tmp_path):
    out_path = tmp_path / "magnet_geometry_analysis.txt"
    text = run_magnet_geometry_analysis(out_path=str(out_path), verbose=False)
    assert out_path.exists()
    written = out_path.read_text()
    assert written == text
    assert "PHASE 19" in text
    assert "Halbach" in text


def test_default_remanence_is_positive_and_reasonable():
    """Sanity bound on the module-level default so a future accidental
    edit (e.g. an order-of-magnitude typo) fails loudly."""
    assert 0.5 < DEFAULT_REMANENCE_T < 2.5