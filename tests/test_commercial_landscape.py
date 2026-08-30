"""Phase 31 addition: core/commercial_landscape.py had no dedicated test
file even though it is wired into main.py's pipeline (see README's Tier-1
test-coverage gap)."""

from core.commercial_landscape import (
    COMMERCIAL_SYSTEMS,
    disambiguation_note,
    model_prediction_at_cooltech_point,
)


def test_commercial_systems_has_both_magnetocaloric_and_naming_collision_entries():
    magnetocaloric = [c for c in COMMERCIAL_SYSTEMS if c.is_magnetocaloric]
    collisions = [c for c in COMMERCIAL_SYSTEMS if not c.is_magnetocaloric]
    assert len(magnetocaloric) >= 2   # Magnotherm Stellar + Cooltech, at minimum
    assert len(collisions) >= 2       # documents at least the 2 known naming collisions


def test_naming_collision_entries_are_flagged_not_magnetocaloric_in_technology_class():
    for c in COMMERCIAL_SYSTEMS:
        if not c.is_magnetocaloric:
            assert "not-magnetocaloric" in c.technology_class


def test_disambiguation_note_mentions_both_collision_technologies():
    note = disambiguation_note()
    assert "magnetic-bearing" in note.lower() or "york" in note.lower() or "johnson controls" in note.lower()
    assert "desiccant" in note.lower() or "circlemiser" in note.lower()
    assert "active magnetic regenerator" in note.lower() or "amr" in note.lower()


def test_cooltech_entry_has_a_claimed_cop_and_range():
    cooltech = next(c for c in COMMERCIAL_SYSTEMS if c.name.startswith("Cooltech"))
    assert cooltech.claimed_COP is not None
    assert cooltech.claimed_COP_range is not None
    assert cooltech.is_magnetocaloric is True


def test_model_prediction_at_cooltech_point_returns_expected_keys():
    result = model_prediction_at_cooltech_point(verbose=False)
    for key in ("model_COP_electrical", "model_Qc_W",
                "cooltech_claimed_COP_range", "gap_pct_vs_claimed_midpoint"):
        assert key in result


def test_model_prediction_at_cooltech_point_cop_is_physically_sane():
    result = model_prediction_at_cooltech_point(verbose=False)
    assert result["model_COP_electrical"] > 0
    assert result["model_Qc_W"] > 0


def test_model_prediction_gap_pct_is_computed_against_claimed_midpoint():
    result = model_prediction_at_cooltech_point(verbose=False)
    cooltech = next(c for c in COMMERCIAL_SYSTEMS if c.name.startswith("Cooltech"))
    expected_gap = 100 * (result["model_COP_electrical"] - cooltech.claimed_COP) / cooltech.claimed_COP
    assert round(result["gap_pct_vs_claimed_midpoint"], 6) == round(expected_gap, 6)


def test_model_prediction_higher_field_gives_higher_or_equal_cooling_capacity():
    # Sanity check that the wiring into AMRSystem is real, not a stub --
    # more field at fixed everything-else should not reduce Qc.
    low_field = model_prediction_at_cooltech_point(mu0H_max=1.0, verbose=False)
    high_field = model_prediction_at_cooltech_point(mu0H_max=2.0, verbose=False)
    assert high_field["model_Qc_W"] >= low_field["model_Qc_W"]
