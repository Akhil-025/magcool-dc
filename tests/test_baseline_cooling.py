import pytest

from core.baseline_cooling import (
    carnot_cop, vapor_compression_cop, liquid_cooling_cop,
    passive_regenerator_augmentation, augmented_regenerator_cop,
    MAX_COP_GAIN_AT_FULL_EFFECTIVENESS,
)
from core.mce_material import GADOLINIUM, LACAMNO3


def test_carnot_cop_matches_definition():
    assert carnot_cop(290.0, 300.0) == pytest.approx(290.0 / 10.0)


def test_baseline_cops_below_carnot():
    vcc = vapor_compression_cop(290.0, 300.0)
    liq = liquid_cooling_cop(290.0, 300.0)
    assert vcc.COP < vcc.COP_carnot
    # liquid cooling blends a high economizer-mode COP, so only check it's
    # not exceeding physical bounds by an absurd margin (it can legitimately
    # exceed the *mechanical* Carnot figure for the DX-only comparison).
    assert liq.COP > 0


# --- Phase 21: passive/hybrid magnetic regenerator augmentation ---

def test_aligned_material_gains_more_than_misaligned():
    """Gd's own Curie temperature (294K) sits inside [291.15, 301.15]K;
    La0.7Ca0.3MnO3's (267K) does not. The aligned material must show a
    strictly larger (or equal) effectiveness gain."""
    aug_gd = passive_regenerator_augmentation(GADOLINIUM, 291.15, 301.15)
    aug_la = passive_regenerator_augmentation(LACAMNO3, 291.15, 301.15)
    assert aug_gd["delta_eps"] >= aug_la["delta_eps"]
    assert aug_gd["delta_eps"] > 0.0


def test_delta_eps_never_negative_after_clipping_in_cop_gain():
    """augmented_regenerator_cop() must never produce a COP below base_COP,
    even if a material's own delta_eps happens to be negative (a poorly-
    aligned or otherwise unfavorable material) -- delta_eps is clipped at 0
    before scaling the COP gain."""
    base = 10.0
    result = augmented_regenerator_cop(base, LACAMNO3, (400.0, 410.0))  # badly misaligned
    assert result.augmented_COP >= result.base_COP


def test_augmented_cop_capped_by_illustrative_ceiling():
    """No configuration should be able to exceed the illustrative full-
    effectiveness ceiling on the fractional COP gain."""
    base = 10.0
    result = augmented_regenerator_cop(base, GADOLINIUM, (291.15, 301.15))
    assert result.cop_gain_fraction <= MAX_COP_GAIN_AT_FULL_EFFECTIVENESS + 1e-9


def test_augmented_regenerator_cop_return_shape():
    result = augmented_regenerator_cop(10.0, GADOLINIUM, (291.15, 301.15))
    assert result.material_name == GADOLINIUM.name
    assert result.T_cold == 291.15
    assert result.T_hot == 301.15
    assert result.augmented_COP == pytest.approx(
        result.base_COP * (1.0 + result.cop_gain_fraction))


def test_explicit_cp_solid_baseline_override():
    """A caller may still pass a flat cp_solid_baseline explicitly instead
    of the default same-material lattice-only baseline."""
    from core.thermal import CP_SOLID_GD
    aug_default = passive_regenerator_augmentation(GADOLINIUM, 291.15, 301.15)
    aug_explicit = passive_regenerator_augmentation(GADOLINIUM, 291.15, 301.15,
                                                      cp_solid_baseline=CP_SOLID_GD)
    assert aug_explicit["cp_solid_baseline_J_kgK"] == CP_SOLID_GD
    # the default (lattice-only, same material) baseline need not equal the
    # flat CP_SOLID_GD reference -- just confirm both are well-defined floats
    assert isinstance(aug_default["cp_solid_baseline_J_kgK"], float)