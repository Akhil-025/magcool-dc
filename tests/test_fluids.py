"""
Tests for core/fluids.py -- the heat-transfer-fluid property library.

This module was, until now, exercised only indirectly (every AMRSystem
call site defaults fluid="water" and threads it through fluid_properties()),
with no direct test coverage of its own -- despite core/fluid_selection_
optimization.py (a wired-in, every-pipeline-pass production step, see
main.py step 18) depending on every entry in FLUID_LIBRARY being
internally consistent and fluid_properties() raising on bad input rather
than silently degrading to water.
"""
import pytest

from core.fluids import FLUID_LIBRARY, FLUID_NAMES, DEFAULT_FLUID, fluid_properties


REQUIRED_KEYS = {"rho", "cp", "mu", "k", "freeze_C", "description"}


def test_fluid_library_entries_have_all_required_fields():
    for name, entry in FLUID_LIBRARY.items():
        missing = REQUIRED_KEYS - set(entry)
        assert not missing, f"{name} missing fields: {missing}"


def test_fluid_library_property_values_are_physically_positive():
    """rho, cp, mu, k must all be positive for every fluid -- a zero or
    negative value here would silently corrupt every downstream Qc/COP
    calculation that divides by or multiplies through these properties."""
    for name, entry in FLUID_LIBRARY.items():
        assert entry["rho"] > 0.0, name
        assert entry["cp"] > 0.0, name
        assert entry["mu"] > 0.0, name
        assert entry["k"] > 0.0, name


def test_fluid_names_matches_library_keys():
    assert set(FLUID_NAMES) == set(FLUID_LIBRARY)
    assert len(FLUID_NAMES) == len(FLUID_LIBRARY)


def test_default_fluid_is_a_valid_library_entry():
    """DEFAULT_FLUID is referenced by name (not value) elsewhere in the
    repo (core/fluid_selection_optimization.py's summary line, main.py's
    step-18 synthesis) -- if it ever drifted to a typo'd/removed key,
    those call sites would raise deep in a pipeline run rather than at
    import time."""
    assert DEFAULT_FLUID in FLUID_LIBRARY


def test_water_is_the_original_idealized_baseline():
    """water must stay in the library (every existing AMRSystem call site
    pins fluid="water" as its own function-level default -- see fluids.py's
    module docstring) with the same property values every calibrated
    Qc/COP number in the repo was fit against."""
    water = FLUID_LIBRARY["water"]
    assert water["rho"] == pytest.approx(997.0)
    assert water["cp"] == pytest.approx(4186.0)
    assert water["mu"] == pytest.approx(8.9e-4)
    assert water["k"] == pytest.approx(0.606)


def test_fluid_properties_returns_expected_keys_only():
    props = fluid_properties("water")
    assert set(props) == {"rho", "cp", "mu", "k"}


def test_fluid_properties_matches_library_for_every_fluid():
    for name in FLUID_LIBRARY:
        props = fluid_properties(name)
        entry = FLUID_LIBRARY[name]
        assert props["rho"] == entry["rho"]
        assert props["cp"] == entry["cp"]
        assert props["mu"] == entry["mu"]
        assert props["k"] == entry["k"]


def test_fluid_properties_default_argument_is_water():
    assert fluid_properties() == fluid_properties("water")


def test_fluid_properties_accepts_but_does_not_require_T_K():
    """T_K is accepted for interface compatibility with the original
    water_properties(T_K) (see module docstring) but does not currently
    vary the returned values -- both call styles must agree."""
    assert fluid_properties("water", T_K=300.0) == fluid_properties("water")
    assert fluid_properties("water", T_K=250.0) == fluid_properties("water")


def test_fluid_properties_unknown_fluid_raises_value_error():
    with pytest.raises(ValueError):
        fluid_properties("liquid_nitrogen")


def test_fluid_properties_unknown_fluid_error_lists_valid_options():
    """The error message is meant to be actionable at a REPL/pipeline
    traceback -- it should name every real option, not just say
    'invalid'."""
    with pytest.raises(ValueError) as excinfo:
        fluid_properties("not_a_real_fluid")
    message = str(excinfo.value)
    for name in FLUID_LIBRARY:
        assert name in message


def test_fluid_properties_never_silently_falls_back_to_water():
    """Regression guard for the exact failure mode the module docstring
    calls out: an unrecognized fluid must raise, not silently return
    water's properties."""
    with pytest.raises(ValueError):
        fluid_properties("wter")  # typo of "water"


def test_ethanol_has_lowest_freeze_point_and_weakest_thermal_properties():
    """Cross-checks the module docstring's own characterization of
    ethanol ('poor cp/k make it a weak choice at room temperature',
    'only relevant for sub-zero-span AMR variants')."""
    ethanol = FLUID_LIBRARY["ethanol"]
    others = [e for name, e in FLUID_LIBRARY.items() if name != "ethanol"]
    assert all(ethanol["freeze_C"] < o["freeze_C"] for o in others)
    assert all(ethanol["cp"] < o["cp"] for o in others)
    assert all(ethanol["k"] < o["k"] for o in others)


def test_glycol_fraction_lowers_freeze_point_relative_to_pure_water():
    """Physical sanity check on the water/ethylene-glycol family: more
    glycol -> lower freeze point (the entire reason these mixtures are
    used over plain water in real hardware)."""
    water = FLUID_LIBRARY["water"]["freeze_C"]
    eg10 = FLUID_LIBRARY["water_eg10"]["freeze_C"]
    eg20 = FLUID_LIBRARY["water_eg20"]["freeze_C"]
    assert water > eg10 > eg20
