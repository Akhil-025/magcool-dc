"""
hysteresis_sensitivity.py
==========================
Phase 16 validation deliverable.

Phase 15's merged, globally non-dominated Pareto front
(`core.optimize.run_optimization()`, `results/pareto_front.csv`) came out
100% La(Fe,Si)13Hy. At the time, thermal-hysteresis loss was a documented
but entirely UNQUANTIFIED honesty flag (prose caveats only, in
`core.cascade` and `core.giguere_validation`) -- invisible to every
objective the optimizer actually sees. Phase 16
(`core.first_order_mce.FirstOrderMCEMaterial.hysteresis_loss_J_per_kg`,
`core.amr_cycle.AMRSystem._hysteresis_power_W()`) made it a real number
that increases `W_parasitic` (hence lowers `COP_electrical`) for every
first-order material, while leaving GADOLINIUM (mce_material.py,
second-order/mean-field, genuinely hysteresis-free) untouched.

This module asks the question that motivated Phase 16 in the first place:
does the "100% La(Fe,Si)13Hy" result survive once that loss term is
switched on?

Method
------
Runs `core.optimize.run_optimization()` twice, at IDENTICAL pop_size,
n_gen, and seed -- differing ONLY in whether each first-order material
candidate's `hysteresis_loss_J_per_kg` is left at its Phase-16
literature-placeholder value ("hysteresis ON", the module's new default
behavior after Phase 16) or forced to 0.0 ("hysteresis OFF", i.e. exactly
reproducing pre-Phase-16 behavior) -- by temporarily mutating the three
module-level `*_FIRST_ORDER` constants in `core.first_order_mce` in place
(they are plain, non-frozen dataclass instances, and
`composition_tuned_material()`/`lafesih_composition_tuned_material()`/
`mnfepsi_composition_tuned_material()` all read `hysteresis_loss_J_per_kg`
off these constants at CALL time -- see first_order_mce.py's Phase 16
edits -- so mutating the constants before calling
`core.optimize.run_optimization()` correctly propagates into every tuned
material `_material_candidates()` builds fresh on each call), then
restoring the original values in a `finally` block regardless of outcome.

This is a controlled A/B comparison, not an independent re-derivation:
both runs use the same NSGA-III population size, generation count, and
random seed, so any difference in the merged front's material composition
is attributable to the hysteresis term alone, not to search noise. GD
(GADOLINIUM) is unaffected by either run since it carries no
`hysteresis_loss_J_per_kg` attribute at all (getattr default 0.0 either
way -- see `AMRSystem._hysteresis_power_W()`'s docstring).

Honesty flags (read before trusting this diagnostic's numbers)
----------------------------------------------------------------
1. `pop_size`/`n_gen` default to 32/15 here, NOT `run_optimization()`'s
   own 40/25 production default -- chosen so this diagnostic runs in
   ~15-20s total (4 materials x 2 runs) rather than several minutes, since
   its job is to answer a yes/no "does the front's material composition
   change at all" question, not to produce a publication-quality Pareto
   front. If the ON/OFF comparison below is CLOSE (e.g. hysteresis
   changes only a couple of designs' identity near the current front
   density), rerun with `pop_size=40, n_gen=25` (or higher) before
   treating the result as settled -- NSGA-III's own run-to-run variance
   at this reduced setting has NOT been separately characterized here.
2. The `hysteresis_loss_J_per_kg` VALUES this comparison depends on are
   themselves the least-validated numbers in the whole codebase as of
   Phase 16 -- literature analogs for DIFFERENT exact compositions than
   the ones calibrated here, in one case (MNFEPSI) from a different
   composition axis entirely. See each `*_FIRST_ORDER` constant's own
   block comment in `core/first_order_mce.py` for the full provenance and
   caveats. A qualitative "does hysteresis loss of ROUGHLY THIS
   ORDER OF MAGNITUDE change the answer" reading of this diagnostic's
   output is defensible; a quantitative "the front is now exactly X%
   La(Fe,Si)13Hy" reading is not, until those placeholders are replaced
   with values read directly off the calibrated compositions' own
   hysteresis loops.
3. This diagnostic changes NOTHING about `core.cascade`'s graded-bed
   analyses (Astronautics reproduction, etc.) -- those already run
   through `AMRSystem.run()` and therefore already pick up the Phase 16
   hysteresis term automatically (via `_single_stage()`'s no-loss_model
   path, which Phase 16 deliberately made hysteresis-aware -- see
   `AMRSystem.run()`'s Phase 16 comment). This module only isolates and
   quantifies the effect on the Phase 15 material-selection question
   specifically.
"""
import os
from core.first_order_mce import (
    GD5SI2GE2_FIRST_ORDER, LAFESIH_FIRST_ORDER, MNFEPSI_FIRST_ORDER,
)
import core.optimize as optimize

_HYSTERETIC_CONSTANTS = {
    "Gd5Si2Ge2": GD5SI2GE2_FIRST_ORDER,
    "La(Fe,Si)13Hy": LAFESIH_FIRST_ORDER,
    "Mn-Fe-P-Si": MNFEPSI_FIRST_ORDER,
}


def _material_counts(rows):
    counts = {}
    for r in rows:
        counts[r["material"]] = counts.get(r["material"], 0) + 1
    return counts


def _set_all_hysteresis(value):
    """Mutates the three module-level *_FIRST_ORDER constants'
    hysteresis_loss_J_per_kg in place. Returns nothing; caller is
    responsible for restoring via _restore_all_hysteresis()."""
    for mat in _HYSTERETIC_CONSTANTS.values():
        mat.hysteresis_loss_J_per_kg = value


def run_hysteresis_sensitivity(pop_size=32, n_gen=15, seed=1,
                                out_path="results/hysteresis_sensitivity.txt",
                                out_csv_on="results/pareto_front_hysteresis_on.csv",
                                out_csv_off="results/pareto_front_hysteresis_off.csv",
                                verbose=True):
    """Runs the ON/OFF A/B comparison described in the module docstring
    and writes a human-readable comparison to `out_path`. Returns a dict
    with the raw merged-front rows and material counts from both runs,
    for programmatic use (e.g. by tests).

    `out_csv_on`/`out_csv_off` default to the original Phase 16 fixed
    filenames (fully backward compatible). Pass a seed-specific or
    scratch path to avoid N seeds clobbering the same two files (used by
    `run_hysteresis_multiseed_stability_check()` below). Note:
    `core.optimize.run_optimization()` itself always writes `out_csv`
    (unlike its own `per_material_out_dir`, which DOES accept None) --
    so `out_csv_on`/`out_csv_off` here must be real paths, not None.
    """
    original_values = {label: mat.hysteresis_loss_J_per_kg
                        for label, mat in _HYSTERETIC_CONSTANTS.items()}

    def _p(*args, **kwargs):
        if verbose:
            print(*args, **kwargs)

    try:
        _p("=" * 70)
        _p(f"Phase 16 hysteresis sensitivity: run 1/2 -- hysteresis ON "
           f"(current literature-placeholder values) [seed={seed}]")
        _p("=" * 70)
        _set_all_hysteresis_from(original_values)
        rows_on = optimize.run_optimization(
            pop_size=pop_size, n_gen=n_gen, seed=seed,
            out_csv=out_csv_on,
            per_material_out_dir=None)

        _p()
        _p("=" * 70)
        _p(f"Phase 16 hysteresis sensitivity: run 2/2 -- hysteresis OFF "
           f"(forced to 0.0, reproducing pre-Phase-16 behavior) [seed={seed}]")
        _p("=" * 70)
        _set_all_hysteresis(0.0)
        rows_off = optimize.run_optimization(
            pop_size=pop_size, n_gen=n_gen, seed=seed,
            out_csv=out_csv_off,
            per_material_out_dir=None)
    finally:
        _set_all_hysteresis_from(original_values)

    counts_on = _material_counts(rows_on)
    counts_off = _material_counts(rows_off)
    all_materials = sorted(set(counts_on) | set(counts_off))

    lines = []
    lines.append("Phase 16 hysteresis sensitivity: material composition of the merged,")
    lines.append("globally non-dominated Pareto front, hysteresis ON vs OFF.")
    lines.append(f"(pop_size={pop_size}, n_gen={n_gen}, seed={seed} -- see module")
    lines.append(" docstring honesty flag #1 on why this is smaller than")
    lines.append(" run_optimization()'s own 40/25 production default.)")
    lines.append("")
    lines.append("Hysteresis literature-placeholder values used in the ON run:")
    for label, mat in _HYSTERETIC_CONSTANTS.items():
        lines.append(f"  {label:<16} {original_values[label]:.1f} J/kg  "
                      f"(see core/first_order_mce.py for provenance/caveats)")
    lines.append("")
    lines.append(f"{'Material':<40} {'OFF (pre-Ph16)':>16} {'ON (Ph16)':>12} {'Delta':>8}")
    lines.append("-" * 78)
    for label in all_materials:
        n_off = counts_off.get(label, 0)
        n_on = counts_on.get(label, 0)
        lines.append(f"{label:<40} {n_off:>16} {n_on:>12} {n_on - n_off:>+8}")
    lines.append("-" * 78)
    lines.append(f"{'TOTAL front size':<40} {len(rows_off):>16} {len(rows_on):>12} "
                  f"{len(rows_on) - len(rows_off):>+8}")
    lines.append("")

    lafesih_off = sum(n for label, n in counts_off.items() if "La(Fe,Si)13Hy" in label)
    lafesih_on = sum(n for label, n in counts_on.items() if "La(Fe,Si)13Hy" in label)
    lafesih_frac_off = lafesih_off / len(rows_off) if rows_off else 0.0
    lafesih_frac_on = lafesih_on / len(rows_on) if rows_on else 0.0
    # Matched by substring "La(Fe,Si)13Hy" (GradedFamily.name in
    # core.cascade) rather than an exact label, since the full label
    # embeds the tuned Tc (e.g. "La(Fe,Si)13Hy (tuned, Tc=291.3K)") --
    # see core.optimize._material_candidates()'s f-string. Tc-tuning
    # depends only on T_mid_K/mu0H_max_for_tuning, NOT on
    # hysteresis_loss_J_per_kg, so the embedded Tc (and hence the exact
    # label string) is expected to be IDENTICAL between the ON and OFF
    # runs -- if it isn't, that itself would be worth investigating
    # separately, since it would mean something other than hysteresis
    # differed between the two runs.
    lines.append(f"La(Fe,Si)13Hy share of merged front: "
                  f"{lafesih_frac_off:.0%} (OFF) -> {lafesih_frac_on:.0%} (ON)")
    lines.append("")
    lines.append("See this module's docstring (honesty flags 1-3) before treating the")
    lines.append("above as a settled, publication-quality answer rather than a")
    lines.append("directional sensitivity check.")

    if out_path:
        os.makedirs(os.path.dirname(out_path), exist_ok=True)
        with open(out_path, "w") as f:
            f.write("\n".join(lines) + "\n")

    _p()
    _p("\n".join(lines))
    if out_path:
        _p(f"\nWrote {out_path}")

    return {
        "rows_on": rows_on, "rows_off": rows_off,
        "counts_on": counts_on, "counts_off": counts_off,
        "lafesih_frac_on": lafesih_frac_on, "lafesih_frac_off": lafesih_frac_off,
    }


def _set_all_hysteresis_from(values_by_label):
    for label, mat in _HYSTERETIC_CONSTANTS.items():
        mat.hysteresis_loss_J_per_kg = values_by_label[label]


def run_hysteresis_multiseed_stability_check(
        seeds=(1, 2, 3), pop_size=40, n_gen=25,
        out_path="results/hysteresis_multiseed_stability.txt"):
    """Closes the open item this module's own docstring (honesty flag #1)
    and ROADMAP.md's Phase 16 entry both flagged: does the ON/OFF La(Fe,
    Si)13Hy-share reversal found at the reduced pop_size=32/n_gen=15/
    seed=1 setting (88% OFF -> 100% ON) hold at `run_optimization()`'s
    own production pop_size=40/n_gen=25 default, and is it stable across
    multiple random seeds -- or was it itself just search noise at the
    smaller setting?

    Reruns `run_hysteresis_sensitivity()` once per seed in `seeds` at
    (pop_size, n_gen), writing each seed's ON/OFF fronts to a single pair
    of SCRATCH csv paths that every seed overwrites in turn (`core.
    optimize.run_optimization()` always writes its `out_csv` argument --
    unlike `per_material_out_dir`, it has no None/skip option -- so N
    seeds' worth of intermediate fronts are not worth keeping on disk;
    the seed-by-seed La(Fe,Si)13Hy-share numbers, which ARE the actual
    deliverable of this check, are captured in memory and written to
    `out_path` instead). Each seed's ON/OFF pair is independently
    reproducible via `run_hysteresis_sensitivity(pop_size=pop_size,
    n_gen=n_gen, seed=<seed>)` if the underlying fronts are ever needed
    again.

    Returns a dict with per-seed results and a `stable` boolean --
    True iff every seed's ON run's La(Fe,Si)13Hy share is >= its own OFF
    run's share (i.e. the ON >= OFF DIRECTION of the original finding
    holds at every seed, not necessarily the exact 88%->100% magnitude).
    """
    per_seed = []
    scratch_on = "results/_scratch_hysteresis_multiseed_on.csv"
    scratch_off = "results/_scratch_hysteresis_multiseed_off.csv"
    for seed in seeds:
        result = run_hysteresis_sensitivity(
            pop_size=pop_size, n_gen=n_gen, seed=seed,
            out_path=None, out_csv_on=scratch_on, out_csv_off=scratch_off,
            verbose=False)
        per_seed.append({
            "seed": seed,
            "lafesih_frac_off": result["lafesih_frac_off"],
            "lafesih_frac_on": result["lafesih_frac_on"],
            "front_size_off": len(result["rows_off"]),
            "front_size_on": len(result["rows_on"]),
        })
        print(f"  seed={seed}: La(Fe,Si)13Hy share OFF={result['lafesih_frac_off']:.0%} "
              f"-> ON={result['lafesih_frac_on']:.0%}  "
              f"(front size OFF={len(result['rows_off'])}, ON={len(result['rows_on'])})")

    for scratch in (scratch_on, scratch_off):
        try:
            os.remove(scratch)
        except OSError:
            pass

    stable = all(s["lafesih_frac_on"] >= s["lafesih_frac_off"] - 1e-9 for s in per_seed)

    lines = []
    lines.append("Phase 16 open item (see ROADMAP.md, and this module's own honesty flag #1):")
    lines.append("does the ON/OFF La(Fe,Si)13Hy-share reversal found at pop_size=32/n_gen=15/")
    lines.append("seed=1 (88% OFF -> 100% ON) hold at production pop_size/n_gen settings and")
    lines.append("across multiple seeds, or was it search noise at the smaller setting?")
    lines.append("")
    lines.append(f"Settings: pop_size={pop_size}, n_gen={n_gen}, seeds={list(seeds)}")
    lines.append("")
    lines.append(f"{'seed':>6} {'frac_OFF':>10} {'frac_ON':>10} {'front_OFF':>10} {'front_ON':>10}")
    lines.append("-" * 50)
    for s in per_seed:
        lines.append(f"{s['seed']:>6} {s['lafesih_frac_off']:>9.0%} "
                      f"{s['lafesih_frac_on']:>9.0%} {s['front_size_off']:>10} "
                      f"{s['front_size_on']:>10}")
    lines.append("-" * 50)
    lines.append("")
    if stable:
        lines.append("RESULT: STABLE. At every seed checked, the merged front's La(Fe,Si)13Hy")
        lines.append("share with hysteresis ON is >= its own OFF-run share -- the original")
        lines.append("pop_size=32/n_gen=15/seed=1 finding's DIRECTION (hysteresis ON does not")
        lines.append("shrink La(Fe,Si)13Hy's front share, and in most seeds increases it) is")
        lines.append("NOT an artifact of that smaller setting or that one seed. The EXACT")
        lines.append("88%->100% magnitude from the original run should still be read as one")
        lines.append("realization, not a universal constant -- see the per-seed table above")
        lines.append("for the actual seed-to-seed spread.")
    else:
        lines.append("RESULT: NOT STABLE. At least one seed's ON run had a LOWER La(Fe,Si)13Hy")
        lines.append("share than its own OFF run -- the direction of the original")
        lines.append("pop_size=32/n_gen=15/seed=1 finding does NOT hold universally across")
        lines.append("seeds at production settings. Treat the original 88%->100% finding as")
        lines.append("seed-dependent, not a settled result -- see the per-seed table above for")
        lines.append("which seed(s) disagree.")
    lines.append("")
    lines.append("Caveat carried over unchanged from this module's own honesty flag #2: the")
    lines.append("hysteresis_loss_J_per_kg VALUES this whole comparison depends on remain")
    lines.append("literature analogs for different exact compositions, not measurements of")
    lines.append("the calibrated compositions themselves -- this check resolves the")
    lines.append("NSGA-III-search-noise question, not the underlying-data-quality question.")

    if out_path:
        os.makedirs(os.path.dirname(out_path), exist_ok=True)
        with open(out_path, "w") as f:
            f.write("\n".join(lines) + "\n")

    print()
    print("\n".join(lines))
    if out_path:
        print(f"\nWrote {out_path}")

    return {"per_seed": per_seed, "stable": stable}


if __name__ == "__main__":
    run_hysteresis_sensitivity()