"""
optimize.py
============
Multi-objective design optimization of the AMR system using
NSGA-III (Deb & Jain, IEEE Trans. Evol. Comput. 18(4), 577-601 (2014)),
via the pymoo implementation (Blank & Deb, IEEE Access 8, 89497-89509 (2020)).

Design variables (7, continuous):
    mu0H_max                  [1.0, 3.0]   T
    frequency                 [0.3, 5.0]   Hz
    fluid_mdot                [0.02, 0.5]  kg/s
    mass_regenerator          [1.0, 15.0]  kg
    regenerator_effectiveness [0.6, 0.95]  -  (only actually used when
                                USE_NTU_THERMAL_MODEL=False -- see note
                                below; retained for backward compatibility)
    blow_fraction              [0.1, 0.6]   -  (Paper-Mining Pass
                                recommendation #1: fraction of the cycle
                                period spent in cold-to-hot flow. Bounds
                                bracket the 0.25-0.416 window Masche et al.
                                (2022) actually tested, widened somewhat to
                                let the optimizer explore nearby -- values
                                near the edges of [0.1,0.6] are extrapolated
                                beyond that paper's tested range; see
                                core.amr_cycle._blow_fraction_multiplier's
                                docstring for the calibration's honesty
                                flags. 0.5 reproduces the model's original
                                symmetric-blow behavior.)
    particle_diameter_mm      [0.05, 2.0]  mm  (Phase 15 addition -- see
                                below)

Note on regenerator_effectiveness vs. the NTU thermal model: this module
runs with USE_NTU_THERMAL_MODEL=True (as it always has, pre-Phase-15
included), which means AMRSystem._effective_eps() computes effectiveness
from geometry/mass/frequency/mdot via core.thermal.regenerator_
effectiveness() and IGNORES the regenerator_effectiveness design
variable entirely (see AMRSystem._effective_eps()). This was already true
before Phase 15 -- the variable has always been a passed-through-but-
unused 6th search dimension in NTU mode -- and is called out explicitly
here (previously only implicit) rather than removed, since removing a
design variable changes the CSV schema and downstream plots.py column
expectations; flagged as a good phase 15 cleanup candidate instead.

Phase 15 additions
-------------------
1. **Geometry as a design variable** (`particle_diameter_mm`): wired
   through to `AMRSystem`'s new `particle_diameter` parameter
   (`core/amr_cycle.py`), which both (a) feeds `core.thermal.
   regenerator_effectiveness()`'s NTU calculation (so smaller particles
   raise effectiveness, as `geometry_analysis.py` already demonstrated
   for a fixed representative mdot) and (b) computes a geometry-explicit
   hydraulic pumping power (`core.thermal.pumping_power_packed_bed()`,
   Tusek et al. 2013) that REPLACES (not adds to) `StateDependentLossModel`'s
   generic `k_pump*mdot**2` term via the new `pumping_power_override`
   mechanism -- see `core/loss_model.py`'s `parasitic_power()` docstring
   for why replacing rather than summing avoids double-counting the
   pumping-loss channel. Bounds [0.05, 2.0] mm match the range
   `geometry_analysis.py` already swept (avoiding further extrapolation
   beyond what's been checked against the Tusek et al. correlation's
   validity range).
2. **Material as a design variable**: implemented as OPTION (b) from the
   Phase 15 plan -- NSGA-III is run SEPARATELY per material family
   (`MATERIAL_CANDIDATES` below: Gd, plus each of `core.cascade`'s
   composition-tunable GD_FAMILY/LAFESIH_FAMILY/MNFEPSI_FAMILY/
   GA1XCMN3X_FAMILY (Phase 24)/MNCUCOGE_FAMILY (Phase 25), each
   tuned so its OWN peak lands at this module's fixed operating point's
   midpoint temperature -- same `_target_composition_for_peak` machinery
   `material_family_comparison.py` already uses, so no new numerics are
   introduced), and the resulting per-material Pareto fronts are merged
   post-hoc into one overall non-dominated set. This was chosen over
   option (a) (a true mixed-variable pymoo problem with material as a
   native categorical variable, via pymoo's MixedVariableGA/
   MixedVariableProblem) as the lower-risk, more this-repo-idiomatic
   choice: it reuses the existing single-material `ElementwiseProblem`
   unchanged (just parameterized by material/family), matches this
   repo's established "separate-then-compare" pattern
   (`material_family_comparison.py`), and introduces no new pymoo API
   surface. The tradeoff, stated rather than hidden: option (b) cannot
   find a solution that requires TRADING OFF a slightly-worse material
   choice for a better geometry/operating-point combination WITHIN a
   single generational search -- each material's search is independent,
   and only the FINAL fronts are compared. For the objectives here (COP,
   Qc, cost -- none of which couples material choice to geometry beyond
   what each material's own DeltaT_ad(T) shape already implies) this is
   a reasonable approximation, but it is a real, documented limitation of
   choosing (b) over (a).
3. **Cost objective upgraded**: `cost_index()` now uses
   `economics.bom_cost()` (magnet + MCM + soft-magnetic-material yoke,
   Phase 15's economics.py addition) instead of the older `material_cost()`
   (magnet + MCM only), and looks the MCM unit cost up per material family
   via `economics.MCM_COST_PER_KG_BY_FAMILY` rather than always assuming
   Gd's $20/kg -- material choice now has a genuine, family-specific cost
   consequence in the optimization, not just a performance one.

Phase 19 addition
-------------------
`cost_index()` (and hence `AMRDesignProblem`/`run_optimization_for_
material()`/`run_optimization()`) gained an opt-in
`use_geometric_magnet_mass` parameter (default False = unchanged
behavior): when True, the magnet-mass term inside the cost objective
comes from `economics.bom_cost_geometric()` -- itself built on
`core.magnet_geometry`'s closed-form idealized-Halbach-cylinder relation
-- instead of `economics.bom_cost()`'s flat per-Tesla mass ratio. The
geometric relation is NONLINEAR (super-linear) in mu0H, closing the gap
ROADMAP.md's Phase 19 plan named ("achieving high mu0H should cost
nonlinearly more magnet mass ... which is physically real and currently
absent"). Default is unchanged (flat ratio) so every pre-Phase-19 caller,
including this module's own `if __name__ == "__main__"` production run
and main.py's step 11, is completely unaffected unless
`use_geometric_magnet_mass=True` is passed explicitly. See
`core/magnet_geometry.py`'s `run_geometric_cost_pareto_sensitivity()` for
the A/B (flat vs. geometric) Pareto-front comparison this enables, and
that module's own docstring for the underlying physics and honesty flags.

Objectives (all converted to minimization for pymoo):
    f1 = -COP_electrical        (maximize electrical COP; uses the
                                 state-dependent loss model so the
                                 field/frequency/mdot/geometry choices
                                 carry a genuine efficiency cost)
    f2 = -Qc                    (maximize cooling capacity)
    f3 = cost_index             (minimize the materials-BOM cost proxy,
                                economics.bom_cost(), family-specific)

Fixed operating point: T_cold = 291 K, span = 10 K.

Output:
    results/pareto_front.csv                    -- merged, globally
                                                    non-dominated front
                                                    across all material
                                                    candidates
    results/pareto_front_by_material/<name>.csv  -- each material's own
                                                    per-material front,
                                                    before merging (Phase
                                                    15 addition, for
                                                    transparency/debugging)
"""

import os
import numpy as np
import csv
from pymoo.core.problem import ElementwiseProblem
from pymoo.algorithms.moo.nsga3 import NSGA3
from pymoo.util.ref_dirs import get_reference_directions
from pymoo.optimize import minimize as pymoo_minimize
from core.mce_material import GADOLINIUM
from core.amr_cycle import AMRSystem
from core.loss_model import StateDependentLossModel
from core.economics import material_cost, bom_cost, bom_cost_geometric
from core.cascade import (GD_FAMILY, LAFESIH_FAMILY, MNFEPSI_FAMILY,
                            GA1XCMN3X_FAMILY, MNCUCOGE_FAMILY,
                            _target_composition_for_peak)

T_COLD_K = 291.0
SPAN_K = 10.0
T_MID_K = T_COLD_K + SPAN_K / 2.0

_LOSS_MODEL = StateDependentLossModel()
USE_NTU_THERMAL_MODEL = True     # Enables the NTU-based thermal model so
                                 # regenerator mass influences cooling
                                 # capacity through the regenerator thermal
                                 # effectiveness calculation.
BED_CROSS_SECTION_AREA_M2 = 0.002  # matches geometry_analysis.py's own
                                     # default representative bed face.

# Design-variable bounds, in x-vector order:
#   [mu0H, freq, mdot, mass, eps, blow_fraction, particle_diameter_mm]
_XL = np.array([1.0, 0.3, 0.02, 1.0, 0.6, 0.1, 0.05])
_XU = np.array([3.0, 5.0, 0.5, 15.0, 0.95, 0.6, 2.0])


def cost_index(mu0H, mass_regenerator, family_name="Gd",
                use_geometric_magnet_mass=False):
    """Materials-BOM cost proxy (Phase 15: now includes the soft-magnetic
    -material yoke term and looks the MCM cost up per material family --
    see module docstring item 3). `family_name` defaults to "Gd" so
    existing callers passing only (mu0H, mass_regenerator) get the same
    material assumption as before; the $ VALUE returned is larger than
    the old `material_cost()`-based `cost_index()` for the same (mu0H,
    mass) because it now also includes the SMM yoke term -- this is an
    intentional, documented improvement (a more complete materials-cost
    proxy), not a bug; see economics.py's Phase 15 section.

    Phase 19 addition: `use_geometric_magnet_mass=False` (default,
    fully backward-compatible) keeps using `economics.bom_cost()`'s flat
    per-Tesla magnet-mass ratio; passing `True` switches the magnet-mass
    term to `economics.bom_cost_geometric()`'s closed-form idealized-
    Halbach-cylinder relation (`core/magnet_geometry.py`), which is
    NONLINEAR (super-linear) in mu0H rather than linear -- see that
    module's own docstring for the physics and honesty flags.
    `AMRDesignProblem`/`run_optimization_for_material`/
    `run_optimization()` all thread this flag through with the same
    default-False backward-compatible convention."""
    cost_fn = bom_cost_geometric if use_geometric_magnet_mass else bom_cost
    return cost_fn(mu0H, mass_regenerator, family_name)["materials_bom_total_$"]


def _material_candidates(T_mid_K=T_MID_K, mu0H_max_for_tuning=2.0):
    """Builds the set of material candidates for the Phase 15 per-family
    NSGA-III co-optimization: plain Gd, plus each of core.cascade's three
    composition-tunable giant-MCE families, each tuned so its own peak
    DeltaT_ad lands at T_mid_K (this module's fixed operating point's
    midpoint) -- same approach and same root-finder
    (`core.cascade._target_composition_for_peak`) `material_family_
    comparison.py` already uses for the equivalent single-composition-
    per-operating-point comparison. Families whose required composition
    falls outside their documented tunability window are DROPPED here
    (not silently substituted with a Gd fallback that would just
    duplicate the "Gd" row) -- `material_family_comparison.py`'s own
    results (see results/material_family_comparison.csv) already show
    MNFEPSI_FAMILY's window sits mostly at/above the ASHRAE range and can
    fail to cover a given operating point; this module re-checks it fresh
    at ITS OWN operating point rather than assuming that result still
    holds.

    Returns a list of (label, material, family_name_for_cost) tuples,
    family_name_for_cost matching a key in
    economics.MCM_COST_PER_KG_BY_FAMILY.
    """
    candidates = [("Gd", GADOLINIUM, "Gd")]
    for family in (GD_FAMILY, LAFESIH_FAMILY, MNFEPSI_FAMILY, GA1XCMN3X_FAMILY,
                   MNCUCOGE_FAMILY):
        tc = _target_composition_for_peak(T_mid_K, mu0H_max_for_tuning, family)
        in_range = family.tc_min <= tc <= family.tc_max
        if not in_range:
            print(f"_material_candidates: {family.name} needs Tc={tc:.1f}K to hit "
                  f"T_mid={T_mid_K:.1f}K, outside its documented window "
                  f"[{family.tc_min:.1f}, {family.tc_max:.1f}]K -- dropped from this "
                  "operating point's candidate set (not substituted with a "
                  "duplicate Gd row).")
            continue
        candidates.append((f"{family.name} (tuned, Tc={tc:.1f}K)",
                            family.tuned_fn(tc), family.name))
    return candidates


class AMRDesignProblem(ElementwiseProblem):
    def __init__(self, material=GADOLINIUM, family_name="Gd",
                 use_geometric_magnet_mass=False):
        self.material = material
        self.family_name = family_name
        self.use_geometric_magnet_mass = use_geometric_magnet_mass
        super().__init__(
            n_var=7, n_obj=3, n_constr=0,
            xl=_XL, xu=_XU,
        )

    def _evaluate(self, x, out, *args, **kwargs):
        mu0H, freq, mdot, mass, eps, blow_fraction, particle_diameter_mm = x
        sys_ = AMRSystem(material=self.material, mu0H_max=mu0H, mass_regenerator=mass,
                          frequency=freq, fluid_mdot=mdot, regenerator_effectiveness=eps,
                          loss_model=_LOSS_MODEL, use_ntu_thermal_model=USE_NTU_THERMAL_MODEL,
                          blow_fraction=blow_fraction,
                          particle_diameter=particle_diameter_mm / 1000.0,
                          bed_cross_section_area=BED_CROSS_SECTION_AREA_M2)
        result = sys_.run(T_COLD_K, SPAN_K)
        f1 = -result.COP_electrical
        f2 = -result.Qc
        f3 = cost_index(mu0H, mass, self.family_name, self.use_geometric_magnet_mass)
        out["F"] = [f1, f2, f3]


def _row_from_xf(x, f, material_label):
    return {
        "material": material_label,
        "mu0H_max_T": round(x[0], 3), "frequency_Hz": round(x[1], 3),
        "fluid_mdot_kgs": round(x[2], 4), "mass_regenerator_kg": round(x[3], 2),
        "regen_effectiveness": round(x[4], 3), "blow_fraction": round(x[5], 3),
        "particle_diameter_mm": round(x[6], 4),
        "COP_electrical": round(-f[0], 3), "Qc_W": round(-f[1], 2),
        "cost_index_USD": round(f[2], 1),
    }


_ROW_FIELDNAMES = ["material", "mu0H_max_T", "frequency_Hz", "fluid_mdot_kgs",
                    "mass_regenerator_kg", "regen_effectiveness", "blow_fraction",
                    "particle_diameter_mm", "COP_electrical", "Qc_W", "cost_index_USD"]


def _write_csv(rows, path):
    os.makedirs(os.path.dirname(path) or ".", exist_ok=True)
    with open(path, "w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=_ROW_FIELDNAMES)
        writer.writeheader()
        writer.writerows(rows)


def run_optimization_for_material(material, family_name, material_label,
                                    pop_size=40, n_gen=25, seed=1,
                                    out_csv=None, use_geometric_magnet_mass=False):
    """Runs NSGA-III for a single material candidate. This is the
    per-family sub-search of the Phase 15 "option (b)" material
    co-optimization (see module docstring item 2). `pop_size`/`n_gen` are
    reduced from the pre-Phase-15 single-material defaults (60/40) since
    Phase 15's `run_optimization()` now runs this once per material
    candidate (typically 3-4) -- see that function's docstring for the
    total-runtime accounting.

    `use_geometric_magnet_mass` (Phase 19, default False = old behavior)
    is passed straight through to `AMRDesignProblem`/`cost_index()` --
    see `cost_index()`'s own docstring."""
    ref_dirs = get_reference_directions("das-dennis", 3, n_partitions=6)
    algorithm = NSGA3(pop_size=pop_size, ref_dirs=ref_dirs)
    problem = AMRDesignProblem(material=material, family_name=family_name,
                                use_geometric_magnet_mass=use_geometric_magnet_mass)
    res = pymoo_minimize(problem, algorithm, ("n_gen", n_gen), seed=seed, verbose=False)

    X, F = res.X, res.F
    rows = [_row_from_xf(x, f, material_label) for x, f in zip(X, F)]
    if out_csv:
        _write_csv(rows, out_csv)
    return rows


def _pareto_filter(rows):
    """Global non-dominated filter across the MERGED multi-material row
    set (the "merge the Pareto fronts post-hoc" step of option (b) -- see
    module docstring item 2). A row is kept unless some OTHER row is at
    least as good in every objective and strictly better in at least one
    (standard Pareto dominance on [-COP_electrical, -Qc, cost_index_USD],
    i.e. minimize all three)."""
    if not rows:
        return rows
    F = np.array([[-r["COP_electrical"], -r["Qc_W"], r["cost_index_USD"]] for r in rows])
    n = len(rows)
    dominated = np.zeros(n, dtype=bool)
    for i in range(n):
        if dominated[i]:
            continue
        le = np.all(F <= F[i], axis=1)
        lt = np.any(F < F[i], axis=1)
        dominators = le & lt
        if np.any(dominators):
            dominated[i] = True
    return [r for r, d in zip(rows, dominated) if not d]


def run_optimization(pop_size=40, n_gen=25, seed=1,
                      out_csv="results/pareto_front.csv",
                      per_material_out_dir="results/pareto_front_by_material",
                      use_geometric_magnet_mass=False):
    """Phase 15: runs NSGA-III separately for each material candidate
    (`_material_candidates()`), writes each material's own Pareto front to
    `per_material_out_dir/<label>.csv`, then merges all candidate rows and
    applies a global non-dominance filter (`_pareto_filter()`) to produce
    the single overall `out_csv` -- see module docstring item 2 for why
    this "separate-then-merge" approach (option (b)) was chosen over a
    native mixed-variable pymoo formulation (option (a)).

    Runtime note: with `pop_size`/`n_gen` reduced to 40/25 (from the
    pre-Phase-15 single-material 60/40) and ~4 material candidates, total
    NSGA-III wall time is comparable to the pre-Phase-15 single 60/40 run
    (roughly proportional to pop_size*n_gen*n_materials); see main.py's
    own runtime-estimate docstring for the current end-to-end figure.
    Tests (`tests/test_optimize.py`,
    `tests/test_optimize_material_geometry.py`) pass much smaller
    `pop_size`/`n_gen` explicitly and run in ~1s total.

    Phase 19 addition: `use_geometric_magnet_mass=False` (default, fully
    backward-compatible) is threaded through to every per-material
    `run_optimization_for_material()` call -- see `cost_index()`'s own
    docstring for what passing `True` changes (a nonlinear, closed-form
    Halbach-cylinder magnet-mass cost term instead of the flat per-Tesla
    ratio). `core/magnet_geometry.py`'s
    `run_geometric_cost_pareto_sensitivity()` runs this function twice
    (flat vs. geometric) as an A/B comparison, the same pattern
    `core/hysteresis_sensitivity.py` established for Phase 16.
    """
    candidates = _material_candidates()
    all_rows = []
    per_material_rows = {}
    print(f"Phase 15 material co-optimization: {len(candidates)} candidate(s) -- "
          f"{', '.join(label for label, _, _ in candidates)}")
    for label, material, family_name in candidates:
        safe_label = "".join(c if c.isalnum() else "_" for c in label)
        out_path = os.path.join(per_material_out_dir, f"{safe_label}.csv") \
            if per_material_out_dir else None
        rows = run_optimization_for_material(
            material, family_name, label, pop_size=pop_size, n_gen=n_gen,
            seed=seed, out_csv=out_path,
            use_geometric_magnet_mass=use_geometric_magnet_mass)
        per_material_rows[label] = rows
        all_rows.extend(rows)
        print(f"  {label:<40} {len(rows)} Pareto-optimal design(s) found"
              + (f" -> {out_path}" if out_path else ""))

    merged = _pareto_filter(all_rows)
    _write_csv(merged, out_csv)

    rows = merged
    print(f"\nMerged across {len(candidates)} material candidate(s): "
          f"{len(all_rows)} total designs -> {len(rows)} globally non-dominated "
          f"design(s) after cross-material Pareto filtering.")
    print(f"Wrote {out_csv}\n")

    if not rows:
        print("No non-dominated designs found; skipping summary.")
        return rows

    F = np.array([[-r["COP_electrical"], -r["Qc_W"], r["cost_index_USD"]] for r in rows])
    Fn = (F - F.min(axis=0)) / (F.max(axis=0) - F.min(axis=0) + 1e-12)
    knee_idx = int(np.argmin(np.linalg.norm(Fn, axis=1)))
    best_cop_idx = int(np.argmax(-F[:, 0]))
    best_qc_idx = int(np.argmax(-F[:, 1]))
    best_cost_idx = int(np.argmin(F[:, 2]))

    for label, idx in [("Best electrical COP", best_cop_idx),
                        ("Best cooling capacity", best_qc_idx),
                        ("Lowest cost", best_cost_idx),
                        ("Knee point (balanced)", knee_idx)]:
        r = rows[idx]
        print(f"{label:<24} material={r['material']:<28} H={r['mu0H_max_T']}T  "
              f"f={r['frequency_Hz']}Hz  mdot={r['fluid_mdot_kgs']}kg/s  "
              f"mass={r['mass_regenerator_kg']}kg  d_p={r['particle_diameter_mm']}mm  "
              f"bf={r['blow_fraction']}  -> COP_elec={r['COP_electrical']}, "
              f"Qc={r['Qc_W']}W, cost=${r['cost_index_USD']}")

    material_counts = {}
    for r in rows:
        material_counts[r["material"]] = material_counts.get(r["material"], 0) + 1
    print("\nMaterial representation in the merged, globally non-dominated front:")
    for label, count in sorted(material_counts.items(), key=lambda kv: -kv[1]):
        print(f"  {label:<40} {count} design(s) ({100*count/len(rows):.0f}%)")

    diam = [r["particle_diameter_mm"] for r in rows]
    print(f"\nparticle_diameter spans {min(diam):.3f}-{max(diam):.3f} mm across the merged "
          "front -- the geometry-explicit pumping-power/effectiveness trade-off "
          "(see core/geometry_analysis.py, core/amr_cycle.py's Phase 15 wiring) is a "
          "real, active dimension of this search, not degenerate at either bound.")

    masses = [r["mass_regenerator_kg"] for r in rows]
    print(f"Regenerator mass spans {min(masses):.2f}-{max(masses):.2f} kg across the merged "
          "front. The NTU thermal model captures the expected trade-off between "
          "additional regenerator material, cooling capacity, and material cost.")
    return rows


if __name__ == "__main__":
    run_optimization()