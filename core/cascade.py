"""
cascade.py
===========
Multi-stage (cascade) AMR model for extending the operating temperature
span beyond the limits of a single-stage system.

Cascade concept: N single-stage AMR modules in series, each handling an
equal share of the total span (Th_total - Tc_total)/N, analogous to
cascade vapor-compression refrigeration. Stage 1 (coldest) absorbs the
data-center heat load Qc at T_cold; its heat rejection becomes the input
to Stage 2, and so on until the final stage rejects heat at T_hot.
Because heat is transferred through the stages in series, each stage must
reject the same cooling load (steady state, neglecting inter-stage losses):

    W_total = Σ W_i(Qc, span/N)
    COP_cascade = Qc / W_total

run_cascade()/compare_staging() below assume identical regenerator stages
(all Gd, or all Gd5Si2Ge2). run_graded_cascade()/compare_graded_cascade()
implement the more advanced Curie-graded variant (ROADMAP.md Phase 7 open
item): each stage uses a hypothetical composition-tuned Gd5(SixGe1-x)4(-Ga)
material whose own peak MCE effect is matched to that stage's local
operating temperature, checked against the literature-documented
composition-tunability range and against an independent direct-measurement
validation (see core/first_order_mce.py and core/giguere_validation.py).

Phase 9 addendum: the Curie-graded machinery was generalized to a pluggable
`GradedFamily` (see below) rather than being hardcoded to the Gd5(SixGe1-x)4
family, so the SAME grading mechanism can be applied to the La(Fe,Si)13Hy
family added in Phase 9 (core.first_order_mce.lafesih_composition_tuned_material)
-- specifically to test whether a 6-layer Curie-graded La(Fe,Si)13Hy bed can
reproduce the REAL Astronautics_rotary_2014 benchmark device (which is
exactly such a bed), closing the gap validation_system.py's Phase 9
single-layer approximation flagged as a "natural next step, not done here".
GD_FAMILY below reproduces the original Gd5(SixGe1-x)4 behavior exactly
(default argument, so existing calls are unaffected); LAFESIH_FAMILY is new.
See validate_astronautics_graded_bed() at the bottom of this module for the
actual comparison against Jacobs et al. (2014)'s reported numbers.
"""

import numpy as np
import csv
import os
import time
import logging
import threading
import contextlib
from concurrent.futures import ProcessPoolExecutor
from dataclasses import dataclass, field
from typing import Callable
from scipy.optimize import brentq
from core.mce_material import GADOLINIUM
from core.amr_cycle import AMRSystem
from core.loss_model import StateDependentLossModel
from core.baseline_cooling import vapor_compression_cop, liquid_cooling_cop
from core.first_order_mce import (composition_tuned_material,
                                    GIANT_MCE_TC_MIN_K, GIANT_MCE_TC_MAX_K,
                                    lafesih_composition_tuned_material,
                                    LAFESIH_TC_MIN_K, LAFESIH_TC_MAX_K,
                                    mnfepsi_composition_tuned_material,
                                    MNFEPSI_TC_MIN_K, MNFEPSI_TC_MAX_K,
                                    mncucoge_composition_tuned_material,
                                    MNCUCOGE_TC_MIN_K, MNCUCOGE_TC_MAX_K,
                                    GD5SI2GE2_FIRST_ORDER, LAFESIH_FIRST_ORDER,
                                    MNFEPSI_FIRST_ORDER, MNCUCOGE_FIRST_ORDER)
from core.antiperovskite_material import (ga1xcmn3x_composition_tuned_material,
                                            GA1XCMN3X_TC_MIN_K, GA1XCMN3X_TC_MAX_K,
                                            GA1XCMN3X_REF)

_LOSS_MODEL = StateDependentLossModel()
USE_NTU_THERMAL_MODEL = True


@dataclass
class GradedFamily:
    """Everything run_graded_cascade() needs to grade a bed from a single
    composition-tunable material family. reference_material is only used to
    seed _target_composition_for_peak()'s fixed-point iteration (its own
    peak-vs-Tc offset is the starting guess); tuned_fn(Tc_target_K) must
    return a FirstOrderMCEMaterial and raise ValueError outside
    [tc_min, tc_max] (both composition_tuned_material() and
    lafesih_composition_tuned_material() already do this)."""
    name: str
    tuned_fn: Callable[[float], object]
    tc_min: float
    tc_max: float
    reference_material: object
    fallback_material: object = field(default_factory=lambda: GADOLINIUM)


GD_FAMILY = GradedFamily(
    name="Gd5(SixGe1-x)4(-Ga)",
    tuned_fn=lambda Tc: composition_tuned_material(Tc, apply_giguere_correction=True),
    tc_min=GIANT_MCE_TC_MIN_K, tc_max=GIANT_MCE_TC_MAX_K,
    reference_material=GD5SI2GE2_FIRST_ORDER, fallback_material=GADOLINIUM,
)

LAFESIH_FAMILY = GradedFamily(
    name="La(Fe,Si)13Hy",
    tuned_fn=lafesih_composition_tuned_material,
    tc_min=LAFESIH_TC_MIN_K, tc_max=LAFESIH_TC_MAX_K,
    reference_material=LAFESIH_FIRST_ORDER, fallback_material=GADOLINIUM,
)

# Paper-Mining Pass recommendation #3: (Mn,Fe)2(P,Si), a third pluggable
# giant-MCE family alongside GD_FAMILY/LAFESIH_FAMILY. Unlike those two
# families' Tc windows, this one (295.3-331.2K, see MNFEPSI_TC_MIN_K/_MAX_K
# in core/first_order_mce.py) is directly measured across five real
# compositions and sits almost entirely AT OR ABOVE the ASHRAE data-center
# supply range -- the opposite tension from GD_FAMILY, whose documented
# giant-MCE ceiling sits just below it.
MNFEPSI_FAMILY = GradedFamily(
    name="(Mn,Fe)2(P,Si)",
    tuned_fn=mnfepsi_composition_tuned_material,
    tc_min=MNFEPSI_TC_MIN_K, tc_max=MNFEPSI_TC_MAX_K,
    reference_material=MNFEPSI_FIRST_ORDER, fallback_material=GADOLINIUM,
)

# Phase 24: fourth pluggable family, Ga1-xCMn3+x (see
# core/antiperovskite_material.py for the full literature-verification
# writeup -- this is NOT the material a user-supplied "10 materials"
# document named ("GaCMn3" at ~296K), which turned out on checking to be
# a wrong Tc for the stoichiometric compound; this is the real,
# literature-measured, second-order, hysteresis-free, Tc-tunable
# composition series that document was gesturing toward).
# reference_material's Tc/peak-vs-Tc offset is only used to SEED
# _target_composition_for_peak()'s root-find below -- for a second-order
# mean-field material built via with_Tc(), that offset is expected to be
# ~0 (unlike the first-order Landau families' own +10-11K offset), and
# the root-finder recomputes the true offset numerically regardless, so
# using GA1XCMN3X_REF here (rather than needing a separate calibrated
# "first_order"-style reference object) is not a special case.
GA1XCMN3X_FAMILY = GradedFamily(
    name="Ga1-xCMn3+x",
    tuned_fn=ga1xcmn3x_composition_tuned_material,
    tc_min=GA1XCMN3X_TC_MIN_K, tc_max=GA1XCMN3X_TC_MAX_K,
    reference_material=GA1XCMN3X_REF, fallback_material=GADOLINIUM,
)

# Phase 25: fifth pluggable family, Mn1-xCuxCoGe (see the block comment
# above core/first_order_mce.py's MNCUCOGE_FIRST_ORDER for the full
# literature-verification writeup). Unlike GA1XCMN3X_FAMILY (a SECOND-
# order family reusing mce_material.py), this is a first-order Landau
# family using the SAME machinery as GD_FAMILY/LAFESIH_FAMILY/
# MNFEPSI_FAMILY below -- so, unlike GA1XCMN3X_FAMILY, its
# reference_material really is the same kind of calibrated-to-a-digitized-
# DeltaS_M-target object the other three first-order families use.
MNCUCOGE_FAMILY = GradedFamily(
    name="Mn1-xCuxCoGe",
    tuned_fn=mncucoge_composition_tuned_material,
    tc_min=MNCUCOGE_TC_MIN_K, tc_max=MNCUCOGE_TC_MAX_K,
    reference_material=MNCUCOGE_FIRST_ORDER, fallback_material=GADOLINIUM,
)


# --------------------------------------------------------------------------
# Phase 16: process-pool parallelization helpers.
#
# GradedFamily objects themselves are never sent across a process boundary
# (GD_FAMILY's tuned_fn is a lambda, and validate_astronautics_graded_bed()'s
# apply_correction branch builds an even-less-picklable closure -- neither
# survives pickle). Instead, only a short family NAME crosses the boundary;
# each worker process re-imports core.cascade normally (as
# ProcessPoolExecutor's own bootstrapping already requires) and so already
# has its own GD_FAMILY/LAFESIH_FAMILY/MNFEPSI_FAMILY globals to resolve the
# name against. _family_name() returns None for any family this module
# doesn't recognize by identity (e.g. a custom closure-based family) --
# every parallel entry point below checks for that and falls back to the
# original sequential behavior in that case, so nothing that used to work
# stops working; it just doesn't get the parallel speedup.
# --------------------------------------------------------------------------

def _pool_worker_init():
    """ProcessPoolExecutor `initializer`, run once per worker process
    BEFORE that process handles any task. This is a defensive second layer
    only -- see _single_threaded_blas_env() below for the fix that
    actually matters and why setting these vars here alone is too late
    (merely unpickling a reference to THIS function already re-imports
    core.cascade, and therefore numpy, in the child first)."""
    for var in ("OPENBLAS_NUM_THREADS", "OMP_NUM_THREADS", "MKL_NUM_THREADS",
                "NUMEXPR_NUM_THREADS", "VECLIB_MAXIMUM_THREADS"):
        os.environ[var] = "1"


@contextlib.contextmanager
def _single_threaded_blas_env():
    """Temporarily sets OPENBLAS_NUM_THREADS / OMP_NUM_THREADS / etc. to
    '1' in THIS (parent) process's environment, for exactly as long as a
    ProcessPoolExecutor's child processes are being spawned/used below.

    This is the actual fix for 'OpenBLAS error: Memory allocation still
    failed after 10 retries, giving up.' Without it, every worker process
    independently spins up its own OpenBLAS thread pool sized to the
    machine's full logical core count; with max_workers itself already
    close to that core count, you get (workers) x (cores) live BLAS
    threads all contending for the same physical cores and OpenBLAS's
    internal thread-local memory arenas -- CPU/memory oversubscription,
    not a bug in the cascade math itself.

    Why the PARENT's environment, and not something set from inside a
    function this module defines: with the 'spawn' start method (the
    default, and the only option on Windows -- where this was actually
    observed), a new worker process's OS-level environment block is a copy
    of the PARENT's os.environ taken at the moment the child process is
    created, before that child's own Python interpreter -- let alone this
    module's `import numpy as np` -- ever runs. Setting these vars here,
    before ProcessPoolExecutor(...) spawns any child, is what actually
    reaches OpenBLAS's library-load-time thread-pool sizing in each child.
    Restored on exit so later pipeline stages (RSM, NSGA-III, ...) keep
    this process's normal BLAS threading."""
    _vars = ("OPENBLAS_NUM_THREADS", "OMP_NUM_THREADS", "MKL_NUM_THREADS",
             "NUMEXPR_NUM_THREADS", "VECLIB_MAXIMUM_THREADS")
    _prev = {v: os.environ.get(v) for v in _vars}
    for v in _vars:
        os.environ[v] = "1"
    try:
        yield
    finally:
        for v, val in _prev.items():
            if val is None:
                os.environ.pop(v, None)
            else:
                os.environ[v] = val


logger = logging.getLogger(__name__)

# --------------------------------------------------------------------------
# Phase 31: process-pool HANG robustness fix.
#
# What was actually wrong (found by directly running this project's own
# test suite in a sandboxed/restricted container, not by inspection alone):
# `tests/test_cascade.py::test_magqueen_mass_sensitivity_parallel_matches_sequential`
# and its `test_cooltech_mass_sensitivity_parallel_matches_sequential`
# analog hung INDEFINITELY (confirmed via manual per-test wall-clock
# bisection, >60s with no output and no exception) in an environment where
# ProcessPoolExecutor's worker processes are created but never actually
# complete a submitted task. Every existing `except Exception:` guard
# around a pool block in this module is a real, working fallback for pool
# CREATION failing outright (which raises immediately, e.g.
# PermissionError) -- but none of them protect against a pool that is
# created successfully and then simply never returns a result, because
# NONE of `pool.map(...)`, `future.result()`, or `pool.shutdown(wait=True)`
# were ever called with a timeout anywhere in this module. In that failure
# mode nothing raises, so the bare `except Exception:` never fires -- the
# calling process (or `pytest`, or a live `main.py` run, or, worst case, a
# grading/demo run under time pressure) just hangs forever. This is a
# genuinely different failure mode from the ones this module's own Phase 16
# docstring above already anticipated (unrecognized family objects,
# single-cell sweeps, outright pool-creation failure) -- it was not
# previously covered by any fallback.
#
# The three helpers below are used at every remaining
# ProcessPoolExecutor call site in this module and share one rule: ANY
# pool failure -- an exception, OR simply exceeding a wall-clock timeout --
# is treated identically, and always results in returning `None` (for the
# two "do the work" helpers) or in giving up on waiting (for the shutdown
# helper), so the CALLING process itself can never hang past
# `timeout_s` (+ a small constant shutdown grace period) no matter what the
# worker processes are doing. Every call site that uses these already has
# its own pre-existing, already-tested sequential fallback path, so a
# `None` return here changes nothing about correctness -- only about
# whether the parallel speedup was actually obtained this time.
# --------------------------------------------------------------------------

_DEFAULT_POOL_TIMEOUT_S = 120

# `run_graded_cascade()` (via _pool_map_or_none below) is called from
# INSIDE a `scipy.optimize.brentq` root-finding loop by every
# validate_*_graded_bed() calibration function (each brentq call re-runs
# the whole cascade, including its Tc-target pool.map(), once per
# iteration -- commonly 10-40+ times to converge). A naive per-call
# timeout would cost up to `timeout_s` seconds on EVERY iteration if the
# pool is chronically broken/slow, i.e. it would make a stuck pool cost
# `timeout_s * n_brentq_iterations` instead of a single bounded delay --
# the opposite of what this fix is supposed to guarantee. So once a given
# executor instance is observed to fail or time out, it is marked
# "poisoned" (a plain attribute set directly on the executor object) and
# every subsequent _pool_map_or_none() call against that SAME instance
# short-circuits to `None` immediately, without retrying the pool at all,
# for the rest of that instance's lifetime. A fresh pool built later
# (e.g. by the next call to validate_astronautics_graded_bed()) starts
# unpoisoned and gets its own single chance.
_MAP_TIMEOUT_S = 20


def _pool_map_or_none(executor, worker_fn, args_list, timeout_s=_MAP_TIMEOUT_S):
    """`list(executor.map(worker_fn, args_list, timeout=timeout_s))`, but
    ANY exception (including `concurrent.futures.TimeoutError` once the
    batch's wall-clock time is exceeded) is caught and turned into a
    logged warning + a `None` return, instead of propagating or hanging.
    See the "poisoned executor" comment directly above this function for
    why a broken executor is only ever retried once, not on every call.
    Does not shut the executor down -- callers that own the executor's
    lifecycle (rather than borrowing one built elsewhere) do that
    themselves via `_safe_pool_shutdown()` below."""
    if getattr(executor, "_phase31_poisoned", False):
        return None
    try:
        return list(executor.map(worker_fn, args_list, timeout=timeout_s))
    except Exception as exc:
        logger.warning(
            "Phase 31: ProcessPoolExecutor.map() failed or exceeded its "
            "%ss timeout (%s: %s) -- falling back to the sequential path "
            "for the rest of this executor's lifetime.",
            timeout_s, type(exc).__name__, exc)
        try:
            executor._phase31_poisoned = True
        except Exception:
            pass  # best-effort marker only; worst case we retry once more
        return None


def _pool_submit_all_or_none(pool, worker_fn, items, timeout_s=_DEFAULT_POOL_TIMEOUT_S):
    """Submits `worker_fn(item)` for every `item` in `items` to `pool`, then
    waits up to `timeout_s` seconds TOTAL (not per-future) for every result.
    Returns the list of results in the same order as `items`, or `None` if
    any future raised or the overall timeout was exceeded. This is the
    direct fix for the two mass-sensitivity hangs described above, which
    used bare `future.result()` (no timeout at all) on exactly this
    submit-then-collect pattern."""
    deadline = time.monotonic() + timeout_s
    futures = [pool.submit(worker_fn, item) for item in items]
    try:
        results = []
        for fut in futures:
            remaining = max(0.0, deadline - time.monotonic())
            results.append(fut.result(timeout=remaining))
        return results
    except Exception as exc:
        logger.warning(
            "Phase 31: a pool.submit()/future.result() call failed or "
            "exceeded its overall %ss timeout (%s: %s) -- falling back to "
            "the sequential path.", timeout_s, type(exc).__name__, exc)
        for fut in futures:
            fut.cancel()  # best-effort only; no effect on already-running tasks
        return None


def _safe_pool_shutdown(pool, wait_s=5):
    """Shuts a ProcessPoolExecutor down without ever blocking the calling
    process indefinitely. `pool.shutdown(wait=True)` has no timeout
    parameter of its own and can itself hang forever if any worker process
    never completes/exits -- the same failure mode the two helpers above
    guard against for in-flight work, but for the cleanup step instead.
    This does a best-effort wait of `wait_s` seconds on a background daemon
    thread, then gives up and returns either way. It does not guarantee
    the pool's OS-level child processes are reaped, only that THIS process
    does not hang -- which is the actual bug being fixed (see this
    section's module-level comment above)."""
    if pool is None:
        return
    t = threading.Thread(target=pool.shutdown, kwargs={"wait": True}, daemon=True)
    t.start()
    t.join(timeout=wait_s)
    if t.is_alive():
        logger.warning(
            "Phase 31: pool.shutdown(wait=True) did not complete within "
            "%ss -- abandoning the wait (worker processes may still be "
            "running in the background; this process itself will not "
            "hang).", wait_s)
        try:
            pool.shutdown(wait=False, cancel_futures=True)
        except Exception:
            pass


def _family_name(family):
    if family is None or family is GD_FAMILY:
        return "GD"
    if family is LAFESIH_FAMILY:
        return "LAFESIH"
    if family is MNFEPSI_FAMILY:
        return "MNFEPSI"
    if family is GA1XCMN3X_FAMILY:
        return "GA1XCMN3X"
    if family is MNCUCOGE_FAMILY:
        return "MNCUCOGE"
    return None


def _resolve_family(family_name):
    if family_name == "GD":
        return GD_FAMILY
    if family_name == "LAFESIH":
        return LAFESIH_FAMILY
    if family_name == "MNFEPSI":
        return MNFEPSI_FAMILY
    if family_name == "GA1XCMN3X":
        return GA1XCMN3X_FAMILY
    if family_name == "MNCUCOGE":
        return MNCUCOGE_FAMILY
    return None


def _stage_target_worker(args):
    """Top-level (picklable) worker for a single stage's Curie-target
    root-search. Used only when run_graded_cascade() is given a live
    `executor` -- see that function's docstring. Each stage's target
    depends only on its own T_mid, so running these in parallel changes
    nothing about the result, only how long it takes to get it."""
    T_mid_stage, mu0H_max, family_name = args
    family = _resolve_family(family_name) or GD_FAMILY
    return _target_composition_for_peak(T_mid_stage, mu0H_max, family)


def _compare_graded_cascade_cell(args):
    """Top-level (picklable) worker for compare_graded_cascade()'s
    per-(span, n_stages) sweep cells. Must stay at module level so
    ProcessPoolExecutor can pickle a reference to it (a nested/closure
    function cannot be pickled this way)."""
    T_cold_K, span, n, mass_per_stage, family_name = args
    family = _resolve_family(family_name)
    res = run_graded_cascade(T_cold_K, span, n, mass_per_stage=mass_per_stage,
                              family=family)
    return span, n, res


def run_cascade(T_cold_K, total_span_K, n_stages, material=None, mu0H_max=2.0,
                 mass_per_stage=2.0, frequency=1.0, fluid_mdot=0.08,
                 regenerator_effectiveness=0.85):
    """Runs n_stages identical AMR modules in series, each covering
    total_span_K/n_stages, all passing the same Qc through in steady state
    (Qc is set by the coldest/first stage's capacity at its local span)."""
    if material is None:
        material = GADOLINIUM
    span_per_stage = total_span_K / n_stages
    T_local = T_cold_K
    # First stage sets the deliverable Qc (bottleneck of the chain)
    stage1 = AMRSystem(material=material, mu0H_max=mu0H_max,
                        mass_regenerator=mass_per_stage, frequency=frequency,
                        fluid_mdot=fluid_mdot, regenerator_effectiveness=regenerator_effectiveness,
                        loss_model=_LOSS_MODEL, use_ntu_thermal_model=USE_NTU_THERMAL_MODEL)
    r1 = stage1.run(T_local, span_per_stage)
    Qc_target = r1.Qc
    if Qc_target <= 0:
        return {"n_stages": n_stages, "span_per_stage_K": span_per_stage,
                "Qc_W": 0.0, "W_total_W": np.nan, "COP_cascade": 0.0,
                "feasible": False}

    W_total = 0.0
    for _i in range(n_stages):
        stage = AMRSystem(material=material, mu0H_max=mu0H_max,
                           mass_regenerator=mass_per_stage, frequency=frequency,
                           fluid_mdot=fluid_mdot, regenerator_effectiveness=regenerator_effectiveness,
                           loss_model=_LOSS_MODEL, use_ntu_thermal_model=USE_NTU_THERMAL_MODEL)
        # each stage handles the same Qc_target at its local span; back out
        # the required work by re-running at span_per_stage and scaling mdot
        # if needed so Qc matches Qc_target (steady-state series constraint)
        r_i = stage.run(T_local, span_per_stage)
        if r_i.Qc > 0:
            scale = Qc_target / r_i.Qc
            W_i = (r_i.W_mag + r_i.W_parasitic) * scale
        else:
            W_i = np.inf
        W_total += W_i
        T_local += span_per_stage

    COP_cascade = Qc_target / W_total if W_total > 0 else 0.0
    return {"n_stages": n_stages, "span_per_stage_K": span_per_stage,
            "Qc_W": round(Qc_target, 1), "W_total_W": round(W_total, 1),
            "COP_cascade": round(COP_cascade, 3), "feasible": True}


def run_explicit_material_cascade(T_cold_K, total_span_K, materials, mu0H_max=2.0,
                                   mass_per_stage=2.0, frequency=1.0, fluid_mdot=0.08,
                                   regenerator_effectiveness=0.85, cycle_type="brayton",
                                   shared_hardware=False):
    """Same series-cascade mechanics as run_cascade()/run_graded_cascade()
    above, but for the case where the REAL per-stage material composition
    is already known from a paper (not searched for, unlike
    run_graded_cascade()'s composition-tuning against a hypothetical
    GradedFamily). `materials` is an explicit, ordered list (coldest stage
    first) of MagnetocaloricMaterial instances, one per stage -- e.g.
    [GADOLINIUM.with_Tc(272.15), GADOLINIUM.with_Tc(282.65), ...] for a
    reported Gd/Gd-Y multi-layer bed. n_stages = len(materials).
    mass_per_stage may be a single float (split evenly across all stages)
    or a list/tuple of per-stage masses matching len(materials), for the

    `shared_hardware` (default False, reproduces every pre-existing call's
    behavior exactly): when True, changes how the N per-stage
    AMRCycleResult.W_parasitic terms are combined into W_total. The
    default (False) sums each stage's OWN independently-computed
    W_eddy+W_pump+W_base (scaled to the bottleneck Qc), which is the right
    model when the "stages" really are N separate physical machines (as in
    run_cascade()/compare_staging()'s cascade-vapor-compression framing --
    see this module's own docstring). It is the WRONG model for a single
    physical multi-layer regenerator bed with one set of magnets, one
    pump, and one motor driving all layers at once (e.g. the real
    Astronautics/DTU_MagQueen/MAGGIE devices this function's callers below
    reproduce) -- W_eddy = k_eddy*f^2*mu0H^2 and W_base = base_frac*Qc do
    not multiply by the number of MATERIAL LAYERS inside one shared
    housing, and `fluid_mdot` is already passed identically (unscaled per
    stage) to every stage below, i.e. it already represents ONE flow
    rate through the whole bed, not N independent flows. Summing N
    independently-scaled copies of a loss term that is physically
    incurred ONCE overcounts it by roughly a factor of N (exactly N for
    W_base after scaling; a Qc-distribution-weighted factor for W_eddy/
    W_pump) -- this was found by inspecting validate_maggie_real_graded_bed
    (n=4, COP error -69.3%), validate_astronautics_graded_bed (n=6, COP
    error -81.1%), and validate_magqueen_graded_bed (n=10, COP error
    -92.2%): the error magnitude increases monotonically with n_stages,
    the fingerprint of an N-times-overcounted shared term rather than a
    genuine per-device modeling gap. When shared_hardware=True, W_mag is
    still summed per stage (thermodynamically legitimate -- each layer's
    own entropy change/work genuinely differs), but W_eddy/W_pump/W_base
    are each computed ONCE from the aggregate (frequency, mu0H_max,
    fluid_mdot, Qc_target) using the same loss_model every stage already
    shares (_LOSS_MODEL), matching how cycle_type is already shared
    unscaled across stages (see run_graded_cascade()'s cycle_type
    docstring paragraph -- the same "one physical field-change mechanism"
    reasoning applies to the loss terms and was not previously carried
    through to them).

    (mass_per_stage, continued: it may be a single float -- the more
    common case where a paper reports a total MCM mass but not a
    per-layer breakdown -- split evenly, or a list/tuple of per-stage
    masses matching len(materials).)

    See validate_maggie_real_graded_bed() below for the motivating use:
    reproducing DTU_Eriksen_MAGGIE_2016's actual reported 4-composition
    Gd/Gd-Y bed using the paper's own measured Curie temperatures,
    rather than a hypothetical composition-tuned family."""
    n_stages = len(materials)
    if n_stages == 0:
        raise ValueError("materials must be a non-empty list")
    if isinstance(mass_per_stage, (list, tuple)):
        if len(mass_per_stage) != n_stages:
            raise ValueError(f"mass_per_stage has {len(mass_per_stage)} entries, "
                              f"expected {n_stages} (one per material)")
        masses = list(mass_per_stage)
    else:
        masses = [mass_per_stage] * n_stages

    span_per_stage = total_span_K / n_stages
    T_local = T_cold_K
    stage1 = AMRSystem(material=materials[0], mu0H_max=mu0H_max,
                        mass_regenerator=masses[0], frequency=frequency,
                        fluid_mdot=fluid_mdot, regenerator_effectiveness=regenerator_effectiveness,
                        loss_model=_LOSS_MODEL, use_ntu_thermal_model=USE_NTU_THERMAL_MODEL,
                        cycle_type=cycle_type)
    r1 = stage1.run(T_local, span_per_stage)
    Qc_target = r1.Qc
    if Qc_target <= 0:
        return {"n_stages": n_stages, "span_per_stage_K": span_per_stage,
                "Qc_W": 0.0, "W_total_W": np.nan, "COP_cascade": 0.0, "feasible": False}

    W_total = 0.0
    W_mag_total = 0.0
    T_local = T_cold_K
    stage_info = []
    for i in range(n_stages):
        stage = AMRSystem(material=materials[i], mu0H_max=mu0H_max,
                           mass_regenerator=masses[i], frequency=frequency,
                           fluid_mdot=fluid_mdot, regenerator_effectiveness=regenerator_effectiveness,
                           loss_model=_LOSS_MODEL, use_ntu_thermal_model=USE_NTU_THERMAL_MODEL,
                           cycle_type=cycle_type)
        r_i = stage.run(T_local, span_per_stage)
        if r_i.Qc > 0:
            scale = Qc_target / r_i.Qc
            W_mag_total += r_i.W_mag * scale
            if not shared_hardware:
                W_total += (r_i.W_mag + r_i.W_parasitic) * scale
        else:
            W_mag_total = np.inf
            if not shared_hardware:
                W_total = np.inf
        stage_info.append({"stage": i + 1, "T_mid_K": round(T_local + span_per_stage / 2.0, 2),
                            "material": materials[i].name, "Tc_K": round(materials[i].Tc, 2),
                            "mass_kg": round(masses[i], 4), "own_Qc_W": round(r_i.Qc, 2)})
        T_local += span_per_stage

    if shared_hardware:
        W_parasitic_shared = _LOSS_MODEL.parasitic_power(frequency, mu0H_max, fluid_mdot, Qc_target)
        W_total = W_mag_total + W_parasitic_shared

    COP_cascade = Qc_target / W_total if W_total > 0 else 0.0
    return {"n_stages": n_stages, "span_per_stage_K": span_per_stage,
            "Qc_W": round(Qc_target, 1), "W_total_W": round(W_total, 1),
            "COP_cascade": round(COP_cascade, 3), "feasible": True, "stage_info": stage_info,
            "shared_hardware": shared_hardware}


def validate_maggie_real_graded_bed(cycle_type="brayton"):
    """REAL (not hypothetical, not composition-searched) reproduction of
    DTU_Eriksen_MAGGIE_2016, using the actual reported per-layer alloy
    compositions and measured Curie temperatures of the physical
    prototype -- unlike validate_risoe_dtu_graded_bed()/
    validate_cooltech_graded_bed() above, which explicitly test a
    HYPOTHETICAL graded redesign because those two devices' real hardware
    is not graded. This device's real hardware IS graded, and the paper
    reports exactly what it's graded with.

    Sources (both now in this repo's Papers/):
      - D. Eriksen, K. Engelbrecht, C.R.H. Bahl, R. Bjork, K.K. Nielsen,
        A.R. Insinga, N. Pryds, "Design and experimental tests of a
        rotary active magnetic regenerator prototype," Int. J.
        Refrigeration (2015) -- Fig. 3 (a labeled photo/diagram, not a
        text table) shows the regenerator's 11 compartments filled with
        four alloys, each labeled with its own measured Curie
        temperature: Gd (Tc=18C), Gd97.5Y2.5 (Tc=14C), Gd95Y5 (Tc=9.5C),
        Gd90Y10 (Tc=-1C), arranged hot-to-cold in that order. Read
        directly off the rendered figure (pdftotext cannot extract text
        embedded in a figure image) -- not inferred or estimated.
      - D. Eriksen, "Active magnetic regenerator refrigeration with
        rotary multi-bed technology" (PhD thesis, DTU, 2016), Ch. 3.2.2/
        Fig. 3.3 -- same figure, same prototype, confirms the same four
        Curie temperatures and the total MCM mass of 1.7kg. Sec. 6.5.2-
        6.5.3 additionally reports the exact operating point this row's
        Qc=81.5W/COP=3.6 figures come from: "a maximum second-law
        efficiency of 18% was obtained at a cooling load of 81.5 W,
        resulting in a temperature span of 15.5 K" at fAMR=0.61Hz, with
        "TH = 19.5C and delta-T = 15.5C" -- giving T_cold=277.15K
        (=19.5C-15.5C=4.0C), NOT the generic T_COLD_ASSUMED_K=289.0
        default validation_system.py falls back to for plain-Gd rows
        without a directly-reported cold-side temperature.

    Modeling choice: Gd-Y solid solutions are, like pure Gd itself, SECOND-
    ORDER (mean-field/Heisenberg) ferromagnets -- diluting Gd with
    non-magnetic Y lowers Tc without introducing first-order/
    magnetostructural character. This is therefore modeled with
    GADOLINIUM.with_Tc() (core/mce_material.py -- already built for
    exactly this purpose by core/inhomogeneous_broadening.py, Phase 22
    item 1: same J/g/M_molar/theta_D as pure Gd, only Tc shifted), NOT
    cascade.py's GD_FAMILY/LAFESIH_FAMILY machinery, which composition-
    tunes a FIRST-ORDER Landau model (core/first_order_mce.py) -- the
    wrong physics class for a Gd-Y solid solution. Because the real Tc's
    are given directly by the paper, no composition-target search
    (run_graded_cascade()'s brentq-per-stage machinery) is needed at all
    -- this calls run_explicit_material_cascade() with the four materials
    already built at their reported Tc's.

    mass_per_stage: the paper reports only the 1.7kg TOTAL MCM mass, not
    a per-layer breakdown across the four compositions (Fig. 3 gives
    compartment WIDTHS in mm for a physical drawing, not per-layer
    masses) -- so this splits 1.7kg evenly across the 4 compositions
    (0.425kg each), flagged here as an assumption the paper does not
    directly support, same honesty standard as every other
    mass-per-stage assumption elsewhere in this module.

    fluid_mdot is calibrated (brentq) to reproduce the reported Qc=81.5W
    exactly, then COP_cascade is compared to the reported COP=3.6 --
    mirroring validate_astronautics_graded_bed()'s own methodology.

    Calls run_explicit_material_cascade() with shared_hardware=True: the
    real MAGGIE prototype is ONE rotary assembly with ONE set of magnets
    and ONE pump driving all 4 Gd/Gd-Y compartments at once (Eriksen et
    al. 2015/2016, same source as the Curie temperatures above), not 4
    separate machines -- see run_explicit_material_cascade()'s own
    docstring for why the previous shared_hardware=False (default)
    behavior overcounted W_eddy/W_pump/W_base by roughly a factor of
    n_stages=4 here."""
    T_cold_K = 277.15   # 4.0C = 19.5C(T_hot) - 15.5K(span), thesis Sec. 6.5.3
    span_K = 15.5
    mu0H = 1.13
    freq = 0.61
    mass_total_kg = 1.7
    Qc_lit = 81.5
    cop_lit = 3.6

    # Hot-to-cold as drawn in Fig. 3; run_explicit_material_cascade wants
    # coldest-stage-first (matches its T_local increment from T_cold_K up).
    layer_Tc_C_hot_to_cold = [18.0, 14.0, 9.5, -1.0]
    materials = [GADOLINIUM.with_Tc(Tc_C + 273.15) for Tc_C in reversed(layer_Tc_C_hot_to_cold)]
    n_stages = len(materials)
    mass_per_stage = mass_total_kg / n_stages

    def qc_residual(mdot):
        r = run_explicit_material_cascade(T_cold_K, span_K, materials, mu0H_max=mu0H,
                                           mass_per_stage=mass_per_stage, frequency=freq,
                                           fluid_mdot=max(mdot, 1e-6), cycle_type=cycle_type,
                                           shared_hardware=True)
        return (r["Qc_W"] if r["feasible"] else 0.0) - Qc_lit

    try:
        mdot_cal = brentq(qc_residual, 1e-6, 5.0, xtol=1e-6)
    except ValueError:
        return {"feasible": False,
                "status": f"no calibration found (reported Qc={Qc_lit}W unreachable within "
                f"mdot in [1e-6, 5.0] kg/s for the real 4-layer Gd/Gd-Y bed at "
                f"T_cold={T_cold_K}K, mu0H={mu0H}T, f={freq}Hz)"}

    result = run_explicit_material_cascade(T_cold_K, span_K, materials, mu0H_max=mu0H,
                                            mass_per_stage=mass_per_stage, frequency=freq,
                                            fluid_mdot=mdot_cal, cycle_type=cycle_type,
                                            shared_hardware=True)
    result["mdot_calibrated_kg_s"] = round(mdot_cal, 5)
    result["Qc_lit_W"] = Qc_lit
    result["COP_lit"] = cop_lit
    if result["feasible"]:
        result["COP_error_pct"] = round(100 * (result["COP_cascade"] - cop_lit) / cop_lit, 1)
    return result


def run_maggie_span_sensitivity(cycle_type="brayton", verbose=True):
    """DTU_Eriksen_rotary_Gd_2015 (10.2K span, same physical prototype,
    earlier/higher-frequency paper) already calibrates cleanly under the
    single-Tc Gd approximation (step 2, err=-2.1%) -- so the real 4-layer
    graded-bed model built for the 15.5K MAGGIE row above should, at
    minimum, not make that already-working companion point worse when
    checked under the SAME real-materials treatment. This is a direct
    consistency check between the two rows sharing one physical device,
    not a new independent validation."""
    maggie = validate_maggie_real_graded_bed(cycle_type=cycle_type)

    T_cold_2015_K, span_2015_K = 285.75, 10.2   # thesis Sec. 5: best result
    # at Fig. 7's own operating point; ambient/ref T not separately
    # re-derived here -- reuses MAGGIE's own T_hot-span convention
    # (291.15K hot side implied) purely as a same-device consistency
    # check, not a claim this is the 2015 paper's own reported T_cold.
    layer_Tc_C_hot_to_cold = [18.0, 14.0, 9.5, -1.0]
    materials = [GADOLINIUM.with_Tc(Tc_C + 273.15) for Tc_C in reversed(layer_Tc_C_hot_to_cold)]
    n_stages = len(materials)
    mass_per_stage = 1.7 / n_stages
    mu0H, freq = 1.13, 0.75
    Qc_lit_2015, cop_lit_2015 = 102.8, 3.1

    def qc_residual(mdot):
        r = run_explicit_material_cascade(T_cold_2015_K, span_2015_K, materials, mu0H_max=mu0H,
                                           mass_per_stage=mass_per_stage, frequency=freq,
                                           fluid_mdot=max(mdot, 1e-6), cycle_type=cycle_type,
                                           shared_hardware=True)
        return (r["Qc_W"] if r["feasible"] else 0.0) - Qc_lit_2015

    try:
        mdot_cal = brentq(qc_residual, 1e-6, 5.0, xtol=1e-6)
        companion = run_explicit_material_cascade(T_cold_2015_K, span_2015_K, materials,
                                                    mu0H_max=mu0H, mass_per_stage=mass_per_stage,
                                                    frequency=freq, fluid_mdot=mdot_cal,
                                                    cycle_type=cycle_type, shared_hardware=True)
        companion["feasible"] = True
        companion["mdot_calibrated_kg_s"] = round(mdot_cal, 5)
        companion["Qc_lit_W"] = Qc_lit_2015
        companion["COP_lit"] = cop_lit_2015
        companion["COP_error_pct"] = round(
            100 * (companion["COP_cascade"] - cop_lit_2015) / cop_lit_2015, 1)
    except ValueError:
        companion = {"feasible": False, "status": "no calibration found for the "
                     "companion 10.2K/0.75Hz point under the real 4-layer bed"}

    if verbose:
        if maggie.get("feasible"):
            print(f"  MAGGIE (15.5K, 0.61Hz):  Qc={maggie['Qc_W']}W (target {maggie['Qc_lit_W']}W)  "
                  f"COP_cascade={maggie['COP_cascade']}  vs. lit COP={maggie['COP_lit']}  "
                  f"({maggie['COP_error_pct']:+.1f}% error)")
        else:
            print(f"  MAGGIE (15.5K, 0.61Hz):  {maggie.get('status')}")
        if companion.get("feasible"):
            print(f"  Companion (10.2K, 0.75Hz):  Qc={companion['Qc_W']}W "
                  f"(target {companion['Qc_lit_W']}W)  COP_cascade={companion['COP_cascade']}  "
                  f"vs. lit COP={companion['COP_lit']}  ({companion['COP_error_pct']:+.1f}% error)")
        else:
            print(f"  Companion (10.2K, 0.75Hz):  {companion.get('status')}")
    return {"maggie": maggie, "companion_2015": companion}


@dataclass
class StagedBaselineResult:
    """Same fields as amr_cycle.AMRCycleResult, plus n_stages so callers
    can tell whether a value came from the plain single-stage cycle or
    from the automatic multi-stage fallback below."""
    n_stages: int
    T_span: float
    Qc: float
    W_mag: float
    W_parasitic: float
    COP: float
    COP_electrical: float
    exergy_eff: float


def staged_baseline_result(T_cold_K, span, material=None, mu0H_max=2.0,
                            mass_regenerator=5.0, frequency=2.0,
                            fluid_cp=4186.0, fluid_mdot=0.08,
                            regenerator_effectiveness=0.85, max_stages=4):
    """Single-stage AMR result if the span is within that material's own
    no-load DeltaT_ad at this operating point; otherwise, the minimum
    number of IDENTICAL stages (up to max_stages) in series -- each
    covering span/n_stages, same mass/frequency/field/flow per stage --
    that reaches a positive Qc.

    Why this exists: AMRSystem.cooling_capacity() (amr_cycle.py) correctly
    returns Qc=0 once a single stage is asked to span more than ~2x its
    own no-load DeltaT_ad (see that function's MODEL LIMITATION docstring)
    -- that zero is honest single-stage physics, not a bug. But a
    "Magnetic (AMR) vs. vapor-compression vs. liquid cooling" comparison
    (main.py's run_baseline_sweep / plots.py's plot_amr_vs_baselines) that
    silently reports COP=0 for those spans reads as "AMR stops working
    above ~16K span," which overstates the limitation: core/cascade.py's
    own run_cascade()/compare_staging() (used for step 7 / fig19-20) show
    a 2-4 stage Gd cascade at the SAME mass/field/flow easily reaches
    17-20K spans with COP well above zero (see cascade_comparison.csv).
    This function brings that same "just add stages" fallback into the
    baseline/vs-baselines comparison so it degrades gracefully instead of
    hard-cutting to zero, while keeping every other modeling choice (no
    NTU thermal model, constant parasitic_fraction=0.15, i.e. NOT
    core/cascade.py's richer state-dependent-loss/NTU settings) identical
    to the single-stage baseline it's falling back from, so single- and
    multi-stage rows in the same table are still comparing the same cycle
    model -- only the stage count differs.

    Returns a StagedBaselineResult with n_stages=1 whenever the plain
    single-stage cycle already works (i.e. this is a no-op change for
    every span where the pre-existing single-stage number was already
    nonzero)."""
    if material is None:
        material = GADOLINIUM

    def _single_stage(T_cold, span_i, mass):
        sys_ = AMRSystem(material=material, mu0H_max=mu0H_max,
                          mass_regenerator=mass, frequency=frequency,
                          fluid_cp=fluid_cp, fluid_mdot=fluid_mdot,
                          regenerator_effectiveness=regenerator_effectiveness)
        return sys_.run(T_cold, span_i)

    r1 = _single_stage(T_cold_K, span, mass_regenerator)
    if r1.Qc > 0 or max_stages <= 1:
        return StagedBaselineResult(n_stages=1, T_span=span, Qc=r1.Qc,
                                     W_mag=r1.W_mag, W_parasitic=r1.W_parasitic,
                                     COP=r1.COP, COP_electrical=r1.COP_electrical,
                                     exergy_eff=r1.exergy_eff)

    for n_stages in range(2, max_stages + 1):
        span_per_stage = span / n_stages
        r_probe = _single_stage(T_cold_K, span_per_stage, mass_regenerator)
        if r_probe.Qc <= 0:
            continue
        Qc_target = r_probe.Qc
        T_local = T_cold_K
        W_mag_total = 0.0
        W_par_total = 0.0
        feasible = True
        for _ in range(n_stages):
            r_i = _single_stage(T_local, span_per_stage, mass_regenerator)
            if r_i.Qc <= 0:
                feasible = False
                break
            scale = Qc_target / r_i.Qc
            W_mag_total += r_i.W_mag * scale
            W_par_total += r_i.W_parasitic * scale
            T_local += span_per_stage
        if not feasible:
            continue
        W_total = W_mag_total + W_par_total
        COP = Qc_target / W_mag_total if W_mag_total > 0 else 0.0
        COP_electrical = Qc_target / W_total if W_total > 0 else 0.0
        T_hot = T_cold_K + span
        COP_carnot = T_cold_K / (T_hot - T_cold_K) if T_hot > T_cold_K else np.inf
        exergy_eff = COP / COP_carnot if np.isfinite(COP_carnot) and COP_carnot > 0 else 0.0
        return StagedBaselineResult(n_stages=n_stages, T_span=span, Qc=Qc_target,
                                     W_mag=W_mag_total, W_parasitic=W_par_total,
                                     COP=COP, COP_electrical=COP_electrical,
                                     exergy_eff=exergy_eff)

    # Not feasible even at max_stages: report the honest zero (same as
    # single-stage), tagged with the stage count that was actually tried.
    return StagedBaselineResult(n_stages=max_stages, T_span=span, Qc=0.0,
                                 W_mag=0.0, W_parasitic=0.0, COP=0.0,
                                 COP_electrical=0.0, exergy_eff=0.0)


def _peak_temperature(material, mu0H_max, T_range=(200.0, 340.0), n=1401):
    """Finds where THIS material's own DeltaT_ad is maximized (same
    approach as giant_mce_analysis.py's find_peak_temperature()). Needed
    because the Landau model's peak does not sit exactly at the nominal
    Tc, and -- discovered while building this cascade -- that offset is
    NOT simply translation-invariant as Tc is shifted (Debye C_lattice(T)
    depends on absolute T, not on T-Tc), so a single global offset applied
    to every stage is not accurate enough; see _target_composition_for_peak
    below for why this matters quantitatively."""
    mu0 = 4 * np.pi * 1e-7
    H = mu0H_max / mu0
    # Coarse pass over the full range, then a fine pass zoomed into the
    # coarse peak's neighborhood -- gets the same resolution near the peak
    # as a single n=1401 pass at a fraction of the evaluations (this
    # function is called inside a root-finder in
    # _target_composition_for_peak, so its cost multiplies fast).
    n_coarse = max(51, n // 10)
    Ts_coarse = np.linspace(*T_range, n_coarse)
    dT_coarse = material.delta_T_adiabatic(Ts_coarse, H)
    i0 = int(np.argmax(dT_coarse))
    lo = Ts_coarse[max(0, i0 - 2)]
    hi = Ts_coarse[min(n_coarse - 1, i0 + 2)]
    if hi <= lo:
        return float(Ts_coarse[i0])
    Ts_fine = np.linspace(lo, hi, max(21, n // 20))
    dT_fine = material.delta_T_adiabatic(Ts_fine, H)
    return float(Ts_fine[int(np.argmax(dT_fine))])


def _target_composition_for_peak(T_target_K, mu0H_max, family, max_iter=6, tol_K=0.02):
    """Solves for the composition Tc whose OWN peak DeltaT_ad lands at
    T_target_K, for the given GradedFamily. peak_T(Tc) is monotonic
    increasing in Tc (verified numerically for both GD_FAMILY and
    LAFESIH_FAMILY), so this is a straightforward bracketed root-find
    (scipy.optimize.brentq) on peak_T(Tc) - T_target_K = 0.

    This turned out to matter more than expected. The original (Phase 7)
    implementation used a simple fixed-point update (Tc += err from a
    single global offset), which worked fine for GD_FAMILY but was found
    (Phase 9, while adding LAFESIH_FAMILY) to visibly FAIL for it: this
    Landau model's transition is narrower for La(Fe,Si)13Hy than for
    Gd5Si2Ge2 (DeltaT_ad can fall by more than an order of magnitude
    within ~0.05K of the true peak -- narrower still than the ~0.2-0.5K
    scale already flagged for Gd5Si2Ge2's own transition, itself already
    narrower than the real, hysteresis/inhomogeneity-broadened transition
    Giguere et al.'s Fig. 3 shows spanning ~10-15K). At that sharpness the
    fixed-point update overshot and oscillated between iterations instead
    of converging (verified: 6-20 fixed-point iterations at tol_K=0.001-
    0.02 landed on Tc values up to ~0.15K apart, several of which left a
    stage's own dTad_noload at 0.6K instead of ~21K -- collapsing that
    stage's Qc to zero via cooling_capacity()'s span_fraction clamp, NOT
    because the material itself is bad). That narrowness is itself a
    genuine physical/numerical limitation of this idealized 6th-order
    Landau fit (flagged here, not smoothed over) -- what changed here is
    only the ROOT-FINDING METHOD used to hit each stage's own true peak as
    precisely as that narrow peak demands, so the graded-bed MECHANISM can
    be evaluated on its own terms without this numerical artifact adding
    noise on top.
    """
    ref = family.reference_material
    # Floor at 100K, well clear of the low-T region where this Landau
    # model's DeltaT_ad numerically diverges (lattice heat capacity ->0 as
    # T->0K, a pre-existing model artifact -- see e.g. T=1K giving
    # DeltaT_ad~30K for GD5SI2GE2_FIRST_ORDER, checked while debugging this
    # search). GD_FAMILY's tc_min=20K would otherwise put the search range
    # at (-20, 330)K and let brentq wander into that spurious peak.
    search_range = (max(100.0, family.tc_min - 40.0), family.tc_max + 40.0)
    offset_guess = _peak_temperature(ref, mu0H_max, T_range=search_range) - ref.Tc
    Tc_guess = T_target_K - offset_guess

    def f(Tc):
        mat = family.tuned_fn(Tc)
        return _peak_temperature(mat, mu0H_max, T_range=search_range) - T_target_K

    lo = max(family.tc_min, Tc_guess - 15.0)
    hi = min(family.tc_max, Tc_guess + 15.0)
    if lo >= hi:
        return Tc_guess  # nothing sane to bracket; let the caller's range check handle it
    try:
        if f(lo) * f(hi) > 0:
            # local bracket doesn't straddle a root -- widen to the family's
            # full documented range once before giving up
            lo, hi = family.tc_min, family.tc_max
            if f(lo) * f(hi) > 0:
                return Tc_guess  # give up gracefully, let range/feasibility checks handle it
        return brentq(f, lo, hi, xtol=min(tol_K, 0.005), maxiter=max(max_iter, 50))
    except ValueError:
        return Tc_guess


def run_graded_cascade(T_cold_K, total_span_K, n_stages, mu0H_max=2.0,
                        mass_per_stage=2.0, frequency=1.0, fluid_mdot=0.08,
                        regenerator_effectiveness=0.85,
                        apply_giguere_correction=True, family=None,
                        executor=None, cycle_type="brayton",
                        particle_diameter=None, blow_fraction=0.5,
                        pump_motor_efficiency=1.0, shared_hardware=False):
    """Curie-graded cascade (ROADMAP.md Phase 7 open item; generalized in
    Phase 9): rather than identical stages of one material (run_cascade
    above), each stage is assigned a hypothetical composition-tuned material
    from `family` (a GradedFamily -- see GD_FAMILY/LAFESIH_FAMILY above)
    whose Curie temperature matches THAT stage's own local midpoint
    temperature, on the Curie-matching principle confirmed in
    giant_mce_analysis.py (a first-order giant-MCE material performs
    strongly at its own Tc and collapses to ~zero capacity away from it).

    family defaults to GD_FAMILY (Gd5(SixGe1-x)4(-Ga), with the Giguere et
    al. (1999) empirical DeltaT_ad correction applied iff
    apply_giguere_correction=True), reproducing this function's original
    Phase 7 behavior exactly -- apply_giguere_correction is IGNORED if you
    pass an explicit `family` (build the correction into family.tuned_fn
    instead, as GD_FAMILY itself does).

    Each stage's target Tc is checked against family's documented
    tunability window (family.tc_min/tc_max). If a stage's target Tc falls
    outside that window, this function does NOT silently extrapolate a
    fictitious material -- it falls back to family.fallback_material (Gd by
    default) for that stage and records the fallback, so the returned
    result honestly reflects what is and is not supported by the
    composition-tunability literature at the requested operating point.

    Phase 16: `executor`, if given a live concurrent.futures.ProcessPoolExecutor,
    fans the n_stages independent Curie-target root-searches below out
    across that pool instead of running them one at a time in this process
    -- each stage's target composition depends only on ITS OWN local T_mid,
    not on any other stage, so parallelizing this loop changes nothing
    about the result, only how long it takes to compute (this is the
    dominant cost inside this function; see _target_composition_for_peak's
    docstring). Falls back to the original sequential loop whenever
    executor is None (the default -- no behavior change for existing
    callers) or when `family` isn't one of this module's own named
    families (a custom/closure-based family can't be resolved inside a
    separate worker process; see _family_name() above).

    `cycle_type`: ROADMAP.md Phase 17 addition, threaded through here as
    a follow-up closing that phase's own "did NOT do" item (cycle_type
    was NOT threaded through core/cascade.py's multi-stage/graded-bed
    helpers). Default "brayton" reproduces every pre-existing call's
    behavior exactly (see core.amr_cycle.AMRSystem's own
    CYCLE_TYPE_FACTORS honesty flag for what "ericsson"/"carnot" do and
    do not claim). Passed unchanged to every per-stage AMRSystem below --
    every stage in a graded cascade shares one physical field-change
    mechanism, so there is no physical reason for cycle_type to vary
    stage-to-stage.

    Phase 29 addition: `particle_diameter` (m, default None),
    `blow_fraction` (default 0.5), and `pump_motor_efficiency` (default
    1.0) are threaded through to every per-stage AMRSystem unchanged from
    their own individual defaults -- i.e. omitting all three reproduces
    every pre-Phase-29 call's behavior exactly (same "opt-in, default
    preserves old behavior" discipline as Phase 15/28's own additions to
    AMRSystem itself). This exists so `core.optimize.LayeredAMRDesignProblem`
    can expose the SAME geometry/blow-fraction/pump-efficiency design
    dimensions to a multi-layer graded-bed search that
    `AMRDesignProblem` already exposes to a single-material search --
    without this, a layered-bed NSGA-III search would be silently
    restricted to a strict subset of the single-stage search's own design
    space, which would bias any single-vs-layered comparison built on top
    of both.

    `shared_hardware` (default False, no change to any existing caller's
    numbers): see run_explicit_material_cascade()'s docstring for the
    full reasoning -- same fix, same justification (one shared magnet/
    pump/motor system driving all graded layers in a single physical
    bed), applied here because validate_astronautics_graded_bed() and
    validate_magqueen_graded_bed() call THIS function, not
    run_explicit_material_cascade(). compare_graded_cascade() below
    still defaults to shared_hardware=False, so graded_cascade_comparison.csv
    and design_recommendations.txt's lever-3 numbers are unchanged by
    this addition -- they describe a hypothetical from-scratch design
    exploration where whether the graded bed is one shared-hardware unit
    or N separate modules is a design choice, not yet a fact about a
    built device, and re-deriving that whole comparison table was judged
    out of scope for this pass; see results/regenerative_amplification_diagnostic.txt
    for a worked example of the resulting magnitude on one operating point."""
    if family is None:
        family = GD_FAMILY if apply_giguere_correction else GradedFamily(
            name=GD_FAMILY.name,
            tuned_fn=lambda Tc: composition_tuned_material(Tc, apply_giguere_correction=False),
            tc_min=GD_FAMILY.tc_min, tc_max=GD_FAMILY.tc_max,
            reference_material=GD_FAMILY.reference_material,
            fallback_material=GD_FAMILY.fallback_material)

    span_per_stage = total_span_K / n_stages

    # Compute each stage's own T_mid the same way (incremental addition, so
    # this is bit-for-bit identical to the pre-Phase-16 code) whether or not
    # the target-composition search below ends up running in parallel.
    mid_temps = []
    T_local = T_cold_K
    for _i in range(n_stages):
        mid_temps.append(T_local + span_per_stage / 2.0)
        T_local += span_per_stage

    family_name = _family_name(family)
    Tc_targets = None
    if executor is not None and family_name is not None and n_stages > 1:
        Tc_targets = _pool_map_or_none(
            executor, _stage_target_worker,
            [(t, mu0H_max, family_name) for t in mid_temps])
    if Tc_targets is None:
        # Phase 31: also the path taken on any pool timeout/failure above,
        # not just when executor is None/unsuitable -- bit-for-bit
        # identical to the pre-Phase-31 sequential result either way.
        Tc_targets = [_target_composition_for_peak(t, mu0H_max, family) for t in mid_temps]

    stage_materials = []
    stage_info = []
    for i, (T_mid_stage, Tc_target) in enumerate(zip(mid_temps, Tc_targets)):
        if family.tc_min <= Tc_target <= family.tc_max:
            mat = family.tuned_fn(Tc_target)
            stage_info.append({"stage": i + 1, "T_mid_K": round(T_mid_stage, 1),
                                "Tc_target_K": round(Tc_target, 1),
                                "material": mat.name, "in_range": True})
        else:
            mat = family.fallback_material
            stage_info.append({"stage": i + 1, "T_mid_K": round(T_mid_stage, 1),
                                "Tc_target_K": round(Tc_target, 1),
                                "material": f"{mat.name} (fallback -- Tc target outside "
                                            f"{family.tc_min:.0f}-{family.tc_max:.0f}K "
                                            f"documented {family.name} range)",
                                "in_range": False})
        stage_materials.append(mat)

    n_fallback = sum(1 for s in stage_info if not s["in_range"])

    T_local = T_cold_K
    stage1 = AMRSystem(material=stage_materials[0], mu0H_max=mu0H_max,
                        mass_regenerator=mass_per_stage, frequency=frequency,
                        fluid_mdot=fluid_mdot, regenerator_effectiveness=regenerator_effectiveness,
                        loss_model=_LOSS_MODEL, use_ntu_thermal_model=USE_NTU_THERMAL_MODEL,
                        cycle_type=cycle_type, particle_diameter=particle_diameter,
                        blow_fraction=blow_fraction, pump_motor_efficiency=pump_motor_efficiency)
    r1 = stage1.run(T_local, span_per_stage)
    Qc_target = r1.Qc
    if Qc_target <= 0:
        return {"n_stages": n_stages, "span_per_stage_K": span_per_stage,
                "Qc_W": 0.0, "W_total_W": np.nan, "COP_cascade": 0.0,
                "feasible": False, "n_stages_out_of_range": n_fallback,
                "stage_info": stage_info}

    W_total = 0.0
    W_mag_total = 0.0
    T_local = T_cold_K
    for i in range(n_stages):
        stage = AMRSystem(material=stage_materials[i], mu0H_max=mu0H_max,
                           mass_regenerator=mass_per_stage, frequency=frequency,
                           fluid_mdot=fluid_mdot, regenerator_effectiveness=regenerator_effectiveness,
                           loss_model=_LOSS_MODEL, use_ntu_thermal_model=USE_NTU_THERMAL_MODEL,
                           cycle_type=cycle_type, particle_diameter=particle_diameter,
                           blow_fraction=blow_fraction, pump_motor_efficiency=pump_motor_efficiency)
        r_i = stage.run(T_local, span_per_stage)
        if r_i.Qc > 0:
            scale = Qc_target / r_i.Qc
            W_mag_total += r_i.W_mag * scale
            if not shared_hardware:
                W_total += (r_i.W_mag + r_i.W_parasitic) * scale
        else:
            W_mag_total = np.inf
            if not shared_hardware:
                W_total = np.inf
        T_local += span_per_stage

    if shared_hardware:
        W_parasitic_shared = _LOSS_MODEL.parasitic_power(frequency, mu0H_max, fluid_mdot, Qc_target)
        W_total = W_mag_total + W_parasitic_shared

    COP_cascade = Qc_target / W_total if W_total > 0 else 0.0
    return {"n_stages": n_stages, "span_per_stage_K": span_per_stage,
            "Qc_W": round(Qc_target, 1), "W_total_W": round(W_total, 1),
            "COP_cascade": round(COP_cascade, 3), "feasible": True,
            "n_stages_out_of_range": n_fallback, "stage_info": stage_info,
            "shared_hardware": shared_hardware}


def compare_graded_cascade(T_cold_C=18.0, spans=range(5, 21), stage_counts=(1, 2, 3, 4),
                            mass_per_stage=2.0, family=None,
                            out_csv="results/graded_cascade_comparison.csv",
                            parallel=True, max_workers=None):
    """Same sweep as compare_staging(), but using run_graded_cascade()
    instead of identical-stage run_cascade(). family is passed straight
    through to run_graded_cascade() (default: GD_FAMILY, i.e. the original
    Phase 7 Gd5(SixGe1-x)4(-Ga) behavior). At the ASHRAE data-center range
    (T_cold_C=18 -> T_cold_K=291.15K) with the default GD_FAMILY, each
    stage's needed composition Tc is checked against GIANT_MCE_TC_MAX_K=
    290K: for small spans/stage-counts every stage stays within that
    documented range and the cascade is fully buildable; for larger spans
    and/or more stages the hottest stage(s) push above 290K and fall back
    to plain Gd for that stage only. See the __main__ block below for the
    actual breakdown across the sweep (computed, not assumed).

    Phase 16: this span x stage_count sweep is embarrassingly parallel --
    every cell is an independent run_graded_cascade() call with no shared
    mutable state -- so by default (parallel=True) it now fans out across a
    ProcessPoolExecutor (max_workers defaults to min(#cells, cpu_count())).
    This is what previously made this function (and the 6-layer
    Astronautics reproduction that reuses it transitively) the single
    slowest stage in the full pipeline. Falls back automatically to the
    original sequential loop if: `family` is a custom object this module
    doesn't recognize by identity (can't be safely rebuilt inside a worker
    process -- see _family_name()); there's only one cell to compute; or
    process-pool creation itself fails (e.g. a sandboxed environment
    without subprocess spawn rights) -- in every case the returned rows and
    CSV are bit-for-bit identical to the sequential path, since each cell
    is computed by the exact same run_graded_cascade() call either way.
    Pass parallel=False to force the old sequential behavior explicitly."""
    T_cold_K = T_cold_C + 273.15
    cells = [(span, n) for span in spans for n in stage_counts]
    family_name = _family_name(family)

    cell_results = None
    if parallel and family_name is not None and len(cells) > 1:
        workers = max_workers or min(len(cells), os.cpu_count() or 1)
        args = [(T_cold_K, span, n, mass_per_stage, family_name) for span, n in cells]
        pool = None
        try:
            with _single_threaded_blas_env():
                pool = ProcessPoolExecutor(max_workers=workers, initializer=_pool_worker_init)
                raw = _pool_map_or_none(pool, _compare_graded_cascade_cell, args)
            if raw is not None:
                cell_results = {(span, n): res for span, n, res in raw}
        except Exception as exc:
            # Pool CREATION itself failed outright (e.g. a sandboxed
            # environment without subprocess spawn rights) -- distinct
            # from the map()-level timeout/failure _pool_map_or_none
            # already handles internally, but the same fallback applies.
            logger.warning("Phase 31: ProcessPoolExecutor creation failed "
                            "(%s: %s) -- falling back to sequential.",
                            type(exc).__name__, exc)
            cell_results = None
        finally:
            # Phase 31: never `with ProcessPoolExecutor(...) as pool:` here --
            # that context manager's own __exit__ calls the blocking
            # pool.shutdown(wait=True), which can hang exactly like the
            # bare future.result()/pool.map() calls this fix targets.
            _safe_pool_shutdown(pool)

    if cell_results is None:
        cell_results = {}
        for span, n in cells:
            cell_results[(span, n)] = run_graded_cascade(
                T_cold_K, span, n, mass_per_stage=mass_per_stage, family=family)

    rows = []
    all_stage_info = []
    for span in spans:
        row = {"span_K": span}
        for n in stage_counts:
            res = cell_results[(span, n)]
            row[f"Graded_{n}stage_COP"] = res["COP_cascade"] if res["feasible"] else None
            row[f"Graded_{n}stage_Qc_W"] = res["Qc_W"] if res["feasible"] else None
            row[f"Graded_{n}stage_n_fallback_to_Gd"] = res["n_stages_out_of_range"]
            all_stage_info.append({"span_K": span, "n_stages": n, "stage_info": res["stage_info"]})
        rows.append(row)

    with open(out_csv, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=rows[0].keys())
        writer.writeheader()
        writer.writerows(rows)
    return rows, all_stage_info


def compare_staging(T_cold_C=18.0, spans=range(5, 21), stage_counts=(1, 2, 3, 4),
                     material=None, mass_per_stage=2.0,
                     out_csv="results/cascade_comparison.csv"):
    T_cold_K = T_cold_C + 273.15
    rows = []
    for span in spans:
        T_hot_K = T_cold_K + span
        vcc = vapor_compression_cop(T_cold_K, T_hot_K)
        liq = liquid_cooling_cop(T_cold_K, T_hot_K)
        row = {"span_K": span, "VaporCompression_COP": round(vcc.COP, 2),
               "LiquidCooling_COP": round(liq.COP, 2)}
        for n in stage_counts:
            res = run_cascade(T_cold_K, span, n, material=material, mass_per_stage=mass_per_stage)
            row[f"AMR_{n}stage_COP"] = res["COP_cascade"] if res["feasible"] else None
            row[f"AMR_{n}stage_Qc_W"] = res["Qc_W"] if res["feasible"] else None
        rows.append(row)

    with open(out_csv, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=rows[0].keys())
        writer.writeheader()
        writer.writerows(rows)
    return rows


def validate_astronautics_graded_bed(apply_correction=None, cycle_type="brayton"):
    """Phase 9 follow-up: builds a 6-layer Curie-graded La(Fe,Si)13Hy bed
    (LAFESIH_FAMILY) at the REAL Astronautics_rotary_2014 operating point
    (Jacobs et al., Int. J. Refrig. 37 (2014) 84-91: mu0H=1.44T, 1.52kg
    total MCM, f=4Hz, T_cold=305K/32C, span=11K, reported Qc=2502W,
    COP=1.9) and runs it through the SAME calibrate-then-validate
    methodology core.validation_system.py uses for every other device:
    fluid_mdot is calibrated (brentq) to reproduce the reported Qc exactly,
    then the resulting COP is compared against the reported COP -- so this
    tests the model's predicted EFFICIENCY, not its ability to predict Qc
    from first principles (same caveat validation_system.py's own
    docstring states for its methodology).

    This directly tests the hypothesis raised in ROADMAP.md Phase 9: that
    validation_system.py's single-Tc=287K LAFESIH_FIRST_ORDER material
    failed to calibrate against this device because the real device is SIX
    Curie-graded layers (~303.6-316.2K), not because the general modeling
    approach is wrong. Layer Tc targets are spread evenly across that
    reported range (linspace(303.6, 316.2, 6)) as the best available
    approximation -- Jacobs et al. (2014) does not tabulate the individual
    layer compositions/Tc values themselves, only the range. mass_per_stage
    = 1.52/6 kg (equal split; per-layer masses are also not individually
    reported).

    apply_correction defaults to LAFESIH_FIRST_ORDER's own dTad_correction
    (1.0, i.e. uncorrected -- see that material's honesty flags for why no
    Giguere-style correction is available for this family) unless
    explicitly overridden.

    Phase 16: brentq calls run_graded_cascade() repeatedly to calibrate
    mdot, and each of those calls does 6 independent per-stage Curie-target
    searches (run_graded_cascade's dominant cost -- see that function's own
    Phase 16 note) -- so this function opens ONE ProcessPoolExecutor up
    front and reuses it across every brentq iteration plus the final
    result call, rather than paying process-startup cost repeatedly. This
    only activates for the default (apply_correction=None) LAFESIH_FAMILY
    case: the apply_correction override below builds a closure-based family
    that can't be sent to a worker process, so that path automatically
    falls back to run_graded_cascade's own sequential loop instead (see
    _family_name()) -- same numeric result either way, just slower.

    Returns a dict with the calibrated mdot, predicted Qc (should equal
    Qc_lit_W by construction), predicted vs. literature COP and its error,
    and the per-stage breakdown -- or a "no calibration found" status dict
    if no mdot in [1e-6, 1.0] kg/s reproduces the reported Qc.

    `cycle_type` (ROADMAP.md Phase 17 addition, follow-up closing that
    phase's "did NOT do" item on cascade.py): default "brayton"
    reproduces pre-existing behavior exactly. Astronautics_rotary_2014 is
    itself the one device `core.validation_system.infer_cycle_type_for_device()`
    flags as "rotary" (continuous-field) via its naming-convention proxy
    -- and it's ALSO the device with the largest single-Tc-approximation
    COP error already on record (single-Tc: -81.1%, see run 8/step 2's
    docstring). `run_astronautics_cycle_type_sensitivity()` below
    directly checks whether "ericsson" narrows that error for the
    graded-bed reproduction the way it narrowed DTU_Eriksen_rotary_Gd_2015's
    (-2.1% -> +0.6%, Phase 17's own single comparable-device result).
    """
    from scipy.optimize import brentq

    T_cold_K = 305.0
    span_K = 11.0
    mu0H = 1.44
    mass_total = 1.52
    n_stages = 6
    freq = 4.0
    Qc_lit = 2502.0
    cop_lit = 1.9

    family = LAFESIH_FAMILY
    if apply_correction is not None:
        base = LAFESIH_FAMILY

        def tuned_fn(Tc, _corr=float(apply_correction)):
            mat = lafesih_composition_tuned_material(Tc)
            mat.dTad_correction = _corr
            return mat
        family = GradedFamily(name=base.name, tuned_fn=tuned_fn, tc_min=base.tc_min,
                               tc_max=base.tc_max, reference_material=base.reference_material,
                               fallback_material=base.fallback_material)

    pool = None
    _blas_env_cm = None
    if _family_name(family) is not None:
        try:
            _blas_env_cm = _single_threaded_blas_env()
            _blas_env_cm.__enter__()
            pool = ProcessPoolExecutor(max_workers=min(n_stages, os.cpu_count() or 1),
                                        initializer=_pool_worker_init)
        except Exception:
            pool = None
            if _blas_env_cm is not None:
                _blas_env_cm.__exit__(None, None, None)
                _blas_env_cm = None

    try:
        def qc_residual(mdot):
            r = run_graded_cascade(T_cold_K, span_K, n_stages, mu0H_max=mu0H,
                                    mass_per_stage=mass_total / n_stages, frequency=freq,
                                    fluid_mdot=max(mdot, 1e-6), family=family, executor=pool,
                                    cycle_type=cycle_type,
                                    shared_hardware=True)
            return (r["Qc_W"] if r["feasible"] else 0.0) - Qc_lit

        try:
            mdot_cal = brentq(qc_residual, 1e-6, 1.0, xtol=1e-6)
        except ValueError:
            return {"feasible": False, "status": "no calibration found "
                    "(reported Qc unreachable within mdot in [1e-6, 1.0] kg/s "
                    "for the 6-layer graded La(Fe,Si)13Hy bed)"}

        result = run_graded_cascade(T_cold_K, span_K, n_stages, mu0H_max=mu0H,
                                     mass_per_stage=mass_total / n_stages, frequency=freq,
                                     fluid_mdot=mdot_cal, family=family, executor=pool,
                                     cycle_type=cycle_type, shared_hardware=True)

    finally:
        # Phase 31: bounded, non-hanging shutdown -- see _safe_pool_shutdown()'s
        # own docstring for why plain pool.shutdown(wait=True) is unsafe here.
        _safe_pool_shutdown(pool)
        if _blas_env_cm is not None:
            _blas_env_cm.__exit__(None, None, None)

    result["mdot_calibrated_kg_s"] = round(mdot_cal, 5)
    result["Qc_lit_W"] = Qc_lit
    result["COP_lit"] = cop_lit
    result["COP_error_pct"] = round(100 * (result["COP_cascade"] - cop_lit) / cop_lit, 1)
    return result


def validate_magqueen_graded_bed(mass_total_kg=1.0, n_stages=10,
                                  apply_correction=None, cycle_type="brayton",
                                  use_internal_pool=True):
    """Extends the validate_astronautics_graded_bed() pattern (ROADMAP.md
    Phase 9) to DTU_MagQueen_2018, the calibration_failure_diagnostics.txt
    step-2c diagnostic's OTHER LAFESIH_FAMILY device (margin=-24.97K,
    "structural"). Like Astronautics_rotary_2014, MagQueen's own reported
    hardware genuinely IS Curie-graded -- data/amr_experimental_benchmarks.csv's
    own row note: "13 packed beds x 10 layers La(Fe,Mn,Si)13Hz alloy each"
    -- so n_stages=10 (one bed's own axial layer count) is a real,
    paper-derived structural choice here, unlike the arbitrary n_stages
    used below for Risoe_DTU/Cooltech (which have no reported layers to
    match). T_cold_K=305.0 reuses validation_system.T_COLD_LAFESIH_K, the
    same La(Fe,Si)13Hy-specific assumption already applied to this row by
    the single-Tc validation in step 2, for an apples-to-apples comparison
    -- NOT independently confirmed for MagQueen's own reported operating
    point (Johra et al. 2019 / Dall'Olio et al. 2018 are not in this
    repo's Papers/; see the CSV row's own note).

    mass_total_kg is NOT reported for this device (no volume/density given
    -- see the CSV row's own note, same gap flagged for Cooltech_2013_rotary)
    and is left as a required argument rather than a silently-assumed
    default; run_magqueen_mass_sensitivity() below sweeps it instead of
    picking one value, per the same reasoning validation_system.py already
    applies (mass=1.0kg fallback, flagged as illustrative) to this row's
    single-Tc check.

    Qc_lit=1200.0W and cop_lit=4.0 are the CSV row's own DERIVED (not
    directly reported) Qc/COP_cooling values -- see that row's note for
    the Qh=Qc+W identity used to get them from the paper's own reported
    heating-mode numbers.

    use_internal_pool=True (default) opens its own ProcessPoolExecutor to
    parallelize the n_stages=10 per-stage Curie-target searches, same as
    validate_astronautics_graded_bed(). Set False when this function is
    ITSELF being called from inside a worker process (e.g. by
    run_magqueen_mass_sensitivity()'s outer, mass-level pool below) --
    ProcessPoolExecutor workers are daemonic and cannot spawn their own
    child processes, so nesting pools would raise
    "daemonic processes are not allowed to have children" rather than
    silently work. With use_internal_pool=False this falls back to
    run_graded_cascade's own sequential per-stage loop -- same numeric
    result either way, just computed without a nested pool.
    """
    from scipy.optimize import brentq

    T_cold_K = 305.0
    span_K = 25.0
    mu0H = 1.6
    freq = 2.0
    Qc_lit = 1200.0
    cop_lit = 4.0

    family = LAFESIH_FAMILY
    if apply_correction is not None:
        base = LAFESIH_FAMILY

        def tuned_fn(Tc, _corr=float(apply_correction)):
            mat = lafesih_composition_tuned_material(Tc)
            mat.dTad_correction = _corr
            return mat
        family = GradedFamily(name=base.name, tuned_fn=tuned_fn, tc_min=base.tc_min,
                               tc_max=base.tc_max, reference_material=base.reference_material,
                               fallback_material=base.fallback_material)

    pool = None
    _blas_env_cm = None
    if use_internal_pool and _family_name(family) is not None:
        try:
            _blas_env_cm = _single_threaded_blas_env()
            _blas_env_cm.__enter__()
            pool = ProcessPoolExecutor(max_workers=min(n_stages, os.cpu_count() or 1),
                                        initializer=_pool_worker_init)
        except Exception:
            pool = None
            if _blas_env_cm is not None:
                _blas_env_cm.__exit__(None, None, None)
                _blas_env_cm = None

    try:
        def qc_residual(mdot):
            r = run_graded_cascade(T_cold_K, span_K, n_stages, mu0H_max=mu0H,
                                    mass_per_stage=mass_total_kg / n_stages, frequency=freq,
                                    fluid_mdot=max(mdot, 1e-6), family=family, executor=pool,
                                    cycle_type=cycle_type,
                                    shared_hardware=True)
            return (r["Qc_W"] if r["feasible"] else 0.0) - Qc_lit

        try:
            mdot_cal = brentq(qc_residual, 1e-6, 1.0, xtol=1e-6)
        except ValueError:
            return {"feasible": False, "mass_total_kg": mass_total_kg,
                    "status": "no calibration found (reported Qc unreachable within "
                    f"mdot in [1e-6, 1.0] kg/s for the {n_stages}-layer graded "
                    f"La(Fe,Si)13Hy bed at mass_total_kg={mass_total_kg})"}

        result = run_graded_cascade(T_cold_K, span_K, n_stages, mu0H_max=mu0H,
                                     mass_per_stage=mass_total_kg / n_stages, frequency=freq,
                                     fluid_mdot=mdot_cal, family=family, executor=pool,
                                     cycle_type=cycle_type, shared_hardware=True)

    finally:
        # Phase 31: bounded, non-hanging shutdown -- see _safe_pool_shutdown()'s
        # own docstring for why plain pool.shutdown(wait=True) is unsafe here.
        _safe_pool_shutdown(pool)
        if _blas_env_cm is not None:
            _blas_env_cm.__exit__(None, None, None)

    result["mdot_calibrated_kg_s"] = round(mdot_cal, 5)
    result["mass_total_kg"] = mass_total_kg
    result["Qc_lit_W"] = Qc_lit
    result["COP_lit"] = cop_lit
    result["COP_error_pct"] = round(100 * (result["COP_cascade"] - cop_lit) / cop_lit, 1)
    return result


def _magqueen_mass_worker(mass_total_kg):
    """Top-level (picklable) worker for run_magqueen_mass_sensitivity()'s
    outer, mass-level pool. use_internal_pool=False because this already
    runs inside a worker process -- see validate_magqueen_graded_bed()'s
    own docstring for why nesting pools would fail outright rather than
    just be slow."""
    return validate_magqueen_graded_bed(mass_total_kg=mass_total_kg, use_internal_pool=False)


def run_magqueen_mass_sensitivity(masses_kg=(0.5, 1.0, 2.0, 5.0, 10.0), verbose=True,
                                   parallel=True):
    """DTU_MagQueen_2018's mass_MCM_kg is not reported (see
    validate_magqueen_graded_bed()'s own docstring), so rather than
    silently picking one placeholder value, this sweeps a small range of
    plausible total-MCM masses (0.5-10 kg -- bracketing Astronautics'
    1.52kg and this repo's other packed-bed devices, e.g.
    DTU_Eriksen_rotary_Gd_2015's 1.7kg, without asserting a specific
    figure this repo's corpus does not contain) and reports how the
    10-layer graded-bed calibration/COP-error result moves with mass, so
    any conclusion about whether this device's structural gap closes is
    read as a sensitivity range rather than a single, unsupported number.

    parallel=True (default) runs the swept masses THEMSELVES in parallel
    across one shared ProcessPoolExecutor (max_workers=min(len(masses_kg),
    cpu_count())), rather than each mass sequentially opening and closing
    its own internal 10-worker pool for its per-stage searches (the
    original behavior, still used when parallel=False or if pool startup
    fails for any reason). Each worker computes one mass's full
    calibration single-threaded (use_internal_pool=False in
    _magqueen_mass_worker) -- ProcessPoolExecutor workers can't spawn
    their own children, so parallelizing across masses and parallelizing
    across stages-within-a-mass are mutually exclusive, not additive;
    given 5 masses vs. 10 stages, parallelizing across masses uses the
    available cores about as well while avoiding 5 separate pool
    startup/teardown cycles."""
    results = None
    if parallel and len(masses_kg) > 1:
        pool = None
        _blas_env_cm = None
        try:
            _blas_env_cm = _single_threaded_blas_env()
            _blas_env_cm.__enter__()
            pool = ProcessPoolExecutor(max_workers=min(len(masses_kg), os.cpu_count() or 1),
                                        initializer=_pool_worker_init)
            # Phase 31 fix: this exact line -- bare `f.result()` with no
            # timeout on a submit-then-collect pattern -- is the reproduced
            # cause of test_magqueen_mass_sensitivity_parallel_matches_sequential
            # hanging indefinitely in a sandboxed environment. See this
            # module's "Phase 31" comment block above _pool_submit_all_or_none
            # for the full diagnosis.
            results = _pool_submit_all_or_none(pool, _magqueen_mass_worker, masses_kg)
        except Exception as exc:
            logger.warning("Phase 31: ProcessPoolExecutor creation failed "
                            "(%s: %s) -- falling back to sequential.",
                            type(exc).__name__, exc)
            results = None  # fall through to the sequential path below
        finally:
            _safe_pool_shutdown(pool)
            if _blas_env_cm is not None:
                _blas_env_cm.__exit__(None, None, None)

    if results is None:
        results = [validate_magqueen_graded_bed(mass_total_kg=m) for m in masses_kg]

    if verbose:
        for m, r in zip(masses_kg, results):
            if r.get("feasible"):
                print(f"  mass_total={m:5.2f}kg  mdot_cal={r['mdot_calibrated_kg_s']:.5f}kg/s  "
                      f"Qc={r['Qc_W']:.1f}W  COP_cascade={r['COP_cascade']:.3f}  "
                      f"COP_error={r['COP_error_pct']:+.1f}%  "
                      f"n_fallback_to_Gd={r['n_stages_out_of_range']}/10")
            else:
                print(f"  mass_total={m:5.2f}kg  {r.get('status', 'infeasible')}")
    n_feasible = sum(1 for r in results if r.get("feasible"))
    if verbose:
        print(f"CONCLUSION: {n_feasible}/{len(masses_kg)} swept masses calibrate at all -- "
              "since mass_MCM_kg is genuinely unreported for this device, treat any single "
              "mass's COP_error_pct above as illustrative of the SENSITIVITY, not as this "
              "repo's calibrated answer for MagQueen (unlike Astronautics_rotary_2014, whose "
              "1.52kg total mass IS directly reported by Jacobs et al. 2014).")
    return results


def validate_risoe_dtu_graded_bed(n_stages=6, apply_correction=None,
                                   cycle_type="brayton", family=None):
    """Extends the graded-bed structural-fix pattern (ROADMAP.md Phase 9,
    validate_astronautics_graded_bed()) to Risoe_DTU_Gd_2011, one of the
    calibration_failure_diagnostics.txt step-2c devices flagged
    STRUCTURAL (margin=2*dTad_noload-span=-23.80K -- the largest-span
    COP-bearing benchmark row in this repo's set, 30K).

    IMPORTANT DIFFERENCE from Astronautics/MagQueen: Risoe_DTU_Gd_2011's
    own reported hardware (Engelbrecht et al., IRC Purdue 2010 / part
    I-II ScienceDirect 2016) is a SINGLE packed bed of plain Gd, not a
    Curie-graded multi-layer bed -- data/amr_experimental_benchmarks.csv's
    own row has no layer/grading note the way Astronautics' and MagQueen's
    rows do. There is therefore NO literature basis for claiming this
    specific device IS graded. What this function tests instead is the
    genuinely different question core/validation_system.py's
    calibration_failure_diagnostics.txt item 3 raises: could a
    HYPOTHETICAL Curie-graded redesign, at this device's same reported
    field/mass/frequency and reported Qc/COP target, close the structural
    gap a single-Tc material cannot? It reuses GD_FAMILY (the same
    Gd5(SixGe1-x)4(-Ga) composition-tunable family compare_graded_cascade()
    already uses for its own non-device-specific sweep, Giguere-corrected
    by default) rather than LAFESIH_FAMILY, since GD_FAMILY's documented
    tc window (20-290K) brackets plain Gd's own ~294K Tc -- the material
    class this device's own reported hardware actually uses -- unlike
    LAFESIH_FAMILY's 190-340K window, centered on an unrelated,
    much-hotter material class this device was never reported to use.

    n_stages=6 is an arbitrary, non-paper-derived choice mirroring
    validate_astronautics_graded_bed()'s own 6-layer count, for a
    like-for-like structural comparison -- Engelbrecht et al. report no
    per-layer breakdown to match against, because the real device has no
    layers to match. T_cold_K=289.0 reuses
    validation_system.T_COLD_ASSUMED_K, the same Gd-centered default
    already applied to this row by the single-Tc validation in step 2."""
    from scipy.optimize import brentq

    T_cold_K = 289.0
    span_K = 30.0
    mu0H = 1.1
    mass_total = 0.1955
    freq = 1.0
    Qc_lit = 35.0
    cop_lit = 5.0

    if family is None:
        family = GD_FAMILY
        if apply_correction is not None:
            base = GD_FAMILY

            def tuned_fn(Tc, _corr=float(apply_correction)):
                mat = composition_tuned_material(Tc, apply_giguere_correction=False)
                mat.dTad_correction = _corr
                return mat
            family = GradedFamily(name=base.name, tuned_fn=tuned_fn, tc_min=base.tc_min,
                                   tc_max=base.tc_max, reference_material=base.reference_material,
                                   fallback_material=base.fallback_material)

    pool = None
    _blas_env_cm = None
    if _family_name(family) is not None:
        try:
            _blas_env_cm = _single_threaded_blas_env()
            _blas_env_cm.__enter__()
            pool = ProcessPoolExecutor(max_workers=min(n_stages, os.cpu_count() or 1),
                                        initializer=_pool_worker_init)
        except Exception:
            pool = None
            if _blas_env_cm is not None:
                _blas_env_cm.__exit__(None, None, None)
                _blas_env_cm = None

    try:
        def qc_residual(mdot):
            r = run_graded_cascade(T_cold_K, span_K, n_stages, mu0H_max=mu0H,
                                    mass_per_stage=mass_total / n_stages, frequency=freq,
                                    fluid_mdot=max(mdot, 1e-6), family=family, executor=pool,
                                    cycle_type=cycle_type)
            return (r["Qc_W"] if r["feasible"] else 0.0) - Qc_lit

        try:
            mdot_cal = brentq(qc_residual, 1e-6, 5.0, xtol=1e-6)
        except ValueError:
            return {"feasible": False, "status": "no calibration found "
                    "(reported Qc unreachable within mdot in [1e-6, 5.0] kg/s "
                    f"for the {n_stages}-stage hypothetical graded Gd-alloy bed)"}

        result = run_graded_cascade(T_cold_K, span_K, n_stages, mu0H_max=mu0H,
                                     mass_per_stage=mass_total / n_stages, frequency=freq,
                                     fluid_mdot=mdot_cal, family=family, executor=pool,
                                     cycle_type=cycle_type)
    finally:
        # Phase 31: bounded, non-hanging shutdown -- see _safe_pool_shutdown()'s
        # own docstring for why plain pool.shutdown(wait=True) is unsafe here.
        _safe_pool_shutdown(pool)
        if _blas_env_cm is not None:
            _blas_env_cm.__exit__(None, None, None)

    result["mdot_calibrated_kg_s"] = round(mdot_cal, 5)
    result["Qc_lit_W"] = Qc_lit
    result["COP_lit"] = cop_lit
    result["COP_error_pct"] = round(100 * (result["COP_cascade"] - cop_lit) / cop_lit, 1)
    return result


def validate_cooltech_graded_bed(n_stages=6, mass_total_kg=1.0,
                                  apply_correction=None, cycle_type="brayton", family=None):
    """Extends the graded-bed structural-fix question to
    Cooltech_2013_rotary, calibration_failure_diagnostics.txt's other
    STRUCTURAL row and this repo's largest-span benchmark device overall
    (42K, margin=-38.99K). Cooltech_2013_rotary is a CAPACITY-ONLY row
    (no COP reported -- Greco et al. 2019's secondary-source review gives
    span/Qc only), so this mirrors
    core.validation_system.run_capacity_only_calibration_check()'s own
    feasibility-only treatment (does the reported Qc become reachable at
    all, for SOME mdot) rather than validate_risoe_dtu_graded_bed()'s
    full Qc-then-COP methodology -- there is no COP_lit here to compare
    against.

    Same caveat as validate_risoe_dtu_graded_bed(): Cooltech's own
    reported hardware (Greco et al. 2019, citing Kitanovski et al. 2015c/
    Chaudron et al. 2018/Lionte et al. 2018b -- none in this repo's
    Papers/) is described only as "packed bed regenerators (Gd or
    Gd-Tb)", with no grading/layering reported -- this tests the same
    HYPOTHETICAL graded-redesign question, using GD_FAMILY for the same
    material-class reason given there. n_stages=6 is the same arbitrary,
    non-paper-derived choice. mass_total_kg is NOT reported for this
    device either (data/amr_experimental_benchmarks.csv's own row note:
    "no volume/mass column ... left blank, so core/validation_system.py's
    calibrate_and_check() falls back to its own existing mass=1.0kg
    default") -- run_cooltech_mass_sensitivity() below sweeps it rather
    than treating any one value as calibrated. T_cold_K=289.0 reuses
    validation_system.T_COLD_ASSUMED_K, same Gd-centered default already
    applied to this row's own single-Tc feasibility check in step 2."""
    T_cold_K = 289.0
    span_K = 42.0
    mu0H = 0.98
    freq = 4.0
    Qc_lit = 120.0

    if family is None:
        family = GD_FAMILY
        if apply_correction is not None:
            base = GD_FAMILY

            def tuned_fn(Tc, _corr=float(apply_correction)):
                mat = composition_tuned_material(Tc, apply_giguere_correction=False)
                mat.dTad_correction = _corr
                return mat
            family = GradedFamily(name=base.name, tuned_fn=tuned_fn, tc_min=base.tc_min,
                                   tc_max=base.tc_max, reference_material=base.reference_material,
                                   fallback_material=base.fallback_material)

    r_probe = run_graded_cascade(T_cold_K, span_K, n_stages, mu0H_max=mu0H,
                                  mass_per_stage=mass_total_kg / n_stages, frequency=freq,
                                  fluid_mdot=1.0, family=family, cycle_type=cycle_type)
    if not r_probe["feasible"] or r_probe["Qc_W"] <= 0:
        return {"feasible": False, "mass_total_kg": mass_total_kg,
                "status": f"no calibration found (the {n_stages}-stage hypothetical graded "
                f"Gd-alloy bed's own Qc(mdot=1.0kg/s)={r_probe.get('Qc_W', 0.0)}W is already "
                "0 -- span still exceeds this bed's own achievable no-load dTad even after "
                "graded splitting)"}

    from scipy.optimize import brentq

    def qc_residual(mdot):
        r = run_graded_cascade(T_cold_K, span_K, n_stages, mu0H_max=mu0H,
                                mass_per_stage=mass_total_kg / n_stages, frequency=freq,
                                fluid_mdot=max(mdot, 1e-9), family=family, cycle_type=cycle_type)
        return (r["Qc_W"] if r["feasible"] else 0.0) - Qc_lit

    try:
        mdot_cal = brentq(qc_residual, 1e-9, 5.0, xtol=1e-9)
    except ValueError:
        return {"feasible": False, "mass_total_kg": mass_total_kg,
                "status": "no calibration found (span_fraction>0 for the graded bed, but "
                f"reported Qc={Qc_lit}W still unreachable within mdot in [1e-9, 5.0] kg/s)"}

    result = run_graded_cascade(T_cold_K, span_K, n_stages, mu0H_max=mu0H,
                                 mass_per_stage=mass_total_kg / n_stages, frequency=freq,
                                 fluid_mdot=mdot_cal, family=family, cycle_type=cycle_type)
    result["mdot_calibrated_kg_s"] = round(mdot_cal, 6)
    result["mass_total_kg"] = mass_total_kg
    result["Qc_lit_W"] = Qc_lit
    result["COP_lit"] = None
    return result


def _cooltech_mass_worker(mass_total_kg):
    """Top-level (picklable) worker for run_cooltech_mass_sensitivity()'s
    mass-level pool. validate_cooltech_graded_bed() never opens its own
    internal pool (unlike validate_magqueen_graded_bed()), so there is no
    nested-pool hazard here -- this wrapper exists only so the function
    passed to ProcessPoolExecutor is a plain top-level name, which is a
    pickling requirement, not a correctness one."""
    return validate_cooltech_graded_bed(mass_total_kg=mass_total_kg)


def run_cooltech_mass_sensitivity(masses_kg=(0.5, 1.0, 2.0, 5.0, 10.0), verbose=True,
                                   parallel=True):
    """Same reasoning as run_magqueen_mass_sensitivity(): Cooltech_2013_rotary's
    mass_MCM_kg is not reported, so this sweeps a small range rather than
    asserting one placeholder as this repo's answer. No COP_lit exists for
    this capacity-only row, so "success" here means the graded bed
    reaches positive Qc feasibility at the reported 42K span, not a COP
    match.

    parallel=True (default) runs the swept masses in parallel across one
    shared ProcessPoolExecutor, the same top-level-parallelism approach
    run_magqueen_mass_sensitivity() uses -- see that function's own
    docstring for why mass-level parallelism (rather than per-stage
    parallelism within each mass) is the right level to parallelize at
    here."""
    results = None
    if parallel and len(masses_kg) > 1:
        pool = None
        _blas_env_cm = None
        try:
            _blas_env_cm = _single_threaded_blas_env()
            _blas_env_cm.__enter__()
            pool = ProcessPoolExecutor(max_workers=min(len(masses_kg), os.cpu_count() or 1),
                                        initializer=_pool_worker_init)
            # Phase 31 fix: same reproduced hang as
            # run_magqueen_mass_sensitivity() above (bare `f.result()` with
            # no timeout) -- see this module's "Phase 31" comment block
            # above _pool_submit_all_or_none for the full diagnosis.
            results = _pool_submit_all_or_none(pool, _cooltech_mass_worker, masses_kg)
        except Exception as exc:
            logger.warning("Phase 31: ProcessPoolExecutor creation failed "
                            "(%s: %s) -- falling back to sequential.",
                            type(exc).__name__, exc)
            results = None
        finally:
            _safe_pool_shutdown(pool)
            if _blas_env_cm is not None:
                _blas_env_cm.__exit__(None, None, None)

    if results is None:
        results = [validate_cooltech_graded_bed(mass_total_kg=m) for m in masses_kg]

    if verbose:
        for m, r in zip(masses_kg, results):
            if r.get("feasible"):
                print(f"  mass_total={m:5.2f}kg  mdot_cal={r['mdot_calibrated_kg_s']:.6f}kg/s  "
                      f"Qc={r['Qc_W']:.1f}W (target {r['Qc_lit_W']:.1f}W)  "
                      f"n_fallback_to_Gd={r['n_stages_out_of_range']}/6")
            else:
                print(f"  mass_total={m:5.2f}kg  {r.get('status', 'infeasible')}")
    n_feasible = sum(1 for r in results if r.get("feasible"))
    if verbose:
        print(f"CONCLUSION: {n_feasible}/{len(masses_kg)} swept masses reach positive Qc "
              "feasibility at this device's reported 42K span (this repo's largest benchmark "
              "span) for a 6-stage hypothetical graded Gd-alloy bed -- mass_MCM_kg is genuinely "
              "unreported for this device, so treat this as a feasibility sensitivity, not a "
              "calibrated result.")
    return results


def run_astronautics_cycle_type_sensitivity(apply_correction=None, verbose=True):
    """ROADMAP.md Phase 17 follow-up, closing that phase's own "did NOT
    do" item: cycle_type was validated system-wide in step 2b
    (core.validation_system.run_cycle_type_validation()) but NEVER
    threaded through core/cascade.py's graded-bed helpers, so
    Astronautics_rotary_2014 -- the one device this repo's own step-2
    single-Tc validation flags with the largest COP error on record
    (-81.1%) -- never got the same "ericsson" check its rotary sibling
    DTU_Eriksen_rotary_Gd_2015 got in step 2b (where ericsson narrowed
    -2.1% to +0.6%). This function runs `validate_astronautics_graded_bed()`
    under both cycle_type="brayton" (baseline, unchanged) and
    cycle_type="ericsson" (the naming-convention-heuristic classification
    `core.validation_system.infer_cycle_type_for_device()` would assign
    this device, since "rotary" is in its name) and reports whether the
    graded-bed reproduction's COP error shrinks.

    Returns a dict with both results and the resulting finding. Does NOT
    modify CYCLE_TYPE_FACTORS or claim this settles whether "ericsson" is
    the objectively correct classification for this device (that
    classification is itself only a naming-convention proxy -- see
    core.validation_system.infer_cycle_type_for_device()'s own honesty
    flag, unchanged by this function)."""
    brayton_result = validate_astronautics_graded_bed(
        apply_correction=apply_correction, cycle_type="brayton")
    ericsson_result = validate_astronautics_graded_bed(
        apply_correction=apply_correction, cycle_type="ericsson")

    both_feasible = bool(brayton_result.get("feasible") and ericsson_result.get("feasible"))
    improved = None
    if both_feasible:
        improved = bool(abs(ericsson_result["COP_error_pct"])
                         < abs(brayton_result["COP_error_pct"]))

    if verbose:
        print("Astronautics_rotary_2014 6-layer graded-bed reproduction: "
              "brayton (baseline) vs. ericsson (ROADMAP.md Phase 17 follow-up)")
        if both_feasible:
            print(f"  brayton:  COP_cascade={brayton_result['COP_cascade']}  "
                  f"COP_error={brayton_result['COP_error_pct']}%")
            print(f"  ericsson: COP_cascade={ericsson_result['COP_cascade']}  "
                  f"COP_error={ericsson_result['COP_error_pct']}%")
            if improved:
                print(f"  FINDING: ericsson narrows the error "
                      f"({brayton_result['COP_error_pct']}% -> "
                      f"{ericsson_result['COP_error_pct']}%), consistent with "
                      f"DTU_Eriksen_rotary_Gd_2015's step-2b result -- a second, "
                      f"independent (graded-bed, not single-Tc) data point in the "
                      f"same direction.")
            else:
                print(f"  FINDING: ericsson does NOT narrow the error "
                      f"({brayton_result['COP_error_pct']}% -> "
                      f"{ericsson_result['COP_error_pct']}%) for this device's "
                      f"graded-bed reproduction -- unlike DTU_Eriksen_rotary_Gd_2015's "
                      f"step-2b result. This is a genuine, single-device disagreement, "
                      f"not smoothed over: the naming-convention 'rotary -> ericsson' "
                      f"proxy does not generalize cleanly across both rotary devices "
                      f"checked so far, at least not for the graded-bed's much larger "
                      f"remaining error (-81.1% at baseline, dominated by other "
                      f"documented gaps -- the single-Tc-approximation-vs-6-real-layers "
                      f"issue this function's own docstring describes, and the model's "
                      f"~2.4x DeltaT_ad overestimate documented in "
                      f"giguere_validation.py -- not by cycle topology).")
        else:
            print("  Not comparable: one or both cycle_type runs did not calibrate.")

    return {"brayton": brayton_result, "ericsson": ericsson_result,
            "both_feasible": both_feasible, "ericsson_improves": improved}


def run_astronautics_giguere_correction_sensitivity(cycle_type="brayton", verbose=True):
    """Paper-Mining Pass review item 2: "the Astronautics graded-bed
    reproduction (-81.1%) is the single worst number in the suite, and a
    likely fix already exists in the repo" -- giguere_validation.py found
    the first-order Landau model overestimates DeltaT_ad by ~2.42x at 7T
    and ships DTAD_CORRECTION_FACTOR (~0.41) for exactly this. The review
    asked whether validate_astronautics_graded_bed() actually applies that
    correction to its 6 La(Fe,Si)13Hy stages.

    UPDATE (later pass): validate_astronautics_graded_bed()'s baseline
    error above (-81.1%) has since been superseded -- that function now
    defaults to shared_hardware=True (see run_explicit_material_cascade()/
    run_graded_cascade()'s own docstrings), which fixed an unrelated
    N-times loss-overcounting bug and moved the UNCORRECTED baseline to
    +0.9%. The Giguere-correction experiment below still runs and is still
    informative (it now narrows +0.9% to -0.3%, still an open, unvalidated
    cross-family experiment, not an adopted correction) but the original
    "-81.1%, single worst number in the suite" framing describing why this
    function exists is no longer current -- kept here for the review-item
    paper trail, not as a live number to quote.

    Checked directly (not assumed): it does NOT, by design --
    validate_astronautics_graded_bed()'s own docstring and
    lafesih_composition_tuned_material()'s own docstring
    (core/first_order_mce.py) both already state DTAD_CORRECTION_FACTOR is
    a Gd5Si2Ge2-specific empirical factor "NOT ... shown to transfer to
    this different first-order compound family -- applying it here would
    be fabricating a validation that doesn't exist." That is a documented
    design decision, not an oversight -- but validate_astronautics_graded_bed()
    DOES already expose an `apply_correction` override for exactly this
    experiment, unexercised anywhere in main.py's pipeline. This function
    runs it: apply_correction=None (baseline, current -81.1%) vs.
    apply_correction=DTAD_CORRECTION_FACTOR (the Gd5Si2Ge2-fit value,
    applied here to LAFESIH_FAMILY as a deliberately-labeled EXPERIMENT,
    not a validated correction), and reports whether COP error shrinks --
    while keeping the existing honesty flag about cross-family transfer
    front and center, not quietly dropped just because the number itself
    looks better.
    """
    from core.giguere_validation import DTAD_CORRECTION_FACTOR

    baseline = validate_astronautics_graded_bed(apply_correction=None, cycle_type=cycle_type)
    corrected = validate_astronautics_graded_bed(
        apply_correction=DTAD_CORRECTION_FACTOR, cycle_type=cycle_type)

    both_feasible = bool(baseline.get("feasible") and corrected.get("feasible"))
    improved = None
    if both_feasible:
        improved = bool(abs(corrected["COP_error_pct"]) < abs(baseline["COP_error_pct"]))

    if verbose:
        print("Astronautics_rotary_2014 6-layer graded-bed reproduction: "
              "uncorrected (baseline) vs. Giguere-corrected "
              f"(DTAD_CORRECTION_FACTOR={DTAD_CORRECTION_FACTOR:.4f}, "
              f"cycle_type={cycle_type})")
        if both_feasible:
            print(f"  uncorrected: COP_cascade={baseline['COP_cascade']}  "
                  f"COP_error={baseline['COP_error_pct']}%")
            print(f"  corrected:   COP_cascade={corrected['COP_cascade']}  "
                  f"COP_error={corrected['COP_error_pct']}%")
            if improved:
                print(f"  FINDING: applying the Giguere correction narrows the error "
                      f"({baseline['COP_error_pct']}% -> {corrected['COP_error_pct']}%). "
                      f"This does NOT establish that DTAD_CORRECTION_FACTOR transfers to "
                      f"La(Fe,Si)13Hy -- it was fit to a Gd5Si2Ge2-specific direct-"
                      f"measurement comparison (core/giguere_validation.py) and no "
                      f"equivalent direct-measurement dataset for La(Fe,Si)13Hy was "
                      f"located in this repo's corpus (see "
                      f"lafesih_composition_tuned_material()'s own honesty flag, "
                      f"core/first_order_mce.py) -- but a smaller error IS a smaller "
                      f"error, and this experiment (unlike leaving the option unexercised) "
                      f"actually shows what applying it would do, as an explicitly-labeled "
                      f"sensitivity check rather than an adopted correction.")
            else:
                print(f"  FINDING: applying the Giguere correction does NOT narrow the "
                      f"error ({baseline['COP_error_pct']}% -> {corrected['COP_error_pct']}%) "
                      f"-- consistent with the existing honesty flag that this factor was "
                      f"never shown to transfer to La(Fe,Si)13Hy. The -81.1% baseline error "
                      f"is not explained by the same DeltaT_ad-overestimate mechanism "
                      f"giguere_validation.py found for Gd5Si2Ge2; the single-Tc-per-stage "
                      f"approximation of the real 6-real-layer device (see this function's "
                      f"own docstring) remains the more likely dominant cause.")
        else:
            print("  Not comparable: one or both apply_correction runs did not calibrate.")

    return {"baseline": baseline, "corrected": corrected,
            "dtad_correction_factor": DTAD_CORRECTION_FACTOR,
            "both_feasible": both_feasible, "correction_improves": improved}


if __name__ == "__main__":
    from core.first_order_mce import GD5SI2GE2_FIRST_ORDER

    print("Cascade AMR staging vs. baselines, ASHRAE 5-20K span sweep")
    print("(mu0H=2T per stage, 5kg regenerator per stage, f=1Hz, mdot=0.08kg/s, NTU thermal model on)")
    print("=" * 100)
    print("\n--- Material: Gd (baseline) ---")
    rows_gd = compare_staging(material=GADOLINIUM, mass_per_stage=5.0,
                                out_csv="results/cascade_comparison.csv")
    header = f"{'span':>5} {'1-stage':>9} {'2-stage':>9} {'3-stage':>9} {'4-stage':>9} {'VCC':>7} {'Liquid':>7}"
    print(header)
    for r in rows_gd:
        def fmt(v):
            return f"{v:9.2f}" if v is not None else f"{'--':>9}"
        print(f"{r['span_K']:>5} {fmt(r['AMR_1stage_COP'])} {fmt(r['AMR_2stage_COP'])} "
              f"{fmt(r['AMR_3stage_COP'])} {fmt(r['AMR_4stage_COP'])} "
              f"{r['VaporCompression_COP']:>7} {r['LiquidCooling_COP']:>7}")
    print("Wrote results/cascade_comparison.csv")

    print("\n--- Material: Gd5Si2Ge2 (giant MCE) ---")
    rows_giant = compare_staging(material=GD5SI2GE2_FIRST_ORDER, mass_per_stage=5.0,
                                   out_csv="results/cascade_comparison_giant_mce.csv")
    print(header)
    for r in rows_giant:
        print(f"{r['span_K']:>5} {fmt(r['AMR_1stage_COP'])} {fmt(r['AMR_2stage_COP'])} "
              f"{fmt(r['AMR_3stage_COP'])} {fmt(r['AMR_4stage_COP'])} "
              f"{r['VaporCompression_COP']:>7} {r['LiquidCooling_COP']:>7}")
    print("Wrote results/cascade_comparison_giant_mce.csv")

    gd_10K = next(r for r in rows_gd if r["span_K"] == 10)
    giant_10K = next(r for r in rows_giant if r["span_K"] == 10)
    print(f"\nAt 10K span: Gd 1-stage COP={gd_10K['AMR_1stage_COP']} vs. "
          f"Gd5Si2Ge2 1-stage COP={giant_10K['AMR_1stage_COP']} "
          f"(VCC={gd_10K['VaporCompression_COP']}, Liquid={gd_10K['LiquidCooling_COP']})")

    print("\n" + "=" * 100)
    print("--- Curie-graded cascade (ROADMAP.md Phase 7 open item) ---")
    print("=" * 100)
    print("Each stage uses a hypothetical composition-tuned Gd5(SixGe1-x)4(-Ga) material")
    print("whose own peak MCE effect is matched (via iterative peak search, see")
    print("_target_composition_for_peak) to that stage's local operating temperature,")
    print("checked against the literature-documented Tc range (20-290K) and scaled by the")
    print("Giguere et al. (1999) empirical correction (core.giguere_validation), so these")
    print("numbers are not built on the raw model's ~2.4x-optimistic DeltaT_ad.\n")

    rows_graded, stage_info_all = compare_graded_cascade(
        T_cold_C=18.0, spans=range(5, 21), mass_per_stage=5.0,
        out_csv="results/graded_cascade_comparison.csv")

    example = next(s for s in stage_info_all if s["span_K"] == 10 and s["n_stages"] == 3)
    print("Example (10K span, 3 stages):")
    for s in example["stage_info"]:
        print(f"    stage {s['stage']}: T_mid={s['T_mid_K']}K, needed composition Tc="
              f"{s['Tc_target_K']}K -> {s['material']}")

    graded_10K_3 = next(r for r in rows_graded if r["span_K"] == 10)
    gd_10K_3 = next(r for r in rows_gd if r["span_K"] == 10)
    print(f"\nAt this point: Graded 3-stage Qc={graded_10K_3['Graded_3stage_Qc_W']}W, "
          f"COP_elec={graded_10K_3['Graded_3stage_COP']}  vs.  plain-Gd 3-stage "
          f"Qc={gd_10K_3['AMR_3stage_Qc_W']}W, COP_elec={gd_10K_3['AMR_3stage_COP']}")
    print("-> consistent with giant_mce_analysis.py's earlier finding: a bigger MCE mostly")
    print("   buys more Qc per kg (here: substantially more), not a materially better COP")
    print("   (loss_model.py's field/frequency/flow-dependent parasitics dominate COP either way).")

    n_cells = sum(1 for row in rows_graded for n in (1, 2, 3, 4))
    n_full_range = sum(1 for row in rows_graded for n in (1, 2, 3, 4)
                        if row[f"Graded_{n}stage_n_fallback_to_Gd"] == 0)
    n_some_fallback = sum(1 for row in rows_graded for n in (1, 2, 3, 4)
                           if 0 < row[f"Graded_{n}stage_n_fallback_to_Gd"] < n)
    n_all_fallback_or_infeasible = n_cells - n_full_range - n_some_fallback
    print(f"\nAcross the full 5-20K span x 1-4 stage-count sweep ({n_cells} cells):")
    print(f"  {n_full_range} cells: every stage's needed composition is within the "
          f"documented 20-290K giant-MCE range")
    print(f"  {n_some_fallback} cells: SOME stages exceed 290K and fall back to plain Gd "
          f"for that stage only (larger spans and/or more stages push the hottest stage's "
          f"needed Tc above the Ga-alloyed ceiling)")
    print(f"  {n_all_fallback_or_infeasible} cells: fully infeasible (all stages fell back, "
          f"or Qc collapsed to ~0)")
    print("\nHONEST CAVEAT: this Landau model's transition turns out to be numerically much")
    print("narrower (DeltaT_ad falls off within a few tenths of a K of its peak at mu0H=2T)")
    print("than the real, hysteresis/inhomogeneity-broadened transition Giguere et al.'s own")
    print("Fig. 3 shows (spread over ~10-15K) -- an idealized-model limitation on top of the")
    print("dTad-magnitude one already flagged, not smoothed over here. A small number of")
    print("individual span/stage-count cells show a stage's Qc collapsing to ~0 despite its")
    print("composition nominally being in-range, from residual peak-alignment error after the")
    print("iterative search (see run_graded_cascade's COP_cascade=0.0 with Qc_W>0 rows in the")
    print("output CSV) -- flagged here as a real numerical fragility of this idealized 6th-")
    print("order Landau fit, not hidden by rounding it away.")
    print("\nWrote results/graded_cascade_comparison.csv")
    print("\nBOTTOM LINE: unlike the fixed Gd5Si2Ge2 comparison above (which collapses to zero")
    print("everywhere in the ASHRAE range because Gd5Si2Ge2's own Tc=276K is fixed and far from")
    print("the operating point), a CURIE-GRADED cascade -- built from literature-documented")
    print("composition tunability and validated against Giguere et al.'s direct measurement --")
    print("genuinely delivers several times the cooling capacity of plain Gd at comparable COP")
    print("for smaller spans/stage-counts, but the same giant-MCE family's documented ~290K")
    print("composition ceiling still constrains larger spans and higher stage counts, echoing")
    print("giant_mce_analysis.py's conclusion that the ASHRAE range sits right at the edge of")
    print("what this material family is documented to reach.")

    print("\n" + "=" * 100)
    print("--- Phase 9: does a 6-layer Curie-graded La(Fe,Si)13Hy bed reproduce the REAL")
    print("    Astronautics_rotary_2014 device? (validation_system.py's single-Tc=287K")
    print("    material could not -- see ROADMAP.md Phase 9) ---")
    print("=" * 100)
    astro = validate_astronautics_graded_bed()
    if astro.get("feasible"):
        print("Layer Tc targets (evenly spread across the device's reported 303.6-316.2K "
              "layer range):")
        for s in astro["stage_info"]:
            print(f"    stage {s['stage']}: T_mid={s['T_mid_K']}K, needed composition Tc="
                  f"{s['Tc_target_K']}K -> {s['material']}")
        print(f"\nmdot calibrated to reproduce reported Qc={astro['Qc_lit_W']}W: "
              f"{astro['mdot_calibrated_kg_s']} kg/s")
        print(f"Predicted COP={astro['COP_cascade']}  vs.  reported COP={astro['COP_lit']} "
              f"({astro['COP_error_pct']:+.1f}% error)")
        print("\n-> a comparable-magnitude error to the Gd devices in validation_system.py's own")
        print("   point-wise validation, and a very different outcome from the flat 'no")
        print("   calibration found' the single-layer LAFESIH_FIRST_ORDER material gave this")
        print("   same device: the graded-bed STRUCTURE, not just the material, was the missing")
        print("   piece. Still an approximation -- layer Tc's are evenly spread across the")
        print("   reported range, not the paper's actual (unpublished here) per-layer values,")
        print("   and mdot is calibrated rather than predicted, same caveat as")
        print("   validation_system.py's own methodology throughout.")
    else:
        print(astro.get("status", "infeasible"))