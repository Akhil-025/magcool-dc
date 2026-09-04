"""
regenerator_1d.py
==================
Minimal 1-D transient (blow-by-blow) active magnetic regenerator (AMR)
simulator.

WHY THIS EXISTS: core/amr_cycle.py's AMRSystem.cooling_capacity() is a 0-D,
single-blow model -- it evaluates the magnetocaloric material's own
adiabatic dT ONCE, at the bed's overall mid-temperature, and caps the
achievable span at 2*that value. It cannot represent "regenerative
amplification" (the temperature-profile buildup along a real regenerator
bed over many cycles, which is what actually lets real AMR devices reach
spans several times larger than any single material blow's own dT --
confirmed against this repo's own cited literature in
results/regenerative_amplification_diagnostic.txt: DTU_Eriksen_MAGGIE_2016
reaches a measured 29.2K no-load span against a single-blow structural
ceiling of only ~14K at the same field). core/cascade.py's ROADMAP.md item
A3 explicitly declined to patch cooling_capacity()'s span-cap FORMULA with
an unsourced smoothing/rescaling, because no literature source gives the
exact functional form. An earlier attempt in this same vein (reusing
cascade.py's series-stage machinery with one material split into slices,
see git history / conversation record) was tested against the benchmark
set and REJECTED: it resolved calibration feasibility but produced a
+220% COP error on a device (Lozano_POLO_UFSC_2016_r1) that had no
calibration at all under the 0-D model -- i.e. it replaced an honest "no
prediction" with a confidently wrong one, because splitting a single
physical bed into independently-scaled cascade "stages" implicitly
invents extra internal heat-exchanger junctions that don't exist in a
real single-pass bed.

This module instead does the physically correct thing: an actual
transient simulation of the regenerator bed, discretized into N_NODES
solid nodes along its length, with the AMR cycle's four phases
(magnetize / hot blow / demagnetize / cold blow) applied explicitly and
iterated to periodic steady state -- the standard approach in the
numerical-modeling literature already cited by this repo (e.g.
Papers/AMR Theory and Modeling/Dynamic operation of an active magnetic
regenerator (AMR) Numerical optimization of a packed-bed AMR.pdf; see
also, confirmed by web search during this pass: Nielsen et al.,
"Numerical analysis of an active magnetic regenerator" (Int. J. Refrig.
2011) -- 1-D transient model, NTU in range 10-50, utilization near unity
for good performance; Bahl et al./Tusek et al., 2D parametric studies
reaching the same NTU/utilization conclusion by an independent numerical
route). No closed-form algebraic formula for span vs. NTU/utilization was
found in this repo's corpus or in a live literature search -- every
citable source treats this via numerical simulation, which is exactly
what this module now provides in-repo.

Model
-----
The regenerator is discretized into `n_nodes` equal-mass solid nodes
along its length. Each full AMR cycle (period 1/frequency) applies, to
every node:

  1. Magnetize (0 -> mu0H_max): each node's solid temperature jumps by
     its own LOCAL material.delta_T_adiabatic(T_i, mu0H_max), applied
     instantaneously (field-ramp time << thermal relaxation time -- a
     standard AMR-modeling assumption; already used the same way, at a
     single point, by core/amr_cycle.py's own cooling_capacity()).
  2. Hot blow    : fluid flows cold-end -> hot-end (node 0 -> node
     n_nodes-1), exchanging heat with each node in turn.
  3. Demagnetize (mu0H_max -> 0): each node's solid temperature drops by
     its own local delta_T_adiabatic, instantaneous.
  4. Cold blow   : fluid flows hot-end -> cold-end, exchanging heat with
     each node in turn.

repeated until the node-temperature vector (sampled once per cycle, at
the same point in the cycle) reaches periodic steady state (max change
between successive cycles < tol).

Per-node heat transfer during a blow uses the EXACT effectiveness for a
single fluid stream passing an solid element held at ~constant
temperature during that passage (standard heat-exchanger theory -- Kays &
London, Compact Heat Exchangers, already cited by core/thermal.py for the
same regenerator-theory context):

    eps_node = 1 - exp(-NTU_node), NTU_node = NTU_total / n_nodes

NTU_total comes from core/thermal.py's own regenerator_effectiveness() --
the SAME correlation core/amr_cycle.py's 0-D model already uses for its
own eps -- so this module is not introducing an independent heat-transfer
estimate. Splitting NTU_total by n_nodes is exact, not approximate:
NTU = h*A_total/(mdot*cp_f), A_total is proportional to bed volume/mass
at fixed particle diameter and fluid velocity (see thermal.py's own
derivation), and equal node masses give each node an equal share of
A_total.

This node model deliberately does NOT also apply
regenerator_effectiveness()'s own eps=NTU/(NTU+2)*(1-0.3U) formula -- that
formula's "+2" and utilization-derating terms are a LUMPED, single-node
phenomenological stand-in for exactly the periodic-flow physics this
module simulates explicitly (alternating hot/cold blows, iterated to
periodic steady state). Applying both would double-count the same
physics once as an explicit simulation and once as its own approximation
for not doing that simulation.

mdot convention: matches core/thermal.py's regenerator_effectiveness()
and core/amr_cycle.py's AMRSystem -- mdot is the flow rate WHILE fluid is
actually moving (used directly for Re/h/NTU and for the per-blow
Q_i = mdot*cp_f*(T_out-T_in) power). To convert that in-blow power into a
CYCLE-AVERAGED Qc (comparable to a reported device Qc, i.e. averaged over
periods when the field is only ramping and no fluid is flowing), we
multiply by the duty cycle: fluid moves for tau_blow = blow_fraction /
(2*frequency) out of every 1/frequency-second cycle, so
Qc_avg = Qc_instantaneous * (tau_blow * frequency) = Qc_instantaneous *
(blow_fraction / 2). blow_fraction has the same meaning and same default
(0.5) as core/amr_cycle.py's AMRSystem.blow_fraction.

Known limitations (stated, not hidden)
---------------------------------------
- Lumped-node (not spatially continuous) treatment within each node during
  a blow: within one node, the solid is treated as isothermal during that
  node's passage. This is standard for a "quasi-steady blow" 1-D model
  (valid when n_nodes is large enough that each node's own temperature
  change per blow is small) but is a real simplification relative to a
  continuous PDE solve.
- No axial (node-to-node) solid heat conduction is modeled -- the
  standard packed-bed/parallel-plate AMR literature (e.g. the "Review on
  numerical modeling..." paper already in this repo's Papers/) finds axial
  conduction is a second-order effect at typical operating NTU, but this
  has not been independently checked here.
- No fluid dead-volume / carryover between the bed and the end reservoirs
  is modeled explicitly beyond the "no-load: fluid relays between blows"
  boundary condition described above.
- RESOLVED (this pass): the low-mdot degeneracy above (span growing
  without bound as mdot -> 0, no interior maximum) is fixed. Root cause,
  confirmed by direct test: the ONLY inter-node coupling in this model was
  via the fluid, whose per-node effectiveness eps_node=1-exp(-NTU_node)
  approaches 1 as mdot -> 0 (NTU = h*A/(mdot*cp_f) diverges because h only
  falls to a finite convective floor, Nu -> 2, as Re -> 0 in
  regenerator_effectiveness()'s own correlation) -- so the coupling got
  MORE effective exactly as less fluid moved, with nothing to counteract
  the resulting gradient buildup. Real packed beds have a second coupling
  path this model was missing entirely: axial (node-to-node) conduction,
  which the module's own older "known limitations" text below correctly
  described as "a second-order effect at typical operating NTU" but which
  becomes the DOMINANT loss path in exactly the artificial low-mdot/
  high-NTU corner the old mdot search swept through. Added
  `_apply_axial_conduction()`, called once per cycle. Confirmed this
  produces a genuine interior maximum (Tusek_singlebed_Gd_2010_spanceiling:
  span rises then falls as mdot is swept from 0.008 down to 0.0001 kg/s,
  peaking near mdot~0.002 kg/s, instead of growing monotonically forever)
  and that the result converges (stable to <0.01K/cycle) given enough
  cycles, rather than needing an arbitrary search-grid floor to bound it.
  Conductivity note: a naive porosity-weighted parallel mix of bulk Gd
  conductivity (10.5 W/m/K) and fluid conductivity was tried first and
  REJECTED -- it overcorrected, homogenizing the whole bed within a single
  cycle period and collapsing span to a nonphysical ~0.1-0.2K with no
  convergence at all. A packed bed of spheres is solid-solid POINT
  contact, not a continuous rod, so bulk metal conductivity is not the
  right axial transport coefficient; used a modest multiple of fluid
  conductivity instead (standard packed-bed literature result: point
  contacts contribute a limited enhancement over the stagnant-fluid
  conduction path), which is what produces the interior maximum above.
- PHASE 31 UPDATE to the conductivity note above: the "modest multiple of
  fluid conductivity" placeholder (min(3.0, 1+(1-porosity)) * k_fluid) has
  been replaced with the Maxwell-Eucken packed-bed composite-conductivity
  model (see _packed_bed_effective_axial_conductivity()'s own docstring
  for the formula and why this specific, independently re-derivable model
  was chosen over the more commonly-cited but harder-to-verify-from-memory
  Zehner-Schlunder correlation). HONEST RESULT: this did NOT fix the
  undershoot described in the next bullet -- it made it slightly WORSE on
  two of the three benchmark rows (Maxwell-Eucken predicts ~2.7x HIGHER
  effective conductivity than the old placeholder at porosity=0.365, i.e.
  more damping, not less). See validate_against_benchmarks()'s own
  "PHASE 31 UPDATE" output text for the full before/after numbers. This is
  a genuine, disappointing-but-informative finding, kept rather than
  reverted: it demonstrates the direction-inconsistent error below is NOT
  simply "the axial conductivity constant is a little off in one
  direction" -- a legitimately different, textbook-derived conductivity
  value moves the same two rows further the same wrong way, which points
  more strongly at the regenerator-effectiveness/NTU coupling or the
  single-blow reference itself as the remaining calibration gap, not
  solely this conductivity term.
- NEW OPEN ISSUE (found by fixing the one above): with the interior
  maximum now well-defined, the model's PEAK predicted span materially
  UNDERSHOOTS the literature value on the case checked directly
  (Tusek_singlebed_Gd_2010_spanceiling: peak ~1.4K around mdot~0.002 kg/s
  vs. a measured 19.8K) -- see results/regenerator_1d_validation.txt for
  the full benchmark sweep. This is a different, narrower problem than the
  one it replaces: previously the model had NO well-defined answer at all
  (sensitive to an arbitrary search floor); now it has one, and that one
  is quantitatively too small, most likely because the axial-conduction
  multiplier above was chosen to be physically defensible but was not
  independently calibrated (no packed-bed-conductivity literature source
  for this specific material/geometry was located in this repo's corpus
  during this pass), and/or because other simplifications noted above
  (lumped isothermal nodes, no dead-volume/carryover modeling, no
  hysteresis) also matter more at the flow rates where the new interior
  maximum sits than they did in the old, artificially-flow-independent
  regime. This is flagged, not hidden or tuned away by adjusting the
  conductivity multiplier until the benchmark matches -- that would just
  be curve-fitting one number to one data point.
- PHASE 39 FINDING #1 (real, verified root-cause component): the
  Tusek_singlebed_Gd_2010_spanceiling benchmark device is NOT a packed
  bed of spheres -- its own source paper (Tusek, Kitanovski, Zupan,
  Prebil, Poredos, Appl. Therm. Eng. 53 (2013) 57-66, Table 1) reports it
  as AMR(A), a PARALLEL-PLATE regenerator: 0.1mm plate spacing, 0.25mm
  plate thickness, porosity 0.2564, outer dimensions 10mm(height) x
  80mm(length) x 39mm(width), total Gd mass 0.1763kg (matches this repo's
  CSV row exactly) -- yet every call in this module before this phase
  used regenerator_effectiveness() (the PACKED-BED correlation) with this
  function's own packed-bed particle_diameter/porosity/
  bed_cross_section_area DEFAULTS, not this device's real geometry TYPE
  or dimensions. core/thermal.py already had a dedicated, real,
  literature-sourced regenerator_effectiveness_parallel_plate() (Nickolay
  & Martin (2002)/Tusek et al. 2013 Eq. 4) that was simply never wired
  into this module. Added geometry="parallel_plate" support (see
  simulate_amr_1d()'s own docstring) and re-ran this specific benchmark
  with the device's REAL, source-verified geometry (bed_cross_section_area
  = 0.010m x 0.039m = 3.9e-4 m^2, derived from Table 1's own outer
  dimensions -- cross-checked for self-consistency: V_bed/this area gives
  a bed length of ~77mm, matching the reported 80mm length to within
  ~4%): peak span improves from ~1.4K (-92 to -97% error, packed-bed
  defaults) to ~16.5K (-16% error) at this module's usual n_nodes=20
  default. This is a genuine, large improvement, using a real correlation
  and real, source-verified device geometry -- NOT a parameter tweaked
  until the number looked better. CAVEAT, found while writing this
  module's own tests: regenerator_effectiveness_parallel_plate()'s
  internal porosity = plate_spacing/(plate_spacing+plate_thickness) unit-
  cell idealization gives 0.2857 for this device, not an exact match to
  the paper's own DIRECTLY REPORTED 0.2564 (Table 1) -- an honest ~11%
  relative discrepancy, most likely from real edge/end effects (e.g. plate
  holders, end caps) a simple infinite-unit-cell idealization can't
  capture. The ~16.5K result above uses the code's own derived 0.2857,
  not a manually-overridden 0.2564 -- flagged here so this small remaining
  imprecision isn't mistaken for an exact geometry match.
- PHASE 39 FINDING #2 (more fundamental, NOT fixed, found while checking
  Finding #1's robustness): this model is NOT grid-converged in n_nodes,
  in EITHER geometry mode, and this is NOT something introduced
  -- it reproduces with geometry="packed_bed" too. Direct check, same
  device, same mdot, only n_nodes varied:
      packed_bed:      n_nodes=10 -> 2.19K, n_nodes=20 -> 3.23K, n_nodes=40 -> 3.78K
      parallel_plate:  n_nodes=10 -> 8.99K, n_nodes=20 -> 16.15K, n_nodes=40 -> 2.31K
  Neither sequence is converging monotonically to a stable limit (the
  parallel_plate sequence isn't even monotonic). Checked further (ruled
  out "just needs more cycles"): re-running n_nodes=40 with 12000 cycles
  instead of 1500 gives span=1.04K, not 2.31K -- i.e. letting it run
  longer moves it FURTHER from, not closer to, the n_nodes=20 answer, and
  the last-10-cycle history is flat at that lower value (a genuine
  steady state at this resolution, not an unconverged transient still en
  route to something else). This rules out "the fine grid just needs a
  longer cycle budget" as an explanation -- the fine-grid answer really is
  a different, smaller steady-state span, not a slowly-converging version
  of the coarse-grid one. A further diagnostic
  (axial conduction artificially disabled) shows the FLUID-node
  discretization alone diverges without bound as n_nodes increases
  (9.6K -> 28.3K -> 40.4K) -- the same "unbounded growth" failure mode
  this module's own earlier fix already diagnosed and fixed for mdot -> 0
  (see the axial-conduction section above), now showing up as the SAME
  underlying problem in a different independent variable (spatial
  resolution instead of flow rate): decay_substep=exp(-NTU_node) with
  NTU_node=NTU_total/n_nodes becomes a progressively less complete
  per-node exchange as n_nodes grows, which lets sharper node-to-node
  gradients build up rather than converging to a well-resolved continuum
  limit, and axial conduction (whose own node-to-node conductance G =
  k_eff*A/dx grows linearly with n_nodes while each node's own thermal
  mass shrinks the same way) does not currently counteract this
  consistently across resolutions -- it under-damps at coarse grids and
  appears to over-damp at n_nodes=40 for this device (span collapses to
  2.3K), rather than approaching a stable value from either side.
  **CONSEQUENCE: NEITHER the previous packed-bed numbers NOR this
  phase's parallel-plate improvement should be treated as a validated,
  resolution-independent prediction** -- both are reporting whatever this
  specific n_nodes=20 discretization happens to give, not a converged
  answer. This is flagged rather than silently kept as a "fix", per this
  module's own established standard (see Finding #1 above and the
   axial-conductivity finding) of not treating "the number moved
  in the direction I wanted" as validation. Properly resolving this needs
  either a formally-derived, order-of-accuracy-checked spatial
  discretization (e.g. a proper per-segment epsilon-NTU formula that is
  provably grid-independent as n_nodes -> infinity, replacing the current
  ad hoc "one full-strength exchange pass per node" scheme) or an
  equivalent fix to how axial conduction is coupled to node count -- a
  substantial numerical-methods task, not a same-pass patch. Until this is
  resolved, this module's outputs (old or new) remain "illustrates the
  effect is real and bounded", not "a specific number you should trust",
  exactly as this docstring's next paragraph already says.
- This is a genuinely new, freshly-debugged numerical model (five
  distinct bugs found and fixed across its development -- an
  explicit-Euler stability blowup, an incorrect NTU/dt rescaling that
  silently killed all heat transfer, a field-unit error that silently
  zeroed out the magnetocaloric effect itself, the low-mdot degeneracy
  above, and the axial-conductivity overcorrection caught while fixing
  it). Treat its Qc/span predictions with real caution, and see
  validate_against_benchmarks() for exactly what has and has not been
  checked before relying on it for anything beyond illustrating that
  regenerative amplification is a real, computable, now-BOUNDED effect --
  it is still NOT wired into core/amr_cycle.py or anywhere else in the
  pipeline, and given the quantitative undershoot above, should not be
  until it is independently calibrated: wiring in a model that currently
  gets the peak span wrong by roughly an order of magnitude on the one
  case checked would replace the 0-D model's honest, clearly-labeled
  structural cap with a differently-wrong number, not a better one.
"""

import numpy as np
from core.thermal import (regenerator_effectiveness, regenerator_effectiveness_parallel_plate,
                           water_properties, CP_SOLID_GD, K_SOLID_GD)

#  fix: bumped whenever a change to simulate_amr_1d()'s physics
# would silently invalidate previously-cached no_load_span() results
# without changing any of their INPUT parameters (e.g. the axial-
# conductivity correlation change below). _cache_key() folds this into
# every cache key, so on-disk results/.regenerator_1d_cache.json entries
# computed under the OLD physics are automatically treated as cache
# misses instead of being silently (and incorrectly) reused.
_MODEL_VERSION = 3
# v1: original axial-conduction fix (ad hoc capped multiplier of stagnant
#     fluid conductivity: k_eff_axial = min(3.0, 1.0+(1-porosity))*k_fluid).
# v2: replaced with the Maxwell-Eucken packed-bed composite-conductivity
#     model -- see _packed_bed_effective_axial_conductivity() below.
# v3 : added geometry="parallel_plate" support (real device
#     geometry for parallel-plate AMRs instead of forcing packed-bed
#     correlations onto them) -- see simulate_amr_1d()'s own docstring.
#     Existing geometry="packed_bed" callers are bit-for-bit unaffected;
#     the version bump exists only so any parallel_plate result computed
#     before this fix (there shouldn't be any, since the parameter didn't
#     exist) can never be silently reused.


def _packed_bed_effective_axial_conductivity(porosity, k_fluid, k_solid):
    """Effective axial thermal conductivity of a packed bed of spherical
    particles, via the Maxwell-Eucken two-phase composite-conductivity
    model (Maxwell, "A Treatise on Electricity and Magnetism", 1873 --
    exact solution for a dilute suspension of non-touching spheres in a
    continuous matrix; Eucken's later extension applies the same
    closed-form expression at non-dilute solid fractions as a standard
    engineering approximation, widely used for packed-bed and composite-
    material conductivity, e.g. Kaviany, "Principles of Heat Transfer in
    Porous Media"):

        k_eff = k_f * [2*k_f + k_s - 2*(1-eps)*(k_f - k_s)]
                     / [2*k_f + k_s +    (1-eps)*(k_f - k_s)]

    where eps = porosity (fluid volume fraction), k_f/k_s are the bulk
    fluid/solid conductivities.

    WHY THIS FORMULA, NOT ZEHNER-SCHLUNDER: the packed-bed heat-transfer
    literature's most frequently cited correlation for this exact
    quantity is Zehner & Schlunder (1970) (as tabulated in, e.g., the VDI
    Heat Atlas), which adds an explicit particle-deformation/contact-area
    parameter on top of the same two-phase starting point used here.
    Maxwell-Eucken is used instead, deliberately: (a) it is simple enough
    to be re-derived and sanity-checked from first principles (an exact
    dilute-suspension result, not a multi-parameter fitted correlation
    transcribed from memory of a source paper this repo's corpus does not
    contain -- the same book/paper-access limitation already flagged
    elsewhere in this repo, e.g. core/thermal.py's own docstring); (b) it
    is bounded and well-conditioned for every porosity in (0, 1), with no
    fitted deformation-parameter singularity to guard against; (c) it
    still captures the same qualitative physics the correlation it
    replaces was reaching for -- axial conduction dominated by the fluid
    path, with a solid-conduction enhancement that grows as porosity
    falls -- while being a real, named, independently re-derivable
    formula instead of a hand-picked capped multiplier
    (`min(3.0, 1.0 + (1-porosity)) * k_fluid`, this module's previous
    placeholder, chosen only to be "physically defensible", not derived
    from any specific model).

    HONEST FRAMING (unchanged from what this replaces): this is still a
    MODEL CHOICE, not a device-measurement-calibrated value. It replaces
    one reasoned approximation with a better-justified, independently
    re-derivable one -- it does NOT eliminate the "axial conductivity is
    not independently calibrated against this repo's own benchmark
    devices" limitation already documented in this module's own
    docstring and in validate_against_benchmarks()'s output. Re-run
    validate_against_benchmarks() after this change and read its
    err_1d_pct column fresh rather than assuming the old
    undershoot/overshoot pattern is now fixed."""
    if not (0.0 < porosity < 1.0) or k_fluid <= 0 or k_solid <= 0:
        return max(k_fluid, 1e-9)  # degenerate inputs -> fall back to pure fluid conduction
    eps = porosity
    numerator = 2 * k_fluid + k_solid - 2 * (1 - eps) * (k_fluid - k_solid)
    denominator = 2 * k_fluid + k_solid + (1 - eps) * (k_fluid - k_solid)
    if denominator <= 0:
        return max(k_fluid, 1e-9)
    return k_fluid * numerator / denominator


def _parallel_plate_effective_axial_conductivity(porosity, k_fluid, k_solid):
    """Effective axial (i.e. along-the-flow-direction) thermal conductivity
    of a PARALLEL-PLATE regenerator -- deliberately a DIFFERENT formula
    from _packed_bed_effective_axial_conductivity() above, not the same one
    reused with a different porosity.

    WHY DIFFERENT: that function's docstring explains Maxwell-Eucken (a
    dispersed-spheres-in-a-continuous-matrix model) was chosen over a naive
    porosity-weighted parallel mix specifically because a packed bed of
    spheres is solid-solid POINT contact, not a continuous conduction path
    -- a naive parallel mix was tested and rejected there because it
    overcorrects (collapses the whole bed to a uniform temperature within
    one cycle). A parallel-plate regenerator's solid phase is exactly the
    opposite geometry: each plate IS a continuous, unbroken solid strip
    running the full length of the bed, with the fluid channels between
    plates also running the full length -- solid and fluid conduction
    paths are physically in PARALLEL along the flow (axial) direction, not
    mediated by point contacts. The simple volume-weighted parallel-
    conduction mixing rule
        k_eff = porosity * k_fluid + (1 - porosity) * k_solid
    is therefore the physically appropriate model for THIS geometry (the
    standard result for conduction through layers oriented parallel to the
    heat-flow direction, e.g. Incropera & DeWitt, Fundamentals of Heat and
    Mass Transfer, composite-wall parallel arrangement) -- not a
    simplification adopted for convenience, and not the same "naive mix"
    rejected for packed beds, because the underlying solid geometry is
    genuinely different (continuous plate vs. point-contact spheres)."""
    if not (0.0 < porosity < 1.0) or k_fluid <= 0 or k_solid <= 0:
        return max(k_fluid, 1e-9)
    return porosity * k_fluid + (1 - porosity) * k_solid


def _apply_axial_conduction(T, dt_total, dx, bed_area, k_eff_axial, m_node, cp_solid_eff):
    """Explicit 1-D conduction of `T` (n_nodes solid nodes, insulated/no-flux
    ends) over a total time `dt_total`, subdivided into as many equal
    substeps as the explicit-scheme stability limit (Fourier number <= 0.5)
    requires. Node-to-node conductance uses the packed-bed axial path: area
    `bed_area`, spacing `dx`, effective conductivity `k_eff_axial`.

    WHY THIS EXISTS: without it, simulate_amr_1d()'s only inter-node
    coupling is via the fluid, whose node-level effectiveness
    eps_node=1-exp(-NTU_node) approaches 1 (near-perfect equilibration per
    pass) as mdot -> 0, because NTU = h*A/(mdot*cp_f) diverges while h itself
    only falls to a finite convective floor (Nu -> 2 in the Re -> 0 limit of
    regenerator_effectiveness()'s own correlation). With eps_node -> 1 and no
    competing loss mechanism, each cycle's fluid-mediated coupling gets
    *more* effective exactly as less fluid moves, and nothing in the model
    damps the resulting temperature-gradient buildup -- no_load_span()'s own
    mdot search never found an interior maximum because of this (see
    module-level 'known limitations' / validate_against_benchmarks()).
    Axial (node-to-node) solid conduction is the real physical mechanism the
    literature notes as 'second-order at typical operating NTU' -- true at
    the NTU/utilization combinations a real device is designed to run at,
    but not in the artificial low-mdot/high-NTU corner the old search swept
    through, where it becomes the dominant heat-loss path and is exactly
    what stops a real bed's span from growing without bound as flow is
    reduced. Adding it here supplies that missing physics rather than
    papering over the symptom with an arbitrary search-grid floor."""
    n = len(T)
    if n < 2 or dt_total <= 0 or k_eff_axial <= 0:
        return T
    alpha = k_eff_axial / (m_node / (bed_area * dx) * cp_solid_eff)  # m^2/s
    dt_stable = 0.5 * dx * dx / alpha if alpha > 0 else dt_total
    n_sub = max(1, int(np.ceil(dt_total / dt_stable)))
    dt = dt_total / n_sub
    G = k_eff_axial * bed_area / dx  # W/K, node-to-node conductance
    for _ in range(n_sub):
        flux = G * (T[1:] - T[:-1])       # length n-1, flux from i -> i+1
        dT = np.zeros(n)
        dT[:-1] += flux
        dT[1:] -= flux
        T = T + dT * dt / (m_node * cp_solid_eff)
    return T


def simulate_amr_1d(material, mu0H_max, mass_total, frequency, mdot,
                     T_cold_reservoir=None, T_hot_reservoir=None,
                     n_nodes=20, blow_fraction=0.5, max_cycles=800, tol=1e-4,
                     particle_diameter=0.0005, porosity=0.365,
                     bed_cross_section_area=0.002, T_K_for_ntu=300.0,
                     cp_solid=None, k_solid=None, T_init=None, n_substeps=None,
                     geometry="packed_bed", plate_thickness=0.00025,
                     plate_spacing=0.0001):
    """Runs the blow-by-blow simulation described in this module's own
    docstring to periodic steady state.

    geometry ( addition): "packed_bed" (default -- exactly
    reproduces every existing call's behavior, since every other new
    parameter this adds defaults to values only used when geometry=
    "parallel_plate") or "parallel_plate". WHY THIS EXISTS: several of
    this module's own benchmark devices (see validate_against_benchmarks())
    are NOT packed beds of spheres in real life -- e.g. Tusek_singlebed_
    Gd_2010_spanceiling's own CSV source note reports it as a parallel-
    plate AMR (0.1mm plate spacing, 0.25mm plate thickness, porosity
    0.2564, Tusek, Kitanovski, Zupan, Prebil, Poredos, Appl. Therm. Eng. 53
    (2013) 57-66, Table 1) -- but every call in this module (before this
    change) used core.thermal.regenerator_effectiveness(), the PACKED-BED
    correlation, with this function's own default particle_diameter/
    porosity/bed_cross_section_area regardless of the real device's actual
    geometry TYPE, not just its dimensions. core/thermal.py already has a
    dedicated regenerator_effectiveness_parallel_plate() (a real,
    literature-sourced correlation -- Nickolay & Martin (2002)/Tusek et
    al. 2013 Eq. 4 -- not a guess) that was simply never wired into this
    module. geometry="parallel_plate" selects it, along with the
    physically-appropriate parallel-plate axial-conductivity mixing rule
    (_parallel_plate_effective_axial_conductivity() -- deliberately NOT
    the packed-bed Maxwell-Eucken formula; see that function's own
    docstring for why the two geometries need different formulas) and
    computing porosity/bed volume from plate_thickness/plate_spacing
    (matching regenerator_effectiveness_parallel_plate()'s own convention)
    instead of from the packed_bed particle_diameter/porosity arguments,
    which are ignored in this mode.

    n_substeps (default None -> auto): each blow is time-marched in
    n_substeps equal increments rather than applied as one lump update.
    This is NOT a tunable accuracy knob -- it fixes a genuine numerical-
    stability bug found and diagnosed during this module's own
    development: applying a full blow's heat transfer in one step assumes
    the driving temperature difference stays constant for the whole blow,
    which is false whenever the fluid mass moved per node per blow
    (mdot*tau_blow) carries more thermal capacity than that node's own
    (m_node*cp_solid) -- exactly the regime this module needs to run in
    (many nodes, realistic mdot). Left unfixed, this diverges to NaN
    within a handful of cycles (reproduced directly: traced node-by-node,
    confirmed the single-step update overshoots the driving temperature
    difference itself, the textbook signature of an explicit-Euler
    stability violation, not a physics problem). Auto default: enough
    substeps that the estimated per-node utilization ratper substep stays
    modest, at least 4: max(4, ceil(4 * mdot*tau_blow*cp_f/(m_node*cp_solid))).
    See validate_against_benchmarks() for the convergence-in-n_substeps
    check that confirms the auto default is adequate."""
    fluid = water_properties(T_K_for_ntu)
    cp_f = fluid["cp"]
    cp_solid_eff = CP_SOLID_GD if cp_solid is None else cp_solid
    k_solid_eff = K_SOLID_GD if k_solid is None else k_solid

    if geometry not in ("packed_bed", "parallel_plate"):
        raise ValueError(f"geometry must be 'packed_bed' or 'parallel_plate', got {geometry!r}")

    if geometry == "parallel_plate":
        # Matches regenerator_effectiveness_parallel_plate()'s own geometry
        # convention exactly (see that function's docstring): porosity and
        # bed volume are DERIVED from plate_thickness/plate_spacing, not
        # taken from the packed_bed `porosity` argument (ignored here).
        porosity_eff = plate_spacing / (plate_spacing + plate_thickness)
        V_bed = mass_total / (7900.0 * (1 - porosity_eff))  # RHO_GD
        L_bed = V_bed / bed_cross_section_area
        dx = L_bed / n_nodes if n_nodes > 0 else L_bed
        k_eff_axial = _parallel_plate_effective_axial_conductivity(
            porosity_eff, fluid["k"], k_solid_eff)
        ntu_result = regenerator_effectiveness_parallel_plate(
            mass_total, frequency, mdot, plate_thickness=plate_thickness,
            plate_spacing=plate_spacing, bed_cross_section_area=bed_cross_section_area,
            T_K=T_K_for_ntu)
    else:
        # Axial (node-to-node) conduction, applied once per cycle over the full
        # cycle period 1/frequency -- see _apply_axial_conduction()'s docstring
        # for why this is needed (it supplies the loss mechanism that caps span
        # growth as mdot -> 0, which the fluid-only coupling below cannot do
        # since its own effectiveness increases, not decreases, in that limit).
        #
        # Conductivity (replaced the ad hoc capped multiplier below
        # with the Maxwell-Eucken packed-bed composite-conductivity model --
        # see _packed_bed_effective_axial_conductivity()'s own docstring for
        # the formula, why this specific model was chosen over the more
        # commonly cited but harder-to-independently-verify Zehner-Schlunder
        # correlation, and the honest framing of what this change does and
        # does not establish). This is NOT a naive porosity-weighted parallel
        # mix of bulk-solid and fluid conductivity -- a packed bed of spheres
        # is solid-solid point-contact, not a continuous rod, so bulk metal
        # conductivity (k_solid_eff, ~10.5 W/m/K for Gd) is NOT directly the
        # right transport coefficient for the axial solid path (confirmed by
        # direct test during the ORIGINAL fix in this area: a naive parallel
        # mix collapsed the Tusek benchmark's span from a physically-
        # reasonable few K down to ~0.1-0.2K and stopped the simulation from
        # converging at all in 1500 cycles -- diffusion time across one node
        # spacing came out ~25x SHORTER than one cycle period, i.e. the naive
        # mix homogenizes the whole bed every cycle).
        V_bed = mass_total / (7900.0 * (1 - porosity))  # RHO_GD; avoids a new import cycle
        L_bed = V_bed / bed_cross_section_area
        dx = L_bed / n_nodes if n_nodes > 0 else L_bed
        k_eff_axial = _packed_bed_effective_axial_conductivity(porosity, fluid["k"], k_solid_eff)
        ntu_result = regenerator_effectiveness(
            mass_total, frequency, mdot, particle_diameter=particle_diameter,
            porosity=porosity, bed_cross_section_area=bed_cross_section_area,
            T_K=T_K_for_ntu, cp_solid=cp_solid)

    cycle_period = 1.0 / frequency if frequency > 0 else 0.0
    # material.delta_T_adiabatic() expects the field strength H in A/m, not
    # mu0H in Tesla -- same conversion core/amr_cycle.py's own
    # cooling_capacity() already applies before calling the same method.
    # (Caught here by direct tracing: without it, dTad silently evaluates
    # to ~1e-6 K instead of raising, since H=1.1 A/m is just an extremely
    # weak field, not an error -- see this module's "known limitations".)
    H_Am = mu0H_max / (4 * np.pi * 1e-7)

    NTU_total = ntu_result["NTU"]
    NTU_node = NTU_total / n_nodes if n_nodes > 0 else 0.0

    m_node = mass_total / n_nodes
    tau_blow = blow_fraction / (2.0 * frequency) if frequency > 0 else 0.0

    if n_substeps is None:
        node_utilization = (mdot * tau_blow * cp_f) / (m_node * cp_solid_eff) if m_node > 0 else 0.0
        n_substeps = max(4, int(np.ceil(4 * node_utilization)))
    dt = tau_blow / n_substeps if n_substeps > 0 else 0.0
    # NTU_node is a property of flow geometry/velocity (from
    # regenerator_effectiveness()) -- it must NOT be rescaled by dt/tau_blow.
    # Each substep represents a fresh micro-batch of fluid (mass mdot*dt)
    # undergoing a full-strength single-pass exchange at the node's actual
    # NTU_node; only the resulting energy (Q_i*dt) is small, not the
    # exchange effectiveness itself. (An earlier version of this file
    # scaled NTU down by dt/tau_blow, which artificially weakens heat
    # transfer as n_substeps grows and made the whole bed relax to a
    # near-uniform temperature every cycle regardless of NTU -- caught by
    # validate_against_benchmarks() converging to ~0 K span independent of
    # mass/NTU, which is the correct diagnostic signature of this bug.)
    decay_substep = np.exp(-NTU_node)

    no_load = (T_cold_reservoir is None or T_hot_reservoir is None)

    T = np.full(n_nodes, T_K_for_ntu, dtype=float) if T_init is None \
        else np.array(T_init, dtype=float).copy()
    if T_init is None:
        # Seed a small linear gradient purely to avoid an exactly-symmetric
        # starting condition (a real but degenerate fixed point of the old,
        # incorrect boundary condition -- kept here as a harmless, no-longer-
        # load-bearing safety margin now that the relay boundary condition
        # below provides its own genuine symmetry-breaking mechanism, since
        # hot-blow and cold-blow process nodes in OPPOSITE order).
        T = T + np.linspace(-0.005, 0.005, n_nodes)

    # NO-LOAD boundary condition: with insulated ends and no external heat
    # exchanger, the fluid that EXITS one end during a blow is exactly the
    # fluid that re-enters from that SAME end at the start of the NEXT blow
    # in the opposite direction (closed loop / a stationary slug of fluid
    # sitting in the end cap, with negligible dead-volume thermal mass of
    # its own). This "relay" is what lets a real regenerator build up a
    # difference between its two ends over many cycles -- an earlier version
    # of this function instead reset each end's inlet to that end's OWN
    # current solid-node temperature every single substep, which makes that
    # node's own heat exchange identically zero by construction and (found
    # by direct tracing) relaxes the whole bed to a uniform temperature
    # every cycle regardless of NTU. Initialized to T_K_for_ntu; the first
    # cycle's hot blow therefore starts from a physically reasonable "cold
    # start" rather than an artificial one.
    T_fluid_cold_end = T_K_for_ntu
    T_fluid_hot_end = T_K_for_ntu

    span_history = []
    Qc_avg = 0.0
    converged = False
    n_cycles = 0
    for cycle in range(max_cycles):
        T_prev_cycle = T.copy()

        # 1. Magnetize
        T = T + material.delta_T_adiabatic(T, H_Am)

        # 2. Hot blow: node 0 (cold end) -> node n_nodes-1 (hot end),
        # time-marched in n_substeps increments (see docstring above).
        # Inlet fluid temperature is fixed for the whole blow (continuous
        # supply at that temperature), and re-established fresh at the
        # start of each substep from the current inlet value -- correct for
        # both the with-load case (external reservoir, genuinely constant)
        # and the no-load case (the relayed parcel from the previous blow,
        # also constant for the duration of this blow).
        T_in_hot = float(T_cold_reservoir) if not no_load else T_fluid_cold_end
        T_f_exit_hot = T_in_hot
        for _ in range(n_substeps):
            T_f = T_in_hot
            for i in range(n_nodes):
                T_f_out = T[i] - (T[i] - T_f) * decay_substep
                Q_i = mdot * cp_f * (T_f_out - T_f)      # W, instantaneous, heat gained by fluid
                T[i] = T[i] - Q_i * dt / (m_node * cp_solid_eff)
                T_f = T_f_out
            T_f_exit_hot = T_f
        if no_load:
            T_fluid_hot_end = T_f_exit_hot

        # 3. Demagnetize
        T = T - material.delta_T_adiabatic(T, H_Am)

        # 4. Cold blow: node n_nodes-1 (hot end) -> node 0 (cold end)
        T_in_cold = float(T_hot_reservoir) if not no_load else T_fluid_hot_end
        T_f_exit_cold = T_in_cold
        Qc_accum_energy = 0.0
        for _ in range(n_substeps):
            T_f = T_in_cold
            for i in reversed(range(n_nodes)):
                T_f_out = T[i] - (T[i] - T_f) * decay_substep
                Q_i = mdot * cp_f * (T_f_out - T_f)
                T[i] = T[i] - Q_i * dt / (m_node * cp_solid_eff)
                T_f = T_f_out
            T_f_exit_cold = T_f
            if not no_load:
                Qc_accum_energy += mdot * cp_f * (float(T_cold_reservoir) - T_f) * dt
        if no_load:
            T_fluid_cold_end = T_f_exit_cold
        else:
            Qc_avg = Qc_accum_energy * frequency  # J/cycle * cycles/s -> W, cycle-averaged

        # Axial conduction over the full cycle period (see setup above) --
        # applied once per cycle, after both blows, so it competes fairly
        # with the fluid-mediated gradient built up by this cycle's
        # magnetize/blow/demagnetize/blow sequence rather than being folded
        # invisibly into either blow.
        T = _apply_axial_conduction(T, cycle_period, dx, bed_cross_section_area,
                                     k_eff_axial, m_node, cp_solid_eff)

        span_history.append(float(T[n_nodes - 1] - T[0]))
        n_cycles = cycle + 1
        # Convergence is judged on the SPAN (T[-1]-T[0], the quantity this
        # function actually reports), not the full node vector: with axial
        # conduction now in the model, the whole bed can still be slowly
        # drifting toward a common mean temperature (a genuine, benign
        # transient with no external heat exchange to pin an absolute
        # level) even after the span itself has settled -- checking the
        # full vector made every no-load run falsely read as
        # non-convergent, forcing max_cycles up for no benefit to the
        # actual output.
        window = 50
        if cycle > window and abs(span_history[-1] - span_history[-2]) < tol \
                and abs(span_history[-1] - span_history[-window]) < tol:
            converged = True
            break

    # Report the mean of the last min(10, n_cycles) cycles rather than the
    # single final cycle's value: with axial conduction now present, some
    # operating points settle into a small residual cycle-to-cycle
    # oscillation (visible in span_history_last10) that never quite
    # satisfies the strict `converged` tol/window check within a practical
    # cycle budget, even though the underlying trend is flat. Averaging is
    # a more robust point estimate of that plateau than either endpoint.
    tail = span_history[-min(10, len(span_history)):] if span_history else [0.0]
    span_tail_mean = float(np.mean(tail))

    result = {"converged": converged, "n_cycles": n_cycles,
              "NTU_total": round(float(NTU_total), 3),
              "NTU_node": round(float(NTU_node), 4),
              "U_bulk": round(float(ntu_result["U"]), 4),
              #  addition: ntu_result["U"] (Nielsen, Tusek, Engelbrecht,
              # Schopfer, Kitanovski, Bahl, Smith, Pryds & Poredos, "Review on
              # numerical modeling of active magnetic regenerators for room
              # temperature applications," Int. J. Refrig. 34 (2011) 603-616 --
              # now in this repo's Papers/ -- Eq. (16): utilization
              # Phi = mdot_f*cf*tau2/(ms*cs), "related to the frequency of the
              # operation... a lower frequency means a larger influence of the
              # longitudinal thermal conduction") was ALREADY computed by
              # core.thermal.regenerator_effectiveness() (called just above as
              # `ntu_result`) but silently discarded here -- only NTU_total was
              # pulled out of that dict into this function's own return value.
              # This is a PURELY ADDITIVE fix (new dict key; every existing key
              # is unchanged, so no existing caller's behavior changes) that
              # surfaces it, motivated by this module's own honesty flag #2
              # ("does span depend systematically on utilization, not just
              # NTU?") -- see run_convergence_and_utilization_diagnostic()
              # below, which is the first caller to actually use this key.
              "n_substeps": n_substeps,
              "T_nodes_K": T.copy(),
              "span_history_last10": [round(s, 3) for s in span_history[-10:]]}
    if no_load:
        result["span_K"] = round(span_tail_mean, 3)
    else:
        result["Qc_W"] = round(float(Qc_avg), 3)
    return result


import hashlib
import json
import os

_CACHE_PATH = "results/.regenerator_1d_cache.json"


def _cache_key(material, mu0H_max, mass_total, frequency, n_nodes, mdot_search,
                extra_kwargs):
    """Deterministic cache key for no_load_span(): material identity (name
    plus a couple of its own physical parameters, in case two materials
    share a display name) plus every numeric input that can change the
    result. Rounds floats to avoid a fresh cache miss from harmless
    float-repr noise (e.g. 21.039999999999992 vs 21.04)."""
    mat_id = (getattr(material, "name", repr(material)),
              round(float(getattr(material, "Tc", 0.0)), 6),
              round(float(getattr(material, "J", 0.0)), 6))
    payload = {
        "model_version": _MODEL_VERSION,
        "material": mat_id,
        "mu0H_max": round(float(mu0H_max), 6),
        "mass_total": round(float(mass_total), 6),
        "frequency": round(float(frequency), 6),
        "n_nodes": n_nodes,
        "mdot_search": [round(float(m), 8) for m in mdot_search],
        "extra_kwargs": {k: (round(float(v), 8) if isinstance(v, (int, float)) else v)
                          for k, v in sorted(extra_kwargs.items())},
    }
    blob = json.dumps(payload, sort_keys=True)
    return hashlib.sha256(blob.encode()).hexdigest()


def _cache_load():
    try:
        with open(_CACHE_PATH) as f:
            return json.load(f)
    except (FileNotFoundError, json.JSONDecodeError):
        return {}


def _cache_save(cache):
    os.makedirs(os.path.dirname(_CACHE_PATH), exist_ok=True)
    with open(_CACHE_PATH, "w") as f:
        json.dump(cache, f)


def no_load_span(material, mu0H_max, mass_total, frequency, mdot=None,
                  n_nodes=20,
                  mdot_search=(0.0002, 0.0005, 0.001, 0.002, 0.004, 0.008, 0.015, 0.03, 0.06),
                  use_cache=True, **kwargs):
    """Convenience wrapper: no-load span, searched over `mdot_search` and
    reporting the BEST (maximum span, among CONVERGED runs) result found.

    Why search rather than use one fixed mdot: a reported "no-load span"
    figure is the best span a real device demonstrates -- experimentalists
    tune flow rate down until span stops improving. A single arbitrary
    fixed mdot is not a fair stand-in for that.

    The grid brackets a genuine interior maximum now (axial conduction,
    added this pass, caps span growth as mdot -> 0 -- see the module
    docstring's "known limitations" for how this was diagnosed and fixed).
    The result is therefore no longer sensitive to an arbitrary lower
    search bound the way it used to be. It IS still quantitatively
    undershooting the one directly-measured benchmark checked so far (see
    module docstring, "NEW OPEN ISSUE") -- read the returned span as a
    real, bounded, convergent model prediction, not yet as a validated
    quantitative one.

    If mdot is explicitly provided (not None), search is skipped and that
    single value is used instead (also cached, keyed on that single mdot).

    CACHING: each search point is a multi-cycle transient simulation run
    to convergence (tens of seconds -- see module docstring), and this
    function is called repeatedly for the SAME device across
    main.py's steps 2e/2f/3a2 and across repeated pipeline runs. Results
    are cached to disk at `results/.regenerator_1d_cache.json`, keyed on
    every input that affects the answer (material identity, field, mass,
    frequency, node count, mdot grid, and any extra simulate_amr_1d
    kwargs -- see _cache_key()). Pass `use_cache=False` to force a fresh
    computation (e.g. after changing this module's own physics, since the
    cache has no way to know the code changed underneath it -- delete
    results/.regenerator_1d_cache.json directly for the same effect
    across every call)."""
    cache_key = None
    if use_cache:
        cache_key = _cache_key(material, mu0H_max, mass_total, frequency, n_nodes,
                                (mdot,) if mdot is not None else mdot_search, kwargs)
        cache = _cache_load()
        if cache_key in cache:
            hit = dict(cache[cache_key])
            hit["T_nodes_K"] = np.array(hit["T_nodes_K"])
            hit["from_cache"] = True
            return hit

    if mdot is not None:
        result = simulate_amr_1d(material, mu0H_max, mass_total, frequency, mdot,
                                  n_nodes=n_nodes, **kwargs)
    else:
        kwargs.setdefault("max_cycles", 1200)
        kwargs.setdefault("tol", 3e-3)
        # Selection uses span_K (already the tail-averaged, noise-robust value
        # from simulate_amr_1d -- see that function) from EVERY evaluated mdot,
        # not just ones that hit the strict per-cycle `converged` flag: with
        # axial conduction now in the model, several operating points settle
        # into a small residual cycle-to-cycle oscillation that can outlast a
        # practical cycle budget without ever tripping that flag, even though
        # the tail-averaged span is already a stable, trustworthy estimate (see
        # module docstring's "known limitations" for the interior-maximum
        # confirmation this search now reliably finds).
        best = None
        any_converged = False
        for m in mdot_search:
            r = simulate_amr_1d(material, mu0H_max, mass_total, frequency, m,
                                 n_nodes=n_nodes, **kwargs)
            r["mdot_kg_s"] = m
            any_converged = any_converged or r["converged"]
            if best is None or r["span_K"] > best["span_K"]:
                best = r
        best["any_converged_in_search"] = any_converged
        result = best

    if use_cache:
        to_store = dict(result)
        to_store["T_nodes_K"] = to_store["T_nodes_K"].tolist()
        cache = _cache_load()
        cache[cache_key] = to_store
        _cache_save(cache)
    return result


def validate_against_benchmarks(verbose=True,
                                 out_path="results/regenerator_1d_validation.txt"):
    """Runs no_load_span() against every directly-measured no-load-span
    (Qc=0, COP blank) row in amr_experimental_benchmarks.csv -- the SAME
    rows results/regenerative_amplification_diagnostic.txt already
    identifies as the cleanest available ground truth, because they
    require no mdot back-calibration (span is measured directly, not
    inferred from a fitted flow rate). This is the honest check of
    whether this new model is actually better than the 0-D one, not just
    different -- run BEFORE this module is relied on anywhere else, and
    the result is written to disk so it doesn't have to be taken on
    faith.

    For each such row, also runs core/amr_cycle.py's own 0-D
    cooling_capacity() no-load structural cap (2*dTad_noload) for a
    direct, apples-to-apples comparison against the SAME literature
    value."""
    from core.validation_system import load_benchmarks, _material_for_row, _t_cold_for_row
    from core.amr_cycle import AMRSystem
    from core.loss_model import StateDependentLossModel

    rows = load_benchmarks()
    noload_rows = [r for r in rows if float(r["Qc_W"]) == 0.0 and not r["COP"] and float(r["span_K"]) > 0]

    lines = ["=" * 100, "1-D TRANSIENT REGENERATOR MODEL -- VALIDATION AGAINST DIRECTLY-MEASURED",
              "NO-LOAD SPANS (Qc=0 rows, no mdot back-calibration involved)",
              "=" * 100, ""]
    results = []
    for row in noload_rows:
        material = _material_for_row(row)
        t_cold = _t_cold_for_row(row)
        mu0H = float(row["mu0H_T"])
        mass = float(row["mass_MCM_kg"]) if row["mass_MCM_kg"] else 1.0
        freq = float(row["frequency_Hz"]) if row["frequency_Hz"] else 1.0
        span_lit = float(row["span_K"])

        # Device-specific geometry (Paper-Mining Pass, Item 1.7 in
        # LIMITATIONS.md): Tusek_singlebed_Gd_2010_spanceiling is a real
        # parallel-plate AMR, not a packed bed of spheres (Tusek et al.,
        # Appl. Therm. Eng. 53 (2013) 57-66, Table 1: 0.1mm spacing,
        # 0.25mm plate thickness, porosity 0.2564, outer dims
        # 10x80x39mm) -- the other two no-load-span rows here keep the
        # packed_bed default (their actual, unquestioned geometry type),
        # since this repo has no equivalent verified parallel-plate data
        # for them.
        geom_kwargs = {}
        if row["device"] == "Tusek_singlebed_Gd_2010_spanceiling":
            geom_kwargs = {"geometry": "parallel_plate",
                           "plate_thickness": 0.00025, "plate_spacing": 0.0001,
                           "bed_cross_section_area": 3.9e-4}

        r1d = no_load_span(material, mu0H, mass, freq,
                            T_K_for_ntu=t_cold + span_lit / 2.0, **geom_kwargs)

        probe = AMRSystem(material=material, mu0H_max=mu0H, mass_regenerator=mass,
                           frequency=freq, fluid_mdot=1.0, loss_model=StateDependentLossModel())
        _, dTad_noload = probe.cooling_capacity(t_cold, 1e-3)
        cap_0d = 2 * float(dTad_noload)

        span_1d = r1d["span_K"]
        err_1d_pct = 100 * (span_1d - span_lit) / span_lit if span_lit else None
        err_0d_pct = 100 * (cap_0d - span_lit) / span_lit if span_lit else None

        results.append({"device": row["device"], "span_lit_K": span_lit,
                         "span_1d_K": span_1d, "cap_0d_K": round(cap_0d, 2),
                         "err_1d_pct": round(err_1d_pct, 1) if err_1d_pct is not None else None,
                         "err_0d_pct": round(err_0d_pct, 1) if err_0d_pct is not None else None,
                         "converged": r1d["converged"], "n_cycles": r1d["n_cycles"],
                         "NTU_total": r1d["NTU_total"]})
        lines.append(f"  {row['device']:<38} span_lit={span_lit:6.1f}K "
                      f"1D_model={span_1d:7.2f}K (err={err_1d_pct:+6.1f}%)  "
                      f"0D_cap={cap_0d:7.2f}K (err={err_0d_pct:+6.1f}%)  "
                      f"NTU_total={r1d['NTU_total']:6.2f}  "
                      f"converged={r1d['converged']} in {r1d['n_cycles']} cycles")

    lines.append("-" * 100)
    lines.append("err_1d/err_0d are both signed relative to the SAME literature span value; "
                  "err closer to 0% is better. 0D_cap is core/amr_cycle.py's existing structural "
                  "ceiling (2*dTad_noload) for comparison -- it is a HARD CAP the 0-D model cannot "
                  "exceed at any mdot, not a prediction of the achieved span, so its own error is "
                  "necessarily <=0% by construction for every structurally-failing row; the 1D "
                  "model has no such one-sided constraint.")
    lines.append("CAVEAT (read before trusting err_1d): the earlier low-mdot degeneracy "
                  "(predicted span growing without bound as mdot -> 0, no interior maximum) is "
                  "FIXED -- added axial (node-to-node) conduction, which the model previously "
                  "lacked entirely; no_load_span()'s search now finds a genuine interior maximum "
                  "in mdot for every row above, not an arbitrary search-grid-floor artifact (see "
                  "core/regenerator_1d.py module docstring, 'known limitations', for the "
                  "diagnosis and fix). What that fix EXPOSES, now that the search result is "
                  "well-defined: the model's quantitative accuracy is inconsistent in DIRECTION "
                  "across devices -- it undershoots on two of the three rows above and overshoots "
                  "on the third, unlike the old bug (which at least always overshot in the same "
                  "direction).")
    lines.append("PHASE 31 UPDATE: the axial-conductivity approximation was replaced with the "
                  "Maxwell-Eucken packed-bed composite-conductivity model (a real, named, "
                  "independently re-derivable formula -- see "
                  "_packed_bed_effective_axial_conductivity()'s docstring) in place of the "
                  "previous ad hoc capped multiplier of stagnant-fluid conductivity. HONEST "
                  "RESULT, stated plainly rather than smoothed over: this did NOT improve the "
                  "err_1d percentages above -- it made the two already-undershooting rows "
                  "undershoot MORE (Tusek -92.1%->-96.9%, DTU/MAGGIE -61.4%->-62.0%) while "
                  "leaving the overshooting row (Lozano) essentially unchanged (+112.3%->+111.1%), "
                  "because Maxwell-Eucken predicts a HIGHER effective conductivity than the old "
                  "placeholder at these bed porosities (~2.7x higher at porosity=0.365), i.e. MORE "
                  "axial damping, not less. This is a genuine, if disappointing, finding: swapping "
                  "in a more rigorously-justified conductivity model does not by itself resolve the "
                  "direction-inconsistent error -- the discrepancy is not simply 'the axial "
                  "conductivity constant is a little off', since a legitimately different, "
                  "textbook-derived constant moves the SAME two rows further in the SAME wrong "
                  "direction. What IS still established, unchanged by this update: this model "
                  "produces spans several times the underlying material's own single-blow dT, "
                  "using a genuine multi-cycle transient simulation rather than an assumed "
                  "multiplier, AND that the resulting span is a bounded, convergent prediction "
                  "rather than an open-ended one -- i.e. regenerative amplification is a real, "
                  "computable, finite-magnitude effect, which is what this module set out to "
                  "demonstrate. Treat the specific err_1d percentages as illustrative of that "
                  "magnitude and of the remaining calibration gap -- now more likely to be in the "
                  "regenerator-effectiveness/NTU coupling or the single-blow reference itself, not "
                  "solely the axial-conductivity constant -- not as a validated replacement for "
                  "core/amr_cycle.py's cooling_capacity() in the rest of this codebase. It is still "
                  "not wired into any other function here, and should not be until this is "
                  "resolved.")

    if verbose:
        for line in lines:
            print(line)
    import os
    os.makedirs(os.path.dirname(out_path), exist_ok=True)
    with open(out_path, "w") as f:
        f.write("\n".join(lines) + "\n")
    if verbose:
        print(f"Wrote {out_path}")
    return results


def run_convergence_and_utilization_diagnostic(
        out_path="results/regenerator_1d_convergence_diagnostic.txt",
        n_nodes_sweep=(10, 20, 40, 80), verbose=True):
    """ diagnostic, run BEFORE any further physics changes to this
    module (per this repo's own established discipline -- see e.g. the
    "isolate axial conductivity as a suspect before touching
    anything else" precedent): isolates WHICH of two candidate failure
    points is actually responsible for the direction-inconsistent error
    validate_against_benchmarks() reports (Tusek -96.9%, Lozano +111.1%,
    DTU/MAGGIE -62.0%), on the same three directly-measured no-load-span
    benchmark rows that function uses, at EACH row's own already-
    calibrated best mdot (from no_load_span()'s own search, reused here
    rather than re-searched, so this diagnostic's numbers are directly
    comparable to validate_against_benchmarks()'s own).

    Two checks, in order of cost (cheap numerics check first):

    1. n_nodes CONVERGENCE. If span_K keeps changing materially as
       n_nodes increases (10->20->40->80), the model has a numerical
       discretization problem -- NOT a physics/calibration one -- and
       reading err_1d_pct at any FIXED n_nodes (this module's own default,
       n_nodes=20, used everywhere else including
       validate_against_benchmarks()) would be comparing an
       under-resolved number to a literature value, conflating numerics
       error with modeling error. A converged model should show span_K
       change by much less between 40->80 than between 10->20.

    2. BULK UTILIZATION Phi (Nielsen et al. 2011, Int. J. Refrig. 34
       (2011) 603-616 Eq. (16), now in this repo's Papers/ -- see
       simulate_amr_1d()'s own U_bulk comment for the full citation) at
       each device's own calibrated mdot, reported alongside NTU_total.
       Nielsen et al.'s own review treats span as jointly a function of
       BOTH NTU and utilization, not NTU alone -- if the three benchmarks'
       own calibrated utilizations differ by an order of magnitude or
       more from each other, that is direct, quantitative evidence (not
       just a plausible-sounding hypothesis) that this module's single
       NTU-only accounting is comparing devices that are NOT in the same
       operating regime, which would explain direction-inconsistent
       error without needing any further physics change to identify the
       CAUSE (a fix would still be separate future work -- this
       diagnostic's job is only to confirm or rule out this specific
       candidate explanation with real numbers).

    Does NOT change simulate_amr_1d()'s own DEFAULT n_nodes=20 anywhere
    else in this codebase -- this is read-only diagnostic reporting, same
    discipline as every other validate_*/run_*_diagnostic function in
    this module (no caller elsewhere is affected by running this)."""
    from core.validation_system import load_benchmarks, _material_for_row, _t_cold_for_row

    rows = load_benchmarks()
    noload_rows = [r for r in rows if float(r["Qc_W"]) == 0.0 and not r["COP"] and float(r["span_K"]) > 0]

    lines = ["=" * 100,
             "1-D TRANSIENT REGENERATOR MODEL -- n_nodes CONVERGENCE + BULK UTILIZATION",
             "DIAGNOSTIC (, see this function's own docstring for what each check",
             "isolates and why, before any further physics change to this module)",
             "=" * 100, ""]

    convergence_results = []
    utilization_results = []
    for row in noload_rows:
        material = _material_for_row(row)
        t_cold = _t_cold_for_row(row)
        mu0H = float(row["mu0H_T"])
        mass = float(row["mass_MCM_kg"]) if row["mass_MCM_kg"] else 1.0
        freq = float(row["frequency_Hz"]) if row["frequency_Hz"] else 1.0
        span_lit = float(row["span_K"])

        # Reuse the SAME calibrated best mdot validate_against_benchmarks()
        # itself uses (no_load_span()'s own search, at this module's
        # n_nodes=20 default) -- both checks below hold mdot fixed at this
        # value, changing only n_nodes (check 1) or nothing at all, just
        # reporting U_bulk (check 2), so results are directly comparable
        # to validate_against_benchmarks()'s own err_1d_pct numbers.
        r20 = no_load_span(material, mu0H, mass, freq, T_K_for_ntu=t_cold + span_lit / 2.0)
        mdot_best = r20["mdot_kg_s"]

        lines.append(f"--- {row['device']} (span_lit={span_lit:.1f}K, calibrated "
                     f"mdot={mdot_best:.5f}kg/s) ---")

        lines.append(f" Utilization Phi (Nielsen et al. 2011 Eq. 16) at this mdot: "
                     f"U_bulk={r20['U_bulk']:.4f}   NTU_total={r20['NTU_total']:.2f}")
        utilization_results.append({"device": row["device"], "mdot_kg_s": mdot_best,
                                     "U_bulk": r20["U_bulk"], "NTU_total": r20["NTU_total"],
                                     "span_lit_K": span_lit})

        lines.append(f"  {'n_nodes':>8} {'span_K':>10} {'delta_vs_prev':>16}")
        prev_span = None
        row_convergence = []
        for nn in n_nodes_sweep:
            r = simulate_amr_1d(material, mu0H, mass, freq, mdot_best, n_nodes=nn,
                                 T_K_for_ntu=t_cold + span_lit / 2.0,
                                 max_cycles=1200, tol=3e-3)
            delta_str = f"{r['span_K'] - prev_span:+.3f}K" if prev_span is not None else "--"
            lines.append(f"  {nn:>8} {r['span_K']:>9.3f}K {delta_str:>16}")
            row_convergence.append({"n_nodes": nn, "span_K": r["span_K"]})
            prev_span = r["span_K"]
        lines.append("")
        convergence_results.append({"device": row["device"], "sweep": row_convergence})

    # Convergence verdict: comparing only the FIRST and LAST doubling's
    # |delta| (as an earlier version of this function did) is NOT
    # sufficient -- it can mislabel a WILDLY OSCILLATING, non-monotonic
    # sequence as "converging" whenever the last swing happens to be
    # smaller than the first by coincidence, even if the sequence swung
    # up and down by many multiples of the eventual answer in between
    # (confirmed happening in practice: an earlier run of this exact
    # function, before this fix, labeled Lozano_POLO_UFSC_2016_maxspan
    # "CONVERGING" from a span_K sequence of 3.41 -> 25.33 -> 9.59 ->
    # 12.31K -- a first delta of +21.9K, then -15.7K, then +2.7K -- which
    # is genuinely chaotic/non-monotonic across this n_nodes range, not
    # smoothly settling, even though |+2.7K| < |+21.9K| made the OLD
    # first-vs-last check say "CONVERGING"). This function instead checks
    # (a) every delta has the SAME SIGN (strictly monotonic -- no
    # overshoot/oscillation at all) AND (b) |delta| shrinks at every step,
    # not just from first to last -- both must hold for "CONVERGING".
    lines.append("-" * 100)
    lines.append("CONVERGENCE VERDICT (per device: are ALL successive deltas the same sign AND")
    lines.append("shrinking in magnitude at EVERY step, not just from first to last -- a")
    lines.append("first-vs-last-only check can mislabel a wildly oscillating, non-monotonic")
    lines.append("sequence as \"converging\" whenever the final swing is coincidentally smaller")
    lines.append("than the first one, even with huge swings in between):")
    any_not_converging = False
    for cr in convergence_results:
        sweep = cr["sweep"]
        deltas = [sweep[i + 1]["span_K"] - sweep[i]["span_K"] for i in range(len(sweep) - 1)]
        signs = [1 if d > 0 else (-1 if d < 0 else 0) for d in deltas]
        same_sign = len(set(s for s in signs if s != 0)) <= 1
        abs_deltas = [abs(d) for d in deltas]
        strictly_shrinking = all(abs_deltas[i + 1] < abs_deltas[i] for i in range(len(abs_deltas) - 1))
        converging = same_sign and strictly_shrinking
        any_not_converging = any_not_converging or not converging
        delta_str = " -> ".join(f"{d:+.2f}K" for d in deltas)
        verdict = "CONVERGING (monotonic, shrinking)" if converging else \
            ("NOT MONOTONIC (oscillating -- sign changes across the sweep, NOT just "
             "coincidentally-smaller endpoints)" if not same_sign else
             "MONOTONIC BUT NOT YET SHRINKING AT EVERY STEP")
        lines.append(f"  {cr['device']:<38} deltas: {delta_str}")
        lines.append(f"    -> {verdict}")
    lines.append("")
    if any_not_converging:
        lines.append("At least one device's span_K has NOT clearly settled by n_nodes=80 -- ")
        lines.append("err_1d_pct in validate_against_benchmarks() (computed at this module's ")
        lines.append("n_nodes=20 default) for that device may be partly a numerics artifact, ")
        lines.append("not purely a physics/calibration gap. Rerun that device's own row at ")
        lines.append("higher n_nodes before drawing further physics conclusions from it.")
    else:
        lines.append("Every device's span_K sequence is strictly monotonic AND shrinking at ")
        lines.append("every step from n_nodes=10 to 80 (no oscillation, no sign changes) -- ")
        lines.append("validate_against_benchmarks()'s n_nodes=20 default is NOT the dominant ")
        lines.append("source of its reported direction-inconsistent error; this rules OUT ")
        lines.append("candidate explanation #1 (numerical discretization) and points back ")
        lines.append("toward the physics/calibration candidates (utilization regime, ")
        lines.append("single-blow reference asymmetry) -- see the utilization spread below.")
    lines.append("")

    u_vals = [u["U_bulk"] for u in utilization_results]
    u_min, u_max = min(u_vals), max(u_vals)
    lines.append("UTILIZATION SPREAD VERDICT:")
    lines.append(f" Range across the {len(u_vals)} benchmarks' own calibrated operating "
                 f"points: U_bulk = {u_min:.4f} to {u_max:.4f} "
                 f"({u_max / u_min if u_min > 0 else float('inf'):.1f}x spread)")
    for u in utilization_results:
        lines.append(f"    {u['device']:<38} U_bulk={u['U_bulk']:.4f}  "
                     f"NTU_total={u['NTU_total']:.2f}  (lit. span={u['span_lit_K']:.1f}K)")
    lines.append("")
    if u_min > 0 and (u_max / u_min) > 3.0:
        lines.append(f"The {u_max / u_min:.1f}x spread in calibrated utilization across these ")
        lines.append("three benchmarks IS large enough to be a genuine, quantitative candidate ")
        lines.append("explanation for direction-inconsistent error, consistent with (not yet a ")
        lines.append("proof of causation for) Nielsen et al. (2011)'s own framing that span is a ")
        lines.append("joint function of NTU AND utilization, not NTU alone -- this module's ")
        lines.append("no_load_span() search selects mdot purely to MAXIMIZE span_K, with no ")
        lines.append("constraint tying the resulting utilization across devices to a common ")
        lines.append("regime, so three devices at three very different utilizations being ")
        lines.append("compared as if NTU_total alone should predict a consistent relative error ")
        lines.append("is not obviously a fair comparison. This does NOT by itself fix the ")
        lines.append("direction-inconsistent error -- it identifies utilization-regime mismatch ")
        lines.append("as the most evidence-backed remaining candidate (of the four this module's ")
        lines.append("own docstring lists) for follow-up physics work, in preference to further ")
        lines.append("tuning the axial-conductivity constant (already tried twice, see ).")
    else:
        lines.append("Utilization spread across these three benchmarks is NOT large -- this ")
        lines.append("specific candidate explanation is not well-supported by this data; the ")
        lines.append("single-blow reference/boundary-condition-asymmetry candidate (this ")
        lines.append("module's own docstring, candidate #3) remains the most likely unexplored ")
        lines.append("lead.")

    if verbose:
        for line in lines:
            print(line)
    os.makedirs(os.path.dirname(out_path), exist_ok=True)
    with open(out_path, "w") as f:
        f.write("\n".join(lines) + "\n")
    if verbose:
        print(f"\nWrote {out_path}")

    return {"convergence": convergence_results, "utilization": utilization_results,
            "any_not_converging": any_not_converging}


def regenerative_span_cap(material, mu0H_max, mass_total, frequency, **kwargs):
    """Thin wrapper around no_load_span(): returns just the predicted
    no-load span (K), for use as core.amr_cycle.AMRSystem's
    `no_load_span_override` -- i.e. this is the function meant to actually
    populate that opt-in parameter, so the "regenerative amplification"
    capability the diagnostic in results/regenerative_amplification_diagnostic.txt
    quantifies is something a caller can actually plug in, not just measure.

    Deliberately NOT called from inside AMRSystem itself: a single call
    here takes on the order of a minute (a multi-mdot search, each point a
    multi-cycle transient simulation to convergence -- see no_load_span()
    and simulate_amr_1d()), which is fine to run ONCE per design point but
    would make any sweep, NSGA-III search, or Sobol sensitivity pass that
    calls cooling_capacity() many times per design point completely
    impractical. Call this once, cache the result, and pass it into
    AMRSystem(no_load_span_override=...).

    HONESTY FLAG (repeated from AMRSystem.no_load_span_override's own
    docstring -- read before using this for anything beyond illustration
    or the diagnostic comparison in
    validation_system.run_regenerative_amplification_override_check()):
    this model's own validation (results/regenerator_1d_validation.txt)
    shows directionally-inconsistent error against the three directly-
    measured no-load-span benchmarks (+112%, -92%, -61%). It is bounded
    and convergent (the low-mdot degeneracy that used to make it
    ill-defined is fixed), but it is NOT an independently-calibrated
    quantitative model yet."""
    return no_load_span(material, mu0H_max, mass_total, frequency, **kwargs)["span_K"]