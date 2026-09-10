"""
plots.py — Publication-quality figures for magcool-dc
========================================================
Generates every figure supported by the magcool-dc physics/economics/
optimization suite: magnetocaloric material validation, first-order
(Landau) giant-MCE modelling, AMR cycle characteristic curves, thermal /
geometry regenerator modelling, loss-model calibration, system-level
validation against published prototypes, cascade and Curie-graded
cascade staging, multi-objective (NSGA-III) design optimization, Sobol
global sensitivity analysis, RSM surrogate fitting, economics/TCO, and
GWP/emissions comparisons.

Run:    python plots.py
Output: results/figures/  (PNG + PDF, one pair per figure)

Figures 27-34 cover analyses added in earlier rounds that previously had no
figure: Tc-broadening (27), nanocomposite off-design robustness (28),
thermal-diode actuation cost (29), magnetocaloric-fluid volume-fraction
sweep (30), passive-regenerator alignment effect (31), rotary-device
cycle-type validation (32), hysteresis-loss Pareto-front sensitivity (33),
and Halbach-cylinder magnet-mass Pareto-front sensitivity (34). Figures 33
and 34 re-run NSGA-III twice each (same as fig18) and fall back to
pre-computed results/pareto_front_*.csv files if pymoo is unavailable.

Figures 36-45 cover later additions: the MCE_Value_Proposition.md figure
set (36-39: COP crossover, emissions breakdown, refrigerant GWP, and
real-deployment comparisons), speculative MnFePSi-doped hysteresis
modelling and its k-fit quality (40-41), a hysteresis-exploiting actuator
cost estimate and the underlying hysteresis-loss landscape (42-43), and
the Phase 37 calibration-drift fix together with the Lozano row-by-row
calibration-status check it motivated (44-45).

Figures 46-56 cover the final rounds of previously-computed-but-unplotted
analyses: cross-technology caloric comparison vs. this repo's own VCC COP
(46), annual water-usage/WUE comparison (47), a Monte Carlo calibration-
uncertainty band on COP_electrical (48), NSGA-III Pareto-front
seed-to-seed stability (49), a two-search regime-crossover null result
(50), the Ames Lab heat-pump architecture check (51), Hypereg parallel-
hydraulic pumping-power sensitivity (52), the corpus-wide regenerative-
amplification gap across the full benchmark set (53), the hybrid
solid-state regenerator's ideal-to-real COP loss funnel vs. VCC (54),
whether the single-design-point COP conclusion survives a full annual
climate profile (55), and the heat-transfer fluid selection trade-off
behind core/fluids.py's own DEFAULT_FLUID choice (56).

Figures 57-61 close out the remaining computed-but-unplotted analyses:
a commercial-landscape reality check against real vendor-claimed
magnetocaloric products (57), the measured system-level cost of not
using the Gd grain-Tc-broadening physics fix as the default material
(58), the 1-D transient regenerator model's own validation against
directly-measured no-load spans (59), a field-dependent Tc-broadening
fit checked against held-out data (60), and the multiseed robustness
check on fig34's single-seed magnet-geometry Pareto finding (61).

Notes
-----
Two analyses in this repository (Sobol sensitivity via SALib, NSGA-III
optimization via pymoo) depend on optional third-party packages. If
those packages are not installed, this script transparently falls back
to the pre-computed results already checked into results/ (sobol_results
*.txt, pareto_front.csv) so that `python plots.py` always produces a
complete figure set regardless of which optional dependencies are
available. Every other figure is computed fresh from the physics/
economics models in core/.
"""

import os
import sys
import csv
import re
import textwrap
from pathlib import Path

import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

# Make sure `core.*` imports resolve regardless of the working directory
# this script is invoked from.
_PROJECT_ROOT = os.path.dirname(os.path.abspath(__file__))
if _PROJECT_ROOT not in sys.path:
    sys.path.insert(0, _PROJECT_ROOT)

from core.mce_material import GADOLINIUM
from core.first_order_mce import GD5SI2GE2_FIRST_ORDER, LAFESIH_FIRST_ORDER
from core.amr_cycle import AMRSystem
from core.baseline_cooling import (vapor_compression_cop, liquid_cooling_cop,
                                    elastocaloric_reference_cop)
from core.thermal import regenerator_effectiveness
from core import validation
from core import validation_system
from core import loss_model as loss_model_mod
from core.loss_model import (StateDependentLossModel, CALIBRATION_POINTS_CORE,
                              CALIBRATION_POINTS_EXTENDED,
                              CALIBRATION_POINTS_FURTHER_EXTENDED)
from core import economics
from core import emissions
from core import cascade
from core import giant_mce_analysis
from core import material_family_comparison
from core import giguere_validation
from core import geometry_analysis
from core import rsm as rsm_mod
from core import inhomogeneous_broadening
from core import nanocomposite_material
from core import thermal_diode_analysis
from core import fluid_mce_analysis
from core import passive_regenerator_analysis
from core import hysteresis_sensitivity
from core import magnet_geometry
from core import beverage_cooler_validation
from core import alternative_caloric_comparison
from core import water_usage
from core import uncertainty_propagation
from core import heat_pump_validation
from core import hypereg_analysis
from core import regime_crossover_analysis
from core import hybrid_solid_state_regenerator
from core import pue_annualized
from core import fluid_selection_optimization
from core import commercial_landscape
from core import regenerator_1d

try:
    from core import sensitivity as sensitivity_mod
    HAVE_SALIB = True
except ImportError:
    HAVE_SALIB = False

try:
    from core import optimize as optimize_mod
    HAVE_PYMOO = True
except ImportError:
    HAVE_PYMOO = False

# ─── Output directories ────────────────────────────────────────────────────
RESULTS_DIR = Path('results')
FIG_DIR = RESULTS_DIR / 'figures'
FIG_DIR.mkdir(parents=True, exist_ok=True)
RESULTS_DIR.mkdir(parents=True, exist_ok=True)

mu0 = 4 * np.pi * 1e-7

# ─── Global style ──────────────────────────────────────────────────────────
plt.rcParams.update({
    'font.family'      : 'DejaVu Sans',
    'font.size'        : 11,
    'axes.labelsize'   : 12,
    'axes.titlesize'   : 13,
    'legend.fontsize'  : 10,
    'figure.dpi'       : 150,
    'lines.linewidth'  : 2.0,
    'axes.grid'        : True,
    'grid.alpha'       : 0.35,
    'axes.spines.top'  : False,
    'axes.spines.right': False,
})

COLOR_MAIN   = '#1f4e79'
COLOR_POWER  = '#c00000'
COLOR_CYCLE  = plt.cm.viridis(np.linspace(0.15, 0.85, 4))
COLOR_CYCLE8 = plt.cm.viridis(np.linspace(0.15, 0.85, 8))


def _atomic_savefig(fig, final_path, **kwargs):
    """Saves to a temp file in the same directory, then atomically renames
    into place. Prevents a mid-write OS-level failure (e.g. an Errno 22
    seen intermittently on Windows, likely AV/file-lock interference) from
    leaving a truncated/corrupted figure file at final_path -- previously,
    fig05's .pdf pass failed partway through and left a malformed PDF
    (missing 'trailer <<...>> startxref' in its final bytes) on disk even
    though the pipeline correctly logged the stage as failed."""
    tmp_path = final_path.with_name(final_path.stem + '.tmp' + final_path.suffix)
    try:
        fig.savefig(tmp_path, **kwargs)
        os.replace(tmp_path, final_path)
    finally:
        if tmp_path.exists():
            tmp_path.unlink()


def save(fig, name):
    path = FIG_DIR / name
    _atomic_savefig(fig, path.with_suffix('.png'), bbox_inches='tight', dpi=150)
    _atomic_savefig(fig, path.with_suffix('.pdf'), bbox_inches='tight')
    plt.close(fig)
    print(f" Saved -> {path}.png / .pdf")


def _read_csv_rows(path):
    """Reads a CSV into a list of dicts, coercing numeric-looking fields
    to float (leaving blanks/non-numeric strings untouched)."""
    rows = []
    with open(path, newline='', encoding='utf-8', errors='replace') as f:
        for r in csv.DictReader(f):
            out = {}
            for k, v in r.items():
                if v is None or v == '':
                    out[k] = None
                    continue
                try:
                    out[k] = float(v)
                except ValueError:
                    out[k] = v
            rows.append(out)
    return rows


def _parse_sobol_txt(path):
    """Fallback parser for a results/sobol_results*.txt file (used only
    when SALib is unavailable and run_sobol() can't be called live)."""
    d = {}
    started = False
    with open(path, encoding='utf-8', errors='replace') as f:
        for line in f:
            if line.strip().startswith('parameter'):
                started = True
                continue
            if not started:
                continue
            parts = line.split()
            if len(parts) == 5:
                name, s1, s1c, st, stc = parts
                try:
                    d[name] = float(st)
                except ValueError:
                    break
            else:
                break
    return d


# ══════════════════════════════════════════════════════════════════════════
# Shared helpers
# ══════════════════════════════════════════════════════════════════════════

def _baseline_amr_sweep(spans=None, T_cold_K=291.15, **kwargs):
    """Single-stage AMR span sweep at the repo's standard baseline
    operating point (matches main.py's run_baseline_sweep defaults):
    Gd, 2T, 5kg, 2Hz, mdot=0.08kg/s, eps=0.85, constant parasitic_fraction."""
    if spans is None:
        spans = np.arange(5, 21, 1)
    defaults = dict(material=GADOLINIUM, mu0H_max=2.0, mass_regenerator=5.0,
                     frequency=2.0, fluid_mdot=0.08,
                     regenerator_effectiveness=0.85)
    defaults.update(kwargs)
    amr = AMRSystem(**defaults)
    results = amr.characteristic_curve(T_cold_K, spans)
    return spans, results


# ══════════════════════════════════════════════════════════════════════════
# FIG 01 — Gd mean-field model validation vs. Dan'kov et al. (1998)
# ══════════════════════════════════════════════════════════════════════════

def plot_gd_validation():
    rows = validation.run_validation(verbose=False)
    Bs = [r[0] for r in rows]
    lit = [r[1] for r in rows]
    model = [r[2] for r in rows]
    err = [r[3] for r in rows]

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(11, 4.5))
    x = np.arange(len(Bs))
    w = 0.35
    ax1.bar(x - w / 2, lit, w, label="Literature (Dan'kov et al. 1998)",
            color=COLOR_MAIN, alpha=0.85, edgecolor='white')
    ax1.bar(x + w / 2, model, w, label='Mean-field (Brillouin) model',
            color=COLOR_POWER, alpha=0.85, edgecolor='white')
    for xi, l, m in zip(x, lit, model):
        ax1.text(xi - w / 2, l + 0.2, f'{l:.1f}', ha='center', fontsize=9)
        ax1.text(xi + w / 2, m + 0.2, f'{m:.1f}', ha='center', fontsize=9)
    ax1.set_xticks(x)
    ax1.set_xticklabels([f'{b:.0f} T' for b in Bs])
    ax1.set_ylabel(r'$\Delta T_{ad}$ [K]  (at T$\approx$294 K)')
    ax1.set_title('Gd Adiabatic Temperature Change:\nModel vs. Literature')
    ax1.legend(fontsize=9)

    colors = [COLOR_POWER if e > 0 else COLOR_MAIN for e in err]
    ax2.bar(x, err, color=colors, alpha=0.85, edgecolor='white')
    ax2.axhline(0, color='k', linewidth=0.8)
    ax2.set_xticks(x)
    ax2.set_xticklabels([f'{b:.0f} T' for b in Bs])
    ax2.set_ylabel('Error vs. literature [%]')
    ax2.set_title('Model Error\n(mean-field theory overpredicts near Tc)')
    err_span = max(err) - min(err) if max(err) != min(err) else 1.0
    label_gap = max(err_span * 0.08, 0.6)
    pad = max(err_span * 0.35, 3.0)
    ax2.set_ylim(min(err) - pad, max(err) + pad)
    for xi, e in zip(x, err):
        va = 'bottom' if e >= 0 else 'top'
        ax2.text(xi, e + (label_gap if e >= 0 else -label_gap), f'{e:+.1f}%',
                  ha='center', va=va, fontsize=9)

    fig.suptitle('Mean-Field MCE Model Validation — Gadolinium (Tc = 294 K)',
                 fontsize=13)
    fig.tight_layout()
    save(fig, 'fig01_gd_mce_validation')


# ══════════════════════════════════════════════════════════════════════════
# FIG 02 — Gd entropy change and DeltaT_ad vs. temperature
# ══════════════════════════════════════════════════════════════════════════

def plot_gd_entropy_dTad():
    Ts = np.linspace(250.0, 340.0, 300)
    fields_T = [1, 2, 5]

    fig, axes = plt.subplots(1, 2, figsize=(12, 4.5))
    for B, color in zip(fields_T, COLOR_CYCLE):
        H = B / mu0
        dS = GADOLINIUM.delta_S_isothermal(Ts, H)
        dT = GADOLINIUM.delta_T_adiabatic(Ts, H)
        axes[0].plot(Ts, -dS, color=color, label=f'{B} T')
        axes[1].plot(Ts, dT, color=color, label=f'{B} T')

    axes[0].axvline(GADOLINIUM.Tc, color='grey', linestyle=':', linewidth=1)
    axes[1].axvline(GADOLINIUM.Tc, color='grey', linestyle=':', linewidth=1,
                     label=f'Tc={GADOLINIUM.Tc:.0f}K')
    axes[0].set_xlabel('Temperature [K]')
    axes[0].set_ylabel(r'$-\Delta S_M$ [J/(kg$\cdot$K)]')
    axes[0].set_title('Isothermal Entropy Change')
    axes[1].set_xlabel('Temperature [K]')
    axes[1].set_ylabel(r'$\Delta T_{ad}$ [K]')
    axes[1].set_title('Adiabatic Temperature Change')
    for ax in axes:
        ax.legend(fontsize=9)

    fig.suptitle('Gadolinium Magnetocaloric Effect vs. Temperature '
                 '(Mean-Field / Brillouin Model)', fontsize=12)
    fig.tight_layout()
    save(fig, 'fig02_gd_entropy_and_dTad_vs_T')


# ══════════════════════════════════════════════════════════════════════════
# FIG 03 — Gd5Si2Ge2 first-order Landau model calibration
# ══════════════════════════════════════════════════════════════════════════

def plot_landau_giant_mce():
    Ts = np.linspace(260.0, 300.0, 401)
    fields_T = [1, 2, 5]

    fig, axes = plt.subplots(1, 2, figsize=(12, 4.5))
    for B, color in zip(fields_T, COLOR_CYCLE):
        H = B / mu0
        dS = GD5SI2GE2_FIRST_ORDER.delta_S_isothermal(Ts, H)
        dT = GD5SI2GE2_FIRST_ORDER.delta_T_adiabatic(Ts, H)
        axes[0].plot(Ts, -dS, color=color, label=f'{B} T')
        axes[1].plot(Ts, dT, color=color, label=f'{B} T')
        i_peak = int(np.argmin(dS))
        axes[0].plot(Ts[i_peak], -dS[i_peak], 'o', color=color, ms=5, zorder=4)

    axes[0].axvline(GD5SI2GE2_FIRST_ORDER.Tc, color='grey', linestyle=':',
                     label=f'nominal Tc={GD5SI2GE2_FIRST_ORDER.Tc:.0f}K')
    axes[1].axvline(GD5SI2GE2_FIRST_ORDER.Tc, color='grey', linestyle=':')
    axes[0].set_xlabel('Temperature [K]')
    axes[0].set_ylabel(r'$-\Delta S_M$ [J/(kg$\cdot$K)]')
    axes[0].set_title('Isothermal Entropy Change\n(first-order Landau model)')
    axes[1].set_xlabel('Temperature [K]')
    axes[1].set_ylabel(r'$\Delta T_{ad}$ [K]')
    axes[1].set_title('Adiabatic Temperature Change\n(uncorrected — see Giguère cross-check, Fig. 04)')
    for ax in axes:
        ax.legend(fontsize=9)

    fig.suptitle('Gd5Si2Ge2 "Giant" MCE — Extended (6th-order) Landau Model, '
                 '(A,B,C)=(10,-4,8)', fontsize=12)
    fig.tight_layout()
    save(fig, 'fig03_landau_giant_mce_calibration')


# ══════════════════════════════════════════════════════════════════════════
# FIG 04 — Giguere et al. (1999) direct-measurement cross-check
# ══════════════════════════════════════════════════════════════════════════

def plot_giguere_validation():
    res = giguere_validation.run_validation(verbose=False)

    labels = ['Direct\n(measured)', 'Clausius-\nClapeyron', 'Indirect\n(Maxwell)',
              'This repo\'s model']
    vals = [giguere_validation.GIGUERE_DIRECT_DTAD_7T,
            giguere_validation.GIGUERE_CLAUSIUS_CLAPEYRON_DTAD_7T,
            giguere_validation.GIGUERE_INDIRECT_MAXWELL_DTAD_7T,
            res['model_dTad_7T_K']]
    colors = [COLOR_MAIN, '#5b8db8', '#e07b54', COLOR_POWER]

    fig, ax = plt.subplots(figsize=(7.5, 5.5))
    bars = ax.bar(labels, vals, color=colors, alpha=0.85, edgecolor='white')
    for b, v in zip(bars, vals):
        ax.text(b.get_x() + b.get_width() / 2, v + 0.3, f'{v:.1f} K',
                 ha='center', fontsize=10)
    ax.set_ylabel(r'$\Delta T_{ad}$ at 7 T [K]')
    ax.set_title("Gd5Si2Ge2 $\\Delta T_{ad}$ at 7T: Model vs. Giguère et al. (1999)\n"
                 f"Model overestimates DIRECT measurement by "
                 f"{res['overestimate_factor_vs_direct']:.2f}$\\times$  "
                 f"(paper's own indirect/direct gap: "
                 f"{res['papers_own_overestimate_factor']:.2f}$\\times$)",
                 fontsize=11)
    fig.tight_layout()
    save(fig, 'fig04_giguere_direct_vs_indirect_validation')


# ══════════════════════════════════════════════════════════════════════════
# FIG 05 — Material comparison: Gd vs. Gd5Si2Ge2 vs. La(Fe,Si)13Hy
# ══════════════════════════════════════════════════════════════════════════

def plot_material_comparison():
    fig, ax = plt.subplots(figsize=(8.5, 5.5))
    H2 = 2.0 / mu0

    Ts_gd = np.linspace(250.0, 340.0, 300)
    dT_gd = GADOLINIUM.delta_T_adiabatic(Ts_gd, H2)
    ax.plot(Ts_gd, dT_gd, color=COLOR_MAIN, label=f'Gd (Tc={GADOLINIUM.Tc:.0f}K)')

    Ts_giant = np.linspace(260.0, 300.0, 401)
    dT_giant = GD5SI2GE2_FIRST_ORDER.delta_T_adiabatic(Ts_giant, H2)
    ax.plot(Ts_giant, dT_giant, color=COLOR_POWER,
            label=f'Gd5Si2Ge2 (Tc={GD5SI2GE2_FIRST_ORDER.Tc:.0f}K, first-order, uncorrected)')

    Ts_la = np.linspace(272.0, 312.0, 401)
    dT_la = LAFESIH_FIRST_ORDER.delta_T_adiabatic(Ts_la, H2)
    ax.plot(Ts_la, dT_la, color='#85bb65',
            label=f'La(Fe,Si)13Hy (Tc={LAFESIH_FIRST_ORDER.Tc:.0f}K, first-order)')

    ax.axvspan(291.15, 300.15, color='grey', alpha=0.15,
               label='ASHRAE 18-27°C supply range')
    ax.set_xlabel('Temperature [K]')
    ax.set_ylabel(r'$\Delta T_{ad}$ at 2 T [K]')
    ax.set_title('Material Comparison: Adiabatic Temperature Change vs. Data-Center Range')
    # The red (Gd5Si2Ge2) and green (La(Fe,Si)13Hy) peaks sit at high y
    # in DIFFERENT x-ranges (~286-289K and ~298-300K respectively), so
    # every in-axes corner wide enough to hold this 4-entry legend
    # collides with one peak or the other -- upper-left hits the red
    # rising edge, upper-right hits the green one. Placing the legend
    # below the axes entirely sidesteps this regardless of curve shape.
    ax.legend(fontsize=9, loc='upper center', bbox_to_anchor=(0.5, -0.15), ncol=2)
    fig.tight_layout()
    save(fig, 'fig05_material_comparison_dTad')


# ══════════════════════════════════════════════════════════════════════════
# FIG 06 — Single-stage AMR characteristic curve
# ══════════════════════════════════════════════════════════════════════════

def plot_amr_characteristic_curve():
    spans, results = _baseline_amr_sweep()
    Qc = [r.Qc for r in results]
    cop_e = [r.COP_electrical for r in results]
    cop_i = [r.COP for r in results]

    fig, ax1 = plt.subplots(figsize=(7.5, 5))
    ax2 = ax1.twinx()
    ax1.plot(spans, Qc, color=COLOR_MAIN, marker='o', ms=4, label='Cooling capacity Qc')
    ax2.plot(spans, cop_e, color=COLOR_POWER, marker='s', ms=4, linestyle='--',
             label='COP (electrical)')
    ax2.plot(spans, cop_i, color=COLOR_POWER, marker='^', ms=4, linestyle=':',
             alpha=0.5, label='COP (ideal, magnetic-cycle only)')

    ax1.set_xlabel('Temperature Span [K]')
    ax1.set_ylabel('Qc [W]', color=COLOR_MAIN)
    ax2.set_ylabel('COP', color=COLOR_POWER)
    lines1, labels1 = ax1.get_legend_handles_labels()
    lines2, labels2 = ax2.get_legend_handles_labels()
    ax1.legend(lines1 + lines2, labels1 + labels2, fontsize=9, loc='upper right')
    ax1.set_title('Single-Stage AMR Characteristic Curve\n'
                  '(Gd, 2T, 5kg, 2Hz, mdot=0.08kg/s, T_cold=291.15K)')
    fig.tight_layout()
    save(fig, 'fig06_amr_characteristic_curve')


# ══════════════════════════════════════════════════════════════════════════
# FIG 07 — AMR energy balance decomposition vs. span
# ══════════════════════════════════════════════════════════════════════════

def plot_amr_energy_balance():
    spans, results = _baseline_amr_sweep()
    Qc = np.array([r.Qc for r in results])
    Wmag = np.array([r.W_mag for r in results])
    Wpar = np.array([r.W_parasitic for r in results])

    fig, ax = plt.subplots(figsize=(7.5, 5))
    ax.stackplot(spans, Wmag, Wpar,
                 labels=['Magnetic work W_mag', 'Parasitic power W_parasitic'],
                 colors=['#5b8db8', '#e07b54'], alpha=0.85)
    ax.plot(spans, Qc, 'k-', linewidth=2.0, label='Cooling capacity Qc')
    ax.set_xlabel('Temperature Span [K]')
    ax.set_ylabel('Power [W]')
    ax.legend(fontsize=9)
    ax.set_title('AMR Energy Balance vs. Span\n(constant parasitic_fraction=0.15 default)')
    fig.tight_layout()
    save(fig, 'fig07_amr_energy_balance_vs_span')


# ══════════════════════════════════════════════════════════════════════════
# FIG 08 — AMR electrical COP vs. vapor-compression / liquid cooling / Carnot
# ══════════════════════════════════════════════════════════════════════════

def plot_amr_vs_baselines(precomputed=None):
    """Unlike fig06/fig07 (explicitly "Single-Stage" curves, left as pure
    single-stage on purpose to show that limitation directly), this figure
    claims to be the overall AMR-vs-baselines comparison, so it uses
    core.cascade.staged_baseline_result() -- same fallback as main.py's
    run_baseline_sweep() -- instead of a bare single-stage sweep. Without
    it, spans past ~16K plotted as a hard 0 here (invisible on a log
    y-axis, but still wrong: a real system would just add stages, not
    stop cooling).

    precomputed, if given, may supply 'baseline_rows' -- the exact row
    list step 4 (run_baseline_sweep()) already computed with these same
    T_cold/span/material/field/mass/frequency/mdot/effectiveness values --
    so this figure reuses it instead of re-running all 16 staged_baseline_result()
    + vapor_compression_cop()/liquid_cooling_cop() calls from scratch a
    second time in the same pipeline invocation."""
    precomputed = precomputed or {}
    baseline_rows = precomputed.get('baseline_rows')

    # elastocaloric reference line. A SINGLE static literature
    # value (core/baseline_cooling.py's own honesty flag explains why it
    # is not a span-dependent simulation the way the other three curves
    # are), repeated across every span so it plots as a flat horizontal
    # line -- the same "reference line, not a simulated system" treatment
    # already used for the Carnot curve, just flat instead of span-varying
    # since no elastocaloric COP(span) relation is available to evaluate.
    elasto = elastocaloric_reference_cop()

    T_cold_K = 291.15
    if baseline_rows is not None:
        spans = [r['span_K'] for r in baseline_rows]
        cop_e = [r['AMR_COP_electrical'] for r in baseline_rows]
        vcc_l = [r['VaporCompression_COP'] for r in baseline_rows]
        liq_l = [r['LiquidCooling_COP'] for r in baseline_rows]
        carnot_l = [r['Carnot_COP'] for r in baseline_rows]
        elasto_l = [r.get('Elastocaloric_COP_ref', elasto.COP_representative)
                    for r in baseline_rows]
    else:
        spans = np.arange(5, 21, 1)
        results = [cascade.staged_baseline_result(
            T_cold_K, float(s), material=GADOLINIUM, mu0H_max=2.0,
            mass_regenerator=5.0, frequency=2.0, fluid_mdot=0.08,
            regenerator_effectiveness=0.85) for s in spans]
        cop_e = [r.COP_electrical for r in results]
        vcc_l, liq_l, carnot_l = [], [], []
        for span in spans:
            Th = T_cold_K + span
            v = vapor_compression_cop(T_cold_K, Th)
            l = liquid_cooling_cop(T_cold_K, Th)
            vcc_l.append(v.COP)
            liq_l.append(l.COP)
            carnot_l.append(v.COP_carnot)
        elasto_l = [elasto.COP_representative] * len(spans)

    fig, ax = plt.subplots(figsize=(8, 5.5))
    ax.plot(spans, cop_e, color=COLOR_POWER, marker='o', label='Magnetic (AMR) — electrical COP')
    ax.plot(spans, vcc_l, color=COLOR_MAIN, marker='s', label='Vapor-compression')
    ax.plot(spans, liq_l, color='#85bb65', marker='^', label='Liquid cooling (blended)')
    ax.plot(spans, elasto_l, color='#8e44ad', linestyle='--',
             label='Elastocaloric (literature ref., static — see honesty flag)')
    ax.axhspan(elasto.COP_low, elasto.COP_high, color='#8e44ad', alpha=0.08, lw=0)
    ax.plot(spans, carnot_l, color='grey', linestyle=':', label='Carnot limit')
    ax.set_yscale('log')
    ax.set_xlabel('Temperature Span [K]')
    ax.set_ylabel('COP (log scale)')
    ax.legend(fontsize=9)
    ax.set_title('Data-Center Cooling COP Comparison\nASHRAE 5-20K Span Range, T_cold=18°C')
    fig.tight_layout()
    save(fig, 'fig08_amr_vs_baselines_cop')


# ══════════════════════════════════════════════════════════════════════════
# FIG 09 — NTU regenerator effectiveness vs. mass and frequency
# ══════════════════════════════════════════════════════════════════════════

def plot_regenerator_effectiveness():
    masses = np.array([0.5, 1, 2, 5, 10, 15])
    eps_mass = [regenerator_effectiveness(m, frequency=1.0, mdot=0.08)['eps'] for m in masses]
    freqs = np.array([0.25, 0.5, 1, 2, 4])
    eps_freq = [regenerator_effectiveness(2.0, frequency=f, mdot=0.08)['eps'] for f in freqs]

    fig, axes = plt.subplots(1, 2, figsize=(11, 4.5))
    axes[0].plot(masses, eps_mass, color=COLOR_MAIN, marker='o')
    axes[0].set_xlabel('Regenerator mass [kg]')
    axes[0].set_ylabel(r'Effectiveness $\varepsilon$')
    axes[0].set_title('ε vs. Regenerator Mass\n(f=1Hz, mdot=0.08kg/s)')

    axes[1].plot(freqs, eps_freq, color=COLOR_POWER, marker='s')
    axes[1].set_xlabel('Frequency [Hz]')
    axes[1].set_ylabel(r'Effectiveness $\varepsilon$')
    axes[1].set_title('ε vs. Frequency\n(mass=2kg, mdot=0.08kg/s)')

    fig.suptitle('NTU Packed-Bed Regenerator Effectiveness (Wakao-Kaguei correlation)',
                 fontsize=12)
    fig.tight_layout()
    save(fig, 'fig09_regenerator_effectiveness_ntu')


# ══════════════════════════════════════════════════════════════════════════
# FIG 10 / 11 — Geometry-dependent pumping power trade-offs
# ══════════════════════════════════════════════════════════════════════════

def plot_geometry_packed_bed():
    rows, best_qc, best_cop = geometry_analysis.sweep_packed_bed_diameter(verbose=False)
    d = [r[0] for r in rows]
    qc = [r[1] for r in rows]
    cop = [r[2] for r in rows]

    fig, ax1 = plt.subplots(figsize=(7.5, 5.5))
    ax2 = ax1.twinx()
    ax1.plot(d, qc, color=COLOR_MAIN, marker='o', label='Qc')
    ax2.plot(d, cop, color=COLOR_POWER, marker='s', linestyle='--', label='COP_aug (hydraulic pumping incl.)')
    ax2.axvline(best_cop[0], color=COLOR_POWER, linestyle=':', linewidth=1)
    ax2.annotate(f'COP optimum\n{best_cop[0]}mm', xy=(best_cop[0], best_cop[2]),
                 xytext=(10, -25), textcoords='offset points', fontsize=8, color=COLOR_POWER,
                 arrowprops=dict(arrowstyle='->', color=COLOR_POWER))
    ax1.set_xscale('log')
    ax1.set_xlabel('Sphere particle diameter [mm] (log)')
    ax1.set_ylabel('Qc [W]', color=COLOR_MAIN)
    ax2.set_ylabel('COP_aug (idealized, no pump/motor efficiency)', color=COLOR_POWER)
    ax1.set_title('Packed-Bed Geometry Trade-off\n'
                  '(Tušek, Kitanovski, Poredoš 2013 friction correlation; fixed mdot=0.08kg/s)')
    fig.tight_layout()
    save(fig, 'fig10_geometry_optimum_packed_bed')


def plot_geometry_parallel_plate():
    rows, best_qc, best_cop = geometry_analysis.sweep_parallel_plate_spacing(verbose=False)
    s = [r[0] for r in rows]
    qc = [r[1] for r in rows]
    cop = [r[2] for r in rows]

    fig, ax1 = plt.subplots(figsize=(7.5, 5.5))
    ax2 = ax1.twinx()
    ax1.plot(s, qc, color=COLOR_MAIN, marker='o', label='Qc')
    ax2.plot(s, cop, color=COLOR_POWER, marker='s', linestyle='--', label='COP_aug')
    ax2.axvline(best_cop[0], color=COLOR_POWER, linestyle=':', linewidth=1)
    ax2.annotate(f'COP optimum\n{best_cop[0]}mm', xy=(best_cop[0], best_cop[2]),
                 xytext=(10, -25), textcoords='offset points', fontsize=8, color=COLOR_POWER,
                 arrowprops=dict(arrowstyle='->', color=COLOR_POWER))
    ax1.set_xscale('log')
    ax1.set_xlabel('Plate spacing [mm] (log)')
    ax1.set_ylabel('Qc [W]', color=COLOR_MAIN)
    ax2.set_ylabel('COP_aug', color=COLOR_POWER)
    ax1.set_title('Parallel-Plate Geometry Trade-off\n'
                  '(fixed plate thickness=0.25mm, mdot=0.08kg/s)')
    fig.tight_layout()
    save(fig, 'fig11_geometry_optimum_parallel_plate')


# ══════════════════════════════════════════════════════════════════════════
# FIG 12 — State-dependent loss-model calibration fit + leave-one-out CV
# ══════════════════════════════════════════════════════════════════════════

def plot_loss_model_calibration():
    cal = loss_model_mod.calibrate_loss_coefficients(verbose=False)
    A, b = loss_model_mod._build_system(CALIBRATION_POINTS_CORE)
    pred = A @ cal['raw']
    names = [p[0].replace('_', ' ') for p in CALIBRATION_POINTS_CORE]
    true = [p[5] for p in CALIBRATION_POINTS_CORE]

    loo = loss_model_mod.leave_one_out_cv(CALIBRATION_POINTS_EXTENDED, verbose=False)
    loo_names = [r[0].replace('_', ' ') for r in loo]
    loo_err = [r[3] for r in loo]

    fig, axes = plt.subplots(1, 2, figsize=(13, 5.5))
    x = np.arange(len(names))
    w = 0.35
    axes[0].bar(x - w / 2, true, w, label='Required (from device COP)',
                color=COLOR_MAIN, alpha=0.85, edgecolor='white')
    axes[0].bar(x + w / 2, pred, w, label='NNLS fit',
                color=COLOR_POWER, alpha=0.85, edgecolor='white')
    axes[0].set_xticks(x)
    axes[0].set_xticklabels(names, rotation=15, ha='right', fontsize=8)
    axes[0].set_ylabel('W_parasitic [W]')
    axes[0].set_title('CORE 3-Point Fit\n(exactly-determined, NNLS)')
    axes[0].legend(fontsize=9)

    colors = [COLOR_POWER if abs(e) > 100 else COLOR_MAIN for e in loo_err]
    axes[1].barh(loo_names, loo_err, color=colors, alpha=0.85, edgecolor='white')
    axes[1].axvline(0, color='k', linewidth=0.8)
    axes[1].set_xlabel('Leave-one-out error [%]')
    axes[1].set_title('Leave-One-Out CV (EXTENDED 4-pt set)\n'
                       'Model cannot generalize across 6.5W-2502W devices')

    fig.suptitle('State-Dependent Loss Model Calibration (loss_model.py)', fontsize=13)
    fig.tight_layout()
    save(fig, 'fig12_loss_model_calibration_fit')


# ══════════════════════════════════════════════════════════════════════════
# FIG 13 — Parasitic-fraction scaling with device size (non-monotonicity)
# ══════════════════════════════════════════════════════════════════════════

def plot_parasitic_fraction_scaling():
    rows = loss_model_mod.analyze_parasitic_fraction_scaling(
        CALIBRATION_POINTS_FURTHER_EXTENDED, verbose=False)
    names = [r[0].replace('_', ' ') for r in rows]
    Qc = [r[1] for r in rows]
    frac = [r[2] for r in rows]

    fig, ax = plt.subplots(figsize=(9.5, 6))
    ax.scatter(Qc, frac, c=np.arange(len(Qc)), cmap='viridis', s=90,
                     edgecolor='white', zorder=3)
    # Points that sit close together in (log Qc, frac) space -- e.g. the
    # four Lozano POLO UFSC rows -- would otherwise get directly
    # overlapping text labels with one fixed offset. Fan the offset out
    # further for each point that lands near an already-placed one.
    log_qc = np.log10(Qc)
    x_span = (max(log_qc) - min(log_qc)) or 1.0
    y_span = (max(frac) - min(frac)) or 1.0
    placed = []
    for qc_, f_, n_, lq in zip(Qc, frac, names, log_qc):
        nearby = sum(1 for plx, ply in placed
                     if abs(plx - lq) / x_span < 0.05 and abs(ply - f_) / y_span < 0.05)
        placed.append((lq, f_))
        dx, dy = 5 + 8 * nearby, 4 + 13 * nearby
        ax.annotate(n_, (qc_, f_), fontsize=7, xytext=(dx, dy), textcoords='offset points')
    ax.set_xscale('log')
    ax.set_xlabel('Device cooling capacity Qc [W] (log)')
    ax.set_ylabel('Parasitic fraction W_parasitic / Qc')
    ax.set_title('Parasitic Fraction vs. Device Scale\n'
                 '(no monotonic size trend — Astronautics outlier attributed to '
                 '"mediocre" electrical-component efficiency, not scale)')
    fig.tight_layout()
    save(fig, 'fig13_parasitic_fraction_scaling')


# ══════════════════════════════════════════════════════════════════════════
# FIG 14 — System-level validation vs. published AMR prototypes
# ══════════════════════════════════════════════════════════════════════════

def plot_system_validation(precomputed=None):
    """precomputed, if given, may supply 'system_validation_results' --
    the exact return value of validation_system.run_system_validation()
    already computed by step 2 in main.py -- so this figure reuses it
    instead of re-running that deterministic calibration search (16
    brentq calibrations against data/amr_experimental_benchmarks.csv)
    a second time in the same pipeline invocation."""
    precomputed = precomputed or {}
    results = precomputed.get('system_validation_results')
    if results is None:
        results = validation_system.run_system_validation()
    ok = [r for r in results if 'COP_error_pct' in r]
    names = [r['device'].replace('_', ' ') for r in ok]
    cop_lit = [r['COP_lit'] for r in ok]
    cop_model = [r['COP_model_electrical'] for r in ok]
    err = [r['COP_error_pct'] for r in ok]

    fig, axes = plt.subplots(2, 2, figsize=(13, 11))
    ax1, ax2, ax3, ax4 = axes[0, 0], axes[0, 1], axes[1, 0], axes[1, 1]

    good_idx = [i for i, e in enumerate(err) if abs(e) <= 30]
    bad_idx = [i for i, e in enumerate(err) if abs(e) > 30]

    def scatter_panel(ax, idx, color, title):
        if not idx:
            ax.set_visible(False)
            return
        x_ = [cop_lit[i] for i in idx]
        y_ = [cop_model[i] for i in idx]
        n_ = [names[i] for i in idx]
        ax.scatter(x_, y_, s=90, color=color, zorder=3, edgecolor='white')
        lims = [0, max(max(x_), max(y_)) * 1.15]
        ax.plot(lims, lims, 'k--', linewidth=1, label='Perfect agreement')
        # Points that sit close together produce overlapping text if
        # labeled in place, so group nearby points into clusters (in
        # axis-fraction units): isolated points keep inline labels, and
        # any cluster of 2+ points is labeled with small numbered
        # markers plus a legend box, which never overlaps.
        span = max(lims[1], 1e-9)
        clusters = []
        for xx, yy, nn in zip(x_, y_, n_):
            for c in clusters:
                if abs(c['x'][0] - xx) / span < 0.06 and abs(c['y'][0] - yy) / span < 0.06:
                    c['x'].append(xx); c['y'].append(yy); c['n'].append(nn)
                    break
            else:
                clusters.append({'x': [xx], 'y': [yy], 'n': [nn]})

        legend_lines = []
        counter = 1
        for c in clusters:
            if len(c['n']) == 1:
                ax.annotate(c['n'][0], (c['x'][0], c['y'][0]), fontsize=7,
                            xytext=(5, 5), textcoords='offset points')
            else:
                for xx, yy, nn in zip(c['x'], c['y'], c['n']):
                    ax.annotate(str(counter), (xx, yy), fontsize=7, fontweight='bold',
                                xytext=(0, 0), textcoords='offset points',
                                ha='center', va='center', color='white')
                    legend_lines.append(f"{counter}: {nn}")
                    counter += 1
        if legend_lines:
            ax.text(0.98, 0.02, '\n'.join(legend_lines), transform=ax.transAxes,
                    fontsize=6.5, ha='right', va='bottom',
                    bbox=dict(boxstyle='round', facecolor='white', alpha=0.85, edgecolor='0.7'))
        ax.set_xlim(lims)
        ax.set_ylim(lims)
        ax.set_xlabel('Literature COP (electrical)')
        ax.set_ylabel('Model COP (electrical)')
        ax.set_title(title)
        ax.legend(fontsize=8)

    def bar_panel(ax, idx, color, title):
        if not idx:
            ax.set_visible(False)
            return
        n_ = [names[i] for i in idx]
        e_ = [err[i] for i in idx]
        ax.barh(n_, e_, color=color, alpha=0.85, edgecolor='white')
        ax.axvline(0, color='k', linewidth=0.8)
        ax.set_xlabel('COP error [%]')
        ax.set_title(title)

    scatter_panel(ax1, good_idx, COLOR_MAIN, 'Well-Modeled Devices')
    scatter_panel(ax2, bad_idx, COLOR_POWER, 'Outlier Devices')
    bar_panel(ax3, good_idx, COLOR_MAIN, 'COP Error: Well-Modeled Devices')
    bar_panel(ax4, bad_idx, COLOR_POWER, 'COP Error: Outlier Devices')

    for ax in axes.flat:
        ax.tick_params(labelsize=8)
    fig.suptitle('System-Level Validation: Model vs. Literature\n'
                 '(fluid mdot calibrated to reproduce reported Qc)', fontsize=13)
    fig.tight_layout(rect=[0, 0, 1, 0.95])
    save(fig, 'fig14_system_validation_scatter')


# ══════════════════════════════════════════════════════════════════════════
# FIG 15 — Curve-level (2-point) Qc(span) shape validation
# ══════════════════════════════════════════════════════════════════════════

def plot_curve_validation():
    results = validation_system.run_curve_validation()
    ok = [r for r in results if 'companion_Qc_model_W' in r]
    names = [r['device_group'].replace('_', ' ') for r in ok]
    lit = [r['companion_Qc_lit_W'] for r in ok]
    model = [r['companion_Qc_model_W'] for r in ok]

    fig, ax = plt.subplots(figsize=(8.5, 5.5))
    x = np.arange(len(names))
    w = 0.35
    bars_lit = ax.bar(x - w / 2, lit, w, label='Literature (companion span point)',
                       color=COLOR_MAIN, alpha=0.85, edgecolor='white')
    bars_model = ax.bar(x + w / 2, model, w, label='Model (predicted from anchor-point calibration)',
                         color=COLOR_POWER, alpha=0.85, edgecolor='white')
    ax.set_xticks(x)
    ax.set_xticklabels(names, rotation=10, ha='right', fontsize=9)
    ax.set_ylabel('Qc at companion span [W]')

    # Explicit value labels: without these, a genuinely-zero (or near-zero)
    # result is visually indistinguishable from an empty/broken chart --
    # which is exactly what happens for the n=1 case below (Tusek's own
    # companion point is a reported zero-capacity max-span point, so both
    # bars are legitimately 0.0W, not missing data).
    ax.bar_label(bars_lit, fmt='%.1f W', fontsize=8, padding=3)
    ax.bar_label(bars_model, fmt='%.1f W', fontsize=8, padding=3)
    ax.set_ylim(bottom=0)  # avoid matplotlib auto-scaling to a +/-0.04 W
                           # range around an all-zero series, which hides
                           # the (correct) 0-vs-0 result behind empty axes

    if len(ok) == 0:
        ax.set_title('Curve-Level (2-Point) Validation:\n'
                      'No device currently has a usable companion point '
                      '(see console log for why each was skipped)')
    elif len(ok) == 1:
        ax.set_title('Curve-Level (2-Point) Validation: Qc(span) Shape Check\n'
                      f'n=1 device has a usable second span point ({names[0]})')
        note = ('Companion point is this device\'s own reported zero-capacity\n'
                '(max-span) point -- model correctly predicts 0W, matching literature.'
                if abs(lit[0]) < 1e-6 and abs(model[0]) < 1e-6 else
                'Only one benchmark device currently has a second, independent\n'
                'span point usable for this shape check; see run_curve_validation().')
        ax.annotate(note, xy=(0.5, 0.92), xycoords='axes fraction', ha='center', va='top',
                    fontsize=8.5, style='italic', color='dimgray')
    else:
        ax.set_title('Curve-Level (2-Point) Validation:\n'
                      'Qc(span) Shape Check — Companion Point NOT Used in Calibration')
    ax.legend(fontsize=9)
    fig.tight_layout()
    save(fig, 'fig15_curve_validation_companion')


# ══════════════════════════════════════════════════════════════════════════
# FIG 16 — Sobol sensitivity: constant-loss vs. state-dependent loss model
# ══════════════════════════════════════════════════════════════════════════

def plot_sobol_sensitivity(precomputed=None):
    """precomputed, if given, may supply 'sobol_const_Si' and/or
    'sobol_state_Si' (the Si dicts already returned by steps 9/9b in
    main.py) so this figure reuses them instead of re-running the ~3s
    SALib Sobol analysis a second time in the same pipeline invocation."""
    names_display = {
        'mu0H_max_T': 'Field\n(mu0H_max)', 'frequency_Hz': 'Frequency',
        'fluid_mdot_kgs': 'Flow rate\n(mdot)', 'regen_effectiveness': 'Regen.\neffectiveness',
        'parasitic_fraction': 'Parasitic\nfraction',
    }

    precomputed = precomputed or {}
    Si_const = precomputed.get('sobol_const_Si')
    Si_state = precomputed.get('sobol_state_Si')

    if Si_const is not None and Si_state is not None:
        names = sensitivity_mod.PROBLEM['names']
        st_const = dict(zip(names, Si_const['ST']))
        st_state = dict(zip(names, Si_state['ST']))
    elif HAVE_SALIB:
        Si_const = sensitivity_mod.run_sobol(
            out_path=str(RESULTS_DIR / 'sobol_results_constant_losses.txt'),
            use_state_dependent_losses=False)
        Si_state = sensitivity_mod.run_sobol(
            out_path=str(RESULTS_DIR / 'sobol_results.txt'),
            use_state_dependent_losses=True)
        names = sensitivity_mod.PROBLEM['names']
        st_const = dict(zip(names, Si_const['ST']))
        st_state = dict(zip(names, Si_state['ST']))
    else:
        print("  [SALib unavailable — falling back to pre-computed results/sobol_results*.txt]")
        st_const = _parse_sobol_txt(RESULTS_DIR / 'sobol_results_constant_losses.txt')
        st_state = _parse_sobol_txt(RESULTS_DIR / 'sobol_results.txt')
        names = list(st_state.keys()) or list(st_const.keys())

    labels = [names_display.get(n, n) for n in names]
    const_vals = [st_const.get(n, 0.0) for n in names]
    state_vals = [st_state.get(n, 0.0) for n in names]

    fig, ax = plt.subplots(figsize=(9.5, 5.5))
    x = np.arange(len(names))
    w = 0.35
    ax.bar(x - w / 2, const_vals, w, label='Constant parasitic_fraction',
           color=COLOR_MAIN, alpha=0.85, edgecolor='white')
    ax.bar(x + w / 2, state_vals, w, label='State-dependent loss model',
           color=COLOR_POWER, alpha=0.85, edgecolor='white')
    ax.set_xticks(x)
    ax.set_xticklabels(labels, fontsize=8)
    ax.set_ylabel('Total-order Sobol index  $S_T$')
    ax.set_title('Sobol Sensitivity of Electrical COP\n'
                 '(T_cold=291K, span=10K, N_base=64 Saltelli samples)')
    ax.legend(fontsize=9)
    fig.tight_layout()
    save(fig, 'fig16_sobol_sensitivity_comparison')


# ══════════════════════════════════════════════════════════════════════════
# FIG 17 — Quadratic RSM surrogate for cooling capacity
# ══════════════════════════════════════════════════════════════════════════

def plot_rsm_surrogate():
    X_train = rsm_mod.sample_design(300, seed=7)
    Y_train = rsm_mod.evaluate(X_train)
    X_test = rsm_mod.sample_design(100, seed=99)
    Y_test = rsm_mod.evaluate(X_test)

    Phi_train, feat_names = rsm_mod.build_quadratic_features(X_train, rsm_mod.VAR_NAMES)
    coeffs, *_ = np.linalg.lstsq(Phi_train, Y_train, rcond=None)
    Phi_test, _ = rsm_mod.build_quadratic_features(X_test, rsm_mod.VAR_NAMES)
    Y_pred = Phi_test @ coeffs
    ss_res = np.sum((Y_test - Y_pred) ** 2)
    ss_tot = np.sum((Y_test - np.mean(Y_test)) ** 2)
    r2 = 1 - ss_res / ss_tot if ss_tot > 0 else float('nan')

    fig, axes = plt.subplots(1, 2, figsize=(13, 5.5))
    axes[0].scatter(Y_test, Y_pred, s=28, color=COLOR_MAIN, alpha=0.7, zorder=3)
    lims = [0, max(Y_test.max(), Y_pred.max()) * 1.05]
    axes[0].plot(lims, lims, 'k--', linewidth=1, label='Perfect agreement')
    axes[0].set_xlim(lims)
    axes[0].set_ylim(lims)
    axes[0].set_xlabel('Actual Qc [W]  (AMRSystem.run())')
    axes[0].set_ylabel('RSM-predicted Qc [W]')
    axes[0].set_title(f'Quadratic RSM Surrogate: Held-out Parity\n'
                       f'R\u00b2={r2:.3f}, n_train=300, n_test=100')
    axes[0].legend(fontsize=9)

    top_idx = np.argsort(-np.abs(coeffs))[:10]
    top_names = [feat_names[i] for i in top_idx][::-1]
    top_vals = [coeffs[i] for i in top_idx][::-1]
    axes[1].barh(top_names, top_vals, color=COLOR_POWER, alpha=0.85, edgecolor='white')
    axes[1].set_xlabel('Coefficient value')
    axes[1].set_title('Top 10 RSM Coefficients by Magnitude')
    axes[1].tick_params(labelsize=8)
    fig.tight_layout()
    save(fig, 'fig17_rsm_surrogate_parity')


# ══════════════════════════════════════════════════════════════════════════
# FIG 18 — NSGA-III multi-objective Pareto front
# ══════════════════════════════════════════════════════════════════════════

def plot_nsga3_pareto(precomputed=None):
    """precomputed, if given, may supply 'pareto_rows' (already produced by
    step 11 in main.py) so this figure reuses that front instead of
    re-running the ~7s NSGA-III optimization a second time in the same
    pipeline invocation."""
    precomputed = precomputed or {}
    rows = precomputed.get('pareto_rows')
    if rows is None:
        if HAVE_PYMOO:
            rows = optimize_mod.run_optimization(out_csv=str(RESULTS_DIR / 'pareto_front.csv'))
        else:
            print("  [pymoo unavailable — falling back to pre-computed results/pareto_front.csv]")
            rows = _read_csv_rows(RESULTS_DIR / 'pareto_front.csv')

    cop = [r['COP_electrical'] for r in rows]
    qc = [r['Qc_W'] for r in rows]
    cost = [r['cost_index_USD'] for r in rows]

    fig, ax = plt.subplots(figsize=(8.5, 6.5))
    materials = [r.get('material', 'Gd') for r in rows]
    unique_materials = sorted(set(materials))
    # phase 15: material is now a co-optimized design variable, not fixed
    # at Gd -- if more than one material actually appears on the merged
    # global front, distinguish them by marker shape (color still encodes
    # cost, as before this pass) so a reader can see at a glance which
    # material wins in which region of the COP/Qc trade-off.
    if len(unique_materials) > 1:
        marker_cycle = ['o', '^', 's', 'D', 'P', 'X']
        sc = None
        for i, mat in enumerate(unique_materials):
            idx = [j for j, m in enumerate(materials) if m == mat]
            sc = ax.scatter([qc[j] for j in idx], [cop[j] for j in idx],
                             c=[cost[j] for j in idx], cmap='viridis',
                             vmin=min(cost), vmax=max(cost),
                             s=55, alpha=0.85, edgecolor='white', linewidth=0.4,
                             marker=marker_cycle[i % len(marker_cycle)], label=mat)
        ax.legend(title='Material (phase 15)', fontsize=8, loc='best')
    else:
        sc = ax.scatter(qc, cop, c=cost, cmap='viridis', s=45, alpha=0.85,
                         edgecolor='white', linewidth=0.3)
    cbar = fig.colorbar(sc, ax=ax)
    cbar.set_label('Materials cost index [$]')
    ax.set_xlabel('Cooling capacity Qc [W]')
    ax.set_ylabel('Electrical COP')
    ax.set_title('NSGA-III Pareto-Optimal AMR Designs\n'
                 '(T_cold=291K, span=10K; state-dependent loss + NTU thermal model;\n'
                 'material + regenerator geometry co-optimized, phase 15)')
    fig.tight_layout()
    save(fig, 'fig18_nsga3_pareto_front')


# ══════════════════════════════════════════════════════════════════════════
# FIG 19 — Multi-stage cascade AMR staging vs. baselines (Gd)
# ══════════════════════════════════════════════════════════════════════════

def plot_cascade_staging_gd(precomputed=None):
    """precomputed, if given, may supply 'cascade_rows_gd' (already
    computed and written to results/cascade_comparison.csv by step 7 in
    main.py) so this figure reuses it instead of re-running the Gd cascade
    sweep a second time in the same pipeline invocation."""
    precomputed = precomputed or {}
    rows = precomputed.get('cascade_rows_gd')
    if rows is None:
        rows = cascade.compare_staging(material=GADOLINIUM, mass_per_stage=5.0,
                                        out_csv=str(RESULTS_DIR / 'cascade_comparison.csv'))
    spans = [r['span_K'] for r in rows]

    fig, ax = plt.subplots(figsize=(8.5, 6))
    for n, color in zip([1, 2, 3, 4], COLOR_CYCLE):
        y = [r[f'AMR_{n}stage_COP'] for r in rows]
        ax.plot(spans, y, marker='o', ms=4, color=color, label=f'{n}-stage AMR')
    ax.plot(spans, [r['VaporCompression_COP'] for r in rows], 'k--', label='Vapor-compression')
    ax.plot(spans, [r['LiquidCooling_COP'] for r in rows], 'k:', label='Liquid cooling')
    ax.set_xlabel('Temperature Span [K]')
    ax.set_ylabel('Electrical COP')
    ax.set_title('Multi-Stage Cascade AMR vs. Baselines\n(Gd, 2T, 5kg/stage)')
    ax.legend(fontsize=8)
    fig.tight_layout()
    save(fig, 'fig19_cascade_staging_gd')


# ══════════════════════════════════════════════════════════════════════════
# FIG 20 — Cascade: Gd vs. Gd5Si2Ge2 (fixed composition)
# ══════════════════════════════════════════════════════════════════════════

def plot_cascade_giant_vs_gd(precomputed=None):
    """Fixed-composition Gd5Si2Ge2 (Tc=276K, fixed) vs. Gd, PLUS a third
    series: the same Gd5(SixGe1-x)4(-Ga) family re-tuned per span so its
    own peak lands at that span's T_mid (core.cascade.GD_FAMILY /
    _target_composition_for_peak -- the same machinery
    material_family_comparison.py already uses for fig26). The
    fixed-composition curve is kept exactly as before (still collapses to
    ~0, honestly, since a fixed Tc really is that far from this range) --
    the tuned curve alongside it shows that the collapse is a
    fixed-composition artifact, not a ceiling on the giant-MCE effect
    itself: composition-tuning the same family recovers real performance
    at every span its documented Tc window covers.

    precomputed, if given, may supply 'cascade_rows_gd' and
    'cascade_rows_giant' (both already computed by step 7 in main.py, which
    runs both materials) so this figure reuses them instead of re-running
    both cascade sweeps a second time in the same pipeline invocation."""
    precomputed = precomputed or {}
    rows_gd = precomputed.get('cascade_rows_gd')
    rows_giant = precomputed.get('cascade_rows_giant')
    if rows_gd is None:
        rows_gd = cascade.compare_staging(material=GADOLINIUM, mass_per_stage=5.0,
                                           out_csv=str(RESULTS_DIR / 'cascade_comparison.csv'))
    if rows_giant is None:
        rows_giant = cascade.compare_staging(material=GD5SI2GE2_FIRST_ORDER, mass_per_stage=5.0,
                                              out_csv=str(RESULTS_DIR / 'cascade_comparison_giant_mce.csv'))
    spans = [r['span_K'] for r in rows_gd]
    y_gd = [r['AMR_1stage_COP'] or 0.0 for r in rows_gd]
    y_giant = [r['AMR_1stage_COP'] or 0.0 for r in rows_giant]

    T_cold_K = 18.0 + 273.15
    y_tuned, tuned_in_range = [], []
    for span in spans:
        T_mid = T_cold_K + span / 2.0
        tc = cascade._target_composition_for_peak(T_mid, 2.0, cascade.GD_FAMILY)
        in_range = cascade.GD_FAMILY.tc_min <= tc <= cascade.GD_FAMILY.tc_max
        tuned_in_range.append(in_range)
        if in_range:
            material = cascade.GD_FAMILY.tuned_fn(tc)
            res = cascade.run_cascade(T_cold_K, span, 1, material=material, mass_per_stage=5.0)
            y_tuned.append(res['COP_cascade'] if res['feasible'] else 0.0)
        else:
            y_tuned.append(np.nan)  # outside the family's documented Tc window -- not plotted

    fig, ax = plt.subplots(figsize=(8.5, 5.5))
    ax.plot(spans, y_gd, color=COLOR_MAIN, marker='o', label='Gd (1-stage)')
    ax.plot(spans, y_giant, color=COLOR_POWER, marker='s',
            label='Gd5Si2Ge2 (1-stage, fixed Tc=276K)')
    ax.plot(spans, y_tuned, color='#85bb65', marker='^',
            label='Gd5(SixGe1-x)4(-Ga) (1-stage, tuned per span)')
    if not all(tuned_in_range):
        first_oor = spans[tuned_in_range.index(False)] if False in tuned_in_range else None
        if first_oor is not None:
            ax.axvline(first_oor, color='#85bb65', linestyle=':', alpha=0.5, linewidth=1,
                       label="Tuned family's documented Tc window ends")
    ax.set_xlabel('Temperature Span [K]')
    ax.set_ylabel('Electrical COP')
    ax.set_title('Gd vs. Gd5Si2Ge2 in the ASHRAE Range\n'
                 '(fixed-Tc Gd5Si2Ge2 collapses to ~0; the same family re-tuned per span does not)')
    ax.legend(fontsize=8.5)
    fig.tight_layout()
    save(fig, 'fig20_cascade_giant_mce_vs_gd')


# ══════════════════════════════════════════════════════════════════════════
# FIG 21 — Curie-graded cascade performance
# ══════════════════════════════════════════════════════════════════════════

def plot_graded_cascade(precomputed=None):
    """precomputed, if given, may supply 'graded_rows' (already computed by
    step 7b in main.py) so this figure reuses it instead of re-running the
    ~90-120s Curie-graded cascade sweep a second time in the same pipeline
    invocation."""
    precomputed = precomputed or {}
    rows_graded = precomputed.get('graded_rows')
    if rows_graded is None:
        rows_graded, _stage_info_all = cascade.compare_graded_cascade(
            T_cold_C=18.0, spans=range(5, 21), mass_per_stage=5.0,
            out_csv=str(RESULTS_DIR / 'graded_cascade_comparison.csv'))
    spans = [r['span_K'] for r in rows_graded]

    fig, axes = plt.subplots(1, 2, figsize=(13, 5.5))
    for n, color in zip([1, 2, 3, 4], COLOR_CYCLE):
        cop = [r[f'Graded_{n}stage_COP'] for r in rows_graded]
        qc = [r[f'Graded_{n}stage_Qc_W'] for r in rows_graded]
        axes[0].plot(spans, cop, marker='o', ms=4, color=color, label=f'{n}-stage')
        axes[1].plot(spans, qc, marker='o', ms=4, color=color, label=f'{n}-stage')
    axes[0].set_xlabel('Span [K]')
    axes[0].set_ylabel('Electrical COP')
    axes[0].set_title('Curie-Graded Cascade COP')
    axes[1].set_xlabel('Span [K]')
    axes[1].set_ylabel('Qc [W]')
    axes[1].set_title('Curie-Graded Cascade Cooling Capacity')
    for ax in axes:
        ax.legend(fontsize=8)
    fig.suptitle('Composition-Tuned Curie-Graded Cascade '
                 '(Gd5(SixGe1-x)4(-Ga) family, Giguère-corrected)', fontsize=12)
    fig.tight_layout()
    save(fig, 'fig21_graded_cascade_performance')


# ══════════════════════════════════════════════════════════════════════════
# FIG 22 — Economics: simplified TCO and lifetime cost breakdown
# ══════════════════════════════════════════════════════════════════════════

def plot_economics():
    capacity_kW = 1.2876  # AMR_Qc_W at the 10K-span baseline point (results/comparison_table.csv)
    cop_electrical = 4.63

    rows_simple = [economics.simple_tco(tco, capacity_kW, annual_hours=8760)
                   for tco in (economics.AMR_MAGNETIC, economics.VAPOR_COMPRESSION,
                               economics.LIQUID_COOLING)]
    lc = economics.lifetime_cost(mu0H_max=2.0, mass_regenerator=5.0,
                                  Qc_avg_W=capacity_kW * 1000.0,
                                  COP_electrical=cop_electrical,
                                  device_lifetime_years=15.0)

    fig, axes = plt.subplots(1, 2, figsize=(13, 5.5))
    names = [r['technology'] for r in rows_simple]
    capex = [r['capex_$'] for r in rows_simple]
    opex = [r['annual_opex_$'] for r in rows_simple]
    x = np.arange(len(names))
    axes[0].bar(x - 0.2, capex, width=0.4, color=COLOR_MAIN, alpha=0.85, label='CAPEX')
    ax0b = axes[0].twinx()
    ax0b.bar(x + 0.2, opex, width=0.4, color=COLOR_POWER, alpha=0.85, label='Annual OPEX')
    axes[0].set_xticks(x)
    axes[0].set_xticklabels(names, fontsize=8, rotation=10, ha='right')
    axes[0].set_ylabel('CAPEX [$]', color=COLOR_MAIN)
    ax0b.set_ylabel('Annual OPEX [$]', color=COLOR_POWER)
    axes[0].set_title(f'Simplified TCO ({capacity_kW:.2f} kW capacity)')

    labels2 = ['Materials floor\n(magnet + MCM)', 'Lifetime electricity\n(15yr, $0.10/kWh)']
    vals2 = [lc['materials_floor_$'], lc['lifetime_electricity_$']]
    axes[1].bar(labels2, vals2, color=[COLOR_MAIN, COLOR_POWER], alpha=0.85, edgecolor='white')
    for i, v in enumerate(vals2):
        axes[1].text(i, v, f'${v:,.0f}', ha='center', va='bottom', fontsize=9)
    axes[1].set_ylabel('$')
    axes[1].set_title('AMR Lifetime Cost Breakdown\n(Bjørk, Bahl & Nielsen 2016 methodology)')

    fig.suptitle('Economics: CAPEX/OPEX and Lifetime Cost', fontsize=13)
    fig.tight_layout()
    save(fig, 'fig22_economics_tco_comparison')


# ══════════════════════════════════════════════════════════════════════════
# FIG 23 — Refrigerant-free emissions comparison
# ══════════════════════════════════════════════════════════════════════════

def plot_emissions():
    results = emissions.compare_emissions(100.0, amr_cop=4.63, vcc_cop=12.23, liquid_cop=19.89)
    names = [r.technology for r in results]
    refrig = [r.refrigerant_GWP_tCO2e_per_year for r in results]
    op = [r.operational_CO2_tCO2e_per_year for r in results]

    fig, ax = plt.subplots(figsize=(8.5, 5.5))
    x = np.arange(len(names))
    ax.bar(x, refrig, label='Refrigerant leakage', color='#e07b54', alpha=0.85, edgecolor='white')
    ax.bar(x, op, bottom=refrig, label='Operational (electricity)', color=COLOR_MAIN,
           alpha=0.85, edgecolor='white')
    ax.set_xticks(x)
    ax.set_xticklabels(names, fontsize=9, rotation=10, ha='right')
    ax.set_ylabel(r'tCO$_2$e / year')
    ax.set_title('Annual Emissions Comparison, 100kW Cooling Capacity\n'
                 '(COPs from the ASHRAE 10K-span operating point)')
    ax.legend(fontsize=9)
    fig.tight_layout()
    save(fig, 'fig23_emissions_comparison')


# ══════════════════════════════════════════════════════════════════════════
# FIG 24 — Giant-MCE targeting: material must match the operating point
# ══════════════════════════════════════════════════════════════════════════

def plot_giant_mce_targeting():
    """Adds a 4th bar to the original 3-bar comparison: the Gd5(SixGe1-x)4
    (-Ga) family composition-tuned so its OWN peak lands at the ASHRAE
    point's T_mid (same core.cascade.GD_FAMILY machinery as fig20/fig26),
    evaluated with the same eval_at() cycle settings as the other three
    bars. This is the direct fix for "Gd5Si2Ge2 @ ASHRAE" reading as
    ~0: that bar is a FIXED composition (Tc=276K) evaluated far from its
    own transition on purpose (see giant_mce_analysis.py's module
    docstring); the new bar shows the same family, tuned, performing
    normally at the same operating point."""
    _loss = StateDependentLossModel()
    mu0H = 2.0
    peak_T_giant = giant_mce_analysis.find_peak_temperature(GD5SI2GE2_FIRST_ORDER, mu0H)
    peak_T_gd = giant_mce_analysis.find_peak_temperature(GADOLINIUM, mu0H)
    span = 10.0
    T_ashrae = 291.0

    def eval_at(material, T_cold, mass=5.0):
        sys_ = AMRSystem(material=material, mu0H_max=mu0H, mass_regenerator=mass,
                          frequency=1.0, fluid_mdot=0.08, loss_model=_loss,
                          use_ntu_thermal_model=True)
        return sys_.run(T_cold, span)

    r_gd_ashrae = eval_at(GADOLINIUM, T_ashrae)
    r_giant_ashrae = eval_at(GD5SI2GE2_FIRST_ORDER, T_ashrae)
    T_cold_favorable = peak_T_giant - span / 2
    r_giant_own = eval_at(GD5SI2GE2_FIRST_ORDER, T_cold_favorable)

    T_mid_ashrae = T_ashrae + span / 2.0
    tc_tuned = cascade._target_composition_for_peak(T_mid_ashrae, mu0H, cascade.GD_FAMILY)
    tuned_in_range = cascade.GD_FAMILY.tc_min <= tc_tuned <= cascade.GD_FAMILY.tc_max
    if tuned_in_range:
        material_tuned = cascade.GD_FAMILY.tuned_fn(tc_tuned)
        r_giant_tuned = eval_at(material_tuned, T_ashrae)
    else:
        r_giant_tuned = None

    labels = ['Gd\n@ ASHRAE (291K)', 'Gd5Si2Ge2\n@ ASHRAE (291K)', 'Gd5Si2Ge2\n@ own peak',
              'Gd5(SixGe1-x)4(-Ga)\n@ ASHRAE, tuned']
    colors4 = [COLOR_MAIN, COLOR_POWER, '#85bb65', '#c9a227']
    qc_vals = [r_gd_ashrae.Qc, r_giant_ashrae.Qc, r_giant_own.Qc,
               r_giant_tuned.Qc if r_giant_tuned else 0.0]
    cop_vals = [r_gd_ashrae.COP_electrical, r_giant_ashrae.COP_electrical, r_giant_own.COP_electrical,
                r_giant_tuned.COP_electrical if r_giant_tuned else 0.0]

    fig, axes = plt.subplots(1, 2, figsize=(12, 5.5))
    axes[0].bar(labels, qc_vals, color=colors4, alpha=0.85, edgecolor='white')
    axes[0].set_ylabel('Qc [W]')
    axes[0].set_title('Cooling Capacity')
    axes[1].bar(labels, cop_vals, color=colors4, alpha=0.85, edgecolor='white')
    axes[1].set_ylabel('Electrical COP')
    axes[1].set_title('Electrical COP')
    for ax in axes:
        ax.tick_params(axis='x', labelsize=7.5)

    fig.suptitle('Giant-MCE Targeting: Material Must Match the Operating Point\n'
                 f'(Gd5Si2Ge2 own peak: {peak_T_giant:.1f}K, Gd own peak: {peak_T_gd:.1f}K; '
                 f'tuned composition Tc={tc_tuned:.0f}K)',
                 fontsize=12)
    fig.tight_layout()
    save(fig, 'fig24_giant_mce_targeting_comparison')


# ══════════════════════════════════════════════════════════════════════════
# FIG 25 — Astronautics 6-layer Curie-graded La(Fe,Si)13Hy bed validation
# ══════════════════════════════════════════════════════════════════════════

def plot_astronautics_validation(precomputed=None):
    """precomputed, if given, may supply 'astro_result' (already computed
    by step 7c in main.py) so this figure reuses it instead of re-running
    the ~33s six-layer Curie-graded Astronautics validation a second time
    in the same pipeline invocation."""
    precomputed = precomputed or {}
    astro = precomputed.get('astro_result')
    if astro is None:
        astro = cascade.validate_astronautics_graded_bed()

    fig, axes = plt.subplots(1, 2, figsize=(12, 5.5))
    if astro.get('feasible'):
        stages = [s['stage'] for s in astro['stage_info']]
        tcs = [s['Tc_target_K'] for s in astro['stage_info']]
        tmids = [s['T_mid_K'] for s in astro['stage_info']]
        axes[0].plot(stages, tmids, 'o-', color=COLOR_MAIN, label='Stage T_mid')
        axes[0].plot(stages, tcs, 's--', color=COLOR_POWER, label='Needed composition Tc')
        axes[0].set_xlabel('Stage (coldest -> hottest)')
        axes[0].set_ylabel('Temperature [K]')
        axes[0].set_title('6-Layer Curie-Graded Bed:\nStage Temperature vs. Needed Composition')
        axes[0].legend(fontsize=9)

        labels = ['Qc [W]', 'COP']
        lit = [astro['Qc_lit_W'], astro['COP_lit']]
        model = [astro['Qc_W'], astro['COP_cascade']]
        x = np.arange(2)
        w = 0.35
        ax2 = axes[1]
        ax2b = ax2.twinx()
        ax2.bar(x[0] - w / 2, lit[0], w, color=COLOR_MAIN, alpha=0.85, edgecolor='white',
                label='Literature (Jacobs et al. 2014)')
        ax2.bar(x[0] + w / 2, model[0], w, color=COLOR_POWER, alpha=0.85, edgecolor='white',
                label='Model (this repo)')
        ax2b.bar(x[1] - w / 2, lit[1], w, color=COLOR_MAIN, alpha=0.85, edgecolor='white')
        ax2b.bar(x[1] + w / 2, model[1], w, color=COLOR_POWER, alpha=0.85, edgecolor='white')
        ax2.set_xticks(x)
        ax2.set_xticklabels(labels)
        ax2.set_ylabel('Qc [W]')
        ax2b.set_ylabel('COP')
        ax2.set_title(f"Astronautics_rotary_2014 Validation\nCOP error: {astro['COP_error_pct']:+.1f}%")
        # loc='best' on a twinned-axis bar chart doesn't see ax2b's bars at
        # all (matplotlib's auto-placement only considers artists on the
        # axis .legend() is called from), so it can and does pick a corner
        # that overlaps the COP bar on ax2b -- confirmed visually. Placed
        # explicitly below the axes instead, which is robust to whatever
        # the Qc/COP bar heights happen to be for any future device.
        handles = ax2.get_legend_handles_labels()[0]
        ax2.legend(handles, ['Literature (Jacobs et al. 2014)', 'Model (this repo)'],
                   fontsize=9, loc='upper center', bbox_to_anchor=(0.5, -0.12), ncol=2)
    else:
        for ax in axes:
            ax.text(0.5, 0.5, astro.get('status', 'infeasible'), ha='center', va='center',
                    wrap=True, fontsize=10, transform=ax.transAxes)
            ax.set_xticks([])
            ax.set_yticks([])

    fig.suptitle('6-Layer Curie-Graded La(Fe,Si)13Hy Bed vs. '
                 'Real Astronautics_rotary_2014 Device', fontsize=12)
    fig.tight_layout()
    save(fig, 'fig25_astronautics_graded_bed_validation')


# ══════════════════════════════════════════════════════════════════════════
# FIG 26 — Eight-way material family comparison (Track A2 item +  + )
# ══════════════════════════════════════════════════════════════════════════

def plot_material_family_comparison(precomputed=None):
    """precomputed, if given, may supply 'material_rows' (already computed
    by step 8d in main.py) so this figure reuses it instead of re-running
    the ~8s eight-way material-family sweep a second time in the same
    pipeline invocation."""
    precomputed = precomputed or {}
    rows = precomputed.get('material_rows')
    if rows is None:
        rows = material_family_comparison.build_comparison_table()
    rep = [r for r in rows if r['span_K'] == material_family_comparison.REPRESENTATIVE_SPAN_K]

    labels = [r['candidate'].replace(' (', '\n(') for r in rep]
    colors8 = [COLOR_MAIN, COLOR_POWER, '#85bb65', '#e8a33d', '#7b52ab', '#c2585d',
               '#2f9599', '#d64545']
    #  added a 6th candidate (nanocomposite);  added a
    # 7th (Ga1-xCMn3+x antiperovskite);  added an 8th (Mn1-xCuxCoGe).
    # colors8 must cover exactly as many entries as `rep` -- checked in
    # tests/test_plots.py::test_plot_material_family_comparison_color_count_matches_candidates.
    colors8 = colors8[:len(rep)] if len(rep) <= len(colors8) else (
        colors8 * (len(rep) // len(colors8) + 1))[:len(rep)]
    qc_vals = [r['1stage_Qc_W'] or 0.0 for r in rep]
    cop_vals = [r['1stage_COP'] or 0.0 for r in rep]
    hatches = ['' if r['in_range'] else '//' for r in rep]

    fig, axes = plt.subplots(1, 2, figsize=(12, 5.5))
    bars0 = axes[0].bar(labels, qc_vals, color=colors8, alpha=0.85, edgecolor='white')
    bars1 = axes[1].bar(labels, cop_vals, color=colors8, alpha=0.85, edgecolor='white')
    for bars in (bars0, bars1):
        for bar, hatch in zip(bars, hatches):
            bar.set_hatch(hatch)
    axes[0].set_ylabel('Qc [W]')
    axes[0].set_title('Cooling Capacity (1-stage)')
    axes[1].set_ylabel('Electrical COP')
    axes[1].set_title('Electrical COP (1-stage)')
    for ax in axes:
        ax.tick_params(axis='x', labelsize=7)
        ax.set_xticklabels(labels, rotation=30, ha='right', fontsize=7)

    span = material_family_comparison.REPRESENTATIVE_SPAN_K
    T_cold = material_family_comparison.T_COLD_K
    fig.suptitle('Material Family Comparison at the ASHRAE Point '
                 f'(T_cold={T_cold:.0f}K, span={span:.0f}K)\n'
                 'Hatched bars: family\'s documented Tc window does not cover this point '
                 '(fell back to Gd)', fontsize=11)
    fig.tight_layout()
    save(fig, 'fig26_material_family_comparison')


# ══════════════════════════════════════════════════════════════════════════
# FIG 27 — Inhomogeneous/polycrystalline Tc-broadening sensitivity
# ══════════════════════════════════════════════════════════════════════════

def plot_inhomogeneous_broadening():
    broad_rows = inhomogeneous_broadening.run_broadening_sweep(verbose=False)
    err_rows = inhomogeneous_broadening.run_dankov_error_sensitivity(verbose=False)

    sigmas = sorted(set(r['sigma_Tc_K'] for r in broad_rows))
    fields = sorted(set(r['mu0H_T'] for r in broad_rows))

    fig, axes = plt.subplots(1, 2, figsize=(12, 5))
    for B, color in zip(fields, COLOR_CYCLE):
        fwhm = [r['fwhm_K'] for r in broad_rows if r['mu0H_T'] == B]
        axes[0].plot(sigmas, fwhm, marker='o', color=color, label=f'{B:.1f} T')
    axes[0].set_xlabel(r'$\sigma_{Tc}$ [K]  (Gaussian Tc spread)')
    axes[0].set_ylabel(r'$\Delta T_{ad}(T)$ FWHM [K]')
    axes[0].set_title('Peak Broadening vs. Tc Inhomogeneity')
    axes[0].legend(title='Field', fontsize=8)

    worst = {}
    for r in err_rows:
        s = r['sigma_Tc_K']
        worst[s] = max(worst.get(s, 0.0), abs(r['err_pct']))
    xs = sorted(worst)
    ys = [worst[s] for s in xs]
    axes[1].plot(xs, ys, marker='s', color=COLOR_POWER)
    axes[1].set_xlabel(r'$\sigma_{Tc}$ [K]')
    axes[1].set_ylabel("Worst-field |error| vs. Dan'kov et al. (1998) [%]")
    axes[1].set_title('Validation Error vs. Tc-Broadening')

    fig.suptitle('Gaussian Inhomogeneous/Polycrystalline Tc-Broadening Sensitivity\n'
                 '(standard literature broadening treatment, '
                 'not digitized book content)', fontsize=11)
    fig.tight_layout()
    save(fig, 'fig27_inhomogeneous_tc_broadening')


# ══════════════════════════════════════════════════════════════════════════
# FIG 28 — Nanocomposite off-design robustness check ( follow-up)
# ══════════════════════════════════════════════════════════════════════════

def plot_nanocomposite_robustness():
    result = nanocomposite_material.run_robustness_check(
        out_path=str(RESULTS_DIR / 'nanocomposite_robustness.txt'), verbose=False)
    rows = result['rows']
    spans = [r['span_K'] for r in rows]
    nano_qc = [r['nanocomposite_Qc_W'] for r in rows]
    single_qc = [r['single_phase_Qc_W'] for r in rows]
    design = [r['is_design_span'] for r in rows]

    fig, ax = plt.subplots(figsize=(8.5, 5.5))
    x = np.arange(len(spans))
    w = 0.35
    ax.bar(x - w / 2, nano_qc, w, label='Nanocomposite (3-phase blend)',
           color=COLOR_MAIN, alpha=0.85, edgecolor='white')
    ax.bar(x + w / 2, single_qc, w, label='Single phase (sharply tuned)',
           color=COLOR_POWER, alpha=0.85, edgecolor='white')
    labels = [f'{s:.0f}K' + ('\n(design)' if d else '') for s, d in zip(spans, design)]
    ax.set_xticks(x)
    ax.set_xticklabels(labels)
    ax.set_xlabel('Evaluated span [K]  (both tuned once at the 10K design span)')
    ax.set_ylabel('Qc [W]   (0 = infeasible at this span)')
    ax.set_title('Nanocomposite Off-Design Robustness\n'
                 'Broad blend keeps Qc>0 where a sharply-tuned single phase collapses')

    # A zero-height bar is otherwise indistinguishable from a missing data
    # point -- mark genuinely-infeasible (Qc=0) bars explicitly.
    y_span = max(max(nano_qc, default=0), max(single_qc, default=0)) or 1.0
    for xi, v in zip(x - w / 2, nano_qc):
        if v == 0:
            ax.text(xi, 0.015 * y_span, '\u2715 infeasible', ha='center', va='bottom',
                    fontsize=7.5, color='gray', rotation=90)
    for xi, v in zip(x + w / 2, single_qc):
        if v == 0:
            ax.text(xi, 0.015 * y_span, '\u2715 infeasible', ha='center', va='bottom',
                    fontsize=7.5, color='gray', rotation=90)

    ax.legend()
    fig.tight_layout()
    save(fig, 'fig28_nanocomposite_offdesign_robustness')


# ══════════════════════════════════════════════════════════════════════════
# FIG 29 — Mechanical-contact thermal-diode actuation-cost sensitivity
# ══════════════════════════════════════════════════════════════════════════

def plot_thermal_diode_sensitivity():
    rows = thermal_diode_analysis.sweep_frequency_with_and_without_diode(verbose=False)
    freqs = [r[0] for r in rows]
    cop_no_diode = [r[1] for r in rows]
    cop_diode = [r[2] for r in rows]
    delta_pct = [r[3] for r in rows]

    fig, axes = plt.subplots(1, 2, figsize=(12, 5))
    axes[0].plot(freqs, cop_no_diode, marker='o', color=COLOR_MAIN, label='No diode')
    axes[0].plot(freqs, cop_diode, marker='s', color=COLOR_POWER, linestyle='--',
                 label='Diode-assisted (actuation cost only)')
    axes[0].set_xlabel('Frequency [Hz]')
    axes[0].set_ylabel('COP_electrical')
    axes[0].set_title('COP vs. Frequency, With/Without Diode')
    axes[0].legend(fontsize=8)

    axes[1].bar([f'{f:.2f}Hz' for f in freqs], delta_pct, color=COLOR_POWER, alpha=0.85,
                edgecolor='white')
    axes[1].axhline(0, color='k', linewidth=0.8)
    axes[1].set_ylabel('COP_electrical change [%]')
    axes[1].set_title('Diode Actuation-Cost Penalty\n(cost-only accounting, no rectification benefit modeled)')

    fig.suptitle('Mechanical-Contact Active Thermal Diode — Sensitivity Study ',
                 fontsize=12)
    fig.tight_layout()
    save(fig, 'fig29_thermal_diode_sensitivity')


# ══════════════════════════════════════════════════════════════════════════
# FIG 30 — Magnetocaloric-fluid (ferrofluid/MR suspension) volume-fraction sweep
# ══════════════════════════════════════════════════════════════════════════

def plot_fluid_mce_sweep():
    sweep = fluid_mce_analysis.volume_fraction_sweep()
    rows = sweep['rows']
    best = sweep['best_row']
    phi = [r['phi'] for r in rows]
    cop = [r['COP_electrical'] for r in rows]
    qc = [r['Qc_W'] for r in rows]

    fig, ax1 = plt.subplots(figsize=(8, 5.5))
    ax2 = ax1.twinx()
    ax1.plot(phi, cop, color=COLOR_MAIN, marker='o', label='COP_electrical')
    ax2.plot(phi, qc, color=COLOR_POWER, marker='s', linestyle='--', label='Qc')
    ax1.axvline(best['phi'], color=COLOR_MAIN, linestyle=':', linewidth=1)
    # Offset well clear of the curve (was sitting directly on the peak,
    # colliding with the line) -- push down-and-right with a visible arrow
    # so the label and the data stay legible independently.
    ax1.annotate(f"interior optimum\nphi={best['phi']:.2f}", xy=(best['phi'], best['COP_electrical']),
                 xytext=(45, -35), textcoords='offset points', fontsize=8, color=COLOR_MAIN,
                 ha='left',
                 arrowprops=dict(arrowstyle='->', color=COLOR_MAIN, lw=1))
    ax1.set_xlabel('Particle volume fraction $\\phi$')
    ax1.set_ylabel('COP_electrical', color=COLOR_MAIN)
    ax2.set_ylabel('Qc [W]  (each $\\phi$ at its own favorable span)', color=COLOR_POWER)
    ax1.set_title('Ferrofluid/MR-Suspension MCE — Volume-Fraction Trade-off \n'
                 'Viscosity-vs-MCE-intensity trade-off produces a genuine interior optimum')
    fig.tight_layout()
    save(fig, 'fig30_fluid_mce_volume_fraction')


# ══════════════════════════════════════════════════════════════════════════
# FIG 31 — Passive/hybrid magnetic-regenerator augmentation alignment effect
# ══════════════════════════════════════════════════════════════════════════

def plot_passive_regenerator_alignment():
    base, results = passive_regenerator_analysis.compare_candidate_materials(verbose=False)
    names = [r.material_name for r in results]
    gain_pct = [r.cop_gain_fraction * 100 for r in results]

    span_rows = passive_regenerator_analysis.span_sweep(verbose=False)
    spans = [r['span_K'] for r in span_rows]
    best_gain = [r['best_cop_gain_fraction'] * 100 for r in span_rows]

    fig, axes = plt.subplots(1, 2, figsize=(12, 5))
    axes[0].bar(names, gain_pct, color=COLOR_MAIN, alpha=0.85, edgecolor='white')
    axes[0].set_ylabel('Augmented COP gain [%]')
    axes[0].set_title('Candidate-Material Alignment\n(T_cold=291.15K, span=10K)')
    axes[0].tick_params(axis='x', labelsize=7, rotation=20)

    # Zero-gain candidates (no usable alignment at this operating point)
    # otherwise render as an empty gap, indistinguishable from missing data.
    y_span = max(gain_pct, default=0) or 1.0
    for xi, v in enumerate(gain_pct):
        if v == 0:
            axes[0].text(xi, 0.015 * y_span, '\u2715 no alignment', ha='center', va='bottom',
                         fontsize=7, color='gray', rotation=90)

    axes[1].plot(spans, best_gain, marker='o', color=COLOR_POWER)
    axes[1].set_xlabel('Span [K]')
    axes[1].set_ylabel('Best-material COP gain [%]')
    axes[1].set_title('Gain vs. Span\n(alignment window widens/narrows with span)')

    fig.suptitle('Passive Magnetic-Regenerator Augmentation of a Vapor-Compression Cycle '
                 '\nIllustrative literature-range ceiling, not a fitted coefficient',
                 fontsize=11)
    fig.tight_layout()
    save(fig, 'fig31_passive_regenerator_alignment')


# ══════════════════════════════════════════════════════════════════════════
# FIG 32 — Rotary-device cycle-type (Ericsson- vs. Brayton-like) validation
# ══════════════════════════════════════════════════════════════════════════

def plot_cycle_type_validation():
    results = validation_system.run_cycle_type_validation(verbose=False, out_path=None)
    comparable = [r for r in results if 'COP_error_pct_baseline_brayton' in r]

    fig, ax = plt.subplots(figsize=(7.5, 5.5))
    if comparable:
        names = [r['device'] for r in comparable]
        brayton_err = [abs(r['COP_error_pct_baseline_brayton']) for r in comparable]
        inferred_err = [abs(r['COP_error_pct_cycle_inferred']) for r in comparable]
        x = np.arange(len(names))
        w = 0.35
        ax.bar(x - w / 2, brayton_err, w, label='brayton (default)',
               color=COLOR_MAIN, alpha=0.85, edgecolor='white')
        ax.bar(x + w / 2, inferred_err, w, label='ericsson (rotary-inferred)',
               color=COLOR_POWER, alpha=0.85, edgecolor='white')
        ax.set_xticks(x)
        ax.set_xticklabels(names, rotation=15, ha='right', fontsize=8)
        # With only n=1 comparable device (the common case here), x has a
        # single tick at 0 and matplotlib's default 5%-of-data-range
        # margin barely pads beyond the bars themselves (width 0.35 each,
        # so the data range is only 0.35 wide) -- the two bars end up
        # stretched to fill almost the entire axes edge-to-edge. Pad the
        # x-limits explicitly so this chart keeps the same breathing room
        # a multi-device version of it would have.
        ax.set_xlim(x[0] - 1.0, x[-1] + 1.0)
    ax.set_ylabel('|COP error| vs. published value [%]')
    ax.set_title('Cycle-Topology Sensitivity: Rotary Benchmark Devices\n'
                 '(naming-convention proxy, not a literature-confirmed classification)')
    ax.legend()
    fig.tight_layout()
    save(fig, 'fig32_cycle_type_validation')


# ══════════════════════════════════════════════════════════════════════════
# FIG 33 — Hysteresis-loss Pareto-front material-selection sensitivity
# ══════════════════════════════════════════════════════════════════════════

def plot_hysteresis_sensitivity(precomputed=None):
    """precomputed, if given, may supply 'hysteresis_result' (the dict
    already returned by step 11b's run_hysteresis_sensitivity() call in
    main.py) so this figure reuses that ON/OFF A/B comparison instead of
    re-running the ~9s NSGA-III optimization a second time in the same
    pipeline invocation."""
    precomputed = precomputed or {}
    result = precomputed.get('hysteresis_result')
    if result is not None:
        counts_on, counts_off = result['counts_on'], result['counts_off']
    elif HAVE_PYMOO:
        result = hysteresis_sensitivity.run_hysteresis_sensitivity(verbose=False)
        counts_on, counts_off = result['counts_on'], result['counts_off']
    else:
        print("  [pymoo unavailable — falling back to pre-computed "
              "results/pareto_front_hysteresis_{on,off}.csv]")
        rows_on = _read_csv_rows(RESULTS_DIR / 'pareto_front_hysteresis_on.csv')
        rows_off = _read_csv_rows(RESULTS_DIR / 'pareto_front_hysteresis_off.csv')
        counts_on = hysteresis_sensitivity._material_counts(rows_on)
        counts_off = hysteresis_sensitivity._material_counts(rows_off)

    materials = sorted(set(counts_on) | set(counts_off))
    n_on = [counts_on.get(m, 0) for m in materials]
    n_off = [counts_off.get(m, 0) for m in materials]

    fig, ax = plt.subplots(figsize=(8.5, 5.5))
    x = np.arange(len(materials))
    w = 0.35
    ax.bar(x - w / 2, n_off, w, label='Hysteresis OFF (previous)',
           color=COLOR_MAIN, alpha=0.85, edgecolor='white')
    ax.bar(x + w / 2, n_on, w, label='Hysteresis ON (literature-placeholder)',
           color=COLOR_POWER, alpha=0.85, edgecolor='white')
    ax.set_xticks(x)
    ax.set_xticklabels([m.replace(' (', '\n(') for m in materials], fontsize=7)
    ax.set_ylabel('Designs on merged Pareto front')
    ax.set_title('Thermal-Hysteresis Loss: Material-Selection Sensitivity \n'
                 'A/B NSGA-III comparison, identical pop_size/n_gen/seed')
    ax.legend(fontsize=8)
    fig.tight_layout()
    save(fig, 'fig33_hysteresis_pareto_sensitivity')


# ══════════════════════════════════════════════════════════════════════════
# FIG 34 — Halbach-cylinder magnet-mass Pareto-front sensitivity
# ══════════════════════════════════════════════════════════════════════════

def plot_magnet_geometry_pareto_sensitivity(precomputed=None):
    """precomputed, if given, may supply 'magnet_geometry_result' (the
    dict already returned by step 11d's run_geometric_cost_pareto_
    sensitivity() call in main.py) so this figure reuses that FLAT-vs-
    GEOMETRIC A/B comparison instead of re-running two full NSGA-III
    optimizations (~10s) a second time in the same pipeline invocation."""
    precomputed = precomputed or {}
    result = precomputed.get('magnet_geometry_result')
    if result is not None:
        rows_flat = result['rows_flat']
        rows_geom = result['rows_geometric']
    elif HAVE_PYMOO:
        rows_flat = optimize_mod.run_optimization(
            out_csv=str(RESULTS_DIR / 'pareto_front_magnet_flat.csv'),
            per_material_out_dir=None, use_geometric_magnet_mass=False)
        rows_geom = optimize_mod.run_optimization(
            out_csv=str(RESULTS_DIR / 'pareto_front_magnet_geometric.csv'),
            per_material_out_dir=None, use_geometric_magnet_mass=True)
    else:
        print("  [pymoo unavailable — falling back to pre-computed "
              "results/pareto_front_magnet_{flat,geometric}.csv]")
        rows_flat = _read_csv_rows(RESULTS_DIR / 'pareto_front_magnet_flat.csv')
        rows_geom = _read_csv_rows(RESULTS_DIR / 'pareto_front_magnet_geometric.csv')

    fields_flat = [r['mu0H_max_T'] for r in rows_flat]
    fields_geom = [r['mu0H_max_T'] for r in rows_geom]

    fig, axes = plt.subplots(1, 2, figsize=(12, 5))
    axes[0].hist(fields_flat, bins=8, alpha=0.6, color=COLOR_MAIN, label='FLAT (previous)')
    axes[0].hist(fields_geom, bins=8, alpha=0.6, color=COLOR_POWER, label='GEOMETRIC (Halbach)')
    axes[0].set_xlabel('mu0H_max [T]')
    axes[0].set_ylabel('Designs on merged Pareto front')
    axes[0].set_title('Field-Strength Distribution\nFLAT vs. GEOMETRIC magnet-mass cost term')
    axes[0].legend(fontsize=8)

    counts_flat = magnet_geometry._material_counts(rows_flat)
    counts_geom = magnet_geometry._material_counts(rows_geom)
    materials = sorted(set(counts_flat) | set(counts_geom))
    x = np.arange(len(materials))
    w = 0.35
    axes[1].bar(x - w / 2, [counts_flat.get(m, 0) for m in materials], w,
                label='FLAT', color=COLOR_MAIN, alpha=0.85, edgecolor='white')
    axes[1].bar(x + w / 2, [counts_geom.get(m, 0) for m in materials], w,
                label='GEOMETRIC', color=COLOR_POWER, alpha=0.85, edgecolor='white')
    axes[1].set_xticks(x)
    axes[1].set_xticklabels([m.replace(' (', '\n(') for m in materials], fontsize=7)
    axes[1].set_ylabel('Designs on merged Pareto front')
    axes[1].set_title('Material-Selection Sensitivity')
    axes[1].legend(fontsize=8)

    fig.suptitle('Halbach-Cylinder Magnet-Mass Cost Term: Pareto-Front Sensitivity ',
                 fontsize=12)
    fig.tight_layout()
    save(fig, 'fig34_magnet_geometry_pareto_sensitivity')


# ══════════════════════════════════════════════════════════════════════════
# FIG 35 — Regenerative-amplification OVERRIDE check (this session)
# ══════════════════════════════════════════════════════════════════════════

def plot_regenerative_amplification_override_check(precomputed=None):
    """Visualizes validation_system.run_regenerative_amplification_
    override_check()'s result: for each benchmark device the 2d/2e
    diagnostics flag as STRUCTURALLY infeasible under
    core.amr_cycle.AMRSystem's default 2*dTad_noload span cap (i.e.
    cooling_capacity() returns a hard Qc=0 at that device's real
    reported span, for ANY mdot), does the opt-in
    `no_load_span_override` -- populated from
    core.regenerator_1d.regenerative_span_cap(), a real multi-cycle
    transient simulation -- recover a usable, evaluable COP prediction,
    and how close does it land?

    This is the first figure for this session's regenerative-
    amplification work (core/regenerator_1d.py, AMRSystem's
    no_load_span_override, and this check) -- everything else from that
    work (results/regenerator_1d_validation.txt,
    results/regenerative_amplification_override_check.txt) was
    previously text-only. Two panels:
      left  -- span reach per device: the old hard cap vs. the 1-D
               model's own span cap vs. the device's actual reported
               span (bar at the real span shows whether the override
               even reaches it).
      right -- for devices the override DOES reach, predicted vs.
               literature COP (a Risoe-DTU-style "doesn't help at all"
               device, if present, is called out by name rather than
               silently omitted -- see run_regenerative_amplification_
               override_check()'s own docstring for why one device can
               show this).

    precomputed, if given, may supply 'override_check_result' (the list
    of dicts already returned by step 2f's
    run_regenerative_amplification_override_check(max_devices=3) call in
    main.py) so this figure reuses it instead of re-running that check
    (each device is a multi-mdot transient search, tens of seconds) a
    second time in the same pipeline invocation. When called standalone
    with no precomputed data, this recomputes the same bounded
    (max_devices=3) check main.py itself runs, for a comparable default
    runtime -- pass core.validation_system.run_regenerative_amplification_
    override_check(max_devices=None) results in via `precomputed` for the
    full, every-flagged-device version instead."""
    precomputed = precomputed or {}
    results = precomputed.get('override_check_result')
    if results is None:
        results = validation_system.run_regenerative_amplification_override_check(
            verbose=False, max_devices=3)

    if not results:
        fig, ax = plt.subplots(figsize=(8, 4))
        ax.text(0.5, 0.5, 'No structurally-infeasible, COP-bearing benchmark rows found',
                ha='center', va='center', transform=ax.transAxes)
        ax.set_xticks([])
        ax.set_yticks([])
        fig.tight_layout()
        save(fig, 'fig35_regenerative_amplification_override_check')
        return

    devices = [r['device'] for r in results]
    n = len(devices)
    x = np.arange(n)

    fig, axes = plt.subplots(1, 2, figsize=(13, 5.5))

    w = 0.25
    axes[0].bar(x - w, [r['old_cap_K'] for r in results], w,
                label='Old cap (2*dTad_noload)', color=COLOR_MAIN, alpha=0.85, edgecolor='white')
    axes[0].bar(x, [r['span_cap_K'] for r in results], w,
                label='1-D span cap (this session)', color=COLOR_POWER, alpha=0.85, edgecolor='white')
    axes[0].bar(x + w, [r['span_K'] for r in results], w,
                label='Actual reported span', color='#2ca02c', alpha=0.85, edgecolor='white')
    axes[0].set_xticks(x)
    axes[0].set_xticklabels([d.replace('_', '\n') for d in devices], fontsize=7)
    axes[0].set_ylabel('Span [K]')
    axes[0].set_title('Span Reach: Old Cap vs. 1-D Model vs. Reality\n'
                       '(bar at "actual span" not reached by the 1-D cap = still infeasible)')
    axes[0].legend(fontsize=8)

    recovered = [r for r in results if r['recovers_nonzero'] and r['COP_pred'] is not None]
    unrecovered = [r for r in results if not (r['recovers_nonzero'] and r['COP_pred'] is not None)]
    if recovered:
        rx = np.arange(len(recovered))
        rw = 0.35
        axes[1].bar(rx - rw / 2, [r['COP_lit'] for r in recovered], rw,
                    label='Literature COP', color=COLOR_MAIN, alpha=0.85, edgecolor='white')
        axes[1].bar(rx + rw / 2, [r['COP_pred'] for r in recovered], rw,
                    label='Model COP (with override)', color=COLOR_POWER, alpha=0.85, edgecolor='white')
        for i, r in enumerate(recovered):
            axes[1].annotate(f"{r['err_pct']:+.1f}%", xy=(i + rw / 2, r['COP_pred']),
                              xytext=(0, 4), textcoords='offset points', ha='center', fontsize=8)
        axes[1].set_xticks(rx)
        axes[1].set_xticklabels([r['device'].replace('_', '\n') for r in recovered], fontsize=7)
        axes[1].set_ylabel('COP')
        axes[1].legend(fontsize=8)
    title = 'COP Recovered by the Override\n(% label = model error vs. literature)'
    if unrecovered:
        title += '\nStill infeasible: ' + ', '.join(r['device'] for r in unrecovered)
    axes[1].set_title(title, fontsize=10)

    fig.suptitle('Regenerative-Amplification Override Check: Does no_load_span_override\n'
                 'Recover a Usable Prediction Where the Old Cap Gives a Hard Qc=0?', fontsize=12)
    fig.tight_layout()
    save(fig, 'fig35_regenerative_amplification_override_check')


# ══════════════════════════════════════════════════════════════════════════
# FIG 36 — Value Proposition: COP crossover, AMR vs VCC across spans
# ══════════════════════════════════════════════════════════════════════════

def plot_value_proposition_cop_crossover():
    """First figure for docs/MCE_Value_Proposition.md's Section 1 (the COP
    case does not exist) -- that section previously had only a markdown
    table. Reproduces run_material_and_field_crossover_search()'s own
    4-span result directly (hardcoded from that module's own printed
    table, not re-run here, to keep this figure cheap -- see that
    module's own run function for the live search) rather than
    re-deriving a new number: at every span, AMR's best found COP (best
    of an 1-7T / multi-material / cascade grid) still trails VCC's, and
    by a WIDER margin as span shrinks, not a narrower one."""
    spans = [5.0, 10.0, 15.0, 20.0]
    amr_cop = [9.70, 7.05, 5.43, 4.32]
    vcc_cop = [24.36, 11.97, 7.84, 5.77]

    fig, ax = plt.subplots(figsize=(8, 5.5))
    ax.plot(spans, amr_cop, 'o-', color=COLOR_POWER, label='Best AMR COP found\n'
            '(La(Fe,Si)13Hy, 1-7T grid, cascades)')
    ax.plot(spans, vcc_cop, 's-', color=COLOR_MAIN, label='VCC COP (eta=0.42)')
    ax.fill_between(spans, amr_cop, vcc_cop, color=COLOR_MAIN, alpha=0.08)
    for s, a, v in zip(spans, amr_cop, vcc_cop):
        ax.annotate(f'{v/a:.1f}x', xy=(s, (a + v) / 2), ha='center', fontsize=9,
                    color='#555555')
    ax.set_xlabel('Temperature Span [K]')
    ax.set_ylabel('COP')
    ax.set_title('No Crossover Found: AMR Never Beats VCC\n'
                  '(labels show VCC/AMR ratio -- widens, not narrows, at small spans)')
    ax.legend(fontsize=9)
    fig.tight_layout()
    save(fig, 'fig36_value_proposition_cop_crossover')


# ══════════════════════════════════════════════════════════════════════════
# FIG 37 — Value Proposition: emissions breakdown, AMR vs VCC vs liquid
# ══════════════════════════════════════════════════════════════════════════

def plot_value_proposition_emissions_breakdown():
    """Second figure for MCE_Value_Proposition.md's Section 2 (refrigerant
    elimination does not rescue emissions at this COP gap) -- reuses
    core.emissions.compare_emissions() at the SAME Magnotherm Eclipse
    operating point (AMR COP=1.76, VCC COP=6.66) the doc's own table
    already cites, rather than a fresh/different call, so this figure and
    that table cannot silently drift apart."""
    results = emissions.compare_emissions(0.4, amr_cop=1.76, vcc_cop=6.66, liquid_cop=9.0)
    names = [r.technology for r in results]
    refrig = [r.refrigerant_GWP_tCO2e_per_year for r in results]
    op = [r.operational_CO2_tCO2e_per_year for r in results]
    total = [r + o for r, o in zip(refrig, op)]

    fig, ax = plt.subplots(figsize=(8.5, 5.5))
    x = np.arange(len(names))
    ax.bar(x, refrig, label='Refrigerant leakage', color='#e07b54', alpha=0.85, edgecolor='white')
    ax.bar(x, op, bottom=refrig, label='Operational (electricity)', color=COLOR_MAIN,
           alpha=0.85, edgecolor='white')
    for i, t in enumerate(total):
        ax.annotate(f'{t:.3f}', xy=(i, t), xytext=(0, 4), textcoords='offset points',
                    ha='center', fontsize=9, fontweight='bold')
    ax.set_xticks(x)
    ax.set_xticklabels(names, fontsize=8.5, rotation=10, ha='right')
    ax.set_ylabel(r'tCO$_2$e / year')
    ax.set_title('Emissions at the Magnotherm Eclipse Operating Point (0.4kW)\n'
                 'AMR comes out 3.6x higher total, driven by the operational (COP) term,\n'
                 'not the eliminated refrigerant term')
    ax.legend(fontsize=9)
    fig.tight_layout()
    save(fig, 'fig37_value_proposition_emissions_breakdown')


# ══════════════════════════════════════════════════════════════════════════
# FIG 38 — Value Proposition: refrigerant GWP landscape
# ══════════════════════════════════════════════════════════════════════════

def plot_value_proposition_refrigerant_gwp():
    """Third figure for the Value Proposition doc -- the regulatory-parity
    argument (Section 4) rests on refrigerant GWP figures already in
    core/emissions.py's own REFRIGERANT_GWP dict; this puts them on one
    axis, with AMR's own zero placed for direct visual contrast, rather
    than leaving the comparison to a sentence."""
    gwp = dict(emissions.REFRIGERANT_GWP)
    names = list(gwp.keys()) + ['Magnetocaloric\n(AMR, no refrigerant)']
    values = list(gwp.values()) + [0]
    colors = list(plt.cm.YlOrRd(np.linspace(0.35, 0.85, len(gwp)))) + ['#2ca02c']

    fig, ax = plt.subplots(figsize=(8, 5.5))
    bars = ax.bar(names, values, color=colors, alpha=0.9, edgecolor='white')
    for b, v in zip(bars, values):
        ax.annotate(f'{v:g}', xy=(b.get_x() + b.get_width() / 2, v), xytext=(0, 4),
                    textcoords='offset points', ha='center', fontsize=9)
    ax.set_ylabel('100-year GWP (IPCC AR5)')
    ax.set_title('Refrigerant GWP Landscape\n'
                 '(the regulatory-tailwind argument: real, but does not by itself\n'
                 'close the COP/emissions gap shown in fig37)')
    ax.tick_params(axis='x', labelsize=8.5)
    fig.tight_layout()
    save(fig, 'fig38_value_proposition_refrigerant_gwp')


# ══════════════════════════════════════════════════════════════════════════
# FIG 39 — Value Proposition: real deployments, what they actually show
# ══════════════════════════════════════════════════════════════════════════

def plot_value_proposition_real_deployments():
    """Fourth figure for the Value Proposition doc's Section 3 (what the
    real, deployed systems actually show) -- plots Magnotherm Eclipse's
    reported real-world energy SAVING (a relative, same-duty comparison,
    not a COP) alongside Polaris's directly reported plug-in COP and
    second-law efficiency (from core/beverage_cooler_validation.py's own
    ECLIPSE_REPORTED_ENERGY_SAVING_PCT / POLARIS_PLUGIN_COP /
    POLARIS_SECOND_LAW_EFF_PCT constants -- not re-typed numbers), next
    to VCC's own COP for scale. Two different metrics deliberately kept
    on two panels rather than forced onto one axis, since a plug-in COP
    of 1.0 and a 15% energy saving are not directly comparable
    quantities -- conflating them into one bar chart would misstate what
    either device actually demonstrated."""
    left_labels = ['Magnotherm Eclipse\nvs. incumbent R290\n(same duty)']
    left_vals = [beverage_cooler_validation.ECLIPSE_REPORTED_ENERGY_SAVING_PCT]

    right_labels = ['Polaris\nplug-in COP', 'Reference VCC\nCOP (eta=0.42,\nsame 15K span)']
    polaris_cop = beverage_cooler_validation.POLARIS_PLUGIN_COP
    T_cold_K = beverage_cooler_validation.ECLIPSE_T_COLD_C + 273.15
    T_hot_K = T_cold_K + beverage_cooler_validation.POLARIS_SPAN_K
    vcc_ref_cop = vapor_compression_cop(T_cold_K, T_hot_K, eta_2nd_law=0.42).COP
    right_vals = [polaris_cop, vcc_ref_cop]

    fig, axes = plt.subplots(1, 2, figsize=(11, 5.5))
    axes[0].bar(left_labels, left_vals, color='#2ca02c', alpha=0.85, edgecolor='white', width=0.5)
    axes[0].set_ylabel('Real-world energy saving [%]')
    axes[0].set_title('Eclipse / REWE Pilot\n(11-week in-store, press-reported)')
    axes[0].annotate(f'{left_vals[0]:.0f}%', xy=(0, left_vals[0]), xytext=(0, 4),
                      textcoords='offset points', ha='center', fontsize=10, fontweight='bold')

    axes[1].bar(right_labels, right_vals, color=[COLOR_POWER, COLOR_MAIN], alpha=0.85,
                edgecolor='white', width=0.5)
    axes[1].set_ylabel('COP')
    axes[1].set_title(f'Polaris (peer-reviewed)\nsecond-law eff.='
                       f'{beverage_cooler_validation.POLARIS_SECOND_LAW_EFF_PCT:.1f}% '
                       f'at {beverage_cooler_validation.POLARIS_FIELD_T}T/'
                       f'{beverage_cooler_validation.POLARIS_SPAN_K:.0f}K')
    for i, v in enumerate(right_vals):
        axes[1].annotate(f'{v:.2f}', xy=(i, v), xytext=(0, 4), textcoords='offset points',
                          ha='center', fontsize=10, fontweight='bold')

    fig.suptitle('Real Deployed Magnetocaloric Systems: Neither Claims COP Superiority',
                 fontsize=12)
    fig.tight_layout()
    save(fig, 'fig39_value_proposition_real_deployments')


# ══════════════════════════════════════════════════════════════════════════
# FIG 40 — Doped (Mn,Fe)2(P,Si) hysteresis-loss speculative estimate
# ══════════════════════════════════════════════════════════════════════════

def plot_mnfepsi_doped_hysteresis_speculative():
    """First figure for Phase 36's mnfepsi_doped_hysteresis_speculative.py
    -- puts the two V-doped estimates (mean-k and least-squares-k) next
    to MNFEPSI_FIRST_ORDER's own existing 25.0 J/kg placeholder on a log
    axis (the estimates are ~30-40x smaller, which a linear axis would
    flatten to invisible bars). Explicitly titled/labeled as
    SPECULATIVE, matching the module's own honesty-flag convention --
    this is not a replacement default."""
    from core.mnfepsi_doped_hysteresis_speculative import VONFE_ESTIMATE, VONMN_ESTIMATE
    from core.first_order_mce import MNFEPSI_FIRST_ORDER

    labels = ['MNFEPSI_FIRST_ORDER\n(current placeholder)',
              'V-on-Fe doped\n(mean-k, SPECULATIVE)', 'V-on-Fe doped\n(least-sq-k, SPECULATIVE)',
              'V-on-Mn doped\n(mean-k, SPECULATIVE)', 'V-on-Mn doped\n(least-sq-k, SPECULATIVE)']
    values = [MNFEPSI_FIRST_ORDER.hysteresis_loss_J_per_kg,
              VONFE_ESTIMATE.estimate_J_per_kg_mean_k, VONFE_ESTIMATE.estimate_J_per_kg_ls_k,
              VONMN_ESTIMATE.estimate_J_per_kg_mean_k, VONMN_ESTIMATE.estimate_J_per_kg_ls_k]
    colors = [COLOR_MAIN, '#e07b54', '#e07b54', '#c9a227', '#c9a227']

    fig, ax = plt.subplots(figsize=(9, 5.5))
    bars = ax.bar(labels, values, color=colors, alpha=0.85, edgecolor='white')
    ax.set_yscale('log')
    for b, v in zip(bars, values):
        ax.annotate(f'{v:.2f}', xy=(b.get_x() + b.get_width() / 2, v), xytext=(0, 4),
                    textcoords='offset points', ha='center', fontsize=9)
    ax.set_ylabel('Hysteresis loss [J/kg] (log scale)')
    ax.set_title('SPECULATIVE: V-Doped (Mn,Fe)2(P,Si) Hysteresis-Loss Estimate\n'
                 'vs. Current 25.0 J/kg Placeholder (Phase 36, not a new default)')
    ax.tick_params(axis='x', labelsize=8)
    fig.tight_layout()
    save(fig, 'fig40_mnfepsi_doped_hysteresis_speculative')


# ══════════════════════════════════════════════════════════════════════════
# FIG 41 — Zhang et al. k-fit quality (proxy for the estimate above)
# ══════════════════════════════════════════════════════════════════════════

def plot_mnfepsi_hysteresis_kfit_quality():
    """Second figure for the same Phase 36 module -- shows the actual fit
    quality of k=W_hys/(dS*T_hys) against Zhang et al.'s own 5-point
    table (module docstring), the thing VONFE_ESTIMATE/VONMN_ESTIMATE are
    extrapolated from. Plots measured W_hys against k_mean*dS*T_hys for
    each of the 5 points, with a y=x reference line, so the ~20-40%
    scatter the module's own docstring quotes is visible directly rather
    than only stated as a number."""
    from core.mnfepsi_doped_hysteresis_speculative import _ZHANG_TABLE, K_MEAN

    xs = [ds * t for (_x, t, _w, ds) in _ZHANG_TABLE]
    ys_actual = [w for (_x, _t, w, _ds) in _ZHANG_TABLE]
    ys_fit = [K_MEAN * x for x in xs]
    comp_labels = [f'x={x_comp:.1f}' for (x_comp, _t, _w, _ds) in _ZHANG_TABLE]

    fig, ax = plt.subplots(figsize=(7, 6))
    lims = [0, max(xs) * 1.15]
    ax.plot(lims, [K_MEAN * v for v in lims], '--', color='#888888',
            label=f'k_mean fit (k={K_MEAN:.3f})')
    ax.scatter(xs, ys_actual, s=70, color=COLOR_POWER, zorder=3, label='Zhang et al. measured')
    # Two points (x=0.9, x=1.0) sit close together in both axes and their
    # labels overlap with a fixed offset; alternate the offset direction
    # for any point that is within a small fraction of the x-range of an
    # already-placed label so they fan out instead of colliding.
    x_span = max(xs) - min(xs) if max(xs) != min(xs) else 1.0
    placed = []
    for x, y, label in zip(xs, ys_actual, comp_labels):
        dx, dy = 6, 4
        for px, py in placed:
            if abs(px - x) / x_span < 0.08:
                dy = -14
                dx = 6
                break
        ax.annotate(label, xy=(x, y), xytext=(dx, dy), textcoords='offset points', fontsize=9)
        placed.append((x, y))
    ax.set_xlabel(r'$\Delta S \times T_{hys}$  [J/kg]')
    ax.set_ylabel(r'$W_{hys}$ measured [J/kg]')
    ax.set_title('Fit Quality: the k-Proxy Used to Extrapolate the Doped\n'
                 'V-on-Fe/V-on-Mn Estimates in fig40 (checked against real data)')
    ax.legend(fontsize=9)
    fig.tight_layout()
    save(fig, 'fig41_mnfepsi_hysteresis_kfit_quality')


# ══════════════════════════════════════════════════════════════════════════
# FIG 42 — Hysteresis-exploiting actuator work estimate
# ══════════════════════════════════════════════════════════════════════════

def plot_hysteresis_exploiting_actuator_estimate():
    """Figure for Phase 36's hysteresis_exploiting_actuator_estimate.py --
    the sigma*epsilon/rho specific-work grid (4 epsilon/rho combinations)
    plotted against this repo's own existing hysteresis_loss_J_per_kg
    range across five material families (2-65 J/kg), the module's own
    stated basis for comparison. Explicitly titled SPECULATIVE / FLOOR
    ESTIMATE, matching the module's own honesty flags (excludes actuator
    inefficiency; would only rise further)."""
    from core.hysteresis_exploiting_actuator_estimate import run_estimate
    result = run_estimate(verbose=False)
    estimates = result['estimates']
    labels = [f'eps={e.epsilon*100:.1f}%\nrho={e.rho_kg_m3:.0f}' for e in estimates]
    values = [e.specific_work_J_per_kg for e in estimates]

    fig, ax = plt.subplots(figsize=(8.5, 5.5))
    ax.axhspan(2.0, 65.0, color=COLOR_MAIN, alpha=0.12,
               label='This repo\'s existing hysteresis_loss_J_per_kg\nrange (5 material families)')
    bars = ax.bar(labels, values, color=COLOR_POWER, alpha=0.85, edgecolor='white', width=0.5)
    for b, v in zip(bars, values):
        ax.annotate(f'{v:.1f}', xy=(b.get_x() + b.get_width() / 2, v), xytext=(0, 4),
                    textcoords='offset points', ha='center', fontsize=9)
    ax.set_ylabel('Specific mechanical work [J/kg per stress half-cycle]')
    ax.set_title('SPECULATIVE FLOOR ESTIMATE: Actuator Work for the\n'
                 'Hysteresis-Exploiting Multicaloric Cycle (Phase 36)\n'
                 '(excludes actuator/hydraulic inefficiency -- would only rise further)')
    ax.legend(fontsize=8.5, loc='upper left')
    fig.tight_layout()
    save(fig, 'fig42_hysteresis_exploiting_actuator_estimate')


# ══════════════════════════════════════════════════════════════════════════
# FIG 43 — Combined hysteresis-loss landscape (Phase 36 findings together)
# ══════════════════════════════════════════════════════════════════════════

def plot_hysteresis_loss_landscape():
    """Ties fig40 and fig42 together on one axis: this repo's existing
    hysteresis_loss_J_per_kg range across five material families, next
    to Phase 36's two SPECULATIVE follow-ups -- the doped-material
    estimate (pushes the LOW end down, by roughly an order of magnitude)
    and the actuator-work estimate (a NEW, structurally different
    parasitic channel that lands INSIDE the existing range, not below
    or dramatically above it). Deliberately one summary figure rather
    than repeating fig40/fig42's own bars, so the two findings' relative
    position is visible at a glance."""
    from core.mnfepsi_doped_hysteresis_speculative import VONFE_ESTIMATE, VONMN_ESTIMATE
    from core.hysteresis_exploiting_actuator_estimate import (
        _EXISTING_HYSTERESIS_LOSS_RANGE_J_PER_KG, run_estimate)

    lo, hi = _EXISTING_HYSTERESIS_LOSS_RANGE_J_PER_KG
    actuator_vals = [e.specific_work_J_per_kg for e in run_estimate(verbose=False)['estimates']]
    doped_vals = [VONFE_ESTIMATE.estimate_J_per_kg_mean_k, VONFE_ESTIMATE.estimate_J_per_kg_ls_k,
                  VONMN_ESTIMATE.estimate_J_per_kg_mean_k, VONMN_ESTIMATE.estimate_J_per_kg_ls_k]

    fig, ax = plt.subplots(figsize=(9, 5.5))
    ax.axhspan(lo, hi, color=COLOR_MAIN, alpha=0.15,
               label=f'Existing hysteresis_loss_J_per_kg range\n({lo:.0f}-{hi:.0f} J/kg, 5 families)')
    ax.scatter([0.2] * len(doped_vals), doped_vals, s=70, color='#c9a227', zorder=3,
               label='SPECULATIVE: V-doped MnFePSi estimates\n(pushes the low end down)')
    ax.scatter([0.8] * len(actuator_vals), actuator_vals, s=70, color=COLOR_POWER, zorder=3,
               label='SPECULATIVE: hysteresis-exploiting\nactuator work (new channel, lands inside)')
    ax.set_yscale('log')
    ax.set_xlim(0, 1)
    ax.set_xticks([0.2, 0.8])
    ax.set_xticklabels(['Doped materials\n(Phase 36)', 'Dual-stimulus actuator\n(Phase 36)'])
    ax.set_ylabel('Specific loss/work [J/kg] (log scale)')
    ax.set_title('Hysteresis-Loss Landscape: Where Phase 36\'s Two\n'
                 'Speculative Follow-Ups Sit Relative to Validated Data')
    ax.legend(fontsize=8, loc='upper center', bbox_to_anchor=(0.5, -0.12), ncol=1)
    fig.tight_layout()
    save(fig, 'fig43_hysteresis_loss_landscape')


# ══════════════════════════════════════════════════════════════════════════
# FIG 44 — Phase 37 calibration drift: predicted vs. literature Qc
# ══════════════════════════════════════════════════════════════════════════

def plot_calibration_drift_phase37():
    """Documents the Phase 37 finding directly: for each of the 4
    recalibrated CORE/MAGGIE_HIGHSPAN points, shows literature Qc against
    (a) what the STALE (pre-Phase-37) mdot predicted and (b) what the
    RE-CALIBRATED mdot predicts (should equal literature Qc by
    construction, i.e. an exact match -- shown for direct visual
    contrast with how far off the stale value had drifted)."""
    from core.mce_material import GADOLINIUM
    T_COLD = 289.0
    # (name, freq, field, mass, span, old_mdot, new_mdot, Qc_lit, no_load_override)
    points = [
        ("Astronautics_rotary_2014", 4.0, 1.44, 1.52, 11.0, 0.252999, 0.309029, 2502.0, None),
        ("DTU_Eriksen_rotary_Gd_2015", 0.75, 1.13, 1.7, 10.2, 0.084666, 0.326616, 102.8, None),
        ("Tusek_singlebed_Gd_2010", 0.3, 1.15, 0.1763, 7.26, 0.007351, 0.021068, 5.27, None),
        ("DTU_Eriksen_MAGGIE_2016", 0.61, 1.13, 1.7, 15.5, 0.014650, 0.015606, 81.5, 21.04),
    ]
    names, old_qc, new_qc, lit_qc = [], [], [], []
    for name, f, H, mass, span, old_mdot, new_mdot, Qc_lit, override in points:
        kwargs = dict(material=GADOLINIUM, mu0H_max=H, mass_regenerator=mass, frequency=f)
        if override is not None:
            kwargs['no_load_span_override'] = override
        sys_old = AMRSystem(fluid_mdot=old_mdot, **kwargs)
        sys_new = AMRSystem(fluid_mdot=new_mdot, **kwargs)
        qc_old, _ = sys_old.cooling_capacity(T_COLD, span)
        qc_new, _ = sys_new.cooling_capacity(T_COLD, span)
        names.append(name.replace('_', '\n'))
        old_qc.append(qc_old)
        new_qc.append(qc_new)
        lit_qc.append(Qc_lit)

    x = np.arange(len(names))
    w = 0.27
    fig, ax = plt.subplots(figsize=(10, 5.5))
    ax.bar(x - w, lit_qc, w, label='Literature Qc (target)', color='#2ca02c', alpha=0.85,
           edgecolor='white')
    ax.bar(x, old_qc, w, label='Stale mdot (pre-Phase-37)', color=COLOR_POWER, alpha=0.85,
           edgecolor='white')
    ax.bar(x + w, new_qc, w, label='Recalibrated mdot (Phase 37)', color=COLOR_MAIN, alpha=0.85,
           edgecolor='white')
    for i, (o, l) in enumerate(zip(old_qc, lit_qc)):
        pct = 100 * (l - o) / l
        ax.annotate(f'{pct:.0f}% low', xy=(i, o), xytext=(0, 4), textcoords='offset points',
                    ha='center', fontsize=8, color=COLOR_POWER)
    ax.set_xticks(x)
    ax.set_xticklabels(names, fontsize=7.5)
    ax.set_ylabel('Qc [W]')
    ax.set_yscale('log')
    ax.set_title('Phase 37: CALIBRATION_POINTS_CORE Drift, Found and Fixed\n'
                 '(stale mdot under-predicted every point; recalibrated mdot matches by construction)')
    ax.legend(fontsize=9)
    fig.tight_layout()
    save(fig, 'fig44_calibration_drift_phase37')


# ══════════════════════════════════════════════════════════════════════════
# FIG 45 — Lozano rows: which ones calibrate, which structurally can't
# ══════════════════════════════════════════════════════════════════════════

def plot_lozano_calibration_status():
    """Documents the other Phase 37 finding: Lozano_POLO_UFSC_2016_r4's
    6.1K target span sits ABOVE this repo's model's own no-load span cap
    at 0.88T/0.4Hz, so Qc=0.0W at every mdot -- a structural, not a
    stale-calibration, non-match (excluded from
    test_rotary_drive_loss_model_substantially_improves_lozano_
    predictions as of Phase 37). Plots Qc vs. mdot for all four Lozano
    rows (r4, r6, r7, r8) on one log-log axis: r4's curve should visibly
    flatline at 0 while the other three rise smoothly and cross their
    own Qc_lit target (marked with a horizontal dashed line + star)."""
    from core.mce_material import GADOLINIUM
    T_COLD_ASSUMED_K = 294.0 - 5.0
    rows = [
        ("Lozano r4 (EXCLUDED, non-calibrating)", 0.4, 0.88, 6.1, 62.5),
        ("Lozano r6", 0.8, 0.88, 5.0, 81.2),
        ("Lozano r7", 0.4, 0.88, 3.7, 80.8),
        ("Lozano r8", 0.8, 0.88, 3.7, 120.4),
    ]
    mdots = np.logspace(-4, 0, 60)
    colors4 = ['#c00000', '#1f4e79', '#2ca02c', '#c9a227']

    fig, ax = plt.subplots(figsize=(8.5, 6))
    all_qcs = []
    for (name, f, H, span, Qc_lit) in rows:
        qcs = []
        for mdot in mdots:
            sys_ = AMRSystem(material=GADOLINIUM, mu0H_max=H, mass_regenerator=1.0,
                              frequency=f, fluid_mdot=mdot)
            qc, _ = sys_.cooling_capacity(T_COLD_ASSUMED_K, span)
            qcs.append(max(qc, 1e-4))
        all_qcs.append(qcs)

    # r7 and r8 land on numerically identical Qc(mdot) curves at this
    # operating point (frequency doesn't move the no-load-span-capped
    # result here), so plotting both as solid lines hides one completely
    # underneath the other. Give any curve that nearly coincides with an
    # already-drawn one a dashed style and a slightly thicker/offset
    # width so both remain visible.
    drawn = []
    for (name, f, H, span, Qc_lit), color, qcs in zip(rows, colors4, all_qcs):
        overlaps = any(
            np.allclose(np.log10(np.array(qcs)), np.log10(np.array(prev)), atol=1e-6)
            for prev in drawn
        )
        linestyle = '--' if overlaps else '-'
        linewidth = 3.5 if overlaps else 2.0
        ax.plot(mdots, qcs, color=color, linestyle=linestyle, linewidth=linewidth,
                label=f'{name} (span={span}K)' + (' [overlaps another curve]' if overlaps else ''))
        ax.axhline(Qc_lit, color=color, linestyle=':', alpha=0.6)
        drawn.append(qcs)

    ax.set_xscale('log')
    ax.set_yscale('log')
    ax.set_xlabel('mdot [kg/s]')
    ax.set_ylabel('Qc [W] (floored at 1e-4 for the log axis)')
    ax.set_title('Which Lozano Rows Calibrate? r4 Flatlines at Qc=0\n'
                 'for EVERY mdot (span exceeds this field/frequency\'s own span cap);\n'
                 'r6/r7/r8 rise smoothly and reach their dotted Qc_lit target')
    ax.legend(fontsize=8)
    fig.tight_layout()
    save(fig, 'fig45_lozano_calibration_status')


# ══════════════════════════════════════════════════════════════════════════
# FIG 46 — Alternative solid-state caloric technologies vs. VCC
# ══════════════════════════════════════════════════════════════════════════

def _load_or_build_vcc_cop_by_span():
    """Reads results/comparison_table.csv if it already exists (written by
    main.py's step 4 / this module's own plot_amr_vs_baselines() call
    earlier in a full pipeline run); otherwise builds it fresh via
    main.run_baseline_sweep() so this figure works standalone
    (`python plots.py`) too, exactly like every other figure here.
    Passes out_path=csv_path explicitly (rather than relying on
    run_baseline_sweep()'s own RESULTS_CSV default) so this keeps working
    when plots.RESULTS_DIR has been redirected -- e.g. under this test
    suite's own tmp_path fixture -- instead of silently reading from a
    location run_baseline_sweep() never actually wrote to."""
    csv_path = RESULTS_DIR / 'comparison_table.csv'
    if not csv_path.exists():
        import main as main_mod
        csv_path.parent.mkdir(parents=True, exist_ok=True)
        main_mod.run_baseline_sweep(out_path=str(csv_path))
    return alternative_caloric_comparison.load_vcc_cop_by_span_from_csv(str(csv_path))


def plot_alternative_caloric_comparison():
    """Figure for alternative_caloric_comparison.py: does ANY solid-state
    caloric cooling technology -- not just this repo's own magnetocaloric
    model -- beat this repo's own VCC electrical COP in the data-center
    5-20K span range? Left panel: literature COP claims (elastocaloric,
    barocaloric, electrocaloric) plotted against this repo's own VCC COP
    at each claim's matched span, MEASURED-device claims shown solid,
    SIMULATION/PROJECTION/vendor claims hatched (matching this repo's own
    material_family_comparison.py hatching convention for "outside the
    documented window"). Right panel: the more rigorous physics-model span
    sweep (this repo's own elastocaloric_cycle.py/barocaloric_cycle.py/
    electrocaloric_cycle.py models, calibrated to those same literature
    points, run across the SAME 5-20K grid as comparison_table.csv) vs.
    VCC, so the headline "does this generalize past MCE" finding is
    visible as a chart, not only as a text report."""
    vcc = _load_or_build_vcc_cop_by_span()
    claim_results = alternative_caloric_comparison.compare_all_at_matched_spans(
        vcc, verbose=False)
    sweep_results, _, _ = alternative_caloric_comparison.compare_physics_models_across_spans(
        vcc, verbose=False)

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(15, 6.5))

    # --- Left: literature claims vs. VCC at each claim's matched span ---
    # Wrap (not mid-word truncate) the system label so short names like
    # "Elastocaloric can-cooler (Ehl et al. 2025)" render in full across
    # up to two lines instead of being cut off mid-word.
    labels = [f"{r['technology']}\n" + textwrap.fill(r['system'], width=22, max_lines=2, placeholder='...')
              for r in claim_results]
    claim_cops = [r['claimed_cop'] for r in claim_results]
    vcc_cops_matched = [r['vcc_cop_at_matched_span'] for r in claim_results]
    measured = [r['is_measured_device_cop'] for r in claim_results]
    x = np.arange(len(labels))
    w = 0.35
    bars_claim = ax1.bar(x - w / 2, claim_cops, w, color=COLOR_POWER, alpha=0.85,
                          edgecolor='white', label='Literature caloric-technology COP')
    for b, m in zip(bars_claim, measured):
        if not m:
            b.set_hatch('//')
    ax1.bar(x + w / 2, vcc_cops_matched, w, color=COLOR_MAIN, alpha=0.85,
            edgecolor='white', label="This repo's own VCC COP (matched span)")
    for xi, c in zip(x, claim_cops):
        ax1.annotate(f'{c:.1f}', xy=(xi - w / 2, c), xytext=(0, 4),
                     textcoords='offset points', ha='center', fontsize=8)
    ax1.set_xticks(x)
    ax1.set_xticklabels(labels, fontsize=7, rotation=20, ha='right')
    ax1.set_ylabel('COP')
    ax1.set_title('Literature Claims vs. This Repo\'s VCC\n'
                   + textwrap.fill('(hatched = simulation/projection/vendor claim, not a measured device)',
                                    width=42),
                   fontsize=10)
    ax1.legend(fontsize=8)

    # --- Right: physics-model span sweep ---
    spans = [r['span_K'] for r in sweep_results]
    vcc_line = [r['vcc_cop'] for r in sweep_results]
    eq_line = [r['elastocaloric_cop_electrical'] for r in sweep_results]
    bc_line = [r['barocaloric_cop_electrical'] for r in sweep_results]
    ec_line = [r['electrocaloric_cop_electrical'] for r in sweep_results]
    ec_far = [r['electrocaloric_is_far_extrapolation'] for r in sweep_results]
    ec_extrap = [r['electrocaloric_is_extrapolated'] for r in sweep_results]

    ax2.plot(spans, vcc_line, color='k', marker='o', label="This repo's VCC COP")
    ax2.plot(spans, eq_line, color='#2ca02c', marker='s', label='Elastocaloric (calibrated model)')
    ax2.plot(spans, bc_line, color='#c9a227', marker='^', label='Barocaloric (calibrated model)')
    # Electrocaloric: solid where it's the one measured calibration point,
    # dashed everywhere else (extrapolated), so the single-point basis for
    # this line is visible rather than implied to be validated everywhere.
    for i in range(len(spans) - 1):
        style = '-' if not (ec_extrap[i] or ec_extrap[i + 1]) else '--'
        color = '#c00000' if not (ec_far[i] or ec_far[i + 1]) else '#e88'
        ax2.plot(spans[i:i + 2], ec_line[i:i + 2], color=color, linestyle=style, linewidth=2)
    ax2.plot([], [], color='#c00000', marker='d', label='Electrocaloric (1 measured pt, extrapolated elsewhere)')
    ax2.scatter([s for s, e in zip(spans, ec_extrap) if not e],
                [c for c, e in zip(ec_line, ec_extrap) if not e],
                color='#c00000', marker='d', s=60, zorder=5)
    ax2.scatter([s for s, e in zip(spans, ec_extrap) if e],
                [c for c, e in zip(ec_line, ec_extrap) if e],
                color='#e88', marker='d', s=30, zorder=4)
    ax2.set_xlabel('Temperature Span [K]')
    ax2.set_ylabel('COP')
    ax2.set_title('Physics-Model Span Sweep\n'
                   + textwrap.fill("(this repo's own literature-calibrated models, same grid as "
                                    "comparison_table.csv)", width=42),
                   fontsize=10)
    ax2.legend(fontsize=7.5)

    fig.suptitle('Does ANY Solid-State Caloric Technology Beat This Repo\'s Own VCC COP?\n'
                 'Honest answer: only one measured device (electrocaloric, 20.9K, 2.1W) — everything else is simulation, projection, or loses',
                 fontsize=11)
    fig.tight_layout(rect=[0, 0, 1, 0.94])
    fig.subplots_adjust(wspace=0.3, bottom=0.28)
    save(fig, 'fig46_alternative_caloric_comparison')


# ══════════════════════════════════════════════════════════════════════════
# FIG 47 — Annual water usage (WUE) comparison
# ══════════════════════════════════════════════════════════════════════════

def plot_water_usage_comparison(precomputed=None):
    """Figure for water_usage.py: converts each technology's own already-
    computed electrical COP into an annual water-consumption figure via
    the industry-standard Water Usage Effectiveness (WUE) metric, the
    same three-technology shape as plot_emissions() (fig23). Uses the
    representative 10K-span row from comparison_table.csv (same
    convention as plot_economics()/plot_emissions()) rather than
    water_usage.py's own illustrative __main__ defaults, unless that row
    isn't available yet, in which case those defaults are used and the
    figure says so explicitly."""
    precomputed = precomputed or {}
    rows = precomputed.get('baseline_rows')
    is_default = False
    if rows is None:
        csv_path = RESULTS_DIR / 'comparison_table.csv'
        if csv_path.exists():
            rows = _read_csv_rows(str(csv_path))
    if rows:
        rep = min(rows, key=lambda r: abs(r['span_K'] - 10.0))
        amr_cop = rep['AMR_COP_electrical']
        vcc_cop = rep['VaporCompression_COP']
        liquid_cop = rep['LiquidCooling_COP']
    else:
        is_default = True
        amr_cop, vcc_cop, liquid_cop = 4.63, 3.2, 4.0

    results = water_usage.compare_water_usage(100.0, amr_cop, vcc_cop, liquid_cop)
    labels = [r.technology for r in results]
    liters_per_kw = [r.annual_water_liters_per_kW_IT for r in results]
    wues = [r.WUE_L_per_kWh_IT for r in results]
    colors = [COLOR_MAIN, COLOR_POWER, '#2ca02c']

    fig, ax = plt.subplots(figsize=(9, 6))
    bars = ax.bar(labels, liters_per_kw, color=colors, alpha=0.85, edgecolor='white')
    # Give the tallest bar's annotation room to breathe above the bar so it
    # can't collide with the (multi-line) title sitting above the axes.
    ax.set_ylim(0, max(liters_per_kw) * 1.18)
    for b, v, w_ in zip(bars, liters_per_kw, wues):
        ax.annotate(f'{v:,.0f} L/kW-IT/yr\n(WUE={w_:.2f} L/kWh)', xy=(b.get_x() + b.get_width() / 2, v),
                    xytext=(0, 4), textcoords='offset points', ha='center', fontsize=9)
    ax.set_ylabel('Annual water consumption [L per kW-IT per year]')
    title = ('Annual Water Usage (WUE) Comparison, 100kW-IT Facility\n'
             + textwrap.fill('(AMR assigned dry/air-cooled rejection by default — a design '
                              'choice, not a property of the magnetocaloric cycle itself)', width=70))
    if is_default:
        title += '\n[ILLUSTRATIVE — comparison_table.csv not found, using placeholder COPs]'
    ax.set_title(title, fontsize=10.5)
    fig.tight_layout()
    save(fig, 'fig47_water_usage_comparison')


# ══════════════════════════════════════════════════════════════════════════
# FIG 48 — Monte Carlo calibration-uncertainty band on COP_electrical
# ══════════════════════════════════════════════════════════════════════════

def plot_uncertainty_propagation():
    """Figure for uncertainty_propagation.py: overlays a 90% confidence
    band (Monte Carlo resampling of the CORE loss-model's 3 calibration
    points under an assumed +/-15% measurement-noise model) on top of
    comparison_table.csv's own point-value AMR_COP_electrical curve,
    across the full 5-20K ASHRAE span sweep -- the "single point-value
    predictions with no uncertainty band" gap this module's own docstring
    names explicitly. n_draws is reduced from the module's own default
    (2000 at a single span / 500 across the sweep) to keep this figure's
    own runtime bounded; the band shape is stable well below that."""
    rows = uncertainty_propagation.uncertainty_band_across_spans(
        spans_K=range(5, 21), n_draws=120, verbose=False)
    spans = [r['span_K'] for r in rows]
    mean = [r['COP_electrical_mean'] for r in rows]
    p05 = [r['COP_electrical_p05'] for r in rows]
    p95 = [r['COP_electrical_p95'] for r in rows]

    # NaN rows (every draw infeasible at that span, e.g. above the 0-D
    # model's own no-load span cap) would otherwise silently break the
    # fill_between/plot calls -- mask them out and note it explicitly
    # rather than plotting a misleading gap-free line through them.
    valid = [i for i in range(len(spans)) if np.isfinite(mean[i])]
    n_masked = len(spans) - len(valid)

    fig, ax = plt.subplots(figsize=(9, 6))
    vs = [spans[i] for i in valid]
    vm = [mean[i] for i in valid]
    vlo = [p05[i] for i in valid]
    vhi = [p95[i] for i in valid]
    ax.fill_between(vs, vlo, vhi, color=COLOR_MAIN, alpha=0.25,
                     label='90% CI (assumed +/-15% calibration-input noise)')
    ax.plot(vs, vm, color=COLOR_MAIN, marker='o', label='Monte Carlo mean COP_electrical')
    if n_masked:
        ax.annotate(f'{n_masked} span(s) omitted: all MC draws infeasible\n'
                    '(0-D model\'s own no-load span cap)',
                    xy=(0.02, 0.02), xycoords='axes fraction', fontsize=8,
                    color='#888888', style='italic')
    ax.set_xlabel('Temperature Span [K]')
    ax.set_ylabel('Electrical COP')
    ax.set_title('Calibration-Uncertainty Band on AMR Electrical COP\n'
                 '(Monte Carlo over the CORE loss-model\'s 3 calibration points, '
                 'assumed +/-15% measurement noise — not a source-derived error bar)')
    ax.legend(fontsize=9)
    fig.tight_layout()
    save(fig, 'fig48_uncertainty_propagation_band')


# ══════════════════════════════════════════════════════════════════════════
# FIG 49 — NSGA-III Pareto-front multiseed stability
# ══════════════════════════════════════════════════════════════════════════

def plot_pareto_multiseed_stability():
    """Figure for pareto_multiseed_stability.py: reruns the main material+
    geometry NSGA-III co-optimization (optimize.run_optimization(), the
    engine behind fig18/pareto_front.csv) across several independent
    seeds, at a REDUCED pop_size/n_gen from the module's own production
    default (40/25) to keep this figure's own runtime bounded -- the
    question being tested (does the search's own random seed move the
    headline numbers) does not require production-scale settings to
    demonstrate qualitatively. Left panel: material-family share of the
    merged Pareto front, mean +/- std across seeds (the "100%
    La(Fe,Si)13Hy"-style claim this repo has made before). Right panel:
    knee-point COP_electrical and best COP_electrical, mean +/- std."""
    from core import pareto_multiseed_stability
    result = pareto_multiseed_stability.run_pareto_multiseed_stability_check(
        seeds=(1, 2, 3, 4, 5), pop_size=20, n_gen=10, verbose=False)
    summary = result['summary']
    if summary is None:
        fig, ax = plt.subplots(figsize=(8, 5))
        ax.text(0.5, 0.5, 'All seeds returned an empty Pareto front --\nnothing to plot.',
                ha='center', va='center', fontsize=11)
        ax.axis('off')
        save(fig, 'fig49_pareto_multiseed_stability')
        return

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(13, 6))

    mats = sorted(summary['material_share_stats'].keys(),
                  key=lambda m: -summary['material_share_stats'][m]['mean_pct'])
    means = [summary['material_share_stats'][m]['mean_pct'] for m in mats]
    stds = [summary['material_share_stats'][m]['std_pct'] for m in mats]
    ax1.bar(range(len(mats)), means, yerr=stds, capsize=4, color=COLOR_MAIN,
            alpha=0.85, edgecolor='white')
    ax1.set_xticks(range(len(mats)))
    ax1.set_xticklabels(mats, rotation=30, ha='right', fontsize=7)
    ax1.set_ylabel('Share of merged Pareto front [%]')
    ax1.set_title('Material-Family Share, Mean ± Std\nAcross 5 NSGA-III Seeds')

    labels2 = ['Best COP\n(across front)', 'Knee-point COP\n(balanced design)']
    means2 = [summary['best_COP_electrical_mean'], summary['knee_COP_electrical_mean']]
    stds2 = [summary['best_COP_electrical_std'], summary['knee_COP_electrical_std']]
    ax2.bar(labels2, means2, yerr=stds2, capsize=6, color=COLOR_POWER, alpha=0.85,
            edgecolor='white', width=0.5)
    for i, (m, s) in enumerate(zip(means2, stds2)):
        # Anchor at the top of the error bar (m + s), not the bar top (m),
        # so the label sits above the error-bar cap instead of being
        # bisected by it.
        ax2.annotate(f'{m:.2f} \u00b1 {s:.2f}', xy=(i, m + s), xytext=(0, 8),
                    textcoords='offset points', ha='center', fontsize=9)
    ax2.set_ylim(0, max(m + s for m, s in zip(means2, stds2)) * 1.15)
    ax2.set_ylabel('Electrical COP')
    consistent = summary['knee_material_consistent_across_seeds']
    ax2.set_title(f'Headline COP Numbers, Mean ± Std\n'
                  f'Knee-point material {"consistent" if consistent else "NOT consistent"} '
                  f'across seeds: {summary["knee_materials_seen"]}')

    fig.suptitle('NSGA-III Pareto Front: Seed-to-Seed Stability Check\n'
                 '(reduced pop_size=20/n_gen=10 for this figure — see docstring)',
                 fontsize=12)
    fig.tight_layout()
    save(fig, 'fig49_pareto_multiseed_stability')


# ══════════════════════════════════════════════════════════════════════════
# FIG 50 — Regime-crossover search: does AMR ever beat VCC anywhere?
# ══════════════════════════════════════════════════════════════════════════

def plot_regime_crossover_analysis():
    """Figure for regime_crossover_analysis.py: the direct visual answer
    to "where, if anywhere, does this repo's own model show magnetic
    cooling beating conventional cooling" -- a systematic search, not one
    favorable-looking point. Left panel: best achievable AMR COP found
    across a broad design grid at each span, against VCC's COP at every
    eta_2nd_law tried (including VCC's own worst-case setting) -- the
    band between VCC's best and worst case vs. the single best-AMR line
    makes the gap's size directly visible. Right panel: the more generous
    best-material/best-field graded-cascade search (fields up to 7T),
    annotated with which field the search actually picked at each span
    (the "search never wants the high fields" finding)."""
    cop_search = regime_crossover_analysis.run_cop_crossover_search(verbose=False)
    mat_search = regime_crossover_analysis.run_material_and_field_crossover_search(verbose=False)

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(13, 6))

    spans1 = [r['span_K'] for r in cop_search['rows']]
    best_amr = [r['best_AMR_COP_electrical'] for r in cop_search['rows']]
    vcc_best = [max(r['VCC_COPs_by_eta'].values()) for r in cop_search['rows']]
    vcc_worst = [min(r['VCC_COPs_by_eta'].values()) for r in cop_search['rows']]
    ax1.fill_between(spans1, vcc_worst, vcc_best, color=COLOR_MAIN, alpha=0.2,
                     label='VCC COP range (eta=0.25-0.55)')
    ax1.plot(spans1, vcc_worst, color=COLOR_MAIN, linestyle=':', linewidth=1.2,
             label="VCC's own WORST-case COP (eta=0.25)")
    ax1.plot(spans1, best_amr, color=COLOR_POWER, marker='o',
             label='Best AMR COP found\n(broad grid: mass/freq/mdot/field, fixed Gd)')
    ax1.set_xlabel('Temperature Span [K]')
    ax1.set_ylabel('Electrical COP')
    ax1.set_title('COP Crossover Search\n(NO crossover found — AMR never beats even VCC\'s own worst case)')
    ax1.legend(fontsize=8)

    spans2 = [r['span_K'] for r in mat_search['rows']]
    best_amr2 = [r['best_AMR_COP_electrical'] for r in mat_search['rows']]
    vcc2 = [r['VCC_COP'] for r in mat_search['rows']]
    fields2 = [r['field_at_best_design_T'] for r in mat_search['rows']]
    ax2.plot(spans2, vcc2, color=COLOR_MAIN, marker='s', label='VCC COP (eta=0.42)')
    ax2.plot(spans2, best_amr2, color=COLOR_POWER, marker='o',
             label='Best AMR COP found\n(best material, graded cascade, field up to 7T)')
    for s, c, f in zip(spans2, best_amr2, fields2):
        ax2.annotate(f'{f:.0f}T', xy=(s, c), xytext=(0, -14), textcoords='offset points',
                    ha='center', fontsize=8, color='#555555')
    ax2.set_xlabel('Temperature Span [K]')
    ax2.set_ylabel('Electrical COP')
    ax2.set_title('Best-Material/Best-Field Search\n'
                  '(labels = field the search picked — never the 5-7T ceiling offered)')
    ax2.legend(fontsize=8)

    fig.suptitle('Does This Repo\'s Own Model Show AMR Beating VCC Anywhere?\n'
                 'Two independent searches, both null — see regime_crossover_analysis.py',
                 fontsize=12)
    fig.tight_layout()
    save(fig, 'fig50_regime_crossover_analysis')


# ══════════════════════════════════════════════════════════════════════════
# FIG 51 — Ames Lab heat-pump architecture check
# ══════════════════════════════════════════════════════════════════════════

def plot_heat_pump_validation():
    """Figure for heat_pump_validation.py: this repo's own default Gd
    packed-bed AMR design (the SAME architecture Ames Lab's real device
    uses) at a representative residential-heat-pump operating point,
    shown ALONGSIDE Ames Lab's real reported whole-device specific power
    density (SPD) figures -- deliberately on two SEPARATE panels with no
    computed ratio between them, since the two quantities have different
    mass denominators (this repo's own MCM-only mass vs. Ames Lab's
    whole-device mass including magnet/motor/housing) and are therefore
    not directly comparable -- see the module's own honesty flag."""
    result = heat_pump_validation.run_ames_lab_architecture_check()

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(11, 5.5))

    if result['model_specific_cooling_power_w_per_kg_MCM'] is not None:
        ax1.bar(['This repo\'s model\n(Gd, packed-bed AMR)'],
                [result['model_specific_cooling_power_w_per_kg_MCM']],
                color=COLOR_MAIN, alpha=0.85, edgecolor='white', width=0.5)
        ax1.annotate(f"{result['model_specific_cooling_power_w_per_kg_MCM']:.1f} W/kg-MCM",
                    xy=(0, result['model_specific_cooling_power_w_per_kg_MCM']),
                    xytext=(0, 4), textcoords='offset points', ha='center', fontsize=10)
    ax1.set_ylabel('Specific cooling power [W per kg of MCM ONLY]')
    ax1.set_title(f"This Repo's Model\n(T_cold={result['T_cold_K']-273.15:.0f}°C, "
                  f"span={result['span_K']:.0f}K, {result['AMR_n_stages']} stage(s))")

    labels = ['Baseline', 'Optimized', 'Projected\nceiling']
    vals = [result['ames_baseline_SPD_w_per_kg_whole_device'],
            result['ames_optimized_SPD_w_per_kg_whole_device'],
            result['ames_projected_max_SPD_w_per_kg_whole_device']]
    ax2.bar(labels, vals, color=COLOR_POWER, alpha=0.85, edgecolor='white')
    for i, v in enumerate(vals):
        ax2.annotate(f'{v:.1f}', xy=(i, v), xytext=(0, 4), textcoords='offset points',
                    ha='center', fontsize=9)
    ax2.set_ylabel('Specific power density [W per kg, WHOLE DEVICE]')
    ax2.set_title('Ames Lab\'s Real Device (Applied Energy 377, 2025)\n'
                  '(magnet + motor + housing + everything)')

    fig.suptitle('Ames Lab Heat-Pump Architecture Check — Same Architecture, Different Metric\n'
                 'NOT directly comparable: MCM-only mass (left) vs. whole-device mass (right) — see docstring',
                 fontsize=11)
    fig.tight_layout()
    save(fig, 'fig51_heat_pump_validation')


# ══════════════════════════════════════════════════════════════════════════
# FIG 52 — Hypereg parallel-hydraulic pumping-power reduction
# ══════════════════════════════════════════════════════════════════════════

def plot_hypereg_analysis():
    """Figure for hypereg_analysis.py: does Klinar et al.'s (2024)
    parallel-hydraulic "Hypereg" regenerator split give a genuine,
    non-negligible COP benefit in this repo's own calibrated loss model?
    Left panel: COP_electrical vs. n_parallel at the module's own
    representative operating point, showing the benefit saturate quickly
    (eddy-current and base-overhead losses are untouched by
    parallelization). Right panel: the ROADMAP-motivated follow-up --
    does the benefit become non-negligible at a 4x higher mdot, where
    pumping power is a bigger share of total loss?"""
    n_rows = hypereg_analysis.sweep_n_parallel(verbose=False)
    n_rows_high_mdot = hypereg_analysis.sweep_n_parallel_at_higher_mdot(verbose=False)

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 5.5))

    ns = [r[0] for r in n_rows]
    cops = [r[2] for r in n_rows]
    ax1.plot(ns, cops, color=COLOR_MAIN, marker='o')
    ax1.set_xscale('log', base=2)
    ax1.set_xticks(ns)
    ax1.set_xticklabels([str(n) for n in ns])
    ax1.set_xlabel('n_parallel sub-regenerators (n=1 = conventional series flow)')
    ax1.set_ylabel('Electrical COP')
    ax1.set_title(f'COP vs. Hypereg Split\n(mdot={hypereg_analysis.MDOT_KG_S}kg/s baseline — benefit saturates)')

    conv_cop = n_rows[0][2]
    best_cop = max(n_rows, key=lambda r: r[2])[2]
    conv_cop_high = n_rows_high_mdot[0][2]
    best_cop_high = max(n_rows_high_mdot, key=lambda r: r[2])[2]
    rel_gain_baseline = 100 * (best_cop - conv_cop) / conv_cop if conv_cop > 0 else 0.0
    rel_gain_high = 100 * (best_cop_high - conv_cop_high) / conv_cop_high if conv_cop_high > 0 else 0.0

    labels = [f'Baseline mdot\n({hypereg_analysis.MDOT_KG_S}kg/s)',
              f'Higher mdot\n({0.3}kg/s)']
    gains = [rel_gain_baseline, rel_gain_high]
    ax2.bar(labels, gains, color=COLOR_POWER, alpha=0.85, edgecolor='white', width=0.5)
    for i, g in enumerate(gains):
        ax2.annotate(f'{g:.2f}%', xy=(i, g), xytext=(0, 4), textcoords='offset points',
                    ha='center', fontsize=10)
    ax2.set_ylabel('Best relative COP gain from Hypereg split [%]')
    ax2.set_title('Does a Higher Flow Rate Make the\nBenefit Non-Negligible? (ROADMAP follow-up)')

    fig.suptitle('Hypereg Parallel-Hydraulic Regenerator Split (Klinar et al. 2024)\n'
                 'Real but modest benefit — pumping power is only one of three loss channels',
                 fontsize=12)
    fig.tight_layout()
    save(fig, 'fig52_hypereg_analysis')


# ══════════════════════════════════════════════════════════════════════════
# FIG 53 — Regenerative-amplification gap across the full benchmark corpus
# ══════════════════════════════════════════════════════════════════════════

def plot_regenerative_amplification_gap():
    """Figure for validation_system.analyze_regenerative_amplification_gap():
    the systematic, corpus-wide version of the point-check fig35 already
    shows for 3 devices. Every span>0 benchmark row is evaluated against
    cooling_capacity()'s own structural ceiling (2*dTad_noload, the most
    span the 0-D single-blow model can EVER reach at that field/T_mid, for
    any mdot) and the ratio actual_span/structural_cap is plotted --
    ratio>1 means the real device demonstrably exceeds what the 0-D model
    can represent, a lower bound on the "regenerative amplification"
    (temperature-profile build-up over many cycles) effect a single-blow
    dTad model cannot capture by construction. This is new information no
    existing figure shows in aggregate: fig35 checks 3 specific spans
    against an override fix; this shows the shape and prevalence of the
    underlying gap across all ~19 usable benchmark rows at once, which is
    what actually motivates fig35's fix and regenerator_1d.py's 1-D
    alternative in the first place."""
    entries = validation_system.analyze_regenerative_amplification_gap(verbose=False)
    clean = [e for e in entries if not e['near_zero']]
    clean.sort(key=lambda e: e['amplification_ratio'])
    names = [e['device'].replace('_', ' ') for e in clean]
    ratios = [e['amplification_ratio'] for e in clean]
    exceeds = [r > 1.0 for r in ratios]
    colors = [COLOR_POWER if e else COLOR_MAIN for e in exceeds]

    fig, ax = plt.subplots(figsize=(10, 8))
    y = np.arange(len(names))
    ax.barh(y, ratios, color=colors, alpha=0.85, edgecolor='white')
    ax.axvline(1.0, color='black', linestyle='--', linewidth=1.2,
               label='Model\'s own structural cap (ratio = 1)')
    ax.set_yticks(y)
    ax.set_yticklabels(names, fontsize=8)
    ax.set_xscale('log')
    ax.set_xlabel('Amplification ratio = span_K / (2\u00d7dTad_noload)  [log]')
    n_exceed = sum(exceeds)
    ax.set_title('Regenerative-Amplification Gap Across the Full Benchmark Corpus\n'
                 f'{n_exceed}/{len(clean)} real devices exceed the 0-D model\'s own structural span cap',
                 fontsize=12)
    from matplotlib.patches import Patch
    ax.legend(handles=[Patch(color=COLOR_POWER, label='Exceeds structural cap (ratio > 1)'),
                        Patch(color=COLOR_MAIN, label='Within structural cap (ratio \u2264 1)'),
                        plt.Line2D([0], [0], color='black', linestyle='--', label='ratio = 1')],
              fontsize=8, loc='lower right')
    fig.tight_layout()
    save(fig, 'fig53_regenerative_amplification_gap_corpus')


# ══════════════════════════════════════════════════════════════════════════
# FIG 54 — Hybrid solid-state regenerator (HMR) vs. VCC, ideal -> real
# ══════════════════════════════════════════════════════════════════════════

def plot_hmr_vs_vcc_loss_funnel(T_cold_K=291.15, span_K=10.0):
    """Figure for hybrid_solid_state_regenerator.compare_to_vcc_realistic():
    a genuinely different alternative-cooling architecture (all-solid-state
    magnetic regenerator, Lin et al. 2024) from the fluid/elastocaloric/
    barocaloric/electrocaloric technologies fig46 already covers, at this
    repo's own ASHRAE-representative operating point. Shows, per tested
    frequency, the same "ideal claim collapses once real losses are added"
    funnel this repo's honesty framing already applies elsewhere (fig04,
    fig08, fig46): COP_ideal (frictionless, magnetic-work only, the
    source paper's own headline number) -> COP_ideal_airgap (adds the
    paper's own worst-case contact-friction derating) -> COP_HMR_electrical
    (this repo's new TIER-3 addition: rotary drivetrain + baseline
    overhead), against VCC's own real installed electrical COP at the same
    point. New information: no existing figure compares this architecture
    to anything."""
    result = hybrid_solid_state_regenerator.compare_to_vcc_realistic(
        T_cold_K, span_K, log=lambda *_: None)
    rows = result['rows']
    freqs = [r['frequency_Hz'] for r in rows]
    ideal = [r['COP_ideal_frictionless'] for r in rows]
    airgap = [r['COP_ideal_airgap'] for r in rows]
    real = [r['COP_HMR_electrical'] for r in rows]
    vcc_cop = result['vcc'].COP

    fig, ax = plt.subplots(figsize=(9.5, 6.5))
    x = np.arange(len(freqs))
    w = 0.25
    ax.bar(x - w, ideal, w, color='#a6c8e0', edgecolor='white',
           label='COP$_{ideal}$ (frictionless, magnetic-work only)')
    ax.bar(x, airgap, w, color='#5b9bd5', edgecolor='white',
           label='COP$_{ideal,airgap}$ (+ contact friction)')
    ax.bar(x + w, real, w, color=COLOR_POWER, edgecolor='white',
           label='COP$_{HMR,electrical}$ (+ drivetrain + baseline overhead)')
    ax.axhline(vcc_cop, color='black', linestyle='--', linewidth=1.5,
               label=f"VCC real electrical COP ({vcc_cop:.2f})")
    ax.set_xticks(x)
    ax.set_xticklabels([f'{f:g} Hz' for f in freqs])
    ax.set_ylabel('COP')
    ax.set_xlabel('Cycle frequency')
    best = max(rows, key=lambda r: r['COP_HMR_electrical'])
    verdict = 'still beats' if best['HMRe_over_VCC'] >= 1.0 else 'does NOT beat'
    ax.set_title('Hybrid Solid-State Regenerator (HMR) vs. VCC: Ideal \u2192 Real\n'
                 f"(Lin et al. 2024, Gd/Cu; best-frequency HMR real electrical COP "
                 f"{verdict} VCC's, {best['HMRe_over_VCC']:.2f}\u00d7)",
                 fontsize=11)
    ax.legend(fontsize=8.5)
    fig.tight_layout()
    save(fig, 'fig54_hmr_vs_vcc_loss_funnel')


# ══════════════════════════════════════════════════════════════════════════
# FIG 55 — Annualized, climate-weighted COP comparison
# ══════════════════════════════════════════════════════════════════════════

def plot_annualized_climate_comparison():
    """Figure for pue_annualized.annualized_energy_comparison(): every
    other COP comparison in this repo (comparison_table.csv and everything
    downstream of it, including fig08/fig19/fig36) is evaluated at ONE
    fixed (T_cold, span) design point. This figure tests whether that
    single-point conclusion survives across a full representative annual
    climate profile (ASHRAE zone 4A-like, 6 bins), which is new
    information -- it is the only figure in this repo that varies the
    outdoor condition rather than holding it fixed. Left panel: COP by
    technology at each climate bin (AMR bars absent where the implied
    span exceeds this repo's own validated 5-20K AMR envelope, flagged
    rather than silently zeroed). Right panel: the bin-weighted annual
    effective COP per technology, with AMR's annotated hours-coverage
    caveat and liquid cooling's economizer-credit caveat carried over
    explicitly from the module's own honest-framing text."""
    result = pue_annualized.annualized_energy_comparison(verbose=False)
    rows = result['rows']
    bins = [r['bin'] for r in rows]
    amr = [r['AMR_COP'] if r['AMR_span_feasible'] else np.nan for r in rows]
    vcc = [r['VCC_COP'] for r in rows]
    liq = [r['Liquid_COP'] for r in rows]

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 6.5),
                                    gridspec_kw={'width_ratios': [1.6, 1]})
    x = np.arange(len(bins))
    w = 0.25
    ax1.bar(x - w, amr, w, color=COLOR_MAIN, edgecolor='white', label='AMR')
    ax1.bar(x, vcc, w, color=COLOR_POWER, edgecolor='white', label='Vapor-compression')
    ax1.bar(x + w, liq, w, color='#2ca02c', edgecolor='white', label='Liquid cooling')
    for xi, a in zip(x, amr):
        if np.isnan(a):
            ax1.annotate('n/a\n(span\n>20K)', xy=(xi - w, 0), xytext=(0, 3),
                        textcoords='offset points', ha='center', fontsize=6.5, color='#888')
    ax1.set_xticks(x)
    ax1.set_xticklabels([textwrap.fill(b, width=10) for b in bins], fontsize=8)
    ax1.set_ylabel('COP at that climate bin')
    ax1.set_title('COP by Climate Bin\n(AMR: n/a where implied span exceeds its validated 5-20K envelope)',
                  fontsize=10)
    ax1.legend(fontsize=8)

    labels2 = ['AMR', 'Vapor-\ncompression', 'Liquid\ncooling']
    vals2 = [result['AMR_effective_annual_COP'] or 0, result['VCC_effective_annual_COP'],
             result['Liquid_effective_annual_COP']]
    colors2 = [COLOR_MAIN, COLOR_POWER, '#2ca02c']
    bars = ax2.bar(labels2, vals2, color=colors2, edgecolor='white', alpha=0.9)
    for b, v in zip(bars, vals2):
        ax2.annotate(f'{v:.2f}', xy=(b.get_x() + b.get_width() / 2, v), xytext=(0, 4),
                    textcoords='offset points', ha='center', fontsize=9)
    ax2.set_ylim(0, max(vals2) * 1.25)
    ax2.set_ylabel('Bin-weighted annual effective COP')
    ax2.set_title('Annual Effective COP', fontsize=11)
    cov = result['AMR_annual_hours_fraction_covered'] * 100
    ax2.annotate(textwrap.fill(
        f'AMR computed over only {cov:.0f}% of annual hours; liquid cooling includes an '
        'economizer credit AMR/VCC do not receive', width=34),
        xy=(0.5, 0.98), xycoords='axes fraction', ha='center', va='top',
        fontsize=7.5, color='#666666', style='italic')

    fig.suptitle('Does the Single-Point COP Conclusion Survive a Full Climate Year?', fontsize=12)
    fig.tight_layout(rect=[0, 0, 1, 0.92])
    save(fig, 'fig55_annualized_climate_comparison')


# ══════════════════════════════════════════════════════════════════════════
# FIG 56 — Heat-transfer fluid selection comparison
# ══════════════════════════════════════════════════════════════════════════

def plot_fluid_selection_comparison(precomputed=None):
    """Figure for fluid_selection_optimization.run_fluid_selection_comparison():
    every other figure in this repo silently assumes one fixed working
    fluid. This is the only figure showing the actual trade-off behind
    that choice -- COP_electrical (with mdot independently re-optimized
    per fluid) at two different operating points, ranked, with each
    fluid's %-COP giveup relative to the pure-water ceiling (not a real
    hardware option, since Gd corrodes in plain water). New information:
    justifies core/fluids.py's DEFAULT_FLUID against the alternatives
    quantitatively rather than by assertion, and shows whether that
    ranking is robust to a second, higher-flow/lower-mass operating
    point."""
    precomputed = precomputed or {}
    result = precomputed.get('fluid_selection_result')
    if result is None:
        result = fluid_selection_optimization.run_fluid_selection_comparison(verbose=False)
    baseline = {r['fluid']: r for r in result['baseline_rows']}
    robustness = {r['fluid']: r for r in result['robustness_rows']}
    order = [r['fluid'] for r in result['baseline_rows']]  # baseline ranking order

    fig, ax = plt.subplots(figsize=(9.5, 6))
    x = np.arange(len(order))
    w = 0.35
    b_cop = [baseline[f]['COP_electrical'] for f in order]
    r_cop = [robustness[f]['COP_electrical'] for f in order]
    ax.bar(x - w / 2, b_cop, w, color=COLOR_MAIN, edgecolor='white',
           label='Baseline point (5kg, 1Hz, 1.5T)')
    ax.bar(x + w / 2, r_cop, w, color=COLOR_POWER, edgecolor='white',
           label='Robustness point (2kg, 2Hz, 1.5T)')
    for xi, f in zip(x, order):
        ax.annotate(f"\u2212{baseline[f]['pct_below_water']:.1f}%", xy=(xi - w / 2, baseline[f]['COP_electrical']),
                    xytext=(0, 3), textcoords='offset points', ha='center', fontsize=7, color='#555')
    ax.set_xticks(x)
    ax.set_xticklabels([f.replace('_', ' ') for f in order], rotation=15, ha='right')
    ax.set_ylabel('Electrical COP (mdot re-optimized per fluid)')
    same = result['same_ranking']
    ax.set_title('Heat-Transfer Fluid Selection: COP Trade-off\n'
                 f"(% labels = COP giveup vs. pure water, not a real hardware option; "
                 f"ranking {'IS' if same else 'is NOT'} robust across both operating points)",
                 fontsize=10.5)
    ax.legend(fontsize=8.5)
    fig.tight_layout()
    save(fig, 'fig56_fluid_selection_comparison')


# ══════════════════════════════════════════════════════════════════════════
# FIG 57 — Commercial landscape reality check (model vs. vendor claims)
# ══════════════════════════════════════════════════════════════════════════

def plot_commercial_landscape_check():
    """Figure for commercial_landscape.py: two vendor-claimed magnetocaloric
    products (Magnotherm Stellar, Cooltech Applications' data-center unit)
    now exist in the commercial/trade-press record, alongside this repo's
    own naming-collision documentation (magnetic-bearing chillers,
    desiccant wheels, which dominate a plain-text search for "magnetic
    cooling" but are not magnetocaloric at all). No existing figure shows
    where this repo's own calibrated model sits relative to those public
    claims. Left panel: claimed capacity of the two real magnetocaloric
    commercial systems (Cooltech's own claimed COP range annotated, since
    Magnotherm's is undisclosed). Right panel: this repo's own model,
    run at a scaled-up Cooltech-class operating point using the SAME CORE
    loss-model calibration already flagged elsewhere as not calibrating
    against this device, against Cooltech's claimed COP range -- the gap
    is reported honestly, as an already-documented model limitation
    continuing to show up here, not a new discrepancy."""
    from core.commercial_landscape import COMMERCIAL_SYSTEMS
    model = commercial_landscape.model_prediction_at_cooltech_point(verbose=False)

    magneto = [c for c in COMMERCIAL_SYSTEMS if c.is_magnetocaloric]
    names = [c.name.replace(' Applications', '').replace(' (data-center-oriented unit)',
                                                           '\n(data-center unit)') for c in magneto]
    capacities = [c.claimed_capacity_kW or 0.0 for c in magneto]

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(11.5, 5.5))
    bars = ax1.bar(names, capacities, color=[COLOR_MAIN, COLOR_POWER], alpha=0.85, edgecolor='white')
    for b, c in zip(bars, magneto):
        label = f'{c.claimed_capacity_kW:.0f} kW'
        if c.claimed_COP_range:
            label += f'\nCOP claimed: {c.claimed_COP_range}'
        else:
            label += '\nCOP: undisclosed'
        ax1.annotate(label, xy=(b.get_x() + b.get_width() / 2, b.get_height()), xytext=(0, 4),
                     textcoords='offset points', ha='center', fontsize=9)
    ax1.set_ylabel('Claimed cooling capacity [kW]')
    ax1.set_title('Commercial Magnetocaloric-AMR Claims\n(trade-press/vendor, not independently audited)', fontsize=10.5)
    ax1.set_ylim(0, max(capacities) * 1.3)

    cooltech = next(c for c in magneto if c.name.startswith('Cooltech'))
    lo, hi = (float(x) for x in cooltech.claimed_COP_range.split('-'))
    labels2 = ['Cooltech claimed\n(trade press)', "This repo's model\n(CORE calibration,\nscaled to Cooltech-class point)"]
    vals2 = [(lo + hi) / 2, model['model_COP_electrical']]
    colors2 = [COLOR_MAIN, COLOR_POWER]
    bars2 = ax2.bar(labels2, vals2, color=colors2, alpha=0.85, edgecolor='white', width=0.55)
    ax2.errorbar([0], [(lo + hi) / 2], yerr=[[(lo + hi) / 2 - lo], [hi - (lo + hi) / 2]],
                 fmt='none', color='black', capsize=6, linewidth=1.5)
    # The Cooltech bar's label sits above the ERROR BAR's own top (hi), not
    # the bar height -- anchoring it at the bar height alone (as the model
    # bar's label does) would print the text directly on top of the
    # error-bar's vertical whisker line, since that whisker extends well
    # above the bar itself.
    for i, (b, v) in enumerate(zip(bars2, vals2)):
        anchor_y = hi if i == 0 else v
        ax2.annotate(f'{v:.2f}', xy=(b.get_x() + b.get_width() / 2, anchor_y), xytext=(0, 6),
                     textcoords='offset points', ha='center', fontsize=10, fontweight='bold')
    ax2.set_ylim(0, hi * 1.15)
    ax2.set_ylabel('Electrical COP')
    ax2.set_title(f"Model vs. Cooltech's Own Claim\ngap vs. claimed midpoint: "
                  f"{model['gap_pct_vs_claimed_midpoint']:+.0f}%\n"
                  f"(already-documented calibration limitation, not a new finding)",
                  fontsize=10.5)

    fig.suptitle('Commercial Landscape Reality Check: Public Claims vs. This Repo\'s Own Model',
                 fontsize=12)
    fig.tight_layout(rect=[0, 0, 1, 0.93])
    save(fig, 'fig57_commercial_landscape_check')


# ══════════════════════════════════════════════════════════════════════════
# FIG 58 — System-level cost of the Gd grain-Tc-broadening physics fix
# ══════════════════════════════════════════════════════════════════════════

def plot_calibrated_gd_system_impact(precomputed=None):
    """Figure for validation_system.run_calibrated_gd_system_level_
    comparison(): core/mce_material.py's more physically-complete Gd model
    (GADOLINIUM_CALIBRATED, with fitted grain-Tc broadening) is NOT wired
    in as this repo's system-wide default -- run_validation()'s own
    material-level comparison uses it, but amr_cycle.py/optimize.py/
    cascade.py all still use plain GADOLINIUM. No existing figure shows
    WHY that's a deliberate choice rather than an oversight: this plots
    the measured COP/Qc shift, side by side, for every already-calibrating
    Gd benchmark device -- and calls out the two devices that go from a
    real, positive Qc under plain Gd to a HARD ZERO under the broadened
    model at the identical mdot/span/field, which is concrete evidence
    (not just a performance-cost argument) for keeping the current
    default."""
    precomputed = precomputed or {}
    rows = precomputed.get('calibrated_gd_rows')
    if rows is None:
        result = validation_system.run_calibrated_gd_system_level_comparison(verbose=False)
        rows = result['rows']

    names = [r['device'].replace('_', '\n') for r in rows]
    names_flat = [r['device'].replace('_', ' ') for r in rows]
    cop_plain = [r['COP_plain'] for r in rows]
    cop_cal = [r['COP_calibrated'] for r in rows]
    zeroed = [r['Qc_calibrated_W'] == 0.0 and r['Qc_plain_W'] > 0 for r in rows]

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(13, 5.5))
    x = np.arange(len(names))
    w = 0.35
    ax1.bar(x - w / 2, cop_plain, w, label='Plain GADOLINIUM (current default)',
            color=COLOR_MAIN, alpha=0.85, edgecolor='white')
    bars_cal = ax1.bar(x + w / 2, cop_cal, w, label='GADOLINIUM_CALIBRATED (grain-Tc broadened)',
                        color=COLOR_POWER, alpha=0.85, edgecolor='white')
    for bar, z in zip(bars_cal, zeroed):
        if z:
            bar.set_hatch('xx')
    ax1.set_xticks(x)
    ax1.set_xticklabels(names, fontsize=7.5)
    ax1.set_ylabel('Electrical COP')
    ax1.set_title('COP: Plain vs. Broadened Gd\n(hatched = device goes to a HARD ZERO under the broadened model)',
                  fontsize=10)
    ax1.legend(fontsize=8)

    cop_shift = [r['COP_shift_pct'] for r in rows]
    colors_shift = [COLOR_POWER if z else COLOR_MAIN for z in zeroed]
    ax2.barh(names_flat, cop_shift, color=colors_shift, alpha=0.85, edgecolor='white')
    ax2.axvline(0, color='k', linewidth=0.8)
    ax2.tick_params(axis='y', labelsize=8)
    ax2.set_xlabel('COP shift vs. plain GADOLINIUM [%]  (-100% = zeroed out)')
    n_zeroed = sum(zeroed)
    ax2.set_title(f'System-Level Impact of the Physics Fix\n'
                  f'{n_zeroed}/{len(rows)} devices lose ALL cooling capacity at the same mdot/span/field',
                  fontsize=10)

    fig.suptitle('Why GADOLINIUM_CALIBRATED Is Not the System-Wide Default: the Measured Cost',
                 fontsize=12)
    fig.tight_layout(rect=[0, 0, 1, 0.93])
    save(fig, 'fig58_calibrated_gd_system_impact')


# ══════════════════════════════════════════════════════════════════════════
# FIG 59 — 1-D transient regenerator model vs. literature no-load spans
# ══════════════════════════════════════════════════════════════════════════

def plot_regenerator_1d_validation(precomputed=None):
    """Figure for regenerator_1d.validate_against_benchmarks(): the direct,
    honest validation of the 1-D transient regenerator model (the tool
    behind fig35/fig53's span-cap discussion) against the cleanest
    available ground truth -- directly-measured no-load spans, which
    require no mdot back-calibration. No existing figure shows this
    model's OWN accuracy; fig35 only shows whether its override recovers
    a nonzero prediction, not how close that prediction lands, and fig53
    only uses the older, cruder 0-D structural cap. This plots literature
    span against both the 1-D model's prediction AND the 0-D model's own
    structural ceiling for the same three devices, on one axis, making
    the honest finding directly visible: the 1-D model undershoots on two
    devices and overshoots on the third -- a genuine, direction-
    inconsistent error, not a single systematic bias to correct for."""
    precomputed = precomputed or {}
    rows = precomputed.get('regenerator_1d_rows')
    if rows is None:
        rows = regenerator_1d.validate_against_benchmarks(
            verbose=False, out_path=str(RESULTS_DIR / 'regenerator_1d_validation.txt'))

    names = [r['device'].replace('_', '\n') for r in rows]
    names_flat = [r['device'].replace('_', ' ') for r in rows]
    lit = [r['span_lit_K'] for r in rows]
    m1d = [r['span_1d_K'] for r in rows]
    cap0d = [r['cap_0d_K'] for r in rows]

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(13, 5.5))
    x = np.arange(len(names))
    w = 0.27
    ax1.bar(x - w, lit, w, label='Literature (directly measured)', color='#2ca02c',
            alpha=0.85, edgecolor='white')
    ax1.bar(x, cap0d, w, label='0-D structural cap (2\u00d7dTad_noload)', color=COLOR_MAIN,
            alpha=0.85, edgecolor='white')
    ax1.bar(x + w, m1d, w, label='1-D transient model', color=COLOR_POWER,
            alpha=0.85, edgecolor='white')
    ax1.set_xticks(x)
    ax1.set_xticklabels(names, fontsize=7.5)
    ax1.set_ylabel('No-load span [K]')
    ax1.set_title('Directly-Measured No-Load Span vs. Both Models', fontsize=10.5)
    ax1.legend(fontsize=8)

    err_1d = [r['err_1d_pct'] for r in rows]
    err_0d = [r['err_0d_pct'] for r in rows]
    w2 = 0.35
    ax2.barh(x - w2 / 2, err_0d, w2, label='0-D cap error (\u22640% by construction)',
             color=COLOR_MAIN, alpha=0.85, edgecolor='white')
    ax2.barh(x + w2 / 2, err_1d, w2, label='1-D model error', color=COLOR_POWER,
             alpha=0.85, edgecolor='white')
    ax2.axvline(0, color='k', linewidth=0.8)
    ax2.set_yticks(x)
    ax2.set_yticklabels(names_flat, fontsize=8)
    ax2.set_xlabel('Error vs. literature span [%]')
    ax2.set_title('Signed Error: the 1-D Model Undershoots on Two Devices,\n'
                  'Overshoots on the Third — Not a Single Systematic Bias', fontsize=10.5)
    ax2.legend(fontsize=8)

    fig.suptitle('1-D Transient Regenerator Model vs. Directly-Measured No-Load Spans\n'
                 '(genuine multi-cycle transient simulation — regenerative amplification is real, '
                 'but accuracy is direction-inconsistent)', fontsize=11.5)
    fig.tight_layout(rect=[0, 0, 1, 0.90])
    save(fig, 'fig59_regenerator_1d_validation')


# ══════════════════════════════════════════════════════════════════════════
# FIG 60 — Field-dependent Tc-broadening: fit vs. held-out data
# ══════════════════════════════════════════════════════════════════════════

def plot_field_dependent_broadening_calibration():
    """Figure for inhomogeneous_broadening.run_field_dependent_broadening_
    calibration(): a DIFFERENT broadening question from fig27's general
    Gaussian sigma sweep -- here a specific, literature-motivated
    field-dependent form (sigma_Tc = k*mu0*H, motivated by Dan'kov et
    al.'s own qualitative field-broadening discussion) is fit to the same
    3 calibration points core/validation.py already uses, THEN checked
    against genuinely held-out data (Dan'kov's own 7.5T point and
    Giguere et al.'s independent 5T/7T Gd cross-check) that did not feed
    the fit. No existing figure shows this held-out check. Reported
    honestly: the least-squares fit converges to k=0 (no broadening
    preferred), so the "fitted" and "sharp" curves are the same line here
    -- itself a genuine negative result worth showing rather than hiding
    behind only a paragraph of text."""
    from core.inhomogeneous_broadening import (calibrate_field_dependent_broadening,
                                                 FieldBroadenedMagnetocaloricMaterial)
    from core.validation import GIGUERE_GD_CROSSCHECK

    fit = calibrate_field_dependent_broadening(verbose=False)
    k_fit = fit['k_fit_K_per_T']
    mat_fit = FieldBroadenedMagnetocaloricMaterial(GADOLINIUM, k_fit)

    fitted_rows = fit['rows']
    fitted_B = [r['mu0H_T'] for r in fitted_rows]
    fitted_lit = [r['dT_lit_K'] for r in fitted_rows]
    fitted_sharp = [r['dT_before_K'] for r in fitted_rows]
    fitted_k = [r['dT_after_K'] for r in fitted_rows]

    B75, dT_lit_75 = 7.5, 15.5
    H75 = B75 / mu0
    dT_sharp_75 = float(np.asarray(GADOLINIUM.delta_T_adiabatic(np.array([294.0]), H75)).ravel()[0])
    dT_k_75 = float(np.asarray(mat_fit.delta_T_adiabatic(np.array([294.0]), H75)).ravel()[0])

    held_out_B = [B75]
    held_out_lit_mid = [dT_lit_75]
    held_out_sharp = [dT_sharp_75]
    held_out_k = [dT_k_75]
    held_out_labels = ["Dan'kov 7.5T\n(Fig.10 pixel-read)"]
    for B, ref in sorted(GIGUERE_GD_CROSSCHECK.items()):
        H = B / mu0
        lo, hi = ref['range_K']
        mid = 0.5 * (lo + hi)
        dT_sharp = float(np.asarray(GADOLINIUM.delta_T_adiabatic(np.array([294.0]), H)).ravel()[0])
        dT_k = float(np.asarray(mat_fit.delta_T_adiabatic(np.array([294.0]), H)).ravel()[0])
        held_out_B.append(B)
        held_out_lit_mid.append(mid)
        held_out_sharp.append(dT_sharp)
        held_out_k.append(dT_k)
        held_out_labels.append(f'Gigu\u00e8re {B:.0f}T\n({lo:.1f}-{hi:.1f}K range)')

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(13, 5.5))

    x1 = np.arange(len(fitted_B))
    w = 0.27
    ax1.bar(x1 - w, fitted_lit, w, label='Literature (fitted points)', color='#2ca02c',
            alpha=0.85, edgecolor='white')
    ax1.bar(x1, fitted_sharp, w, label='Sharp (k=0)', color=COLOR_MAIN, alpha=0.85, edgecolor='white')
    ax1.bar(x1 + w, fitted_k, w, label=f'Fitted k={k_fit:.3f} K/T', color=COLOR_POWER,
            alpha=0.85, edgecolor='white')
    ax1.set_xticks(x1)
    ax1.set_xticklabels([f'{b:.0f} T' for b in fitted_B])
    ax1.set_ylabel(r'$\Delta T_{ad}$ [K] at 294K')
    ax1.set_title(f'3-Point Fit (SAME points core/validation.py uses)\n'
                  f'k converges to {k_fit:.3f} K/T — essentially no broadening preferred',
                  fontsize=10)
    ax1.legend(fontsize=8)

    x2 = np.arange(len(held_out_B))
    ax2.bar(x2 - w, held_out_lit_mid, w, label='Literature (held-out, NOT in the fit)',
            color='#2ca02c', alpha=0.85, edgecolor='white')
    ax2.bar(x2, held_out_sharp, w, label='Sharp (k=0)', color=COLOR_MAIN, alpha=0.85, edgecolor='white')
    ax2.bar(x2 + w, held_out_k, w, label='Fitted k', color=COLOR_POWER, alpha=0.85, edgecolor='white')
    ax2.set_xticks(x2)
    ax2.set_xticklabels(held_out_labels, fontsize=8)
    ax2.set_ylabel(r'$\Delta T_{ad}$ [K] at 294K')
    ax2.set_title('Held-Out Cross-Check\n(sharp and fitted-k bars coincide: k\u22480 changes nothing)',
                  fontsize=10)
    # Dan'kov (the tallest category) sits at the LEFT, directly under
    # where a default/'best' legend placement lands, clipping the tops
    # of its Sharp/Fitted-k bars behind the legend box. Add headroom
    # above the tallest bar so the legend has clear whitespace above the
    # data instead of sitting on top of it.
    ax2.set_ylim(0, max(held_out_lit_mid + held_out_sharp + held_out_k) * 1.3)
    ax2.legend(fontsize=8)

    # tight_layout(rect=[...]) reserves a FIXED amount of vertical space
    # for the suptitle regardless of the rect's own top value (confirmed:
    # axes top tracks rect_top - 0.153 exactly), which left a large blank
    # band between the suptitle and the subplot titles that bbox_inches=
    # 'tight' (used by save()) can't crop since it only trims the outer
    # figure margin, not this internal gap. tight_layout() with no rect,
    # followed by an explicit subplots_adjust(top=...), sizes the axes
    # tightly under both title rows with no leftover band.
    fig.tight_layout()
    fig.suptitle('Field-Dependent Tc-Broadening (\u03c3_Tc = k\u00b7\u03bc\u2080H): a Literature-Motivated '
                 'Fit That Converges to No Broadening', fontsize=12, y=0.99)
    fig.subplots_adjust(top=0.80)
    save(fig, 'fig60_field_dependent_broadening_calibration')


# ══════════════════════════════════════════════════════════════════════════
# FIG 61 — Magnet-geometry Pareto sensitivity: multiseed stability check
# ══════════════════════════════════════════════════════════════════════════

def _parse_magnet_geometry_multiseed_txt(path):
    """Fallback parser for the checked-in results/magnet_geometry_
    multiseed_stability.txt table (used when pymoo is unavailable so this
    figure doesn't have to re-run 6 NSGA-III optimizations live)."""
    rows = []
    started = False
    with open(path, encoding='utf-8', errors='replace') as f:
        for line in f:
            s = line.strip()
            if s.startswith('seed') and 'mean_FLAT_T' in s:
                started = True
                continue
            if not started:
                continue
            if s.startswith('---') or not s:
                if rows:
                    break
                continue
            parts = s.split()
            if len(parts) == 5:
                try:
                    seed, flat, geom, ff, fg = parts
                    rows.append({'seed': int(seed), 'mean_flat_T': float(flat),
                                 'mean_geometric_T': float(geom),
                                 'front_size_flat': int(ff), 'front_size_geometric': int(fg)})
                except ValueError:
                    break
    return rows


def plot_magnet_geometry_multiseed_stability(precomputed=None):
    """Figure for magnet_geometry.run_magnet_geometry_multiseed_stability_
    check(): fig34 shows the FLAT-vs-GEOMETRIC magnet-mass Pareto
    comparison at a single seed/reduced NSGA-III setting. This is the
    direct robustness follow-up that single-seed result explicitly
    invited (Paper-Mining Pass review item 4): does the "geometric cost
    pulls the mean field down" direction hold at production pop_size/
    n_gen settings across multiple independent seeds? The honest answer,
    plotted directly rather than left in a text report: NO -- at least
    one seed's GEOMETRIC run has a HIGHER mean field than its own FLAT
    run, the same kind of seed-sensitivity fig49 already found for the
    main material/geometry Pareto front. This is new information no
    existing figure shows: it is the specific reliability check on
    fig34's own headline claim, not a repeat of fig34 itself.

    precomputed, if given, may supply 'magnet_geometry_multiseed_result'
    (the dict returned by run_magnet_geometry_multiseed_stability_check())
    to reuse an already-computed run instead of re-running 6 NSGA-III
    optimizations (production pop_size=40/n_gen=25 settings, ~1 minute)
    a second time. When called standalone with pymoo unavailable, this
    falls back to the pre-computed results/magnet_geometry_multiseed_
    stability.txt table checked into the repo."""
    precomputed = precomputed or {}
    result = precomputed.get('magnet_geometry_multiseed_result')
    if result is not None:
        per_seed = result['per_seed']
    elif HAVE_PYMOO:
        result = magnet_geometry.run_magnet_geometry_multiseed_stability_check(
            out_path=None)
        per_seed = result['per_seed']
    else:
        print("  [pymoo unavailable — falling back to pre-computed "
              "results/magnet_geometry_multiseed_stability.txt]")
        per_seed = _parse_magnet_geometry_multiseed_txt(
            RESULTS_DIR / 'magnet_geometry_multiseed_stability.txt')

    seeds = [s['seed'] for s in per_seed]
    mean_flat = [s['mean_flat_T'] for s in per_seed]
    mean_geom = [s['mean_geometric_T'] for s in per_seed]
    reversed_seeds = [g > f for f, g in zip(mean_flat, mean_geom)]

    fig, ax = plt.subplots(figsize=(8, 5.5))
    x = np.arange(len(seeds))
    w = 0.35
    ax.bar(x - w / 2, mean_flat, w, label='FLAT (previous cost term)', color=COLOR_MAIN,
           alpha=0.85, edgecolor='white')
    bars_geom = ax.bar(x + w / 2, mean_geom, w, label='GEOMETRIC (Halbach magnet-mass cost)',
                        color=COLOR_POWER, alpha=0.85, edgecolor='white')
    for bar, rev in zip(bars_geom, reversed_seeds):
        if rev:
            bar.set_hatch('//')
    ax.set_xticks(x)
    ax.set_xticklabels([f'seed {s}' for s in seeds])
    ax.set_ylabel('Merged Pareto front mean mu0H_max [T]')
    n_rev = sum(reversed_seeds)
    ax.set_title('Does "Geometric Cost Pulls Mean Field Down" Hold at Production Settings?\n'
                 f'{n_rev}/{len(seeds)} seed(s) reverse the expected direction '
                 '(hatched bar = GEOMETRIC mean field HIGHER than FLAT)',
                 fontsize=11)
    # The default 'best' location sat directly on top of the tallest bar
    # (seed 2's GEOMETRIC bar came closest to the field ceiling), clipping
    # its hatching -- add headroom above the tallest bar so the legend has
    # clear whitespace to sit in instead of overlapping data.
    ax.set_ylim(0, max(mean_flat + mean_geom) * 1.22)
    ax.legend(fontsize=9)
    fig.tight_layout()
    save(fig, 'fig61_magnet_geometry_multiseed_stability')


# ══════════════════════════════════════════════════════════════════════════
# Generate all figures
# ══════════════════════════════════════════════════════════════════════════

# ══════════════════════════════════════════════════════════════════════════
# Generate all / selected figures
# ══════════════════════════════════════════════════════════════════════════

def _build_figure_registry(precomputed=None):
    """Builds the full (label, fn) figure registry used by both run_all()
    (which runs every entry) and run_figures()/the CLI (which runs a
    filtered subset). Kept as its own function so a specific figure can be
    regenerated -- e.g. after a fix, or just to check one result -- without
    re-running the other 60, most of which are cheap but a handful of
    which (fig18/33/34/49/61's NSGA-III sweeps) are not.

    precomputed is a dict that may carry results already computed earlier
    in the SAME pipeline run (main.py steps 2/2f/4/7/7b/7c/8d/9/9b/11/
    11b/11d), reused instead of recomputing them:
      - 'system_validation_results' (step 2)
      - 'baseline_rows' (step 4)
      - 'sobol_const_Si', 'sobol_state_Si' (step 9 / step 9b)
      - 'pareto_rows' (step 11)
      - 'cascade_rows_gd', 'cascade_rows_giant' (step 7)
      - 'graded_rows' (step 7b)
      - 'astro_result' (step 7c)
      - 'material_rows' (step 8d)
      - 'hysteresis_result' (step 11b)
      - 'magnet_geometry_result' (step 11d)
      - 'override_check_result' (step 2f)
      - 'fluid_selection_result' (fig56)
      - 'calibrated_gd_rows' (fig58)
      - 'regenerator_1d_rows' (fig59)
      - 'magnet_geometry_multiseed_result' (fig61)
    Figures 8, 14, 16, 18, 19, 20, 21, 25, 26, 33, 34, 35, 47, 56, 58, 59,
    and 61 reuse whichever of these are supplied instead of recomputing
    them, which otherwise adds several minutes of redundant work. When
    called standalone (`python plots.py`, precomputed=None) every figure
    still computes its own data fresh from core/, exactly as before."""
    precomputed = precomputed or {}

    figure_fns = [
        ("01 Gd MCE validation vs. Dan'kov et al. (1998)", plot_gd_validation),
        ("02 Gd entropy change & DeltaT_ad vs. T", plot_gd_entropy_dTad),
        ("03 Gd5Si2Ge2 Landau model calibration", plot_landau_giant_mce),
        ("04 Giguere et al. (1999) direct-measurement cross-check", plot_giguere_validation),
        ("05 Material comparison (Gd / Gd5Si2Ge2 / La(Fe,Si)13Hy)", plot_material_comparison),
        ("06 AMR characteristic curve", plot_amr_characteristic_curve),
        ("07 AMR energy balance vs. span", plot_amr_energy_balance),
        ("08 AMR vs. baselines COP comparison",
         lambda: plot_amr_vs_baselines(precomputed)),
        ("09 NTU regenerator effectiveness", plot_regenerator_effectiveness),
        ("10 Geometry trade-off — packed bed", plot_geometry_packed_bed),
        ("11 Geometry trade-off — parallel plate", plot_geometry_parallel_plate),
        ("12 Loss-model calibration + leave-one-out CV", plot_loss_model_calibration),
        ("13 Parasitic-fraction scaling with device size", plot_parasitic_fraction_scaling),
        ("14 System-level validation vs. published prototypes",
         lambda: plot_system_validation(precomputed)),
        ("15 Curve-level (2-point) Qc(span) validation", plot_curve_validation),
        ("16 Sobol global sensitivity analysis",
         lambda: plot_sobol_sensitivity(precomputed)),
        ("17 RSM surrogate parity + coefficients", plot_rsm_surrogate),
        ("18 NSGA-III Pareto front",
         lambda: plot_nsga3_pareto(precomputed)),
        ("19 Cascade staging vs. baselines (Gd)",
         lambda: plot_cascade_staging_gd(precomputed)),
        ("20 Cascade: Gd vs. Gd5Si2Ge2",
         lambda: plot_cascade_giant_vs_gd(precomputed)),
        ("21 Curie-graded cascade performance",
         lambda: plot_graded_cascade(precomputed)),
        ("22 Economics: TCO and lifetime cost", plot_economics),
        ("23 Emissions comparison", plot_emissions),
        ("24 Giant-MCE targeting comparison", plot_giant_mce_targeting),
        ("25 Astronautics graded-bed validation",
         lambda: plot_astronautics_validation(precomputed)),
        ("26 Material family comparison (Track A2 item)",
         lambda: plot_material_family_comparison(precomputed)),
        ("27 Inhomogeneous/polycrystalline Tc-broadening sensitivity ",
         plot_inhomogeneous_broadening),
        ("28 Nanocomposite off-design robustness ( follow-up)",
         plot_nanocomposite_robustness),
        ("29 Mechanical-contact thermal-diode sensitivity ",
         plot_thermal_diode_sensitivity),
        ("30 Magnetocaloric-fluid volume-fraction sweep ",
         plot_fluid_mce_sweep),
        ("31 Passive magnetic-regenerator alignment effect ",
         plot_passive_regenerator_alignment),
        ("32 Rotary-device cycle-type validation ",
         plot_cycle_type_validation),
        ("33 Hysteresis-loss Pareto-front sensitivity ",
         lambda: plot_hysteresis_sensitivity(precomputed)),
        ("34 Magnet-geometry (Halbach) Pareto-front sensitivity ",
         lambda: plot_magnet_geometry_pareto_sensitivity(precomputed)),
        ("35 Regenerative-amplification override check (this session)",
         lambda: plot_regenerative_amplification_override_check(precomputed)),
        ("36 Value proposition: COP crossover, AMR vs. VCC across spans",
         plot_value_proposition_cop_crossover),
        ("37 Value proposition: emissions breakdown, AMR vs. VCC vs. liquid",
         plot_value_proposition_emissions_breakdown),
        ("38 Value proposition: refrigerant GWP comparison",
         plot_value_proposition_refrigerant_gwp),
        ("39 Value proposition: real-deployment comparison",
         plot_value_proposition_real_deployments),
        ("40 MnFePSi-doped hysteresis (speculative)",
         plot_mnfepsi_doped_hysteresis_speculative),
        ("41 MnFePSi hysteresis k-fit quality",
         plot_mnfepsi_hysteresis_kfit_quality),
        ("42 Hysteresis-exploiting actuator cost estimate",
         plot_hysteresis_exploiting_actuator_estimate),
        ("43 Hysteresis-loss landscape",
         plot_hysteresis_loss_landscape),
        ("44 Phase 37: CALIBRATION_POINTS_CORE drift, found and fixed",
         plot_calibration_drift_phase37),
        ("45 Which Lozano rows calibrate?",
         plot_lozano_calibration_status),
        ("46 Alternative caloric technologies vs. this repo's VCC COP",
         plot_alternative_caloric_comparison),
        ("47 Annual water-usage (WUE) comparison",
         lambda: plot_water_usage_comparison(precomputed)),
        ("48 Monte Carlo calibration-uncertainty band on COP_electrical",
         plot_uncertainty_propagation),
        ("49 NSGA-III Pareto-front seed-to-seed stability",
         plot_pareto_multiseed_stability),
        ("50 Regime-crossover search: does AMR ever beat VCC anywhere?",
         plot_regime_crossover_analysis),
        ("51 Ames Lab heat-pump architecture validation",
         plot_heat_pump_validation),
        ("52 Hypereg parallel-hydraulic pumping-power sensitivity",
         plot_hypereg_analysis),
        ("53 Regenerative-amplification gap across the full benchmark corpus",
         plot_regenerative_amplification_gap),
        ("54 Hybrid solid-state regenerator (HMR) vs. VCC, ideal -> real",
         plot_hmr_vs_vcc_loss_funnel),
        ("55 Annualized, climate-weighted COP comparison",
         plot_annualized_climate_comparison),
        ("56 Heat-transfer fluid selection comparison",
         lambda: plot_fluid_selection_comparison(precomputed)),
        ("57 Commercial landscape reality check",
         plot_commercial_landscape_check),
        ("58 System-level cost of the Gd grain-Tc-broadening physics fix",
         lambda: plot_calibrated_gd_system_impact(precomputed)),
        ("59 1-D transient regenerator model vs. literature no-load spans",
         lambda: plot_regenerator_1d_validation(precomputed)),
        ("60 Field-dependent Tc-broadening: fit vs. held-out data",
         plot_field_dependent_broadening_calibration),
        ("61 Magnet-geometry Pareto sensitivity: multiseed stability",
         lambda: plot_magnet_geometry_multiseed_stability(precomputed)),
    ]

    return figure_fns


def _execute_figures(figure_fns):
    """Runs each (label, fn) pair in figure_fns, printing progress and
    catching exceptions per-figure so one failure doesn't abort the rest.
    Returns the list of labels that failed."""
    failures = []
    for label, fn in figure_fns:
        print(f"[{label}]")
        try:
            fn()
        except Exception as exc:  # noqa: BLE001 - keep generating the rest
            print(f"  !!! FAILED: {label}: {exc}")
            failures.append(label)
        print()
    return failures


def _report_failures(failures):
    print(f"Figures written to: {FIG_DIR.resolve()}")
    if failures:
        print(f"\n{len(failures)} figure(s) failed:")
        for label in failures:
            print(f"  - {label}")
    if not HAVE_SALIB:
        print("\nNote: SALib not installed — fig16 used pre-computed results/sobol_results*.txt.")
    if not HAVE_PYMOO:
        print("Note: pymoo not installed — fig18 used pre-computed results/pareto_front.csv.")


def _figure_number(label):
    """Extracts the leading figure number from a figure_fns label, e.g.
    '07 AMR energy balance vs. span' -> 7. Returns None if the label
    doesn't start with a number (shouldn't happen for any entry in
    _build_figure_registry(), but handled defensively anyway)."""
    m = re.match(r'\s*0*(\d+)', label)
    return int(m.group(1)) if m else None


def list_available_figures(precomputed=None):
    """Prints every figure's number and label, one per line, in registry
    order -- what `python plots.py --list` shows. Useful for finding the
    right number/substring to pass to run_figures() or the CLI without
    reading this whole file."""
    for label, _ in _build_figure_registry(precomputed):
        print(f"  {label}")


def _resolve_selection(tokens, figure_fns):
    """Resolves user-supplied tokens against figure_fns, matching each
    token against, in order of precedence:
      1. an inclusive numeric range, 'N-M' (e.g. '53-56')
      2. a single figure number, with or without a 'fig' prefix or
         leading zeros ('53', 'fig53', '053')
      3. a case-insensitive substring of the figure's own label (e.g.
         'hysteresis' matches every figure whose label mentions it)
    Returns the matching (label, fn) pairs in the registry's own order,
    de-duplicated (a figure matched by two different tokens is only run
    once). Raises ValueError naming any token that matched nothing, so a
    typo doesn't just silently generate zero figures."""
    numbered = [(_figure_number(label), label, fn) for label, fn in figure_fns]
    selected_indices = []
    unmatched = []

    def _add(i):
        if i not in selected_indices:
            selected_indices.append(i)

    for token in tokens:
        t = str(token).strip()
        matched = False

        range_match = re.match(r'^(\d+)\s*-\s*(\d+)$', t)
        if range_match:
            lo, hi = int(range_match.group(1)), int(range_match.group(2))
            for i, (num, _label, _fn) in enumerate(numbered):
                if num is not None and lo <= num <= hi:
                    _add(i)
                    matched = True
            if not matched:
                unmatched.append(token)
            continue

        num_match = re.match(r'^(?:fig)?0*(\d+)$', t, re.IGNORECASE)
        if num_match:
            want = int(num_match.group(1))
            for i, (num, _label, _fn) in enumerate(numbered):
                if num == want:
                    _add(i)
                    matched = True
                    break
            if not matched:
                unmatched.append(token)
            continue

        t_lower = t.lower()
        for i, (_num, label, _fn) in enumerate(numbered):
            if t_lower in label.lower():
                _add(i)
                matched = True
        if not matched:
            unmatched.append(token)

    if unmatched:
        raise ValueError(
            f"No figure matched: {', '.join(repr(u) for u in unmatched)}. "
            f"Run with --list to see every available figure number/label."
        )
    return [(numbered[i][1], numbered[i][2]) for i in selected_indices]


def run_figures(selection, precomputed=None):
    """Generates only the requested figure(s) instead of the full
    run_all() sweep -- useful for regenerating one figure after a fix, or
    checking a single result, without waiting on the NSGA-III-backed
    figures (18/33/34/49/61) that dominate a full run's time.

    `selection` is a list/tuple of tokens; each one may be:
      - a figure number, as an int or string, with or without a 'fig'
        prefix or leading zeros: 53, '53', 'fig53', '053'
      - an inclusive number range: '53-56'
      - a case-insensitive substring of a figure's own label, e.g.
        'hysteresis' or 'commercial'
    A single bare token (not in a list) is also accepted for convenience,
    e.g. run_figures(53) or run_figures('hysteresis').

    precomputed works exactly as in run_all() -- pass the same dict to
    reuse already-computed results for whichever selected figures accept
    it. Returns the list of any selected figure labels that failed (same
    convention as run_all()); an unresolvable selection raises ValueError
    up front, before anything is computed.
    """
    if isinstance(selection, (str, int)):
        selection = [selection]
    figure_fns = _build_figure_registry(precomputed)
    chosen = _resolve_selection(selection, figure_fns)
    print(f"Generating {len(chosen)} selected figure(s):")
    for label, _fn in chosen:
        print(f"  [{label}]")
    print()
    failures = _execute_figures(chosen)
    _report_failures(failures)
    return failures


def run_all(precomputed=None):
    """Generates all 61 figures. See _build_figure_registry()'s docstring
    for what may be passed via `precomputed`. To generate only some of
    the figures, use run_figures(selection, precomputed) instead, or run
    `python plots.py <figure(s)>` from the command line."""
    print("Generating all figures...\n")
    figure_fns = _build_figure_registry(precomputed)
    failures = _execute_figures(figure_fns)
    _report_failures(failures)
    return failures


if __name__ == '__main__':
    import argparse

    parser = argparse.ArgumentParser(
        description="Generate magcool-dc figures. With no arguments, generates the "
                    "full figure set (equivalent to calling run_all()). Pass one or "
                    "more figure numbers, ranges, or name substrings to generate only "
                    "those instead.")
    parser.add_argument(
        'figures', nargs='*',
        help="Figure(s) to generate, e.g. 53, fig53, 53-56, or a label substring "
             "like hysteresis. Space-separated for multiple. Omit to generate every "
             "figure.")
    parser.add_argument(
        '--list', action='store_true',
        help="List every available figure's number and label, then exit without "
             "generating anything.")
    args = parser.parse_args()

    if args.list:
        list_available_figures()
    elif args.figures:
        try:
            run_figures(args.figures)
        except ValueError as exc:
            print(f"Error: {exc}")
            sys.exit(1)
    else:
        run_all()