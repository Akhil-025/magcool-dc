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
from core.baseline_cooling import vapor_compression_cop, liquid_cooling_cop
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
    print(f"  Saved -> {path}.png / .pdf")


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
    for xi, e in zip(x, err):
        ax2.text(xi, e + (1.5 if e >= 0 else -3.0), f'{e:+.1f}%',
                  ha='center', fontsize=9)

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
    ax.legend(fontsize=9)
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

    T_cold_K = 291.15
    if baseline_rows is not None:
        spans = [r['span_K'] for r in baseline_rows]
        cop_e = [r['AMR_COP_electrical'] for r in baseline_rows]
        vcc_l = [r['VaporCompression_COP'] for r in baseline_rows]
        liq_l = [r['LiquidCooling_COP'] for r in baseline_rows]
        carnot_l = [r['Carnot_COP'] for r in baseline_rows]
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

    fig, ax = plt.subplots(figsize=(8, 5.5))
    ax.plot(spans, cop_e, color=COLOR_POWER, marker='o', label='Magnetic (AMR) — electrical COP')
    ax.plot(spans, vcc_l, color=COLOR_MAIN, marker='s', label='Vapor-compression')
    ax.plot(spans, liq_l, color='#85bb65', marker='^', label='Liquid cooling (blended)')
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
# FIG 10 / 11 — Geometry-dependent pumping power trade-offs (Phase 7)
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
    for qc_, f_, n_ in zip(Qc, frac, names):
        ax.annotate(n_, (qc_, f_), fontsize=7, xytext=(5, 4), textcoords='offset points')
    ax.set_xscale('log')
    ax.set_xlabel('Device cooling capacity Qc [W] (log)')
    ax.set_ylabel('Parasitic fraction  W_parasitic / Qc')
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

    fig, axes = plt.subplots(1, 2, figsize=(13, 5.5))
    axes[0].scatter(cop_lit, cop_model, s=90, color=COLOR_MAIN, zorder=3, edgecolor='white')
    lims = [0, max(max(cop_lit), max(cop_model)) * 1.15]
    axes[0].plot(lims, lims, 'k--', linewidth=1, label='Perfect agreement')
    for x_, y_, n_ in zip(cop_lit, cop_model, names):
        axes[0].annotate(n_, (x_, y_), fontsize=7, xytext=(5, 5), textcoords='offset points')
    axes[0].set_xlim(lims)
    axes[0].set_ylim(lims)
    axes[0].set_xlabel('Literature COP (electrical)')
    axes[0].set_ylabel('Model COP (electrical)')
    axes[0].set_title('System-Level Validation: Model vs. Literature\n'
                       '(fluid mdot calibrated to reproduce reported Qc)')
    axes[0].legend(fontsize=9)

    colors = [COLOR_POWER if abs(e) > 30 else COLOR_MAIN for e in err]
    axes[1].barh(names, err, color=colors, alpha=0.85, edgecolor='white')
    axes[1].axvline(0, color='k', linewidth=0.8)
    axes[1].set_xlabel('COP error [%]')
    axes[1].set_title('COP Prediction Error by Device')
    for ax in axes:
        ax.tick_params(labelsize=8)
    fig.tight_layout()
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
    ax.bar(x - w / 2, lit, w, label='Literature (companion span point)',
           color=COLOR_MAIN, alpha=0.85, edgecolor='white')
    ax.bar(x + w / 2, model, w, label='Model (predicted from anchor-point calibration)',
           color=COLOR_POWER, alpha=0.85, edgecolor='white')
    ax.set_xticks(x)
    ax.set_xticklabels(names, rotation=10, ha='right', fontsize=9)
    ax.set_ylabel('Qc at companion span [W]')
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
            out_path=str(RESULTS_DIR / 'sobol_results_phase2_constant.txt'),
            use_state_dependent_losses=False)
        Si_state = sensitivity_mod.run_sobol(
            out_path=str(RESULTS_DIR / 'sobol_results.txt'),
            use_state_dependent_losses=True)
        names = sensitivity_mod.PROBLEM['names']
        st_const = dict(zip(names, Si_const['ST']))
        st_state = dict(zip(names, Si_state['ST']))
    else:
        print("  [SALib unavailable — falling back to pre-computed results/sobol_results*.txt]")
        st_const = _parse_sobol_txt(RESULTS_DIR / 'sobol_results_phase2_constant.txt')
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
            ax.axvline(first_oor, color='#85bb65', linestyle=':', alpha=0.5, linewidth=1)
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
        ax2.legend(fontsize=9)
    else:
        for ax in axes:
            ax.text(0.5, 0.5, astro.get('status', 'infeasible'), ha='center', va='center',
                    wrap=True, fontsize=10, transform=ax.transAxes)
            ax.set_xticks([])
            ax.set_yticks([])

    fig.suptitle('Phase 9: 6-Layer Curie-Graded La(Fe,Si)13Hy Bed vs. '
                 'Real Astronautics_rotary_2014 Device', fontsize=12)
    fig.tight_layout()
    save(fig, 'fig25_astronautics_graded_bed_validation')


# ══════════════════════════════════════════════════════════════════════════
# FIG 26 — Four-way material family comparison (Track A2 item)
# ══════════════════════════════════════════════════════════════════════════

def plot_material_family_comparison(precomputed=None):
    """precomputed, if given, may supply 'material_rows' (already computed
    by step 8d in main.py) so this figure reuses it instead of re-running
    the ~6s four-way material-family sweep a second time in the same
    pipeline invocation."""
    precomputed = precomputed or {}
    rows = precomputed.get('material_rows')
    if rows is None:
        rows = material_family_comparison.build_comparison_table()
    rep = [r for r in rows if r['span_K'] == material_family_comparison.REPRESENTATIVE_SPAN_K]

    labels = [r['candidate'].replace(' (', '\n(') for r in rep]
    colors5 = [COLOR_MAIN, COLOR_POWER, '#85bb65', '#e8a33d', '#7b52ab']
    qc_vals = [r['1stage_Qc_W'] or 0.0 for r in rep]
    cop_vals = [r['1stage_COP'] or 0.0 for r in rep]
    hatches = ['' if r['in_range'] else '//' for r in rep]

    fig, axes = plt.subplots(1, 2, figsize=(12, 5.5))
    bars0 = axes[0].bar(labels, qc_vals, color=colors5, alpha=0.85, edgecolor='white')
    bars1 = axes[1].bar(labels, cop_vals, color=colors5, alpha=0.85, edgecolor='white')
    for bars in (bars0, bars1):
        for bar, hatch in zip(bars, hatches):
            bar.set_hatch(hatch)
    axes[0].set_ylabel('Qc [W]')
    axes[0].set_title('Cooling Capacity (1-stage)')
    axes[1].set_ylabel('Electrical COP')
    axes[1].set_title('Electrical COP (1-stage)')
    for ax in axes:
        ax.tick_params(axis='x', labelsize=8)

    span = material_family_comparison.REPRESENTATIVE_SPAN_K
    T_cold = material_family_comparison.T_COLD_K
    fig.suptitle('Material Family Comparison at the ASHRAE Point '
                 f'(T_cold={T_cold:.0f}K, span={span:.0f}K)\n'
                 'Hatched bars: family\'s documented Tc window does not cover this point '
                 '(fell back to Gd)', fontsize=11)
    fig.tight_layout()
    save(fig, 'fig26_material_family_comparison')


# ══════════════════════════════════════════════════════════════════════════
# Generate all figures
# ══════════════════════════════════════════════════════════════════════════

def run_all(precomputed=None):
    """Generates all 26 figures.

    precomputed, if given, is a dict that may carry results already
    computed earlier in the SAME pipeline run (main.py steps 2/4/7/7b/7c/8d/9/9b/11):
      - 'system_validation_results'         (step 2)
      - 'baseline_rows'                     (step 4)
      - 'sobol_const_Si', 'sobol_state_Si'  (step 9 / step 9b)
      - 'pareto_rows'                       (step 11)
      - 'cascade_rows_gd', 'cascade_rows_giant'  (step 7)
      - 'graded_rows'                       (step 7b)
      - 'astro_result'                      (step 7c)
      - 'material_rows'                     (step 8d)
    Figures 8, 14, 16, 18, 19, 20, 21, 25, and 26 reuse whichever of these
    are supplied instead of recomputing them, which otherwise adds several
    minutes of redundant work (mostly fig21's graded-cascade sweep and
    fig25's Astronautics validation) on top of the identical computation
    already performed earlier in the same run. When called standalone
    (`python plots.py`, precomputed=None) every figure still computes its
    own data fresh from core/, exactly as before.
    """
    print("Generating all figures...\n")
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
    ]

    failures = []
    for label, fn in figure_fns:
        print(f"[{label}]")
        try:
            fn()
        except Exception as exc:  # noqa: BLE001 - keep generating the rest
            print(f"  !!! FAILED: {label}: {exc}")
            failures.append(label)
        print()

    print(f"All figures written to: {FIG_DIR.resolve()}")
    if failures:
        print(f"\n{len(failures)} figure(s) failed:")
        for label in failures:
            print(f"  - {label}")
    if not HAVE_SALIB:
        print("\nNote: SALib not installed — fig16 used pre-computed results/sobol_results*.txt.")
    if not HAVE_PYMOO:
        print("Note: pymoo not installed — fig18 used pre-computed results/pareto_front.csv.")


if __name__ == '__main__':
    run_all()