"""
hypereg_analysis.py
====================
Phase 15 addition: demonstrates the pumping-power effect of a Klinar et
al. (2024) "Hypereg"-style parallel-hydraulic regenerator split (see
`results/hypereg_findings.md` for the full literature findings note this
module implements the code side of, and `core.thermal.
pumping_power_packed_bed_hypereg()` for the underlying model).

Scope, deliberately narrow (read the findings note before trusting any
number here)
--------------------------------------------------------------------------
The Klinar et al. (2024) review gives ONE illustrative example
(n_parallel_subregenerators=4, Fig. 19) and no validated pressure-drop
data -- there is no published pressure-drop-reduction curve to fit or
reproduce. So this module does NOT claim an optimum `n`; it shows, at
this repo's own representative operating point, how COP_electrical and
the achievable no-pumping-penalty frequency shift as `n` is swept,
qualitatively confirming (or not) that the mechanism is worth pursuing
further -- the same "does a genuine trade-off/benefit exist in this
repo's own model" question `geometry_analysis.py` already asks for
particle diameter and plate spacing, applied here to Hypereg's parallel-
split hydraulic idea instead.
"""

from core.mce_material import GADOLINIUM
from core.amr_cycle import AMRSystem
from core.loss_model import StateDependentLossModel

T_COLD_K = 291.0
SPAN_K = 10.0
MU0H_T = 1.5
MASS_KG = 5.0
MDOT_KG_S = 0.08
PARTICLE_DIAMETER_M = 0.0005   # 0.5mm, representative lab-scale sphere size
                                # (matches thermal.py's own default)

_LOSS_MODEL = StateDependentLossModel()


def _run(n_parallel, frequency=1.0):
    sys_ = AMRSystem(GADOLINIUM, mu0H_max=MU0H_T, mass_regenerator=MASS_KG,
                      frequency=frequency, fluid_mdot=MDOT_KG_S,
                      regenerator_effectiveness=0.85, loss_model=_LOSS_MODEL,
                      use_ntu_thermal_model=True, particle_diameter=PARTICLE_DIAMETER_M,
                      hypereg_n_parallel=n_parallel)
    return sys_.run(T_COLD_K, SPAN_K)


def sweep_n_parallel(n_values=(1, 2, 4, 8, 16), frequency=1.0, verbose=True):
    """n_parallel=1 is the conventional series-flow case (identical to
    particle_diameter set with hypereg_n_parallel=None -- see
    AMRSystem._geometry_pumping_power_W()'s docstring). Sweeping n shows
    the pumping-power benefit saturating, since eddy-current and baseline
    ('base_frac') losses (both frequency/field/Qc-driven, not pumping-
    driven) are untouched by n and eventually dominate W_parasitic."""
    rows = []
    for n in n_values:
        res = _run(n, frequency=frequency)
        rows.append((n, res.Qc, res.COP_electrical, res.W_parasitic))
        if verbose:
            print(f"  n_parallel={n:3d}   Qc={res.Qc:8.2f}W   "
                  f"COP_electrical={res.COP_electrical:7.4f}   "
                  f"W_parasitic={res.W_parasitic:7.2f}W")
    return rows


def sweep_frequency_at_fixed_n(frequencies=(0.5, 1.0, 2.0, 4.0, 8.0), n_parallel=4,
                                 verbose=True):
    """Shows how the pumping-power SAVING from Hypereg (vs. conventional
    series flow at the same frequency/mdot) changes with frequency. Since
    this model's pumping power depends on mdot, not frequency, directly
    (see thermal.pumping_power_packed_bed()'s Darcy-flow dP ~ mdot
    dependence), and eddy-current loss DOES scale with f^2, the pumping-
    power saving itself does not grow with frequency in this 0-D model --
    stated plainly rather than implied to be a high-frequency-specific
    benefit beyond what the paper itself claims (the paper's own framing
    is that reduced pumping power is what MAKES higher-frequency, and
    therefore higher-flow, operation practical in the first place, not
    that the per-unit-flow pumping saving itself grows with frequency)."""
    rows = []
    for f in frequencies:
        conv = _run(None, frequency=f)
        hyp = _run(n_parallel, frequency=f)
        saving_pct = 100 * (1 - hyp.W_parasitic / conv.W_parasitic) if conv.W_parasitic > 0 else 0.0
        rows.append((f, conv.COP_electrical, hyp.COP_electrical, saving_pct))
        if verbose:
            print(f"  f={f:5.2f}Hz   COP_conventional={conv.COP_electrical:7.4f}   "
                  f"COP_hypereg(n={n_parallel})={hyp.COP_electrical:7.4f}   "
                  f"W_parasitic saving={saving_pct:5.1f}%")
    return rows


def run_hypereg_analysis(out_path="results/hypereg_analysis.txt"):
    import io, contextlib
    buf = io.StringIO()
    with contextlib.redirect_stdout(buf):
        print("=" * 90)
        print("PHASE 15: Hypereg parallel-hydraulic pumping-power reduction, motivated by")
        print("Klinar et al. (2024), Adv. Energy Mater. 14, 2401739, Figs. 18-21.")
        print("See results/hypereg_findings.md for the full literature findings note.")
        print("=" * 90)
        print(f"\nOperating point: T_cold={T_COLD_K}K, span={SPAN_K}K, mu0H={MU0H_T}T, "
              f"mass={MASS_KG}kg, mdot={MDOT_KG_S}kg/s, d_p={PARTICLE_DIAMETER_M*1000}mm")

        print("\n--- Step 1: sweep n_parallel_subregenerators at f=1.0Hz "
              "(n=1 == conventional series flow) ---")
        n_rows = sweep_n_parallel()

        print("\n--- Step 2: pumping-power saving vs. frequency at a fixed n=4 "
              "(the paper's own illustrative example) ---")
        f_rows = sweep_frequency_at_fixed_n()

        print("\n--- Conclusion ---")
        conv_cop = n_rows[0][2]
        best_cop = max(n_rows, key=lambda r: r[2])
        print(f"At this repo's own representative operating point, splitting the "
              f"regenerator into n={best_cop[0]} parallel sub-beds raises "
              f"COP_electrical from {conv_cop:.3f} (conventional, n=1) to "
              f"{best_cop[2]:.3f} -- a real but modest benefit in THIS model, because "
              "pumping power (W_pump) is only one of three loss channels "
              "(k_eddy*f^2*H^2 + k_pump*mdot^2 + base_frac*Qc) and is not the "
              "dominant one at this lab-scale operating point. The benefit saturates "
              "quickly with n (diminishing returns are visible in Step 1) since the "
              "other two loss channels are untouched by parallelization. This "
              "confirms the qualitative mechanism the paper describes (reduced "
              "pumping power from parallel vs. series hydraulic flow) is real and "
              "representable in this repo's own model, WITHOUT claiming a validated "
              "optimum n or a device-level performance prediction for an as-yet-"
              "unbuilt concept -- see results/hypereg_findings.md's honesty flags.")

    text = buf.getvalue()
    print(text, end="")
    with open(out_path, "w") as fh:
        fh.write(text)
    return text


if __name__ == "__main__":
    run_hypereg_analysis()