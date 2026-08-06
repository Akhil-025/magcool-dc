# COMSOL Step-by-Step Tutorial: 2-D/3-D AMR Regenerator-Bed Model

## How to use this tutorial

This is a click-by-click implementation of `COMSOL_Setup_Guide.md`, written for someone who has never opened COMSOL Multiphysics before. It covers **both** geometry options from the original guide (packed-bed and parallel-plate) and **both** dimensionalities (2-D axisymmetric and full 3-D), and calls out clearly, at every step, where those four tracks differ. Where the source guide didn't specify how to actually operate the software, the step is explicitly labeled **"Software operation (not specified by the guide)"** and a standard COMSOL 6.x workflow is described instead. No equation, assumption, benchmark value, material property, or geometry value from the original guide has been changed.

Menu paths below assume **COMSOL Multiphysics 6.x** (Windows/Linux desktop client, ribbon-style interface). If you're on 5.6 or earlier, the same operations exist under the classic (non-ribbon) menu tree — the model tree structure (Model Wizard → Component → Physics → Geometry → Materials → Mesh → Study) is unchanged across versions.

> ⚠️ **Preserve this warning from the original guide, word for word in spirit:** *This is a setup specification, not a validated result.* Nothing in the original guide was ever built, meshed, or solved in COMSOL. Every geometry, boundary condition, and material property is transcribed from this repo's own 0-D/1-D correlations (`core/thermal.py`, `core/mce_material.py`, `core/first_order_mce.py`, `core/geometry_analysis.py`), so a model built by following this tutorial is *consistent* with the physics this repo already validates — it is **not independent confirmation of anything**. Treat every number this tutorial produces as unvalidated until it passes the degeneracy check in Step 16 (§8 of the original guide).

## 0. Why this model exists (background — no clicking yet)

`core/amr_cycle.py` and `core/thermal.py` are 0-D/lumped models: the regenerator bed is treated as a single effectiveness number (`eps`) fed by an NTU correlation, not as a spatial temperature field. That's adequate for system-level COP comparisons (the `amr_cycle.py` module docstring calls it "good enough for system-level COP comparison... NOT a replacement for a full 2-D/3-D COMSOL regenerator-bed solve"), but it cannot show axial temperature gradients, thermal dispersion, or demagnetization-front propagation within the bed. Those matter if you're optimizing internal bed geometry (channel shape, layer transitions in a Curie-graded bed, wall conduction losses) rather than just system-level mass/frequency/field trade-offs.

This tutorial builds **one AMR unit cell** (one packed-bed or parallel-plate channel, periodic in the cycle) at the same operating point the Python repo already reports numbers for, so you can check the 2-D/3-D result against a real, calibrated benchmark before trusting it for anything new.

## 1. Scope

- **In scope:** a single-material (Gd), single-stage regenerator bed, packed-sphere OR parallel-plate geometry, one full AMR cycle (magnetize → cold-to-hot flow → demagnetize → hot-to-cold flow), reproducing the bed-internal temperature field and resulting Qc/COP at a literature-calibrated operating point.
- **Out of scope:** the 6-layer Curie-graded La(Fe,Si)13Hy bed (`core/cascade.py`), multi-stage cascades, magnet-circuit field mapping (this tutorial assumes a prescribed, spatially uniform `mu0H(t)` waveform, not a solved magnetic field), and structural/mechanical analysis of the regenerator housing. Do not silently extend this geometry to cover those without saying so.

## 2. Prerequisites checklist

Before opening COMSOL, confirm:

| Item | Requirement | Why |
|---|---|---|
| COMSOL license | Includes **Heat Transfer Module** and either **Porous Media and Subsurface Flow Module** or **CFD Module** | The Brinkman Equations interface and Heat Transfer in Porous Media interface both live in these add-on modules, not the base package |
| Lookup CSV file | `gd_dTad_vs_T_1p13T.csv` generated (Step 5 below shows you how) | Both the COMSOL and Fluent guides require the *same* CSV so results stay comparable |
| Python environment | Working `core/mce_material.py` importable (i.e. the repo installed, `numpy` available) | To generate the CSV in Step 5 |
| Values on hand | The Step 3/4/7 tables below (geometry, material, operating point) | You will type every one of these into COMSOL by hand |

**Common mistake:** trying to add the Brinkman Equations interface without the Porous Media module licensed — COMSOL will simply not list it in the Model Wizard. If it's missing, check `File → Help → About COMSOL Multiphysics` module list before troubleshooting anything else.

---

## 3. Creating the model (COMSOL Model Wizard)

*Software operation (not specified by the guide) — standard COMSOL model-creation workflow.*

**Goal:** produce a new, empty COMSOL model file with the correct spatial dimension and the physics interfaces from §2 of the original guide already added.

**Why this step is needed:** COMSOL organizes every simulation as a "Model" containing one or more "Components," each with a fixed spatial dimension and a list of physics interfaces. You must choose the dimension (2-D axisymmetric vs 3-D) and physics up front — dimension can be changed later only by deleting and re-adding the component, so get it right here.

### 3.1 Launch and choose dimension

1. Open **COMSOL Multiphysics** (Windows: Start Menu → COMSOL 6.x; Linux: run `comsol` from the install `bin` directory).
2. On the startup screen, click **Model Wizard**.
3. A dimension-selection screen appears with icons: 0D, 1D, 1D Axisymmetric, 2D, 2D Axisymmetric, 3D.
   - **2-D axisymmetric track:** click **2D Axisymmetric**. Use this if you're modeling a cylindrical packed-bed segment as described in §3.1 of the original guide — it lets you resolve radial + axial gradients while solving only a 2-D `(r, z)` slice, which is far cheaper than full 3-D and is adequate as long as the bed and flow are genuinely axisymmetric (uniform packing, no circumferential asymmetry).
   - **3-D track:** click **3D**. Use this for a true representative-volume packed bed, or for the parallel-plate geometry if you want to resolve the finite plate width/depth rather than assuming translational symmetry.
   - **Parallel-plate 2-D-as-cross-section note:** the parallel-plate geometry is not naturally axisymmetric (it's a stack of flat plates, not a body of revolution). If you choose 2-D for the parallel-plate track, use plain **2D** (not 2D Axisymmetric) and treat it as a cross-section through one flow channel with symmetry/periodic boundaries on the top/bottom of the modeled slice — see Step 4.4.

**Expected result:** the Model Builder tree on the left now shows a single **Component 1** node under a **Global Definitions** node, with the dimension you picked locked in.

**How to verify:** hover over the Component node — its tooltip or the ribbon's Component settings should confirm the chosen dimension (e.g., "2D Axisymmetric").

**Common mistake:** picking plain **2D** when you meant **2D Axisymmetric** (or vice versa) for the packed-bed track — the geometry, boundary conditions, and even the meaning of "radius" differ between the two, and COMSOL will not warn you; it will just solve the wrong problem.

### 3.2 Add physics interfaces

**Goal:** add the two physics interfaces the original guide's §2 specifies — a porous-flow interface and a heat-transfer interface — so COMSOL knows what equations to solve in which domain.

**Why this step is needed:** COMSOL doesn't automatically know you want porous flow coupled to heat transfer with a custom volumetric source; you must explicitly add both interfaces (and, later, a Multiphysics coupling node) to the component.

4. In the **Add Physics** window (opens automatically after dimension selection), expand **Fluid Flow → Porous Media and Subsurface Flow**.
5. Select **Brinkman Equations (br)**. Click the blue **+ (Add)** arrow to add it to Component 1.
   - *Why Brinkman and not "Free and Porous Media Flow":* the original guide (§2.1) recommends Brinkman Equations with a Darcy-Forchheimer drag term for a bed-averaged, 1-unit-cell model. Only switch to **Free and Porous Media Flow** if you specifically need to resolve the explicit void space between individual spheres in a small representative volume — that's a materially different (and much more expensive) model than what this tutorial builds.
6. Expand **Heat Transfer → Heat Transfer in Porous Media (ht)**. Select it and click **+ (Add)**.
   - This interface will host both the fluid-side convective heat transport and the solid-side MCE heat source term (Step 9).
7. Click the blue arrow **Study** at the top of the Add Physics window to proceed (or click it later from the ribbon — see Step 12).

**Screens that should appear:** the Model Builder tree now shows, under Component 1: **Definitions**, **Geometry 1**, **Materials**, **Brinkman Equations (br)**, **Heat Transfer in Porous Media (ht)**, **Multiphysics**, **Mesh 1**.

**Expected result before moving on:** both physics interface nodes are present with no red error icons.

**How to verify:** right-click **Brinkman Equations (br)** → it should list default sub-nodes **Fluid and Matrix Properties 1**, **Initial Values 1**. Same check for **Heat Transfer in Porous Media (ht)**: default sub-nodes **Porous Medium 1**, **Initial Values 1**, **Thermal Insulation 1**.

**Common mistake:** adding **Laminar Flow** instead of **Brinkman Equations** — Laminar Flow solves the full Navier-Stokes equations in open fluid domains and has no porous-drag formulation; it is the wrong interface for this bed-averaged model and will not accept the friction-factor correlation from Step 8.

---

## 4. Geometry creation

The original guide (§3) gives two geometry tracks. Build **whichever track(s) match your research question**; both are documented fully below because you're building both per your instructions.

### 4.1 Packed-bed unit cell — parameter values to enter

| Quantity | Value | Where it's used | Source in original guide |
|---|---|---|---|
| Particle diameter `d_p` | 0.5 mm default (sweep range 0.05–2 mm) | Sphere/cylinder radius in geometry | §3.1 — Phase 7 geometry sweep found COP-optimal packed-bed diameter at 0.5 mm for the 291K/10K-span/2kg/1Hz/1.5T operating point; **not necessarily optimal** at the 1.13T/10.2K-span/1.7kg/0.75Hz benchmark point used for the Step 14 degeneracy check — re-sweep if geometry itself is the research question |
| Porosity `ε` | 0.365 | Domain-averaged void fraction in Brinkman/porous settings | §3.1 — `core/thermal.py`'s `regenerator_effectiveness()` default |
| Bed cross-section | 0.002 m² (~5×4 cm face) | Overall bed footprint | §3.1 — `core/thermal.py` default, representative of lab-scale devices in `data/amr_experimental_benchmarks.csv` |
| Solid density `ρ_Gd` | 7900 kg/m³ | Material property, not geometry — listed here because the original guide lists it in the same table | §3.1 — `core/thermal.py`'s `RHO_GD` |

**Bed length:** the original guide does not give an explicit axial bed length in §3; derive it from mass, density, porosity and cross-section: `L_bed = m_bed / (ρ_Gd × (1-ε) × A_cross)`. Using the Step 13 operating-point mass (1.7 kg) and this section's values: `L_bed = 1.7 / (7900 × (1-0.365) × 0.002) = 0.170 m` (≈17 cm). *Software operation (not specified by the guide): this derivation follows directly from the guide's own stated quantities, but the guide itself does not spell out the arithmetic — shown here for completeness.*

#### 4.1a Building the packed-bed geometry — 2-D axisymmetric track

**Goal:** create a 2-D `(r, z)` rectangle representing the bed as a single homogenized porous domain (not individual resolved spheres — this is a bed-averaged Brinkman model, so you draw one bulk domain, not a sphere packing).

**Why this step is needed:** the Brinkman equations interface treats the bed as a continuum with an effective porosity and drag law; you don't need — and should not build — a literal packed-sphere CAD geometry for this level of fidelity (see §2.1 note: "a bed-averaged 1-unit-cell model" is the target fidelity, not the fully resolved void-space model).

1. In the Model Builder, right-click **Geometry 1 → Rectangle**.
2. In the Rectangle settings panel:
   - **Width**: enter the bed's radius (pick a representative unit-cell radius consistent with the 0.002 m² cross-section, e.g. for a circular cross-section `r = sqrt(A/π) = sqrt(0.002/π) ≈ 0.0252 m`).
   - **Height**: `0.170 m` (bed length, computed above).
   - **Position → r**: `0`. **Position → z**: `0`.
3. Click **Build Selected** (or the **Build All Objects** icon — the checkered cube in the ribbon).

**Screens that should appear:** a rectangle appears in the Graphics window; the geometry sequence shows a green checkmark on **Rectangle 1**.

**Expected result before moving on:** one 2-D domain, axisymmetric about `r = 0`, height 0.170 m, radius ≈0.0252 m.

**How to verify:** click **Rectangle 1**, check the **Measurements** or hover-readout in the Graphics window matches the entered dimensions; confirm the domain sits entirely at `r ≥ 0` (COMSOL will error on Build if any geometry crosses `r < 0` in an axisymmetric component).

**Common mistake:** forgetting axisymmetric geometry is drawn in the **half-plane** `(r ≥ 0, z)` — do not mirror the rectangle across `r = 0`; COMSOL handles the revolution internally.

#### 4.1b Building the packed-bed geometry — 3-D track

**Goal:** create a 3-D cylindrical (or rectangular-prism) representative volume for the same bed.

1. Right-click **Geometry 1 → Cylinder**.
2. In the Cylinder settings:
   - **Radius**: `0.0252 m` (from the 0.002 m² cross-section, as above — or use a rectangular prism with **Block** instead of **Cylinder** if your lab-scale housing is rectangular; the original guide's "~5×4 cm face" note under §3.1 suggests a rectangular cross-section is equally valid — if so, use **Block** with **Width = 0.05 m, Depth = 0.04 m**).
   - **Height**: `0.170 m`.
   - **Axis**: default `z`-axis, **Position**: `(0,0,0)`.
3. Click **Build Selected**.

**Expected result:** a single 3-D solid domain, either a cylinder (r=0.0252 m, h=0.170 m) or a 0.05×0.04×0.170 m rectangular prism.

**How to verify:** rotate the Graphics view (drag with left mouse button) to confirm it's a genuine 3-D solid, not a surface; check **Geometry 1 → Cylinder 1/Block 1** has a green checkmark.

**Common mistake:** using a **Work Plane + Extrude** sequence when a direct **Cylinder**/**Block** primitive would do — not wrong, just unnecessarily complex for this simple unit-cell shape.

### 4.2 Parallel-plate alternative — parameter values to enter

| Quantity | Value | Where it's used | Source in original guide |
|---|---|---|---|
| Plate spacing (fluid gap) | 0.1–0.2 mm typical; Phase-7 COP optimum at **0.1 mm** at the 291K/10K operating point | Channel gap dimension in geometry | §3.2 — `core/geometry_analysis.py` |
| Plate thickness | 0.5 mm default | Solid plate dimension | §3.2 — `core/thermal.py`'s `regenerator_effectiveness_parallel_plate()` |
| Hydraulic diameter | `d_h = 2 × spacing` | Used later in the friction-factor correlation, Step 8 | §3.2 — Eq. 7 slot limit |

#### 4.2a Building the parallel-plate geometry — 2-D track

**Goal:** draw a 2-D cross-section (in the flow-times-gap plane) showing one fluid channel bounded by two solid plates, representing a unit cell of the plate stack.

1. This track uses plain **2D** dimension (not axisymmetric — see Step 3.1's parallel-plate note). If you built the component as 2D Axisymmetric for the packed-bed track, you need a **separate Component** for this track: right-click the root **Global Definitions** area → **Add Component → 2D**, then repeat Step 3.2's physics-interface additions in the new component.
2. Right-click **Geometry 1 → Rectangle** (fluid domain):
   - **Width** (flow direction, length of bed): `0.170 m` (reuse the same derived bed length as §4.1 for consistency, unless the parallel-plate device has its own documented length — the original guide does not give one).
   - **Height** (gap): `0.0001 m` (0.1 mm, the Phase-7 optimum) — or `0.0002 m` if using the "typical" upper end of the 0.1–0.2 mm range.
   - **Position**: `(0, 0)`.
3. Add a second **Rectangle** for the solid plate immediately above the fluid domain:
   - **Width**: `0.170 m`. **Height**: `0.0005 m` (0.5 mm plate thickness).
   - **Position**: `(0, 0.0001)` (stacked directly on top of the fluid domain's height).
4. Click **Build All Objects**.

**Screens that should appear:** two adjacent rectangles — a thin fluid layer (0.1 mm tall) directly below a thicker solid layer (0.5 mm tall), both 0.170 m long.

**Expected result before moving on:** two touching, non-overlapping 2-D domains sharing one boundary edge (the plate/fluid interface).

**How to verify:** click on the shared edge in the Graphics window — it should highlight as a single boundary shared by both domains (check via **Geometry → Measure** or by clicking **Form Union** finalization and confirming no duplicate/overlapping boundaries in the **Errors and Warnings** log).

**Common mistake:** leaving a gap or overlap between the fluid and solid rectangles — COMSOL's default **Form Union** finalization will either silently merge them incorrectly or throw a geometry error; always verify with **View → Wireframe Rendering** that the two domains share a single, contiguous edge.

#### 4.2b Building the parallel-plate geometry — 3-D track

**Goal:** extrude the 2-D cross-section above into a finite-width, finite-depth 3-D channel stack (useful if you want to check for lateral/depth-direction effects the 2-D idealization can't show).

1. Build the 2-D cross-section as in §4.2a within a **Work Plane** (right-click **Geometry 1 → Work Plane**, draw the two rectangles on the work plane).
2. Right-click **Geometry 1 → Extrude**, select the Work Plane, and set **Distance**: `0.04 m` (using the "~4 cm" dimension implied by the original guide's 0.002 m² packed-bed cross-section note, as a reasonable representative depth — the guide does not give a specific plate width/depth, so this value is a *software operation (not specified by the guide)* choice made for dimensional consistency with the packed-bed cross-section; document this choice explicitly in any report).
3. Click **Build All Objects**.

**Expected result:** a 3-D block-shaped fluid channel (0.170 m × 0.0001 m × 0.04 m) sandwiched against a 3-D solid plate slab of the same footprint, 0.5 mm thick.

**How to verify:** rotate the 3-D view; confirm two distinct, correctly proportioned solid domains (the fluid layer will look like an extremely thin sliver next to the plate — this is expected given the 0.1 mm gap vs 0.170 m length; zoom in near one end to confirm the gap geometry is not degenerate/zero-thickness due to a units error).

**Common mistake:** entering plate spacing in mm without converting to meters (COMSOL's default length unit is meters) — a `0.1` typed where `0.0001` was meant produces a channel 1000× too large and silently wrong results. Always check the **Length unit** dropdown in **Geometry 1** settings (top of the Geometry node) reads **m**, and always write the full value in meters as shown in the tables above.


---

## 5. Generating and importing the MCE source-term lookup table

**Goal:** produce a `(T, ΔT_ad, C_total)` table from this repo's own material model and import it into COMSOL as an Interpolation function, so the 2-D model's heat source is provably the same physics as the 0-D model it's being checked against.

**Why this step is needed:** the original guide (§2.2–2.3) is explicit that you must **not** re-derive the mean-field/Brillouin MCE model from scratch inside COMSOL — you export a lookup table from `core/mce_material.py` and import it directly, so any known bias in that model (see the warning below) is inherited faithfully rather than silently re-implemented differently.

### 5.1 Generate the CSV (outside COMSOL, in Python)

Run exactly this script (reproduced unchanged from the original guide's §2.3):

```python
import numpy as np
from core.mce_material import GADOLINIUM

mu0H = 1.13  # T -- matches the DTU_Eriksen_rotary_Gd_2015 benchmark, see §8
Ts = np.linspace(270, 310, 401)
dTad = GADOLINIUM.delta_T_adiabatic(Ts, mu0H / (4 * np.pi * 1e-7))
C = GADOLINIUM.total_heat_capacity(Ts, mu0H / (4 * np.pi * 1e-7))
np.savetxt("gd_dTad_vs_T_1p13T.csv",
           np.column_stack([Ts, dTad, C]), delimiter=",",
           header="T_K,dTad_K,C_total_J_per_kgK", comments="")
```

**Expected result:** a file `gd_dTad_vs_T_1p13T.csv` with 401 rows and 3 columns (`T_K`, `dTad_K`, `C_total_J_per_kgK`), spanning T = 270–310 K.

**How to verify:** open the CSV and confirm the `dTad_K` column peaks somewhere near T = 294 K (the Gd Curie temperature — see Step 6's material table) — a flat or zero column means the import path or field value is wrong before you've even opened COMSOL.

> ⚠️ **Preserve this warning:** this table inherits the mean-field model's own documented **+29% to +49% overprediction** of ΔT_ad relative to Dan'kov et al. (1998) at 1–2 T, improving to **−7.5%** at 5 T (`core/validation.py`, `results/pipeline.log` step 1). **A 2-D model built on this table inherits that same known bias; it does not fix it.**

### 5.2 Import as an Interpolation function

*Software operation (not specified by the guide) — standard COMSOL import workflow.*

3. In the Model Builder, right-click **Global Definitions → Functions → Interpolation**.
4. In the Interpolation settings panel:
   - **Data source**: select **File**.
   - **Filename**: browse to `gd_dTad_vs_T_1p13T.csv`. Click **Import**.
   - **Number of arguments**: `1`.
   - Under **Function names**, you should see two entries auto-populated for the two data columns beyond the first (T is the argument column): rename them (or confirm names) as `dTad_lookup` and `Ctotal_lookup` for clarity matching the ANSYS guide's UDF variable names.
   - **Interpolation**: set to **Linear** (piecewise-linear interpolation — matches the ANSYS guide's §2.3 requirement that both tools use the same interpolation scheme so results stay comparable).
   - **Extrapolation**: set to **Constant, extrapolate using nearest function value** to avoid nonphysical extrapolated source terms if the solver's temperature briefly steps outside 270–310 K during transients.
5. Click **Plot** (in the ribbon, if available for the Interpolation node) to render `dTad_lookup(T)` and confirm visually it peaks near 294 K.

**Screens that should appear:** the Interpolation node shows a green checkmark; the function plot (if generated) shows a peak-shaped curve in T, consistent with a λ-anomaly-like feature near the Curie point.

**Expected result before moving on:** two named, callable global functions `dTad_lookup(T)` and `Ctotal_lookup(T)`, usable anywhere else in the model (e.g., `dTad_lookup(T)` inside a heat-source expression).

**How to verify:** in any later expression field (e.g., a test point in **Global Definitions → Parameters**), type `dTad_lookup(294[K])` and evaluate — it should return a nonzero value with units of K.

**Common mistake:** forgetting to set the interpolation argument's **unit** — if COMSOL treats the imported `T_K` column as dimensionless instead of Kelvin, later expressions like `dTad_lookup(T)` (where `T` carries units of K from the physics interface) will throw a units-mismatch error or silently evaluate at the wrong temperature. Check the **Argument** tab of the Interpolation node and explicitly set the unit to `K`.

### 5.3 (Conditional) Import the full temperature-dependent heat capacity table

The original guide's §4 material table flags: *"If the model geometry spans temperatures near Tc=294K, use `core/mce_material.py`'s full `total_heat_capacity(T)` (lattice + magnetic λ-anomaly) instead of this flat constant [236 J/(kg·K)]."* Since the Step 13 benchmark operating point (T_cold=289 K, span=10.2 K → up to ~299 K) straddles Tc=294 K, **use the full table** — this is exactly the `Ctotal_lookup` function already imported in §5.2, so no additional action is needed here beyond referencing `Ctotal_lookup(T)` instead of a flat 236 J/(kg·K) constant wherever heat capacity appears (Step 9's source term, and Step 6's material Cp field).


---

## 6. Material creation

**Goal:** define two materials — solid Gd and water — with the exact property values from the original guide's §4, and assign them to the correct geometric domains.

**Why this step is needed:** COMSOL's physics interfaces read density, heat capacity, thermal conductivity, and viscosity from Material nodes, not from the geometry itself; every property below must be entered exactly, or the porous-flow and heat-transfer equations will use COMSOL's built-in defaults (usually air/water at 20°C) instead of this repo's calibrated values.

### 6.1 Material property values (reproduced unchanged from the original guide's §4)

| Property | Value | Source |
|---|---|---|
| Gd solid density | 7900 kg/m³ | `core/thermal.py::RHO_GD` |
| Gd solid specific heat (lattice, near-room-T representative) | 236 J/(kg·K) | `core/thermal.py::CP_SOLID_GD` — **use `Ctotal_lookup(T)` from §5.3 instead of this flat value**, because the model spans temperatures near Tc=294K |
| Gd Curie temperature Tc | 294.0 K | `core/mce_material.py::GADOLINIUM` |
| Gd Debye temperature θ_D | 169.0 K | same |
| Gd molar mass | 157.25 g/mol | same |
| Water density | 997 kg/m³ | `core/thermal.py::water_properties()` |
| Water specific heat | 4186 J/(kg·K) | same |
| Water dynamic viscosity | 8.9×10⁻⁴ Pa·s | same |
| Water thermal conductivity | 0.606 W/(m·K) | same |

> ⚠️ **Preserve this note:** these water properties are constant-value simplifications ("adequate for this 0-D estimate; a full model would use IAPWS correlations" — `water_properties()` docstring). A 2-D/3-D model with real temperature gradients is exactly the case where this stops being adequate — if you replace them with COMSOL's built-in IAPWS-97 water material (Step 6.4), **note in any writeup that this is a genuine improvement over the 0-D model, not a discrepancy to explain away.**

### 6.2 Create the Gd solid material

1. Right-click **Component 1 → Materials → Blank Material**.
2. Rename it (double-click, or F2): **Gadolinium (Gd) — solid**.
3. In the **Material Contents** table, add the following properties (click the **+** row or type directly into the **Value** column for each existing default row; add new rows via the **Material property** browse tree on the right, e.g. under **Basic Properties**):
   - **Density (rho)**: `7900[kg/m^3]`
   - **Heat capacity at constant pressure (Cp)**: `Ctotal_lookup(T)` (an expression referencing the imported function from §5.3 — this makes Cp temperature-dependent rather than a flat number)
   - **Thermal conductivity (k)**: the original guide does not give an explicit Gd thermal-conductivity value in its §4 table. *Software operation (not specified by the guide):* COMSOL's Heat Transfer in Porous Media interface requires a conductivity value for the solid matrix; if the source repo doesn't specify one, use a literature value for polycrystalline Gd near room temperature (~10.5 W/(m·K), per standard elemental-Gd references) and flag this explicitly as **not sourced from the original guide** in any report.
4. Under **Geometric Entity Selection**, set **Selection** to the solid domain(s) — for the packed-bed track this is the single bed domain (since the Brinkman porous formulation treats the bed as one homogenized domain with both fluid and solid character — see §6.5 below on how porosity separates the two within one domain); for the parallel-plate track, select specifically the **plate rectangle/block**, not the fluid channel.

**Expected result:** the Materials tree shows **Gadolinium (Gd) — solid** with a filled-in (not default/blank) icon, applied to the correct domain(s).

**How to verify:** click the material node — the **Graphics** window should highlight the assigned domain(s) in blue; hover over any property row to confirm the expression evaluates without a red error underline.

### 6.3 Create the water material

5. Right-click **Materials → Blank Material** again. Rename: **Water (constant properties)**.
6. Enter:
   - **Density (rho)**: `997[kg/m^3]`
   - **Heat capacity at constant pressure (Cp)**: `4186[J/(kg*K)]`
   - **Dynamic viscosity (mu)**: `8.9e-4[Pa*s]`
   - **Thermal conductivity (k)**: `0.606[W/(m*K)]`
7. Assign this material to the fluid domain: for the packed-bed track, the same bed domain as the solid material (Brinkman treats both fluid and solid material data as inputs to one domain via the porosity-weighted formulation — see §6.5); for the parallel-plate track, the **fluid channel rectangle/block only**.

**Common mistake:** entering `8.9e-4` without the `[Pa*s]` unit tag — COMSOL will interpret a bare number using its current default unit system, which may not be Pa·s, silently scaling the viscosity by orders of magnitude. **Always type explicit unit tags in square brackets** for every numeric material property in COMSOL.

### 6.4 (Conditional) Upgrade water to COMSOL's built-in IAPWS-97 property set

*Software operation (not specified by the guide) — optional COMSOL built-in material.*

8. If you decide the temperature-dependent water properties are warranted (per the §6.1 warning), instead of a Blank Material: right-click **Materials → Add Material from Library**, search **"Water, liquid"** in the built-in Material Library, and add it. This replaces the four constant values above with COMSOL's temperature-dependent IAPWS-97 correlations automatically.
9. Document this substitution explicitly wherever you report results — per the original guide, this is "a genuine improvement over the 0-D model, not a discrepancy to explain away."

### 6.5 How porosity connects the two materials in one Brinkman domain (packed-bed track)

**Goal:** tell the Brinkman Equations and Heat Transfer in Porous Media interfaces that the bed domain is `ε = 0.365` fluid by volume and `(1-ε) = 0.635` solid by volume, not two separately-meshed sub-regions.

10. Click **Brinkman Equations (br) → Fluid and Matrix Properties 1**.
11. Set **Fluid material**: **Water (constant properties)** (or the IAPWS material from §6.4).
12. Set **Porosity model**: choose the option that lets you enter a fixed porosity value; enter **Porosity (εp)**: `0.365`.
13. Click **Heat Transfer in Porous Media (ht) → Porous Medium 1**.
14. Under **Porous Matrix**, set **Solid material**: **Gadolinium (Gd) — solid**.
15. Under **Fluid**, set **Fluid material**: **Water (constant properties)**.
16. Confirm the **Porosity** field here matches the value entered in Step 12 (`0.365`) — COMSOL does not automatically synchronize porosity between the Brinkman and Heat Transfer interfaces; you must enter it in both, and a mismatch here is a silent, hard-to-spot bug.

**Expected result before moving on:** both physics interfaces reference the same two materials and the same porosity value for the bed domain.

**How to verify:** open both **Fluid and Matrix Properties 1** and **Porous Medium 1** side by side (or check each in turn) and confirm porosity = 0.365 in both, and that Solid material = Gd, Fluid material = Water, in both.

**Common mistake:** entering porosity in only one of the two interfaces — COMSOL will not error, it will just silently use its default porosity (often 1, i.e., treating the domain as pure fluid with no solid heat capacity) in whichever interface you forgot, producing a bed that stores far too little thermal mass.

*(For the parallel-plate track, porosity does not apply the same way — the fluid and solid are geometrically separate domains, ε=1 in the fluid domain and the solid domain is pure Gd with no porous formulation needed; skip the porosity entries above for that track and instead just assign the materials to their respective, separate domains as described in §6.2–6.3.)*


---

## 7. Porous-media drag law (packed-bed and parallel-plate friction factors)

**Goal:** make the Brinkman equations' drag term reproduce the exact same `f(Re)` correlation `core/thermal.py` uses, rather than falling back to COMSOL's generic default.

**Why this step is needed:** the original guide (§2.1) is explicit: *"reproduce the same f(Re) as a User Defined Function or Analytic function in COMSOL rather than falling back to a generic Ergun-equation default, so the 2-D pressure drop is directly comparable to the 0-D number."*

### 7.1 Correlations to implement (reproduced unchanged from the original guide)

| Geometry | Correlation | Validity | Source |
|---|---|---|---|
| Packed bed | `f = 23.462 · Re^-0.6716` | `10 < Re < 5e5` | Tušek, Kitanovski, Poredoš, *Int. J. Refrig.* 36 (2013) 1456-1464, Eq. 5; `Re = ρ_f u_s d_p / μ_f` (superficial velocity, particle diameter) |
| Parallel plate | `f = 24/Re` (laminar) | `Re < 2300` | same paper, Eq. 6; `d_h = 2·(plate spacing)` |

Hydraulic diameter for the packed bed: `d_h = 4·V_bed·ε/A_total` (Eq. 7 of the same paper), equivalent to `d_h = (2/3)·d_p·ε/(1-ε)`.

### 7.2 Create the Analytic function

*Software operation (not specified by the guide) — standard COMSOL Analytic function workflow, implementing the correlation the guide specifies.*

1. Right-click **Global Definitions → Functions → Analytic**.
2. Rename: **f_friction_packedbed** (packed-bed track) or **f_friction_plate** (parallel-plate track) — build one Analytic function per geometry track you're implementing.
3. **Expression**:
   - Packed bed: `23.462*Re^(-0.6716)`
   - Parallel plate: `24/Re`
4. **Arguments**: `Re`. **Argument units**: dimensionless (`1`).
5. Under **Plot Parameters**, set the Re range to match validity: packed bed `10` to `5e5`; parallel plate `1` (or a small positive number to avoid divide-by-zero) to `2300`.
6. Click **Plot** to confirm the curve is monotonically decreasing (both correlations are power-law decays in Re) and passes through a sanity-check point, e.g. packed bed at Re=100: `f = 23.462 × 100^-0.6716 ≈ 1.02`.

**Expected result before moving on:** a callable function `f_friction_packedbed(Re)` (or `f_friction_plate(Re)`) usable in expressions elsewhere in the model.

**How to verify:** in **Global Definitions → Parameters**, add a test parameter `Re_test = 100` and an expression `f_friction_packedbed(Re_test)`; confirm it evaluates to ≈1.02.

### 7.3 Wire the friction factor into the Brinkman drag term

7. Click **Brinkman Equations (br) → Fluid and Matrix Properties 1**.
8. Under **Porous Matrix Properties**, find the permeability/drag model dropdown. Select **User defined** (not "Kozeny-Carman" or "Ergun" — those are COMSOL's generic defaults the original guide says not to fall back on).
9. For a Darcy-Forchheimer formulation, COMSOL's User-Defined drag option typically asks for a **permeability κ** and, for the Forchheimer term, a **drag coefficient**. Relate these to the friction factor via the standard Darcy-Weisbach identity `ΔP/L = f · (ρ_f u_s²)/(2 d_h) · (1/ε)` (packed-bed) or the plate-flow equivalent, and enter:
   - **Permeability (κ)**: derived from the correlation's Reynolds-independent (Darcy/viscous) term.
   - **Forchheimer/inertial coefficient**: derived from the correlation's Reynolds-dependent (inertial) term, since the packed-bed correlation's exponent (`-0.6716`, not `-1`) implies a mixed viscous-inertial regime over the stated Re range, not pure Darcy flow.
   - *Software operation (not specified by the guide):* the original guide specifies the target `f(Re)` correlation but does not give the exact COMSOL UI field mapping from `f(Re)` to COMSOL's built-in `κ`/Forchheimer parameterization — fit these two parameters against `f_friction_packedbed(Re)` (or `_plate`) over the Reynolds range the Step 13 operating point actually spans, using COMSOL's **Optimization Module** parameter-estimation tool if available, or a manual least-squares fit performed outside COMSOL and then entered as fixed numbers. **Do not accept COMSOL's Ergun-equation default without checking it against this fit.**

**Expected result:** the Brinkman interface's pressure-drop calculation, when evaluated at a representative Re from the operating point, matches `f_friction_packedbed`/`f_friction_plate` to within the fit tolerance.

**How to verify:** after meshing and a first solve (Steps 11–12), use **Derived Values → Line Average** to extract `ΔP` across the bed at the operating-point superficial velocity, and hand-check it against `ΔP = f(Re)·(ρu_s²)/(2d_h)·L/ε` using the same `f` from §7.1.

**Common mistake:** leaving the drag model on COMSOL's default (Ergun equation) — this is a *different* correlation with different constants than Tušek et al. Eq. 5/6, and the resulting pressure drop (and hence flow-field/residence-time distribution) will not match the 0-D model this repo already validates, breaking the whole point of the degeneracy check in Step 16.


---

## 8. Implementing the MCE heat source term

**Goal:** add a volumetric heat source to the solid (Gd) domain that switches on during magnetize/demagnetize phases and reproduces exactly:

```
q_MCE(T) = ρ_solid · C_total(T) · ΔT_ad(T, μ0H_max) / t_mag
```

where `t_mag` is the magnetization ramp duration and `ΔT_ad(T, μ0H_max)`/`C_total(T)` come from the §5 lookup table (not re-derived inside COMSOL).

**Why this step is needed:** the regenerator solid is not a passive conductor — during magnetization and demagnetization it releases/absorbs heat via the magnetocaloric effect. This is the physics that makes AMR cooling work, and it must be modeled as an explicit source term because COMSOL's built-in heat-transfer physics has no native "magnetocaloric material" option.

### 8.1 Compute `t_mag`

`t_mag` is the magnetization ramp duration, derived from the AMR cycle frequency and blow fraction (Step 10 covers this fully). With the Step 13 operating point (frequency = 0.75 Hz, so period = 1/0.75 = 1.333 s) and the default blow fraction `BLOW_FRACTION_MASCHE = 0.5` (symmetric cycle — each of the four phases gets an equal quarter of the period unless a specific device's own blow fraction is documented differently): `t_mag = 0.25 × 1.333 s = 0.333 s`.

*Software operation (not specified by the guide): the guide states the four-phase structure and default blow fraction but does not spell out this exact division-by-four arithmetic for computing t_mag — shown here for completeness, following directly from the guide's own definitions.*

### 8.2 Add a global parameter for the cycle state and t_mag

1. Go to **Global Definitions → Parameters 1**. Add:
   - `t_mag` = `0.333[s]` — description: "Magnetization/demagnetization ramp duration"
   - `rho_solid` = `7900[kg/m^3]` — description: "Gd solid density (RHO_GD)"
   - `mu0H_max` = `1.13[T]` — description: "Peak applied field, DTU_Eriksen benchmark"

### 8.3 Add the Heat Source node

2. Right-click **Heat Transfer in Porous Media (ht) → Heat Source**. (If not directly visible, right-click the **ht** interface node itself and choose **Heat Source** from the physics feature list.)
3. In **Domain Selection**, choose the solid Gd domain(s) — for the packed-bed track, this is the whole bed domain (since porous formulation carries both fluid and solid material behavior in one domain — the heat source physically belongs to the solid fraction, which COMSOL's Heat Transfer in Porous Media interface accounts for automatically via the porosity-weighted energy balance); for the parallel-plate track, select the **plate solid domain only**, not the fluid channel.
4. In the **Heat Source** field, enter the expression:

```
rho_solid * Ctotal_lookup(T) * dTad_lookup(T) / t_mag * cycle_sign
```

where `cycle_sign` is a piecewise indicator variable defined in Step 10 (Events interface) that equals `+1` during magnetization, `-1` during demagnetization, and `0` during the two flow phases — this reproduces the original guide's sign convention ("MCE source term active with the sign flipped" during demagnetize) in one expression rather than four separate source nodes.

**Expected result before moving on:** a Heat Source node in the solid domain(s), referencing `Ctotal_lookup`, `dTad_lookup`, `t_mag`, and `cycle_sign` — all previously defined.

**How to verify:** click into the expression field — COMSOL underlines undefined variables in red; confirm no red underlines. Then, with `cycle_sign` temporarily hardcoded to `1` for a quick check (before building the full Events logic in Step 10), run **Physics → Check Consistency** or a trial **Study → Get Initial Value** to confirm the source term evaluates to a finite, physically reasonable number (order of magnitude: with `Ctotal_lookup(294K)` on the order of a few hundred J/(kg·K), `dTad_lookup(294K, 1.13T)` on the order of a few K, and `t_mag = 0.333 s`, expect `q_MCE` on the order of 10⁶–10⁷ W/m³ — a source term many orders of magnitude off this suggests a units error, most likely a missing `[K]`/`[s]` tag).

**Common mistake:** applying the heat source to the fluid domain instead of (or in addition to) the solid domain — the fluid (water) has no magnetocaloric effect; only the Gd solid releases/absorbs heat this way. Double-check **Domain Selection** on the Heat Source node against the material assignments from Step 6.


---

## 9. Boundary conditions

**Goal:** set the cold-end, hot-end, wall, and inlet boundary conditions exactly as specified in the original guide's §6.

**Why this step is needed:** boundary conditions determine how heat and mass enter/leave the unit cell each cycle phase — get these wrong and the bed never reaches a meaningful periodic state no matter how correct the interior physics is.

### 9.1 Boundary condition values (reproduced unchanged from the original guide's §6)

| Boundary | Condition | Active during | Value |
|---|---|---|---|
| Cold end | Fixed-temperature or fixed-heat-flux reservoir | Hot-to-cold flow (phase 4), when fluid draws heat from this reservoir into the bed | `T_cold = 289 K` |
| Hot end | Fixed-temperature reservoir | Cold-to-hot flow (phase 2) | `T_hot = T_cold + span = 289 + 10.2 = 299.2 K` |
| Bed walls | Adiabatic (no radial heat loss) | Always | matches the 0-D model, which has no wall-loss term at all |
| Fluid inlet | Mass flow rate `mdot` | Only during the two flow phases (2 and 4), OFF during magnetize/demagnetize | `mdot = 0.084666 kg/s` |

> ⚠️ **Preserve this note:** adding a wall-conduction loss term would be a genuine addition beyond what `amr_cycle.py`/`loss_model.py` currently capture — **if added, say so explicitly** rather than silently changing what's being compared.

### 9.2 Set the cold-end and hot-end temperature boundaries

1. Click **Heat Transfer in Porous Media (ht) → Boundaries**. Identify the two end-faces of your bed geometry (bottom face at z=0, top face at z=L_bed, for the packed-bed track; the two flow-direction end-faces for the parallel-plate track).
2. Right-click **ht → Temperature**. Rename: **Cold reservoir BC**.
   - **Boundary Selection**: the cold-end face.
   - **Temperature**: `T_cold * flag_flow4` — using a piecewise flag from Step 10 so this BC is only "active" (imposing 289 K) during phase 4; during other phases, switch this node's activation using the Events interface (§10.4) rather than trying to force it through the expression alone, since a Temperature BC that's always technically present will always clamp the boundary even if its value expression evaluates strangely — the cleaner approach is an Events-driven **enable/disable** of the whole BC node (see §10.4).
3. Right-click **ht → Temperature** again. Rename: **Hot reservoir BC**.
   - **Boundary Selection**: the hot-end face.
   - **Temperature**: `T_hot` = `299.2[K]` (or the parameter expression `T_cold + span` with `span = 10.2[K]` defined in Global Parameters — preferred, since it keeps the dependency explicit and matches how `core/amr_cycle.py` computes `T_hot`).

### 9.3 Set adiabatic walls

4. Confirm **Heat Transfer in Porous Media (ht) → Thermal Insulation 1** (COMSOL's default boundary condition for any boundary not otherwise assigned) covers the bed's lateral/radial walls — for the 2-D axisymmetric packed-bed track this is the outer radial boundary (`r = r_max`); for 3-D, the cylindrical or prismatic side walls. No additional node is usually needed here since Thermal Insulation (zero heat flux) is COMSOL's default — but explicitly verify it, don't just assume.

**How to verify:** click **Thermal Insulation 1**, check its **Boundary Selection** highlights exactly the lateral wall boundaries (blue in the Graphics window) and not accidentally the cold/hot end faces.

### 9.4 Set the fluid inlet mass flow rate

5. Right-click **Brinkman Equations (br) → Inlet**. 
   - **Boundary Selection**: the appropriate end face depending on flow direction (this must switch between phases 2 and 4 — see §10.4 for how Events handles the direction reversal).
   - **Boundary condition type**: **Mass flow**.
   - **Mass flow rate**: `mdot` = `0.084666[kg/s]`, active only during flow phases per §10.4's enable/disable logic; `0` (or the node disabled entirely) during magnetize/demagnetize.
6. Right-click **Brinkman Equations (br) → Outlet** on the opposite end face, with a **Pressure, no viscous stress** condition (`p0 = 0[Pa]`, i.e. a reference/gauge pressure outlet) — this is a standard companion boundary condition to a mass-flow inlet. *Software operation (not specified by the guide): the guide specifies the inlet mass flow rate but does not name a specific outlet condition — a pressure outlet is the standard COMSOL pairing for a mass-flow inlet and is used here for that reason.*

**Expected result before moving on:** four to five boundary-condition nodes (Cold reservoir BC, Hot reservoir BC, Thermal Insulation on walls, Inlet mass flow, Outlet pressure), all with valid (non-red) Boundary Selections and expressions.

**Common mistake:** applying the mass-flow Inlet condition to the *same* face throughout the whole cycle without reversing it for phases 2 vs 4 — the original guide explicitly requires **flow direction to reverse** between cold→hot (phase 2) and hot→cold (phase 4); a single static Inlet/Outlet pair will only correctly model one of those two phases.


---

## 10. Events for the AMR cycle (COMSOL Events interface)

**Goal:** implement the four-step cycle from the original guide's §5 — adiabatic magnetization, cold-to-hot flow, adiabatic demagnetization, hot-to-cold flow — switching the MCE source term, flow on/off, flow direction, and which reservoir BC is active, all automatically as the time-dependent solve progresses.

**Why this step is needed:** without Events (or an equivalent explicit state machine), COMSOL has no way to know that flow should be off during magnetization, on (in one direction) during phase 2, off again during demagnetization, and on (reversed) during phase 4 — this is the software mechanism that turns four separate physical sub-processes into one coherent time-dependent study.

### 10.1 The four-step cycle (reproduced unchanged from the original guide's §5)

1. **Adiabatic magnetization** (`H: 0 → H_max`): MCE source term active, fluid flow OFF.
2. **Cold-to-hot flow**: MCE source term OFF, fluid flow ON (cold→hot direction), duration set by the blow fraction (`BLOW_FRACTION_MASCHE`, default 0.5 — symmetric; check if a specific device documents its own blow fraction, since most benchmark devices in `data/amr_experimental_benchmarks.csv` implicitly assume 0.5).
3. **Adiabatic demagnetization** (`H: H_max → 0`): MCE source term active with sign flipped, fluid flow OFF.
4. **Hot-to-cold flow**: MCE source term OFF, fluid flow ON (hot→cold direction).

Total period = `1/frequency`. With `frequency = 0.75 Hz` (Step 13): period = `1.333 s`. With `BLOW_FRACTION_MASCHE = 0.5` applied symmetrically across the two flow phases and the remaining time split between the two static phases: each phase gets `1.333/4 = 0.333 s` (matching `t_mag` from §8.1).

### 10.2 Add the Events interface

1. Right-click **Component 1 → Add Physics**, or in the ribbon **Physics → Add Physics**. Under **Mathematics → ODE and DAE Interfaces**, or search **"Events"** directly — select **Events (ev)**. Add it to the component.

*Software operation (not specified by the guide) — the guide names the Events interface but doesn't give exact node-by-node settings; the discrete-state / explicit-event configuration below is the standard COMSOL pattern for a phase-switching cyclic problem.*

### 10.3 Define a Discrete State for cycle phase

2. Right-click **Events (ev) → Discrete States → Discrete State**. Name the state variable `phase` (integer-like, values 1–4 corresponding to the four steps above). Set **Initial value**: `1`.
3. Also add a discrete state `cycle_sign` (used in §8's source-term expression): **Initial value**: `1`.

### 10.4 Define Explicit Events to switch phase at the right times

4. Right-click **Events (ev) → Explicit Event**. Create four such nodes, one per phase transition, each firing at a specific time within the period and repeating every period:
   - **Event 1 (start of magnetize)**: **Start time** `0`, **Time interval** `1.333[s]` (i.e., recurs every period). **Reinitialize variable**: set `phase = 1`, `cycle_sign = 1`.
   - **Event 2 (start of cold-to-hot flow)**: **Start time** `0.333[s]`, **Time interval** `1.333[s]`. Set `phase = 2`, `cycle_sign = 0`.
   - **Event 3 (start of demagnetize)**: **Start time** `0.667[s]`, **Time interval** `1.333[s]`. Set `phase = 3`, `cycle_sign = -1`.
   - **Event 4 (start of hot-to-cold flow)**: **Start time** `1.000[s]`, **Time interval** `1.333[s]`. Set `phase = 4`, `cycle_sign = 0`.
5. Back in **Heat Transfer in Porous Media (ht) → Heat Source** (§8.3), the expression `... * cycle_sign` now automatically tracks the correct sign and zero-out per phase via this Discrete State.

### 10.5 Gate the flow boundary conditions on `phase`

6. Go back to **Brinkman Equations (br) → Inlet** (§9.4). Change the **Mass flow rate** expression to:

```
mdot * (phase==2) - mdot * (phase==4)
```

using COMSOL's boolean-to-number coercion (`phase==2` evaluates to `1` when true, `0` when false) so the sign automatically reverses between phase 2 (cold→hot, positive) and phase 4 (hot→cold, negative) at the **same** boundary — this avoids needing two separate Inlet nodes on two different faces with manual enabling/disabling, and is the more robust standard pattern for a reversing-flow BC in COMSOL.
7. Similarly, gate the Cold/Hot reservoir Temperature BCs (§9.2) so each is only meaningfully "in effect" during its correct phase — e.g., **Hot reservoir BC** temperature expression: `T_hot` but with the BC node itself toggled via **Events (ev) → Explicit Event**'s ability to enable/disable physics features (right-click the Explicit Event node → check for an "Enable/Disable" action target, or alternatively use each BC's own **Active** checkbox driven by a boolean expression like `phase==2` if COMSOL's version supports conditional activation on that BC type).

**Expected result before moving on:** running a short test solve (even before full mesh/solver tuning) shows `phase` stepping through 1→2→3→4→1... at exactly `0.333 s` intervals, and `cycle_sign`/the inlet mass flow direction switching in lockstep.

**How to verify:** after a trial time-dependent solve, plot `phase` and `cycle_sign` vs. time using **Results → 1D Plot Group → Global** — confirm a clean square-wave/step pattern with the correct period (1.333 s) and phase durations (0.333 s each).

**Common mistake:** off-by-one timing errors in the four **Start time** values (e.g., using `0.25×period` intervals that don't actually sum to one full period, or forgetting the events must repeat every period via **Time interval**, not just fire once at t=0) — always verify the plotted `phase` signal in the check above before trusting anything downstream.


---

## 11. Mesh generation

*Software operation (not specified by the guide) — standard COMSOL meshing workflow; the guide specifies geometry and physics but not explicit mesh density.*

**Goal:** produce a mesh fine enough that the solution doesn't change appreciably under refinement (mesh independence), while keeping element count manageable for a multi-cycle transient solve.

**Why this step is needed:** an under-resolved mesh is explicitly called out in the original guide's §8 as one of the three most likely explanations (along with a units error or wrong blow-fraction assumption) if the degeneracy check fails — meshing correctly the first time avoids wasted debugging later.

1. Click **Mesh 1** in the Model Builder.
2. **Sequence type**: choose **User-controlled mesh** (gives more control than Physics-controlled for this thin, high-aspect-ratio geometry — especially the parallel-plate track's 0.1 mm gap against a 0.170 m length, an aspect ratio of ~1700:1).
3. For the **2-D axisymmetric packed-bed track**:
   - Right-click **Mesh 1 → Mapped**, apply to the rectangular domain — a mapped (structured, quadrilateral) mesh is natural for this simple rectangular cross-section.
   - Under **Mapped 1 → Distribution**, set **Number of elements** along the axial (z) direction: start with `100`; along the radial (r) direction: `20`.
4. For the **3-D packed-bed track**:
   - Right-click **Mesh 1 → Swept**, with the sweep direction along the bed's axial length — this keeps element count manageable versus a fully unstructured tetrahedral mesh.
   - Set **Distribution**: `100` elements along the sweep (axial) direction. For the cross-section, use **Free Triangular** with a target element size scaled to the bed radius (e.g., **Maximum element size** ≈ radius/10).
5. For the **parallel-plate tracks (2-D and 3-D)**:
   - Because the fluid gap (0.1 mm) is far smaller than the bed length (0.170 m), use a **Mapped** (2-D) or **Swept** (3-D) mesh with **strong boundary-layer-style clustering** in the gap-normal direction: right-click the mesh distribution on the gap-normal edges and set **Element ratio** (a stretching factor, e.g. 5–10) with more elements clustered near the fluid/solid interface, since that's where the steepest temperature and velocity gradients occur.
   - Axial (flow-direction) elements: `100` as a starting point, matching the packed-bed track for consistency.

**Screens that should appear:** after clicking **Build All**, the Graphics window shows a full quadrilateral (2-D) or hexahedral/prism (3-D swept) mesh with no red error markers; the **Mesh Statistics** panel (right-click **Mesh 1 → Statistics**) reports the element count and quality metrics.

**Expected result before moving on:** a **Minimum element quality** (in Mesh Statistics) above ~0.1 (COMSOL's rule-of-thumb threshold for "acceptable" quality; higher is better) with no visibly degenerate (near-zero-area) elements, especially in the thin parallel-plate gap.

**How to verify:** right-click **Mesh 1 → Statistics** — check **Minimum element quality**, **Average element quality**, and total element count. Also visually zoom into the parallel-plate gap region to confirm the boundary-layer clustering produced well-shaped (not excessively skewed) elements.

### 11.1 Mesh independence check

Following the same discipline `core/geometry_analysis.py`'s own sweeps apply at the 0-D level (multiple resolution points per sweep, not a single value):

6. Solve once at the mesh from Steps 3–5 (after completing Steps 12–13 below), record Qc at periodic steady state.
7. Halve the mesh element counts in the flow/axial direction (e.g., 100 → 50), re-solve, record Qc again.
8. Double the original element counts (100 → 200), re-solve, record Qc a third time.

**Expected result:** Qc changes by less than ~2% between the 100-element and 200-element cases (using the ANSYS guide's own mesh-independence tolerance as a reasonable cross-tool standard, since the original COMSOL guide doesn't state a numeric tolerance itself — *software operation (not specified by the guide)*, borrowed here from the companion ANSYS guide's explicit §3 tolerance for consistency between the two tools).

**Common mistake:** declaring "mesh independence" after checking only one refinement level, or refining only in one direction (e.g., only axial) while leaving a genuinely under-resolved cross-stream/radial mesh — always check at least two refinement levels bracketing your working mesh.


---

## 12. Solver settings and time stepping

*Software operation (not specified by the guide) — standard COMSOL time-dependent study configuration; the guide specifies the physics and cycle structure but not explicit solver tolerances/time-step sizes.*

**Goal:** configure a Time Dependent study that can resolve the fast-switching Events (phase changes every 0.333 s) accurately while running efficiently over many periods.

1. In the ribbon, click **Study → Add Study**. Select **Time Dependent**. Click **Add Study** (or **Done**) to insert **Study 1 → Step 1: Time Dependent**.
2. Click **Step 1: Time Dependent**.
   - **Times**: use `range(0, 0.01, N_periods*1.333)` syntax in the Times field (COMSOL's `range(start, step, stop)`), where a `0.01 s` output-storage interval is fine-grained enough to resolve within-phase transients without generating an excessive dataset; `N_periods` is set per Step 13 below (start with `N_periods = 20` for an initial convergence check, per §14).
   - **Tolerance**: switch to **User controlled**; set **Relative tolerance**: `1e-4` (tighter than COMSOL's default `1e-3` — the sharp on/off switching of the MCE source term and flow BCs at each Event benefits from tighter tolerance to avoid the solver skating past a discontinuity).
3. Right-click **Step 1: Time Dependent → Time-Dependent Solver 1** (visible after right-clicking **Study 1** and selecting **Show Default Solver**).
   - Under **Time Stepping**, confirm the solver uses COMSOL's default **BDF** (Backward Differentiation Formula) method with **Free** step selection, but set **Maximum step**: `0.01[s]` — this caps the solver's automatic step size so it can't step over an Event's instantaneous phase change (a much larger free step could otherwise straddle a discontinuity and produce an inaccurate or non-converging step).
   - Confirm **Events** are being respected: the solver settings should show an **Events** node was automatically added when the Events (ev) interface was added to the component (Step 10.2) — this ensures the solver stops exactly at each Event's trigger time rather than stepping past it.

**Expected result before moving on:** the Study tree shows **Step 1: Time Dependent** with the modified Times/Tolerance settings, and the auto-generated solver sequence includes an **Events** solver node.

**How to verify:** right-click **Study 1 → Compute** (or click the **Compute** icon in the ribbon) on a short trial run (`N_periods = 2`) first — confirm it completes without solver-convergence errors and that the **Convergence Plot** (from **Messages** or **Progress** log) shows a reasonable, non-diverging iteration history at each time step.

**Common mistake:** leaving **Maximum step** unset (COMSOL default "Free," letting the solver pick arbitrarily large steps) — this is a very common way for Events-driven models to silently miss or blur a phase transition, since the adaptive step-size algorithm has no inherent reason to slow down near a discontinuity unless you tell it to respect the Events.

---

## 13. Operating point and running multiple cycles

**Goal:** run the model at the exact literature-calibrated operating point from the original guide's §7, for enough periods that the bed reaches periodic steady state.

### 13.1 Operating point (reproduced unchanged from the original guide's §7)

| Parameter | Value |
|---|---|
| Material | Gd (single-Tc approximation of the real Curie-graded 11-layer bed — same simplification `core/validation_system.py` uses for this row) |
| μ0H_max | 1.13 T |
| Regenerator mass | 1.7 kg |
| Frequency | 0.75 Hz |
| T_cold | 289 K (per `loss_model.py`'s calibration comment) |
| Span | 10.2 K |
| Calibrated mdot | 0.084666 kg/s (reproduces Qc=102.8 W exactly under the 0-D model) |

Source: Eriksen, Engelbrecht, Bahl, Bjørk, Nielsen, Insinga, Pryds, "Design and experimental tests of a rotary active magnetic regenerator prototype," *Int. J. Refrigeration* (2015), doi:10.1016/j.ijrefrig.2015.05.004. Reported result: Qc=102.8 W at span=10.2 K, COP=3.1 ("the COP of 3.1 is 11.3% of the Carnot efficiency"). This repo's own 0-D model reproduces this to −2.1% error (`results/pipeline.log`, step 2: `DTU_Eriksen_rotary_Gd_2015 ... err=-2.1% implied_parasitic=0.255`).

Enter every one of these values into **Global Definitions → Parameters** before running, using the exact parameter names referenced throughout this tutorial (`mu0H_max = 1.13[T]`, `frequency = 0.75[Hz]`, `T_cold = 289[K]`, `span = 10.2[K]`, `mdot = 0.084666[kg/s]`, mass used only to cross-check the derived bed length in §4.1).

### 13.2 Run multiple full cycles

The original guide's §5 requires: *"run multiple full cycles and check the bed's periodic-steady-state temperature profile has converged (cycle-to-cycle ΔT at any point below some tolerance, e.g. 0.01K) before reading off Qc — a single-cycle result from a bed starting at uniform T is not the AMR's actual operating point."*

1. Set **Times** in Step 12.2 to run at least `N_periods = 20` full periods initially (`20 × 1.333 s = 26.67 s` total simulated time) — a starting point, not a guarantee of convergence; §14 covers how to check and extend if needed.
2. Click **Compute**. This is a computationally expensive step (Events-driven transient with a fine mesh, run for tens of periods) — expect run times from minutes to hours depending on mesh size and machine; consider running unattended (COMSOL supports batch/background solves via **File → Save and Run in Background**, or the command-line `comsolbatch` utility for headless runs on a workstation/cluster).

**Expected result before moving on:** the solve completes without solver failure, and the Study log shows all `N_periods` fully computed (check the **Progress** window's reported final time equals `N_periods × 1.333 s`).

**Common mistake:** starting the bed at a non-physical initial condition (e.g., uniform T = 294 K everywhere) and reading off Qc from the *first* cycle — the original guide is explicit that this is not the AMR's actual operating point; always run to periodic steady state (§14) before trusting any Qc number.


---

## 14. Convergence checking (periodic steady state)

**Goal:** confirm the bed's cycle-to-cycle temperature profile has stopped changing before reading off Qc, per the original guide's §5 tolerance of cycle-to-cycle ΔT below ~0.01 K at any point.

*Software operation (not specified by the guide) — the guide states the 0.01 K tolerance and the requirement itself but not the exact COMSOL post-processing clicks; standard workflow below.*

1. After the Step 13 solve completes, go to **Results → 1D Plot Group**. Rename: **Cycle Convergence Check**.
2. Right-click **Cycle Convergence Check → Point Graph**.
   - **Selection**: pick 2–3 representative points in the bed (e.g., the cold-end face, mid-bed, hot-end face).
   - **y-axis data**: **Temperature (ht.T)**. 
   - **x-axis data**: **Time**.
3. Click **Plot**. You should see the point temperatures oscillating with the cycle period but with a decaying envelope — early cycles show large swings as the bed moves away from its uniform initial condition, later cycles should show the oscillation settle into a repeating pattern.
4. To check the quantitative tolerance: right-click **Results → Derived Values → Point Evaluation**, evaluate temperature at the same point(s) at the *end* of each period (`t = k × 1.333 s` for `k = 1, 2, ..., N_periods`), export to a table (**Table 1** or Excel export), and compute the cycle-to-cycle difference `ΔT_k = T(k×period) - T((k-1)×period)`.

**Expected result:** `|ΔT_k|` at every monitored point drops below `0.01 K` by the final few periods of the Step 13 run.

**How to verify:** the exported table's last few `ΔT_k` rows should all read below `0.01`. If not, extend `N_periods` in Step 12.2's Times field (e.g., 20 → 40 → 60) and re-run (COMSOL can resume from Step 13's stored solution rather than restarting from t=0, via **Study → Continue** or by extending the Times range and re-computing, depending on version) until the tolerance is met.

**Common mistake:** checking convergence only by eye on the plot (Step 3) without the quantitative table check (Step 4) — a plot that "looks" converged can still hide a slow residual drift larger than 0.01 K, especially near the bed ends where the reservoir BCs impose the largest cycle-to-cycle temperature swings.

---

## 15. Post-processing

**Goal:** extract Qc and COP from the converged periodic-steady-state solution, matching the definitions the 0-D benchmark (Step 13.1) uses, so the comparison in Step 16 is apples-to-apples.

*Software operation (not specified by the guide) — the guide specifies the target quantities (Qc, COP) but not the exact COMSOL derived-value clicks.*

### 15.1 Compute Qc (cooling capacity)

1. Right-click **Results → Derived Values → Surface Integration** (2-D/axisymmetric) or **Volume Integration**/**Surface Integration** as appropriate for a 3-D face integral.
2. **Selection**: the cold-end boundary face.
3. **Expression**: the convective heat flux drawn from the cold reservoir into the bed during phase 4 (hot-to-cold flow) — in COMSOL's Heat Transfer in Porous Media interface, this is available as a built-in boundary heat-flux variable (e.g., `ht.ntflux`, the normal total heat flux) integrated over the cold-end face **and** restricted to phase 4 (multiply the integrand by `(phase==4)` so only the cold-drawing part of the cycle contributes), then divided by the cycle period to get an average (time-averaged) cooling power in Watts — matching how the 0-D model's Qc is a per-cycle-averaged quantity, not an instantaneous one.
4. Evaluate this integral over the **final converged period only** (i.e., integrate from `t = (N_periods-1)×1.333 s` to `t = N_periods×1.333 s`), not over the whole transient including the non-periodic early cycles.

### 15.2 Compute COP

5. COP requires the net work input in addition to Qc. The original guide's §7-8 doesn't specify a COMSOL-side work calculation directly (the 0-D `implied_parasitic=0.255` figure comes from the Python-side `loss_model.py`, not from anything the 2-D model computes) — *software operation (not specified by the guide):* if you want a COMSOL-derived COP rather than just comparing Qc, you would need to separately integrate the magnetic work input (from the prescribed `H(t)` waveform and the material's magnetization curve — genuinely out of scope for this guide per its own §1, since it assumes a prescribed field, not a solved magnetic circuit) and the pumping work (from the pressure drop computed in §7.3 times the volumetric flow rate). **For the Step 16 degeneracy check, Qc alone is the primary comparison target** — the original guide's own §8 checks Qc specifically, not COP, so this is consistent with the guide's own validation approach.

### 15.3 Visualize the temperature field

6. Right-click **Results → 2D Plot Group** (or 3D Plot Group). Add a **Surface** plot of **Temperature (ht.T)**, evaluated at a specific point in the cycle (e.g., end of phase 2, the point of maximum hot-end temperature) — this is the genuinely new information this 2-D/3-D model adds over the 0-D model (real axial temperature gradients within the bed, per the original guide's §9).
7. Add an **Animation** (right-click **Results → Export → Animation**) sweeping through one full converged period, to visualize the demagnetization-front propagation and thermal wave motion through the bed.

**Expected result before moving on:** a single time-averaged Qc value (Watts) from the final converged period, plus a spatial temperature-field visualization showing a genuine axial gradient (not a flat/uniform bed — a flat result here suggests either non-convergence, Step 14, or a boundary-condition error, Step 9).

**How to verify:** the temperature field plot should show T monotonically varying (roughly) from `T_cold` near the cold end to `T_hot` near the hot end at the moment of peak flow, with visible curvature/nonlinearity reflecting the porous thermal-dispersion and MCE source effects — not a straight line (which would suggest the MCE source term isn't actually contributing) and not a flat/uniform field (which would suggest the flow or BCs aren't actually driving the bed).

**Common mistake:** integrating Qc over the *entire* transient (including the non-periodic early cycles where the bed is still warming/cooling from its initial condition) instead of restricting to the final converged period — this systematically biases the reported Qc and makes the Step 16 comparison meaningless.

---

## 16. Validation against the benchmark (degeneracy check)

**Goal — reproduced unchanged from the original guide's §8:** *"Before drawing any conclusion from the 2-D/3-D model that isn't already in this repo, reproduce this one number:"*

> **Qc = 102.8 W at span = 10.2 K, μ0H = 1.13 T, mass = 1.7 kg, f = 0.75 Hz, mdot = 0.084666 kg/s.**

### 16.1 Compare

1. Take the Step 15.1 time-averaged Qc from the final converged period.
2. Compute percent error: `error = (Qc_COMSOL - 102.8) / 102.8 × 100%`.

### 16.2 Interpret the result

> ⚠️ **Preserve this warning, unchanged:** *"If the COMSOL model's Qc at periodic steady state is not within roughly the same error band the 0-D model already achieves (−2.1%, or generously ±10-15% to allow for genuine spatial effects the 0-D model can't capture), do not trust any new geometry/gradient conclusion from the model until the discrepancy is understood — it more likely means a units error, a wrong blow-fraction assumption, or an under-resolved mesh than a genuine new physical finding."*
>
> This mirrors the standard this repo already holds itself to everywhere else (`core/validation.py`, `core/giguere_validation.py`, `core/validation_system.py`) — **a new, more sophisticated model earns trust by first reproducing what the simpler, already-checked model gets right, not by producing an interesting-looking number no one has checked.**

**How to verify pass/fail:** if `|error| ≤ 10–15%`, proceed to draw conclusions about axial gradients, geometry optima, etc. (Step 17/§9 topics below). If `|error| > 15%`, work back through: (a) unit tags on every parameter (Step 9's common mistake), (b) the blow-fraction assumption (Step 10.1 — confirm 0.5 is right for this device, or find its documented value), (c) mesh independence (Step 11.1 — re-check at finer resolution), before treating the mismatch as a genuine physical finding.

### 16.3 What a working 2-D/3-D model would add beyond the 0-D model (reproduced unchanged from the original guide's §9)

- Real axial temperature gradients within the bed (the 0-D model represents the whole bed by one `eps` number).
- Whether the linear `span_fraction = max(0, 1 - T_span/(2·dTad_noload))` approximation in `amr_cycle.py::cooling_capacity()` — which that function's own docstring flags as producing "a sharper, straight-line cutoff... than a real AMR device would show" — is actually a reasonable approximation, by directly resolving the spatial temperature profile near the no-load span limit instead of assuming a shape with no literature source.
- A genuine check on whether the packed-bed/parallel-plate COP optima `core/geometry_analysis.py` finds at a *fixed representative mdot* (documented there as a real methodological simplification, since free mdot optimization is degenerate in this repo's 2nd-law work model) survive when mdot and geometry are optimized jointly with spatially resolved thermal-hydraulics.
- For a Curie-graded bed specifically (out of scope for this guide, but the natural next extension): whether `cascade.py`'s treatment of each graded layer as independently peak-tuned, ignoring inter-layer axial conduction, materially changes the graded-cascade Qc/COP numbers in `results/graded_cascade_comparison.csv`.


---

## 17. Known limitations of this guide itself (reproduced unchanged from the original guide's §10)

- Never built or solved by the original guide's author — see §0.
- The MCE source-term lookup table (Step 5) inherits the mean-field model's own documented ~+30-50% overprediction of ΔT_ad near Tc at low field (`core/validation.py`); a 2-D model built on it will reproduce that bias faithfully, not correct it.
- Water properties are treated as temperature-independent constants unless explicitly replaced per Step 6.4's note.
- No wall-conduction or radiative loss term is specified (adiabatic walls, Step 9) — this matches the 0-D model's scope but is itself a simplification worth flagging in any writeup.
- Magnetic field is assumed spatially uniform and prescribed as a time waveform, not solved from an actual magnet-circuit model.

---

## 18. Complete checklist

- [ ] COMSOL license confirmed to include Heat Transfer Module + Porous Media and Subsurface Flow Module (or CFD Module) (§2)
- [ ] `gd_dTad_vs_T_1p13T.csv` generated from `core/mce_material.py::GADOLINIUM` at `mu0H = 1.13 T`, T = 270–310 K (Step 5.1)
- [ ] Model created with correct dimension: 2D Axisymmetric and/or 3D, per geometry track (Step 3.1)
- [ ] Brinkman Equations (br) and Heat Transfer in Porous Media (ht) physics interfaces added (Step 3.2)
- [ ] Geometry built: packed-bed (2D-axi and/or 3D) and/or parallel-plate (2D and/or 3D), per your chosen track(s) (Step 4)
- [ ] Interpolation function imported (`dTad_lookup`, `Ctotal_lookup`), linear interpolation, constant extrapolation, unit `K` on the argument (Step 5.2)
- [ ] Gd solid material created: ρ=7900 kg/m³, Cp=`Ctotal_lookup(T)`, k (literature value, flagged as not guide-sourced) (Step 6.2)
- [ ] Water material created: ρ=997 kg/m³, Cp=4186 J/(kg·K), μ=8.9e-4 Pa·s, k=0.606 W/(m·K) — or IAPWS-97 built-in library material if upgraded (Step 6.3/6.4)
- [ ] Porosity ε=0.365 entered in **both** Brinkman and Heat Transfer in Porous Media interfaces (packed-bed track) (Step 6.5)
- [ ] Friction-factor Analytic function created and fit into the Brinkman drag model — NOT the default Ergun equation (Step 7)
- [ ] MCE Heat Source node added to solid domain(s), referencing `dTad_lookup`, `Ctotal_lookup`, `t_mag`, `cycle_sign` (Step 8)
- [ ] Cold reservoir BC (T_cold=289 K), Hot reservoir BC (T_hot=299.2 K), adiabatic walls, mass-flow inlet (mdot=0.084666 kg/s) with phase-gated direction reversal (Step 9)
- [ ] Events interface added; Discrete States `phase` and `cycle_sign`; four Explicit Events firing every 0.333 s within a 1.333 s period (Step 10)
- [ ] Mesh built and mesh-independence checked (≤2% Qc change between refinement levels) (Step 11)
- [ ] Time Dependent study configured: relative tolerance 1e-4, maximum step 0.01 s, Events respected (Step 12)
- [ ] Solved for ≥20 full periods at the DTU_Eriksen_rotary_Gd_2015 operating point (Step 13)
- [ ] Cycle-to-cycle ΔT < 0.01 K confirmed at monitored points before reading off Qc (Step 14)
- [ ] Time-averaged Qc extracted from the final converged period only, over the cold-end face during phase 4 (Step 15)
- [ ] Degeneracy check performed: `|error| ≤ 10–15%` vs. Qc=102.8 W target (Step 16)

## 19. Parameter table (single reference, all values used in this tutorial)

| Parameter | Symbol | Value | Unit | Where entered in COMSOL |
|---|---|---|---|---|
| Applied field | `mu0H_max` | 1.13 | T | Global Definitions → Parameters |
| Regenerator mass | — | 1.7 | kg | Used to derive bed length (§4.1); not entered directly |
| Cycle frequency | `frequency` | 0.75 | Hz | Global Definitions → Parameters |
| Cold-side temperature | `T_cold` | 289 | K | Global Definitions → Parameters; Cold reservoir BC |
| Temperature span | `span` | 10.2 | K | Global Definitions → Parameters |
| Hot-side temperature (derived) | `T_hot` | 299.2 | K | `T_cold + span` expression; Hot reservoir BC |
| Calibrated mass flow rate | `mdot` | 0.084666 | kg/s | Global Definitions → Parameters; Brinkman Inlet |
| Particle diameter (packed bed) | `d_p` | 0.5 (default; sweep 0.05–2) | mm | Geometry sphere/cylinder dimension |
| Porosity (packed bed) | `ε` (`eps`) | 0.365 | — | Brinkman Fluid and Matrix Properties; ht Porous Medium |
| Bed cross-section | — | 0.002 | m² | Used to derive geometry radius/footprint (§4.1) |
| Derived bed length | `L_bed` | 0.170 | m | Rectangle/Cylinder height dimension |
| Plate spacing (parallel plate) | — | 0.1–0.2 (0.1 optimum used) | mm | Fluid-channel rectangle height |
| Plate thickness (parallel plate) | — | 0.5 | mm | Solid-plate rectangle height |
| Magnetization ramp duration | `t_mag` | 0.333 (derived: period/4) | s | Global Definitions → Parameters; Heat Source expression |
| Cycle period | — | 1.333 (=1/frequency) | s | Time Dependent study Times field; Events intervals |
| Gd solid density | `rho_solid` | 7900 | kg/m³ | Global Definitions → Parameters; Gd Material |
| Gd lattice Cp (flat, use only away from Tc) | — | 236 | J/(kg·K) | Not used directly — superseded by `Ctotal_lookup(T)` |
| Gd Curie temperature | `Tc` | 294.0 | K | Reference only (embedded in the imported lookup table) |
| Gd Debye temperature | `θ_D` | 169.0 | K | Reference only (embedded in the imported lookup table) |
| Gd molar mass | — | 157.25 | g/mol | Reference only (embedded in the imported lookup table) |
| Water density | — | 997 | kg/m³ | Water Material |
| Water specific heat | — | 4186 | J/(kg·K) | Water Material |
| Water dynamic viscosity | — | 8.9×10⁻⁴ | Pa·s | Water Material |
| Water thermal conductivity | — | 0.606 | W/(m·K) | Water Material |
| Packed-bed friction factor | `f` | `23.462·Re^-0.6716` (10<Re<5e5) | — | Analytic function `f_friction_packedbed` |
| Parallel-plate friction factor | `f` | `24/Re` (Re<2300) | — | Analytic function `f_friction_plate` |
| Solver relative tolerance | — | 1e-4 | — | Time Dependent Step → Tolerance |
| Solver maximum time step | — | 0.01 | s | Time-Dependent Solver 1 → Time Stepping |
| Target Qc (degeneracy check) | — | 102.8 | W | Comparison target, Step 16 |
| Target COP | — | 3.1 | — | Reference only (0-D model comparison) |
| Acceptable Qc error band | — | ±10–15 (0-D model's own error: −2.1%) | % | Step 16 pass/fail criterion |

## 20. Boundary-condition summary

| Boundary | Node type | Active phase(s) | Value |
|---|---|---|---|
| Cold end | Temperature | Phase 4 (hot-to-cold flow) | T_cold = 289 K |
| Hot end | Temperature | Phase 2 (cold-to-hot flow) | T_hot = 299.2 K |
| Bed/plate walls | Thermal Insulation (adiabatic) | Always | q=0 |
| Fluid inlet/outlet | Mass flow (Inlet) / Pressure (Outlet) | Phases 2 & 4 only, direction reverses | mdot = 0.084666 kg/s, sign = `(phase==2) - (phase==4)` |
| Solid domain | MCE Heat Source | Phases 1 & 3 only | `q_MCE = ρ·C_total(T)·ΔT_ad(T)/t_mag`, sign flips between phase 1 (+) and phase 3 (−) |

## 21. Material-property table

| Material | Property | Value | Unit |
|---|---|---|---|
| Gadolinium (solid) | Density | 7900 | kg/m³ |
| Gadolinium (solid) | Specific heat | `Ctotal_lookup(T)` (temperature-dependent; flat 236 valid only away from Tc) | J/(kg·K) |
| Gadolinium (solid) | Curie temperature | 294.0 | K |
| Gadolinium (solid) | Debye temperature | 169.0 | K |
| Gadolinium (solid) | Molar mass | 157.25 | g/mol |
| Gadolinium (solid) | Thermal conductivity | ~10.5 (literature value — **not from the original guide**) | W/(m·K) |
| Water | Density | 997 (or IAPWS-97 if upgraded) | kg/m³ |
| Water | Specific heat | 4186 (or IAPWS-97) | J/(kg·K) |
| Water | Dynamic viscosity | 8.9×10⁻⁴ (or IAPWS-97) | Pa·s |
| Water | Thermal conductivity | 0.606 (or IAPWS-97) | W/(m·K) |

## 22. Solver settings summary

| Setting | Value |
|---|---|
| Study type | Time Dependent |
| Total simulated time | ≥ 20 × period = ≥ 26.67 s (extend if not converged per Step 14) |
| Output time step | 0.01 s |
| Relative tolerance | 1e-4 (User controlled) |
| Maximum solver time step | 0.01 s |
| Time-stepping method | BDF (default), Free step selection with capped maximum |
| Events handling | Automatic Events solver node (from Events (ev) interface, Step 10.2) |
| Mesh independence tolerance | ≤2% Qc change between refinement levels |

## 23. Validation checklist

- [ ] Lookup table (`gd_dTad_vs_T_1p13T.csv`) generated from the exact repo script in Step 5.1, not hand-fit
- [ ] Friction-factor correlation matches Tušek et al. Eq. 5 (packed bed) or Eq. 6 (parallel plate) exactly — not COMSOL's default Ergun equation
- [ ] Porosity ε=0.365 entered identically in both Brinkman and Heat Transfer interfaces
- [ ] Ran ≥20 full periods; cycle-to-cycle ΔT < 0.01 K at all monitored points before reading Qc
- [ ] Qc extracted from the final converged period only, restricted to phase 4 at the cold-end face
- [ ] Degeneracy check performed against Qc = 102.8 W (span=10.2 K, μ0H=1.13T, mass=1.7kg, f=0.75Hz, mdot=0.084666 kg/s)
- [ ] Result falls within ±10–15% (or ideally closer to the 0-D model's own −2.1% error) before trusting any new geometry/gradient conclusion
- [ ] Any new conclusion (axial gradients, geometry optimum shift, span_fraction shape) explicitly flagged as **beyond** what the 0-D model checks, per Step 16.3
- [ ] Every deviation from the 0-D model's assumptions (wall-loss terms, IAPWS water, resolved magnet field) explicitly flagged in any writeup, not silently introduced

## 24. Troubleshooting

| Symptom | Likely cause | Fix |
|---|---|---|
| Qc comes out orders of magnitude too large or small | Missing unit tag on a parameter (e.g., `mdot=0.084666` typed without `[kg/s]`) | Re-check every numeric entry in Global Definitions → Parameters and every Material property has an explicit unit tag |
| Bed temperature field is flat/uniform even after many periods | MCE Heat Source applied to the wrong domain, or `cycle_sign`/`phase` never actually changing | Re-check Step 8.3 Domain Selection; plot `phase` vs. time (Step 10.5) to confirm Events are firing |
| Solver fails to converge or takes excessively small time steps | Maximum step size not capped, solver stepping across an Event discontinuity | Set **Maximum step** = 0.01 s in Time-Dependent Solver settings (Step 12.3) |
| Cycle-to-cycle ΔT never drops below 0.01 K | Not enough periods run, or a genuine periodic limit-cycle with amplitude above tolerance | Extend `N_periods`; if still not converging after e.g. 100 periods, re-check for a modeling error (e.g., a BC that never actually turns off) |
| Pressure drop across the bed doesn't match the hand-calculated `f(Re)` value | Brinkman drag model still on COMSOL's Ergun default instead of the fitted User-Defined values | Re-check Step 7.3 — confirm **User defined** (not Ergun/Kozeny-Carman) is selected and the κ/Forchheimer fit was actually applied |
| Import of the Interpolation function fails or produces zero values | CSV column headers/units mismatched, or wrong number of arguments specified | Re-check Step 5.2 — confirm 1 argument, correct file path, and that COMSOL parsed 3 columns (T, ΔT_ad, C_total) not 1 |
| Geometry Build fails with an error about domain crossing r<0 | Axisymmetric geometry accidentally mirrored or offset into negative r | Re-check Step 4.1a — axisymmetric geometry is drawn only in the r≥0 half-plane |
| Parallel-plate mesh has degenerate/zero-area elements | Plate spacing entered in mm instead of m (off by 1000×), or mesh not clustered in the thin gap | Re-check units on plate spacing (Step 4.2a); apply boundary-layer-style element clustering (Step 11, parallel-plate bullet) |

## 25. Common errors and fixes

| Error message / symptom | Meaning | Fix |
|---|---|---|
| "Failed to find a solution" (Time Dependent solver) | Solver diverged, often from an overly large step near an Event | Cap Maximum step at 0.01 s (Step 12.3); tighten Relative tolerance to 1e-4 |
| "Undefined variable" (red underline in an expression field) | A referenced parameter, function, or discrete state hasn't been defined yet, or is misspelled | Cross-check spelling against Global Definitions → Parameters and Functions; confirm the Interpolation function names match `dTad_lookup`/`Ctotal_lookup` exactly |
| "Inconsistent units" | A numeric literal without a unit tag was combined with a unit-bearing quantity | Add explicit `[unit]` tags to every literal, per Step 6.3's common mistake |
| Geometry "Boolean operation failed" / gaps in Form Union | Overlapping or non-touching domains (common in the parallel-plate two-rectangle geometry) | Re-check Step 4.2a's exact Position values so domains share a single edge with no gap/overlap |
| Brinkman interface not available in Add Physics list | Porous Media and Subsurface Flow Module not licensed | Check `File → Help → About COMSOL Multiphysics`; contact your COMSOL administrator about module licensing (§2 prerequisite) |
| Mesh "Failed to create swept mesh" (3-D packed-bed or parallel-plate track) | Source/target faces for the sweep not correctly identified, common with complex or thin geometries | Re-check the Swept mesh node's **Source Face** / **Destination Face** selections; for very thin domains (parallel-plate gap), consider a Mapped mesh in a Work Plane cross-section extruded instead, per Step 4.2b |