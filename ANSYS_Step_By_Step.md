# ANSYS Fluent Step-by-Step Tutorial: 2-D/3-D AMR Regenerator-Bed Model

## How to use this tutorial

This is a click-by-click implementation of `ANSYS_Setup_Guide.md`, written for someone who has never opened ANSYS Workbench or Fluent before. It covers **both** geometry options from the original guide (packed-bed and parallel-plate) and **both** dimensionalities (2-D and full 3-D), and calls out clearly, at every step, where those four tracks differ. Where the source guide didn't specify how to actually operate the software, the step is explicitly labeled **"Software operation (not specified by the guide)"** and a standard ANSYS 2024 R1+ workflow is described instead. No equation, assumption, benchmark value, material property, or geometry value from the original guide has been changed.

Menu paths below assume **ANSYS Workbench 2024 R1 or later** with **Fluent**, **SpaceClaim**, and **Fluent Meshing** available. Where the original COMSOL guide is referenced for a shared value or definition, that content is inlined directly here rather than sending you to a different document, per the original guide's own instruction: *"Read `COMSOL_SETUP_GUIDE.md` first if you haven't — this guide assumes the same physics and mostly covers where the ANSYS workflow differs, rather than re-deriving everything from scratch."*

> ⚠️ **Preserve this warning from the original guide, unchanged:** *This is a setup specification, not a validated result*, for exactly the same reason as the COMSOL guide: no ANSYS license was available in the environment this repo was developed in, so nothing below has actually been built, meshed, or solved. This guide specifies the same physical model as the COMSOL guide — same geometry, same material properties, same MCE source term, same degeneracy-check target — translated into ANSYS Fluent's porous-media + UDF workflow instead of COMSOL's Brinkman-equations + interpolation-function workflow. **If both guides are ever built out, they should be checked against the same benchmark number (Step 20-21) and against each other, not just against the 0-D model independently** — a genuine agreement between two independently-implemented 2-D/3-D solvers would be much stronger evidence than either one alone.

## 0. Why this model exists / scope (background — no clicking yet)

`core/amr_cycle.py` is a lumped 0-D model; this tutorial builds a spatially-resolved single-cell regenerator model (packed-bed or parallel-plate, one AMR cycle) to check whether that 0-D model's structural approximations — the linear `span_fraction` cutoff, the fixed-mdot geometry sweeps in `core/geometry_analysis.py`, the phenomenological `(1 - 0.3U)` utilization correction in `core/thermal.py` — hold up once real gradients are resolved.

- **In scope:** a single-material (Gd), single-stage regenerator bed, packed-sphere OR parallel-plate geometry, one full AMR cycle (magnetize → cold-to-hot flow → demagnetize → hot-to-cold flow), reproducing the bed-internal temperature field and resulting Qc/COP at a literature-calibrated operating point.
- **Out of scope:** no Curie-graded multi-layer bed, no solved magnet-circuit field (a uniform, prescribed `H(t)` waveform is assumed throughout), no structural/mechanical analysis of the housing.

## 1. Tool choice within ANSYS

**Goal:** decide which ANSYS products to launch before you start clicking, since the wrong choice here (e.g. building this in CFX) makes later steps not match this tutorial.

- **ANSYS Fluent** — the tool this entire tutorial uses. Its **Porous Zone** model (Cell Zone Conditions → Porous Zone) natively supports the viscous/inertial resistance coefficients this problem needs, and its **UDF (User-Defined Function)** framework is the natural place to implement the MCE volumetric source term, since it needs to switch on and off with cycle phase and depend on local cell temperature — not a static material property.
- **ANSYS Mechanical / Maxwell** — out of scope for this guide (no structural or magnetic-circuit analysis specified). Would be the natural next tool if a follow-on step solves the actual permanent-magnet field distribution instead of prescribing a uniform `H(t)` waveform, or checks regenerator-housing thermal stress from the cyclic temperature swing.
- **ANSYS CFX** — a viable alternative to Fluent for the porous-flow part but has a less flexible user-source-term mechanism (CEL expressions vs. Fluent's compiled/interpreted UDFs) for something as cycle-dependent as the MCE source term. Fluent is the more direct match to what this problem needs — **the rest of this tutorial assumes Fluent.**

## 2. Prerequisites checklist

| Item | Requirement | Why |
|---|---|---|
| ANSYS license | Fluent (with Porous Zone capability, included in standard Fluent) + a C compiler for UDF compilation | Fluent's compiled UDFs (§13) require a supported compiler (Visual Studio on Windows, gcc on Linux) registered with Fluent |
| Lookup CSV file | `gd_dTad_vs_T_1p13T.csv` — **the exact same file** the COMSOL guide's §2.3 script produces | Both guides must read the same CSV so COMSOL and Fluent results stay directly comparable — do not regenerate it with different bounds for this guide |
| Python environment | Working `core/mce_material.py` importable | To generate the CSV (Step 4 below) |
| C compiler | Visual Studio (Windows) or gcc (Linux), Fluent-compatible version per your Fluent release notes | Required for **Compiled** UDFs (recommended over Interpreted for this problem's complexity) |
| Values on hand | The Step 6/7/20 tables below (geometry, material, operating point) | You will type or hard-code every one of these by hand |

**Common mistake:** attempting to load a Compiled UDF without a registered, compatible C compiler — Fluent will report a build failure referencing `udf` makefile errors; check `Help → About Fluent` and your compiler's installed version against ANSYS's published compiler-compatibility table for your Fluent release before troubleshooting the UDF code itself.

---

## 3. Workbench project setup

*Software operation (not specified by the guide) — standard ANSYS Workbench project creation, since the guide specifies the physics/geometry/UDF content but not the Workbench-level project scaffolding.*

**Goal:** create a Workbench project with the correct system chain: Geometry → Mesh → Fluent (Setup/Solution/Results), so the geometry and mesh you build feed directly into a Fluent case.

1. Launch **ANSYS Workbench** (Start Menu → ANSYS 2024 R1 → Workbench, or `runwb2` on Linux).
2. From the **Toolbox** on the left, expand **Analysis Systems**.
3. Drag **Fluid Flow (Fluent)** onto the empty **Project Schematic** canvas. This automatically creates a linked chain: **A2 Geometry → A3 Mesh → A4 Setup → A5 Solution → A6 Results**.
4. Rename the system (right-click the header cell → **Rename**): `AMR_PackedBed_2Daxi` (or an equivalent name reflecting your chosen track — repeat this whole Workbench system once per geometry/dimensionality track you're building, e.g. `AMR_PackedBed_3D`, `AMR_ParallelPlate_2D`, `AMR_ParallelPlate_3D`, each as its own **Fluid Flow (Fluent)** system on the same schematic).
5. Save the project (**File → Save As**): `AMR_regenerator_bed.wbpj`.

**Expected result before moving on:** the Project Schematic shows one **Fluid Flow (Fluent)** system per geometry/dimensionality track, each with cells A2–A6, all showing the blue "unfulfilled" question-mark icon (not yet built).

**How to verify:** hover over cell **A2 Geometry** — the tooltip should read "Geometry — Not yet defined starting from Workbench." No red error icons should be present yet (those appear only after you've started editing and left something inconsistent).

**Common mistake:** using a single **Fluid Flow (Fluent)** system and trying to swap geometries in and out of it for each track — this makes it very easy to accidentally carry over stale mesh/setup settings from a previous track. Use a **separate system per track**, as in Step 4 above; Workbench handles multiple parallel systems on one schematic natively.

---

## 4. Generating the MCE source-term lookup table

**Goal:** produce the exact same `(T, ΔT_ad, C_total)` CSV the COMSOL guide uses, so Fluent's UDF-based source term is provably tied to the same repo-validated material model.

**Why this step is needed:** per the original guide's §2.3, this is "identical to `COMSOL_SETUP_GUIDE.md` §2.3 — reuse the same CSV, don't regenerate it with different field/temperature bounds for the two guides, or the COMSOL and Fluent results stop being comparable to each other."

### 4.1 Generate the CSV (outside Fluent, in Python)

Run exactly this script (identical to the COMSOL guide's script, inlined here per the original guide's own instruction not to force the reader to jump elsewhere):

```python
import numpy as np
from core.mce_material import GADOLINIUM

mu0H = 1.13  # T -- matches the DTU_Eriksen_rotary_Gd_2015 benchmark, see Step 20
Ts = np.linspace(270, 310, 401)
dTad = GADOLINIUM.delta_T_adiabatic(Ts, mu0H / (4 * np.pi * 1e-7))
C = GADOLINIUM.total_heat_capacity(Ts, mu0H / (4 * np.pi * 1e-7))
np.savetxt("gd_dTad_vs_T_1p13T.csv",
           np.column_stack([Ts, dTad, C]), delimiter=",",
           header="T_K,dTad_K,C_total_J_per_kgK", comments="")
```

**Expected result:** a file `gd_dTad_vs_T_1p13T.csv` with 401 rows, 3 columns, T = 270–310 K.

**How to verify:** confirm the `dTad_K` column peaks near T = 294 K.

> ⚠️ **Preserve this warning:** this table inherits the mean-field model's own documented **+29% to +49% overprediction** of ΔT_ad relative to Dan'kov et al. (1998) at 1–2 T, improving to **−7.5%** at 5 T (`core/validation.py`). **A Fluent model built on it will reproduce that bias faithfully, not correct it.**

### 4.2 Keep this file accessible to the UDF

Copy `gd_dTad_vs_T_1p13T.csv` into the same working folder Fluent will use for this case (typically the Workbench project's `user_files` directory, or a fixed absolute path you'll hard-code into the UDF's file-open call in Step 13). Note the exact path — you will need it verbatim in the UDF source code.

---

## 5. Geometry creation (SpaceClaim)

The original guide's §3 specifies: *"Same unit-cell geometry as `COMSOL_SETUP_GUIDE.md` §3 (packed bed: d_p=0.5mm default, ε=0.365, 0.002 m² cross-section; parallel plate: 0.1-0.2mm spacing, 0.5mm plate thickness) — build the same CAD geometry in Fluent's meshing tool (Fluent Meshing / Watertight Geometry workflow, or SpaceClaim if extruding a swept 3-D unit cell)."* This tutorial uses **SpaceClaim** for CAD geometry creation (recommended for its more direct sketch/extrude workflow, especially for the thin parallel-plate geometry) and passes the result to **Fluent Meshing**.

### 5.1 Geometry values to enter (identical to the COMSOL guide, inlined here)

**Packed-bed unit cell:**

| Quantity | Value | Source |
|---|---|---|
| Particle diameter `d_p` | 0.5 mm default; sweep 0.05–2 mm | Phase 7 geometry sweep found COP-optimal packed-bed diameter at 0.5mm for the 291K/10K-span/2kg/1Hz/1.5T operating point — **not necessarily optimal** at the 1.13T/10.2K-span/1.7kg/0.75Hz benchmark used for the Step 21 degeneracy check |
| Porosity `ε` | 0.365 | `core/thermal.py`'s `regenerator_effectiveness()` default |
| Bed cross-section | 0.002 m² (~5×4 cm face) | `core/thermal.py` default |
| Solid density `ρ_Gd` | 7900 kg/m³ | `core/thermal.py::RHO_GD` (material property, listed here for completeness — enter in Fluent Materials, Step 11) |

Derived bed length (same arithmetic as the COMSOL tutorial, following directly from the guide's own stated quantities): `L_bed = m_bed / (ρ_Gd × (1-ε) × A_cross) = 1.7 / (7900 × 0.635 × 0.002) = 0.170 m`.

**Parallel-plate alternative:**

| Quantity | Value | Source |
|---|---|---|
| Plate spacing (fluid gap) | 0.1–0.2 mm typical; Phase-7 COP optimum at 0.1 mm | `core/geometry_analysis.py` |
| Plate thickness | 0.5 mm default | `core/thermal.py`'s `regenerator_effectiveness_parallel_plate()` |
| Hydraulic diameter | `d_h = 2 × spacing` | Eq. 7 slot limit |

### 5.2 Build the packed-bed geometry — 2-D track

*Software operation (not specified by the guide) — standard SpaceClaim sketch/extrude workflow.*

1. In the Workbench Project Schematic, double-click **A2 Geometry** to launch SpaceClaim.
2. In SpaceClaim, confirm units: **File → SpaceClaim Options → Units**, set **Length**: **Meter** (or work in mm and let SpaceClaim convert — but be consistent and double-check the final geometry's reported dimensions match the table above).
3. Click **Sketch Mode**, select the **XZ plane** (or XY, any consistent choice — just be consistent through meshing and Fluent setup).
4. Draw a **Rectangle**: one corner at the origin, opposite corner at `(r=0.0252, z=0.170)` meters — reusing the same circular-cross-section-equivalent radius as the COMSOL tutorial (`r = sqrt(0.002/π) ≈ 0.0252 m`).
5. Exit sketch mode. Since Fluent 2-D solves a genuine 2-D planar problem (Fluent does not have a distinct "2-D axisymmetric" geometry type the way COMSOL does — axisymmetry in Fluent is a **solver setting**, not a geometry-drawing choice), leave this sketch as a flat 2-D **surface** (do not extrude) if you plan to use Fluent's **Axisymmetric** solver option (Step 12.3); Fluent will treat the sketch's local x-axis as the axial direction and y-axis as radial when Axisymmetric is enabled.
6. Click **Pull** is not needed for a pure 2-D case — instead, right-click the sketch curve region and confirm it becomes a **Surface** body (not a Solid).

**Expected result before moving on:** a single 2-D surface body, 0.0252 m × 0.170 m, visible in the SpaceClaim **Structure** tree as `Surface1`.

**How to verify:** click the surface, check the **Properties** panel reports Area ≈ `0.0252 × 0.170 = 0.004284 m²` (this is the 2-D "footprint" area of the sketch itself, not the 0.002 m² physical bed cross-section — those are two different quantities: the sketch is a *radial-axial* slice, while the 0.002 m² is the bed's *flow* cross-section).

**Common mistake:** confusing the 0.002 m² bed cross-section (perpendicular to flow) with the 2-D sketch's own in-plane area (which is a radial-by-axial slice, a different plane entirely) — these should not numerically match, and trying to force them to is a units/geometry-interpretation error.

### 5.3 Build the packed-bed geometry — 3-D track

1. Launch SpaceClaim (as above), confirm units.
2. Sketch a circle of radius `0.0252 m` on the XY plane, centered at the origin (or a `0.05 × 0.04 m` rectangle if using the guide's "~5×4 cm face" rectangular-cross-section interpretation — same choice as the COMSOL tutorial's §4.1b).
3. Click **Pull** and extrude to **Height = 0.170 m** along Z, creating a solid cylinder (or rectangular prism).

**Expected result:** a single 3-D solid body, radius 0.0252 m (or 0.05×0.04 m footprint), height 0.170 m.

**How to verify:** **Properties** panel reports Volume ≈ `π × 0.0252² × 0.170 = 3.40×10⁻⁴ m³` (cylinder) or `0.05 × 0.04 × 0.170 = 3.40×10⁻⁴ m³` (rectangular prism — same volume by construction since both use the same 0.002 m²-equivalent cross-section × 0.170 m length).

### 5.4 Build the parallel-plate geometry — 2-D track

1. Sketch on a plane, draw two adjacent rectangles sharing one edge (identical construction to the COMSOL tutorial's §4.2a):
   - Fluid rectangle: `0.170 m` (flow direction) × `0.0001 m` (0.1 mm gap), corner at origin.
   - Solid plate rectangle: `0.170 m` × `0.0005 m` (0.5 mm), stacked directly above the fluid rectangle (corner at `(0, 0.0001)`).
2. Confirm the two rectangles form a single, continuous **Surface** with an internal shared edge (this shared edge becomes the fluid/solid interface boundary in Fluent) — in SpaceClaim, sketching two adjacent, touching rectangles and exiting sketch mode should automatically produce two separate surface bodies sharing one edge; if SpaceClaim instead merges them into a single body with no internal edge, use **Split Body** along the shared line to reintroduce the interface, since Fluent needs to distinguish the fluid zone from the solid zone.

**Expected result before moving on:** two adjacent 2-D surface bodies (fluid, solid) sharing exactly one edge, each with the correct dimensions.

**How to verify:** select each surface individually in the **Structure** tree, check reported Area: fluid = `0.170 × 0.0001 = 1.7×10⁻⁵ m²`; solid = `0.170 × 0.0005 = 8.5×10⁻⁵ m²`.

**Common mistake:** identical to the COMSOL tutorial's — entering plate spacing in millimeters where SpaceClaim expects meters (or vice versa, depending on your Options → Units setting) produces a channel 1000× the intended size. Always re-check the **Properties** panel's reported dimensions against the table in §5.1 before proceeding to meshing.

### 5.5 Build the parallel-plate geometry — 3-D track

1. Sketch the same two-rectangle cross-section as §5.4 on the XY plane.
2. **Pull**/extrude both bodies together (or extrude the sketch, then split) to a depth of `0.04 m` (same representative-depth choice as the COMSOL tutorial's §4.2b — *software operation (not specified by the guide)*, chosen for dimensional consistency with the packed-bed cross-section; flag this choice explicitly in any report since the original guide doesn't specify a plate width/depth).

**Expected result:** two adjacent 3-D solid bodies (fluid channel, solid plate), sharing one face, each `0.170 m × 0.04 m` in footprint, `0.0001 m` (fluid) or `0.0005 m` (solid) thick.

**How to verify:** Volume check — fluid ≈ `0.170 × 0.0001 × 0.04 = 6.8×10⁻⁷ m³`; solid ≈ `0.170 × 0.0005 × 0.04 = 3.4×10⁻⁶ m³`.

### 5.6 Share topology between fluid and solid bodies

*Software operation (not specified by the guide) — required for Fluent to recognize a conjugate (fluid/solid conduction-coupled) interface rather than two disconnected bodies.*

For **both** parallel-plate tracks (§5.4, §5.5): before leaving SpaceClaim, select both bodies in the Structure tree, right-click → **Share** (or use the **Component** grouping + **Share Topology** tool under the **Workbench** ribbon tab in SpaceClaim). This ensures the shared face/edge between fluid and solid is meshed with **conformal** (matching) nodes on both sides, which Fluent needs for a proper conjugate heat transfer interface — without this, Fluent Meshing may create two independent, non-matching meshes at the interface, which then requires a non-conformal interface (more error-prone and not what a clean conjugate setup needs).

**Common mistake:** skipping this step and discovering only after meshing that Fluent treats the fluid/solid boundary as two separate, unconnected walls rather than one coupled interface — always share topology in SpaceClaim (or the equivalent **Multizone**/**Shared Topology** setting in the Fluent Meshing app if building geometry there instead) before meshing.

---

## 6. Meshing (Fluent Meshing / Watertight Geometry workflow)

**Goal:** produce a mesh fine enough for mesh independence (Qc changes < 2% under refinement, per the original guide's explicit §3 tolerance) while keeping element count and solve time manageable for a multi-cycle transient UDF-driven solve.

**Why this step is needed:** the original guide's §3 recommends "a structured/swept mesh along the flow direction where possible (porous-zone flows don't need boundary-layer-resolved cells the way a fully-resolved packed-bed CFD would, since the porous model is already a volume-averaged approximation)."

1. In Workbench, double-click **A3 Mesh**. This launches either the classic ANSYS Meshing app or **Fluent Meshing** depending on your Workbench cell configuration — for this tutorial, use **Fluent Meshing's Watertight Geometry workflow** (right-click **A3 Mesh** → **Edit...** and confirm it's routed to Fluent Meshing, or launch Fluent Meshing directly and import the SpaceClaim geometry via **File → Import → Geometry**).
2. **Watertight Geometry workflow steps** (Fluent Meshing's guided task list, left panel):
   - **Import Geometry**: browse to the SpaceClaim `.scdoc`/exported CAD file. Confirm the units shown match meters.
   - **Add Local Sizing**: for all tracks, add a **Body of Influence** or **Face Sizing** control with a target element size appropriate to the geometry — for the parallel-plate tracks, add a much finer sizing control specifically on the thin fluid-gap region (target size ≈ gap/10, i.e. ~0.00001 m for the 0.1 mm gap) to properly resolve the steep gradients there.
   - **Generate the Surface Mesh**: accept defaults for surface mesh generation, then inspect for any red/flagged bad-quality faces before proceeding.
   - **Describe Geometry**: tell Fluent Meshing whether the geometry represents fluid, solid, or both regions with defined boundaries, so it can auto-detect zones — select **"Solid and fluid regions will be identified"** since this problem has both.
   - **Update Boundaries / Update Regions**: Fluent Meshing auto-names zones based on the geometry topology; rename them here to something meaningful (`bed_fluid`, `bed_solid` for the packed-bed track — recall the Brinkman-equivalent porous approach in Fluent still uses **one combined cell zone** treated as a **Porous Zone** for the packed bed, not physically separate fluid/solid sub-domains, since Fluent's Porous Zone model — like COMSOL's Brinkman formulation — represents the bed as a single homogenized zone with a porosity fraction, not two literal geometric regions; `fluid_channel`, `solid_plate` for the parallel-plate track, which **does** have genuinely separate geometric fluid and solid domains).
   - **Add Boundary Layers**: for the packed-bed porous-zone track, boundary layers are generally not needed (per the original guide's §3 note that porous-zone flows don't need boundary-layer-resolved cells). For the parallel-plate track, the fine sizing control from the Local Sizing step already provides adequate near-wall resolution given the extreme aspect ratio of the gap.
   - **Generate the Volume Mesh**: choose **Poly-Hexcore** (a good general-purpose choice balancing element count and quality) for the 3-D tracks, or accept the automatically-generated triangular/quadrilateral volume mesh for the 2-D tracks (Fluent Meshing solves 2-D problems as an extruded pseudo-3-D single-cell-thick mesh internally, but this is transparent to you).

**Screens that should appear:** each Watertight Geometry task shows a green checkmark as it completes; the final **Volume Mesh** step displays the mesh in the Graphics window with no red error/warning banners in the task list.

**Expected result before moving on:** a complete volume mesh with defined zones (`bed` porous zone for packed-bed tracks; `fluid_channel` + `solid_plate` for parallel-plate tracks), ready to switch to **Solution Mode** (the button to transition from Fluent Meshing into the Fluent solver).

**How to verify:** in the **Mesh → Check** step (automatically run, or manually via the Fluent Meshing console command `mesh/check`), confirm zero negative-volume cells and a reported minimum orthogonal quality above ~0.1 (a common rule-of-thumb acceptable-quality threshold).

**Common mistake:** skipping the fine local sizing on the parallel-plate fluid gap — a mesh sized for the 0.170 m bed length but not specifically refined in the 0.0001 m gap direction will produce only 1-2 cells across the gap, nowhere near enough to resolve the velocity and temperature profile Fluent's Porous Zone (or, for genuinely resolved flow, laminar) model needs.

### 6.1 Mesh independence check

Per the original guide's explicit §3 tolerance: *"Mesh independence check: halve the axial cell count and confirm Qc at the periodic-steady operating point changes by less than ~2% before trusting any geometry-sensitivity conclusion."*

1. After completing Steps 12–19 below and obtaining a converged Qc at your working mesh, return to Fluent Meshing (or re-run the Watertight workflow) with the axial cell count halved; re-solve; record Qc.
2. Repeat with the axial cell count doubled versus the original working mesh; re-solve; record Qc.

**Expected result:** Qc changes by less than ~2% between the halved and doubled cases (bracketing the original working mesh).

---

## 7. Named selections

**Goal:** name every mesh zone/boundary you'll need to reference in the Fluent solver setup (Cell Zone Conditions, Boundary Conditions, UDF domain hooks), so later steps can select them unambiguously by name rather than by clicking geometry each time.

*Software operation (not specified by the guide) — standard Fluent Meshing/Workbench named-selection workflow.*

1. In the Fluent Meshing **Named Selections** panel (or, if meshing was done in classic ANSYS Meshing, right-click faces/bodies in the Outline tree → **Create Named Selection**), create:

| Named selection | Applies to | Purpose |
|---|---|---|
| `cold_end` | The bed/channel end face at z=0 | Cold reservoir boundary condition (Step 15) |
| `hot_end` | The bed/channel end face at z=L_bed | Hot reservoir boundary condition (Step 15) |
| `bed_wall` (packed-bed) or `plate_outer_wall` (parallel-plate) | Lateral/outer walls | Adiabatic wall boundary condition (Step 15) |
| `bed_porous` (packed-bed) | The single homogenized bed cell zone | Porous Zone cell zone condition (Step 10) |
| `fluid_channel` (parallel-plate) | The fluid domain | Fluid cell zone / Inlet-Outlet boundaries |
| `solid_plate` (parallel-plate) | The solid domain | Solid cell zone / MCE source UDF hook |
| `fluid_solid_interface` (parallel-plate) | The shared fluid/solid face | Conjugate heat transfer coupled wall |

**Expected result before moving on:** all named selections appear in the Fluent Meshing (or Mesh app) Outline tree with the names above, each highlighting the correct geometric entity when clicked.

**How to verify:** click each named selection in turn; confirm the Graphics window highlights exactly the intended face(s)/body — a `cold_end` selection that accidentally includes part of `bed_wall`, for instance, will apply the wrong boundary condition type to part of the wall later.

**Common mistake:** naming `cold_end`/`hot_end` inconsistently between different Workbench systems (if building multiple geometry/dimensionality tracks as separate systems per Step 3) — keep the naming convention above identical across all tracks so the Fluent Setup steps below (Steps 8+) can be followed verbatim for each track without renaming anything.

---

## 8. Fluent launcher options

*Software operation (not specified by the guide) — standard Fluent launcher configuration.*

**Goal:** launch Fluent's solver with the correct dimension (2-D vs 3-D) and double precision, matching the mesh you built.

1. In Workbench, double-click **A4 Setup**. The **Fluent Launcher** dialog appears.
2. **Dimension**: select **2D** or **3D** matching the geometry track you built in Step 5.
3. **Options**: check **Double Precision** (recommended for this problem — the parallel-plate geometry's extreme aspect ratio and the sharp MCE source-term switching benefit from the added numerical precision).
4. **Processing**: select **Serial** for an initial test run (faster to set up and debug), or **Parallel (Local Machine)** with your available core count for the full multi-cycle production run (Step 19) — the transient, UDF-heavy solve benefits significantly from parallelization once you're past initial setup/debugging.
5. Click **Start**.

**Expected result before moving on:** the Fluent solver GUI opens, with the mesh from Step 6 already loaded (since launched from the Workbench A4 cell) and the correct 2D/3D dimension shown in the Fluent window's title bar and the **General** task page.

**How to verify:** in the **General** task page (left tree, under **Setup**), confirm the reported mesh statistics (cell/face/node counts) are nonzero and roughly match what you saw in Fluent Meshing's Step 6 check.

**Common mistake:** launching Fluent with **Single Precision** — for a problem with a 0.1 mm gap and a bed spanning tens of centimeters, or temperature differences on the order of 0.01 K (the periodic-convergence tolerance, Step 19), single precision can introduce meaningful round-off error; always use double precision for this model.

---

## 9. Solver settings and the energy equation

**Goal:** configure Fluent's General solver settings (transient, pressure-based) and enable the Energy equation, which is required for any heat-transfer physics including the MCE source term.

### 9.1 General settings

1. **Setup → General**:
   - **Solver → Type**: **Pressure-Based** (standard choice for this low-speed, incompressible-liquid-water flow — not Density-Based, which targets compressible high-speed flows).
   - **Solver → Time**: **Transient** (required — the AMR cycle is inherently time-dependent, with phases switching every `t_mag` seconds).
   - **Velocity Formulation**: **Absolute** (default; no rotating reference frame needed for this stationary unit-cell model).
   - **2D Space** (2-D tracks only): if modeling the packed-bed track as an axisymmetric slice (matching the COMSOL guide's 2-D axisymmetric approach), select **Axisymmetric** here — this is Fluent's equivalent of COMSOL's "2D Axisymmetric" dimension choice, applied as a solver setting rather than a separate geometry type (as noted in Step 5.2). For the parallel-plate 2-D track, leave this as **Planar**.

**Expected result before moving on:** the General task page shows Pressure-Based, Transient, and (for the packed-bed 2-D track) Axisymmetric selected, with no warning icons.

**How to verify:** click **Check** (bottom of the General page, if present in your Fluent version) to run Fluent's built-in mesh/setup consistency check.

### 9.2 Enable the Energy equation

2. **Setup → Models → Energy**: double-click, or select and click **Edit...**, then check the box to **turn Energy on**.

**Why this step is needed:** without the Energy equation enabled, Fluent solves only the flow field (momentum/continuity) and has no temperature field at all — the MCE source term, the reservoir temperature boundary conditions, and Qc itself all require Energy to be on.

**Expected result:** the Models tree shows **Energy: On**.

**How to verify:** the **Models** task page's summary table lists "Energy" with status "On."

**Common mistake:** forgetting to enable Energy and then being confused why temperature-related boundary condition fields (Step 15) don't appear in the Boundary Conditions task page — Fluent only exposes thermal BC fields once Energy is on.

### 9.3 Laminar viscous model

3. **Setup → Models → Viscous**: confirm **Laminar** is selected (not k-epsilon or another turbulence model) — the flow velocities implied by `mdot = 0.084666 kg/s` through this bed's cross-section, combined with the correlation validity ranges (`Re < 5e5` packed bed, `Re < 2300` parallel plate — both dominated by laminar-regime friction-factor correlations per the original guide's §2.1), indicate a laminar-flow regime is the correct and consistent choice, matching the correlations' own stated Reynolds-number validity.

**Common mistake:** leaving Fluent's default **k-epsilon (2 eqn)** turbulence model selected without checking — this would silently apply a wall-function/turbulent-viscosity treatment inconsistent with the laminar friction-factor correlations §10 implements, producing a pressure drop and heat transfer that don't match the correlations you carefully fit in the next step.

---

## 10. Porous zone setup (packed-bed track)

**Goal:** configure the bed as a Fluent **Porous Zone**, with viscous (`1/α`) and inertial (`C2`) resistance coefficients derived from the same friction-factor correlation `core/thermal.py` uses — not Fluent's own Ergun-equation porous-media defaults.

**Why this step is needed:** per the original guide's §2.1, *"do not accept Fluent's own Ergun-equation porous-media defaults without checking they match."*

### 10.1 Correlations to implement (reproduced unchanged from the original guide)

| Geometry | Correlation | Validity | Fluent form |
|---|---|---|---|
| Packed bed | `f = 23.462 · Re^-0.6716` | `10 < Re < 5e5` | Convert to `ΔP = (μ/α)·v + C2·(1/2)ρ·v²` by fitting `α` and `C2` against `f(Re)` over the Reynolds range the operating point spans (§20) — a single-point linear/quadratic fit is not adequate if Re varies significantly across the cycle (flow is zero during magnetize/demagnetize, so this mostly matters within the flow-on steps) |
| Parallel plate | `f = 24/Re` | `Re < 2300` | Pure viscous drag, so `C2 ≈ 0` and `α` follows directly from Poiseuille flow between plates of spacing `2·d_h⁻¹` (see the hydraulic-diameter definition and Eq. 7's slot limit, §5.1) |

Hydraulic diameter, porosity, specific surface area: identical definitions to the COMSOL guide's §2.1/§3 tables (already inlined above in this tutorial's §5.1) — reuse those values rather than re-deriving them.

### 10.2 Fit `α` and `C2` against the correlation

*Software operation (not specified by the guide): the guide specifies the target f(Re) correlation and the general Fluent ΔP form but does not give the exact numeric α/C2 fit — perform this fit outside Fluent (e.g. in Python or a spreadsheet) before entering the two resulting numbers into the Cell Zone Conditions panel.*

1. Over the Reynolds-number range the Step 20 operating point's flow phases actually span (compute `Re = ρ_f u_s d_p / μ_f` using the calibrated `mdot=0.084666 kg/s`, the bed cross-section `0.002 m²`, porosity `0.365`, and water properties from §11, to get the superficial velocity `u_s`), sample `f(Re)` from the packed-bed correlation at several points across that range.
2. Convert each `f(Re)` to a pressure-drop-per-length `ΔP/L = f · (ρ_f u_s²)/(2 d_h) · (1/ε)` (Darcy-Weisbach form, consistent with the COMSOL tutorial's §7.3 derivation).
3. Fit `ΔP/L = (μ_f/α)·u_s + C2·(1/2)ρ_f·u_s²` (Fluent's standard porous-media form) against these sampled `(u_s, ΔP/L)` points via least-squares (linear regression on `u_s` and `u_s²` terms) to obtain `α` (permeability, m²) and `C2` (inertial resistance factor, 1/m).

### 10.3 Enter the Porous Zone settings in Fluent

4. **Setup → Cell Zone Conditions**, select the `bed_porous` zone (named selection from Step 7).
5. Check **Porous Zone**. Click **Edit...** to open the Porous Zone dialog.
6. **Viscous Resistance** tab: enter `1/α` (the inverse permeability from Step 10.2's fit) in the **Direction-1 Viscous Resistance** field (and Direction-2/3 if the porous zone is anisotropic — for this bed-averaged, isotropic model, use the same value in all directions unless you have a specific reason to differentiate).
7. **Inertial Resistance** tab: enter `C2` from the fit.
8. **Porosity** tab: enter **Porosity**: `0.365`.
9. **Fluid Porosity / Solid Material** settings: under this same dialog (Fluent's Porous Zone panel typically includes a **Solid** sub-tab when Energy is on), set the **Solid Material** to the Gadolinium material (Step 11) and confirm the porosity value here matches Step 10.3.8 — Fluent, like COMSOL, requires the porosity to be entered consistently everywhere it's referenced; a mismatch here silently produces wrong thermal mass, exactly as flagged in the COMSOL tutorial's §6.5.

**Expected result before moving on:** the `bed_porous` Cell Zone Condition shows **Porous Zone: enabled**, with viscous/inertial resistance and porosity values matching the fit from Step 10.2 and the 0.365 value from the original guide.

**How to verify:** after a first trial solve (Steps 17–19), extract the pressure drop across the bed at the operating-point flow rate and hand-check it against the `f(Re)` correlation, exactly as the COMSOL tutorial's §7.3 verification step does.

**Common mistake:** leaving Fluent's Porous Zone dialog on its **Power-Law Model** default (which fits `ΔP = C0·v^C1`, Fluent's own generic form, not tied to the Tušek et al. correlation) instead of explicitly entering the α/C2 values fit in Step 10.2 — always use the **Viscous Resistance / Inertial Resistance** (Darcy-Forchheimer) tabs with your own fitted values, not Fluent's built-in defaults or alternate porous-media models.


---

## 11. Material creation

**Goal:** define Gadolinium (solid) and water (fluid) materials in Fluent with the exact property values from the original guide's §4, identical to the COMSOL guide's table.

### 11.1 Material property values (reproduced unchanged from the original guide's §4)

| Property | Value | Source |
|---|---|---|
| Gd solid density | 7900 kg/m³ | `core/thermal.py::RHO_GD` |
| Gd solid specific heat (lattice, near-room-T representative) | 236 J/(kg·K) | `core/thermal.py::CP_SOLID_GD` — **use the full `total_heat_capacity(T)` table** (via UDF lookup, Step 13) instead of this flat constant, since the model spans temperatures near Tc=294K |
| Gd Curie temperature Tc | 294.0 K | `core/mce_material.py::GADOLINIUM` |
| Gd Debye temperature θ_D | 169.0 K | same |
| Gd molar mass | 157.25 g/mol | same |
| Water density | 997 kg/m³ | `core/thermal.py::water_properties()` |
| Water specific heat | 4186 J/(kg·K) | same |
| Water dynamic viscosity | 8.9×10⁻⁴ Pa·s | same |
| Water thermal conductivity | 0.606 W/(m·K) | same |

> ⚠️ **Preserve this note:** "prefer Fluent's built-in temperature-dependent property database over the constant values above if the model spans a temperature range wide enough for it to matter" (same caveat as the COMSOL guide's §4) — if you upgrade water to Fluent's built-in temperature-dependent database, flag this explicitly as a genuine improvement over the 0-D model, not a discrepancy to explain away.

### 11.2 Create the Gadolinium material

1. **Setup → Materials**, right-click **Solid → Create/Edit...**
2. **Name**: `gadolinium-gd`.
3. **Properties**:
   - **Density**: `7900` kg/m³.
   - **Cp (Specific Heat)**: change the dropdown from **Constant** to **User-Defined** (or **polynomial**, but User-Defined Function is the correct choice here since Step 13's UDF will supply `Ctotal_lookup(T)` directly, matching the COMSOL guide's identical treatment) — select the compiled UDF property function once it's built (Step 13 covers compiling and hooking this).
   - **Thermal Conductivity**: the original guide doesn't give an explicit Gd thermal-conductivity value. *Software operation (not specified by the guide):* use a literature value (~10.5 W/(m·K), polycrystalline Gd near room temperature) and flag this explicitly as **not sourced from the original guide**, identical to the COMSOL tutorial's Step 6.2 note — keep both tutorials' values consistent if you built both.
4. Click **Change/Create**.

### 11.3 Create the water material

5. **Materials**, right-click **Fluid → Create/Edit...**
6. **Name**: `water-liquid-constant`.
7. **Properties**:
   - **Density**: `997` kg/m³.
   - **Cp**: `4186` J/(kg·K).
   - **Thermal Conductivity**: `0.606` W/(m·K).
   - **Viscosity**: `8.9e-4` Pa·s (Fluent's default viscosity unit is kg/(m·s), numerically identical to Pa·s — confirm the units dropdown next to the field reads `kg/(m-s)` and enter the same numeric value).
8. Click **Change/Create**.

**Expected result before moving on:** both materials appear in the **Materials** tree with the values above, and are available to assign in Cell Zone Conditions.

**How to verify:** re-open **Create/Edit Material** for each and confirm the entered values persisted; check units dropdowns next to each field match the table above (Fluent lets you change display units per-field, which can silently produce a wrong numeric entry if you don't check).

**Common mistake:** leaving Cp on Fluent's default **Constant** setting for the Gd material — the original guide requires the temperature-dependent `Ctotal_lookup(T)` because the operating point spans temperatures near Tc=294K; a flat 236 J/(kg·K) value here would silently drop the λ-anomaly physics the source repo's own material model captures.

---

## 12. Assigning materials to cell zones

*Software operation (not specified by the guide) — standard Fluent cell-zone material assignment, following naturally from Steps 7/10/11.*

1. **Setup → Cell Zone Conditions**:
   - **Packed-bed track**: select `bed_porous`. Under **Material Name**, set the fluid-phase material to `water-liquid-constant`; under the Porous Zone dialog's **Solid** tab (Step 10.3.9), set to `gadolinium-gd`.
   - **Parallel-plate track**: select `fluid_channel` → **Material Name**: `water-liquid-constant` (Type: **Fluid**). Select `solid_plate` → **Material Name**: `gadolinium-gd` (Type: **Solid**).

**Expected result before moving on:** each cell zone's **Material Name** dropdown shows the correct material, with **Type** correctly set to Fluid or Solid.

**How to verify:** **Setup → Cell Zone Conditions** summary table lists both zones with the correct Type/Material columns filled in — no zones left on Fluent's default `air`/`aluminum` placeholders.

**Common mistake:** leaving a cell zone on Fluent's default material (often `air` for fluid zones, `aluminum` for solid) after only assigning the *other* zone — always check both zones explicitly, since Fluent does not warn you about an unassigned/default zone.

---

## 13. UDF compilation and hooking DEFINE_SOURCE

**Goal:** implement the MCE volumetric heat source as a compiled `DEFINE_SOURCE` UDF, exactly reproducing:

```
q_MCE(T) = ρ_solid · C_total(T) · ΔT_ad(T, μ0H_max) / t_mag
```

with the sign flipping between magnetize (+1) and demagnetize (−1), and zero during the two flow phases — hooked to the Energy equation in the solid (or porous-solid) zone.

### 13.1 The UDF source code (reproduced unchanged from the original guide's §2.2)

```c
#include "udf.h"

/* Lookup table generated by core/mce_material.py -- see Step 4 above.
   Loaded at UDF init from the same CSV the COMSOL guide's Step 5 script
   produces (gd_dTad_vs_T_1p13T.csv), NOT re-derived inside the UDF. */
extern real dTad_lookup(real T);
extern real Ctotal_lookup(real T);

DEFINE_SOURCE(mce_source, cell, thread, dS, eqn)
{
    real T = C_T(cell, thread);
    real rho_solid = 7900.0;          /* core/thermal.py::RHO_GD */
    real q = 0.0;

    /* cycle_phase and t_mag are set by DEFINE_EXECUTE_AT_END or a
       global UDM (User-Defined Memory) driven by a separate cycle-
       timing UDF -- see Step 14. Sign flips between magnetize (+1) and
       demagnetize (-1); zero during the two flow phases. */
    if (cycle_phase == MAGNETIZE)
        q = rho_solid * Ctotal_lookup(T) * dTad_lookup(T) / t_mag;
    else if (cycle_phase == DEMAGNETIZE)
        q = -rho_solid * Ctotal_lookup(T) * dTad_lookup(T) / t_mag;

    dS[eqn] = 0.0;   /* source term has no explicit T-dependence needed
                         for a converged linearization here; revisit if
                         convergence is poor near Tc where dTad_lookup
                         has its sharpest curvature */
    return q;
}
```

`dTad_lookup`/`Ctotal_lookup` should read the same CSV table generated in Step 4 (from `core/mce_material.py::MagnetocaloricMaterial.delta_T_adiabatic()` / `.total_heat_capacity()`), via a simple piecewise-linear interpolation over the tabulated `(T, ΔT_ad, C_total)` triples — loaded once at `DEFINE_INIT` via Fluent's UDF file I/O, **not hardcoded as a polynomial fit**. This is the same requirement as the COMSOL guide's lookup-table treatment: the point is that both FEA tools' source term is provably tied to this repo's own validated material model, not a re-derivation that could silently diverge from it.

### 13.2 Implement the CSV-reading lookup functions

*Software operation (not specified by the guide) — the original guide's code block references `dTad_lookup`/`Ctotal_lookup` and specifies they should be loaded via `DEFINE_INIT` file I/O with piecewise-linear interpolation, but doesn't give the actual C implementation; a standard implementation is provided below, consistent with those stated requirements.*

```c
#include "udf.h"

#define MAX_ROWS 401
static real T_table[MAX_ROWS];
static real dTad_table[MAX_ROWS];
static real Ctotal_table[MAX_ROWS];
static int n_rows = 0;

DEFINE_INIT(load_mce_table, domain)
{
    FILE *fp;
    char line[256];
    int i = 0;

    /* Hard-code the exact path to gd_dTad_vs_T_1p13T.csv on the machine
       running this case -- see Step 4.2's note on keeping this file
       accessible. */
    fp = fopen("gd_dTad_vs_T_1p13T.csv", "r");
    if (fp == NULL)
    {
        Message("ERROR: could not open gd_dTad_vs_T_1p13T.csv\n");
        return;
    }
    fgets(line, sizeof(line), fp);  /* skip header row */
    while (fgets(line, sizeof(line), fp) != NULL && i < MAX_ROWS)
    {
        sscanf(line, "%lf,%lf,%lf", &T_table[i], &dTad_table[i], &Ctotal_table[i]);
        i++;
    }
    n_rows = i;
    fclose(fp);
    Message("Loaded %d rows from gd_dTad_vs_T_1p13T.csv\n", n_rows);
}

real dTad_lookup(real T)
{
    int i;
    if (T <= T_table[0]) return dTad_table[0];
    if (T >= T_table[n_rows - 1]) return dTad_table[n_rows - 1];
    for (i = 0; i < n_rows - 1; i++)
    {
        if (T >= T_table[i] && T <= T_table[i + 1])
        {
            real frac = (T - T_table[i]) / (T_table[i + 1] - T_table[i]);
            return dTad_table[i] + frac * (dTad_table[i + 1] - dTad_table[i]);
        }
    }
    return 0.0;
}

real Ctotal_lookup(real T)
{
    int i;
    if (T <= T_table[0]) return Ctotal_table[0];
    if (T >= T_table[n_rows - 1]) return Ctotal_table[n_rows - 1];
    for (i = 0; i < n_rows - 1; i++)
    {
        if (T >= T_table[i] && T <= T_table[i + 1])
        {
            real frac = (T - T_table[i]) / (T_table[i + 1] - T_table[i]);
            return Ctotal_table[i] + frac * (Ctotal_table[i + 1] - Ctotal_table[i]);
        }
    }
    return 0.0;
}
```

**Common mistake:** using a relative path (`"gd_dTad_vs_T_1p13T.csv"`) that resolves differently depending on which directory Fluent was launched from — prefer a full absolute path in the `fopen()` call, or confirm (every time you launch this case) that Fluent's working directory matches where the CSV lives (`File → Show Working Directory`, or check the console's startup banner).

### 13.3 Compile the UDF

1. **User-Defined → Functions → Compiled...**
2. Click **Add...**, browse to and add both `.c` source files (the `DEFINE_SOURCE` file and the lookup-table file — or combine them into one `.c` file, which is simpler and avoids needing to manage extern linkage across compilation units).
3. Click **Build**. Fluent invokes your registered C compiler (Step 2 prerequisite) to compile a shared library.
4. Once built successfully, click **Load** to load the compiled library into the running Fluent session.

**Expected result before moving on:** the Fluent console reports a successful build with no compiler errors, and `mce_source` becomes selectable as a UDF option in later dropdown menus (Step 13.4).

**How to verify:** **User-Defined → Functions → Compiled...** dialog should list `libudf` (or your chosen library name) with status "Loaded."

**Common mistake:** editing the `.c` file after loading and forgetting to re-**Build** and re-**Load** — Fluent does not auto-detect source file changes; you must explicitly rebuild.

### 13.4 Hook DEFINE_SOURCE to the Energy equation

5. **Setup → Cell Zone Conditions**, select the solid domain (`bed_porous`'s Solid sub-tab for the packed-bed track, or `solid_plate` directly for the parallel-plate track).
6. Under the zone's **Source Terms** tab, check **Source Terms** enabled, then under **Energy**, click **Edit...** and add a source term. From the dropdown, select **udf mce_source** (or however your compiled function is named/exposed in the dropdown — Fluent lists compiled `DEFINE_SOURCE` functions by their registered name).
7. Click **OK**.

**Expected result before moving on:** the solid zone's Cell Zone Conditions summary shows a non-default Energy source term referencing `mce_source`.

**How to verify:** run a short trial solve (a few time steps, before the full Step 19 production run) and monitor a solid-zone temperature point — with the cycle-phase UDM hard-coded temporarily to force `MAGNETIZE` (before Step 14's full state machine is wired in), confirm the temperature visibly rises, consistent with `q_MCE > 0`.

---

## 14. Cycle state machine (Fluent Named Expressions / UDF-driven, implementing the AMR cycle)

**Goal:** implement the same four-step cycle as the COMSOL guide's Events interface, but as an explicit UDF-driven state machine, since Fluent has no direct equivalent to COMSOL's Events interface.

**Why this step is needed:** per the original guide's §5, Fluent needs an explicit mechanism to (a) track elapsed time against the cycle period and blow fraction, (b) update a phase flag, and (c) gate both the `DEFINE_SOURCE` UDF (Step 13) and the inlet boundary condition's mass-flow profile off that flag.

### 14.1 The four-step cycle (reproduced unchanged from the original guide's §5)

1. **Adiabatic magnetization** (`H: 0 → H_max`): MCE source term active, fluid flow OFF.
2. **Cold-to-hot flow**: MCE source term OFF, fluid flow ON (cold→hot direction), duration set by the blow fraction (`BLOW_FRACTION_MASCHE`, default 0.5).
3. **Adiabatic demagnetization** (`H: H_max → 0`): MCE source term active with sign flipped, fluid flow OFF.
4. **Hot-to-cold flow**: MCE source term OFF, fluid flow ON (hot→cold direction).

With `frequency = 0.75 Hz` (period = `1.333 s`) and blow fraction 0.5 applied symmetrically: each phase gets `1.333/4 = 0.333 s` (identical arithmetic to the COMSOL tutorial's Step 8.1 — `t_mag = 0.333 s`).

### 14.2 Register a User-Defined Memory (UDM) location for `cycle_phase`

1. **User-Defined → Memory...** (or **Parameters and Customization → User Defined → Memory** depending on Fluent version). Set **Number of User-Defined Memory Locations**: at least `1` (UDM 0 will hold the integer phase flag; if you also want a separate UDM slot for a real-valued `cycle_sign`, request `2`).

### 14.3 The cycle-timing UDF (reproduced/expanded from the original guide's §5 description)

```c
#include "udf.h"

#define MAGNETIZE 1
#define FLOW_C2H  2
#define DEMAGNETIZE 3
#define FLOW_H2C  4

real t_mag = 0.333;      /* seconds, = period/4 at blow fraction 0.5 */
real period = 1.333;     /* seconds, = 1/frequency, frequency = 0.75 Hz */
int cycle_phase = MAGNETIZE;   /* global, read by mce_source (Step 13) */
real cycle_time = 0.0;         /* time elapsed within the current period */

DEFINE_EXECUTE_AT_END(update_cycle_phase)
{
    real t = CURRENT_TIME;
    cycle_time = fmod(t, period);

    if (cycle_time < t_mag)
        cycle_phase = MAGNETIZE;
    else if (cycle_time < 2.0 * t_mag)
        cycle_phase = FLOW_C2H;
    else if (cycle_time < 3.0 * t_mag)
        cycle_phase = DEMAGNETIZE;
    else
        cycle_phase = FLOW_H2C;

    Message("t=%f  cycle_time=%f  phase=%d\n", t, cycle_time, cycle_phase);
}
```

Compile this alongside the Step 13 UDF files (all in one `libudf` library, since `cycle_phase` needs to be visible to both `mce_source` and the boundary-condition profile UDF in Step 14.4 — declare `extern int cycle_phase;` in whichever `.c` file doesn't define it directly, or place them all in one file to avoid extern-linkage issues).

**Common mistake:** using `CURRENT_TIME` incorrectly across a restarted/resumed solve — `CURRENT_TIME` is the *absolute* simulation time, not time-since-restart; the `fmod(t, period)` in the code above already handles this correctly since it always wraps to within one period regardless of how many periods have already elapsed, but double-check this if you write your own variant.

### 14.4 Gate the inlet mass-flow boundary condition on `cycle_phase`

```c
#include "udf.h"

extern int cycle_phase;  /* defined in the cycle-timing UDF, Step 14.3 */

DEFINE_PROFILE(inlet_massflow_profile, thread, position)
{
    face_t f;
    real mdot = 0.084666;   /* kg/s, calibrated value -- Step 20 */
    real value = 0.0;

    if (cycle_phase == 2)        /* FLOW_C2H: cold-to-hot */
        value = mdot;
    else if (cycle_phase == 4)   /* FLOW_H2C: hot-to-cold */
        value = -mdot;
    else
        value = 0.0;             /* MAGNETIZE or DEMAGNETIZE: flow off */

    begin_f_loop(f, thread)
    {
        F_PROFILE(f, thread, position) = value;
    }
    end_f_loop(f, thread)
}
```

This single profile, applied at one boundary as a **mass-flow-inlet** with sign convention matched to the mesh's face normal, reproduces the direction reversal between phase 2 and phase 4 the same way the COMSOL tutorial's Step 10.5 boolean expression does — avoiding the need for two separate, manually-toggled Inlet/Outlet boundary pairs.

**Expected result before moving on:** compiling and loading this alongside the Step 13 UDFs (§13.3) succeeds with no errors; the compiled-UDF list in Fluent shows `update_cycle_phase` and `inlet_massflow_profile` as available hooks.

**How to verify:** run a short trial solve and check the Fluent console's `Message()` output from `update_cycle_phase` (Step 14.3) — confirm `phase` cycles 1→2→3→4→1... at the expected 0.333 s intervals.

**Common mistake:** forgetting `DEFINE_EXECUTE_AT_END` runs *after* each time step completes, not before — if the very first time step needs the correct initial phase (t=0 should be `MAGNETIZE`), also initialize `cycle_phase = MAGNETIZE` explicitly in `DEFINE_INIT` (Step 13.2's `load_mce_table` function is a convenient place to add this one extra line), so the first step doesn't run with an undefined or stale phase value.

---

## 15. Boundary conditions

**Goal:** set the cold-end, hot-end, wall, and inlet boundary conditions exactly as specified in the original guide's §6 (identical to the COMSOL guide's §6, translated into Fluent's boundary-condition types).

### 15.1 Boundary condition values (reproduced unchanged from the original guide's §6)

| Boundary | Fluent BC type | Active during | Value |
|---|---|---|---|
| Cold end | `pressure-outlet` (during phase 4, when it's the outflow-adjacent reservoir) or `velocity-inlet`/`mass-flow-inlet` depending on flow direction at that instant | Hot-to-cold flow (phase 4), when fluid draws heat from this reservoir into the bed | `T_cold = 289 K` |
| Hot end | `pressure-outlet`/`velocity-inlet` pairing, gated by the same phase flag | Cold-to-hot flow (phase 2) | `T_hot = T_cold + span = 299.2 K` |
| Bed/plate walls | `wall`, adiabatic | Always | matches the 0-D model, which has no wall-loss term at all |
| Fluid inlet | `mass-flow-inlet` with the Step 14.4 profile | Only during the two flow phases | `mdot = 0.084666 kg/s` |

> ⚠️ **Preserve this note:** adiabatic bed walls match the 0-D model's lack of any wall-loss term — flag explicitly if this is later relaxed with a real wall-conduction/radiative loss term.

### 15.2 Set up the inlet/outlet pair with reversing flow

Because the flow direction physically reverses between phase 2 and phase 4, the cleanest Fluent implementation (matching the single-profile approach in Step 14.4) is:

1. **Setup → Boundary Conditions**, select the `cold_end` face. Set **Type**: **mass-flow-inlet**.
2. Under **Momentum**, set **Mass Flow Rate**: select **udf inlet_massflow_profile** (the Step 14.4 hook) instead of a constant value — this makes the sign (and hence effective direction) switch automatically with `cycle_phase`.
3. Under **Thermal**, set the boundary's **Temperature** specification method to a value gated the same way — since a single `mass-flow-inlet` boundary needs a temperature for whichever direction is currently *inflow*: use a small UDF (`DEFINE_PROFILE`, similar structure to Step 14.4) that returns `T_cold` when `cycle_phase==4` (this face is the genuine cold-reservoir inflow) and is simply unused/ignored by Fluent when the mass flow at this face is actually outflow (phase 2), since Fluent only applies an inlet's specified temperature to inflowing mass.
4. Repeat symmetrically at `hot_end`: **mass-flow-inlet** with a **Temperature** UDF profile returning `T_hot` when `cycle_phase==2`.

*Software operation (not specified by the guide): the guide specifies the physical temperature/flow-direction requirements but not this exact dual-mass-flow-inlet Fluent implementation pattern — this is the standard way to implement a reversing-flow-with-reservoir-temperature boundary pair in Fluent without needing to dynamically change boundary condition *types* mid-solve (which Fluent does not support without a case-modification UDF, a much more complex alternative not needed here).*

5. During the magnetize/demagnetize phases (1 and 3), the `inlet_massflow_profile` UDF (Step 14.4) already returns `value = 0.0` at both ends, so no flow enters or leaves — Fluent naturally treats a zero-mass-flow inlet as a no-flow (effectively wall-like, for mass purposes) boundary for that period.

### 15.3 Set adiabatic walls

6. Select the lateral wall named selection (`bed_wall` or `plate_outer_wall`). **Type**: **wall**. Under **Thermal**, confirm **Heat Flux** = `0` (Fluent's wall default is already adiabatic/zero-heat-flux unless you specify otherwise — explicitly verify this rather than assuming).

**Expected result before moving on:** `cold_end` and `hot_end` are both `mass-flow-inlet` type with UDF-driven mass flow and temperature profiles; the lateral walls are `wall` type with zero heat flux.

**How to verify:** **Setup → Boundary Conditions** summary table lists all boundaries with the correct Type column; open each and confirm the UDF hooks (not constant values) are selected where specified.

**Common mistake:** setting the cold/hot end boundaries as plain **pressure-outlet** with a fixed temperature — a `pressure-outlet`'s specified temperature in Fluent only applies to *backflow* (reverse flow into the domain), which is a different physical meaning than what's needed here; the **mass-flow-inlet** with a phase-gated profile (as above) is the correct pattern for this reversing, reservoir-driven flow.

---

## 16. Initialization

*Software operation (not specified by the guide) — standard Fluent solution initialization; the guide requires "a bed starting at uniform T is not the AMR's actual operating point" once past cycle 1, but the very first cycle does need some starting condition.*

1. **Solution → Initialization**. Select **Standard Initialization**.
2. **Compute from**: `all-zones` (or select one boundary, e.g. `cold_end`, as the reference).
3. Set **Temperature**: an initial guess, e.g. the average of `T_cold` and `T_hot` = `(289+299.2)/2 = 294.1 K` — a reasonable midpoint starting condition, though (per the original guide's §5) this initial condition's specific value matters only for how many cycles it takes to reach periodic steady state, not for the final converged Qc.
4. Set **Gauge Pressure**: `0`.
5. Click **Initialize**.
6. Also initialize the UDM `cycle_phase` (Step 14) to `MAGNETIZE` explicitly, and confirm the `load_mce_table` `DEFINE_INIT` UDF (Step 13.2) has run (check the console for its `"Loaded %d rows..."` message) before proceeding to Step 19's production run.

**Expected result before moving on:** the Graphics window (if you plot a contour of Temperature at this point) shows the uniform initial field; the console confirms the UDF table load message.

**How to verify:** **Solution → Report Definitions** or a quick **Contours** plot of Temperature immediately after initialization should show a flat, uniform field at your chosen initial value — not an inherited stale field from a previous, unrelated case.

**Common mistake:** re-using a previous case/data file's field as the "initialization" for a new geometry/mesh without re-initializing — always explicitly click **Initialize** for a fresh production run, especially after any mesh or boundary-condition change.

---

## 17. Time-stepping and transient solver settings

*Software operation (not specified by the guide) — standard Fluent transient run configuration; the guide specifies the physics and cycle structure but not explicit time-step sizes/solver tolerances.*

1. **Solution → Run Calculation**.
2. **Time Step Size**: set to `0.001 s` — small enough to resolve within-phase dynamics and, critically, to land close to each phase transition (every `0.333 s`) without a large step straddling the discontinuity in the MCE source term or the flow-direction reversal (the same concern the COMSOL tutorial's Step 12.3 addresses with a capped Maximum step; Fluent's fixed time-step transient scheme needs this chosen directly since it has no automatic Events-respecting solver the way COMSOL does).
3. **Number of Time Steps**: for a first trial run, `2000` (covering `2000 × 0.001 = 2 s`, about 1.5 periods) — enough to sanity-check the phase-switching logic (Step 14) before committing to the full multi-cycle production run (Step 19, which needs many more steps).
4. **Max Iterations/Time Step**: `20` (a reasonable default for this pressure-based transient solve; increase if residuals aren't dropping sufficiently within each time step — see Step 18).
5. Under **Solution Methods**: confirm **Scheme**: **SIMPLE** (a robust, standard choice for pressure-velocity coupling in this type of problem); **Gradient**: **Least Squares Cell Based**; **Momentum** and **Energy** discretization: **Second Order Upwind** (better accuracy than First Order for resolving the sharp thermal gradients this problem produces, especially in the thin parallel-plate gap).
6. Under **Solution Controls**, leave under-relaxation factors at Fluent's defaults initially; reduce them (e.g. Pressure 0.3, Momentum 0.5) only if the trial run in Step 17.3 shows oscillating/diverging residuals.

**Expected result before moving on:** the Run Calculation page shows the settings above with no red warnings; clicking **Calculate** on the short trial run (Step 17.3's 2000 steps) completes without a reported divergence.

**How to verify:** watch the Fluent console/residual plot during the trial run (see Step 18) — confirm residuals drop at least 2-3 orders of magnitude within most time steps, and that the `cycle_phase` `Message()` output (Step 14.3) shows the expected 1→2→3→4 progression at the correct times.

**Common mistake:** choosing too large a time step (e.g. `0.01 s` or larger) relative to the `0.333 s` phase duration — while this wouldn't literally skip a phase transition (Fluent will still evaluate `cycle_phase` correctly at each `DEFINE_EXECUTE_AT_END` call), a coarse time step relative to the phase duration under-resolves the transient response *within* each phase, particularly the sharp thermal response right after a phase switch; keep the time step at least ~300x smaller than the phase duration as a starting point (as in the `0.001 s` value above, vs. `0.333 s` phases).

---

## 18. Residual monitoring

*Software operation (not specified by the guide) — standard Fluent residual-monitoring workflow.*

**Goal:** confirm the solver is actually converging within each time step, not just running without crashing.

1. **Solution → Monitors → Residuals**: confirm the default residual plots for **continuity**, **x-velocity**, **y-velocity** (and **z-velocity** for 3-D), and **energy** are enabled.
2. Set **Convergence Criterion** for each: `1e-3` for continuity/momentum (Fluent's common default), and a tighter `1e-6` for **energy** — energy residuals need to be much smaller than momentum residuals for a problem where the actual physics of interest (Qc) is a *thermal* quantity; a loosely-converged energy equation directly corrupts the very number you're trying to validate against.
3. Also add a **Report Definition** monitor tracking a representative point temperature (e.g., mid-bed) and the cold-end mass-weighted average temperature, plotted vs. time step, so you can visually track the periodic-steady-state approach in real time (feeding directly into Step 20's convergence check) rather than only inspecting after the full run completes.

**Expected result before moving on:** during the Step 17.3 trial run, the residual plot shows all monitored residuals dropping below their convergence criteria within each time step's inner iterations (not pegged at the iteration-count ceiling, Step 17.4's Max Iterations = 20, without actually converging).

**How to verify:** if residuals are still decreasing right up to iteration 20 at many time steps (i.e., hitting the iteration ceiling rather than the convergence criterion), increase **Max Iterations/Time Step** (Step 17.4) or reduce the time step size (Step 17.2) until residuals reliably converge within the iteration budget.

**Common mistake:** treating "the run finished without an error" as equivalent to "the solution converged" — always check the actual residual trend, not just that Fluent didn't crash; a transient run can complete all requested time steps while never actually converging at each step, silently producing an inaccurate time-integrated result.

---

## 19. Running the full transient simulation and periodic convergence

**Goal:** run enough full AMR periods for the bed's periodic-steady-state temperature profile to converge, per the original guide's §5/§8 requirement (identical convergence criterion to the COMSOL guide — cycle-to-cycle ΔT below a small tolerance, e.g. 0.01 K, at every monitored point) before reading off Qc.

### 19.1 Operating point (reproduced unchanged from the original guide's §7 — identical to the COMSOL guide's benchmark)

| Parameter | Value |
|---|---|
| Material | Gd (single-Tc approximation) |
| μ0H_max | 1.13 T |
| Regenerator mass | 1.7 kg |
| Frequency | 0.75 Hz |
| T_cold | 289 K |
| Span | 10.2 K |
| Calibrated mdot | 0.084666 kg/s |
| Target Qc | 102.8 W |
| Target COP | 3.1 (this repo's 0-D model: −2.1% error, `implied_parasitic=0.255`) |

Source: Eriksen, Engelbrecht, Bahl, Bjørk, Nielsen, Insinga, Pryds, *Int. J. Refrigeration* (2015), doi:10.1016/j.ijrefrig.2015.05.004.

### 19.2 Extend the run to many periods

1. Return to **Solution → Run Calculation**. With the Step 17-18 settings confirmed good from the trial run, set **Number of Time Steps** to cover at least 20 full periods: `20 × 1.333 s / 0.001 s = 26,660` time steps (extend if periodic convergence, checked in §19.3, isn't yet met).
2. Switch **Processing** to **Parallel** if not already (Step 8.4) — this run is computationally significant (tens of thousands of time steps, each with UDF calls and up to 20 inner iterations) and benefits substantially from multiple cores.
3. Click **Calculate**. Consider running this via Fluent's batch/journal-file mode for an unattended, resumable run on a workstation or cluster (`fluent 3ddp -g -i case_setup.jou` from the command line, where `case_setup.jou` is a Fluent journal script automating the read-case, initialize, and run-calculation steps) rather than the interactive GUI, given the run length.

**Expected result before moving on:** the solve completes the requested number of time steps (or you can interrupt periodically to check periodic convergence, per §19.3, and extend as needed) without a reported divergence.

### 19.3 Check periodic convergence

Per the original guide's explicit requirement (§5 and §8, identical to the COMSOL guide's Step 14): *"Run enough full periods for the bed's periodic-steady-state temperature profile to converge (same convergence criterion as the COMSOL guide — cycle-to-cycle ΔT below a small tolerance at every monitored point) before reading off Qc."*

4. Using the Step 18.3 Report Definition monitors (point/average temperatures vs. time), export the monitor history (**File → Export → Solution Data**, or the auto-generated `.out` monitor file) and compute, for each monitored point, the temperature at the *end* of each period `T(k×1.333 s)` for `k = 1, 2, ..., N_periods`, then `ΔT_k = |T(k×period) - T((k-1)×period)|`.

**Expected result:** `ΔT_k` at every monitored point drops below `0.01 K` by the final few periods.

**How to verify:** the exported monitor data's last several `ΔT_k` values should all read below `0.01`. If not, extend the run (Step 19.2's Number of Time Steps, resuming from the existing data file rather than restarting) until the tolerance is met — Fluent's **File → Read → Case & Data** followed by **Run Calculation** with additional time steps continues seamlessly from where the previous run left off.

**Common mistake:** checking periodic convergence only at one monitored point (e.g. mid-bed) — the original guide requires convergence "at every monitored point," and points nearer the reservoir boundaries (cold_end, hot_end) can converge at a different rate than the bed interior; monitor at least 3 points (cold end, mid-bed, hot end) as a minimum.

---

## 20. Post-processing

**Goal:** extract Qc from the converged periodic-steady-state solution.

*Software operation (not specified by the guide) — the guide specifies the target quantity (Qc) but not the exact Fluent post-processing clicks.*

### 20.1 Compute Qc (cooling capacity)

1. **Results → Reports → Surface Integrals** (or set up a **Report Definition** of type **Surface Heat Transfer Rate** during the run itself for a live-updating monitor, recommended so you don't have to reconstruct this from raw field data after the fact).
2. **Report Type**: **Heat Transfer Rate** (or "Total Heat Transfer Rate," depending on Fluent version's exact naming).
3. **Surface**: `cold_end`.
4. Restrict the time window to the **final converged period only** (the period identified as converged in Step 19.3) — if using a live **Report Definition** monitor, this means reading off the monitor's time-history export and averaging (or integrating and dividing by period) only over that final period's time range, matching exactly how the COMSOL tutorial's Step 15.1 restricts to the final converged period.
5. Additionally restrict the integration to the time sub-window when `cycle_phase == 4` (hot-to-cold flow, when the cold-end boundary is genuinely drawing heat from the reservoir into the bed) within that final period, consistent with the original guide's definition of Qc as tied specifically to that phase.
6. Time-average the resulting heat-transfer-rate values over that phase-4 sub-window to get a single Qc value in Watts.

### 20.2 COP

Per the same reasoning as the COMSOL tutorial's Step 15.2: the original guide's §7-8 comparison target is **Qc specifically**, not a Fluent-computed COP (the 0-D `implied_parasitic=0.255` figure comes from the Python-side `loss_model.py`). A full COP calculation would additionally require the magnetic work input (out of scope — a solved magnet-circuit field is explicitly excluded by §0/§1) and the pumping work (derivable from the Step 10 pressure-drop fit times volumetric flow rate, if desired, but not required for the Step 21 degeneracy check).

### 20.3 Visualize the temperature field

7. **Results → Graphics → Contours**, plot **Temperature** on the bed/channel domain at a specific point in the converged cycle (e.g., end of phase 2, the point of maximum hot-end temperature) — the genuinely new information this model adds over the 0-D model.
8. **Results → Animations → Solution Animation Playback**, if you saved intermediate data files during the Step 19.2 run (enable **Autosave** with a reasonable interval, e.g. every 50 time steps, *before* starting the production run if you want a smooth animation afterward) — animate one full converged period to visualize the demagnetization-front propagation and thermal wave motion.

**Expected result before moving on:** a single time-averaged Qc value (Watts) from the final converged period's phase-4 window, plus a temperature-field contour plot showing a genuine axial gradient (a flat/uniform result suggests non-convergence, Step 19, or a boundary-condition error, Step 15).

**How to verify:** the temperature contour should show T varying smoothly and monotonically (roughly) from `T_cold` near the cold end to `T_hot` near the hot end at the moment of peak flow — not flat (source term or flow not actually contributing) and not discontinuous/jagged (mesh under-resolution or a units error in the source term).

**Common mistake:** integrating Qc over the *whole* run including the early, non-periodic cycles instead of restricting to the final converged period's phase-4 window — this systematically biases the reported Qc, exactly the same pitfall flagged in the COMSOL tutorial's Step 15.

---

## 21. Validation against the benchmark (degeneracy check)

**Goal — reproduced unchanged from the original guide's §8:** *"Same standard as `COMSOL_SETUP_GUIDE.md` §8: reproduce Qc≈102.8W at span=10.2K before trusting any new geometry or gradient conclusion from the Fluent model."*

### 21.1 Compare

1. Take the Step 20.1 time-averaged Qc.
2. Compute percent error: `error = (Qc_Fluent - 102.8) / 102.8 × 100%`.

### 21.2 Interpret the result

> ⚠️ **Preserve this warning, unchanged:** *"If it doesn't land within roughly the 0-D model's own -2.1% error band (generously ±10-15% to allow for genuine spatial effects), the most likely explanations are a UDF unit error, a wrong blow-fraction assumption, or an unconverged periodic-steady state — check those before treating a mismatch as a physical finding."*
>
> **"If both the COMSOL and Fluent models are eventually built, cross-checking them against each other, in addition to this benchmark, is the strongest validation step available"** — two independently-coded solvers agreeing with each other and with the literature number is meaningfully stronger evidence than either one alone.

**How to verify pass/fail:** identical criterion to the COMSOL tutorial's Step 16.2 — if `|error| ≤ 10–15%`, proceed to draw conclusions about axial gradients, geometry optima, etc. (§21.3 below). If not, work back through: (a) UDF unit errors (Step 13's `rho_solid`, `mdot` hard-coded values — confirm every constant in the compiled C code matches the table in §19.1/§11.1 exactly, since a typo in compiled C fails silently, unlike COMSOL's unit-tagged fields which at least throw a units-mismatch warning), (b) the blow-fraction assumption (Step 14.1), (c) periodic-steady-state convergence (Step 19.3 — re-check the 0.01 K tolerance was genuinely met, not just assumed).

### 21.3 What a working Fluent model would add beyond the 0-D model (reproduced unchanged from the original guide's §9)

Real axial gradients within the bed; a direct check on whether `amr_cycle.py`'s linear `span_fraction` cutoff near the no-load span limit is a reasonable approximation (that function's own docstring already flags it as producing an artificially sharp corner vs. the gradual rolloff published Qc(span) curves show); whether the fixed-mdot geometry optimum `core/geometry_analysis.py` finds survives joint mdot+geometry optimization once gradients are resolved; and, if extended to a graded bed, whether inter-layer axial conduction (currently ignored by `cascade.py`'s independently-peak-tuned-layer treatment) matters.

---

## 22. Known limitations of this guide itself (reproduced unchanged from the original guide's §10)

- Never built or solved by the original guide's author (§0).
- MCE source term inherits the mean-field model's documented overprediction of ΔT_ad near Tc at low field (`core/validation.py`).
- Water properties as constants unless explicitly upgraded (§11.1's note).
- No wall-conduction/radiative loss term specified (§15).
- Uniform, prescribed magnetic field — no solved magnet-circuit model.
- **The UDF-based cycle state machine (Step 14) is a genuine additional implementation surface (compiled C code, not a built-in solver feature) compared to COMSOL's Events interface — more room for a units/logic bug to hide, which is exactly why the Step 21 degeneracy check matters more here, not less.**

---

## 23. Complete checklist

- [ ] Fluent license confirmed to support Compiled UDFs; C compiler installed and registered (§2)
- [ ] `gd_dTad_vs_T_1p13T.csv` generated from `core/mce_material.py::GADOLINIUM` at `mu0H=1.13T`, T=270–310K — the exact same file used by the COMSOL guide (Step 4.1)
- [ ] Workbench Fluid Flow (Fluent) system created per geometry/dimensionality track (Step 3)
- [ ] Geometry built in SpaceClaim: packed-bed (2D and/or 3D) and/or parallel-plate (2D and/or 3D), with Share Topology applied for parallel-plate conjugate interfaces (Step 5)
- [ ] Mesh generated via Fluent Meshing Watertight Geometry workflow, with fine local sizing on the parallel-plate gap; mesh independence checked (≤2% Qc change) (Step 6)
- [ ] Named selections created: `cold_end`, `hot_end`, wall(s), porous/fluid/solid zones, interface (Step 7)
- [ ] Fluent launched with correct 2D/3D dimension, Double Precision (Step 8)
- [ ] Pressure-Based, Transient, Axisymmetric (if applicable), Laminar, Energy equation ON (Step 9)
- [ ] Porous Zone configured (packed-bed track): viscous/inertial resistance fit against `f=23.462·Re^-0.6716` (packed) or `f=24/Re` (plate), porosity=0.365 entered consistently (Step 10)
- [ ] Materials created: `gadolinium-gd` (ρ=7900, Cp=UDF lookup, k~10.5 literature value) and `water-liquid-constant` (ρ=997, Cp=4186, μ=8.9e-4, k=0.606) (Step 11)
- [ ] Materials assigned to correct cell zones (Step 12)
- [ ] `DEFINE_SOURCE mce_source` UDF compiled and hooked to the solid-zone Energy source term (Step 13)
- [ ] Cycle-timing UDF (`update_cycle_phase`) and inlet mass-flow profile UDF compiled, loaded, and confirmed cycling correctly (Step 14)
- [ ] Boundary conditions set: cold/hot end as UDF-driven mass-flow-inlet pairs, adiabatic walls (Step 15)
- [ ] Solution initialized with a reasonable starting temperature; UDM/UDF state confirmed loaded (Step 16)
- [ ] Time step = 0.001 s, SIMPLE scheme, Second Order Upwind discretization confirmed (Step 17)
- [ ] Residual convergence criteria set (1e-3 momentum/continuity, 1e-6 energy); confirmed converging within iteration budget (Step 18)
- [ ] Ran ≥20 full periods (or extended until periodic convergence met) at the DTU_Eriksen_rotary_Gd_2015 operating point (Step 19)
- [ ] Cycle-to-cycle ΔT < 0.01 K confirmed at ≥3 monitored points before reading off Qc (Step 19.3)
- [ ] Qc extracted from the final converged period's phase-4 window only, at `cold_end` (Step 20)
- [ ] Degeneracy check performed: `|error| ≤ 10–15%` vs. Qc=102.8 W target (Step 21)

## 24. Parameter table (single reference, all values used in this tutorial)

| Parameter | Value | Unit | Where entered in Fluent |
|---|---|---|---|
| Applied field | 1.13 | T | Embedded in the Step 4 CSV generation; hard-coded reference in UDF comments |
| Regenerator mass | 1.7 | kg | Used to derive bed length (§5.1); not entered directly |
| Cycle frequency | 0.75 | Hz | UDF `period = 1.333` (=1/frequency), Step 14.3 |
| Cold-side temperature | 289 | K | UDF temperature-profile hook, Step 15.2 |
| Temperature span | 10.2 | K | Used to derive T_hot |
| Hot-side temperature (derived) | 299.2 | K | UDF temperature-profile hook, Step 15.2 |
| Calibrated mass flow rate | 0.084666 | kg/s | UDF `mdot`, Step 14.4 |
| Particle diameter (packed bed) | 0.5 (default; sweep 0.05–2) | mm | Geometry sketch dimension, SpaceClaim |
| Porosity (packed bed) | 0.365 | — | Porous Zone dialog, Step 10.3 |
| Bed cross-section | 0.002 | m² | Used to derive geometry radius/footprint (§5.1) |
| Derived bed length | 0.170 | m | Sketch/extrude dimension |
| Plate spacing (parallel plate) | 0.1–0.2 (0.1 optimum used) | mm | Fluid-channel rectangle height |
| Plate thickness (parallel plate) | 0.5 | mm | Solid-plate rectangle height |
| Magnetization ramp duration | 0.333 (derived: period/4) | s | UDF `t_mag`, Step 14.3 |
| Cycle period | 1.333 (=1/frequency) | s | UDF `period`, Step 14.3 |
| Gd solid density | 7900 | kg/m³ | Material panel; UDF `rho_solid` constant, Step 11.2/13.1 |
| Gd lattice Cp (flat, use only away from Tc) | 236 | J/(kg·K) | Not used directly — superseded by `Ctotal_lookup(T)` UDF |
| Gd Curie temperature | 294.0 | K | Reference only (embedded in the imported lookup table) |
| Gd Debye temperature | 169.0 | K | Reference only (embedded in the imported lookup table) |
| Gd molar mass | 157.25 | g/mol | Reference only (embedded in the imported lookup table) |
| Water density | 997 | kg/m³ | Material panel |
| Water specific heat | 4186 | J/(kg·K) | Material panel |
| Water dynamic viscosity | 8.9×10⁻⁴ | Pa·s (kg/(m·s)) | Material panel |
| Water thermal conductivity | 0.606 | W/(m·K) | Material panel |
| Packed-bed friction factor | `23.462·Re^-0.6716` (10<Re<5e5) | — | Fit to Porous Zone α/C2, Step 10.2 |
| Parallel-plate friction factor | `24/Re` (Re<2300) | — | Fit to Porous Zone α/C2 (C2≈0), Step 10.2 |
| Time step size | 0.001 | s | Run Calculation, Step 17.2 |
| Residual criterion (momentum/continuity) | 1e-3 | — | Monitors → Residuals, Step 18 |
| Residual criterion (energy) | 1e-6 | — | Monitors → Residuals, Step 18 |
| Target Qc (degeneracy check) | 102.8 | W | Comparison target, Step 21 |
| Target COP | 3.1 | — | Reference only (0-D model comparison) |
| Acceptable Qc error band | ±10–15 (0-D model's own error: −2.1%) | % | Step 21 pass/fail criterion |

## 25. Boundary-condition summary

| Boundary | Fluent BC type | Active phase(s) | Value |
|---|---|---|---|
| Cold end (`cold_end`) | `mass-flow-inlet`, UDF-driven mass flow + temperature | Phase 4 (hot-to-cold flow) is genuine inflow | T_cold = 289 K, mdot = −0.084666 kg/s (sign convention per Step 14.4) |
| Hot end (`hot_end`) | `mass-flow-inlet`, UDF-driven mass flow + temperature | Phase 2 (cold-to-hot flow) is genuine inflow | T_hot = 299.2 K, mdot = +0.084666 kg/s |
| Bed/plate walls | `wall`, adiabatic | Always | Heat Flux = 0 |
| Solid zone (`bed_porous`'s solid tab, or `solid_plate`) | `DEFINE_SOURCE` Energy source | Phases 1 & 3 only | `q_MCE = ρ·C_total(T)·ΔT_ad(T)/t_mag`, sign flips between phase 1 (+) and phase 3 (−) |

## 26. Material-property table

| Material | Property | Value | Unit |
|---|---|---|---|
| Gadolinium (solid) | Density | 7900 | kg/m³ |
| Gadolinium (solid) | Specific heat | UDF `Ctotal_lookup(T)` (temperature-dependent; flat 236 valid only away from Tc) | J/(kg·K) |
| Gadolinium (solid) | Curie temperature | 294.0 | K |
| Gadolinium (solid) | Debye temperature | 169.0 | K |
| Gadolinium (solid) | Molar mass | 157.25 | g/mol |
| Gadolinium (solid) | Thermal conductivity | ~10.5 (literature value — **not from the original guide**) | W/(m·K) |
| Water | Density | 997 (or Fluent's built-in temperature-dependent database if upgraded) | kg/m³ |
| Water | Specific heat | 4186 (or built-in database) | J/(kg·K) |
| Water | Dynamic viscosity | 8.9×10⁻⁴ (or built-in database) | Pa·s |
| Water | Thermal conductivity | 0.606 (or built-in database) | W/(m·K) |

## 27. Solver settings summary

| Setting | Value |
|---|---|
| Solver type | Pressure-Based, Transient |
| Viscous model | Laminar |
| Energy equation | On |
| Axisymmetric | On (packed-bed 2-D track only) |
| Time step size | 0.001 s |
| Total simulated time | ≥ 20 × period = ≥ 26.66 s (extend if not converged per Step 19.3) |
| Max iterations per time step | 20 |
| Pressure-velocity coupling | SIMPLE |
| Momentum/Energy discretization | Second Order Upwind |
| Gradient scheme | Least Squares Cell Based |
| Residual criterion (momentum/continuity) | 1e-3 |
| Residual criterion (energy) | 1e-6 |
| Mesh independence tolerance | ≤2% Qc change between refinement levels |

## 28. Validation checklist

- [ ] Lookup table (`gd_dTad_vs_T_1p13T.csv`) is the exact same file used by the COMSOL model, not regenerated with different bounds
- [ ] Porous Zone α/C2 (or parallel-plate equivalent) fit matches Tušek et al. Eq. 5/6 exactly — not Fluent's default Power-Law or Ergun-style porous model
- [ ] Porosity ε=0.365 entered consistently in both the Porous Zone dialog's flow and solid-material tabs
- [ ] UDF constants (`rho_solid=7900`, `mdot=0.084666`, `t_mag=0.333`, `period=1.333`) match the operating-point table exactly — checked by re-reading the compiled `.c` source, since a UDF typo fails silently
- [ ] Ran ≥20 full periods; cycle-to-cycle ΔT < 0.01 K at ≥3 monitored points before reading Qc
- [ ] Qc extracted from the final converged period only, restricted to the phase-4 window at `cold_end`
- [ ] Degeneracy check performed against Qc = 102.8 W (span=10.2K, μ0H=1.13T, mass=1.7kg, f=0.75Hz, mdot=0.084666 kg/s)
- [ ] Result falls within ±10–15% (or ideally closer to the 0-D model's own −2.1% error) before trusting any new geometry/gradient conclusion
- [ ] If a COMSOL model was also built, cross-checked against it directly (the "strongest validation step available" per the original guide)
- [ ] Every deviation from the 0-D model's assumptions (wall-loss terms, temperature-dependent water, resolved magnet field) explicitly flagged in any writeup, not silently introduced

## 29. Troubleshooting

| Symptom | Likely cause | Fix |
|---|---|---|
| Qc comes out orders of magnitude too large or small | A hard-coded UDF constant (e.g. `mdot`, `rho_solid`) has a units or typo error — no COMSOL-style unit-tag safety net in raw C | Re-read every numeric literal in the compiled `.c` files against §24's parameter table |
| Bed temperature field is flat/uniform even after many periods | `mce_source` UDF not actually hooked to the solid zone's Energy source term, or `cycle_phase` never changing | Re-check Step 13.4's Source Terms hookup; confirm the `Message()` output from Step 14.3 shows phase actually cycling |
| UDF fails to compile | Missing/incompatible C compiler, or a syntax error in the `.c` file | Check `Help → About Fluent`'s compiler-compatibility table against your installed compiler version; review the Fluent console's exact compiler error line |
| Solver diverges partway through the run | Time step too large relative to phase duration, or under-relaxation factors too aggressive for the sharp MCE source-term switching | Reduce time step (Step 17.2); reduce under-relaxation factors (Solution Controls) |
| Cycle-to-cycle ΔT never drops below 0.01 K | Not enough periods run, or a genuine periodic limit-cycle with amplitude above tolerance | Extend the run (resume from the existing data file); if still not converging after e.g. 100 periods, re-check for a modeling error such as a boundary condition that never actually turns off |
| Pressure drop across the bed doesn't match the hand-calculated f(Re) value | Porous Zone still on Fluent's default Power-Law model instead of the fitted Viscous/Inertial Resistance values | Re-check Step 10.3 — confirm the Darcy-Forchheimer (Viscous/Inertial Resistance) tabs, not the default Power-Law tab, hold the fitted values |
| `fopen()` fails in the UDF (table not loaded, `dTad_lookup` returns 0 everywhere) | CSV path wrong relative to Fluent's working directory | Use an absolute path in the `fopen()` call, or verify Fluent's working directory via `File → Show Working Directory` |
| Named selection not visible in Boundary Conditions or Cell Zone Conditions list | Named selection created in SpaceClaim/Meshing wasn't correctly passed through to the solver session (common if geometry was re-imported after naming) | Re-check the Named Selections were created *after* final geometry edits, and confirm they still appear correctly in the **Setup → General → Mesh → Info** zone list once inside Fluent |

## 30. Common errors and fixes

| Error message / symptom | Meaning | Fix |
|---|---|---|
| "UDF Library Build Failed" | Compiler not found, or a syntax error in the `.c` source | Check compiler registration (Step 2 prerequisite); review the exact line number in the Fluent console's build log |
| "Divergence detected in AMG solver" | Numerical instability, often from too large a time step or too coarse a mesh at a sharp gradient (e.g. the parallel-plate gap) | Reduce time step (Step 17.2); refine mesh locally (Step 6's fine sizing on the gap); reduce under-relaxation |
| "Negative cell volume detected" | Meshing error, often from a poorly-formed geometry (overlapping or non-shared-topology bodies at the fluid/solid interface) | Re-check Step 5.6's Share Topology step was applied before meshing |
| Segmentation fault when running with the compiled UDF loaded | Buffer overrun in the lookup-table arrays (Step 13.2), e.g. `MAX_ROWS` smaller than the actual CSV row count, or `n_rows` referenced before `DEFINE_INIT` has run | Confirm `MAX_ROWS=401` matches the Step 4.1 CSV's 401 rows exactly; confirm `DEFINE_INIT` runs (via its console message) before any solve step that calls `dTad_lookup`/`Ctotal_lookup` |
| Boundary condition Temperature field greyed out / not editable | Energy equation not yet enabled | Re-check Step 9.2 — Energy must be On before thermal BC fields become editable |
| Named UDF function not appearing in a dropdown menu (e.g. Source Terms, Boundary Conditions) | UDF library not yet **Loaded** (only **Built**), or the function's `DEFINE_*` macro type doesn't match the dropdown's expected hook type | Re-check Step 13.3's Load step; confirm you're selecting from the correct dropdown category (e.g. `DEFINE_SOURCE` functions only appear in Source Terms dropdowns, `DEFINE_PROFILE` functions only in boundary-condition profile dropdowns) |