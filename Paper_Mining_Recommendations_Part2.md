# Paper-Mining Pass, Part 2 — Deeper Search

Went back through the remaining unmined papers more carefully (tables,
in-text numeric callouts, secondary reviews' own citation tables) rather
than abstract-level skims. Four more findings, ordered by how directly
usable they are.

---

## 1. A comparative table of ~12 AMR prototypes — several usable, several not (source clearly a secondary citation)

**Source:** Kamran, Ahmad, Wang, *"Review on the developments of active
magnetic regenerator refrigerators – Evaluated by performance,"* Renewable
and Sustainable Energy Reviews 133 (2020) 110247 (`Papers/Reviews/`) —
its Table 2, "A comparative summary of few AMR prototypes." This table
didn't extract cleanly as text (merged/multi-value cells), so I rendered
the page as an image and read it directly rather than trust the raw
PDF-to-text output — worth knowing if anyone re-extracts this later, since
`pdftotext`/`pdfplumber` scramble the row alignment here.

Like the existing `Okamura_Hirano_2013` row (already in
`amr_experimental_benchmarks.csv`, itself sourced from a *review's* table
rather than the primary paper), these are **secondary-source citations** —
the review cites the primary papers ([69], [70], [71] etc.) but those
primary PDFs are not in this repo's `Papers/` folder. Same caveat the
existing Okamura row already carries should apply to any of these if added.

**Clean, addable rows** (unambiguous after re-reading the rendered image):

| Lab | MCM | V_reg (cm³) | n_reg | f (Hz) | ΔB (T) | Qc,max (W) | ΔT_max (K) |
|---|---|---|---|---|---|---|---|
| Chubu Electric/Toshiba, Japan [69] | Gd spheres | 484 | 2 | 0.167 | 4 | 100 | 26 |
| — same device, lower field | Gd spheres | 484 | 2 | 0.167 | 2 | 40 | 24 |
| Institute of Tech., Chubu, Japan [70] | GdDy spheres | 844 | 4 | 0.39–0.42 | 1.1 | 540 | 0.2 |
| Nanjing University, China [71] — GdDy | GdDy spheres | 200 | 2 | 0.25 | 1.4 | 40 → 0 | 5 → 25 |
| Riso Lab, Denmark [18,50,62,72–74] | Gd spheres | 23 | 24 | 2.25 | 1.24 | 200 | 18.9 |

Two are worth flagging specifically:
- **Chubu Electric/Toshiba** gives a genuine two-field-point curve at fixed
  geometry (4T→100W/26K, 2T→40W/24K) — a clean field-sensitivity check the
  current benchmark set doesn't have (every existing device is single-field).
- **Institute of Tech., Chubu** (540 W at 0.2 K span) is a near-zero-span
  max-capacity extreme, same shape as the existing `_zerospan`/`_maxcap`
  companion rows, and at 540 W it's a genuinely different operating regime
  (high mass 844 cm³ Gd/GdDy) than anything currently calibrated.
- **Riso Lab** entry here (23 cm³, **24 regenerators**, 2.25 Hz, 1.24 T) does
  *not* match the existing `Risoe_DTU_Gd_2011` row's parameters (0.1955 kg,
  1 Hz, 30 K span) — different regenerator count and frequency strongly
  suggest this is a *different* Risø/DTU paper or device generation, not a
  duplicate of the one already in the CSV. I did not merge them; flagging
  as a second, independent Risø data point pending identification of its
  primary source (refs [18,50,62,72–74] in the review, none of which are
  in this repo's `Papers/`).

**Ambiguous rows I chose not to extract as clean numbers** (Astronautics,
Nanjing's Gd/GdSiGe sub-rows, U. Victoria) genuinely have merged multi-device
cells in the source table even in the rendered image — presenting a single
number from these would be guessing at which sub-value belongs to which
citation. Left out rather than reported with false precision.

**Rows citing labs whose primary papers aren't in this repo at all**
(POLO plates/1.65T — an *earlier*, different POLO device than the
`Lozano_POLO_UFSC_2016` one already in the CSV; U. Salerno; G2E Grenoble;
U. Tokyo; Wroclaw; and Teyber et al.'s superconducting-magnet device at
3.3 T reaching a **100 K span** — noted only because it's a striking
outlier, not because it's usable) are context, not addable data, until
their source papers are obtained.

---

## 2. A second, unused Gd validation point already sitting in a paper you've already exploited

**Source:** Giguère et al., PRL 83, 2262 (1999) — already in this repo,
already used by `giguere_validation.py`, but only for its Gd5Si2Ge2
content. Re-reading the same paper's *pure-Gd* cross-check paragraph (used
there only as the authors' own sanity check before moving to the alloy)
gives numbers `mce_material.py`'s Gd validation hasn't used:

> At 5 T, ΔT_ad ≈ 10.5–11.5 K (this paper's high-purity Gd vs. AMES
> laboratory's independent measurement, agreeing within 1 K). At 7 T,
> ΔT_ad ≈ 12–13 K (industrial- and high-purity Gd respectively), "agrees
> well with that of Brown (14 K)."

`validation.py` currently checks Gd against Dan'kov et al. (1998) at 1 T,
2 T, and 5 T only (per ROADMAP.md Phase 8's own note: "1T error... 5T
improved to -7.5%"). This gives a **free 7 T extension point** from a paper
already in the repo, without digitizing anything — just reading past the
alloy section into the paper's own Gd baseline paragraph.

---

## 3. A precise, quotable, testable field-dependence number for Tc — from the paper you already cite as your primary Gd source

**Source:** Dan'kov, Tishin, Pecharsky, Gschneidner, *"Magnetic phase
transitions and the magnetothermal properties of gadolinium,"* Phys. Rev. B
57, 3478 (1998) — this is the actual PRB paper behind `validation.py`'s Gd
calibration (confirmed from its own header: same four authors, 1998,
PRB 57). It contains a number not currently exercised anywhere in the code:

> Above 2 T, the Curie-point transition temperature increases almost
> linearly with field at a rate of **≈6 K/T**, up to fields of 7.5 T.
> Zero-field T_C = 294 ± 1 K, confirmed by four independent techniques.

`mce_material.py`'s mean-field/Weiss model treats `Tc=294.0` as a fixed
input parameter, but the effective peak-of-ΔS(T) temperature at a given
field is an *emergent* property of the self-consistent M(T,H) solve, not
hardcoded — so this is a genuine, free, zero-fabrication check: does the
model's own emergent peak-shift-with-field come out near 6 K/T when swept
over 2–7.5 T? This is exactly the kind of implicit-prediction check
`validation.py` already does for ΔS and ΔT_ad; it just hasn't been pointed
at this specific number yet. If it doesn't reproduce ~6 K/T, that's a real
finding about the mean-field model's field-dependence, not just a
confirmation exercise.

---

## 4. Confirmed dead ends (so they don't get re-mined later)

- **`Papers/Data center cooling/...`** (Ebrahimi, Jones & Fleischer 2014,
  already in `Literature_Review.md`): checked specifically for VCC/CRAC
  baseline COP or PUE figures that could ground `baseline_cooling.py`.
  Found none — the COP numbers in this paper (0.4–0.9) are for
  *absorption* cooling systems running on captured waste heat, an entirely
  different technology from the vapor-compression baseline
  `baseline_cooling.py` models. Not usable for that module.
- **`Papers/Magnetocaloric effect and materials physics/Materials
  challenges for high performance magnetocaloric refrigeration
  devices.pdf`** (Smith et al., Adv. Energy Mater. 2012): checked for
  quantitative corrosion/fatigue/hysteresis-loss numbers that could feed a
  device-lifetime or maintenance-cost term in `economics.py`. Found only
  qualitative discussion (Gd corrosion is mitigated by automotive
  antifreeze per Engelbrecht et al.; no systematic corrosion-rate studies
  published) — no number to add, but it does independently corroborate
  the existing honesty flag in `economics.py`/ROADMAP.md that
  manufacturing/maintenance/auxiliary costs are a real, still-unquantified
  gap, from a materials-focused angle rather than the Bjørk group's
  costing-focused angle.
- **`Papers/Reviews/Current perspective in magnetocaloric materials
  research.pdf`** (Law et al., J. Appl. Phys. 133, 040903 (2023)): a
  forward-looking perspective piece (hysteresis/reversibility, mechanical
  stability, emerging high-entropy-alloy magnetocalorics). No tabulated
  ΔS/ΔT_ad numbers for a specific composition were found in the pages
  checked — it's a landscape piece, useful for the paper's own "future
  work" framing, not a validation source.

---

## Updated priority list (combining both passes)

1. Blow-fraction asymmetry in `amr_cycle.py` (Part 1, §1) — still the
   single largest reported effect size found across both passes.
2. (Mn,Fe)₂(P,Si) as a third material family (Part 1, §3).
3. **New, cheap addition: extend `validation.py`'s Gd checks to 7 T using
   Giguère et al.'s own Gd paragraph (§2 above) and check the ≈6 K/T
   Curie-shift prediction from Dan'kov et al. 1998 (§3 above)** — both are
   already-in-repo papers, zero new sourcing required, and both are
   genuine held-out checks rather than re-fits.
4. Chubu Electric/Toshiba's two-field-point row (§1) as a new secondary-source
   CSV entry, styled like the existing Okamura row, if a field-sensitivity
   check is wanted without waiting on the parallel-plate digitization.
5. Flag the parallel-plate validation gap (Part 1, §2) — unchanged from
   Part 1, still open.
