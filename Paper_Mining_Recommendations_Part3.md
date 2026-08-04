# Paper-Mining Pass, Part 3 — Remaining Papers + Reference Books

Went through everything not yet touched in Parts 1–2: the economics paper in
full, both "Development of a rotary..." papers (confirmed which is which),
the original 1997 discovery paper, the theory review, the solid-state
caloric cooling review, and both reference books. One of these is the
single richest source found across all three passes.

---

## 1. A ~25-device table, cleaner and larger than the one in Part 2 — several strong new candidates

**Source:** Greco, Aprea, Maiorino, Masselli, *"A review of the state of the
art of solid-state caloric cooling processes at room-temperature before
2019,"* Int. J. Refrigeration (2019), Table 2, "Magnetocaloric devices built
to date... introduced after the year 2010." This table lists essentially
every notable AMR prototype built 2009–2018, with Q̇_ref,max [W],
ΔT_span,max [K], magnet type/field, and MCM+HTF for each.

**Caveat up front:** this PDF has a diagonal "ACCEPTED MANUSCRIPT" watermark
that bleeds single letters/fragments ("US", "AN", "IP", "CR", "T", "M",
"ED", "PT", "CE", "AC") into the extracted text of several rows. I checked
each row I'm reporting below against the surrounding column structure to
make sure the watermark fragments aren't corrupting an actual number, but
this is text-layer extraction, not a rendered-image read like Part 2's
table — treat these as good-confidence, not verified-by-eye.

**Strongest new candidates**, all independent of anything currently in your
benchmark CSV:

| Lab / device | Year | Type | Q̇_max (W) | ΔT_span,max (K) | Field (T) | MCM |
|---|---|---|---|---|---|---|
| Cooltech Application Co., France | 2011 | Reciprocating | 150 | **38** | 1.27 | Gd/Gd-Tb, parallel-plates |
| Cooltech Application Co., France | 2013 | Rotary | 120 | **42** | 0.98 | Gd/Gd-Tb, packed-bed |
| Cooltech Application Co., France | 2014 | Rotary | 300 | 38 | 1.17 | Gd/Gd-Er/La-Fe-Si, parallel-plates |
| DTU + ENOVHEAT partners, "MagQueen" | 2018 | Rotary | 1500 (heat) | 25 | 1.6 | La(Fe,Mn,Si)13Hz spheres |
| Kobe Univ. + Railway Tech. Research Inst. et al. (Hirano et al., 2014) | 2014 | Rotary | 1400 | 21 | 1.5 | Gd/Gd-Y/Gd-Dy flat plates |
| TOBB Univ./ASELSAN/Ankara Univ. (Turkey) | 2018 | Rotary | 500 | 5 | (4-pole permanent) | packed bed |

The **Cooltech 2013 device (42 K span)** is worth calling out specifically —
that's a larger span than anything currently in your benchmark set (your
existing max-span devices top out well below that). If you want a
stress-test case for the model at the extreme end of the achievable-span
range, that's it. The **DTU MagQueen** device is directly relevant to your
existing `LAFESIH_FAMILY` cascade material (same La-Fe-Si-H chemistry,
just Mn-doped) and gives a real large-scale (1500 W-class) operating point
for that family, which your current benchmarks don't have — every
LaFeSiH-material benchmark you have is Astronautics' rotary device; this
would be an independent lab/design cross-check on the same material class.

**One that needs a flag, not a clean add:** the table lists Astronautics'
2014 device (Jacobs et al.) as Q̇_max=3042 W / ΔT_span,max=**18 K**. The
review's own body text (which I read directly, not just the table)
confirms the 3042 W zero-span figure and the 2502 W-at-11-K operating point
that are already in your CSV — but it does *not* explain where "18 K"
comes from; that number appears only in the table, not the surrounding
prose. It may be the device's true zero-capacity max span (a third,
more-extreme point beyond your current two), or it could be a table
transcription issue on the review authors' part. I'd treat it as "worth
checking against the original Jacobs et al. 2013/2014 IJR paper's figures
before using," not as a confirmed third data point — I didn't find it
independently confirmed anywhere else in this corpus.

**Lower-confidence rows** (Univ. of Victoria PM-I/PM-II, TU Denmark
2011/2012/2014 reciprocating and rotary devices, Univ. of Ljubljana 2012/
2018, Kobe 2011, Sanden Co. 2013, Univ. of Salerno/Naples 2014, Univ. of
Coruña 2013, Univ. of Zaragoza 2014, Wroclaw 2014, TU Delft 2016, KTH 2016,
Baotou 2016): all present in the table with plausible-looking numbers, but
I did not independently cross-check each one against a primary source in
this corpus. Several duplicate labs already represented in your CSV
(Ljubljana ≈ Tušek's group, already there via a different paper) — I'm
flagging the table's existence and its cleanest rows rather than
transcribing all ~25 entries, since several would need the primary papers
(mostly not in `Papers/`) to verify before being trusted as calibration
targets.

---

## 2. The original 1997 discovery paper has a ΔT_ad ratio number that's more conservative than the ~2× usually quoted — and it independently supports your existing correction factor

**Source:** Pecharsky & Gschneidner, *"Giant Magnetocaloric Effect in
Gd5Si2Ge2,"* Phys. Rev. Lett. 78, 4494 (1997) — already your cited source
for `GD5SI2GE2_FIRST_ORDER`'s Tc and J, but I hadn't previously pulled a
number out of its own text, only used it for the composition/Tc citation.

The paper states directly, comparing its own heat-capacity-derived
ΔT_ad curves for Gd5Si2Ge2 against pure Gd:

> The ΔT_ad values of Gd5Si2Ge2 are larger than the corresponding ΔT_ad
> values for Gd by about 30%, comparing the peak values, regardless of
> the temperature.

That's a **~1.3× peak ΔT_ad ratio**, not the ~2× enhancement often quoted
for the *entropy* change (ΔS_M). This is independently consistent with
what your own `first_order_mce.py` already documents from Giguère et
al. (1999) — that the indirect/Maxwell-relation ΔS_M route overstates the
*true* ΔT_ad enhancement, and that a correction factor is warranted. This
gives you a **second, independent primary source** (heat-capacity-based,
not pulse-field-thermometry-based like Giguère) making the same
qualitative point with a specific number (~1.3×) rather than just Giguère's
single 7 T cross-check point. If useful, `first_order_mce.py`'s
`dTad_correction` derivation could cite this as corroborating evidence,
or — more rigorously — you could check whether the model's own computed
peak-ΔT_ad(Gd5Si2Ge2)/peak-ΔT_ad(Gd) ratio, after applying the existing
Giguère-derived correction, lands near 1.3, as a second held-out check
using a completely different physical measurement method than the one the
correction was fit to.

---

## 3. Confirmed identity of the two "rotary refrigerator development" papers — one is design-only, no performance numbers

- `Development of a novel rotary magnetic refrigerator.pdf` = **Lozano et
  al. 2016**, IJR, DOI 10.1016/j.ijrefrig.2016.04.005 — this is the primary
  source already behind your `Lozano_POLO_UFSC_2016` CSV row. Confirmed,
  no new content beyond what's presumably already extracted.
- `Development of a rotary magnetic refrigerator.pdf` = **Tušek, Zupan,
  Šarlah, Prebil, Poredoš, 2010**, IJR 33, 294–300 — this is a mechanical/
  magnet-design paper for an *earlier* Ljubljana prototype, not a
  performance-testing paper. It reports zero Qc/span numbers; its content
  is entirely about the permanent-magnet assembly design and a "pros and
  cons" table of the mechanical build (shaft-seal leakage causing
  rotational friction/heat generation, large magnet-structure weight,
  assembly/disassembly complexity). Useful only as qualitative
  engineering-realism context (e.g., for a discussion of real mechanical
  parasitic losses that `amr_cycle.py`'s idealized cycle doesn't capture),
  not as a numeric validation source. Worth knowing so nobody expects
  Qc/span data from this specific file.

---

## 4. Economics paper (Bjørk et al. 2011) — fully already mined, one loose end confirmed closed

Checked the full text, not just the abstract already used in
`economics.py`. The paper's Fig. 9 (minimum system cost vs. operating
frequency, at fixed 20 K span / 100 W) is the only content not already
reflected in your `COST_MCM_PER_KG`/`COST_MAGNET_PER_KG` constants — it's a
qualitative trend ("increasing frequency reduces cost," no simple optimum
found over the frequency range studied) presented only as a figure, not a
digitizable text value. Given your ROADMAP.md already treats full BOM cost
as open pending real HX/pump/motor data, this figure wouldn't move that
needle even if digitized — it's the same magnet+MCM-only cost model you
already use, just re-plotted vs. frequency. Not worth the digitization
effort relative to what it would add.

---

## 5. Reference books: one is front-matter only, the other is untapped but costly

- **Kitanovski et al., *Magnetocaloric Energy Conversion* (2015)**: the
  file in this corpus is only 30 pages — title page, preface,
  acknowledgments, table of contents, and the first ~9 pages of Chapter 1
  (thermodynamic fundamentals: COP and exergy-efficiency definitions). It
  is **not** the full book; Chapters 4 (AMR performance), 7 (prototypes by
  country — which the TOC shows covers French/Canadian/Chinese/Brazilian/
  Slovenian devices), and 9 (costs) are listed in the TOC but their pages
  aren't included in this file. The exergy-efficiency definition given in
  the fragment matches the standard formula (η = (q_R − T_amb·Δs_R)/(q_H −
  q_R)) and doesn't add anything beyond what's already standard. If you
  want the prototype/cost chapters, you'd need a more complete copy of
  this book — worth knowing so nobody expects that data to be extractable
  from this file.
- **Tishin & Spichkin, *The Magnetocaloric Effect and its Applications*
  (2003)**: 486 pages, but **entirely scanned page images with no OCR text
  layer** — I confirmed zero pages return extractable text via
  `pdfplumber`. This is genuinely the deepest, most comprehensive
  materials-property compendium in the whole corpus (rare-earth and alloy
  MCE data tables, historical device reviews), but getting anything out of
  it means running OCR across ~486 pages, which I haven't done here given
  the effort-to-yield tradeoff versus the primary-source papers already
  mined. If you want me to OCR specific chapters (e.g., its materials
  property tables), point me at a page range or topic and I can run
  Tesseract on that subset rather than the whole book.

---

## Updated priority list (all three passes combined)

1. Blow-fraction asymmetry in `amr_cycle.py` (Part 1, §1).
2. (Mn,Fe)₂(P,Si) as a third material family (Part 1, §3).
3. Extend Gd validation to 7 T + check the ≈6 K/T Curie-shift prediction,
   using papers already in the repo (Part 2, §2–3) — cheapest win, zero
   new sourcing.
4. **Cooltech's 42 K-span device and DTU's MagQueen** (Part 3, §1) as new
   secondary-source stress-test/cross-check benchmark rows.
5. Cross-check the Gd5Si2Ge2 ΔT_ad correction factor against the
   Pecharsky & Gschneidner ~1.3× peak-ratio figure (Part 3, §2) as a second
   independent test of `dTad_correction`.
6. Flag the parallel-plate validation gap (Part 1, §2) — still open,
   unchanged.
7. If you want the Tishin & Spichkin book's data tables, tell me which
   topic/chapter and I'll OCR that section specifically rather than the
   whole 486 pages.
