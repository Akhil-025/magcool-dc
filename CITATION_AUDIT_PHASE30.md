# Citation Audit — Phase 30

Direct primary-source verification pass against the PDFs actually present
in this project's `Papers/` corpus, focused on the numbers with the most
downstream weight: the three `CALIBRATION_POINTS_CORE` benchmark devices
in `core/loss_model.py` (every COP_electrical prediction in
`comparison_table.csv`/`pareto_front.csv` traces back to these), plus a
spot-check on the two hysteresis-loss citations touched in the previous
pass. Method: for each claim, the actual PDF text was extracted and
searched for the specific numeric values used in code, not re-typed from
memory.

## 1. DTU_Eriksen_rotary_Gd_2015 — CONFIRMED, exact match

Source in code: `core/loss_model.py` CALIBRATION_POINTS_CORE,
`("DTU_Eriksen_rotary_Gd_2015", 0.75, 1.13, 0.084666, 102.8, 26.18)`.

Checked against: `Papers/Design and experimental tests of a rotary active
magnetic regenerator prototype.pdf` (Eriksen, Engelbrecht, Bahl, Bjørk,
Nielsen, Insinga, Pryds, Int. J. Refrig., accepted manuscript,
DOI 10.1016/j.ijrefrig.2015.05.004).

Direct quote-adjacent text found in the PDF: "the best result is obtained
at a temperature span of 10.2 K at a cooling load of 102.8 W... the COP of
3.1 is 11.3% of the Carnot Efficiency" and "an AMR operational frequency
of 0.75 Hz." Field: "current 1.13 T" (paper notes a future-work path to
1.4 T, confirming 1.13 T was the AS-TESTED field, not a design target).

**Result: f=0.75Hz, H=1.13T, Qc=102.8W, COP=3.1 all verified directly
against the primary source, word-for-word.** This is the strongest of the
three CORE points and needs no correction.

## 2. Astronautics_rotary_2014 — CONFIRMED, exact match (with one existing caveat re-confirmed)

Source in code: `("Astronautics_rotary_2014", 4.0, 1.44, 0.252999, 2502.0,
1133.70)`.

Checked against: `Papers/AMR systems and prototypes/The performance of a
large-scale rotary magnetic refrigerator.pdf` (Jacobs, Auringer, Boeder,
Chell, Komorowski, Leonard, Russek, Zimm — Astronautics/Xerox Fuel Cell
Systems; accepted manuscript, DOI 10.1016/j.ijrefrig.2013.09.025;
"Received Date: 30 April 2013" — the paper is dated 2013/2014, so the
project's "2014" label is reasonable for the published issue year, though
the accepted manuscript itself is dated 2013).

Direct text found: "2502 W over a span of 11.0 °C, corresponding to
[173 W/(T·L)]" and "electrical COP of 1.9. Moreover, this COP was obtained
using electrical components with mediocre efficiency" and "peak field of
1.44 tesla."

**Result: H=1.44T and Qc=2502W verified directly. COP=1.9 (electrical) at
that operating point is also directly stated in the source ("mediocre
efficiency" — the exact phrase this repo's own loss_model.py docstring
already paraphrases correctly).** One thing worth flagging that this
repo's own comments already partially note but is worth stating plainly:
the source device's beds are packed with **LaFeSiH** (six Curie-graded
layers), not Gd — the code's own comment ("Gd-only approximation...
intentionally different from validation_system.py's LAFESIH-based
system-validation row for the same device") already discloses this, so
no correction needed, but this audit independently confirms that
disclosure is accurate and necessary, not overcautious.

## 3. Tusek_singlebed_Gd_2010 — CONFIRMED DISCREPANCY (pre-existing, now independently verified)

Source in code: `("Tusek_singlebed_Gd_2010", 0.25, 1.69, 0.002422, 6.5,
0.76)`.

Checked against: `Papers/AMR Theory and Modeling/A comprehensive
experimental analysis of gadolinium active magnetic regenerators.pdf`
(Tušek, Kitanovski, Zupan, Prebil, Poredoš, Applied Thermal Engineering 53
(2013) 57-66 — this is the paper the code's own comments already point to
for "Figs. 10-11" digitization).

Direct text found for the actual reported operating point used in the
paper's own cooling-capacity/COP-vs-span figures (Figs. 10-11): **field
"1.15 T"** ("the magnet assembly provides a measured magnetic field of
about 1.15T"; "the finite magnetic flux density is 1.15T" in the air
gap), and **frequency "0.3 Hz"** ("obtained for an operating frequency of
0.3 Hz" — the caption text for Figs. 8-11), with total Gd mass **0.1763
kg** for the best-performing AMR (A.) configuration (parallel-plate,
0.1 mm spacing).

**This independently confirms — via direct primary-source text, not
re-derivation — the discrepancy this repo's own `loss_model.py` comment
already flags**: the CORE calibration point's field (1.69T) and frequency
(0.25Hz) do NOT match what this exact paper reports for the operating
point its own Figs. 10-11 curves are built from (1.15T, 0.3Hz). The
paper's own text does mention an "optimum frequency... between 0.22 and
0.3 Hz" range for the parallel-plate AMRs, so 0.25Hz is at least inside
the paper's discussed range (not fabricated), but 1.69T is not found
anywhere in this PDF's text — the paper's own reported field is
consistently 1.15T throughout. No occurrence of "1.69" was found in the
source text at all.

**Recommendation, not yet implemented in code (flagging honestly rather
than silently fixing)**: this repo's own prior comment already explains
why the corrected 1.15T value was tried and rejected ("does NOT calibrate
at all" under the current `cooling_capacity()` model) — this audit does
not change that finding, it only confirms the *source* side of the
discrepancy is real and not a typo introduced by this project. The
underlying gap (this device's real reported field doesn't calibrate
under the model) remains open, exactly as ROADMAP.md already documents,
and is a legitimate limitation to state plainly in the paper rather than
paper over.

## 4. Lozano POLO/UFSC calibration points — RESOLVED (earlier flag was a false alarm)

**Update (same Phase 30 pass, follow-up web search): this item is now
resolved.** The "2016" citation is correct. A dedicated web search located
the real paper: **Lozano, Capovilla, Trevizoli, Bahl, Engelbrecht,
Nielsen, Barbosa, "Development of a novel rotary magnetic refrigerator,"
International Journal of Refrigeration 68 (2016) 187-197** — confirmed
via both DTU's own institutional publication database (orbit.dtu.dk) and
ScienceDirect's listing, independently, as a genuinely 2016-dated paper.

The earlier version of this audit item was a false alarm caused by
checking the code's 2016-cited numbers against the WRONG Lozano paper —
this project's local `Papers/` corpus happens to contain a different,
earlier Lozano et al. paper (Int. J. Refrig. 37 (2014) 92-98, a distinct
device), and a direct text search of that PDF unsurprisingly found no
match for numbers that were never claimed to come from it.

The real 2016 paper's own reported headline figures (device: ~1.7 kg Gd
spheres, 8 regenerator-bed pairs, ~1 T rotor-stator field; results: max
zero-load span 12 K at 1 Hz/150 L·h⁻¹; max zero-span cooling power 150 W
at 0.8 Hz/200 L·h⁻¹; 80.4 W thermal load producing a 7.1 K span) are
consistent in device class, field strength, and order of magnitude with
the four `Lozano_POLO_UFSC_2016_r4/r6/r7/r8` (Qc, COP) pairs already in
`core/loss_model.py`. The full PDF was not directly fetchable in this
pass (paywalled/rate-limited), so the exact per-row digitized values
(62.5/81.2/80.8/120.4 W and their COP multipliers) were not individually
re-verified against the paper's own data tables — a much smaller,
lower-priority residual gap than "wrong citation year," and one that
doesn't block citing the paper itself. `core/loss_model.py`'s own comment
has been updated to reflect this resolution.

**Updated status: citation year CONFIRMED CORRECT (2016); exact per-point
digitized values still not individually re-verified (minor, downgraded
from the earlier "citation year mismatch" finding).**

## 5. MNFEPSI_FIRST_ORDER hysteresis (Hanggai et al. 2026) — CONFIRMED (already updated in prior pass)

Re-confirmed directly against `Papers/Magnetocaloric effect and materials
physics/Impact of F and S Doping on (Mn,Fe)2(P,Si) Giant Magnetocaloric
Materials.pdf` (Hanggai et al., Acta Materialia 302 (2026) 121677): the
exact calibrated composition's row (Mn=0.68) reports |ΔS_max,2T| = 13.3(3)
J/(kg·K) and ΔThys = 5.4(1) K, both used in the prior pass's update to
`core/first_order_mce.py`'s `hysteresis_loss_J_per_kg=71.8`. This audit
independently re-verified the table values directly from the PDF text
(not re-typed from the earlier pass's own notes) — no discrepancy found.

## 6. Bjørk field-vs-cost citation (magnet_geometry.py) — NOT independently re-checked this pass

`core/magnet_geometry.py`'s own comments already document and correct a
citation number (arXiv:1410.6248, not :1410.1987) from an earlier pass.
This audit did not re-derive that correction independently (the arXiv
paper itself is not a PDF present in this project's local `Papers/`
corpus to check against offline) — flagged as still resting on the
earlier pass's own correction, not re-verified against the primary source
in this session.

## 7. Web-search follow-up pass: Okamura & Hirano, Dan'kov et al., Giguere et al., Qian et al., Wu/Li et al.

### Okamura & Hirano (2013) -- secondary-source citation confirmed real, primary paper year uncertain

`core/loss_model.py`'s `EXTENDED` calibration set already documents that
this point is read from a REVIEW's table, not a primary paper. That
review is confirmed real and correctly cited: **Kamran, M.S., Ahmad,
H.O., Wang, H.S., "Review on the developments of active magnetic
regenerator refrigerators – Evaluated by performance," Renewable and
Sustainable Energy Reviews 133 (2020) 110247. DOI: 10.1016/j.rser.2020.110247.**
Confirmed via two independent listings (RePEc/IDEAS and the journal's own
metadata) with an exactly matching volume/page/DOI.

The underlying PRIMARY Okamura device paper could not be pinned to a
specific 2013 publication -- Okamura's own device work in the literature
trail appears as **Okamura, T., Yamada, K., Hirano, N., 2005, "Performance
of a room-temperature rotary magnetic refrigerator," Proc. 1st Int. Conf.
Magnetic Refrigeration at Room Temperature, Montreux, 319-324** and a
follow-up **Okamura, T., Rachi, R., Hirano, N., Nagaya, S., 2007,
"Improvement of 100W class room temperature magnetic refrigerator," Proc.
2nd Int. Conf. Magnetic Refrigeration at Room Temperature, Portoroz,
377-382** -- both non-DOI'd conference proceedings, neither dated 2013.
The "2013" in this repo's `Okamura_Hirano_2013` row label most likely
reflects how the SOURCE REVIEW (Kamran et al. 2020) itself dated or
indexed the entry in its Table 2, not a primary-paper publication year --
this is a device-labeling provenance note, not a numeric-value concern,
since the code already correctly treats this as a secondary-source point
and the review it's drawn from is confirmed real.

### Dan'kov, Tishin, Pecharsky & Gschneidner (1998) -- CONFIRMED real, DOI found

**S. Yu. Dan'kov, A. M. Tishin, V. K. Pecharsky, K. A. Gschneidner Jr.,
"Magnetic phase transitions and the magnetothermal properties of
gadolinium," Phys. Rev. B 57, 3478 (1998). DOI: 10.1103/PhysRevB.57.3478.**
Confirmed as a real, frequently-cited paper via multiple independent
reference-list appearances in unrelated later papers.

### Giguere, Foldeaki, Gopal, Chahine, Bose, Frydman & Barclay (1999) -- CONFIRMED real, DOI found

**A. Giguere, M. Foldeaki, B. R. Gopal, R. Chahine, T. K. Bose, A. Frydman,
J. A. Barclay, "Direct Measurement of the 'Giant' Adiabatic Temperature
Change in Gd5Si2Ge2," Phys. Rev. Lett. 83, 2262 (1999). DOI:
10.1103/PhysRevLett.83.2262.** Confirmed as a real paper, independently
cited in multiple unrelated later Gd5Si2Ge2 papers with matching
volume/page.

### Qian et al. (2023) and "Wu et al." (2023) -- Qian confirmed, "Wu" was a MISATTRIBUTION (now fixed)

**Qian, Catalini, Muehlbauer, Liu, Mevada, Hou, Hwang, Radermacher &
Takeuchi, "High-performance multimode elastocaloric cooling system,"
Science 380 (6646), 722-727 (2023). DOI: 10.1126/science.adg7043.**
Confirmed real via an independent Wikipedia citation cross-reference with
exactly matching volume/issue/pages/DOI/PMID.

**Genuine finding, now fixed in `core/baseline_cooling.py`:** the second
elastocaloric anchor was attributed in this repo's code as "Wu et al.
(2023)," but the actual paper -- confirmed via PubMed Central, matching
exactly on journal/volume/page/year/DOI/the specific COP=3.7 and ~1K-span
figures the code already uses -- is **Li, Hua & Sun, "Continuous and
efficient elastocaloric air cooling by coil-bending," Nature
Communications 14, 7982 (2023). DOI: 10.1038/s41467-023-43611-6.**
Authors: Xueshi Li, Peng Hua, Qingping Sun (Hong Kong University of
Science and Technology) -- no author named Wu appears on this paper. The
journal, volume, page, year, and the specific COP=3.7 (device-level
system COP) and ~1K span figures were all independently re-verified this
pass and are correct; only the author surname was wrong, likely
confused with a different, similarly-themed elastocaloric paper by a
Wu-led group. Fixed in `core/baseline_cooling.py`'s
`ELASTOCALORIC_COP_SOURCE_NOTE` and surrounding comments.

## Summary

| Calibration point | Status |
|---|---|
| DTU_Eriksen_rotary_Gd_2015 | **Confirmed, exact match** |
| Astronautics_rotary_2014 | **Confirmed, exact match** |
| Tusek_singlebed_Gd_2010 | Field/frequency mismatch independently reconfirmed (pre-existing, documented, unresolved) |
| Lozano_POLO_UFSC_2016 (×4, diagnostic-only) | **Resolved: citation year confirmed correct (Int. J. Refrig. 68 (2016) 187-197); exact per-row digits not yet independently re-verified** |
| MNFEPSI_FIRST_ORDER hysteresis (Hanggai 2026) | **Confirmed** (re-verified from prior pass) |
| Bjørk field-vs-cost arXiv correction | Not re-checked (paper not available locally) |
| Okamura_Hirano_2013 (secondary source) | Source review (Kamran et al. 2020, DOI 10.1016/j.rser.2020.110247) confirmed real; primary paper year still unpinned |
| Dan'kov et al. (1998), Phys. Rev. B 57, 3478 | **Confirmed real, DOI 10.1103/PhysRevB.57.3478** |
| Giguere et al. (1999), Phys. Rev. Lett. 83, 2262 | **Confirmed real, DOI 10.1103/PhysRevLett.83.2262** |
| Qian et al. (2023), Science 380, 722-727 | **Confirmed real, DOI 10.1126/science.adg7043** |
| "Wu et al. (2023)" elastocaloric anchor | **Misattribution found and fixed: actual paper is Li, Hua & Sun (2023) Nat. Commun. 14, 7982, DOI 10.1038/s41467-023-43611-6** |

Bottom line for the paper: the production CORE calibration's two cleanest
points (DTU, Astronautics) are independently, directly verified against
primary-source PDF text, word-for-word. The one already-flagged weak
point (Tušek) is confirmed weak by the same method, not newly discovered.
The diagnostic-only FURTHER_EXTENDED set's "Lozano 2016" citation, initially
flagged as a possible year mismatch, was traced to the correct real paper
(Lozano, Capovilla, Trevizoli et al., Int. J. Refrig. 68 (2016) 187-197)
via a follow-up web search and confirmed genuinely 2016-dated with
consistent device-class/order-of-magnitude figures -- the citation itself
is sound; only the exact per-row digitized values remain formally
unverified, a minor residual item, not a provenance failure.
