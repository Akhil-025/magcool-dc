# Tušek et al. (2013) Figs. 10-11 — digitization attempt notes

Source: Tušek, Kitanovski, Zupan, Prebil, Poredoš (2013), "A comprehensive
experimental analysis of gadolinium active magnetic regenerators," Appl.
Therm. Eng. 53, 57-66. `Papers/AMR Theory and Modeling/A comprehensive
experimental analysis of gadolinium active magnetic regenerators.pdf`.

## What was actually done

Extracted the two embedded raster images at native resolution (476x1093 and
476x1095 px — pdfplumber confirms 0 vector paths on this page, only 2
embedded bitmaps, so there is no vector-graphics shortcut; this has to be
genuine pixel work).

OCR (tesseract, 3x upscale) confirmed real structure with high confidence
(90-96%): each figure is **three stacked subplots**, each labeled "cooling
capacity [W]" (Fig. 10) — one subplot per AMR geometry (A, B, F per the
paper's text) — with a shared x-axis labeled "temperature span [K]",
ticked at 5/10/15 K (OCR confidence 95-96% on two of the three subplots'
tick labels; the third subplot's "15" tick did not OCR cleanly).

Row-wise dark-pixel-density analysis locates candidate horizontal
gridlines at ~70px vertical spacing within each subplot band, consistent
with evenly spaced y-axis gridlines.

## Where this stopped

- **Y-axis numeric tick labels did not OCR reliably** at this resolution/
  font size (tried whitelisted-digit OCR on the tick-label column; no
  confident hits). Without OCR'd tick values, the y-axis calibration
  (W per pixel) rests on the gridline *spacing* alone, not a confirmed
  absolute scale — usable as a rough guide, not as a citable calibration.
- **Series separation was not attempted.** The image is grayscale with no
  color channel to distinguish the 3 geometries x 3 flow-ratio V*
  combinations (9 lines per subplot); distinguishing them requires marker
  *shape* recognition (circle/square/triangle etc.) at native resolution,
  which needs sub-pixel template matching or a human doing point-and-click
  digitization (e.g. WebPlotDigitizer), not blind pixel statistics. This
  is exactly the false-precision trap ROADMAP.md already flagged: a
  confident-looking auto-extracted table here would be worse than no
  table, because 9 overlapping series without color cues cannot be
  reliably auto-separated.

## Bottom line

Real, checkable progress was made (confirmed image resolution, confirmed
3-subplot structure and x-axis calibration via OCR, located candidate
y-gridline spacing) but the item is still not digitized to the standard
`validation_system.py` would need for a system-level curve check. The
concrete next step is unchanged from the prior note: a human walking
through Figs. 10-11 series-by-series with a point-and-click digitizer,
using the x-axis calibration above as a starting point.