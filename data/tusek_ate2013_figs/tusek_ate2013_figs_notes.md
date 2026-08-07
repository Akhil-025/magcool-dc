# Tusek et al. (2013) Figs. 10-11 digitization notes

Source: Tusek, Kitanovski, Zupan, Prebil, Poredos, "A comprehensive
experimental analysis of gadolinium active magnetic regenerators,"
Appl. Therm. Eng. 53 (2013) 57-66. PDF supplied directly (page 8 of the
10-page article contains both figures).

This supersedes the earlier `results/tusek_ate2013_figs_notes.md`
"non-authoritative, by-eye" placeholder referenced in
`core/validation_system.py`'s module docstring and in the
`Tusek_singlebed_Gd_2010` CSV row. This is a real pixel-calibrated
digitization, not a by-eye guess.

## Method

1. Both figures are embedded raster images, confirmed via `pdfplumber`
   (`page.images`) to be the only image content on the page -- 0 vector
   paths, 0 rects on that page's drawing layer -- so pixel-level digitization
   was required exactly as ROADMAP.md anticipated. Extracted at native
   resolution with PyMuPDF: Fig. 10 is 476x1093 px, Fig. 11 is 476x1095 px.

2. Each figure is 3 stacked subplots (V*=0.16, 0.42, 0.95), each with 3
   series (AMR A, B, F -- circle/square/star-or-X markers respectively).
   Panel boundaries were found from rows of near-pure-white pixels.

3. **Axis calibration** (per panel, done independently for all 6 panels
   across both figures): the plot's bounding box and horizontal gridlines
   are solid dark rows/columns spanning (almost) the full plot width/height,
   detected by thresholding and counting dark pixels per row/column.
   - Fig. 10: gridlines at y = 8, 6, 4, 2, 0 W, evenly spaced 70 px apart in
     every panel (confirmed independently per panel).
   - Fig. 11: gridlines at y = 40, 30, 20, 10, 0 (COP), same 70 px spacing
     pattern.
   - x-axis: left axis (x=0) and right border (x=20 K) pixel columns were
     used directly (no intermediate tick detection needed -- the plot's
     xlim exactly matches the drawn box, verified against the y=0 axis row
     spanning exactly between those two columns).
   - Resulting scale: ~19.8 px/K (x), 35 px per 2 W or 35 px per 5 COP units
     (y), consistent across all panels to within 1 px.

4. **Marker detection**: binary threshold (pixel < 130) + one iteration of
   3x3 binary erosion to strip the ~1-2 px connecting line strokes while
   preserving the larger (~7x7 px) marker blobs, then connected-component
   labeling (`scipy.ndimage.label`) and centroid extraction. Legend-box
   glyphs were excluded by position (top-right corner of each panel).

5. **Series assignment and disambiguation**: automated detection gives
   precise pixel->data coordinates but does not reliably distinguish marker
   *shape* (open circle = A, filled square = B, star/X = F) where two
   markers sit close together (e.g. near curve crossings) or where a thin
   X-marker eroded to a smaller blob than expected. Every panel's candidate
   point list was cross-checked against a 3-6x zoomed crop of that region
   of the source image and assigned/corrected by eye against the visible
   marker glyph -- this is the "one series at a time, by eye" step
   ROADMAP.md flagged as unavoidable; the pixel calibration above just makes
   each by-eye read quantitatively precise rather than a ruler-on-screen
   guess. A handful of spurious sub-threshold blobs (line self-crossings,
   not real markers) were identified this way and discarded.

## Known uncertainty / open flag

- **Fig. 11, V*=0.95, AMR F**: one candidate point at (span=4.78 K,
  COP=4.39) has no counterpart in Fig. 10's V*=0.95/AMR F curve (which
  only has points at span=7.07 and 12.85 K in that region -- verified by
  re-zooming the corresponding Fig. 10 region, which shows no marker near
  span=4.78 K). It is visually a real X-shaped marker in Fig. 11, not an
  artifact, so it has been kept in `fig11_data.csv` but flagged
  `uncertain_no_fig10_counterpart`. Possible explanations: COP was
  reported at an additional operating point not shown in the Qc plot, or
  a misread. Treat this single point with extra caution; it is not used
  in any validation calibration below.
- The two figures' "shared endpoint" markers (where AMR A and AMR B, or
  A and F, curves converge near their common zero-capacity companion span)
  render as a single overlapping glyph in the source image. Where the two
  series' true endpoints are close but not pixel-identical, a single
  representative (span, value) pair is used for both and flagged
  `shared_endpoint_with_*`.
- **Estimated pixel-to-value uncertainty**: +/-0.15 K on span, +/-0.15 W on
  Qc, +/-0.3 on COP for well-isolated markers (based on the ~1 px
  calibration-line uncertainty and typical ~3 px marker-centroid
  uncertainty). Points near curve crossings or shared endpoints (flagged
  above) should be treated as +/-0.5 K / +/-0.3 W / +/-1 COP.

## Files

- `fig10_data.csv`: span_K, Qc_W for all 9 series (3 AMR geometries x 3
  V* flow ratios) digitized from Fig. 10.
- `fig11_data.csv`: span_K, COP for the same 9 series, digitized from
  Fig. 11.

## Not digitized

Fig. 6's single reported point (19.8 K span at 0 W cooling capacity,
V*=0.365, freq=0.3 Hz, AMR A) was **not** digitized -- the paper states
this number directly in text/caption ("19.8 K", Section 3.2, Fig. 6
caption), so it is used as-is in the benchmark CSV without pixel reading.