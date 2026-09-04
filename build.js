const pptxgen = require("pptxgenjs");
const path = require("path");

const BG = "E3E8F0";
const DARK = "1C4A48";
const PINK = "F3B7C5";
const PINKDEEP = "EDA0B3";
const TEXT = "27504D";
const WHITE = "FFFFFF";
const MUTED = "5C7B79";
const AMBER = "C97A2B";
const AMBERBG = "FBEBD9";

const FIGDIR = "slide_figs";

let pres = new pptxgen();
pres.layout = "LAYOUT_WIDE"; // 13.33 x 7.5
const PW = 13.33, PH = 7.5;

let pageNum = 1;

function corner(slide, dark) {
  const c1 = dark ? PINK : DARK;
  const c2 = dark ? WHITE : PINK;
  slide.addShape("line", { x: 0.35, y: PH - 0.35, w: 0.9, h: 0, line: { color: c1, width: 1.5 } });
  slide.addShape("line", { x: 0.35, y: PH - 1.25, w: 0, h: 0.9, line: { color: c1, width: 1.5 } });
  slide.addShape("line", { x: 0.55, y: PH - 0.55, w: 0.65, h: 0, line: { color: c2, width: 1.5 } });
  slide.addShape("line", { x: 0.55, y: PH - 1.2, w: 0, h: 0.65, line: { color: c2, width: 1.5 } });
  slide.addShape("line", { x: PW - 1.25, y: 0.35, w: 0.9, h: 0, line: { color: c1, width: 1.5 } });
  slide.addShape("line", { x: PW - 0.35, y: 0.35, w: 0, h: 0.9, line: { color: c1, width: 1.5 } });
  slide.addShape("line", { x: PW - 1.2, y: 0.55, w: 0.65, h: 0, line: { color: c2, width: 1.5 } });
  slide.addShape("line", { x: PW - 0.55, y: 0.55, w: 0, h: 0.65, line: { color: c2, width: 1.5 } });
}

function bgSlide(dark) {
  let s = pres.addSlide();
  s.background = { color: dark ? DARK : BG };
  return s;
}

function divider(slide) {
  slide.addShape("line", { x: 1.0, y: 1.55, w: PW - 2.0, h: 0, line: { color: DARK, width: 0.75, dashType: "solid" } });
  slide.addShape("ellipse", { x: PW / 2 - 0.06, y: 1.52, w: 0.12, h: 0.06, fill: { color: PINK }, line: { type: "none" } });
}

function title(slide, kicker, text, opts) {
  opts = opts || {};
  if (kicker) {
    slide.addText(kicker.toUpperCase(), {
      x: 0.7, y: 0.5, w: PW - 1.4, h: 0.35, fontFace: "Calibri", fontSize: 13,
      color: PINKDEEP, bold: true, charSpacing: 2, isTextBox: true, margin: 0
    });
  }
  slide.addText(text, {
    x: 0.7, y: kicker ? 0.82 : 0.55, w: PW - 1.4, h: 0.7, fontFace: "Cambria", fontSize: opts.size || 30,
    color: DARK, bold: true, isTextBox: true, margin: 0
  });
  if (opts.divider !== false) divider(slide);
}

function footer(slide) {
  pageNum++;
  slide.addText(String(pageNum), {
    x: PW - 0.9, y: PH - 0.55, w: 0.5, h: 0.3, fontFace: "Calibri", fontSize: 10,
    color: MUTED, align: "right", isTextBox: true, margin: 0
  });
  slide.addText("MAGCOOL-DC \u00B7 ASHRAE REGION XV \u00B7 SAMUDRA 2026", {
    x: 1.55, y: PH - 0.55, w: 6, h: 0.3, fontFace: "Calibri", fontSize: 9,
    color: MUTED, isTextBox: true, margin: 0
  });
}

function statBadge(slide, x, y, w, h, value, label, color, opts) {
  opts = opts || {};
  slide.addShape("roundRect", { x, y, w, h, rectRadius: 0.08, fill: { color: WHITE }, line: { color: color, width: 1 }, shadow: { type: "outer", color: "3A5A58", opacity: 0.15, blur: 4, offset: 2, angle: 90 } });
  const autoSize = value.length > 14 ? 15 : value.length > 10 ? 18 : 22;
  slide.addText(value, { x: x + 0.06, y: y + 0.12, w: w - 0.12, h: h * 0.55, align: "center", fontFace: "Cambria", fontSize: opts.valueSize || autoSize, bold: true, color: color, isTextBox: true, margin: 0, fit: "shrink" });
  slide.addText(label, { x: x + 0.1, y: y + h * 0.6, w: w - 0.2, h: h * 0.35, align: "center", fontFace: "Calibri", fontSize: 10, color: TEXT, isTextBox: true, margin: 0, fit: "shrink" });
}

function card(slide, x, y, w, h, heading, body, opts) {
  opts = opts || {};
  const headBg = opts.dark ? DARK : PINK;
  const headColor = opts.dark ? WHITE : DARK;
  slide.addShape("roundRect", { x, y, w, h, rectRadius: 0.07, fill: { color: WHITE }, line: { color: "C7D1DE", width: 0.75 } });
  slide.addShape("roundRect", { x, y, w, h: 0.5, rectRadius: 0.07, fill: { color: headBg }, line: { type: "none" } });
  slide.addShape("rect", { x, y: y + 0.25, w, h: 0.25, fill: { color: headBg }, line: { type: "none" } });
  slide.addText(heading, { x: x + 0.18, y, w: w - 0.36, h: 0.5, valign: "middle", fontFace: "Calibri", fontSize: 13, bold: true, color: headColor, isTextBox: true, margin: 0 });
  slide.addText(body, { x: x + 0.18, y: y + 0.62, w: w - 0.36, h: h - 0.8, fontFace: "Calibri", fontSize: 11, color: TEXT, isTextBox: true, margin: 0, valign: "top", lineSpacingMultiple: 1.15 });
}

function bulletList(slide, x, y, w, h, items, opts) {
  opts = opts || {};
  const paras = items.map((t) => ({
    text: t, options: { bullet: { code: "25AA", color: PINKDEEP }, color: TEXT, fontSize: opts.size || 14,
      fontFace: "Calibri", breakLine: true, paraSpaceAfter: 10 }
  }));
  slide.addText(paras, { x, y, w, h, isTextBox: true, margin: 0, valign: "top" });
}

function img(slide, file, x, y, w, h) {
  slide.addShape("roundRect", { x: x - 0.06, y: y - 0.06, w: w + 0.12, h: h + 0.12, rectRadius: 0.06, fill: { color: WHITE }, line: { color: "C7D1DE", width: 0.75 }, shadow: { type: "outer", color: "3A5A58", opacity: 0.18, blur: 6, offset: 2, angle: 90 } });
  slide.addImage({ path: path.join(FIGDIR, file), x, y, w, h, sizing: { type: "contain", w, h } });
}

function chipRow(slide, chips, y) {
  const w = 3.35, h = 0.72, gap = 0.28, perRow = 3;
  chips.forEach((label, i) => {
    const colIdx = i % perRow;
    const rowIdx = Math.floor(i / perRow);
    const rowItems = Math.min(perRow, chips.length - rowIdx * perRow);
    const rowW = rowItems * w + (rowItems - 1) * gap;
    const startX = (PW - rowW) / 2;
    const x = startX + colIdx * (w + gap);
    const cy = y + rowIdx * (h + 0.25);
    const dark = i % 2 === 0;
    slide.addShape("rect", { x, y: cy, w, h, fill: { color: dark ? DARK : PINK }, line: { type: "none" } });
    slide.addText(label, { x, y: cy, w, h, align: "center", valign: "middle", fontFace: "Calibri", fontSize: 13, color: dark ? WHITE : DARK, isTextBox: true, margin: 0 });
  });
}

function sectionDivider(partLabel, titleText, chips) {
  let s = bgSlide(false);
  corner(s, false);
  s.addText(partLabel.toUpperCase(), {
    x: 0.7, y: 1.75, w: PW - 1.4, h: 0.4, align: "center", fontFace: "Calibri", fontSize: 14,
    color: PINKDEEP, bold: true, charSpacing: 2, isTextBox: true, margin: 0
  });
  s.addText(titleText, {
    x: 0.9, y: 2.15, w: PW - 1.8, h: 1.15, align: "center", fontFace: "Cambria", fontSize: 32,
    color: DARK, bold: true, isTextBox: true, margin: 0
  });
  s.addShape("line", { x: PW / 2 - 2.4, y: 3.5, w: 4.8, h: 0, line: { color: DARK, width: 0.75 } });
  s.addShape("ellipse", { x: PW / 2 - 0.06, y: 3.47, w: 0.12, h: 0.06, fill: { color: PINK }, line: { type: "none" } });
  chipRow(s, chips, 4.15);
  footer(s);
  return s;
}

// ---------------- SLIDE: TITLE ----------------
{
  let s = bgSlide(true);
  corner(s, true);
  s.addText("PHYSICS-BASED SIMULATION & MULTI-OBJECTIVE OPTIMISATION", {
    x: 1.0, y: 1.85, w: 11.3, h: 0.4, fontFace: "Calibri", fontSize: 14, color: PINK, bold: true, charSpacing: 2, isTextBox: true, margin: 0, align: "center"
  });
  s.addText("MAGNETOCALORIC COOLING\nFOR DATA CENTERS", {
    x: 0.8, y: 2.2, w: 11.73, h: 1.9, fontFace: "Cambria", fontSize: 44, color: WHITE, bold: true, align: "center", isTextBox: true, margin: 0
  });
  s.addText("Toward Refrigerant-Free, Low-Carbon HVAC&R", {
    x: 1.0, y: 4.05, w: 11.3, h: 0.4, fontFace: "Calibri", fontSize: 15, italic: true, color: "BFE0DD", align: "center", isTextBox: true, margin: 0
  });
  s.addShape("line", { x: PW / 2 - 1.3, y: 4.65, w: 2.6, h: 0, line: { color: PINK, width: 1 } });
  
  // AUTHORS ADDED HERE
  s.addText("Presented by: AKHIL PILLAI   •   AMRUTA PATIL", {
    x: 1.0, y: 4.95, w: 11.3, h: 0.35, fontFace: "Calibri", fontSize: 14, color: WHITE, bold: true, align: "center", isTextBox: true, margin: 0
  });
  s.addText("ASHRAE Region XV — 3rd CRC 2026  |  Student Competition  |  Theme: SAMUDRA", {
    x: 1.0, y: 5.45, w: 11.3, h: 0.35, fontFace: "Calibri", fontSize: 12, color: "BFE0DD", align: "center", isTextBox: true, margin: 0
  });
  s.addText("Department of Mechanical Engineering, SPCE, Mumbai  |  Guide: Kunal Bhavsar", {
    x: 1.0, y: 5.75, w: 11.3, h: 0.35, fontFace: "Calibri", fontSize: 12, color: "BFE0DD", align: "center", isTextBox: true, margin: 0
  });
}

// ================= PART I: MOTIVATION & METHODOLOGY =================
sectionDivider("Part I", "Motivation & Methodology", [
  "Timeline", "Why This Matters", "The AMR Cycle", "Research Gap", "Key Literature", "Objectives", "Methodology Pipeline"
]);

// ---------------- SLIDE: TIMELINE ----------------
{
  let s = bgSlide(false);
  corner(s, false);
  title(s, "Project Journey", "Timeline");
  const timeline = [
    ["1", "Feb", "Literature review; mean-field & Landau material models"],
    ["2", "Mar", "AMR cycle model; Gd / benchmark-device validation"],
    ["3", "Apr", "Loss calibration; Sobol sensitivity; RSM surrogate"],
    ["4", "May", "NSGA-III multi-objective optimisation; cascade design"],
    ["5", "Jun–Jul", "Economics, emissions, extension studies"],
    ["6", "Aug–Sep", "Design recommendations; report & CRC presentation"]
  ];
  const w = 3.6, gap = 0.35, startX = 1.0;
  timeline.forEach((item, i) => {
    const colIdx = i % 3;
    const rowIdx = Math.floor(i / 3);
    const x = startX + colIdx * (w + gap);
    const y = 2.15 + rowIdx * 1.8;
    statBadge(s, x, y, w, 1.4, item[1], item[2], DARK, { valueSize: 18 });
  });
  footer(s);
}

// ---------------- SLIDE: WHY THIS MATTERS ----------------
{
  let s = bgSlide(false);
  corner(s, false);
  title(s, "The Problem", "Why This Matters");
  const items = [
    ["~40%", "of a data center's total energy goes to cooling — the single largest non-IT load"],
    ["High GWP", "HFC/HFO refrigerants in vapor-compression and liquid cooling carry high Global Warming Potential"],
    ["Solid-State", "Magnetocaloric cooling promises a refrigerant-free alternative — no compressor, no leak pathway"]
  ];
  const w = 3.6, gap = 0.35, startX = (PW - (w * 3 + gap * 2)) / 2;
  items.forEach((it, i) => {
    const x = startX + i * (w + gap);
    statBadge(s, x, 2.1, w, 1.5, it[0], "", i === 0 ? DARK : (i === 1 ? "B23A56" : DARK));
    s.addText(it[1], { x, y: 3.75, w, h: 2.2, fontFace: "Calibri", fontSize: 13, color: TEXT, align: "center", isTextBox: true, margin: 0, valign: "top" });
  });
  footer(s);
}

// ---------------- SLIDE: AMR CYCLE ----------------
{
  let s = bgSlide(false);
  corner(s, false);
  title(s, "The Core Physics", "The Active Magnetic Regenerator (AMR) Cycle");
  const stages = [
    ["1", "MAGNETIZE", "Field applied — material heats (adiabatic)"],
    ["2", "HOT BLOW", "Heat rejected to the sink"],
    ["3", "DEMAGNETIZE", "Field removed — material cools (adiabatic)"],
    ["4", "COLD BLOW", "Heat absorbed from the load"]
  ];
  const w = 2.75, gap = 0.28, startX = (PW - (w * 4 + gap * 3)) / 2, y = 2.3;
  stages.forEach((st, i) => {
    const x = startX + i * (w + gap);
    s.addShape("ellipse", { x: x + w / 2 - 0.32, y, w: 0.64, h: 0.64, fill: { color: i % 2 === 0 ? DARK : PINK }, line: { type: "none" } });
    s.addText(st[0], { x: x + w / 2 - 0.32, y, w: 0.64, h: 0.64, align: "center", valign: "middle", fontFace: "Cambria", bold: true, fontSize: 20, color: i % 2 === 0 ? WHITE : DARK, isTextBox: true, margin: 0 });
    s.addText(st[1], { x, y: y + 0.85, w, h: 0.35, align: "center", fontFace: "Calibri", bold: true, fontSize: 13, color: DARK, isTextBox: true, margin: 0 });
    s.addText(st[2], { x, y: y + 1.25, w, h: 1.1, align: "center", fontFace: "Calibri", fontSize: 11, color: TEXT, isTextBox: true, margin: 0, valign: "top" });
    if (i < 3) s.addText("→", { x: x + w, y: y + 0.05, w: gap, h: 0.6, align: "center", fontFace: "Calibri", fontSize: 22, color: MUTED, isTextBox: true, margin: 0 });
  });
  s.addText("Regeneration amplifies one stage's small ΔT into a usable temperature span across the bed.", {
    x: 1.5, y: 5.1, w: 10.3, h: 0.6, align: "center", italic: true, fontFace: "Calibri", fontSize: 13, color: DARK, isTextBox: true, margin: 0
  });
  footer(s);
}

// ---------------- SLIDE: AMR UNIT SCHEMATIC ----------------
{
  let s = bgSlide(false);
  corner(s, false);
  title(s, "Hardware Layout", "Anatomy of the AMR Unit");
  img(s, "amr_unit_schematic.png", 3.9, 2.1, 5.6, 4.55);
  bulletList(s, 0.9, 2.1, 2.85, 4.4, [
    "NdFeB permanent-magnet assembly cycles field across the packed bed",
    "Regenerator bed: packed spheres or parallel plates of the MCM",
    "Hot/cold-side heat exchangers reject/absorb heat each half-cycle",
    "Pump drives working fluid through the bed in sync with the field"
  ], { size: 12.5 });
  footer(s);
}

// ---------------- SLIDE: RESEARCH GAP ----------------
{
  let s = bgSlide(false);
  corner(s, false);
  title(s, "Positioning", "Research Gap");
  const rows = [
    [{ text: "Gap", options: { bold: true, color: WHITE, fill: { color: DARK } } },
     { text: "This Work", options: { bold: true, color: WHITE, fill: { color: DARK } } },
     { text: "Typical Literature", options: { bold: true, color: WHITE, fill: { color: DARK } } }],
    ["Loss modeling", "State-dependent (f, ṁ, Qc-based)", "Constant parasitic-fraction assumption"],
    ["Design scope", "Material + geometry co-optimized", "Fixed material, manual geometry choice"],
    ["Honesty on gap vs. VCC", "Explicit crossover search, reported", "Rarely benchmarked against real HVAC"],
    ["Scale", "Data-center 5–20K span", "Mostly domestic-refrigeration scale"]
  ].map((r, ri) => ri === 0 ? r : r.map((c, ci) => ({ text: c, options: { color: TEXT, fill: { color: ri % 2 === 0 ? "FFFFFF" : "F1F4F9" }, bold: ci === 0 } })));
  s.addTable(rows, {
    x: 0.9, y: 2.15, w: 11.5, h: 3.6, fontFace: "Calibri", fontSize: 13, border: { type: "solid", color: "C7D1DE", pt: 0.75 },
    autoPage: false, colW: [3.0, 4.5, 4.0], valign: "middle", margin: [6, 8, 6, 8]
  });
  footer(s);
}

// ---------------- SLIDE: KEY LITERATURE ----------------
{
  let s = bgSlide(false);
  corner(s, false);
  title(s, "Standing on Six Sources", "Key Literature This Work Builds On");
  const rows = [
    [{ text: "Foundation", options: { bold: true, color: WHITE, fill: { color: DARK } } },
     { text: "Source", options: { bold: true, color: WHITE, fill: { color: DARK } } },
     { text: "What It Anchors in This Work", options: { bold: true, color: WHITE, fill: { color: DARK } } }],
    ["Giant MCE discovery", "Pecharsky & Gschneidner, PRL 78, 4494 (1997)", "First-order Gd\u2085Si\u2082Ge\u2082 transition — motivates the Landau extension"],
    ["Gd benchmark data", "Dan\u2019kov et al., PRB 57, 3478 (1998)", "\u0394T_ad at 1/2/5 T — primary material-model validation target"],
    ["AMR concept", "Barclay, US Patent 4332135 (1982)", "Defines the four-stage cycle this model simulates"],
    ["Prototype benchmarks", "Kamran, Ahmad & Wang (2020); Greco et al. (2019)", "Comparative device tables behind this repo's 16-device benchmark set"],
    ["Loss & geometry physics", "Tu\u0161ek et al. (2013); Klinar et al. (2024)", "Geometry-dependent pumping power and the Hypereg parallel-hydraulic model"],
    ["Mean-field limits", "de Oliveira & von Ranke, Phys. Rep. 489 (2010)", "Documents the near-T\u1d04 bias this model's own validation reproduces"]
  ].map((r, ri) => ri === 0 ? r : r.map((c, ci) => ({ text: c, options: { color: TEXT, fill: { color: ri % 2 === 0 ? "FFFFFF" : "F1F4F9" }, bold: ci === 0, fontSize: 10.5 } })));
  s.addTable(rows, {
    x: 0.9, y: 2.05, w: 11.5, h: 4.35, fontFace: "Calibri", fontSize: 11, border: { type: "solid", color: "C7D1DE", pt: 0.75 },
    autoPage: false, colW: [2.5, 3.6, 5.4], valign: "middle", margin: [6, 8, 6, 8]
  });
  footer(s);
}

// ---------------- SLIDE: OBJECTIVES ----------------
{
  let s = bgSlide(false);
  corner(s, false);
  title(s, "Scope", "Objectives");
  const cards = [
    ["Material & Cycle Modeling", "• Mean-field (Gd) + Landau (giant-MCE)\n• 0-D AMR cooling/COP model\n• Validate vs. literature & real devices"],
    ["Loss Calibration & Sensitivity", "• State-dependent eddy/pumping/base loss model\n• Calibrate to real AMR prototypes\n• Sobol sensitivity analysis"],
    ["Multi-Objective Optimization", "• RSM surrogate + NSGA-III\n• Techno-economic + emissions comparison\n• Assess feasibility vs. conventional HVAC"]
  ];
  const w = 3.75, gap = 0.3, startX = (PW - (w * 3 + gap * 2)) / 2;
  cards.forEach((c, i) => card(s, startX + i * (w + gap), 2.15, w, 3.9, c[0], c[1], { dark: i === 1 }));
  footer(s);
}

// ---------------- SLIDE: METHODOLOGY PIPELINE ----------------
{
  let s = bgSlide(false);
  corner(s, false);
  title(s, "How It All Fits Together", "Methodology Pipeline");
  const stages = ["Material\nPhysics", "AMR Cycle\nModel", "Benchmark\nValidation", "Loss\nCalibration", "Sobol\nSensitivity", "RSM\nSurrogate", "NSGA-III\nOptimization", "Design\nRecs", "Economics &\nEmissions"];
  const w = 1.28, gap = 0.06, startX = (PW - (w * 9 + gap * 8)) / 2, y = 2.5;
  stages.forEach((st, i) => {
    const x = startX + i * (w + gap);
    s.addShape("roundRect", { x, y, w, h: 1.0, rectRadius: 0.06, fill: { color: i % 2 === 0 ? DARK : PINK }, line: { type: "none" } });
    s.addText(st, { x: x + 0.05, y, w: w - 0.1, h: 1.0, align: "center", valign: "middle", fontFace: "Calibri", fontSize: 9.5, bold: true, color: i % 2 === 0 ? WHITE : DARK, isTextBox: true, margin: 0 });
    if (i < 8) s.addText("→", { x: x + w - 0.02, y: y + 0.32, w: gap + 0.04, h: 0.4, align: "center", fontFace: "Calibri", fontSize: 12, color: MUTED, isTextBox: true, margin: 0 });
  });
  s.addText("9 stages  ·  ~30 core modules  ·  452 automated tests  ·  35 auto-generated figures", {
    x: 1.5, y: 3.9, w: 10.3, h: 0.5, align: "center", fontFace: "Cambria", fontSize: 16, bold: true, color: DARK, isTextBox: true, margin: 0
  });
  s.addText("Each stage feeds validated outputs into the next, ensuring traceable, reproducible results.", {
    x: 1.5, y: 4.45, w: 10.3, h: 0.5, align: "center", italic: true, fontFace: "Calibri", fontSize: 12, color: MUTED, isTextBox: true, margin: 0
  });
  footer(s);
}

// ================= PART II: MATERIAL & SYSTEM MODELING =================
sectionDivider("Part II", "Material & System Modeling", [
  "Material Physics Primer", "Material Validation", "Giant-MCE Honesty Flag", "AMR Device Physics", "Loss Modeling", "System Validation"
]);

// ---------------- SLIDE: MATERIAL PHYSICS PRIMER ----------------
{
  let s = bgSlide(false);
  corner(s, false);
  title(s, "Foundations", "How the Model \u201cSees\u201d a Magnetocaloric Material");
  img(s, "slide8.png", 6.55, 2.2, 6.1, 2.7);
  bulletList(s, 0.9, 2.2, 5.3, 4.0, [
    "Gd: mean-field (Brillouin/Weiss) theory, H_eff = H + λM — λ fixed by the known Curie temperature, not fit",
    "Giant-MCE families (Gd₅Si₂Ge₂, La(Fe,Si)₁₃Hy, (Mn,Fe)₂(P,Si)): 6th-order Landau free energy",
    "Needed because these are first-order transitions mean-field theory cannot capture"
  ]);
  s.addText("Entropy change & ΔT_ad vs. temperature (1/2/5 T), peak at Tc = 294 K", { x: 6.55, y: 5.0, w: 6.1, h: 0.4, italic: true, fontFace: "Calibri", fontSize: 10.5, color: MUTED, isTextBox: true, margin: 0 });
  footer(s);
}

// ---------------- SLIDE: MATERIAL VALIDATION ----------------
{
  let s = bgSlide(false);
  corner(s, false);
  title(s, "Validation", "Gd Model Validation vs. Dan\u2019kov et al. (1998)");
  const w = 1.9, gap = 0.25, startX = 0.9;
  const vals = [["+4.2%", "at 1 T", DARK], ["\u22121.1%", "at 2 T", "B08A2B"], ["\u22128.1%", "at 5 T", "B23A56"]];
  vals.forEach((v, i) => statBadge(s, startX + i * (w + gap), 2.2, w, 1.3, v[0], v[1], v[2]));
  bulletList(s, 0.9, 3.85, 5.9, 2.3, [
    "Mean-field theory's known near-Tc bias — physically explained, not hidden",
    "Error trend consistent with textbook mean-field limitations near the Curie point"
  ], { size: 13 });
  img(s, "slide9.png", 7.0, 2.2, 5.65, 3.9);
  footer(s);
}

// ---------------- SLIDE: GIANT-MCE HONESTY FLAG ----------------
{
  let s = bgSlide(false);
  corner(s, false);
  title(s, "Reported Openly", "An Honest Calibration Gap: Giant-MCE \u0394T_ad");
  img(s, "slide10.png", 6.85, 2.2, 5.8, 3.9);
  s.addShape("roundRect", { x: 0.9, y: 2.2, w: 5.6, h: 2.5, rectRadius: 0.08, fill: { color: AMBERBG }, line: { color: AMBER, width: 1.25 } });
  s.addText([
    { text: "Model:  ", options: { bold: true, color: TEXT } }, { text: "24.2 K\n", options: { bold: true, color: AMBER } },
    { text: "Direct measurement (Giguère et al. 1999):  ", options: { bold: true, color: TEXT } }, { text: "10.0 K\n", options: { bold: true, color: AMBER } },
    { text: "Overestimate:  ", options: { bold: true, color: TEXT } }, { text: "2.42×", options: { bold: true, color: AMBER } }
  ], { x: 1.15, y: 2.45, w: 5.1, h: 1.5, fontFace: "Calibri", fontSize: 15, isTextBox: true, margin: 0, lineSpacingMultiple: 1.4 });
  s.addText("Every downstream giant-MCE / cascade number in this talk should be read as an upper bound.", {
    x: 1.15, y: 4.0, w: 5.1, h: 0.6, italic: true, fontFace: "Calibri", fontSize: 12, color: TEXT, isTextBox: true, margin: 0
  });
  footer(s);
}

// ---------------- SLIDE: AMR DEVICE PHYSICS ----------------
{
  let s = bgSlide(false);
  corner(s, false);
  title(s, "System Behavior", "The AMR Characteristic Curve");
  img(s, "slide11.png", 3.9, 2.15, 6.5, 4.35);
  s.addText("Cooling capacity and COP both fall as span widens — and there is a hard structural ceiling near 16 K for a single stage.", {
    x: 1.0, y: 3.6, w: 2.7, h: 2.4, fontFace: "Calibri", fontSize: 13, color: TEXT, isTextBox: true, margin: 0, valign: "top"
  });
  footer(s);
}

// ---------------- SLIDE: LOSS MODELING ----------------
{
  let s = bgSlide(false);
  corner(s, false);
  title(s, "Realism", "Why \u201cIdeal\u201d COP Is the Wrong Number");
  img(s, "slide12.png", 6.7, 2.2, 5.95, 3.5);
  s.addText([
    { text: "W_eddy  \u221D  f \u00B7 H\u00B2\n", options: { fontFace: "Cambria", fontSize: 16, color: DARK, bold: true } },
    { text: "W_pump  \u221D  \u1e41\u00B2\n", options: { fontFace: "Cambria", fontSize: 16, color: DARK, bold: true } },
    { text: "W_base  \u221D  Qc", options: { fontFace: "Cambria", fontSize: 16, color: DARK, bold: true } }
  ], { x: 0.9, y: 2.3, w: 5.4, h: 1.5, isTextBox: true, margin: 0, lineSpacingMultiple: 1.5 });
  s.addText("Calibrated to 3 real devices — an exactly-determined fit, not a regression.", {
    x: 0.9, y: 4.0, w: 5.4, h: 0.8, italic: true, fontFace: "Calibri", fontSize: 13, color: TEXT, isTextBox: true, margin: 0
  });
  footer(s);
}

// ---------------- SLIDE: SYSTEM VALIDATION ----------------
{
  let s = bgSlide(false);
  corner(s, false);
  title(s, "Credibility Check", "System-Level Validation \u2014 6 Real AMR Prototypes");
  img(s, "slide13.png", 3.7, 2.15, 6.9, 3.75);
  s.addText("4 of 6 devices within \u00B115%  \u00B7  2 outliers identified and physically explained", {
    x: 1.0, y: 5.95, w: 11.3, h: 0.4, align: "center", fontFace: "Calibri", fontSize: 14, bold: true, color: DARK, isTextBox: true, margin: 0
  });
  footer(s);
}

// ---------------- SLIDE: 1-D REGENERATOR VS NO-LOAD SPANS ----------------
{
  let s = bgSlide(false);
  corner(s, false);
  title(s, "Going Beyond 0-D", "1-D Transient Model vs. 3 More Prototypes");
  img(s, "graded_bed_validation.png", 3.5, 2.1, 7.4, 3.9);
  s.addShape("roundRect", { x: 0.9, y: 6.05, w: 11.5, h: 0.55, rectRadius: 0.06, fill: { color: AMBERBG }, line: { color: AMBER, width: 1 } });
  s.addText("Direction-inconsistent across devices (undershoots 2 of 3, overshoots 1) \u2014 a genuine, unresolved calibration gap, reported openly rather than smoothed over.", {
    x: 1.15, y: 6.05, w: 11.0, h: 0.55, valign: "middle", italic: true, fontFace: "Calibri", fontSize: 11.5, color: AMBER, isTextBox: true, margin: 0
  });
  footer(s);
}

// ================= PART III: OPTIMIZATION & DESIGN =================
sectionDivider("Part III", "Optimization & Design", [
  "Sobol Discovery", "Surrogate Modeling", "Optimization Setup", "Pareto Front Hero", "Recommended Design Point", "Geometry Optimum", "Multi-Stage Cascades"
]);

// ---------------- SLIDE: SOBOL DISCOVERY ----------------
{
  let s = bgSlide(false);
  corner(s, false);
  title(s, "The #2 Insight", "Frequency, Not Field, Runs the Show");
  img(s, "slide14.png", 6.9, 2.2, 5.75, 3.7);
  const rows = [["Frequency", "0.87"], ["Flow rate", "0.15"], ["Field", "0.01"], ["Regen. effectiveness", "0.01"]];
  let y = 2.3;
  rows.forEach(r => {
    s.addText(r[0], { x: 0.9, y, w: 3.6, h: 0.5, fontFace: "Calibri", fontSize: 14, color: TEXT, isTextBox: true, margin: 0, valign: "middle" });
    s.addText(r[1], { x: 4.6, y, w: 1.5, h: 0.5, fontFace: "Cambria", bold: true, fontSize: 16, color: DARK, isTextBox: true, margin: 0, valign: "middle" });
    y += 0.6;
  });
  s.addText("A naive constant-loss model attributes 100% of sensitivity to a made-up parameter instead.", {
    x: 0.9, y: 4.9, w: 5.4, h: 0.9, italic: true, fontFace: "Calibri", fontSize: 12.5, color: TEXT, isTextBox: true, margin: 0
  });
  footer(s);
}

// ---------------- SLIDE: SURROGATE MODELING ----------------
{
  let s = bgSlide(false);
  corner(s, false);
  title(s, "Scaling Up", "Making 40,000+ Evaluations Possible");
  img(s, "slide15.png", 3.6, 2.15, 7.1, 3.6);
  s.addText("Held-out R\u00B2 = 0.856   \u00B7   RMSE = 365 W   \u00B7   300 training samples   \u00B7   100 held-out test samples   \u00B7   5 design variables", {
    x: 1.0, y: 5.85, w: 11.3, h: 0.5, align: "center", fontFace: "Calibri", fontSize: 13, bold: true, color: DARK, isTextBox: true, margin: 0
  });
  footer(s);
}

// ---------------- SLIDE: OPTIMIZATION SETUP ----------------
{
  let s = bgSlide(false);
  corner(s, false);
  title(s, "NSGA-III", "Optimization Setup");
  card(s, 0.9, 2.15, 5.7, 3.9, "Design Variables (7 + 1 categorical)", "Field \u00B7 frequency \u00B7 flow rate \u00B7 regenerator mass \u00B7 effectiveness \u00B7 blow fraction \u00B7 particle diameter\n\nMaterial family is co-optimized as an 8th, categorical choice.", { dark: true });
  card(s, 6.75, 2.15, 5.7, 3.9, "Objectives", "Maximize COP\n\nMaximize cooling capacity (Qc)\n\nMinimize BOM cost ($)\n\nAcross 3 independent NSGA-III seeds for stability.");
  footer(s);
}

// ---------------- SLIDE: PARETO FRONT HERO ----------------
{
  let s = bgSlide(false);
  corner(s, false);
  title(s, "Hero Finding", "The Pareto Front & a Material Discovery");
  img(s, "slide17.png", 3.55, 2.15, 6.2, 3.5);
  s.addText("La(Fe,Si)\u2081\u2083Hy dominates 67\u2013100% of Pareto-optimal designs across 3 independent seeds \u2014 a genuine, seed-stable finding.", {
    x: 1.5, y: 5.85, w: 10.3, h: 0.6, align: "center", italic: true, fontFace: "Calibri", fontSize: 13, color: DARK, isTextBox: true, margin: 0
  });
  footer(s);
}

// ---------------- SLIDE: RECOMMENDED DESIGN POINT ----------------
{
  let s = bgSlide(false);
  corner(s, false);
  title(s, "Spec Sheet", "Recommended Design Point");
  const specs = [
    ["Material", "La(Fe,Si)\u2081\u2083Hy, tuned, Tc = 285.3 K"],
    ["Field", "1.01 T"], ["Frequency", "2.52 Hz"], ["Flow", "0.49 kg/s"],
    ["Mass", "3.29 kg"], ["Effectiveness", "0.90"], ["Particle diameter", "0.54 mm"]
  ];
  let x = 0.9, y = 2.15, cw = 3.7, ch = 0.85;
  specs.forEach((sp, i) => {
    const col = i % 3, row = Math.floor(i / 3);
    const cx = x + col * (cw + 0.2), cy = y + row * (ch + 0.2);
    s.addShape("roundRect", { x: cx, y: cy, w: cw, h: ch, rectRadius: 0.06, fill: { color: WHITE }, line: { color: "C7D1DE", width: 0.75 } });
    s.addText(sp[0], { x: cx + 0.15, y: cy + 0.08, w: cw - 0.3, h: 0.3, fontFace: "Calibri", fontSize: 10.5, color: MUTED, bold: true, isTextBox: true, margin: 0 });
    s.addText(sp[1], { x: cx + 0.15, y: cy + 0.36, w: cw - 0.3, h: 0.45, fontFace: "Cambria", fontSize: 14, color: DARK, bold: true, isTextBox: true, margin: 0 });
  });
  s.addShape("roundRect", { x: 0.9, y: 5.35, w: 11.55, h: 0.95, rectRadius: 0.08, fill: { color: DARK }, line: { type: "none" } });
  s.addText("COP_electrical = 9.91      \u00B7      Qc = 32.4 kW      \u00B7      Cost = $1,919", {
    x: 0.9, y: 5.35, w: 11.55, h: 0.95, align: "center", valign: "middle", fontFace: "Cambria", fontSize: 19, bold: true, color: WHITE, isTextBox: true, margin: 0
  });
  footer(s);
}

// ---------------- SLIDE: MATERIAL FAMILY COMPARISON ----------------
{
  let s = bgSlide(false);
  corner(s, false);
  title(s, "Eight Candidates, One Winner", "Material Family Comparison at the ASHRAE Point");
  img(s, "material_family_comparison.png", 2.6, 2.05, 8.1, 4.55);
  s.addText("La(Fe,Si)\u2081\u2083Hy leads every span tested (5\u201320K, 2T, 5kg/stage) \u2014 consistent with its Pareto dominance.", {
    x: 1.0, y: 6.6, w: 11.3, h: 0.45, align: "center", italic: true, fontFace: "Calibri", fontSize: 12.5, color: TEXT, isTextBox: true, margin: 0
  });
  footer(s);
}

// ---------------- SLIDE: MAGNET COST MODEL SENSITIVITY ----------------
{
  let s = bgSlide(false);
  corner(s, false);
  title(s, "Robustness Check", "Does the Cost Model Change the Winner?");
  img(s, "magnet_geometry_sensitivity.png", 3.6, 2.1, 7.1, 3.9);
  s.addText("Pareto front size shrinks slightly under a geometric magnet-cost term (21\u219219), but La(Fe,Si)\u2081\u2083Hy still dominates \u2014 the material finding is robust to how magnet cost is modeled.", {
    x: 1.0, y: 6.15, w: 11.3, h: 0.6, align: "center", italic: true, fontFace: "Calibri", fontSize: 12, color: TEXT, isTextBox: true, margin: 0
  });
  footer(s);
}

// ---------------- SLIDE: NANOCOMPOSITE ROBUSTNESS ----------------
{
  let s = bgSlide(false);
  corner(s, false);
  title(s, "Design Exploration", "A Broader Blend Trades Peak Performance for Range");
  img(s, "nanocomposite_robustness.png", 3.5, 2.1, 7.4, 3.6);
  s.addText("At the design span the sharply-tuned single-phase material wins outright \u2014 but off-design, it can collapse to zero cooling while the nanocomposite blend still delivers positive Qc.", {
    x: 1.0, y: 5.9, w: 11.3, h: 0.75, align: "center", italic: true, fontFace: "Calibri", fontSize: 12.5, color: TEXT, isTextBox: true, margin: 0
  });
  footer(s);
}

// ---------------- SLIDE: GEOMETRY OPTIMUM ----------------
{
  let s = bgSlide(false);
  corner(s, false);
  title(s, "Design Tradeoff", "Regenerator Geometry: an Interior Optimum");
  img(s, "slide19.png", 3.4, 2.15, 7.3, 4.15);
  s.addText("Smaller particles improve heat transfer but raise pumping loss. Optimum packed-bed sphere diameter: 0.5 mm", {
    x: 1.0, y: 6.15, w: 11.3, h: 0.5, align: "center", italic: true, fontFace: "Calibri", fontSize: 12.5, color: TEXT, isTextBox: true, margin: 0
  });
  footer(s);
}

// ---------------- SLIDE: MULTI-STAGE CASCADES ----------------
{
  let s = bgSlide(false);
  corner(s, false);
  title(s, "Extending Span", "Staging Beats the Single-Stage Ceiling");
  img(s, "slide20.png", 3.2, 2.15, 7.6, 4.2);
  s.addText("Curie-graded, multi-stage regenerators extend usable span well past the single-stage collapse.", {
    x: 1.0, y: 6.25, w: 11.3, h: 0.4, align: "center", italic: true, fontFace: "Calibri", fontSize: 12.5, color: TEXT, isTextBox: true, margin: 0
  });
  footer(s);
}

// ================= PART IV: RESULTS & HONEST VERDICT =================
sectionDivider("Part IV", "Results & Honest Verdict", [
  "Honest Verdict Hero", "Economics", "Emissions & Water", "Design Roadmap", "Limitations & Open Items", "Conclusion"
]);

// ---------------- SLIDE: HONEST VERDICT HERO ----------------
{
  let s = bgSlide(false);
  corner(s, false);
  title(s, "Hero Finding #2", "The Honest Verdict vs. Conventional Cooling");
  img(s, "slide21.png", 3.6, 1.95, 7.0, 3.35);
  s.addShape("roundRect", { x: 1.6, y: 5.45, w: 10.1, h: 1.55, rectRadius: 0.08, fill: { color: "FCE9EC" }, line: { color: "B23A56", width: 1.25 } });
  s.addText([
    { text: "A systematic search \u2014 span 3\u201330K, vapor-compression \u03B7 = 0.25\u20130.55 \u2014 found ", options: { color: TEXT } },
    { text: "NO combination", options: { bold: true, color: "B23A56" } },
    { text: " where calibrated AMR beats vapor-compression on COP.\n", options: { color: TEXT } },
    { text: "Best AMR COP found: 5.89   vs.   VCC's own worst-case COP: 24.33", options: { bold: true, color: DARK } }
  ], { x: 1.85, y: 5.6, w: 9.6, h: 1.25, fontFace: "Calibri", fontSize: 12.5, isTextBox: true, margin: 0, lineSpacingMultiple: 1.3 });
  footer(s);
}

// ---------------- SLIDE: ECONOMICS ----------------
{
  let s = bgSlide(false);
  corner(s, false);
  title(s, "Techno-Economics", "Economics: Where the Cost Actually Lives");
  img(s, "slide22.png", 6.7, 2.2, 6.0, 3.9);
  statBadge(s, 0.9, 2.2, 5.5, 1.3, "$0.1027/kWh_cooling", "Levelized cost of cooling ($0.0535 capital + $0.0492 electricity, 15-yr, 6% discount)", DARK);
  statBadge(s, 0.9, 3.75, 5.5, 1.3, "$5,825", "Materials BOM (2T, 5kg Gd): Magnet $5,388 (134.7kg) \u00B7 MCM $100 \u00B7 Yoke $337 (67.4kg)", "B08A2B");
  s.addText("Magnet mass is the dominant cost driver \u2014 not electronics, not electricity.", {
    x: 0.9, y: 5.2, w: 5.5, h: 0.9, italic: true, fontFace: "Calibri", fontSize: 12, color: TEXT, isTextBox: true, margin: 0
  });
  footer(s);
}

// ---------------- SLIDE: EMISSIONS & WATER ----------------
{
  let s = bgSlide(false);
  corner(s, false);
  title(s, "Told Honestly", "Emissions & Water: the Refrigerant-Free Case");
  img(s, "slide23.png", 3.5, 2.15, 7.3, 3.6);
  s.addText([
    { text: "Emissions: ", options: { bold: true, color: "B23A56" } },
    { text: "AMR runs ~3.6\u00D7 higher than VCC (0.989 vs 0.275 tCO\u2082e/yr) \u2014 operational electricity dominates.      ", options: { color: TEXT } },
    { text: "Water: ", options: { bold: true, color: DARK } },
    { text: "AMR uses ~36\u00D7 less water than VCC (0.05 vs 1.80 L/kWh) \u2014 a genuine advantage today.", options: { color: TEXT } }
  ], { x: 1.0, y: 5.9, w: 11.3, h: 0.9, align: "center", fontFace: "Calibri", fontSize: 12.5, isTextBox: true, margin: 0, lineSpacingMultiple: 1.3 });
  footer(s);
}

// ---------------- SLIDE: DESIGN ROADMAP ----------------
{
  let s = bgSlide(false);
  corner(s, false);
  title(s, "Ranked, Not Guessed", "What Would Actually Close the Gap");
  const items = [
    ["1", "Operating frequency", "Dominant lever (Sobol S_T = 0.87)", 0.95],
    ["2", "Material / composition choice", "La(Fe,Si)\u2081\u2083Hy: +44% COP vs. plain Gd", 0.75],
    ["3", "Curie-grading / staging", "Extends usable span past the single-stage ceiling", 0.55],
    ["4", "Regenerator geometry", "Interior optimum at 0.5mm particle size", 0.35],
    ["5", "Field / flow balance", "Co-optimized, not maximized independently", 0.2]
  ];
  let y = 2.2;
  items.forEach(it => {
    s.addShape("ellipse", { x: 0.9, y, w: 0.5, h: 0.5, fill: { color: DARK }, line: { type: "none" } });
    s.addText(it[0], { x: 0.9, y, w: 0.5, h: 0.5, align: "center", valign: "middle", fontFace: "Cambria", bold: true, fontSize: 16, color: WHITE, isTextBox: true, margin: 0 });
    s.addText(it[1], { x: 1.6, y: y - 0.02, w: 4.0, h: 0.5, fontFace: "Calibri", bold: true, fontSize: 13.5, color: DARK, isTextBox: true, margin: 0, valign: "middle" });
    s.addText(it[2], { x: 5.7, y: y - 0.02, w: 3.7, h: 0.5, fontFace: "Calibri", fontSize: 11.5, color: TEXT, isTextBox: true, margin: 0, valign: "middle" });
    s.addShape("rect", { x: 9.55, y: y + 0.13, w: 2.9, h: 0.22, fill: { color: "D6DEE8" }, line: { type: "none" } });
    s.addShape("rect", { x: 9.55, y: y + 0.13, w: 2.9 * it[3], h: 0.22, fill: { color: PINK }, line: { type: "none" } });
    y += 0.72;
  });
  footer(s);
}

// ---------------- SLIDE: MAGNETOCALORIC FLUIDS ----------------
{
  let s = bgSlide(false);
  corner(s, false);
  title(s, "Design Exploration", "Magnetocaloric Fluids: an Alternative Working Body");
  img(s, "fluid_mce_volume_fraction.png", 6.7, 2.1, 5.95, 3.9);
  s.addText([
    { text: "Ferrofluid / MR-suspension bed, swept over volume fraction \u03C6.\n\n", options: { color: TEXT } },
    { text: "COP_elec peaks near \u03C6 \u2248 0.10", options: { bold: true, color: DARK } },
    { text: " \u2014 too little MCE mass hurts capacity, too much raises viscosity/pumping loss.\n\n", options: { color: TEXT } },
    { text: "Design-exploration tool, not yet validated against a real device.", options: { italic: true, color: MUTED } }
  ], { x: 0.9, y: 2.3, w: 5.4, h: 3.6, fontFace: "Calibri", fontSize: 12.5, isTextBox: true, margin: 0, lineSpacingMultiple: 1.3 });
  footer(s);
}

// ---------------- SLIDE: THERMAL DIODE SENSITIVITY ----------------
{
  let s = bgSlide(false);
  corner(s, false);
  title(s, "Cost-Only Study", "Mechanical Thermal Diode: an Upper-Bound Check");
  img(s, "thermal_diode_sensitivity.png", 3.6, 2.1, 7.1, 3.6);
  s.addShape("roundRect", { x: 0.9, y: 5.9, w: 11.55, h: 0.85, rectRadius: 0.07, fill: { color: AMBERBG }, line: { color: AMBER, width: 1 } });
  s.addText("Illustrative actuation cost reduces COP_electrical by \u2264 0.03% across 0.5\u20138 Hz \u2014 a small, cost-only effect; no offsetting heat-transfer benefit is modeled here.", {
    x: 1.15, y: 5.9, w: 11.05, h: 0.85, valign: "middle", italic: true, fontFace: "Calibri", fontSize: 11.5, color: AMBER, isTextBox: true, margin: 0
  });
  footer(s);
}

// ---------------- SLIDE: LIMITATIONS & OPEN ITEMS (UPDATED) ----------------
{
  let s = bgSlide(false);
  corner(s, false);
  title(s, "Told Straight", "Limitations, Future Scope & Open Items");
  s.addShape("roundRect", { x: 0.9, y: 1.95, w: 11.5, h: 0.55, rectRadius: 0.06, fill: { color: AMBERBG }, line: { color: AMBER, width: 1 } });
  s.addText("Future scope and unresolved items deliberately documented — a roadmap for the next hardware build.", {
    x: 1.15, y: 1.95, w: 11.0, h: 0.55, valign: "middle", italic: true, fontFace: "Calibri", fontSize: 12, color: AMBER, isTextBox: true, margin: 0
  });
  bulletList(s, 0.9, 2.7, 11.5, 4.4, [
    "Purpose-built experimental prototype and hardware-in-the-loop validation needed for data-center-relevant spans and loads.",
    "Resolve the open 1-D regenerator calibration gap (directionally-inconsistent axial-conduction fit across benchmark devices).",
    "Magnetocaloric-fluid and passive/hybrid-regenerator architectures are mature design-exploration tools, but need to become validated features.",
    "Full bottom-up manufactured-system BOM (pumps, motors, controls) needed to replace the current order-of-magnitude cost multiplier.",
    "Targeted OCR/access required to mine remaining chapters in Halbach design, AMR cycle topology, and thermal-diode heat switching.",
    "AI-assisted or digital-twin real-time optimisation for adaptive data-center cooling control."
  ], { size: 12.5 });
  footer(s);
}

// ---------------- SLIDE: CONCLUSION ----------------
{
  let s = bgSlide(true);
  corner(s, true);
  s.addText("CONCLUSION", { x: 0.9, y: 0.6, w: 11.3, h: 0.7, fontFace: "Cambria", fontSize: 32, bold: true, color: WHITE, isTextBox: true, margin: 0 });
  const paras = [
    "Gd model validated to +4.2% / \u22121.1% / \u22128.1% vs. Dan\u2019kov et al. (1998)",
    "System model validated against 6 real AMR prototypes \u2014 4 of 6 within \u00B115%",
    "Sobol sensitivity: electrical COP dominated by frequency (S_T = 0.87), not field",
    "NSGA-III discovered a La(Fe,Si)\u2081\u2083Hy-dominated Pareto front (67\u2013100% share, seed-stable)",
    "No design point found where calibrated AMR beats vapor-compression on COP or emissions",
    "Refrigerant-free operation + ~36\u00D7 lower water use are real advantages today"
  ].map(t => ({ text: t, options: { bullet: { code: "25AA", color: PINK }, color: "E7F2F0", fontSize: 15, fontFace: "Calibri", breakLine: true, paraSpaceAfter: 14 } }));
  s.addText(paras, { x: 0.9, y: 1.65, w: 11.3, h: 4.3, isTextBox: true, margin: 0, valign: "top" });
  s.addText("\u201cA validated, self-checking engineering framework that tells you exactly where it stands today, why, and what would have to change.\u201d", {
    x: 0.9, y: 6.2, w: 11.3, h: 0.8, italic: true, fontFace: "Cambria", fontSize: 13, color: PINK, isTextBox: true, margin: 0
  });
}

pres.writeFile({ fileName: "MagCool_DC.pptx" }).then(() => console.log("done"));