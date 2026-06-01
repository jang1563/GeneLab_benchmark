# Visual Premium Quality Standard

Date: 2026-05-31

Purpose: define the quality bar for upgrading GeneLab Benchmark /
SpaceBio-Bench v1-v9 visual assets into a premium presentation and
paper-ready visual system.

## North Star

Premium visuals for this project should feel like a high-end scientific
resource/benchmark paper, not like a collection of exploratory analysis plots.

The visual system must communicate three things at once:

1. Scientific credibility: every number has a traceable source.
2. Benchmark clarity: tasks, splits, metrics, and claims are visually obvious.
3. Editorial polish: figures are calm, legible, consistent, and presentation
   ready without live interactivity.

## Non-Negotiable Standard

A figure is premium-ready only if it passes all of these:

- The main message is understandable in 5 seconds.
- The figure remains legible when projected and when printed in a manuscript.
- Every plotted number has a source-of-truth file and captionable context.
- Panel hierarchy is deliberate: one lead panel, supporting panels, no filler.
- Typography, color, spacing, and stroke weights are consistent across figures.
- Color is semantically stable across the deck/paper.
- The figure is static-exportable without internet access.
- The figure does not rely on hover, buttons, or hidden interactive state.
- Claim boundaries are visible where needed, especially v8/v9.
- It can survive reviewer scrutiny without sounding broader than the evidence.

## Quality Levels

| Level | Meaning | Use |
|---|---|---|
| L0 raw | Exploratory output, not presentation-ready | Internal only |
| L1 readable | Correct and viewable, but visually inconsistent | Lab update / backup |
| L2 deck-ready | Legible, clean, and captioned | Talks and progress decks |
| L3 manuscript-ready | Static, source-verified, polished, journal-safe | Paper main/supplement |
| L4 premium | Editorially cohesive, beautiful, source-verified, reusable across deck and paper | Main story figures |

Target:

- Main deck figures: L3-L4.
- Main manuscript figures: L4.
- Supplemental figures: L2-L3.
- Diagnostic v8/v9 status visuals: L2-L3, with explicit claim boundaries.

## Review Rubric

Score each asset from 0 to 3 in each category.

| Category | 0 | 1 | 2 | 3 |
|---|---|---|---|---|
| Message clarity | unclear | visible but cluttered | clear with caption | obvious without explanation |
| Data trust | unknown source | source implied | source listed | canonical source verified |
| Static reproducibility | no export path | screenshot only | static export possible | scripted/static export reproducible |
| Typography | inconsistent/tiny | readable only large | mostly consistent | publication-grade |
| Layout | crowded/empty/imbalanced | usable | balanced | editorially composed |
| Color semantics | arbitrary | partly meaningful | consistent within figure | consistent across project |
| Accessibility | poor contrast/color-only | partially accessible | readable grayscale-ish | colorblind/print robust |
| Label economy | over-labeled or cryptic | acceptable | clean | minimal but self-contained |
| Claim boundary | missing | caption only | visible note | integrated visual status language |
| Venue fit | exploratory | slide only | slide + supplement | main-figure ready |

Pass thresholds:

- Main figure candidate: average >= 2.6 and no category below 2.
- Slide figure candidate: average >= 2.2 and no critical category below 1.
- Anything with claim-boundary score 0 cannot be used for v8/v9 main slides.

## Design System

### Typography

Use one sans-serif family throughout. Recommended hierarchy:

- Figure title: 18-24 pt in slides, 10-12 pt in manuscript panels.
- Panel title: 13-16 pt in slides, 8-10 pt in manuscript panels.
- Axis labels: 11-13 pt in slides, 7-9 pt in manuscript panels.
- Tick labels: 9-11 pt in slides, 6-8 pt in manuscript panels.
- Callouts: short, bold, never paragraph-like.

Avoid:

- Tiny rotated labels unless unavoidable.
- Multiple font families.
- Long explanatory text inside panels.

### Text and Jargon

Figure-visible text must be reader-facing first and method-facing second.

Rules:

- Titles should state the biological or benchmark message in plain language.
- Panel titles should describe the comparison, not the internal pipeline step.
- Acronyms are not allowed in titles unless they are globally recognized by the
  target venue or are part of the project name.
- Metric acronyms such as AUROC, FDR, and macro-F1 may appear in axis labels or
  legends when they are the plotted quantity.
- Internal shorthand such as LOMO, NES, FM, RRRM, payload, diagnostic-only,
  feature-effect, and top-k should be moved to captions, manifests, or source
  tables unless there is no accurate plain-language substitute.
- A first-pass plain-language replacement should be attempted before accepting
  any jargon.

Preferred replacements:

| Internal term | Figure-facing term |
|---|---|
| LOMO | held-out mission |
| NES concordance | pathway activity agreement |
| FM | foundation model |
| DE rows | gene rows or significant gene rows |
| top-k | Top 50 / Top 100 / selected top genes |
| diagnostic-only | biology check or secondary check |
| feature-effect | model gene effects |
| provenance | evidence trail |
| payload | data files or local data copy |
| metadata alpha | metadata-only preview, with exact release term in caption |

### Color Semantics

Use color as a contract, not decoration.

Recommended semantic palette:

- Flight / spaceflight: high-saturation coral or red-orange.
- Ground / control: neutral gray or muted blue-gray.
- Pathway / biology abstraction: teal/green.
- Classical ML: deep blue.
- Foundation models / LLMs: violet or magenta, used sparingly.
- v8 hypothesis-only: amber.
- v9 platform/provenance: slate/teal.
- Blocker / not claimed: muted red.
- Ready / passed: green.
- Diagnostic / draft: amber.

Rules:

- Do not use rainbow palettes for categories unless there is no semantic
  ordering.
- Use diverging palettes only for signed quantities.
- Use sequential palettes only for magnitude.
- Keep tissue colors stable when tissues appear across figures.
- Make status colors consistent: pass, diagnostic, blocker, not claimed.

Current global color contract:

| Meaning | Color family | Notes |
|---|---|---|
| Tissue highlight / lead result | muted burnt orange | Used sparingly for a biological lead such as thymus |
| Classical baseline / evidence | deep blue | Do not reuse for status unless clearly separated |
| Gene-level input | blue | Within pathway figures only |
| Pathway / biology summary | teal/green | Within pathway figures and biology checks |
| Ready / passed | green | Status meaning only |
| Draft / needs update | amber | Status meaning only |
| Blocked / not claimed | muted red | Status meaning only |
| Foundation / LLM model family | violet/magenta | Model comparison only, restrained use |

### Layout

Premium figures should use:

- A clear lead panel.
- Fewer panels with stronger hierarchy.
- Consistent outer margins.
- Shared legends where possible.
- Direct labels instead of long legends when space allows.
- Journal-style dot plots, lollipop plots, thin-rule tables, heatmaps,
  schematics, and axis-based summaries rather than UI-like cards.
- No empty panels.
- No decorative backgrounds.
- No unneeded gridlines.
- No UI buttons in exported figures.
- No dashboard-style card boxes for scientific figures; they read as low-end
  product UI rather than high-profile journal visuals.

### Figure/Table Boundary

Do not make a table pretend to be a figure.

Use a figure when the visual encodes:

- spatial, temporal, workflow, or causal structure;
- quantitative patterns that benefit from axes or geometry;
- model/tissue/task comparisons where shape, rank, or uncertainty matters;
- schematic relationships that are easier to see than read.

Use a table when the content is primarily:

- lane readiness rows;
- release gates;
- source inventories;
- allowed/prohibited language;
- decision matrices;
- exact summary rows.

Decks may contain status slides with table-like content, but manuscript figure
sets should split these into a schematic figure plus a real table or
supplementary table.

For slides:

- Prefer 16:9 landscape figures.
- One major point per slide.
- Use "evidence + claim boundary" instead of dense multi-panel overload.
- Treat the slide as a layered scientific scene, not a pasted chart.

For manuscript:

- Use fewer colors and more whitespace.
- Captions carry methodological detail; panels carry evidence.
- Extended data can hold dense support panels.

### Layered Scene Standard For Premium Slides

The slide-deck version of a figure should use depth deliberately. The goal is
not decoration; the goal is to help the audience see evidence, interpretation,
and trust boundaries in the right order.

Production rule:

> Layered PNG scene + editable scientific interpretation.

This replaces the weaker rule of "PNG proof plate plus text shell." A premium
slide should have a source-verified proof object embedded inside a composed
scene, with interpretation and caveats remaining editable.

Use this z-stack for major 2026 ISGP-grade slide figures:

| Layer | Role | Recommended assets | Rules |
|---|---|---|---|
| Z0 canvas / atmosphere | Set the scientific context and visual depth | full-canvas dark field, source-image crop, cellular texture, evidence-paper texture | Quiet, low contrast, never compete with the data |
| Z1 measurement layer | Show that the slide is about measured biology, not generic design | faint grid, orbital arc, mission timeline tick, assay lane, tissue/sample lane | Thin, geometric, aligned to the evidence object |
| Z2 evidence surface | Carry the source-verified result | high-resolution PNG/PDF proof object, source-derived chart, microscope/object crop, backplate/shadow | Must be numerically audited or traceable to a source manifest |
| Z3 interpretation layer | Tell the audience where to look | editable callout, focus ring, highlight path, short plain-language phrase | Editable in the slide tool; no paragraph text |
| Z4 trust/status layer | State scope without killing the visual | small caveat, data source label, preliminary/internal tag when appropriate | Reader-facing wording only; no raw file paths or internal operating notes |
| Z5 motion/focus layer | Create one narrative movement | trajectory line, crop window, transfer break, before/after reveal path | One dominant motion idea per slide |

Layering constraints:

- Z0-Z1 should create depth but remain nearly silent.
- Z2 is the proof plate; it can be rasterized if it is source-verified.
- Z3-Z5 should stay editable for slide iteration and venue-specific wording.
- Do not use decorative blobs, generic gradients, or UI cards as fake depth.
- Do not let Z4 reintroduce internal language such as release-decision notes,
  source file paths, or review instructions.
- For manuscript exports, collapse the slide scene back into a clean static
  figure or move atmospheric/context layers into the cover/graphical abstract.

Depth cues that are acceptable:

- low-opacity tissue/cell texture behind the data;
- thin orbital or timeline arcs tied to the mission-held-out story;
- subtle backplate/shadow around a proof object;
- masked crop windows that reveal a local result;
- one focus ring or highlight path that matches the verbal claim.

Depth cues to avoid:

- boxed dashboard cards;
- unmotivated shadows around every object;
- cinematic backgrounds that make the plotted data feel secondary;
- more than one competing trajectory line;
- small text baked into a PNG when it should be editable.

## Project-Specific Figure Rules

### v1-v7 Benchmark Core

Main figure style:

- Calm benchmark paper aesthetic.
- AUROC and metric plots should use consistent 0.5 chance reference line.
- Confidence intervals should be visible but not visually heavy.
- Tissue ranking should be explicit.
- FM/LLM comparisons must label subset coverage.

Required boundary:

- Do not imply all FM rows are the same 8-tissue surface.
- Keep v4 canonical table and raw best-row table separate until reconciled.

### v8 SpaceMed

Main figure style:

- Translational incubator / hypothesis map.
- Use amber "hypothesis-only" labels on all intervention/Mars slides.
- Avoid clinical/countermeasure styling that implies recommendation.

Required boundary:

- No intervention plot can appear without "hypothesis-only" language.
- Mars extrapolation should be backup unless the visual is reframed as
  sensitivity/exploratory flagging.

### v9 SpaceBio-Bench Platform

Main figure style:

- Platform architecture and release-boundary schematics.
- Use manifest/registry/evaluator/provenance flow diagrams, thin-rule status
  tables, compact heatmaps, and lollipop evidence summaries.
- Avoid card-box dashboards. Status can be shown with row rules, dots, glyphs,
  and restrained color, not tiled UI panels.

Required boundary:

- Public bulk alpha is metadata-only, not payload-frozen.
- Organoid and multispecies tracks are draft/diagnostic unless explicitly
  promoted later.
- Single-cell track is not runnable until canonical payload staging or
  regeneration is complete.

## Preferred Upgrade Strategy

Do not try to make every existing HTML figure premium.

Decision after JK review:

> If an existing visual asset is weak, cluttered, technically brittle, or not
> aligned with the premium narrative, rebuild it from scratch. Existing figures
> are evidence/design references, not obligations.

Preferred route:

1. Select 8-10 main-story figures.
2. Regenerate them in a unified static system from canonical source files.
3. Keep older HTML figures as reference/backup/supplement candidates.
4. Export each final figure as PNG, PDF, and source-data-linked manifest.

Suggested premium main-story figure set:

1. Mission-held-out benchmark design.
2. Tissue hierarchy: thymus vs liver.
3. Pathway rescue and batch resistance.
4. Model tier comparison.
5. v4 hardening grid.
6. Temporal/single-cell failure-mode extension.
7. v8 translational incubator, hypothesis-only.
8. v9 platform architecture and status matrix.

## Session Operating Protocol

Each visual-premiumization session should do one bounded block:

1. Pick 1-3 assets or one figure family.
2. Confirm source-of-truth numbers.
3. Decide keep / regenerate / retire.
4. Draft or regenerate the visual.
5. Render and visually inspect.
6. Record the QA result and next action.

Never mix:

- figure design;
- numeric reconciliation;
- claim-boundary decisions;
- slide-deck formatting;

unless the scope is deliberately small. The project is big enough that mixing
all four causes drift.

## Output Convention

Recommended future output layout:

- `output/visual_premium_qa/`: screenshots, render checks, visual QA notes.
- `output/premium_figures/`: regenerated static figure exports.
- `docs/VISUAL_*`: decisions, standards, and audits.

Each final figure should have:

- `.png` for slides.
- `.pdf` or `.svg` for vector/paper use.
- `.json` or `.csv` source manifest.
- one-line claim boundary.
- one-line source-of-truth statement.
