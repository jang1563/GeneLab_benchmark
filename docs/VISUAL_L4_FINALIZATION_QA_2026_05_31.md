# Visual L4 Finalization QA

Date: 2026-05-31

Purpose: record the first L4-oriented QA pass after P0 figure restructuring,
table separation, and jargon reduction.

## What Was Added

New audit/render scripts:

- `scripts/audit_premium_visual_sources.py`
- `scripts/render_premium_manuscript_previews.py`

New QA outputs:

- `output/premium_audit/premium_visual_numeric_audit.csv`
- `output/premium_audit/premium_visual_numeric_audit.json`
- `output/premium_figures/manuscript_previews/manuscript_preview_manifest.csv`
- `output/premium_figures/manuscript_previews/manuscript_preview_manifest.json`
- `output/premium_figures/manuscript_previews/manuscript_preview_contact_sheet.png`

## Numeric Source Audit

Command:

```bash
python3 scripts/audit_premium_visual_sources.py
```

Result:

- 117 checks run.
- 117 checks passed.
- 0 failures.

One true rounding issue was found and fixed:

- `Mouse-Geneformer:skin` delta changed from `-0.265` to `-0.264`, matching
  `0.557 - 0.821` in `evaluation/RESULTS_SUMMARY.md`.

Interpretation:

- Current regenerated figure/table numbers match their cited source files.
- This does not mean the source files are publication-frozen; it only verifies
  that the visual layer is numerically faithful to current sources.

## Manuscript-Width Preview

Command:

```bash
MPLCONFIGDIR=output/.matplotlib MPLBACKEND=Agg \
  python3 scripts/render_premium_manuscript_previews.py
```

Preview setting:

- two-column width: 7.2 inches;
- render DPI: 300;
- output width: 2160 px.

Contact sheet:

- `output/premium_figures/manuscript_previews/manuscript_preview_contact_sheet.png`

## Figure-Level Readiness

| Figure | Numeric audit | 2-column preview | Current status |
|---|---|---|---|
| Fig1 tissue transfer | PASS | PASS | strongest clean main-figure candidate |
| Fig2 pathway | PASS | PASS/BORDERLINE | strong main candidate; may benefit from manuscript-specific panel sizing |
| Fig3 model | PASS | BORDERLINE | deck-ready; manuscript variant needed |
| Fig4 platform architecture | PASS | BORDERLINE | resource overview/schematic; pair with tables |
| Fig5 release boundary | PASS | SUPPLEMENT/STATUS | not a main science figure |
| Fig6 organoid | PASS | PASS/BORDERLINE | strong extension figure; caption must keep biology checks secondary |

## Strict Interpretation

### Fig1

This is the cleanest current manuscript candidate. It is now mostly one plot,
has direct tissue ranking, and no longer embeds a table-like readout panel.

Remaining work:

- caption definition for AUROC and held-out mission;
- final source freeze.

### Fig2

Scientific story remains strong. The two-column preview is usable, but Panel A
and Panel C are dense enough that a manuscript-specific version may be better.

Remaining work:

- caption must clarify that mission, hardware, and gravity label prediction is
  a diagnostic check with coupled factors;
- consider splitting into two panels for a paper version:
  - unwanted mission-label signal;
  - pathway task-rescue and activity agreement.

### Fig3

The current figure is still deck-first. It survives as a two-column preview but
Panel C remains too text-heavy for a high-profile paper main figure.

Required before manuscript use:

- create a manuscript variant with a dot plot rather than wide bars;
- convert coverage warning into a compact table or caption;
- consider an axis starting at zero or explicitly mark the zoomed range.

### Fig4

The table separation helped, but this remains a resource/platform overview.
The lane text is small at manuscript width.

Recommended use:

- resource paper overview figure or deck slide;
- pair with `table4_v9_lane_readiness` and `table4_v9_evidence_footprint`.

### Fig5

This is correctly demoted from main scientific figure to status/release
schematic. It is useful, but not a main result.

Recommended use:

- release deck;
- supplement;
- Methods/Data availability governance graphic.

### Fig6

The decision table removal helped. The figure now works better as a cautious
extension figure, although Panel A/B text is still dense.

Remaining work:

- caption must state that flight-response and model-gene-effect checks are
  secondary biology checks;
- consider moving Panel A footprint details into the supplement if a shorter
  manuscript figure is needed.

## Current L4 Verdict

Not all figures are L4 yet.

Current best L4 path:

1. Promote Fig1 after caption/source freeze.
2. Build a manuscript-specific Fig2 variant.
3. Build a manuscript-specific Fig3 variant.
4. Keep Fig4 as resource overview, not primary result.
5. Keep Fig5 outside main scientific figure set.
6. Promote Fig6 after caption hardening and possible layout tightening.

## Next Recommended Figure Work

Priority order:

1. Fig2 manuscript variant.
2. Fig3 manuscript variant.
3. Caption pack for Fig1/Fig2/Fig3/Fig6.
4. Grayscale/colorblind preview.
5. Final source freeze table.

