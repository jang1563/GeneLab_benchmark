# Visual Text and Jargon Audit

Date: 2026-05-31

Scope: figure-visible text in the regenerated premium figure family under
`output/premium_figures/`.

Purpose: make the figures understandable to a broad computational biology,
space-biology, and benchmark-paper audience without hiding the technical
methods needed for reviewer trust.

## Audit Rule

Figure titles and panel titles should state the biological or benchmark message
in plain language. Technical terms may remain in axis labels, legends, footers,
source manifests, and captions when they are the actual metric or release term.

Working rule:

- Title: reader-facing claim.
- Panel title: reader-facing comparison.
- Axis label: metric name is allowed.
- Footer/caption: method boundary and source detail.
- Manifest: exact internal terminology and source provenance.

## P0 Changes Applied

### Figure 1

Current primary asset:

- `output/premium_figures/premium_fig1_tissue_transfer_hierarchy.png`
- `output/premium_figures/premium_fig1_tissue_transfer_hierarchy.pdf`

Legacy path also regenerated for continuity:

- `output/premium_figures/premium_fig1_core_tissue_pathway.png`
- `output/premium_figures/premium_fig1_core_tissue_pathway.pdf`

Changes:

- Split Figure 1 into a tissue-transfer-only figure.
- Removed detection and pathway-rescue panels from the core visual.
- Replaced the multi-claim bar layout with a journal-style dot/interval plot
  plus a thin-rule summary table.
- Replaced `cross-mission` and `LOMO` language in visible figure text with
  `held-out mission`.
- Kept `AUROC` because it is the plotted metric, but moved scope details into
  the subtitle/footer.

Status:

- P0 structural issue resolved.
- Current figure is a much cleaner main-deck/manuscript candidate.
- Final numeric audit still needed before L4.

### Figure 2

Current asset:

- `output/premium_figures/premium_fig2_pathway_artifact_rescue.png`
- `output/premium_figures/premium_fig2_pathway_artifact_rescue.pdf`

Changes:

- Removed the fitted trend line from Panel C.
- Replaced the risky `Rank r=0.9 excluding gastrocnemius outlier` callout with
  a plain exploratory note.
- Replaced `NES concordance` in the axis label with `Pathway activity
  agreement`; the exact NES source remains in the manifest.
- Shortened panel titles to prevent overlap and reduce jargon.
- Replaced broad `artifact abstraction` phrasing with reader-facing language
  about unwanted mission signals and selected weak tasks.

Status:

- P0 statistical-visualization risk resolved for deck use.
- Manuscript caption still needs to state that mission, hardware, and gravity
  labels are coupled diagnostic checks, not biological endpoint metrics.

### Figure 6

Current asset:

- `output/premium_figures/premium_fig6_organoid_diagnostic_surface.png`
- `output/premium_figures/premium_fig6_organoid_diagnostic_surface.pdf`

Changes:

- Rebuilt Panel B as three separate mini-axes:
  - primary task;
  - flight-response pattern;
  - model gene effects.
- Removed the single shared bar axis that made AUROC, macro-F1, direction
  match, and rank correlation look like one leaderboard.
- Replaced `diagnostic-only`, `feature-effect`, `response signature`, `DE`,
  and `top-k` where possible with simpler figure-facing language.
- Kept metric names such as `AUROC`, `Macro-F1`, and `rank correlation` where
  the figure is explicitly reporting those metrics.

Status:

- P0 metric-family comparison issue resolved for deck use.
- For a manuscript main figure, Panel D may still move to supplement.

## Term-Level Decisions

| Term | Problem in figure text | Current decision |
|---|---|---|
| LOMO | Internal shorthand; opaque to many readers | Replace with `held-out mission` in figure-visible text |
| NES | Method-specific pathway statistic | Replace with `pathway activity agreement` in axes/titles; retain in manifest/caption |
| FM | Abbreviation not universally clear | Replace with `foundation model` or `single-cell model` |
| DE | Common but still abrupt in titles | Use `gene contrasts` or `significant gene rows` in figures; retain DE in captions/source files |
| top-k | ML jargon in visible label | Use `Top 50`, `Top 100`, etc. |
| diagnostic-only | Reads as internal gating language | Use `biology check`, `secondary check`, or `show as secondary` |
| feature-effect | Method jargon | Use `model gene effects` |
| response signature | Can be obscure in title | Use `flight-response pattern` |
| provenance | Resource-paper term but not always intuitive | Use `evidence trail` in figure tables; retain provenance in docs |
| payload | Engineering/release jargon | Use `data files` or `local data copy` in figure text |
| metadata alpha | Release jargon | Use `metadata-only preview` in most figure text; retain exact alpha wording where claim boundary needs it |
| leaderboard | Retain only when explicitly warning against leaderboard interpretation | Allowed in claim-boundary warnings |

## Current Figure-Visible Text Risk

| Figure | Current risk | Remaining action |
|---|---|---|
| Figure 1 tissue transfer | Low | Keep AUROC but explain once in caption/deck notes |
| Figure 2 pathway | Medium-low | Caption must clarify diagnostic label tasks and coupled mission/hardware/gravity factors |
| Figure 3 model tiers | Medium | Manuscript variant should reduce bar-chart ranking feel and possibly start axis at zero |
| Figure 4 v9 platform | Medium | Good for deck/resource overview; too operational for a pure result figure |
| Figure 5 public bulk boundary | Medium | Strong governance figure; do not use as main scientific result |
| Figure 6 organoid | Medium-low | Keep biology-check wording; consider moving decision table to supplement |

## Figure Text Contract Going Forward

Required before a figure can be marked L4:

1. No unexplained acronym in title or panel title.
2. No internal pipeline shorthand in a main visual unless it is the study name.
3. No one-axis comparison of different metric families.
4. No figure text that implies final release, frozen benchmark, clinical
   recommendation, or leaderboard status beyond the actual evidence.
5. Captions must carry method detail; panels must carry message and evidence.

