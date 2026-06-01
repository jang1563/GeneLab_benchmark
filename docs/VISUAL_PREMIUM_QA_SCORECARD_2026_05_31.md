# Visual Premium QA Scorecard

Date: 2026-05-31

Purpose: score regenerated premium visual prototypes against
`docs/VISUAL_PREMIUM_QUALITY_STANDARD_2026_05_31.md`.

Strict follow-up review:

- `docs/VISUAL_EXPERT_PANEL_REVIEW_2026_05_31.md`

Use the expert-panel review as the stricter gate for future L4 decisions. This
scorecard is an internal QA record; it should not override expert-panel
major-revision or reject-for-main-figure verdicts.

P0 update:

- Figure 1 has been replaced by the cleaner tissue-transfer hierarchy asset.
- Figure 2 Panel C has been hardened by removing the fitted trend line.
- Figure 6 Panel B now separates metric families into mini-axes.
- Figure-visible jargon audit:
  `docs/VISUAL_TEXT_JARGON_AUDIT_2026_05_31.md`.
- Figure/table separation audit:
  `docs/VISUAL_FIGURE_TABLE_SEPARATION_2026_05_31.md`.
- L4 finalization QA:
  `docs/VISUAL_L4_FINALIZATION_QA_2026_05_31.md`.
- Numeric source audit result: 117/117 checks passed.

## Scoring Scale

0 = fail, 1 = weak, 2 = acceptable, 3 = premium.

Main-figure target:

- Average >= 2.6.
- No critical category below 2.

## Assets Scored

| Asset | Current level | Mean score | Verdict |
|---|---:|---:|---|
| `output/premium_figures/premium_fig1_tissue_transfer_hierarchy.png` | L3 prototype | 2.7 | Cleaner after P0 split; numeric audit still needed |
| `output/premium_figures/premium_fig2_pathway_artifact_rescue.png` | L3 prototype | 2.7 | Stronger premium candidate |
| `output/premium_figures/premium_fig3_model_tier_comparison.png` | L3 prototype | 2.7 | Strong story figure if mixed-surface boundary stays explicit |
| `output/premium_figures/premium_fig4_v9_platform_architecture.png` | L3 prototype | 2.8 | Stronger as schematic; readiness rows now tables |
| `output/premium_figures/premium_fig5_public_bulk_release_boundary_schematic.png` | L2/L3 status schematic | 2.5 | Release/status schematic; detailed gates now tables |
| `output/premium_figures/premium_fig6_organoid_diagnostic_surface.png` | L3 prototype | 2.8 | Strong organoid extension figure after moving decision table out |

## Figure 1: Tissue-Transfer Hierarchy Prototype

Asset:

- `output/premium_figures/premium_fig1_tissue_transfer_hierarchy.png`

Strengths:

- Static PNG/PDF with source manifest.
- Clean lead panel showing held-out mission tissue transfer.
- Good use of chance reference line and confidence intervals.
- Claim boundary included.
- Much more coherent than the previous three-claim composite.

Weaknesses:

- Final numeric audit still needed.
- Summary panel is deck-friendly; manuscript version may use a smaller inset or
  caption-only interpretation.
- AUROC remains in the axis label and should be defined once in the caption.

Scores:

| Category | Score | Note |
|---|---:|---|
| Message clarity | 3 | Lead panel is immediately clear |
| Data trust | 2 | Manifest exists; final numeric audit still needed |
| Static reproducibility | 3 | PNG/PDF/script/manifest all exist |
| Typography | 2 | Mostly good, title long |
| Layout | 3 | P0 split removed competing right panels |
| Color semantics | 2 | Good start; needs project-wide tissue/status palette |
| Accessibility | 2 | Contrast OK; color-only differences remain |
| Label economy | 2 | Some labels could be simplified |
| Claim boundary | 3 | Footer boundary present |
| Venue fit | 2 | Good deck prototype, not yet main-paper final |

Decision:

- Keep as the current Figure 1 candidate.
- Use the legacy `premium_fig1_core_tissue_pathway` path only for continuity;
  the preferred asset is `premium_fig1_tissue_transfer_hierarchy`.
- The main-readout summary is now exported as
  `output/premium_tables/table1_tissue_transfer_summary.csv`.

## Figure 2: Pathway Artifact / Rescue Prototype

Asset:

- `output/premium_figures/premium_fig2_pathway_artifact_rescue.png`

Strengths:

- Strong standalone message.
- Panel A gives a clear "gene artifacts vs pathway abstraction" argument.
- Panel B is concise and visually honest about task-specific rescue.
- Panel C connects biology to transfer success.
- Static PNG/PDF with source manifest.

Weaknesses:

- Title is long for slide use.
- Panel C label placement should be rechecked at manuscript scale.
- Gastrocnemius outlier needs clearer visual encoding if central to caption.
- Confounder tasks mix mission, hardware, and gravity; caption must explain
  collinearity and exploratory scope.

Scores:

| Category | Score | Note |
|---|---:|---|
| Message clarity | 3 | Very clear central message |
| Data trust | 2 | Manifest exists; final numeric audit still needed |
| Static reproducibility | 3 | PNG/PDF/script/manifest all exist |
| Typography | 3 | Good readability |
| Layout | 3 | Balanced and cleaner than Figure 1 |
| Color semantics | 2 | Gene/pathway colors work; project palette still needs finalization |
| Accessibility | 2 | Good contrast but needs colorblind check |
| Label economy | 3 | Labels are mostly efficient |
| Claim boundary | 3 | Footer boundary present |
| Venue fit | 3 | Candidate for main deck/paper after numeric audit |

Decision:

- Keep as premium P0/P1 candidate.
- This is likely a better home for pathway rescue than Figure 1 panel C.

## Cross-Figure Design Decision

The first six regenerated figures already show a better direction than legacy
HTML assets:

- static;
- reproducible;
- source-manifested;
- visually consistent;
- claim-boundary aware.

JK style correction applied:

- Avoid card-box/dashboard visual language in scientific figures.
- Prefer high-profile journal patterns: thin-rule tables, dot/lollipop
  summaries, axis-based plots, schematic spines, and restrained whitespace.
- Figure 4, Figure 5, and Figure 6 were revised accordingly after their first
  card-like prototypes.

Next design decision:

> Should Figure 1 become a cleaner tissue-hierarchy-only figure, with all
> pathway rescue/batch-resistance content moved into Figure 2?

Recommendation:

- Yes for manuscript.
- Maybe no for a compact talk, where Figure 1 can remain a combined "core
  result" slide.

## Figure 3: Model-Tier Comparison Prototype

Asset:

- `output/premium_figures/premium_fig3_model_tier_comparison.png`

Strengths:

- Directly supports the manuscript/deck thesis that model scale does not erase
  small-n spaceflight domain shift.
- Panel A is immediately readable and shows the tuned classical reference above
  current FM/LLM rows.
- Panel B prevents over-aggregation by showing tissue-specific FM deltas.
- Panel C makes the benchmark-surface boundary visible inside the figure.
- Static PNG/PDF with source manifest.

Weaknesses:

- The comparison intentionally mixes surfaces: 6-tissue v1 means, zero-shot
  text tasks, and v3 extension best rows.
- Panel A uses a zoomed x-axis starting at 0.4, which is good for slide
  readability but should be reviewed for manuscript conservatism.
- Panel C is text-heavy and may need a shorter manuscript version.
- Final caption must say this is not a universal 8-tissue FM leaderboard.

Scores:

| Category | Score | Note |
|---|---:|---|
| Message clarity | 3 | Strong and immediate |
| Data trust | 2 | Manifest exists; final numeric audit still needed |
| Static reproducibility | 3 | PNG/PDF/script/manifest all exist |
| Typography | 3 | Readable at deck scale |
| Layout | 3 | Balanced A/B/C hierarchy |
| Color semantics | 2 | Family colors work; palette still needs project-wide lock |
| Accessibility | 2 | Good contrast, needs colorblind check |
| Label economy | 2 | Panel C is dense |
| Claim boundary | 3 | Boundary is explicit in panel and footer |
| Venue fit | 3 | Strong main-story candidate after caption hardening |

Decision:

- Keep as premium P0/P1 candidate.
- Use this as the model-family result slide, with scFoundation/UCE treated as
  extension rows rather than as direct 6-tissue leaderboard entries.

## Figure 4: v9 Platform Architecture Prototype

Asset:

- `output/premium_figures/premium_fig4_v9_platform_architecture.png`

Strengths:

- Makes v9 legible as a benchmark platform rather than a single score table.
- Shows the source-to-release architecture as a schematic.
- Keeps extension-lane readiness as a schematic, while exact rows are exported
  as tables.
- Static PNG/PDF with source manifest.

Weaknesses:

- More deck-native than manuscript-native.
- Some status phrases are manually curated from multiple local reports; final
  manuscript use should pair it with a formal status table.
- The lane schematic still carries operational status language, so it should be
  used as a resource/platform figure rather than a biological result figure.

Scores:

| Category | Score | Note |
|---|---:|---|
| Message clarity | 3 | Clear v9 status/story |
| Data trust | 2 | Manifest exists; status rows need final table freeze |
| Static reproducibility | 3 | PNG/PDF/script/manifest all exist |
| Typography | 3 | Good deck-scale readability |
| Layout | 3 | Architecture plus lane schematic is cleaner than a matrix figure |
| Color semantics | 3 | Ready/draft/blocked semantics are explicit |
| Accessibility | 2 | Needs colorblind check, but text also carries status |
| Label economy | 3 | Dense but controlled |
| Claim boundary | 3 | Boundary visible in title/subtitle/footer |
| Venue fit | 3 | Excellent deck overview; manuscript needs adaptation |

Decision:

- Keep as the main v9 overview schematic.
- Pair with `output/premium_tables/table4_v9_lane_readiness.csv` and
  `output/premium_tables/table4_v9_evidence_footprint.csv`.

## Figure 5: Public Bulk Release Boundary Schematic

Asset:

- `output/premium_figures/premium_fig5_public_bulk_release_boundary_schematic.png`

Strengths:

- Clearly separates metadata-only preview from frozen-data release.
- No longer disguises gate and language tables as figure panels.
- Strong as a release-boundary schematic or deck divider.
- Static PNG/PDF with source manifest.

Weaknesses:

- This is governance/release communication, not a scientific result figure.
- It is intentionally sparse and should not be promoted as a main result.
- Two gates are still marked `needs update`, so it should be refreshed after
  the dataset card and Data Package are updated.

Scores:

| Category | Score | Note |
|---|---:|---|
| Message clarity | 3 | Boundary is immediate |
| Data trust | 2 | Manifest exists; display text is curated from claim tables |
| Static reproducibility | 3 | PNG/PDF/script/manifest all exist |
| Typography | 3 | Readable after overflow fix |
| Layout | 2 | Clean schematic, but intentionally sparse |
| Color semantics | 3 | Pass/update/blocked is clear |
| Accessibility | 2 | Text carries status; colorblind check still needed |
| Label economy | 3 | Dense but controlled after shortening |
| Claim boundary | 3 | This is the figure's main purpose |
| Venue fit | 3 | Strong deck/release slide; paper needs formal table support |

Decision:

- Keep as a v9 public-bulk alpha governance slide.
- Do not present it as a benchmark-performance figure.
- Pair with the three table outputs in `output/premium_tables/`.

## Figure 6: Human Organoid Diagnostic Surface Prototype

Asset:

- `output/premium_figures/premium_fig6_organoid_diagnostic_surface.png`

Strengths:

- Gives the organoid extension a clear visual role without overpromoting it.
- Panel A quickly establishes the footprint: 2 sources, 42 samples, 8 matched
  DE contrasts, 242.7k DE-reference rows, and draft status.
- Panel B separates primary prediction metrics from flight-response pattern and
  model-gene-effect biology checks.
- Panel C shows enrichment among top-ranked model genes while preserving
  heterogeneity and non-primary interpretation.
- Default, secondary, and negative biology-check decisions are now exported as
  a separate table.
- Static PNG/PDF with source manifest.

Weaknesses:

- Panel B now separates metric families, but caption language must still
  prevent leaderboard interpretation.
- Organoid surface remains draft, small-n, and coupled by source/disease/donor
  factors.

Scores:

| Category | Score | Note |
|---|---:|---|
| Message clarity | 3 | Clear diagnostic-extension story |
| Data trust | 2 | Manifest exists; final table freeze still needed |
| Static reproducibility | 3 | PNG/PDF/script/manifest all exist |
| Typography | 3 | Good after overlap fixes |
| Layout | 3 | Three-panel figure is cleaner after moving the decision table out |
| Color semantics | 2 | Draft/diagnostic/blocked colors work, but palette needs global lock |
| Accessibility | 2 | Text carries status; colorblind check still needed |
| Label economy | 3 | Dense but controlled |
| Claim boundary | 3 | Boundary is explicit throughout |
| Venue fit | 3 | Strong deck slide; manuscript needs simplification |

Decision:

- Keep as the main human organoid extension slide.
- Use it to show breadth and diagnostic maturity, not leaderboard performance.
- Pair with `output/premium_tables/table6_organoid_biology_check_decisions.csv`.

## Next Figure Family

Continue v9-native visuals:

- multispecies extension-task figure;
- single-cell blocker/status diagram if the RRRM lane stays in the talk.
