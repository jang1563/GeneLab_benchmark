# Visual Expert Panel Review

Date: 2026-05-31

Scope: strict review of newly generated premium figures under a high-profile
journal standard. This review intentionally treats the current figures as if
they were being screened by a scientific graphics editor, statistical reviewer,
domain reviewer, and accessibility/layout reviewer.

Reviewed assets:

- `output/premium_figures/premium_fig1_tissue_transfer_hierarchy.png`
- `output/premium_figures/premium_fig2_pathway_artifact_rescue.png`
- `output/premium_figures/premium_fig3_model_tier_comparison.png`
- `output/premium_figures/premium_fig4_v9_platform_architecture.png`
- `output/premium_figures/premium_fig5_public_bulk_release_boundary_schematic.png`
- `output/premium_figures/premium_fig6_organoid_diagnostic_surface.png`

P0 revision note:

- The detailed review below records the strict pre-P0 critique and should be
  read as the design rationale for the current revisions.
- P0 changes have now been applied in `scripts/build_premium_visuals.py`:
  Figure 1 is tissue-transfer-only; Figure 2 Panel C no longer uses a fitted
  trend or outlier-exclusion claim; Figure 6 Panel B now separates primary
  prediction, flight-response pattern, and model-gene-effect checks into
  separate mini-axes.
- Figure-visible jargon was audited in
  `docs/VISUAL_TEXT_JARGON_AUDIT_2026_05_31.md`.
- Table-like content was split out in
  `docs/VISUAL_FIGURE_TABLE_SEPARATION_2026_05_31.md`.

## Executive Verdict

None of the six figures should be treated as final high-profile journal main
figures yet.

Current best candidates:

| Rank | Figure | Verdict | Use now |
|---:|---|---|---|
| 1 | Figure 2 pathway artifact/rescue | strongest scientific figure | main deck / manuscript candidate after caption hardening |
| 2 | Figure 6 organoid diagnostic surface | strong extension figure | deck figure; decision table now separate |
| 3 | Figure 3 model-tier comparison | strong story, high overclaim risk | deck figure; manuscript after axis and coverage redesign |
| 4 | Figure 1 tissue-transfer hierarchy | cleaner after P0 split | main deck candidate; manuscript after numeric audit |
| 5 | Figure 4 v9 platform architecture | useful platform overview | deck/resource-paper figure; readiness table separate |
| 6 | Figure 5 public bulk boundary schematic | governance/status schematic | release deck or supplement only |

Strict global issue:

- Figure 4-6 are improved after removal of card-box UI, but Figure 4 and Figure
  5 still feel more like release/status communication than scientific evidence.
  That is acceptable for a deck, but not for a main research figure.

## Expert Panel Criteria

The review applies four expert lenses:

1. Scientific graphics editor: does it look like a polished paper figure rather
   than dashboard output?
2. Statistical visualization reviewer: are axes, baselines, comparisons, and
   exclusions visually honest?
3. Biological domain reviewer: does the visual match the scientific claim and
   avoid overinterpretation?
4. Accessibility/layout reviewer: does it survive grayscale, projection, and
   manuscript reduction?

Ratings:

- `Accept`: usable with minor caption edits.
- `Major revision`: scientifically promising but figure-level redesign needed.
- `Reject for main figure`: useful internally or in supplement/deck, but not
  main-paper material in current form.

## Figure 1: Core Tissue / Pathway

P0 update:

- Rebuilt as `premium_fig1_tissue_transfer_hierarchy`.
- Current version focuses on held-out mission tissue transfer only.
- Detection and pathway-rescue panels were removed from Figure 1 to avoid a
  three-claim composite.
- The pre-P0 critique below is retained as rationale for the split.

Panel summary:

- A: cross-mission tissue transfer ranking.
- B: Category A detection.
- C: pathway abstraction rescue.

Verdict: major revision.

Main strengths:

- Panel A is a strong lead panel.
- Tissue hierarchy is clear.
- Confidence intervals and chance line are helpful.
- Claim boundary is visible.

Major concerns:

- The figure tries to carry three different claims: tissue transfer, detection
  significance, and pathway rescue. This weakens the main message.
- Panel C duplicates Figure 2's conceptual space and should not compete with
  the core tissue hierarchy.
- The title promises a broad conclusion but the figure reads as a composite
  slide rather than a refined main paper figure.
- The `FDR<0.05` marker in Panel B is visually ambiguous unless the caption
  explicitly says it refers to skin only.
- Color semantics are not fully locked: gold means thymus in Panel A and Panel
  B, red means skin significance, blue means default tissue. This is workable
  but not yet journal-clean.

Required fix:

- Split into a tissue-transfer main figure and move pathway rescue fully to
  Figure 2.
- Keep Figure 1 focused on mission-held-out benchmark behavior.
- Consider a manuscript version with only Panel A plus a small fold/source
  schematic or compact significance inset.

Use decision:

- Deck: acceptable after title shortening.
- Manuscript main: not yet.

## Figure 2: Pathway Artifact / Rescue

P0 update:

- Panel C no longer contains a fitted trend line or a rank-correlation callout
  based on excluding gastrocnemius.
- Visible `NES concordance` language was replaced with `pathway activity
  agreement`; NES remains in source manifests/captions.

Panel summary:

- A: gene versus pathway confounder prediction.
- B: pathway-minus-gene rescue delta.
- C: NES concordance versus transfer success.

Verdict: strongest current candidate; P0 statistical hardening applied, still
needs careful caption boundaries.

Main strengths:

- The scientific claim is coherent: pathway abstraction can reduce artifacts
  and rescue selected weak gene-level tasks.
- Panel A is memorable and high-value.
- Panel B is honest about task specificity.
- Panel C creates a biologically interpretable bridge.

Major concerns:

- Panel A includes confounded tasks with gene bars near 1.0. This is visually
  dramatic but can be interpreted as "pathways solve batch effects" unless the
  caption states that mission, hardware, and gravity labels are coupled.
- Pre-P0 concern now resolved: Panel C used a fitted trend with very small n
  and an explicit outlier exclusion. The current version removes the fitted
  line and uses an exploratory label.
- The gastrocnemius outlier needs clearer treatment. Right now it is gray, but
  the exclusion logic is embedded in text rather than visually formalized.
- Title is strong for a talk but slightly too declarative for a paper unless
  the scope is carefully bounded.

Required fix:

- Done: convert Panel C to an exploratory association panel with no fitted
  line.
- Add a small note that confounder classifiers are diagnostic artifacts, not
  biological endpoint metrics.
- Keep this as the main pathway-rescue figure after removing Panel C from
  Figure 1.

Use decision:

- Deck: strong.
- Manuscript main: likely, after Panel C hardening.

## Figure 3: Model-Tier Comparison

Panel summary:

- A: mean score comparison across classical ML, single-cell FMs, and text LLMs.
- B: tissue-level FM deltas.
- C: coverage boundary.

Verdict: major revision for manuscript, strong for deck.

Main strengths:

- The central message is clear: model scale does not automatically solve
  small-n spaceflight domain shift.
- Panel B prevents a naive mean-only interpretation.
- The figure explicitly warns against a uniform 8-tissue FM leaderboard.

Major concerns:

- Panel A uses a truncated x-axis starting at 0.4. This is defensible for deck
  readability but will draw reviewer/editor scrutiny in a paper.
- Panel A mixes metric surfaces: 6-tissue v1 means, zero-shot text tasks, and
  v3 extension best rows. Panel C flags this, but the visual still invites
  direct ranking.
- Panel C remains too text-heavy and looks like an annotation block rather than
  a scientific panel.
- The purple/magenta family colors are visually loud and less journal-like than
  restrained point/interval encodings.

Required fix:

- Make a manuscript variant using dot plots or slope/delta plots rather than
  horizontal bars.
- Either start the Panel A axis at 0 or explicitly mark it as a zoomed score
  range with subdued bars.
- Convert Panel C into a small thin-rule coverage table.

Use decision:

- Deck: good.
- Manuscript main: not yet.

## Figure 4: v9 Platform Status

Panel summary:

- A: v9 lane readiness table.
- B: evidence footprint list.
- C: platform schematic spine.

Verdict: useful overview, reject as main scientific result figure.

Main strengths:

- Much improved after removal of card/tile UI.
- Thin-rule table style is closer to journal language.
- The lane distinction prevents scope drift.
- Panel C is a clean platform schematic.

Major concerns:

- This is a platform/release overview, not a scientific result.
- Panel B mixes unlike quantities: sources, folds, rows, samples, blockers.
  It works as evidence footprint, not as quantitative comparison.
- The figure may be too operational for a manuscript main figure unless the
  paper is framed as a resource/platform paper.
- Some wording still feels internal: `payload blockers`, `metadata alpha`.

Required fix:

- For paper, use as a resource overview/supplement or convert to a cleaner
  platform schematic plus small status table.
- Replace internal project labels with reader-facing terms.

Use decision:

- Deck: good.
- Manuscript main: only if paper is explicitly a benchmark resource paper.

## Figure 5: Public Bulk Alpha Boundary

Panel summary:

- A: metadata-alpha gates.
- B: selected release path.
- C: allowed/prohibited language.

Verdict: reject for main scientific figure; useful release/governance figure.

Main strengths:

- Clear and honest claim boundary.
- Figure now uses rule-table/spine style rather than cards.
- Good for preventing overclaiming in public materials.

Major concerns:

- This is not a scientific evidence figure.
- The title is operational and release-process oriented.
- It should not appear in a high-profile journal main figure set unless the
  manuscript includes a data-resource governance section.
- Panel C is valuable internally but reads like editorial policy, not research
  result.

Required fix:

- Keep for deck, appendix, release package, or supplement.
- For paper, convert to prose/table in Methods or Supplementary Data.

Use decision:

- Deck/release: strong.
- Manuscript main: no.

## Figure 6: Human Organoid Diagnostic Surface

P0 update:

- Panel B now uses separate mini-axes for primary task, flight-response
  pattern, and model gene effects.
- The figure no longer visually ranks AUROC, macro-F1, direction match, and
  rank correlation on a single score axis.
- Figure text now uses `biology check`, `secondary check`, and
  `model gene effects` instead of internal diagnostic labels where possible.

Panel summary:

- A: provenance/evidence spine.
- B: primary classification versus diagnostic-only biological metrics.
- C: enrichment among top-ranked model genes.

Verdict: promising extension figure; P0 metric-family separation applied,
decision table moved out for manuscript/deck cleanliness.

Main strengths:

- Panel A is much improved after removing card tiles.
- The figure gives organoids a credible role without claiming a leaderboard.
- Panel C is scientifically useful: top-k enrichment is shown with expected
  values and significance pattern.
- Default, secondary, and negative biology-check decisions are preserved as a
  separate table rather than as a figure panel.

Major concerns:

- Pre-P0 concern now resolved: Panel B placed AUROC, macro-F1, direction
  match, and rank correlation on one axis, inviting inappropriate comparison.
  The current version separates these into metric-family mini-axes.
- The flight-response pattern values still need caption language clarifying
  that they are biology checks, not the same endpoint as prediction metrics.
- The small-n/coupled-source limitation is in title/subtitle but could be more
  visually integrated.

Required fix:

- Done: split Panel B into grouped mini-axes for primary task,
  flight-response pattern, and model gene effects.
- Done: replace `diagnostic-only` text with biology-check and secondary-check
  language.
- Done: move the decision table out of the figure and into
  `output/premium_tables/table6_organoid_biology_check_decisions.csv`.

Use decision:

- Deck: strong.
- Manuscript main: possible after Panel B redesign.

## Cross-Figure System Review

What is working:

- Static PNG/PDF generation is solid.
- Source manifests are a major strength.
- Claim boundaries are consistently present.
- The visual tone improved substantially after removing dashboard/card boxes.
- Figure 2, 3, and 6 contain publishable scientific stories.

What still blocks L4:

- Too many figures carry deck-style explanatory titles.
- Several panels mix metric families or evidence types.
- Color semantics are not yet fully project-wide:
  - gold = thymus / draft / warning in different contexts;
  - blue = classical ML / gene / single-cell evidence in different contexts;
  - green = pathway / ready / enrichment.
- Some status figures are operational rather than scientific.
- Final numeric audit is still pending.

## Priority Revision Plan

P0 before any serious deck/paper packaging:

1. Done: split Figure 1 into a cleaner tissue-transfer figure.
2. Done: hardened Figure 2 Panel C by removing the fitted trend line.
3. Done: redesigned Figure 6 Panel B so different metric families are not
   compared as one score axis.
4. In progress: global color and text contract created; needs colorblind and
   manuscript-size verification before L4.

P1 for manuscript-grade polish:

1. Create manuscript variants for Figures 2, 3, and 6.
2. Move Figure 5 to supplement/release appendix.
3. Make Figure 4 a platform schematic plus concise status table if retained.
4. Add exact source-data tables for each final figure.

P2 for final premium system:

1. Run grayscale/colorblind checks.
2. Render manuscript-size panels, not only 16:9 deck figures.
3. Add figure captions with scope, fold definitions, and claim boundaries.
4. Freeze final source-of-truth numbers before assigning L4.

## Recommended Main Story Set After Revision

Main figures:

1. Benchmark/task design and mission-held-out tissue hierarchy.
2. Pathway artifact suppression and task-specific rescue.
3. Model-family comparison under small-n domain shift.
4. v9 benchmark platform architecture and extension lanes.
5. Human organoid diagnostic extension.
6. Multispecies extension readiness/performance once built.

Supplement/release figures:

- Public bulk metadata-alpha boundary.
- Detailed v9 lane readiness table.
- Single-cell payload blocker status.
- Negative/secondary organoid diagnostics.
