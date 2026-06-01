# Visual Methods Storyboard

Date: 2026-06-01

Purpose: prepare the methods section of the presentation/manuscript so an
unfamiliar viewer can follow data collection, processing, benchmark creation,
analysis, and release boundaries visually.

## Methods Section Goal

The viewer should be able to say:

> The project starts from public GeneLab/OSDR studies, turns them into
> mission-held-out tasks, builds gene/pathway feature layers inside each fold,
> compares models on held-out missions, audits results, and separates
> metadata-ready resources from pending data-file release gates.

## Proposed Slide Sequence

| Slide | Working Title | One-Sentence Claim | Visual Form | Notes |
|---:|---|---|---|---|
| B1 | Why a benchmark layer is needed | Public space omics needs an auditable evaluation layer | scattered sources compressed into benchmark surface | Content brief required before build |
| B2 | Studies become tasks | A task is a source-traceable sample/label contract | source -> mission -> samples -> labels -> task record | Content brief required before build |
| B3 | The test set is a mission | The benchmark hides an entire mission, not random samples | mission field with one mission behind boundary | Build before result slides |
| B4 | Train-only guard | Feature processing is learned from training missions only | two-lane process-control schematic | Build before feature/result slides |
| M1 | From public studies to held-out tasks | Public studies become auditable benchmark tasks | process rail with source, sample table, hidden mission, features, models, audit | Build first |
| M2 | One mission stays hidden | Mission is the independence unit | train/test split diagram with one mission behind a wall | Needed before Fig1 |
| M3 | Features are built inside each fold | Train-only processing prevents leakage | two-lane leakage guard schematic | Critical for methods trust |
| M4 | Genes and pathways are different views | Pathway summaries can reduce unwanted labels and expose biology | gene cloud to pathway bands | Bridge into Fig2 |
| M5 | Model comparisons have scope | Shared rows are direct comparisons; other rows are contextual | family lanes into shared task table | Bridge into Fig3 |
| M6 | Extension datasets are diagnostic | Organoid/multispecies/single-cell are different readiness lanes | source footprint plus caution layer | Bridge into Fig6/v9 |
| M7 | Release boundary is explicit | Metadata can be ready while data-file mirroring remains pending | provenance-document release gate | Use v9 light grammar |

## First-Time-Viewer Copy Rules

Use visible words like:

- `public studies`;
- `clean sample table`;
- `hide one mission`;
- `gene/pathway features`;
- `train models`;
- `score held-out mission`;
- `audit and release boundary`.

Avoid visible words like:

- `LOMO` without expansion;
- `payload`;
- `NES`;
- `artifact`;
- `macro-F1` unless the metric is the visual topic;
- `alpha` unless paired with `metadata-only`.

## Speaker-Note Anchors

For M1:

- OSDR/GeneLab provides public sample tables and RNA-seq matrices.
- The benchmark turns studies into task manifests with labels, folds, features,
  metrics, and provenance.
- The critical idea is hiding a whole mission, not just random samples.

For M2/M3:

- Each fold trains on all but one mission and scores the hidden mission.
- Any feature selection or scaling must be learned from training missions only.
- This prevents the test mission from influencing feature construction.

For M4:

- Gene-level features can carry strong mission or hardware signatures.
- Pathway summaries can dampen some unwanted labels, but do not universally
  improve every task.

For M5:

- Direct model comparisons require the same task rows.
- Some foundation-model rows are contextual and must be captioned that way.

For M6/M7:

- Extension datasets are valuable but not equivalent to the frozen mouse bulk
  benchmark surface.
- Public bulk v9 is currently metadata-ready, not payload-frozen.

## Current Build Status

Created as visual prototype:

- `output/premium_methods_scenes/methods_data_to_evaluation_overview.png`
- `output/premium_methods_scenes/methods_data_to_evaluation_overview.pdf`
- `output/premium_methods_scenes/methods_data_to_evaluation_overview_manifest.json`

This is the first methods explainer candidate. It should be followed by a
dedicated `one mission stays hidden` slide before final deck production.

## Production Gate

Do not render additional polished methods slides until the slide has a content
brief covering:

- audience question;
- decision headline;
- evidence object;
- visual move;
- bridge role;
- visible caveat;
- terms allowed on the slide;
- terms pushed to notes;
- premium constraint.

Bridge/content planning files:

- `docs/VISUAL_METHODS_BRIDGE_AND_CONSULTING_BRIEF_2026_06_01.md`
- `docs/VISUAL_BRIDGE_CONTENT_BRIEFS_B1_B4_2026_06_01.md`
- `docs/VISUAL_TECHNICAL_PRODUCTION_PROTOCOL_2026_06_01.md`

Technical gate:

- choose a stable `slide_id`;
- declare evidence sources before rendering;
- write `scene_plate.png`, `rendered_preview.png`, `rendered_preview.pdf`,
  `overlay_spec.json`, `manifest.json`, and `qa.json` or equivalent files;
- inspect the rendered PNG at full size and thumbnail size;
- record forbidden-term and visible-overlap checks.

Recommended next production order:

1. B3 `The test set is a mission, not a random sample`.
2. B4 `Feature processing must stay on the training side`.
3. B2 `Studies become tasks through missions, samples, and labels`.
4. B1 `Public space omics needs an evaluation layer`.
