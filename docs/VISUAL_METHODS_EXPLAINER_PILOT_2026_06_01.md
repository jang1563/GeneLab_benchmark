# Visual Methods Explainer Pilot

Date: 2026-06-01

Purpose: create and review the first methods slide that explains the benchmark
workflow to a viewer who has not followed the project.

## Output

- `output/premium_methods_scenes/methods_data_to_evaluation_overview.png`
- `output/premium_methods_scenes/methods_data_to_evaluation_overview.pdf`
- `output/premium_methods_scenes/methods_data_to_evaluation_overview_manifest.json`
- generator: `scripts/build_methods_explanation_scenes.py`

## Current Visual Claim

Public studies become benchmark tasks by hiding one mission, training only on
the other missions, building features inside that fold, scoring the hidden
mission, and auditing the resulting evidence.

## Manual Visual QA

Pass status: draft deck candidate.

Checks performed:

- Rendered the slide as a 16:9 high-resolution PNG.
- Opened the rendered PNG and inspected layout at full-slide scale.
- Removed metric and model shorthand from visible step labels.
- Removed the overlapping red leakage annotation from the main process rail.
- Kept the leakage concept visible in the lower mission-holdout inset.
- Confirmed that the figure reads left-to-right without requiring a table.
- Confirmed that the slide does not expose internal planning terms such as
  `artifact rescue`, `payload freeze`, or `alpha snapshot`.

## Before/After Copy Decision

Visible text that was too technical for a first-time viewer:

- `FM`;
- `LLM baselines`;
- `AUROC`;
- `macro-F1`;
- `provenance` as a primary step label.

Replaced with:

- `baseline and foundation models`;
- `held-out scores`;
- `audit trail`;
- `clean sample table`;
- `other missions train / hidden mission tests`.

The technical terms are still available for captions and speaker notes, but the
first visual read stays plain.

## Explanation Hotspots To Prepare

These are the places where the slide deck should add either a speaker-note
sentence, a small caption, or a later dedicated slide:

| Hotspot | Why It Needs Explanation | Prepared Handling |
|---|---|---|
| One mission is hidden | Viewers may expect random train/test splits | Build a dedicated mission-holdout slide next |
| Train-only feature processing | Leakage prevention is invisible without a schematic | Build a leakage-guard slide after the overview |
| Gene versus pathway features | Viewers may not know why both views are tested | Use a compression visual, not a pathway acronym list |
| Foundation model comparison | Some comparisons are direct only on shared task rows | Keep shared-row scope in captions |
| Scores and audit trail | Metrics alone do not explain trust | Pair score readouts with source/provenance checks |
| Release boundary | Metadata readiness and payload readiness differ | Use the v9 light provenance-document grammar |

## Design Assessment

The slide uses the light methods/provenance grammar rather than the dark result
grammar. This is appropriate because it teaches process rather than claiming a
scientific result.

Layer mapping:

- Z0 canvas: warm paper field with faint measurement rules.
- Z1 measurement layer: horizontal structure and process rail.
- Z2 evidence surface: mission nodes, feature compression, and train-only
  operation marks.
- Z3 interpretation layer: plain-language step titles.
- Z4 trust/status layer: hidden-mission and no-test-leakage cues.

## Remaining Work

Next methods visuals:

1. `One mission stays hidden`: make the mission-held-out split impossible to
   misread as random cross-validation.
2. `Train-only guard`: show feature selection, scaling, and model fitting on
   train missions only.
3. `Genes and pathways`: explain why pathway views are biological summaries,
   not a separate dataset.
4. `Model comparison scope`: distinguish direct shared-row comparison from
   contextual model rows.

Do not make these as table-heavy slides. Tables can be appendix or manuscript
supporting tables; methods figures should teach the workflow visually.
