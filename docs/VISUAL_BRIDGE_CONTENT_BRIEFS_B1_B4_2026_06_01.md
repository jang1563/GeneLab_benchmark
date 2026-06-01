# Visual Bridge Content Briefs B1-B4

Date: 2026-06-01

Purpose: define the first four bridge slides before any additional slide
rendering. These briefs are the gate between content strategy and visual
production.

Companion:

- `docs/VISUAL_METHODS_BRIDGE_AND_CONSULTING_BRIEF_2026_06_01.md`
- `docs/VISUAL_METHODS_EXPLANATION_GAP_MAP_2026_06_01.md`
- `docs/VISUAL_METHODS_STORYBOARD_2026_06_01.md`
- `docs/V1_V9_PRESENTATION_AND_MANUSCRIPT_MASTER_OUTLINE_2026_05_31.md`

## Shared Design Rule For B1-B4

These are explanation slides, not result slides.

Use the light methods/provenance grammar:

- warm off-white canvas;
- thin rules;
- sparse evidence strips;
- no dashboard cards;
- no decorative space wallpaper;
- one dominant visual movement;
- one compact status/caveat label.

Do not use the dark result-slide grammar until the deck reaches the main result
block.

## B1. Public Space Omics Needs An Evaluation Layer

### Audience Question

Why is this project more than another reanalysis of public GeneLab/OSDR data?

### Decision Headline

Public space omics is valuable but fragmented; SpaceBio-Bench adds an auditable
evaluation layer.

### Bridge Role

- Previous: title/thesis slide.
- Next: study/mission/sample/task hierarchy.

B1 tells the viewer what problem the project solves before asking them to care
about mission-held-out splitting.

### What To Put On The Slide

Visible content:

- `Public studies` as the source field.
- `Missions, tissues, samples, labels` as the data-structure field.
- `Benchmark tasks` as the organized evaluation field.
- `Audited scores` as the output field.
- One small status label: `public-data resource; source-traceable`.

Evidence/support:

- `docs/PROJECT_SLIDE_CONTENT_INVENTORY_V1_TO_V9_2026_05_31.md`
- `docs/PROJECT_RESULTS_LOCATION_INVENTORY_2026_05_31.md`
- `docs/SIMILAR_PROJECTS_AND_POSITIONING_SCAN_2026_05_31.md`
- v9 task/source inventory docs and manifests where available.

### What Not To Put

- Full source inventory table.
- All version history.
- Every species/modality lane.
- Model family names.
- `payload`, `RRRM`, `alpha`, or other internal status terms.

### Visual Move

Scattered source dots compress into a single horizontal evaluation surface.

Layout:

- left third: loosely arranged public-study source marks;
- center: alignment rail labeled with mission/tissue/sample/label;
- right third: one clean benchmark surface with score/audit symbols;
- bottom right: source-traceability status label.

### Premium Consulting Feel

Make it feel like a strategic problem-definition slide:

- headline is a conclusion, not a label;
- left side is intentionally fragmented, right side intentionally ordered;
- no more than four visual nouns;
- use thin-rule annotation rather than boxes.

### Speaker Note

NASA OSDR/GeneLab already contains rich public studies. The contribution here is
to make them evaluable under a repeatable mission-held-out benchmark contract:
source records, task definitions, models, metrics, and audit trails.

### Acceptance Test

A new viewer should be able to say:

> The project organizes public space omics into a benchmark, rather than just
> plotting one dataset.

## B2. Studies Become Tasks Through Missions, Samples, And Labels

### Audience Question

What is the unit of the benchmark, and how does a public study become a task?

### Decision Headline

A benchmark task is a source-traceable sample/label contract, not just an
expression matrix.

### Bridge Role

- Previous: B1 explains why an evaluation layer is needed.
- Next: B3 explains why a mission is hidden during evaluation.

B2 prevents the common confusion between study, mission, sample, tissue, and
task.

### What To Put On The Slide

Visible content:

- `Public source`.
- `Mission context`.
- `Samples`.
- `Labels`.
- `Tissue / assay`.
- `Task record`.

One simplified example row:

- source/study mark;
- mission mark;
- sample count placeholder or sample icons;
- flight/ground label chips;
- tissue label;
- task record output.

Evidence/support:

- v9 public bulk source inventory and task manifest materials.
- `docs/PROJECT_RESULTS_LOCATION_INVENTORY_2026_05_31.md`
- `docs/V1_V9_SLIDE_ASSET_MANIFEST_2026_05_31.md`

### What Not To Put

- Do not show all task rows.
- Do not show raw accession soup.
- Do not make `task manifest` the only explanatory phrase.
- Do not include code/package structure.

### Visual Move

Nested containment opens into a task record:

`source -> mission -> samples -> labels/tissue -> task record`

The task record should look like a thin scientific document strip, not a UI
card.

### Premium Consulting Feel

Make the conceptual hierarchy obvious:

- one nested stack on the left;
- one task-record strip on the right;
- a single arrow between them;
- definitions are tiny margin notes, not body paragraphs.

### Speaker Note

The task is not simply "classify rows." Each task records where the samples came
from, what labels are being predicted, what tissue and assay context define the
comparison, and how the fold and metric should be evaluated.

### Acceptance Test

A new viewer should be able to say:

> A task is an auditable contract connecting source data, mission context,
> samples, labels, and evaluation.

## B3. The Test Set Is A Mission, Not A Random Sample

### Audience Question

Why is mission-held-out evaluation central to this benchmark?

### Decision Headline

The benchmark hides an entire mission so models must generalize across mission
conditions.

### Bridge Role

- Previous: B2 defines what a task is.
- Next: B4 explains how train-only processing prevents leakage.

B3 is the highest-priority methods bridge because it changes how the viewer
interprets every result figure that follows.

### What To Put On The Slide

Visible content:

- three or four training mission marks;
- one hidden test mission mark;
- boundary line between train and test;
- repeated fold hint, but only as a small rotational/timeline cue;
- label: `hidden mission scored after training`.

Evidence/support:

- `docs/VISUAL_METHODS_EXPLANATION_GAP_MAP_2026_06_01.md`
- `docs/VISUAL_METHODS_STORYBOARD_2026_06_01.md`
- v1/v7 result summaries for mission-held-out framing.

### What Not To Put

- `LOMO` as the headline.
- Random CV notation.
- Fold tables.
- Dataset-specific mission IDs unless the slide intentionally uses one
  concrete example.
- Model names.

### Visual Move

One mission moves behind a clean boundary while the remaining missions feed the
training lane.

Layout:

- left: training mission field;
- center: model training/evaluation rail;
- right: hidden mission behind a vertical trust boundary;
- bottom: one sentence caveat: `No hidden-mission samples are used to choose
  features or fit models`.

### Premium Consulting Feel

This slide should feel like a diligence control diagram:

- calm, spacious layout;
- no dramatic warning graphics;
- red used only for the boundary/test mission;
- all verbs plain: `train`, `hide`, `score`.

### Speaker Note

The independence unit is the mission. The model cannot see any sample from the
hidden mission during training, which makes the score a cross-mission
generalization test rather than a random-sample performance estimate.

### Acceptance Test

A new viewer should be able to say:

> The test set is a whole mission, so the benchmark is asking about
> cross-mission generalization.

## B4. Feature Processing Must Stay On The Training Side

### Audience Question

How does the benchmark avoid learning from the mission it is supposed to test?

### Decision Headline

Feature selection, scaling, and model fitting are learned from training
missions only.

### Bridge Role

- Previous: B3 shows that one mission is hidden.
- Next: gene/pathway feature explanation and then the Fig1/Fig2 result block.

B4 makes the methodology trustworthy before any performance number is shown.

### What To Put On The Slide

Visible content:

- train mission lane;
- hidden mission lane;
- operations on train lane only:
  - `choose features`;
  - `scale`;
  - `fit model`;
- hidden mission bypass arrow to:
  - `score`;
- one status label: `train-only processing`.

Evidence/support:

- pipeline/task manifest documentation where available;
- `docs/VISUAL_METHODS_EXPLANATION_GAP_MAP_2026_06_01.md`;
- current methods overview output:
  `output/premium_methods_scenes/methods_data_to_evaluation_overview.png`.

### What Not To Put

- function names;
- package names;
- scaler/model implementation details;
- every preprocessing option;
- long caveat paragraph.

### Visual Move

Two lanes move in parallel until the hidden mission joins only at scoring.

Layout:

- upper lane: training missions -> choose features -> scale -> fit model;
- lower lane: hidden mission -> wait outside fitting operations -> score;
- vertical guard before scoring;
- status label beside guard.

### Premium Consulting Feel

This should feel like a process-control slide:

- disciplined lines and gates;
- no cartoon lock icon;
- minimal labels;
- operations written as verbs;
- caveat is precise but visually quiet.

### Speaker Note

If feature selection or scaling used the hidden mission, the benchmark would
leak information from the test set. The fold must learn those operations on the
training missions and apply them to the hidden mission only at scoring.

### Acceptance Test

A new viewer should be able to say:

> The hidden mission is not used to choose features, scale data, or fit the
> model.

## B1-B4 Production Gate

Before rendering any of B1-B4:

1. Confirm the exact evidence object or manifest reference.
2. Write the headline as a conclusion.
3. Limit visible text to fewer than 45 words unless the slide is a deliberate
   appendix/table slide.
4. Identify one visual move.
5. Identify one caveat/status label.
6. Check that no internal language appears in visible slide text.
7. Sketch thumbnail structure before producing a full-resolution scene.

## Recommended Production Order

1. B3: strongest methodological misunderstanding risk.
2. B4: necessary trust bridge before results.
3. B2: clarifies task construction.
4. B1: anchors opening narrative.

This order is better for design development because B3/B4 establish the visual
grammar for the methods trust layer. B1/B2 can then be built to match that
grammar without becoming generic overview slides.
