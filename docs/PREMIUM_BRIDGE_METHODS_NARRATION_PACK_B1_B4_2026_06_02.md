# Premium Bridge Methods Narration Pack: B1-B4

Date: 2026-06-02

Purpose:

- prepare B1-B4 for deck assembly without rushing into slide production;
- make the data collection, task construction, evaluation split, and train-only
  processing understandable to first-time viewers;
- keep visible figures clean while moving definitions into captions and
  speaker notes.

Companion assets:

- B1: `output/premium_bridge_scenes/b1_evaluation_layer/rendered_preview.png`
- B2: `output/premium_bridge_rebuild_scenes/b2_study_to_task_premium/rendered_preview.png`
- B3: `output/premium_bridge_rebuild_scenes/b3_mission_held_out_premium/rendered_preview.png`
- B4: `output/premium_bridge_rebuild_scenes/b4_train_only_guard_premium/rendered_preview.png`
- Family QA: `output/premium_bridge_family_review/b1_b4_premium_family_contact_sheet.png`

## Deck-Level Bridge

Recommended sequence:

1. B1: public studies need an evaluation layer.
2. B2: source records become benchmark tasks.
3. B3: one whole mission is held out.
4. B4: feature processing and model fitting stay on the training side.

Narrative promise:

> We are not just collecting public space omics studies. We are turning them
> into auditable benchmark tasks where source context is preserved, one mission
> is hidden for testing, and all processing choices are learned before the
> hidden mission is scored.

Use this block before any result figure. It earns trust before performance
numbers appear.

## B1. Evaluation Layer

Slide caption:

> Public space-omics studies are fragmented across source records, mission
> metadata, sample annotations, labels, tissues, and assays. SpaceBio-Bench
> organizes those public studies into benchmark tasks and audited score records
> so each result can be followed back to source evidence.

20-second speaker note:

> The public resource is already valuable, but it is not automatically a
> benchmark. The first step is to align studies, missions, tissues, samples, and
> labels into a repeatable evaluation layer. That lets us report scores with an
> audit trail rather than as isolated reanalyses.

Bridge to B2:

> The next question is what exactly counts as a benchmark task.

Definition to hold in notes:

- Evaluation layer: the source records, task definitions, splits, model runs,
  scores, and audit trail that make a public study comparable across methods.

Avoid saying:

- "We solved all public space-omics heterogeneity."
- "The public data are fully standardized."
- "This is a frozen payload release."

## B2. Benchmark Task

Slide caption:

> A benchmark task preserves source context. The task record keeps the public
> source, mission context, sample set, label contrast, tissue/assay context, and
> evaluation scope together instead of treating the expression matrix as a
> context-free table.

20-second speaker note:

> The unit of evaluation is not just rows in a matrix. A task records where the
> samples came from, what label is being predicted, what tissue and assay define
> the comparison, and how that task should be evaluated. This is why source rows
> remain auditable outside the figure.

Bridge to B3:

> Once a task has mission context, the split can be defined at the mission
> level rather than by random rows.

Definition to hold in notes:

- Benchmark task: a source-linked task record that binds samples, labels,
  biological context, split rule, metric, and source evidence.

Avoid saying:

- "Task contract" as the main visible metaphor.
- "Loose matrix" unless explaining why that older draft language was removed.
- "Task manifest" without first explaining that it is a machine-readable task
  definition.

## B3. Mission-Held-Out Evaluation

Slide caption:

> The test set is a whole mission. Training missions are used to fit the model,
> while the hidden mission stays outside training until the final score is
> computed. This makes the score a cross-mission generalization test rather
> than a random-sample performance estimate.

20-second speaker note:

> The independence unit is the mission. If the model sees samples from the
> hidden mission while learning, the score can overstate generalization. Holding
> out an entire mission asks whether the model can transfer across mission
> conditions rather than simply recognize nearby samples.

Bridge to B4:

> But hiding a mission is not enough; every preprocessing choice must also avoid
> using that mission.

Definition to hold in notes:

- Mission-held-out evaluation: train on all eligible missions except one, then
  score on the hidden mission.
- LOMO, if used later: leave-one-mission-out.

Avoid saying:

- "Random cross-validation."
- "Validation mission" if the hidden mission is being used as the test set.
- "Independent validation" without specifying that the independence unit is the
  mission.

## B4. Train-Only Guard

Slide caption:

> Feature selection, scaling, and model fitting are learned from training
> missions only. The hidden mission bypasses those fitting steps and is touched
> only when the trained pipeline produces a score.

20-second speaker note:

> Leakage can happen before the model is trained. If feature selection or
> scaling uses the hidden mission, the test set has already influenced the
> pipeline. The guardrail is simple: learn choices on training missions, then
> apply the learned transformation to the hidden mission only at scoring.

Bridge to next methods/result slide:

> With the split and guardrail defined, we can now explain what feature layers
> the models see and then interpret the score figures.

Definition to hold in notes:

- Train-only processing: feature selection, normalization or scaling, model
  fitting, and threshold choices are learned only on training missions.

Avoid saying:

- "No leakage" as an absolute claim unless every implemented pipeline has been
  audited at that level.
- Package or function names.
- Every preprocessing variant; those belong in methods tables or appendix.

## Caption Pack For Manuscript-Style Figure

If B1-B4 become one multi-panel methods figure:

> **Figure X. SpaceBio-Bench task construction and mission-held-out evaluation.**
> (A) Public space-omics studies are organized into an auditable evaluation
> layer linking source records, mission/sample annotations, benchmark tasks, and
> score records. (B) Each benchmark task preserves source context by binding the
> public source, mission context, sample set, label contrast, tissue/assay
> context, and evaluation scope. (C) Evaluation holds out an entire mission so
> the test set measures cross-mission generalization rather than random-sample
> performance. (D) Feature selection, scaling, and model fitting are learned
> from training missions only; the hidden mission is used only for scoring.

Short deck caption:

> Public studies become auditable benchmark tasks through source-linked task
> records, mission-held-out splits, and train-only processing.

## Remaining Explanation Needs

High priority before final deck export:

- a one-slide feature-layer bridge explaining genes versus pathway summaries;
- a compact model/metric caption explaining when model rows are direct
  shared-row comparisons;
- a release-boundary slide separating metadata/task readiness from local data
  payload readiness.

Do not solve these by adding more text inside B1-B4. Keep B1-B4 visually quiet
and use the next bridge slides or speaker notes to carry the definitions.
