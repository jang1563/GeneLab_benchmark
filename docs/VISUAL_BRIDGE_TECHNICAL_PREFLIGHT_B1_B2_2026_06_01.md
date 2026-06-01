# Visual Bridge Technical Preflight: B1-B2

Date: 2026-06-01

Purpose: translate the B1/B2 content briefs into production constraints before
rendering. These two slides introduce the evaluation-layer logic, so they must
feel explanatory rather than decorative.

Companion files:

- `docs/VISUAL_BRIDGE_CONTENT_BRIEFS_B1_B4_2026_06_01.md`
- `docs/VISUAL_METHODS_BRIDGE_AND_CONSULTING_BRIEF_2026_06_01.md`
- `docs/VISUAL_TECHNICAL_PRODUCTION_PROTOCOL_2026_06_01.md`
- `scripts/build_visual_bridge_scene_contracts.py`
- `scripts/audit_visual_scene_contracts.py`

## Why B1 And B2 Matter

B1/B2 are the deck's first methods bridge.

They should answer two questions before any result figure appears:

- why public space-omics studies need an evaluation layer;
- what concrete object becomes a benchmark task.

If these slides are vague, later mission-held-out and train-only processing
slides can look like arbitrary modeling choices. If they are too technical, a
first-time viewer can get stuck on accessions, manifests, and package structure
before understanding the benchmark unit.

## B1 Contract

Slide ID:

- `b1_evaluation_layer`

Output family:

- `output/premium_bridge_scenes/b1_evaluation_layer/`

Decision headline:

> Public space omics is valuable but fragmented; SpaceBio-Bench adds an
> auditable evaluation layer.

Primary visual move:

Scattered public-study marks compress into one organized evaluation surface.

Required visible elements:

- `Public studies`
- `Missions, tissues, samples, labels`
- `Benchmark tasks`
- `Audited scores`
- `public-data resource; source-traceable`

Forbidden visible terms:

- internal release status;
- package names;
- raw accession lists;
- model family names;
- raw local paths.

Technical risk:

- The slide can become a generic overview diagram. Preserve one movement from
  fragmentation to evaluation, and keep text sparse.

## B2 Contract

Slide ID:

- `b2_study_to_task`

Output family:

- `output/premium_bridge_scenes/b2_study_to_task/`

Decision headline:

> A benchmark task is a source-traceable sample/label contract.

Audience question:

> What is the unit of the benchmark, and how does a public study become a task?

Primary visual move:

Nested source context opens into a thin task-record strip:

`public source -> mission context -> samples -> labels/tissue -> task record`

Required visible elements:

- `Public source`
- `Mission context`
- `Samples`
- `Labels`
- `Tissue / assay`
- `Task record`
- `source-traceable task record`

Required caveat:

- `Exact rows stay in source tables.`

Evidence sources:

- `docs/VISUAL_BRIDGE_CONTENT_BRIEFS_B1_B4_2026_06_01.md`
- `docs/PROJECT_RESULTS_LOCATION_INVENTORY_2026_05_31.md`
- `docs/V1_V9_SLIDE_ASSET_MANIFEST_2026_05_31.md`
- `v9/source_inventory.csv`
- `v9/task_manifest_index.csv`

Forbidden visible terms:

- raw accession lists;
- `task manifest` as the only explanatory phrase;
- package or code structure;
- `LOMO`;
- `payload`;
- `artifact`;
- `RRRM`;
- `alpha`;
- raw local paths.

Technical risk:

- A table-shaped slide would make the task look like a spreadsheet row rather
  than an auditable scientific contract. Use a document strip and a single
  transformation path instead of a full table.

## Shared Pre-Render Checks

The generated contract must pass:

1. `manifest.json` exists.
2. `overlay_spec.json` exists.
3. `qa.json` exists.
4. slide IDs match across all contract files.
5. canvas is `3840 x 2160`, 16:9.
6. content brief and technical preflight paths exist.
7. all evidence source paths exist.
8. visible text is under the 45-word bridge-slide budget.
9. forbidden visible terms are absent.
10. coordinates stay inside the safe area.
11. output image paths are declared before rendering.

## What This Does Not Approve

Passing preflight only means the slide is ready for a prototype render.

The post-render gate still requires:

- full-size visual inspection;
- thumbnail/contact-sheet inspection;
- text overlap check;
- caveat visibility check;
- source/claim-boundary review;
- family review against B3/B4.
