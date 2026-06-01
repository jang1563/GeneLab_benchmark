# Visual Methods Bridge And Consulting Brief

Date: 2026-06-01

Purpose: pause slide production and define what must be explained, how the
slides should bridge between scientific claims, and how the final deck should
feel like a premium scientific consulting deliverable rather than a set of
loosely connected figures.

## Current Correction

Do not create additional slides until each candidate has a content brief.

The current visual assets are useful, but the next quality jump will not come
from more rendering. It will come from controlling the audience's mental model:

1. what problem this project solves;
2. how public space omics data are converted into benchmark tasks;
3. why the split design is stricter than random cross-validation;
4. how feature construction avoids leakage;
5. what the major result figures prove;
6. what the extension/resource slides do and do not claim.

The deck should feel calm, expensive, and self-evident. Every slide should have
one sentence of decision-grade meaning, one evidence object, and one visible
boundary where the evidence could be overread.

## Audience Mental Model

A first-time viewer needs to build this model in order:

| Step | Viewer Question | What They Must Learn | What Should Not Happen |
|---:|---|---|---|
| 1 | What is this project? | A benchmark/resource for evaluating models on public space-biology omics | They should not think it is only a mouse tissue analysis |
| 2 | What is the data unit? | Public studies contain samples, labels, tissues, missions, and expression matrices | They should not confuse study, mission, sample, and task |
| 3 | What is held out? | A whole mission is hidden during model training | They should not assume random sample cross-validation |
| 4 | Why is leakage a risk? | Feature selection/scaling can accidentally learn from the hidden mission | They should not trust scores without seeing train-only processing |
| 5 | What do models see? | Gene-level and pathway-level feature views, plus model family-specific inputs | They should not read pathway summaries as a new dataset |
| 6 | What is compared? | Direct comparisons require shared task rows and aligned evaluation scope | They should not read contextual rows as a strict leaderboard |
| 7 | What are the biological findings? | Tissue transfer is uneven; pathway summaries can reduce coupled-label signals; model tiers behave differently | They should not expect one universal winner |
| 8 | What are extensions? | Organoids/species/single-cell are diagnostic/readiness lanes, not all equally frozen | They should not read extensions as validated external cohorts |
| 9 | What is released? | Metadata, source inventory, checks, and baseline artifacts can be separated from local mirrored data files | They should not confuse metadata readiness with payload readiness |

## Narrative Spine

The presentation should follow an assertion-evidence spine:

1. Public space omics has high-value data, but the evaluation surface is
   fragmented.
2. SpaceBio-Bench turns public studies into mission-held-out benchmark tasks.
3. The methodological trust hinge is train-only processing.
4. The first result: biological generalization is tissue-dependent.
5. The second result: pathway summaries can change what the benchmark sees.
6. The third result: model-family comparisons must be scoped carefully.
7. The extension result: organoids and other species are useful diagnostic
   lanes, but not a replacement for the core benchmark.
8. The platform result: v9 is becoming a resource layer with explicit release
   boundaries.

This spine should appear as short chapter claims, not as an agenda slide full
of section labels.

## Bridge Slides Needed Before Final Deck

### B1. Why A Benchmark Is Needed

Audience question:

> Why isn't this just another reanalysis of GeneLab studies?

Must include:

- public space omics is fragmented across studies, missions, tissues, and assay
  contexts;
- model evaluation needs a repeatable task surface;
- the project contributes a benchmark/resource layer, not only individual
  biological observations.

Visual content:

- left: scattered public studies/sources;
- middle: alignment into task records;
- right: auditable evaluation surface.

Premium consulting implementation:

- use an executive thesis headline, such as `Public space omics needs an
  evaluation layer`;
- avoid database screenshots or busy accession tables;
- show source diversity as a controlled source field, then compress into one
  clean benchmark rail;
- use a small bottom-right status note: `public-data resource; source-traceable`.

Do not include:

- raw file paths;
- full source inventory table;
- detailed model names;
- version-by-version internal history.

### B2. Study, Mission, Sample, Task

Audience question:

> What exactly becomes a benchmark task?

Must include:

- study/source accession;
- mission or flight context;
- samples with labels;
- tissue/assay context;
- task manifest as the auditable output.

Visual content:

- nested object diagram: public source -> mission -> samples -> labels ->
  task record;
- one simple example row, but not a full table.

Premium consulting implementation:

- use thin lines and nested containment, not boxes inside boxes;
- a single canonical example can sit as a small evidence strip;
- put definitions in an appendix note or speaker note, not as paragraphs on the
  slide.

Do not include:

- `task manifest` as the only visible phrase without explanation;
- every metadata column;
- accession soup that overwhelms the conceptual hierarchy.

### B3. One Mission Stays Hidden

Audience question:

> Why is mission-held-out evaluation stricter than random splitting?

Must include:

- training missions;
- hidden mission;
- repeated folds;
- no sample from the hidden mission is used in training.

Visual content:

- one mission physically separated behind a clean vertical boundary;
- a subtle repeat-fold indicator, not a dense fold table.

Premium consulting implementation:

- make the hidden mission the one dominant visual move;
- use a restrained red boundary line and muted training missions;
- headline should be direct: `The test set is a mission, not a random sample`.

Do not include:

- `LOMO` as the main title;
- cross-validation notation;
- fold IDs as the main visual object.

### B4. Train-Only Guard

Audience question:

> How do we know the hidden mission did not influence feature construction?

Must include:

- train missions feed feature selection/scaling/model fitting;
- hidden mission bypasses all fitting steps;
- hidden mission enters only at final scoring.

Visual content:

- two-lane schematic: train lane and hidden-test lane;
- a literal guard/wall before scoring;
- operation marks: choose features, scale, fit model, score.

Premium consulting implementation:

- use process architecture, not software architecture;
- keep operation labels as verbs;
- use one thin trust annotation: `learned from train missions only`.

Do not include:

- package names;
- class/function names;
- code-like arrows;
- a large table of preprocessing variants.

### B5. Genes And Pathways Are Views

Audience question:

> Why introduce pathway summaries, and what do they change?

Must include:

- gene-level matrix;
- pathway-level biological program summary;
- claim boundary that pathways help selected checks, not every task.

Visual content:

- many gene signals compress into a smaller set of pathway bands;
- a small paired-result inset can show the diagnostic change.

Premium consulting implementation:

- use a visual compression metaphor with scientific restraint;
- color genes and pathways consistently with existing semantics;
- headline should say `Pathways summarize biology; they do not guarantee
  better scores`.

Do not include:

- raw pathway acronym lists;
- `NES` unless the slide is explicitly a metric explanation;
- a pathway table masquerading as a figure.

### B6. Model Comparison Scope

Audience question:

> Are these model comparisons directly fair?

Must include:

- model families;
- shared task rows versus contextual rows;
- score metrics only where the figure needs them;
- direct comparison boundary.

Visual content:

- model family lanes feeding a common evaluation surface;
- a visible `shared rows` bracket;
- contextual rows in a lower-emphasis layer.

Premium consulting implementation:

- avoid leaderboard aesthetics;
- use a sober comparison surface;
- headline should say `Direct comparisons require the same task rows`.

Do not include:

- winner badges;
- oversized model logos;
- unbounded claims about foundation models.

### B7. Extension Lanes Are Diagnostic

Audience question:

> How do organoids, other species, and single-cell data fit into a mouse bulk benchmark?

Must include:

- core mouse bulk benchmark remains the anchored surface;
- human organoids are a biology-check extension;
- other species can be indexed as future/resource lanes;
- readiness differs by data type.

Visual content:

- core benchmark lane at center;
- extension lanes branching with explicit readiness/status markers;
- small status layer: `diagnostic`, `indexed`, `pending payload`, or similar.

Premium consulting implementation:

- use a portfolio/readiness map, not a table-heavy slide;
- keep status language conservative and visible;
- headline should say `Extensions widen coverage, but not all lanes carry the
  same claim`.

Do not include:

- `organoids validate the benchmark`;
- species lists as the primary visual;
- unsupported translational language.

### B8. Release Boundary

Audience question:

> What is ready now, and what remains gated?

Must include:

- metadata/source/task/checksum/baseline readiness;
- local mirrored data-file copy pending where applicable;
- what is safe to claim publicly;
- what is not yet a frozen payload release.

Visual content:

- light provenance-document grammar;
- release gate with ready and pending sides;
- exact claim language in a small, calm status layer.

Premium consulting implementation:

- look like a diligence memo, not a product dashboard;
- use thin rules, document surfaces, and traceability cues;
- headline should say `The resource boundary is explicit`.

Do not include:

- internal operating notes;
- raw local paths;
- release-decision backchannel text;
- UI cards.

## Slide-Level Content Template

Every new slide must be briefed with this template before rendering:

| Field | Required Answer |
|---|---|
| Audience question | What question does this slide answer for a first-time viewer? |
| Decision headline | What should the viewer believe after 5 seconds? |
| Evidence object | What source-verified figure/table/manifest supports the claim? |
| Visual move | What single spatial or narrative movement explains the point? |
| Bridge role | What previous slide does it connect from, and what next slide does it set up? |
| Visible caveat | What boundary prevents overclaiming? |
| Terms allowed on slide | Which technical terms are necessary and already explained? |
| Terms pushed to notes | Which accurate details belong in notes/caption instead? |
| Premium constraint | What would make this look cheap, busy, or dashboard-like? |

No slide should be generated unless these fields are filled.

## Consulting-Grade Style Principles

This deck should borrow the discipline of high-end consulting decks without
looking corporate or generic.

### 1. Assertion-Evidence, Not Decorative Title Plus Figure

Use a headline that states the conclusion, then make the figure prove it.

Good:

- `The test set is a mission, not a random sample`
- `Pathways summarize biology, but gains are task-specific`
- `The v9 resource boundary is explicit`

Weak:

- `Methods`
- `Model Results`
- `Organoid Analysis`

### 2. One Dominant Visual Move Per Slide

Each slide should have one movement:

- scattered sources compress into tasks;
- one mission moves behind a boundary;
- train operations stop before the hidden mission;
- genes compress into pathway bands;
- model lanes enter a shared-row evaluation surface;
- extension lanes branch away from the core benchmark;
- ready metadata stops at a data-file gate.

Multiple competing arrows, rings, or badges will make the slide feel cheap.

### 3. Thin-Rule Editorial Layout

Use:

- generous margins;
- thin divider rules;
- quiet grids;
- small status captions;
- direct labels;
- disciplined alignment.

Avoid:

- card boxes;
- rounded UI panels;
- heavy shadows;
- decorative icon rows;
- gradient blobs;
- dense labels around every object.

### 4. Evidence First, Decoration Last

The slide should contain a proof object or source-derived schematic before any
atmospheric treatment is added.

Allowed depth:

- faint measurement grid;
- tissue/cell texture at very low contrast;
- orbital/timeline cue when tied to mission holdout;
- document texture for provenance/release slides.

Not allowed:

- stock space wallpaper;
- generic sci-fi glow;
- background that competes with data;
- dark theme for a slide whose job is to explain a release boundary.

### 5. Status Language Is A Design Layer

Claim boundaries should look intentional, not apologetic.

Use small, precise status labels:

- `held-out mission evaluation`;
- `train-only feature processing`;
- `shared task rows only`;
- `biology-check extension`;
- `metadata-ready; data-file mirror pending`;
- `hypothesis-only`.

Avoid:

- `preliminary` as a vague shield;
- internal decision text;
- status paragraphs.

### 6. Tables Belong To Tables

Do not convert readiness matrices, source inventories, or claim-language
matrices into fake figures.

For the deck:

- use a schematic to teach the relationship;
- place exact rows in a clean appendix table if needed.

For the manuscript:

- main figure should show structure or quantitative pattern;
- exact release gates and source lists should become tables or supplement.

## Preferred Deck Architecture

### Opening

1. Executive thesis: SpaceBio-Bench creates an auditable evaluation layer for
   public space omics.
2. Problem bridge: public studies are valuable but fragmented.
3. Method bridge: public studies become mission-held-out tasks.

### Methods Trust

4. One mission stays hidden.
5. Train-only guard.
6. Genes and pathways are feature views.

### Results

7. Tissue transfer is uneven.
8. Pathway summaries change coupled-label signals and selected task behavior.
9. Model comparison is scope-dependent.

### Extensions And Resource

10. Organoids and species extensions are diagnostic/readiness lanes.
11. v9 platform/resource boundary is explicit.
12. What is ready, what remains gated, and what the next release will unlock.

This is the minimum coherent deck. Extra figures should go into appendix unless
they strengthen one of these transitions.

## Immediate Next Step

Do not build another polished slide yet.

Next work item:

1. Write content briefs for B1-B4.
2. Decide which existing figures/tables support each bridge.
3. Only then prototype one slide, starting with B3 or B4 because these are the
   highest-risk methods misunderstandings.

Recommended first prototype after briefing:

- B3 `The test set is a mission, not a random sample`.

Reason:

- it explains the core benchmark discipline;
- it bridges from data collection into result interpretation;
- it prevents a common reviewer misunderstanding before Fig1 appears.
