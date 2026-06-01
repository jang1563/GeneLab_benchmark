# Visual Methods Explanation Gap Map

Date: 2026-06-01

Purpose: identify where a first-time viewer may fail to understand the
SpaceBio-Bench methods, and prepare slide-level explanations before building
the final deck.

## Audience Problem

An expert in genomics may know RNA-seq and AUROC but not the spaceflight data
structure. A space-biology viewer may know OSDR/GeneLab but not benchmark split
discipline. A general AI/ML viewer may know benchmarks but not why mission is
the independence unit.

Therefore the deck needs methods visuals that show:

1. where the data come from;
2. what is held out;
3. how leakage is prevented;
4. what models see;
5. what scores mean;
6. what is released versus still pending.

## Explanation Gaps

| Gap | Why It Confuses New Viewers | Visual Fix |
|---|---|---|
| Mission-held-out evaluation | `LOMO` sounds like internal jargon; viewers may assume random CV | Show one mission physically hidden while all other missions train |
| Study/sample/source relationship | OSDR/GeneLab accessions, missions, samples, and tissues are nested | Use a source-to-sample-to-task flow, not a table-only explanation |
| Leakage prevention | Train-only filtering and feature selection are invisible in result plots | Add a guard rail: transformations learned only from training missions |
| Gene versus pathway features | Viewers may not see why pathway summaries matter | Show genes as many noisy signals compressed into pathway-level biology |
| Model comparison boundary | Model rows use shared task subsets and sometimes different task sets | Mark `shared rows only` and keep subset caveats near model slides |
| Confounder diagnostics | Mission/hardware/gravity label checks can sound like the main task | Label these as diagnostic labels, separate from flight prediction |
| Organoid extension | Human organoids can look like validation if not bounded | Keep `biology checks`, `small public set`, and `source factors coupled` visible |
| v8 intervention layer | Perturbation candidates can be misread as countermeasures | Use `hypothesis-only`, not treatment or recommendation language |
| v9 metadata alpha | Metadata readiness can be mistaken for a frozen data payload | Separate metadata records from data-file mirror in a release-boundary visual |
| Single-cell track | Indexed assets may look runnable | Show blockers and missing canonical payload staging explicitly |

## Required Methods Explainer Slides

### 1. Data-To-Task Map

Question it answers:

> How do public GeneLab/OSDR studies become benchmark tasks?

Must show:

- public studies;
- sample tables and expression matrices;
- tissue/mission labels;
- task manifests;
- one held-out mission split.

Primary visual:

- full-width process rail with one hidden mission as the key visual move.

### 2. Leakage Guard

Question it answers:

> How do we avoid accidentally learning from the test mission?

Must show:

- train missions only;
- variance filtering or feature selection inside each fold;
- scaler/model fit on training only;
- test mission touched only at scoring.

Primary visual:

- two-lane train/test split with a visible wall between training operations and
  held-out scoring.

### 3. Feature-Layer Map

Question it answers:

> What does gene versus pathway mean in this benchmark?

Must show:

- gene expression matrix;
- pathway activity summary;
- diagnostic label suppression;
- selected tissue/task gains.

Primary visual:

- gene cloud compressed into pathway bands, with a small note that pathway
  summaries do not universally improve every task.

### 4. Model And Metric Surface

Question it answers:

> What is being compared, and when are comparisons direct?

Must show:

- classical baseline;
- single-cell/gene-expression foundation models;
- text LLMs;
- shared-row caveat;
- AUROC/macro-F1 as score readouts.

Primary visual:

- model families as lanes feeding a shared evaluation table, not as a product
  leaderboard.

### 5. Platform And Release Boundary

Question it answers:

> What exists now, and what is still pending?

Must show:

- task/source/checksum/baseline metadata present;
- local data-file mirror pending;
- extension lanes diagnostic;
- single-cell blockers.

Primary visual:

- light provenance-document grammar, not dark result-slide grammar.

## Caption/Notes Prep

Terms that need first-use expansion:

- `LOMO`: one mission held out.
- `task manifest`: machine-readable task definition.
- `source inventory`: public data source table.
- `checksum evidence`: source-side file identity evidence.
- `payload`: local mirrored data files.
- `pathway summary`: sample-level biological program activity.
- `diagnostic label`: a label used to test confounding, not the main prediction
  endpoint.

Terms to avoid in visible slide text unless necessary:

- `artifact rescue`;
- `payload freeze`;
- `RRRM scaffold`;
- `NES`;
- `macro-F1` without nearby context;
- `alpha snapshot` without `metadata-only`.

## Immediate Visual Build

Build one first-time-viewer overview slide:

> Public studies -> clean sample table -> hide one mission -> build
> gene/pathway features -> train models -> score and audit.

The slide should use plain verbs and keep technical terms in smaller secondary
labels. It should not look like a dashboard or a software architecture diagram.
