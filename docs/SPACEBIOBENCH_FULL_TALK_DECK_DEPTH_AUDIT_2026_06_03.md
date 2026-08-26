# SpaceBio-Bench Full Talk Deck Depth Audit

Date: 2026-06-03

## Bottom Line

The v0.3 full-talk deck is strong on claim discipline but still too thin as a
teaching artifact. It repeatedly tells the audience what not to overclaim, but
it does not yet give a first-time viewer enough scaffolding to understand the
data units, evaluation recipe, metric reading, and result plots without a strong
live narration.

The current deck is not scientifically shallow in intent. The problem is that
too much methodological depth is parked in speaker notes, footers, or small
proof-object labels. For a specialist audience this may be acceptable. For a
mixed biology/AI/NASA-style audience, slides 4-6 and the first result slides
need more explicit on-slide pedagogy.

## Audience Assumption

The next pass should assume a smart first-time viewer who may know one of these
areas, but not all three:

- space biology data repositories such as GeneLab/OSDR
- ML benchmark construction and held-out evaluation
- transcriptomic feature engineering and pathway summaries

That viewer needs the deck to answer four questions early:

- What is one data record in this benchmark?
- What exactly is hidden from the model?
- What feature representation is being compared?
- How should I read AUROC, tissue transfer, and negative results?

## Current Strengths

- The claim boundaries are unusually clear and should be preserved.
- The opening story is coherent: benchmark niche, methods bridge, result spine,
  platform/extension boundary, roadmap.
- The visual rhythm works at contact-sheet scale.
- Slides 10-14 have a stronger result narrative than the earlier methods
  slides.
- The v8/v9 extension section is careful about alpha, hypothesis-only, and
  blocker status.

## Core Weakness

The deck currently treats methods as visual atmosphere plus presenter narration.
That creates a depth gap:

- Slide 4 says public studies become auditable tasks, but the audience does not
  see a concrete task record with fields.
- Slide 5 says the hidden test set is a mission, but does not show the leakage
  rule in operational form.
- Slide 6 says feature views can change, but does not show what gene-level,
  pathway-level, PCA, or model inputs mean for an actual sample.
- Slide 7 introduces AUROC/tissue transfer before the audience has been taught
  how to read the axis, confidence interval, chance line, sample count, or
  mission-held-out split.
- Later platform/extension slides are boundary-correct but assume the viewer
  already understands why manifests, hashes, h5ad, STARsolo, and payload freeze
  matter.

## Slide-Level Audit

| Slides | Current Role | Depth Risk | Recommended Fix |
| --- | --- | --- | --- |
| 1-3 | Opening thesis and positioning | Good macro story, but GeneLab/OSDR and `v1-v9` remain undefined for newcomers. | Add a compact glossary strip or one-line definitions in slide 2 or 3. |
| 4 | Task construction | Too schematic. It looks like documents turning into boxes, not a benchmark data contract. | Replace or augment with a concrete task card: source study, tissue, assay, samples, labels, train/test rule, metric, provenance fields. |
| 5 | Mission-held-out guard | Strong title, but the actual anti-leakage protocol is not visible enough. | Show four-step protocol: choose held-out mission, freeze preprocessing/model selection, train on other missions, score hidden mission once. |
| 6 | Feature layers | Visually clear but under-explained. | Add labeled lanes: gene expression matrix, pathway scores, model input, classifier, AUROC. Define that the same samples are evaluated under different feature views. |
| 7 | Tissue transfer | Result arrives before the audience has a metric primer. Labels are too small for novice reading. | Add a `How to read this plot` callout: dot is AUROC, line is uncertainty, 0.5 is chance, each row is a tissue. |
| 8-9 | Pathway and model comparison | Interesting but likely too compressed. | Add a visible interpretation ladder: observed pattern, what it supports, what it does not support. |
| 10 | Hardening grid | Stronger because it shows scale, counts, and grid logic. | Keep, but define `8 tissues x 8 classifiers x 4 feature views = 256 evaluations` more explicitly as a method/result bridge. |
| 11 | Timepoint/RRRM | Dense domain-specific terms. | Add a small decoder: RR-8, RR-6, RRRM, preservation, recovery projection. |
| 12 | Negative boundary | Good boundary slide, but novice may not know what was tested. | Add a small panel: tested models, tested task, observed negative pattern. |
| 13 | Biological interpretation | The hypothesis boundary is good. | Add source-to-hypothesis chain: signal -> pathway/target tier -> follow-up hypothesis. |
| 14 | Human translation | Good caution, but may need a species/feature bridge. | Show mouse signal, human perturbation signal, overlap layer, and non-overlap layer as a four-part visual. |
| 15-17 | v7/v8 boundary | Boundary clear, but v8 SpaceMed may sound sudden. | Add a transition note: "Why translation appears after benchmark evidence." |
| 18-19 | v9 platform/public bulk | Too infrastructure-heavy for newcomers. | Add definitions for manifest, evaluator, run record, checksum, payload hash. |
| 20-22 | Extension tracks | Honest and useful, but terms like h5ad/STARsolo appear late without decoding. | Add a small `What is blocked?` explanation: raw data exists, processed analysis-ready matrix does not. |
| 23-24 | Close and roadmap | Clear closing guardrails. | Keep, but make roadmap map to three audiences: data release, benchmark paper, extension work. |

## Needed Conceptual Inserts

The deck needs at least three explicit teaching moments. They can be inserted
as new slides in a longer seminar deck or folded into existing slides in a
24-slide conference deck.

### 1. Data Orientation

Add a beginner-friendly "what data are we using?" slide or panel.

Recommended visual:

`Repository -> Study -> Mission -> Tissue -> Sample -> Expression Matrix -> Task`

On-slide definitions:

- GeneLab/OSDR: public NASA space biology data repositories.
- Study/mission: biological experiment context that must stay attached to the
  samples.
- Tissue: organ or sample context; transfer differs by tissue.
- Expression matrix: samples by molecular features.
- Label: condition being predicted, such as flight versus ground or timepoint.

### 2. Evaluation Recipe

Add a concrete task recipe before the first result.

Recommended visual:

`Train missions -> frozen choices -> held-out mission -> score -> source record`

Required on-slide labels:

- Train only: feature selection, model choice, preprocessing.
- Hidden mission: no tuning after boundary.
- Score: AUROC or task-specific metric.
- Audit trail: source study, sample IDs, task ID, model, feature view, result.

### 3. Metric/Plot Primer

Add one small primer before slide 7 or build it into slide 7.

Recommended visual:

- AUROC 0.5 = chance.
- Higher AUROC = better separation for the tested task.
- Dot = score.
- Horizontal line = uncertainty or interval.
- Each row = one tissue or task family.
- A result can be strong in one tissue and weak in another.

## Methodology Visualization Recommendations

The methods section should become less abstract and more operational.

Recommended replacement sequence:

1. `Public data becomes a task record`
   - Show a native editable table/card with fields: source, assay, tissue,
     samples, labels, split, features, metric, provenance.
2. `Held-out means mission-held-out`
   - Show train missions on the left, a vertical freeze boundary, hidden
     mission on the right, and a single score output.
3. `Same samples, different feature views`
   - Show gene matrix and pathway summary as parallel inputs feeding the same
     evaluator.
4. `How to read the first result`
   - Either an inserted primer slide or a redesigned slide 7 with annotated
     axis, chance line, score dot, and uncertainty line.

This would make the result spine feel earned instead of abrupt.

## Data Explanation Recommendations

Each data/result slide should carry a small "data card" in a consistent place.

Suggested fields:

- Dataset or source family
- Tissue or biological context
- Task label
- Split rule
- Feature view
- Metric
- Claim boundary

This can be very compact, but it should be visible on the slide, not only in
speaker notes.

## Extension Section Recommendations

Slides 18-22 are honest but need one level of translation for mixed audiences.

Add a compact infrastructure glossary across slides 18-22:

- Manifest: the list of expected source files and task records.
- Evaluator: the code path that computes benchmark metrics.
- Run record: stored output from a completed evaluation.
- Checksum/hash: proof that a local file matches the expected payload.
- h5ad/STARsolo: processed single-cell output formats needed before scoring.
- Payload blocker: a missing or unverifiable analysis-ready data object.

The key teaching point should be:

`metadata readiness is not the same as payload readiness`

Slide 19 already says this, but slide 18 should prepare the audience to
understand why it matters.

## Recommended Next Pass

Create a `v0.4 depth pass` rather than only applying institutional branding.
The pass should preserve the current visual system but deepen the explanatory
surface.

Priority order:

1. Rewrite slides 4-6 as native editable method explainers, not atmospheric
   schematics.
2. Add a plot-reading scaffold to slide 7.
3. Add data cards to slides 7-14.
4. Add infrastructure definitions to slides 18-22.
5. Expand speaker notes into audience-facing mini-scripts for slides 4-8 and
   18-22.

## Suggested Deck Variants

### 24-Slide Conference Version

Keep 24 slides. Replace the least explanatory method visuals with clearer
native diagrams. Do not add many new slides.

Best for:

- 20-25 minute talk
- already technical audience
- conference presentation

### 30-Slide Seminar Version

Add six teaching slides:

- data orientation
- task record anatomy
- mission-held-out protocol
- feature-view primer
- AUROC/plot-reading primer
- platform/payload glossary

Best for:

- lab seminar
- thesis/job-style talk
- mixed audience across biology and AI

## Final Judgment

The current deck is not too shallow for a fast expert talk, but it is too
shallow for a first-time mixed audience. It shows the boundary of each claim
better than it teaches the machinery behind the claim. The next high-value
iteration is therefore not more polish. It is a depth pass: make the data units,
evaluation protocol, metric reading, and platform blockers visible on the
slides themselves.
