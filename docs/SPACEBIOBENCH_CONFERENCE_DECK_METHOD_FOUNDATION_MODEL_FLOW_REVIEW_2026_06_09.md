# SpaceBio-Bench Conference Deck Method and Foundation-Model Flow Review

Date: 2026-06-09

## Scope

This pass reviewed whether the conference deck explains the methodology,
analysis workflow, ML comparison, and foundation-model comparison clearly enough
for an audience that may not know ML or LLM terminology.

Reviewed and edited deck:

- `output/spacebiobench_conference_deck_v0_5/spacebiobench_conference_deck_v0_5.pptx`

## Finding Before Edit

Slides 4-7 already made the benchmark mechanics understandable:

- Slide 4: public repository records become a task contract.
- Slide 5: held-out means a hidden mission, not random sample validation.
- Slide 6: models receive numerical views, not biology directly.
- Slide 7: feature views feed one held-out AUROC score.

The gap was between slide 7 and the later model-comparison result. The deck did
not yet show how `compressed input` becomes a comparison among classical ML,
foundation models, and text LLM checks. That made slide 11's foundation-model
result too dependent on verbal explanation.

## Change Made

Added a new native methods slide:

- Slide 8: `Model families share one task`

The new slide inserts a visual bridge before the result section:

`fixed split/view/metric -> model-family tiers -> one result surface`

The model-family tiers are:

- `CLASSICAL ML`: task-trained models; PCA-LR / RF / LR use gene or pathway features
- `FOUNDATION MODEL`: pretrained, then adapted; scGPT / Geneformer start from expression pretraining
- `TEXT LLM CHECK`: prompt-only language model; Gemini / Llama / DeepSeek see prompts, not matrices

The slide also adds the explicit reading rule:

`this is a controlled model-family stress test, not a universal intelligence ranking`

June 10 glossary update:

- Slide 8 was revised from model-name-first labels to definition-first labels.
- Speaker notes now define the three model families before naming the fixed
  split/view/metric comparison surface.
- Slide 11's model guide was revised to echo the same vocabulary:
  `families: classical / FM / text LLM`.

## Visual Review

Pass after iteration.

- The slide is readable at full-slide scale.
- The first render had the bottom reading rule too close to the `TEXT LLM CHECK`
  card; this was corrected.
- The model-tier body text initially ran too close to the outgoing arrows; it was
  shortened and moved under each tier title.
- The final glossary version puts the definition in the large text and the model
  examples in the smaller supporting line.
- The final visual hierarchy is clear: left setup, middle model tiers, right
  result surface, bottom reading rule.

## Narrative Flow

Pass.

The updated methods sequence is now:

1. Data becomes a task contract.
2. Held-out evaluation hides a mission.
3. Models see numerical views of samples.
4. Feature views become one held-out score.
5. Model families are compared only after split, view, and metric are fixed.
6. Result slides can then be read as local benchmark evidence.

This makes the later `Scale alone does not transfer` slide easier to understand
without needing a long spoken aside about foundation models.

## Claim Safety

Pass.

- The new slide is explicitly a method scaffold.
- It does not claim foundation-model failure in general.
- It separates prompt-only text LLM checks from expression-matrix model training.
- It frames the comparison as tested-setting evidence, not an AI ranking.

## Residual Notes

- Slide 11's embedded model-comparison figure remains small, but the guide card
  and new slide 8 reduce the burden on that figure. The slide 11 guide now
  explicitly links the result back to the slide 8 model-family glossary.
- No additional method slide is recommended right now. Adding more would likely
  slow the conference rhythm more than it helps.

## QA

PPTX XML smoke:

- Slide count: 26
- Notes count: 26
- New slide text present:
  - `Model families share one task`
  - `task-trained models`
  - `pretrained, then adapted`
  - `prompt-only language model`
  - `Gemini / Llama / DeepSeek`
  - `classical / FM / text LLM`
- Visible forbidden phrase hit for `single-cell leaderboard`: none
- Extension term hits:
  - slide 23: `organoid`
  - slide 24: `OSD-120`

PDF export:

- Pages: 26
- Page size: 960.009 x 540 pt
- PDF version: 1.7
- Tagged: yes
- Encrypted: no
- JavaScript: no
- File size: 6,292,695 bytes

Output sizes:

- PPTX: 28,915,825 bytes
- PDF: 6,292,695 bytes
- Contact sheet: 3,788,646 bytes
- Speaker notes markdown: 15,999 bytes

## Decision

Keep the new slide 8.

The methods flow now explains the ML and foundation-model comparison visually
enough for a mixed audience before results begin.
