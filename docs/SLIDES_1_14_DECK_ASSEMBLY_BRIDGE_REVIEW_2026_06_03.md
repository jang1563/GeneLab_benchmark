# Slides 1-14 Deck Assembly Bridge Review

Date: 2026-06-03

## Scope

This review locks the first production bridge for slides 1-14 before building the actual presentation deck. The goal is to connect the opening thesis, methods bridge, feature view, early result family, and core result spine without drifting into the v8/v9 extension material too early.

## Outputs

- Assembly board: `output/slides_1_14_deck_assembly_bridge_v0_1/slides_1_14_deck_assembly_bridge.png`
- Grayscale QA board: `output/slides_1_14_deck_assembly_bridge_v0_1/qa/slides_1_14_deck_assembly_bridge_grayscale.png`
- Caption pack: `output/slides_1_14_deck_assembly_bridge_v0_1/slides_1_14_deck_assembly_caption_pack.md`
- Caption pack JSON: `output/slides_1_14_deck_assembly_bridge_v0_1/slides_1_14_deck_assembly_caption_pack.json`
- Overlay rules: `output/slides_1_14_deck_assembly_bridge_v0_1/slides_1_14_deck_overlay_rules.json`
- Manifest: `output/slides_1_14_deck_assembly_bridge_v0_1/slides_1_14_deck_assembly_manifest.json`
- QA: `output/slides_1_14_deck_assembly_bridge_v0_1/slides_1_14_deck_assembly_qa.json`

## Production Decision

The 14-slide default path should compress B2 into notes for slide 4 and B4 into notes for slide 5. This keeps the talk moving from "why this benchmark exists" to "how public studies become auditable tasks" to "what the model sees" before the first result block.

An expanded methods option remains valid if the audience is methods-heavy. In that case, B2 can become its own task-unit slide and B4 can become its own train-only guardrail slide. For the main deck, those details should live in notes, appendix, or a backup slide, not in the core 1-14 spine.

## Slide Block Findings

Slides 1-3 are production briefs, not final visuals. They need a title visual, external positioning matrix, and project evolution timeline. These should be original opening visuals rather than decorative reuse of result thumbnails.

Slides 4-6 have usable methods assets. Slide 4 carries the evaluation-layer bridge, slide 5 carries mission-held-out testing, and slide 6 defines feature layers before quantitative results appear.

Slides 7-9 can serve as the first result family: tissue hierarchy, pathway rescue, and model comparison. Their placement is correct because the audience has already seen the task and feature logic.

Slides 10-14 can serve as the hardened core result family: v4 hardening, temporal/RRRM guardrails, negative boundary, biological interpretation, and human translation. These slides should retain their claim-boundary lines.

## Style Caution

Slides 4-5 are light proof-stage scenes, while slide 6 onward uses a darker result grammar. This can work if the final deck treats the contrast as intentional: methods proof on a light paper-like surface, then results on dark analytical canvases. If tighter continuity is preferred, slides 4-5 should be rebuilt as dark variants before PPTX assembly.

## Overlay Rules

Use the existing rendered PNGs as Z2 evidence plates only. The final deck should rebuild the slide title, headline, bridge line, source, caveat, and presenter-facing text as editable objects above the scene.

Do not place table-like panels inside figure cards. If a dense table is needed, move it to a manuscript table, backup slide, or appendix table.

Do not introduce organoid, OSD-120, or other v8/v9 extension visuals before slide 20. The first 14 slides should keep the v1-v7 benchmark result core legible and bounded.

## First-Time-Viewer Bridge

The first-time viewer still needs explicit bridge language at four transition points:

- Slide 1 to 3: why a benchmark is needed and how the project evolved.
- Slide 3 to 4: how versions become an auditable evaluation workflow.
- Slide 6 to 7: why feature layers matter before comparing transfer results.
- Slide 14 onward: why partial human translation motivates, but does not prove, downstream v8/v9 extensions.

## QA

The first render had clipped opening-card headlines. The render script was revised to wrap long placeholder headlines and card headlines.

The regenerated assembly board was inspected visually in color and grayscale. Automatic QA reports 14 slides, no missing required assets, existing caption and overlay files, and 3840 x 2160 PNG outputs.

Current QA status: pass for assembly planning. The board is not a final slide deck; it is a production map for building the actual PPTX.

## Next Step

Before PPTX assembly, choose one of two production paths:

- Build original premium opening visuals for slides 1-3 first.
- Decide whether slides 4-5 remain light proof-stage scenes or are rebuilt as dark methods scenes for tighter visual continuity.
