# Visual Layered Scene Blueprint

Date: 2026-06-01

Purpose: define how to convert the current source-verified figures into
premium slide scenes without losing scientific trust. This is the slide-deck
companion to `docs/VISUAL_PREMIUM_QUALITY_STANDARD_2026_05_31.md`.

## Bottom Line

For the 2026 ISGP-grade deck, the target is:

> Layered PNG scene + editable scientific interpretation.

The current manuscript variants are clean proof objects. The next slide pass
should not simply paste those figures onto white slides. Each major slide
should add depth through a controlled z-stack:

- Z0-Z1: context and measurement atmosphere;
- Z2: source-verified evidence object;
- Z3-Z5: editable interpretation, trust boundary, and focus/motion.

## Deliverable Contract

Each premium slide should have three files or file-equivalents:

1. `scene_plate.png`
   - raster scene including Z0-Z2;
   - may include background texture, assay lanes, grid, proof object, and
     subtle backplate;
   - should not include long explanatory text.
2. `editable_overlay`
   - slide-native text, callouts, focus rings, labels, and caveats;
   - must remain editable until the final deck export.
3. `layer_manifest.json`
   - source figure path;
   - source table or manifest path;
   - layer list;
   - claim boundary;
   - color/accessibility QA status.

The manuscript figure remains a separate static PNG/PDF. The slide scene is a
communication layer around that proof object, not a replacement for it.

## Layer Definitions

| Layer | Slide role | What belongs here | What does not belong here |
|---|---|---|---|
| Z0 canvas / atmosphere | Establish domain and depth | dark field, subdued tissue/cell texture, faint evidence-paper texture | stock-looking space wallpaper, decorative blobs, busy photography |
| Z1 measurement layer | Reveal the experimental frame | orbital arc, held-out mission tick, assay lane, sample/tissue lane, faint grid | dense methods diagram, decorative icons |
| Z2 evidence surface | Carry the audited result | high-res figure crop, chart proof plate, source-derived image, shaded backplate | unsourced screenshots, editable text baked into the PNG |
| Z3 interpretation layer | Guide attention | focus ring, highlighted transfer gap, arrow, two- to six-word callout | paragraphs, method details, reviewer notes |
| Z4 trust/status layer | Set scope | short caveat, public/draft/preliminary status, source family label | raw file paths, allowed/prohibited language, internal decision text |
| Z5 motion/focus layer | Create one visual move | trajectory line, crop window, reveal mask, highlighted break | multiple simultaneous animation ideas |

## Figure-Specific Blueprints

### Fig1 Tissue Transfer

Evidence object:

- `output/premium_figures/premium_fig1_tissue_transfer_hierarchy.png`

Slide scene:

- Z0: dark, quiet transcriptomic/cellular field with very low contrast.
- Z1: mission-held-out timeline ticks or orbital arc behind the plot.
- Z2: lollipop/interval proof plate cropped to the ranking region.
- Z3: focus ring around thymus and liver, with one short callout:
  `Some tissues transfer; others do not`.
- Z4: small scope label:
  `Mouse bulk RNA-seq; held-out mission evaluation`.
- Z5: one highlighted transfer break from thymus/gastrocnemius to liver/kidney.

Risk:

- Do not over-space-theme this slide. The audience should read tissue
  generalization first, spaceflight atmosphere second.

### Fig2 Pathway Summaries

Evidence object:

- `output/premium_figures/manuscript_variants/premium_fig2_pathway_rescue_manuscript.png`

Slide scene:

- Z0: faint gene-to-pathway texture or pathway-node field.
- Z1: assay lane showing `gene-level input` to `pathway summary`.
- Z2: proof plate centered on Panel A and the detection-task change panel.
- Z3: highlight the paired drop in coupled-label signal and the selected task
  rescue region.
- Z4: scope label:
  `Coupled experimental labels; diagnostic check`.
- Z5: one transition line from noisy gene-level labels to cleaner pathway
  summary.

Risk:

- Avoid implying that pathways universally improve every task. The visual move
  should say `reduce unwanted label signals` plus `selected task gains`.

### Fig3 Model Comparison

Evidence object:

- `output/premium_figures/manuscript_variants/premium_fig3_model_tier_comparison_manuscript.png`

Slide scene:

- Z0: dark benchmark grid, not a product-dashboard background.
- Z1: three model-family lanes: classical baseline, single-cell model,
  text LLM.
- Z2: proof plate with the dot/lollipop model comparison.
- Z3: focus ring on PCA-LR and the negative tissue deltas.
- Z4: scope label:
  `Shared 6-task rows are direct comparisons`.
- Z5: one break line at the chance/reference region or one reveal from
  aggregate score to tissue-specific deltas.

Risk:

- Do not make model family colors decorative. The point is benchmark surface
  fairness and small-sample domain shift, not a leaderboard aesthetic.

### Fig6 Human Organoid Biology Check

Evidence object:

- `output/premium_figures/manuscript_variants/premium_fig6_organoid_biology_check_manuscript.png`

Slide scene:

- Z0: subdued organoid/cellular texture. Use generated or open-license visual
  texture unless a source image has clear reuse permission.
- Z1: RNA-seq assay lane and small source/sample tick marks.
- Z2: proof plate with the dataset footprint, metric strip, and enrichment
  panel.
- Z3: focus ring around the separation between primary prediction and secondary
  biology checks.
- Z4: scope label:
  `Small public human organoid dataset; source factors coupled`.
- Z5: crop window that moves from primary AUROC/macro-F1 to gene-enrichment
  check.

Risk:

- Do not visually over-promote organoids as a finished validation cohort. It is
  an extension biology-check dataset.

### v8 Translational Incubator

Evidence object:

- final v8 intervention or bridge figure to be selected in the deck build pass.

Slide scene:

- Z0: muted biological pathway/stressor field.
- Z1: mission or stressor lane.
- Z2: source-verified pathway/stressor/intervention proof object.
- Z3: hypothesis focus ring.
- Z4: mandatory scope label:
  `Hypothesis-only; not a countermeasure recommendation`.
- Z5: one transition from benchmark signal to candidate hypothesis.

Risk:

- This slide needs the strongest trust/status layer. The design must look
  exciting without sounding validated.

### v9 Platform / Public Bulk Alpha

Evidence object:

- `output/premium_figures/premium_fig4_v9_platform_architecture.png`
- `output/premium_figures/premium_fig5_public_bulk_release_boundary_schematic.png`

Slide scene:

- Z0: subdued evidence-paper or manifest texture.
- Z1: registry/evaluator/provenance lane.
- Z2: architecture or release-boundary proof plate.
- Z3: callout around the public-bulk alpha path.
- Z4: scope label:
  `Metadata indexed; local data-file copy pending`.
- Z5: one path from source inventory to task manifest to baseline check.

Risk:

- Avoid product UI. This should look like a scientific resource figure with a
  trust trail, not a SaaS dashboard.

## Slide Build Workflow

For each slide:

1. Choose the one-sentence claim.
2. Select the proof object and verify its source manifest.
3. Sketch Z0-Z5 in a short layer brief.
4. Render or assemble the scene plate.
5. Add editable overlay labels and caveats.
6. Export a full-slide PNG/PDF.
7. Inspect at full screen and thumbnail size.
8. Run grayscale/color-vision QA if the slide carries quantitative meaning.
9. Record the layer manifest and final verdict.

## Acceptance Criteria

A slide reaches premium layered-scene status only if:

- the evidence object remains readable without the background;
- the background improves context but does not steal attention;
- the main claim is visible in under five seconds;
- the caveat/status layer is present but not visually dominant;
- the figure can be explained from the slide without reading speaker notes;
- the editable overlay can be revised without regenerating the proof plate;
- no internal process language appears in visible slide text;
- the proof object still links to source data and audit results.

## Pilot Record

Completed layered-scene pilots:

- `docs/VISUAL_LAYERED_SCENE_PILOT_FIG1_2026_06_01.md`
- `docs/VISUAL_LAYERED_SCENE_PILOT_FIG2_2026_06_01.md`
- `docs/VISUAL_LAYERED_SCENE_PILOT_FIG3_2026_06_01.md`
- `docs/VISUAL_LAYERED_SCENE_PILOT_FIG6_2026_06_01.md`
- `docs/VISUAL_LAYERED_SCENE_FAMILY_QA_FIG1_FIG6_2026_06_01.md`
- `docs/VISUAL_V9_PROVENANCE_DOCUMENT_GRAMMAR_2026_06_01.md`
- `docs/VISUAL_PREMIUM_DECK_ORDER_QA_2026_06_01.md`

Current lesson:

- Fig1, Fig2, Fig3, and Fig6 can share one result-slide grammar: dark Z0
  field, quiet Z1 measurement lane, source-derived Z2 proof crop, and editable
  Z3-Z5 interpretation.
- Fig2 showed that large rings quickly look coarse on quantitative panels.
  Corner focus windows and a single motion path are a better default.
- Fig3 showed that crop repair is a premium-quality issue, not a mechanical
  export issue. Long model names, panel titles, and value labels must remain
  readable before any overlay is judged.
- Fig6 can use the same result-slide grammar if the caution layer is stronger.
  `Organoids add biology checks` is acceptable; `organoids validate` is not.
- The Fig1-Fig6 family contact sheet passes as a draft main-deck sequence.
  The right-side evidence layout is repetitive but acceptable for the result
  block; v9 platform slides should decide whether to switch grammars.
- v9 platform/resource slides should switch grammars. Use a light
  provenance-document canvas for Fig4/Fig5 so release boundaries and pending
  data-file gates are visually primary.
- The deck-order contact sheet confirms the grammar transition: dark result
  slides first, light v9 provenance-document slides second.
- Jargon should move down into the proof panel or caption. The visible slide
  claim should be reader-facing.

Next move:

- before final deck export, rebuild Z3-Z5 as editable slide-native objects and
  re-run full-screen, thumbnail, grayscale, and color-vision QA.
- optional next visual block: build one Fig1 corner-window alternate and one
  v9 public-bulk zoomed variant, then compare against the current candidates.
