# Biological Symbol Module Pack v0.3 Review

Date: 2026-06-02

Purpose:

- continue v0.3 biological visual-system work beyond slide mockups;
- create reusable symbols, icons, badges, and small modules for premium
  consulting-style materials;
- make biological method explanations easier to understand without turning
  every slide into a dense table or generic box-arrow diagram.

## Output

Asset root:

- `assets/biovis_symbol_module_pack_v0_3/`

Generator:

- `scripts/build_biovis_symbol_module_pack_v3.py`

Generated files:

- `assets/biovis_symbol_module_pack_v0_3/symbols/svg/`
- `assets/biovis_symbol_module_pack_v0_3/symbols/png/`
- `assets/biovis_symbol_module_pack_v0_3/badges/svg/`
- `assets/biovis_symbol_module_pack_v0_3/badges/png/`
- `assets/biovis_symbol_module_pack_v0_3/modules/svg/`
- `assets/biovis_symbol_module_pack_v0_3/modules/png/`
- `assets/biovis_symbol_module_pack_v0_3/modules_dark/svg/`
- `assets/biovis_symbol_module_pack_v0_3/modules_dark/png/`
- `assets/biovis_symbol_module_pack_v0_3/preview/biovis_symbol_module_pack_v0_3_symbols.png`
- `assets/biovis_symbol_module_pack_v0_3/preview/biovis_symbol_module_pack_v0_3_badges_modules.png`
- `assets/biovis_symbol_module_pack_v0_3/preview/biovis_symbol_module_pack_v0_3_dark_modules.png`
- `assets/biovis_symbol_module_pack_v0_3/qa/biovis_symbol_module_pack_v0_3_symbols_grayscale.png`
- `assets/biovis_symbol_module_pack_v0_3/qa/biovis_symbol_module_pack_v0_3_badges_modules_grayscale.png`
- `assets/biovis_symbol_module_pack_v0_3/qa/biovis_symbol_module_pack_v0_3_micro_symbol_qa.png`
- `assets/biovis_symbol_module_pack_v0_3/qa/biovis_symbol_module_pack_v0_3_qa.json`
- `assets/biovis_symbol_module_pack_v0_3/manifest.json`
- `assets/biovis_symbol_module_pack_v0_3/manifest.csv`

## Asset Count

Generated:

- 36 biological and trust symbols;
- 8 evidence/status badges;
- 7 reusable method/scope/trust modules;
- 7 dark-field module variants;
- 58 total editable SVG assets, each with PNG preview.

Automatic QA:

- SVG assets created: pass;
- PNG previews created: pass;
- manifest JSON/CSV: pass;
- grayscale preview sheets: pass;
- dark module contact sheet: pass;
- micro-size symbol QA sheet: pass;
- contact sheet rendering: pass after switching preview PNG generation from
  nested SVG references to Pillow-composited PNG sheets.

## v0.3 Hardening Pass

After the initial symbol/module pack, three weaknesses were addressed:

- dark-field modules were added because transparent light-module PNGs had dark
  text and were not suitable on deep slide canvases;
- species/model-system symbols were simplified toward restrained silhouettes,
  especially mouse, rat, and fly, to reduce playful/cartoon reading;
- micro-size QA was added because normal contact sheets overestimate how well
  symbols survive in figure footers, legends, and compact method bars.

This hardening pass does not change the purpose of the pack. It makes the pack
safer to use in actual premium decks.

## Why This Pack Was Needed

The previous symbol pack v0.1 was useful but too close to a simple icon set.
It lacked:

- species/model-system coverage;
- status badges for source/proof/caveat semantics;
- reusable modules for explaining methods and scope;
- explicit allowed/prohibited-use metadata;
- a way to distinguish biological signposts from proof objects.

The v0.3 motif pack solved biological texture and atmosphere. This v0.3 symbol
module pack solves the smaller but equally important problem: what goes in
headers, lanes, badges, method bridges, figure annotations, and recurring deck
controls.

## Design System

Three layers are now separated:

| Layer | Purpose | Example Use |
|---|---|---|
| Symbols | small biological or trust signposts | tissue, RNA, organoid, mouse, split guard |
| Badges | evidence/status tags | source, processed, validated, caveat, train-only |
| Modules | reusable explanatory components | sample-to-feature, species coverage, claim boundary |

Design rule:

> Use these as semantic visual vocabulary, not as source-derived evidence.

## Symbol Coverage

Molecular:

- `gene_locus`
- `rna_readout`
- `protein_program`
- `pathway_network`

Cellular:

- `mitochondrial_stress`
- `cell_state`
- `cell_population`
- `organoid_rosette`

Tissue and organ:

- `tissue_section`
- `spatial_spots`
- `vascular_context`
- `epithelial_barrier`
- `brain_organ`
- `lung_organ`
- `kidney_organ`
- `liver_organ`

Assay/modality:

- `sample_tube`
- `bulk_rna_assay`
- `single_cell_droplet`
- `proteomics_assay`
- `microscopy_frame`
- `omics_matrix`

Species/model systems:

- `human_model`
- `mouse_model`
- `rat_model`
- `fly_model`
- `worm_model`
- `fish_model`
- `plant_model`
- `microbe_model`

Space-biology context:

- `orbit_microgravity`
- `radiation_stressor`

Trust/provenance:

- `source_record`
- `checksum_record`
- `split_guard`
- `caveat_flag`

## Badges

Generated status badges:

- `generated_context`: schematic visual context only;
- `source_proof`: real source-derived proof object;
- `processed_result`: computed or processed output;
- `validated_result`: reviewed/validated claim surface;
- `preliminary`: exploratory or draft analysis;
- `caveat`: explicit limitation;
- `train_only`: train-only processing;
- `held_out`: held-out evaluation split.

Recommended deck use:

- use one or two badges per panel, not a row of many badges;
- keep caveat/source/proof badges near the figure they qualify;
- do not let badges become decorative tags detached from claims.

## Modules

Generated reusable modules:

- `biological_scale_ladder`: molecule to cell to tissue to organ to model;
- `sample_to_feature_stack`: source material to assay to matrix to program to
  model task;
- `species_coverage_strip`: human/organoid/mouse core plus extension systems;
- `modality_lane_set`: bulk, single-cell, spatial, proteomics, imaging,
  organoid lanes;
- `trust_status_ribbon`: source, freeze, split, result, caveat;
- `claim_boundary_bar`: schematic, source, processed, validated;
- `space_bio_assay_lane`: mission/stressor to tissue to assay to feature layer
  to benchmark guard.

Each module now has:

- a light-surface SVG/PNG in `modules/`;
- a dark-field SVG/PNG in `modules_dark/`.

Immediate high-value modules:

- `sample_to_feature_stack`: best for methods explanation;
- `claim_boundary_bar`: best for preventing overclaiming;
- `species_coverage_strip`: useful for mouse/human core plus organoid/species
  expansion framing;
- `space_bio_assay_lane`: useful for SpaceBio-Bench-specific bridge slides.

## Visual QA

Contact sheet:

- pass: symbols are visible and grouped by semantic family;
- pass: molecular, cell, assay, organ, and trust symbols read clearly;
- pass: preview contact sheet correctly displays icons after switching to
  Pillow-composited PNG preview generation;
- conditional pass: species icons are useful but can read slightly playful if
  enlarged; use them as small scope markers, not hero visuals.

Grayscale:

- pass: most symbols retain shape identity without color;
- pass: trust/provenance symbols remain legible;
- conditional pass: organ-specific symbols need labels when used small because
  liver/kidney/lung silhouettes can be ambiguous at icon scale.

Modules:

- pass: `sample_to_feature_stack`, `claim_boundary_bar`, and
  `space_bio_assay_lane` explain method flow better than generic boxes;
- pass: `trust_status_ribbon` makes provenance/split/caveat logic visible;
- pass: dark-field variants now work on deep canvases and preserve text
  contrast better than transparent light modules;
- conditional pass: dark modules should still be used as restrained method
  components, not enlarged into full-slide diagrams.

Micro-size QA:

- pass: 48 px and 64 px versions are generally stable;
- conditional pass: 32 px versions need nearby text labels, especially organ,
  species, and matrix-style symbols;
- pass: trust symbols such as `source_record`, `checksum_record`,
  `split_guard`, and `caveat_flag` remain readable enough for footer use.

## Claim Boundary

Allowed:

- use symbols as biological signposts;
- use badges to mark status, source, caveat, split, or validation semantics;
- use modules as compact explanations of method, scope, or evidence state;
- use symbols in PPT/Keynote/Illustrator after editing colors or labels.

Avoid:

- using any symbol as proof of microscopy, tissue morphology, spatial data, or
  measured expression;
- enlarging species symbols into hero graphics;
- using `source`, `validated`, or `held-out` badges before the underlying
  source/result status is actually true;
- placing many badges on a slide just to make it look more technical.
- relying on 32 px icons without labels for organ/species categories.

## Consulting-Quality Assessment

Verdict:

- pass as a v0.3 reusable asset layer;
- stronger than v0.1 because it includes scope, status, and method modules;
- stronger after hardening because it now supports light and dark deck surfaces
  and includes micro-size QA;
- suitable for controlled use in consulting-style decks, especially for method
  explanation and claim-boundary signaling;
- not yet a final source-proof layout kit.

Most useful near-term applications:

- methods bridge slides;
- source-to-feature explanations;
- data collection and modality summaries;
- species/organoid expansion maps;
- figure footers and claim-boundary labels;
- executive decks where the audience needs to know what is context, proof,
  processed result, or validation.

## Remaining Gaps for v0.4

Priority upgrades:

1. Native PPTX export or SVG-to-PPT insertion test.
2. 24 px micro-icon QA for extremely small footer/status use.
3. Colorblind simulation beyond grayscale.
4. Source-proof placeholder modules that reserve space for real microscopy,
   dataset screenshots, or result panels.
5. Additional restrained species/model variants when a deck needs more formal
   manuscript tone.

Recommended next step:

- test the dark/light modules inside the actual SpaceBio-Bench methods bridge,
  not only as standalone contact sheets.
