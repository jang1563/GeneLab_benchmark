# Biological Symbol Module Application and PPTX Smoke Review

Date: 2026-06-02

Purpose:

- follow the recommended next step after the v0.3 symbol/module pack;
- preserve the reusable asset pack in git;
- test the modules inside real bridge-style layouts;
- run a PPTX/SVG smoke test before using these assets in a production deck;
- define the next v0.4 source-proof placeholder direction.

## Commit

Committed the v0.3 biological symbol/module pack:

- commit: `1f6380a Add biological symbol module asset pack`

This commit preserves:

- 36 biological/trust symbols;
- 8 evidence/status badges;
- 7 light modules;
- 7 dark-field modules;
- micro-size QA;
- grayscale QA;
- manifest and review documentation.

## Bridge Application Test

Generator:

- `scripts/build_biovis_symbol_module_bridge_application.py`

Output root:

- `output/biovis_symbol_module_bridge_application/`

Generated outputs:

- `output/biovis_symbol_module_bridge_application/01_light_bridge_application.png`
- `output/biovis_symbol_module_bridge_application/02_dark_bridge_application.png`
- `output/biovis_symbol_module_bridge_application/03_light_dark_ab_contact_sheet.png`
- `output/biovis_symbol_module_bridge_application/biovis_symbol_module_bridge_application_manifest.json`
- `output/biovis_symbol_module_bridge_application/biovis_symbol_module_bridge_application_qa.json`

## Bridge Verdict

Dark bridge:

- pass;
- strongest near-term direction for premium SpaceBio-Bench bridge slides;
- dark modules preserve contrast and give the canvas more depth;
- source/status/proof layers are easier to separate.

Light bridge:

- conditional pass;
- the sample-to-feature module explains the method better than generic boxes;
- lower half has too much open space for a final slide;
- final production slide should either tighten the lower evidence modules or
  replace open space with real source-proof objects.

Most useful modules after actual application:

- `sample_to_feature_stack`;
- `space_bio_assay_lane`;
- `claim_boundary_bar`;
- `trust_status_ribbon`.

Use with care:

- `species_coverage_strip` is useful as compact scope grammar;
- it should not become a hero visual or main thesis.

## PPTX/SVG Smoke Test

Generator:

- `scripts/build_biovis_symbol_module_pptx_smoke.mjs`

Output root:

- `output/biovis_symbol_module_bridge_application/pptx_smoke/`

Generated outputs:

- `output/biovis_symbol_module_bridge_application/pptx_smoke/biovis_symbol_module_svg_editability_smoke.pptx`
- `output/biovis_symbol_module_bridge_application/pptx_smoke/biovis_symbol_module_pptx_smoke_manifest.json`
- `output/biovis_symbol_module_bridge_application/pptx_smoke/rendered_pdf/biovis_symbol_module_svg_editability_smoke.pdf`
- `output/biovis_symbol_module_bridge_application/pptx_smoke/rendered_pdf/biovis_symbol_module_svg_editability_smoke-1.png`

Smoke result:

- pass: PPTX was generated;
- pass: representative SVG assets are present inside `ppt/media`;
- pass: local LibreOffice opened and exported the PPTX to PDF;
- pass: PDF-to-PNG render shows all representative SVG modules/symbols/badges.

Technical note:

- `@oai/artifact-tool` import was blocked in this local session by a macOS
  native-library code-signing issue;
- this was not treated as a production deck build;
- fallback used `pptxgenjs` only for a smoke test;
- production deck generation should still use the presentation runtime when it
  is available.

## What This Proves

The v0.3 symbol/module pack is no longer only a static asset gallery.

It can now support:

- real bridge-style method explanation;
- dark and light canvas variants;
- compact scope/status markers;
- PPTX embedding with SVG media preserved;
- local office export/render without obvious asset loss.

## Remaining Risk

These tests do not prove full PowerPoint editability at the level of editing
individual vector paths inside PowerPoint.

What passed:

- SVG media are embedded;
- PowerPoint-compatible package structure exists;
- local office tooling can render/export the deck.

What still needs human app-level checking:

- whether PowerPoint allows "Convert to Shape" or direct vector editing for
  each SVG;
- whether Keynote preserves the same SVG appearance;
- whether brand-color retinting is easier via native PPT controls or by editing
  SVG source before insertion.

## v0.4 Source-Proof Placeholder Direction

The next asset-system step should not be more icons.

The next step should be source-proof placeholder modules: reusable frames that
reserve space for real evidence objects while keeping claim boundaries visible.

Recommended v0.4 modules:

1. `source_dataset_record_plate`
   - source screenshot, dataset accession, sample count, source badge;
   - use for OSDR/GeneLab/GEO record proof.

2. `expression_matrix_proof_plate`
   - real matrix heatmap or QC snapshot;
   - source/proceeded badges and checksum line.

3. `single_cell_embedding_proof_plate`
   - real embedding or placeholder for real embedding;
   - include species, tissue, cell-count, and caveat controls.

4. `organoid_evidence_plate`
   - organoid source image or source-derived assay panel;
   - explicit extension/preliminary badge.

5. `held_out_task_proof_plate`
   - train/test or mission-held-out boundary;
   - use with split guard and held-out badge.

6. `result_claim_plate`
   - result plot plus claim-boundary footer;
   - separates processed result from validated claim.

Design rule for v0.4:

> The module should reserve the proof object, status badge, source line, and
> caveat position before the final figure is available.

This prevents the slide from being designed around decorative biology assets
and then trying to retrofit evidence later.

## Recommended Next Execution Step

Build a small v0.4 source-proof placeholder pack with 4 to 6 modules, then test
those modules inside the same bridge application output.

Priority order:

1. source dataset record plate;
2. expression matrix proof plate;
3. single-cell embedding proof plate;
4. held-out task proof plate;
5. organoid evidence plate.
