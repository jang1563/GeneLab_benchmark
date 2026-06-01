# Biological Visual Symbol Asset Pack

Date: 2026-06-02

Purpose:

- collect reusable biological symbols for premium SpaceBio-Bench slides;
- add cell, molecular, tissue, organ, and assay-specific visual vocabulary;
- avoid generic workflow diagrams by giving method/result slides biological
  signposts that match the project visual identity.

## Local Asset Pack

Output root:

- `assets/biovis_symbol_pack_v0_1/`

Generated files:

- `assets/biovis_symbol_pack_v0_1/svg/`
- `assets/biovis_symbol_pack_v0_1/png/`
- `assets/biovis_symbol_pack_v0_1/manifest.json`
- `assets/biovis_symbol_pack_v0_1/manifest.csv`
- `assets/biovis_symbol_pack_v0_1/biovis_symbol_pack_contact_sheet.svg`
- `assets/biovis_symbol_pack_v0_1/biovis_symbol_pack_contact_sheet.png`

Generator:

- `scripts/build_biovis_symbol_asset_pack.py`

Implementation note:

- SVG is the source of truth.
- PNG previews are generated from SVG using `rsvg-convert`.
- These are original project assets and follow the repository MIT license.

## Included Symbols

| Category | Assets | Best Use |
|---|---|---|
| Cellular | `cell_nucleus`, `mitochondrion`, `organoid_sphere` | sample/cell state, organelle biology, human organoid extension |
| Tissue | `tissue_section`, `epithelium_layer`, `vascular_slice` | tissue context, barrier/lining tissue, vascular/perfusion hints |
| Molecular | `dna_helix`, `rna_strand`, `protein_complex`, `pathway_nodes` | gene/RNA/protein/pathway feature explanation |
| Organ | `brain_outline`, `lung_outline`, `kidney_outline`, `liver_lobule` | tissue-specific result slides and organ-level navigation |
| Assay | `sample_tube`, `expression_matrix` | sample collection, expression matrix/task input |

## Design Verdict

Pass for v0.1 asset library.

Strengths:

- style matches the current visual identity tokens;
- SVG assets are editable and scale cleanly;
- biological specificity is stronger than generic circles, arrows, and boxes;
- the set covers the immediate deck needs: cell, RNA, pathway, tissue, organ,
  sample, and matrix.

Use with restraint:

- these are scientific signposts, not source-derived proof images;
- do not use a symbol to imply an assay, tissue, or modality that the data do
  not support;
- do not enlarge these into hero artwork unless they are paired with a real
  evidence surface or generated/curated biological texture.

## Recommended Slide Uses

B1-B4 bridge slides:

- B1: use `cell_nucleus`, `tissue_section`, or `expression_matrix` as faint
  source-field texture only if the source field feels too abstract.
- B2: use `sample_tube`, `rna_strand`, and `expression_matrix` as small
  context markers around the task record.
- B3: avoid adding biology symbols unless they clarify a specific mission or
  tissue; the mission-held-out boundary should remain the visual thesis.
- B4: use `expression_matrix` sparingly to identify the feature matrix before
  train-only processing.

Feature-layer bridge:

- primary symbols: `rna_strand`, `dna_helix`, `pathway_nodes`,
  `protein_complex`;
- likely visual move: expression matrix or RNA signal compresses into pathway
  nodes, with a caveat that pathway summaries are biological features, not
  universally better features.

Organoid extension:

- primary symbols: `organoid_sphere`, `brain_outline`, `rna_strand`;
- caveat required: organoid extension is a bounded biology check or draft
  extension, not a completed validation benchmark unless the underlying task
  and payload are ready.

Tissue/result slides:

- primary symbols: `kidney_outline`, `liver_lobule`, `lung_outline`,
  `brain_outline`;
- use as small panel anchors beside actual result plots, not as substitutes for
  measured data.

## External Source Reserve

Use these only when the local original pack is insufficient or when a slide
needs a recognizable biological object that would be costly to redraw.

BioIcons:

- useful for broad biology/chemistry SVG icons;
- BioIcons is a free open-source SVG icon library, but individual icons can
  carry icon-specific license and citation metadata;
- deck rule: import only icons with explicit per-icon license metadata, record
  attribution in the slide manifest, and restyle to SpaceBio-Bench tokens.

Servier Medical Art:

- useful for anatomical and medical illustrations;
- Servier states the library has more than 3,000 free professional medical
  illustrations under CC BY 4.0;
- deck rule: safe as attributed medical/anatomy illustration reserve, but do
  not let Servier's stock style dominate the brand.

SciDraw:

- useful for scientific SVG drawings with DOI-backed records;
- SciDraw states content is shared under Creative Commons CC-BY unless stated
  otherwise;
- deck rule: use only with creator/source/DOI attribution preserved in the
  manifest.

Allen Cell Explorer:

- useful for real cell-image proof objects, not generic decorative assets;
- citation and terms language should be checked per use because some content is
  framed for noncommercial research and citation compliance;
- deck rule: use only as source-derived evidence with citation, not as generic
  background texture.

## Next Additions

The v0.1 pack is missing several likely useful biology/space-biology symbols:

- species: mouse, rat, fly, fish, Arabidopsis/plant;
- assay/modality: bulk RNA-seq, scRNA-seq droplet, spatial spot grid,
  proteomics/mass-spec;
- space-biology context: mission patch placeholder, orbit arc, microgravity
  cue, radiation/stressor cue;
- data-trust context: source record, checksum, task manifest, train/test split
  guard.

Recommended v0.2:

1. Add species/modality symbols for organoid/multispecies/single-cell slides.
2. Create a flatter manuscript-safe variant with black/gray strokes only.
3. Build a slide test sheet that places symbols on B1-B4 and the feature-layer
   bridge to check whether they improve comprehension or merely decorate.
