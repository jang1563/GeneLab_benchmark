# Biological Motif v0.3 Deep Research

Date: 2026-06-02

Purpose:

- review why the current biological motif pack still feels underpowered;
- compare it against high-profile journal figure and graphical abstract
  expectations;
- define a v0.3 upgrade path before building more slides.

## Local Diagnosis

Current assets:

- `assets/biovis_motif_pack_v0_2/`
- `output/premium_feature_layer_bridge/`

v0.2 is a useful improvement over v0.1, but it is still mostly a symbolic motif
pack. The assets read as polished icons or soft diagrams, not as biological
evidence surfaces.

Observed weaknesses:

- too many motifs use a similar oval field, thin line, and soft glow grammar;
- tissue, organoid, organ, and molecular motifs differ by label more than by
  biological morphology;
- no scale bars, channel legends, source tags, or microscopy-like visual
  metadata;
- no species grammar for human, mouse, rat, fly, worm, fish, plant, yeast, or
  microbe;
- no modality grammar for bulk RNA-seq, scRNA-seq, spatial transcriptomics,
  proteomics, histology, immunofluorescence, phenotyping, or bioimaging;
- premium depth exists, but the z-stack is atmospheric rather than
  evidence-driven;
- motifs can support a bridge slide, but they cannot yet carry a high-profile
  figure opener.

## External Design Evidence

### Nature Figure Guide

Source:

- https://research-figure-guide.nature.com/
- https://research-figure-guide.nature.com/figures/building-and-exporting-figure-panels/

Relevant lessons:

- Nature treats figures as core research communication, not decoration.
- Main figures should preserve editable vector layers where possible.
- Labels should use standard fonts and remain legible after reduction.
- Accessibility is part of figure preparation.
- Decorative icons should be removed unless they add real context.
- Text labels are often clearer than ambiguous icons, especially for species.

Implication for v0.3:

- motifs must be useful scientific context, not decorative fillers;
- every motif needs an editable vector layer plus a rendered preview;
- species should be labeled explicitly, not only represented by silhouettes.

### Cell Press Graphical Abstract Guidelines

Source:

- https://crosstalk.cell.com/hubfs/Files/GA_guide.pdf

Relevant lessons:

- a single-panel scientific visual should give immediate understanding of the
  take-home message;
- it should read left-to-right or top-to-bottom;
- biological context should be visually indicated, including subcellular
  location, tissue or cell type, and species;
- text should be sparse and labels simple;
- heavily saturated primary colors are distracting;
- one process or one point should be made clear.

Implication for v0.3:

- every motif must have a specific communication job;
- saturated colors should be reserved for focus, not used as the base mood;
- context markers should encode cell/tissue/species/modality without clutter.

### Journal of Cell Biology Graphical Abstract Guidelines

Source:

- https://rupress.org/jcb/pages/graphical-abstract

Relevant lessons:

- graphical abstracts should be simple, informative, self-explanatory, and
  serious;
- they should point to the most important take-home finding;
- they should avoid excess detail, speculative conclusions, and data;
- they allow very little text and require high-resolution export.

Implication for v0.3:

- v0.3 should separate "explanation motifs" from "evidence panels";
- speculative or future-work motifs must visibly carry draft/future status;
- data-like elements must not imply unperformed analyses.

### Image-Based Figure Quality

Source:

- Jambor et al., "Creating clear and informative image-based figures for
  scientific publications", PLOS Biology 2021:
  https://pmc.ncbi.nlm.nih.gov/articles/PMC8041175/

Relevant lessons:

- common failures include missing scale bars, unclear insets, color choices
  inaccessible to colorblind readers, and insufficient explanations of species,
  tissue, labels, annotations, or colors;
- image-based figures need self-contained interpretation aids.

Implication for v0.3:

- source-like biological texture plates need scale bars or scale disclaimers;
- color channels need legends;
- tissue/species/object identity should be explicit;
- inset/focus rings should be standardized.

### Color Accessibility

Source:

- Bang Wong, "Points of view: Color blindness", Nature Methods 2011:
  https://www.nature.com/articles/nmeth.1618

Relevant lessons:

- scientific colors must be distinguishable for color-vision-deficient readers;
- color should not be the only encoding.

Implication for v0.3:

- every motif should pass grayscale and colorblind-simulation checks;
- line style, texture, pattern, and label should back up color.

## External Asset Landscape

### Open Biological Icons

Sources:

- BioIcons: https://bioicons.com/
- BioIcons GitHub: https://github.com/duerrsimon/bioicons
- NIH BioArt Source announcement:
  https://pubweb-prod.niaid.nih.gov/news-events/nih-bioart-source-free-biomedical-art-resource-library

Interpretation:

- these are useful references for semantic coverage and biological object
  vocabulary;
- they should not define the premium style by themselves;
- direct reuse needs license/citation tracking and style harmonization.

### Open Biological Image / Texture Sources

Sources:

- Allen Cell Explorer data download:
  https://www.allencell.org/data-downloading.html
- Human Protein Atlas license:
  https://www.proteinatlas.org/about/licence
- Image Data Resource:
  https://www.nature.com/articles/nmeth.4326
- EMPIAR:
  https://www.ebi.ac.uk/empiar/

Interpretation:

- these are better references for source-derived texture, microscopy-like
  material, and biological morphology;
- license and attribution differ by source and must be tracked per asset;
- raw source use should be kept separate from original schematic artwork.

### Space Biology Scope

Sources:

- NASA GeneLab overview:
  https://www.nasa.gov/centers-and-facilities/ames/what-is-nasas-genelab/
- NASA OSDR overview:
  https://science.nasa.gov/biological-physical/data/osdr/
- NASA OSDR data processing overview:
  https://science.nasa.gov/biological-physical/data/osdr/osdr-data-processing/
- OSD-863 human cortical organoid study:
  https://osdr.nasa.gov/bio/repo/data/studies/OSD-863

Interpretation:

- the benchmark visual language should not look mouse-only;
- OSDR includes omics, environmental, phenotypic, behavioral, bioimaging, and
  physiological data;
- GeneLab/OSDR scope includes rodents, fruit flies, plants, microbes, fish,
  human cells, and organoids;
- v0.3 should include species and modality axes as first-class visual grammar.

## v0.3 Design Direction

v0.3 should not be another icon pack.

It should become a biological visual system with three tiers:

1. Proof Texture Plates
   - microscopy-like cell field;
   - tissue-section field;
   - organoid rosette/spheroid field;
   - spatial spot tissue field;
   - omics matrix/heatmap field;
   - single-cell embedding/cloud field;
   - proteomics spectrum/complex field.

2. Editable Semantic Motifs
   - cell, nucleus, mitochondria, membrane, organoid lumen, tissue contour;
   - pathway node-link, protein complex, RNA trace, DNA/gene marker;
   - species/model-system labels: human, mouse, rat, fly, worm, fish, plant,
     yeast, microbe;
   - organ/tissue labels: brain, muscle, bone, heart, kidney, liver, immune,
     gut, plant tissue.

3. Trust / Annotation Overlays
   - scale bar or scale-disclaimer layer;
   - source tag: original schematic, source-derived, generated texture,
     internal draft;
   - channel legend;
   - inset/focus ring;
   - train/test or hidden-mission boundary tag;
   - species/modality pill labels.

## Required Visual Grammar

Each v0.3 asset should declare:

- biological level: molecule, organelle, cell, tissue, organ, organism,
  experiment, dataset;
- modality: RNA-seq, scRNA-seq, spatial, proteomics, microscopy, phenotype,
  behavior, environmental telemetry;
- species/model system;
- evidence status: schematic, source-derived, generated texture, draft/future;
- allowed slide roles: hero texture, bridge scene, result anchor, callout,
  caption-only support;
- prohibited uses.

## Style Changes

Move away from:

- same oval backing field repeated across all motifs;
- generic line-art organ shapes;
- high similarity between tissue, cell, and molecular assets;
- motifs that need labels to be biologically identifiable;
- decorative cards or icon boxes.

Move toward:

- irregular biological morphology;
- restrained, source-like textures;
- spatial density differences;
- scale/context marks;
- limited but meaningful z-depth;
- fewer colors per asset;
- grayscale-readable contours;
- label-adjacent, not label-dependent, recognition.

## v0.3 Asset Proposal

P0 assets:

- `cell_micrograph_texture_field`
- `tissue_section_texture_field`
- `organoid_rosette_texture_field`
- `spatial_spot_section_field`
- `single_cell_embedding_field`
- `omics_matrix_texture_field`
- `pathway_network_summary_field`
- `species_model_strip`

P1 assets:

- `mitochondrial_ultrastructure_field`
- `protein_complex_orbit_field`
- `immune_cell_population_field`
- `muscle_fiber_section_field`
- `bone_remodeling_context`
- `plant_root_or_leaf_section_field`
- `microbe_colony_field`
- `behavior_phenotype_trace_field`

## QA Gate For v0.3

Every new asset sheet should include:

- full-size visual inspection;
- fit-to-screen thumbnail inspection;
- grayscale export;
- colorblind-safe simulation or palette check;
- line-weight and text-size check;
- source/license manifest;
- claim-boundary manifest;
- "decorative icon" check;
- slide-in-context test, not only contact sheet.

## Recommendation

Build v0.3 as a smaller but stricter pack.

Do not expand from 16 motifs to 40 motifs yet. First create 8 P0 assets that
solve the main premium gap:

> source-like biological texture plus editable scientific interpretation.

The immediate next build should produce:

- a v0.3 contact sheet;
- a v0.3 usage sheet;
- a grayscale QA sheet;
- a manifest with biological level, modality, species, evidence status, and
  allowed/prohibited uses;
- one direct replacement test inside the feature-layer bridge scene.

## Working Verdict

v0.2 is acceptable as a restrained context layer.

v0.3 is necessary before more premium deck production because the deck needs a
visual language that feels biological before it feels decorative.
