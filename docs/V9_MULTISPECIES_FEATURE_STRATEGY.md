# v9 Multispecies Feature Strategy

Status: initial strategy
Date: 2026-05-21
Depends on: `docs/V9_ORGANOID_AND_SPECIES_EXTENSION_REVIEW_2026_05_21.md`

## Decision

Do not build a single raw-gene leaderboard across mouse, human, fly, plant, and
microbial data.

v9 multispecies work should use two layers:

1. Species-native tasks, where each species keeps its own gene identifiers and
   only within-species splits are evaluated.
2. Cross-species bridge tasks, where the shared feature space is pathway,
   signature, or ortholog-group based.

The first cross-species bridge should use pathway/NES-style features because
this repository already contains historical pathway and NES transfer outputs,
and because pathway-level features reduce avoidable gene-id and orthology
friction.

## Feature namespaces

| Namespace | Use | Allowed task family | First status |
|---|---|---|---|
| `mouse_gene` | Current mouse bulk LOMO | `bulk_lomo` | Active |
| `human_gene` | Human organoid RNA-seq pilot | `human_organoid_spaceflight` | Draft |
| `species_gene_ids` | Within-species fly or plant tasks | `multispecies_spaceflight` | Draft |
| `pathway_nes` | Cross-species pathway/signature transfer | `bridge_cross_species` | Recommended first bridge |
| `ortholog_group` | Gene-level cross-species comparisons after ortholog mapping | `bridge_cross_species` | Later |
| `curated_stress_signature` | Hand-curated conserved stress axes | `stressor_radiation_quality`, `bridge_cross_species` | Later |

## Species-native task rule

Species-native tasks may use raw genes only when train and test samples are from
the same species and the same feature namespace.

Examples:

- OSD-207 Drosophila flight-versus-ground, with genotype blocking.
- OSD-37 Arabidopsis seedling-pool flight-versus-ground, with ecotype blocking.
- OSD-120 Arabidopsis root response, with genotype and light-condition blocking.

These tasks can share evaluation code, but they should not share a leaderboard
row with mouse raw-gene LOMO unless the task family, organism, and feature
namespace are explicit.

## Cross-species bridge rule

Cross-species tasks must declare a non-raw-gene feature strategy before any
baseline is run.

Recommended first bridge:

- Convert species-native differential signatures to pathway or NES-like
  summaries.
- Evaluate rank correlation, direction concordance, and pathway-level
  classification/transfer metrics.
- Keep animal, plant, and human cell-model branches separately reported until
  there is enough evidence that one bridge metric is fair across them.

Later bridge:

- Add ortholog groups when the mapping, dropped-feature rules, one-to-many
  handling, and species coverage are explicit.
- Compare ortholog-group results against the pathway/NES bridge as a sensitivity
  analysis, not as the default first result.

## Manifest requirements

Any multispecies or bridge manifest must include:

- `organism`
- `taxon_id`
- `species_common_name`
- `material_type`
- `model_system`
- `biospecimen_type`
- `assay_modality`
- `feature_namespace`
- `orthology_strategy`
- `donor_or_strain_block`
- split-blocking fields relevant to genotype, ecotype, light condition, donor,
  organoid type, or microglia condition

For cross-species bridge manifests, also include:

- pathway database or signature collection
- gene-id namespace before mapping
- mapping version or source
- one-to-many ortholog handling policy, if applicable
- dropped-feature count and reason

## Initial implementation order

1. Keep current mouse bulk LOMO as `mouse_gene`.
2. Keep OSD-863/OSD-871 human organoid RNA-seq as `human_gene`.
3. Add OSD-207/OSD-37/OSD-120 species-native draft task cards using
   `species_gene_ids`.
4. Draft a `pathway_nes` bridge task after species-native source integrity is
   checked.
5. Add `ortholog_group` only after pathway/NES baselines establish a comparison
   point.

## Claim boundary

Allowed early claim:

> v9 now separates mouse, human organoid, and non-mouse species sources by
> organism, model system, modality, and feature namespace, enabling future
> multispecies task cards without collapsing incompatible raw feature spaces.

Not allowed yet:

- A cross-species model ranking.
- A claim that fly, plant, mouse, and human organoid responses are directly
  comparable at raw gene level.
- A claim that pathway/NES bridge metrics are biologically universal.

