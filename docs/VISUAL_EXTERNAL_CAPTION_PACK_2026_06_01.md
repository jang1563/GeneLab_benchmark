# External-Facing Caption Pack

Date: 2026-06-01

Purpose: provide manuscript-style captions for the current Fig1/Fig2/Fig3/Fig6
set without internal operating language. Source paths and implementation notes
belong in provenance manifests, not inside the figures.

## Figure 1. Tissue-Specific Cross-Mission Generalization

Held-out mission performance varies strongly by tissue in mouse bulk RNA-seq.
Points show mean AUROC for one-mission-held-out evaluation; marker shape marks
high, mid-range, and near-chance transfer groups. Horizontal intervals show
bootstrap uncertainty, and the dashed vertical line marks chance performance.
Thymus and gastrocnemius show the strongest cross-mission generalization,
whereas liver and kidney are close to chance.

## Figure 2. Pathway Summaries Reduce Recoverable Experimental Labels

Pathway-level summaries reduce recoverable mission, hardware, and gravity label
signals relative to gene-level inputs while improving selected tissue-level
flight-detection tasks. Panel A shows macro-F1 for coupled-label prediction
using gene-level and pathway-level inputs. Panel B shows the pathway-minus-gene
change in held-out mission AUROC for flight detection. Panel C compares
tissue-level pathway activity agreement with held-out mission AUROC as a
descriptive association across six tissues; no fitted trend is used.

## Figure 3. Classical Baselines Remain Competitive on the Shared Task Set

Model scale alone does not resolve the small-sample held-out mission prediction
problem. Panel A summarizes mean AUROC or task score across model families,
with the shared six-task rows providing the direct comparison among PCA-LR,
scGPT, and Mouse-Geneformer. Panel B shows tissue-specific AUROC changes for
scGPT and Mouse-Geneformer relative to matched classical baselines. Task-scope
definitions for text-only rows and single-tissue extension rows are reported
in the accompanying model-comparison table.

## Figure 6. Human Neural Organoids as a Small Biology-Check Dataset

Public human neural-organoid RNA-seq data provide a small extension dataset for
checking whether model signals align with flight-response biology. Panel A
summarizes the available public sources, samples, and gene-level contrast
footprint. Panel B separates primary flight/ground prediction from secondary
biology checks. Panel C shows enrichment of significant flight-response genes
among top-ranked model genes. Interpretation should account for the small
sample count and coupled source, disease, donor, and microglia factors.

## Figure Text Policy

Use inside figures:

- scientific readouts, axis labels, panel titles, sample counts, and direct
  model/tissue names;
- concise method qualifiers such as "held-out mission" or "descriptive
  association" when they are needed to read the plot.

Keep outside figures:

- source-file paths;
- implementation status labels;
- package/release decisions;
- internal wording about allowed/prohibited claims;
- review instructions about how a figure should or should not be used.
