# MAQC 2026 Abstract Draft v1

Status: working draft, not submitted  
Deadline: August 20, 2026  
Format: fewer than 400 words; current presentation availability is poster or lightning talk  
Proposed authors: JangKeun Kim; Christopher E. Mason, pending author review and approval  
Proposed affiliation: Weill Cornell Medicine, New York, NY, USA

## Title

**SpaceBio-Bench: How Evaluation Design Changes Apparent Foundation-Model Gains in Spaceflight Omics**

## Abstract

Foundation models are increasingly applied to small, heterogeneous omics studies, but technical reproducibility does not guarantee a valid comparison under mission or cohort shift. We audited two public spaceflight-omics benchmark surfaces to determine how model selection, score aggregation, and outcome construction affect apparent gains over conventional methods.

SpaceBio-Bench comprises 549 mouse bulk RNA-seq tissue profiles from six tissues and 22 leave-one-mission-out folds. We compared stored results for a fixed PCA-logistic regression pipeline, scGPT, and Mouse-Geneformer, auditing when the held-out mission was evaluated and how performance was aggregated across missions. In SpaceOmicsBench, we traced a cell-free RNA differential-expression prediction task from its 466 gene labels to its 29 input features and source statistical rule.

The initial public SpaceBio summary selected each foundation model's best held-out test epoch. Using a fixed final-epoch report reduced the six-tissue macro AUROC from 0.666 to 0.599 for scGPT and from 0.476 to 0.458 for Mouse-Geneformer. A fixed PCA-logistic regression pipeline achieved 0.730. Thus, the conclusion that the tested foundation models did not outperform the classical baseline remained, but the estimated gap increased. Aggregation also materially changed estimates. For thymus, the mean of mission-level AUROCs was 0.923, whereas pooled out-of-fold AUROC was 0.631. In SpaceOmicsBench, the label-defining effect-size conditions selected 803 genes and contained all 466 positive labels. The resulting TabPFN AUPRC of 0.957 therefore could not be interpreted as an independent test of biological generalization.

In small-sample omics, apparent model gains can be driven by evaluation choices even when code and data are available. We propose a reporting contract that separates train-side model selection, mission-level and pooled estimands, label-construction provenance, and versioned run records. These controls preserve useful negative results while preventing a reproducible pipeline from becoming an unsupported scientific claim.

## Internal Submission Boundaries

- The abstract reports a retrospective benchmark audit, not a new foundation-model leaderboard.
- `0.730` is the six-tissue macro mean from one fixed PCA-LR pipeline on the v1 mission-held-out task surface.
- `0.599` and `0.458` are fixed epoch-10 summaries computed from stored scGPT and Mouse-Geneformer histories. They are audit estimates, not substitutes for future train-side validation-based model selection.
- `0.923` and `0.631` are distinct thymus aggregation estimands. The pooled value is not presented as the uniquely correct estimate because predictions come from different fold-specific models and may not share a calibrated probability scale.
- The SpaceOmicsBench result is a construct-validity audit of one cfRNA task. It does not invalidate every SpaceOmicsBench task or support a broad ranking claim.
- Do not add pathway-rescue, v4 best-row, text-LLM, countermeasure, or clinical claims to this abstract.
- Confirm author list, title, and permission to name Christopher E. Mason before submission.
