# MAQC 2026 Abstract Draft v2

Status: working draft, not submitted  
Deadline: August 20, 2026  
Format: fewer than 400 words; current presentation availability is poster or lightning talk  
Proposed authors: JangKeun Kim; Christopher E. Mason, pending author review and approval  
Proposed affiliation: Weill Cornell Medicine, New York, NY, USA

## Recommended Title

**From Reproducible Pipelines to Reproducible Claims: An Audit of Model Selection and Aggregation in Spaceflight Omics**

## Abstract

Technical reproducibility is necessary but not sufficient for a reproducible scientific claim. Spaceflight omics is a demanding test case because small cohorts and mission-specific hardware, protocols, age, tissue handling, and processing are entangled with biological variation. We examined how model selection and score aggregation affect conclusions about model performance under cross-study shift.

SpaceBio-Bench contains 549 mouse bulk RNA-seq profiles from six tissues organized into 22 leave-one-mission-out task folds. We audited a fixed PCA-logistic regression pipeline and two widely used first-generation single-cell-pretrained models, scGPT and Mouse-Geneformer, selected because their fold-level training histories were available across the audited benchmark surface. This was a retrospective audit rather than a comprehensive comparison of current models. It separated two choices that had been combined in the initial summary: selecting an epoch using the held-out mission and defining performance as either the mean of mission-level AUROCs or a pooled out-of-fold AUROC. Result coverage was retained explicitly: 21 folds for scGPT and 22 for Mouse-Geneformer.

The initial summary selected each foundation model's best held-out test epoch. In a fixed epoch-10 sensitivity analysis, the six-tissue macro AUROC decreased from 0.666 to 0.599 for scGPT and from 0.476 to 0.458 for Mouse-Geneformer; a fixed PCA-logistic regression pipeline achieved 0.730. Within this audited benchmark surface, neither tested model outperformed the conventional pipeline, and the estimated gaps increased under the fixed-epoch analysis. Aggregation also changed the reported effect. For thymus, the PCA-logistic regression mean mission-level AUROC was 0.923, whereas its pooled out-of-fold AUROC was 0.631. These values answer different questions, and neither is interpretable without naming the evaluation unit and aggregation rule.

Mission-held-out evaluation alone therefore did not make the resulting model-comparison claim fully specified. We propose a minimum reporting profile that records task and fold coverage, train-side model-selection rules, per-mission and pooled estimands, model and preprocessing versions, and prediction-level run records. Spaceflight omics provides a public stress test for moving from rerunnable pipelines to traceable, reproducible claims about AI-enabled science.

## Submission Boundaries

- This is a retrospective evaluation audit, not a new foundation-model leaderboard.
- scGPT and Mouse-Geneformer are included as widely used first-generation reference models with near-complete fold-level training histories, not as a comprehensive 2026 model panel.
- `0.599` and `0.458` are post hoc fixed epoch-10 sensitivity estimates. They are not substitutes for a future train-side validation or nested mission-selection analysis.
- `0.730` is the six-tissue macro AUROC from one fixed PCA-LR workflow on the audited task surface.
- `0.923` and `0.631` are different thymus estimands. The pooled estimate is not presented as uniquely correct because fold-specific models may not share a calibrated probability scale.
- The qualitative comparison is restricted to the tested pipeline configurations. Do not generalize it to all foundation models or to space biology as a domain.
- Mission labels may encode hardware, protocol, age, tissue handling, and processing differences. Classification performance is not evidence for a pure microgravity mechanism.
- Confirm the author list and permission to name Christopher E. Mason before submission.

## Manuscript Expansion After MAQC

The conference abstract is the first presentation of the audit, not the final publication. A subsequent methods manuscript should add:

1. train-side or nested mission model selection rather than test-epoch selection;
2. a fully matched task and result surface across conventional and foundation-model pipelines;
3. a version-locked current-model panel, including a modality-aligned bulk model where species and feature mapping can be made defensible;
4. reconciliation of canonical summaries with raw run records;
5. uncertainty for differences between evaluation surfaces;
6. an independent rerun or external analyst replication;
7. a public evaluation profile, claim register, and versioned prediction records.
