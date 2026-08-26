# MAQC 2026 Submission-Ready Abstract

Status: submitted through EasyChair on August 19, 2026  
Deadline: August 20, 2026; no deadline time or timezone published  
Available presentation formats: poster or lightning talk  
Body length: 339 words by whitespace count  

## Title

**From Reproducible Pipelines to Reproducible Claims: An Audit of Model Selection and Aggregation in Spaceflight Omics**

## Authors

JangKeun Kim; Christopher E. Mason

## Affiliation

Department of Physiology and Biophysics, Weill Cornell Medicine, New York, NY, USA

## Abstract

An analysis can be rerun exactly even when its scientific conclusion depends on unreported choices. That distinction matters in spaceflight omics, where small cohorts and differences in hardware, protocols, animal age, tissue handling, and processing can be confounded with biological variation. We examined how model selection and score aggregation changed conclusions about model performance under cross-study shift.

SpaceBio-Bench includes 549 mouse bulk RNA-seq profiles from six tissues, evaluated in 22 leave-one-mission-out task folds. We audited a fixed PCA-logistic regression pipeline and two established single-cell-pretrained expression models, scGPT and Mouse-Geneformer. Because the aim was to audit the evaluation process rather than rank current models, we selected models with near-complete fold-level training records. The audit separated two choices that had been combined in the original summary: selecting an epoch using the held-out mission and calculating performance as either the mean of mission-level AUROCs or a pooled out-of-fold AUROC. We report coverage explicitly: 21 folds for scGPT and 22 for Mouse-Geneformer.

The original summary retained the epoch with the highest held-out test score for each pretrained model. In a sensitivity analysis using epoch 10 for every task, the six-tissue macro AUROC decreased from 0.666 to 0.599 for scGPT and from 0.476 to 0.458 for Mouse-Geneformer. The fixed PCA-logistic regression pipeline achieved 0.730. Across these six tissues, neither pretrained model outperformed the conventional pipeline, and the performance gaps widened under the fixed-epoch analysis. The aggregation rule also changed the reported performance substantially. For thymus, the mean mission-level AUROC for PCA-logistic regression was 0.923, whereas the pooled out-of-fold AUROC was 0.631. Because these estimands answer different questions, the evaluation unit and aggregation rule must be reported with the result.

Holding out entire missions was therefore not enough to make the model comparison unambiguous. We propose a minimum reporting profile that records task and fold coverage, model-selection rules and the data used to apply them, per-mission and pooled estimands, model and preprocessing versions, and prediction-level outputs. Together, these records would make it possible to reproduce not only the pipeline, but the reported model comparison itself.

## Suggested Keywords

spaceflight omics; reproducibility; foundation models; mission-held-out evaluation; model selection

## Submission Boundaries

- This is a retrospective evaluation audit, not a current-model leaderboard.
- The epoch-10 analysis is a post hoc fixed-epoch sensitivity analysis, not nested validation.
- The PCA-logistic regression comparison is limited to the same six-tissue audited task surface.
- Mean mission-level and pooled out-of-fold AUROC are distinct estimands; the pooled value is not privileged as uniquely correct.
- Mission labels may encode hardware, protocol, age, tissue handling, and processing differences. Classification performance is not evidence for a pure microgravity mechanism.
- Christopher E. Mason approved the abstract for submission in Slack on August 15, 2026.
