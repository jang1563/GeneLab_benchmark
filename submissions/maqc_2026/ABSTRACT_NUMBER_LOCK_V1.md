# MAQC 2026 Abstract Number Lock v1

Status: draft audit ledger; numbers must be rechecked before submission  
Companion abstract: `ABSTRACT_DRAFT_V2.md`

| Item | Reported value | Interpretation | Primary local source | Boundary |
|---|---:|---|---|---|
| Profiles | 549 | Total mouse bulk RNA-seq profiles on the six-tissue audited SpaceBio-Bench surface | task inputs and model-result `n_test` records | scGPT has results for 533 profiles because one 16-sample fold is absent |
| Task folds | 22 | Six-tissue leave-one-mission-out task inventory | `evaluation/geneformer_mouse_gf_all_tissues_summary.json` and tissue result files | This is task inventory, not equal model coverage |
| scGPT coverage | 21 folds | Stored scGPT result coverage | `evaluation/scgpt/scgpt_whole_human_all_tissues_summary.json` | One eye fold is absent |
| Mouse-Geneformer coverage | 22 folds | Stored Geneformer result coverage | `evaluation/geneformer_mouse_gf_all_tissues_summary.json` and `evaluation/geneformer_mouse_gf_A*_lomo_results.json` | Best-test-epoch reporting in original summaries |
| scGPT initial macro AUROC | 0.666 | Mean of six tissue-level summaries after selecting each fold's best held-out test epoch | `evaluation/scgpt/scgpt_whole_human_all_tissues_summary.json` | Historical initial surface |
| scGPT epoch-10 macro AUROC | 0.599 | Mean of six tissue-level fixed epoch-10 fold means | stored `epoch_aurocs` arrays in the scGPT summary | Post hoc sensitivity, not corrected nested validation |
| Mouse-Geneformer initial macro AUROC | 0.476 | Mean of six tissue-level summaries after selecting each fold's best held-out test epoch | `evaluation/geneformer_mouse_gf_all_tissues_summary.json` | Historical initial surface |
| Mouse-Geneformer epoch-10 macro AUROC | 0.458 | Mean of six tissue-level fixed epoch-10 fold means | `history[-1].test_auroc` in `evaluation/geneformer_mouse_gf_A*_lomo_results.json` | Post hoc sensitivity, not corrected nested validation |
| Fixed PCA-LR macro AUROC | 0.730 | Mean of the six tissue-level PCA-LR mean AUROCs (`0.729715210092` before rounding) | `.pca_lr.mean_auroc` in `evaluation/A1_baseline_results.json` through `evaluation/A6_baseline_results.json` | Do not substitute the historical mixed best-classical value of 0.758 |
| Thymus PCA-LR mean mission AUROC | 0.923 | Unweighted mean of four fold-level AUROCs | `evaluation/submission_PCA-LR_baseline_v1_A4_eval.json` | Mission-level estimand |
| Thymus PCA-LR pooled OOF AUROC | 0.631 | AUROC after pooling predictions across four held-out missions | `evaluation/submission_PCA-LR_baseline_v1_A4_eval.json` | Different estimand; calibration across fold-specific models is not guaranteed |

## Required Pre-submission Checks

- Save a machine-readable recomputation artifact for the `0.730` fixed PCA-LR macro value.
- Recompute all four macro values from source files with one script and save the output.
- Confirm the 549-profile denominator and document the absent 16-sample scGPT eye fold.
- Confirm that every abstract number is rounded from the locked artifact, not copied from prose.
- Keep the full scheduled task or profile count as the denominator when describing coverage.
- Do not introduce v4 best-row, SpaceOmicsBench, pathway, countermeasure, or clinical results into this abstract.
