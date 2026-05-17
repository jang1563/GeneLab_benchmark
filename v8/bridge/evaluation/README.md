# BRIDGE Result Map

Pillar question: which mouse pathway signals transfer to human spaceflight data?

## Claims And Artifacts

| Claim | Primary output | Script | Status |
|---|---|---|---|
| Aggregate mouse pathway NES correlates with human I4/Twins compartments. | `species_transfer_nes.json`, `species_transfer_nes.csv` | `v8/bridge/link_spaceomicsbench.py` | hpc_validated; manifest recorded |
| Per-tissue mouse NES identifies the tissue-specific direction and strength of human transfer. | `tissue_nes_spearman.json`, `tissue_nes_spearman.csv` | `v8/bridge/tissue_nes_bridge.py` | hpc_validated; manifest recorded |
| Mouse tissue NES features improve the I4/Twins conserved-pathway classifier over SpaceOmics-only features. | `supervised_conservation.json` | `v8/bridge/supervised_conservation.py` | hpc_validated; manifest recorded |
| Supervised BRIDGE feature construction excludes the target label from model features and records deterministic stratified folds. | `bridge_leakage_audit.json` | `v8/bridge/leakage_audit.py` | hpc_validated audit; upstream SpaceOmicsBench feature-builder freeze still required before beta |
| Human analog isolation betas are exploratory external-data checks, not primary v8 claims. | `human_analog_betas_*.csv`, `rodent_to_human_analog_transfer.json` | `v8/bridge/human_analog_isolation.py` | exploratory |
| Per-mission decomposition of the pooled mouse-tissue proteostasis-11 signature (paired with PRs #7-#10). Only gastrocnemius shows both missions consistently translation-suppressed (RR-1 -2.04 / RR-9 -2.47, range 0.42) with significant positive Spearman r against the human i4_skin anchor (p = 0.001 and 0.008). Skin shows striking within-tissue heterogeneity (range 5.14): MHU-2_dorsal alone (-2.23) reproduces the gastrocnemius-like suppression while MHU-2_femoral (+1.85) and RR-6 (+2.91) are oppositely directed. Companion supplementary figure visualizes per-(tissue, mission) mean NES + sign of Spearman r against i4_skin. | `proteostasis_mission_stability.csv`, `proteostasis_mission_stability.json`, `proteostasis_mission_stability.draft_manifest.json`, `v8/figures/SupplementaryFigure_Proteostasis_Mission_Stability.{png,pdf}` | `v8/bridge/proteostasis_mission_stability.py` + `v8/bridge/proteostasis_mission_stability_figure.py` | exploratory; per-mission n=11 pathways so test is not power-limited at n |

## Promotion Requirements

- `SPACEOMICS_ROOT` must resolve to the SpaceOmicsBench checkout used for the run.
- Manifest must record the SpaceOmicsBench commit or release tag.
- Supervised runs must record bootstrap seed, fold strategy, feature list, and model hyperparameters.
- Claims should be written as pathway-level transfer, not raw gene-level prediction.
