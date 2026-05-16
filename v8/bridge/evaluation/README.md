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
| Cross-tissue generalization of the 11-pathway translation/proteostasis suppression signature (paired with `skin_transfer_driver_decomp.py`): 23 of 25 valid human compartments (all 9 Inspiration4 PBMC subsets, 8 of 10 NASA_Twins CD4/CD8 contrasts, the i4_skin anchor) have mean proteostasis-11 NES below -2.0, matching the mouse gastrocnemius anchor mean of -2.25. A within-overlap control on the 14 non-proteostasis pathways shows a weaker, less uniform suppression pattern, supporting ribosomal-protein-driven translation suppression as a specifically conserved cross-tissue spaceflight signal. | `proteostasis_generalization.csv`, `proteostasis_generalization.json`, `proteostasis_generalization.draft_manifest.json` | `v8/bridge/proteostasis_generalization.py` | exploratory; per-compartment Spearman r noisy at n=4 pathways for Twins compartments, conclusion anchored on the mean-NES distribution |

## Promotion Requirements

- `SPACEOMICS_ROOT` must resolve to the SpaceOmicsBench checkout used for the run.
- Manifest must record the SpaceOmicsBench commit or release tag.
- Supervised runs must record bootstrap seed, fold strategy, feature list, and model hyperparameters.
- Claims should be written as pathway-level transfer, not raw gene-level prediction.
