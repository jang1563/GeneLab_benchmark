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
| Mouse tissue NES vs Inspiration4 skin post-flight NES does not support a strict tissue-matching hypothesis on the locally auditable (Hallmark-dominated, n=25 shared pathways) overlap. Mouse skin r=-0.665 (perm p=0.001), mouse gastrocnemius r=+0.826 (the only positive), mouse kidney/eye/liver r in -0.87 to -0.91. | `human_skin_transfer.csv`, `human_skin_transfer.json`, `human_skin_transfer.draft_manifest.json` | `v8/bridge/human_skin_transfer.py` | exploratory; awaiting HPC c2cgp/c5bp re-run to lift n |

## Promotion Requirements

- `SPACEOMICS_ROOT` must resolve to the SpaceOmicsBench checkout used for the run.
- Manifest must record the SpaceOmicsBench commit or release tag.
- Supervised runs must record bootstrap seed, fold strategy, feature list, and model hyperparameters.
- Claims should be written as pathway-level transfer, not raw gene-level prediction.
