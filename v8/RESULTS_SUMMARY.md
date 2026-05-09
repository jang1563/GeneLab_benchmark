# SpaceMed v8 Results Summary

## Key Findings Across Three Pillars

### Pillar 1: Species Transfer (BRIDGE)
- RF supervised AUROC (mouse NES → I4 post-flight): **0.888**
- Baseline AUROC (SpaceOmics features only): 0.712
- Mean AUROC improvement: **+0.175 [0.134, 0.219]** (95% CI)

### Pillar 2: Stressor Decomposition
- Top ICP-stability stressor: **T** (ICP=0.5401)
- Interaction dominance: 44–61% variance in top-responsive genes
- Radiation quality discovery: Low-LET γ (r=+0.36) vs high-LET HZE (r=–0.22) opposite signs

### Pillar 3: Countermeasure-Hypothesis Triage
- Multi-tissue perturbation reversals: 2 on Pareto front
- Top multi-tissue signature-reversal hypothesis: **CGP-60474**
- CRISPR orthogonal support: CDK9 KO supports CDK-inhibitor target plausibility

### Integrated Mars Regime Flagging
- Linear extrapolation breaks at >5× dose amplification
- Top Mars-sensitive genes in the linear stress test: WNT10B (~1052×), KRTAP19-2 (~414×), ENSMUSG00000092534 (~182×)
- Bootstrap CIs propagated from factorial β ± SE; outputs are exploratory regime flags, not point predictions

---

## Detailed Tissue-Level Results

```
              species_transfer_NES_r species_transfer_n_compartments interaction_variance_frac_spleen n_sig_genes_spleen interaction_variance_frac_hze_endocrine n_sig_genes_hze_endocrine radiation_quality_gamma_r        l1000cds2_top_reverser l1000cds2_top_score crispr_supported_drugs mars_vs_flight_spearman_r mars_n_genes interaction_variance_frac_skin_analog n_sig_genes_skin_analog interaction_variance_frac_brain n_sig_genes_brain crispr_supported_drugs_list
thymus                        -0.449                               9                           61.79%               4556                                  44.00%                      1347                    +0.162                   ALBENDAZOLE               0.077                      0                   +0.0593        28111                                   NaN                     NaN                             NaN               NaN                         NaN
gastrocnemius                  0.469                               9                              NaN                NaN                                     NaN                       NaN                       NaN  Dorsomorphin dihydrochloride               0.076                      0                       NaN          NaN                                   NaN                     NaN                             NaN               NaN                         NaN
skin                          -0.384                               9                              NaN                NaN                                     NaN                       NaN                    -0.128                       MLN2238               0.114                      0                   -0.1263        28262                                39.92%                    7791                             NaN               NaN                         NaN
eye                           -0.642                               9                              NaN                NaN                                     NaN                       NaN                    +0.004                   penfluridol               0.059                      0                   -0.0164        27801                                   NaN                     NaN                          21.83%             21249                         NaN
liver                         -0.310                               9                              NaN                NaN                                     NaN                       NaN                       NaN  Dorsomorphin dihydrochloride               0.060                      0                       NaN          NaN                                   NaN                     NaN                             NaN               NaN                         NaN
kidney                        -0.204                               9                              NaN                NaN                                     NaN                       NaN                       NaN                        EI-156               0.091                      2                       NaN          NaN                                   NaN                     NaN                             NaN               NaN           AZD-5438, AT-7519
```