# INTERVENE Result Map

Pillar question: which perturbations reverse tissue signatures enough to justify
pre-clinical validation?

## Claims And Artifacts

| Claim | Primary output | Script | Status |
|---|---|---|---|
| Tissue-level flight signatures can be exported as LINCS-ready up/down gene sets. | `signatures/signatures_manifest.json`, `signatures/*_{up,dn}.txt`, `signatures/*_ranked.csv` | `v8/intervene/export_signatures.py` | hpc_validated; manifest recorded |
| L1000CDS2 returns tissue-specific chemical reversal hypotheses. | `lincs_summary.json`, `lincs_*_top50.csv` | `v8/intervene/lincs_query.py` | hpc_validated snapshot; API payload archive still required before beta freeze |
| Multi-tissue Pareto filtering prioritizes candidates with broader signature reversal. | `pareto_front.csv`, `multi_tissue_drug_matrix.csv`, `multi_tissue_drug_scores.json` | `v8/intervene/pareto_multi_tissue.py` | hpc_validated snapshot; exploratory |
| Enrichr CRISPR KO signatures provide orthogonal target plausibility checks. | `crispr_orthog_summary.json`, `crispr_orthog_*_{up,dn}.csv` | `v8/intervene/perturb_seq_orthog.py` | hpc_validated snapshot; API payload archive still required before beta freeze |
| Offline DGIdb reversal scoring is a lightweight internal sanity check only. | `offline_reversal_summary.json`, `offline_reversal_*.csv` | `v8/intervene/offline_reversal_scorer.py` | exploratory |

## Promotion Requirements

- Manifest must record query timestamp, query payload hash, API endpoint, `db-version`, top-N gene list size, and gene-symbol normalization.
- Drug language must remain hypothesis-generating until independent safety and efficacy validation.
- Raw API dumps should stay outside Git unless compact and explicitly needed for reproducibility.
