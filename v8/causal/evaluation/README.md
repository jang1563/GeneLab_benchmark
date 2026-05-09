# Causal Result Map

Integration question: which stressor-to-tissue edges are stable enough to become
causal candidates for later counterfactual modeling?

## Claims And Artifacts

| Claim | Primary output | Script | Status |
|---|---|---|---|
| ICP-style stability ranks stressor effects across analog and flight environments. | `icp_stressor_pathway_scores.csv`, `icp_dag_summary.json` | `v8/causal/build_dag.py` | hpc_validated; manifest recorded |
| DAG edges are evidence links, not validated intervention effects. | `dag_edges.csv` | `v8/causal/build_dag.py` | exploratory |

## Promotion Requirements

- Manifest must record all environments used, minimum environment count, edge-selection thresholds, and input result files.
- Manuscript language must separate observational stability from causal intervention claims.
- Counterfactual examples require explicit assumptions and falsification checks before promotion.
