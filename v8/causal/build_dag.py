"""v8 Causal DAG — Invariant Causal Prediction across spaceflight environments.

Uses stressor environments as ICP "settings": each (analog, tissue) pair provides
a distinct interventional environment. Identifies causal parents of each tissue's
pathway NES by checking which stressor betas remain stable (low variance across
settings) — the ICP stability criterion.

Environments defined:
  E1: ISS-like (µg + low-LET GCR)     → flight Wald-z
  E2: HLU-only (ground, OSD-211/237)  → beta_HLU from factorial
  E3: IR-only (ground, 0.04 Gy 57Co)  → beta_IR from factorial
  E4: HLU+IR interaction              → beta_HLUxIR from factorial
  E5: Time (isolation proxy, OSD-202) → beta_T from brain 3-way

Causal stability score for (stressor S → tissue T pathway P):
  ICP_score(S,P) = 1 - CV(beta_S across environments sharing P)
  A high score (low CV) means S → P is environment-invariant = causal candidate.

Outputs
-------
  v8/causal/evaluation/icp_stressor_pathway_scores.csv
  v8/causal/evaluation/icp_dag_summary.json
  v8/causal/evaluation/dag_edges.csv   (node_from, node_to, icp_score, evidence)
"""
from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import variation

REPO_ROOT = Path(__file__).resolve().parents[2]
DECOMPOSE_DIR = REPO_ROOT / "v8" / "decompose" / "evaluation"
BRIDGE_DIR = REPO_ROOT / "v8" / "bridge" / "evaluation"
SIG_DIR = REPO_ROOT / "v8" / "intervene" / "signatures"
OUT_DIR = Path(__file__).resolve().parent / "evaluation"

TISSUES = ["thymus", "gastrocnemius", "skin", "eye", "liver", "kidney"]
STRESSORS = ["HLU", "IR", "HLUxIR", "T"]


def load_factorial_betas() -> dict[str, pd.DataFrame]:
    """Load per-gene factorial β tables for each analog."""
    analogs = {
        "spleen": DECOMPOSE_DIR / "factorial_spleen_betas.csv",
        "skin_analog": DECOMPOSE_DIR / "factorial_skin_analog_betas.csv",
        "brain": DECOMPOSE_DIR / "factorial_brain_betas.csv",
    }
    out = {}
    for name, path in analogs.items():
        if path.exists():
            out[name] = pd.read_csv(path)
    return out


def load_flight_signatures() -> dict[str, pd.DataFrame]:
    out = {}
    for tissue in TISSUES:
        p = SIG_DIR / f"{tissue}_ranked.csv"
        if p.exists():
            df = pd.read_csv(p)
            df["gene"] = df["gene"].astype(str).str.upper()
            out[tissue] = df
    return out


def build_stressor_environment_matrix(
    betas: dict[str, pd.DataFrame], flight: dict[str, pd.DataFrame]
) -> pd.DataFrame:
    """
    Build long DataFrame: (gene, environment, stressor, beta_value).

    Environments:
      - flight_{tissue}: Wald-z from ISS (surrogate for combined µg+GCR response)
      - analog_{analog}_{stressor}: β from factorial OLS (isolated stressor)
    """
    rows = []
    # Analog environments
    for analog_name, df in betas.items():
        for col in [c for c in df.columns if c.startswith("beta_") and c != "beta_intercept"]:
            stressor = col.replace("beta_", "")
            for _, row in df.iterrows():
                rows.append({
                    "gene": str(row["gene"]).upper(),
                    "environment": f"analog_{analog_name}",
                    "stressor": stressor,
                    "beta": float(row[col]),
                })

    # Flight environments (treat flight Wald-z as combined stressor = HLU+IR)
    for tissue, df in flight.items():
        for _, row in df.iterrows():
            rows.append({
                "gene": str(row["gene"]).upper(),
                "environment": f"flight_{tissue}",
                "stressor": "combined_flight",
                "beta": float(row["mean_z"]),
            })

    return pd.DataFrame(rows)


def icp_stability(env_matrix: pd.DataFrame, min_envs: int = 2) -> pd.DataFrame:
    """
    Per (stressor, gene) compute coefficient of variation across environments.
    Low CV → invariant → causal candidate.
    """
    records = []
    for (stressor, gene), grp in env_matrix.groupby(["stressor", "gene"]):
        vals = grp["beta"].dropna().values
        envs = grp["environment"].tolist()
        if len(vals) < min_envs:
            continue
        mean_b = float(np.mean(vals))
        std_b = float(np.std(vals))
        cv = std_b / (abs(mean_b) + 1e-9)
        icp_score = float(1.0 / (1.0 + cv))
        records.append({
            "stressor": stressor,
            "gene": gene,
            "n_environments": len(vals),
            "mean_beta": mean_b,
            "std_beta": std_b,
            "cv": cv,
            "icp_score": icp_score,
            "environments": "|".join(envs),
        })
    return pd.DataFrame(records).sort_values("icp_score", ascending=False)


def build_dag_edges(icp_df: pd.DataFrame, flight: dict[str, pd.DataFrame],
                    top_n: int = 50) -> pd.DataFrame:
    """
    For each tissue, find which stressor → tissue pathway/gene edges are
    most ICP-invariant. Cross-reference with flight signature magnitude.
    """
    edges = []
    for tissue, sig_df in flight.items():
        sig_top = set(
            sig_df.nlargest(top_n, "mean_z")["gene"].tolist() +
            sig_df.nsmallest(top_n, "mean_z")["gene"].tolist()
        )
        tissue_icp = icp_df[icp_df["gene"].isin(sig_top)].copy()
        if tissue_icp.empty:
            continue
        by_stressor = tissue_icp.groupby("stressor").agg(
            mean_icp=("icp_score", "mean"),
            n_genes=("gene", "count"),
            top_gene=("gene", lambda x: x.iloc[0]),
        ).reset_index()
        for _, row in by_stressor.iterrows():
            edges.append({
                "node_from": f"stressor:{row['stressor']}",
                "node_to": f"tissue:{tissue}",
                "icp_score": float(row["mean_icp"]),
                "n_genes": int(row["n_genes"]),
                "evidence": f"ICP across {row['n_genes']} genes; top={row['top_gene']}",
            })

    df_edges = pd.DataFrame(edges).sort_values("icp_score", ascending=False)
    return df_edges


def main() -> None:
    OUT_DIR.mkdir(parents=True, exist_ok=True)

    betas = load_factorial_betas()
    flight = load_flight_signatures()
    print(f"Loaded {len(betas)} analog factorial tables, {len(flight)} flight tissue signatures")

    env_matrix = build_stressor_environment_matrix(betas, flight)
    print(f"Environment matrix: {len(env_matrix)} rows, "
          f"{env_matrix['environment'].nunique()} environments, "
          f"{env_matrix['gene'].nunique()} genes")

    # ICP stability per (stressor, gene)
    icp_df = icp_stability(env_matrix, min_envs=2)
    out_path = OUT_DIR / "icp_stressor_pathway_scores.csv"
    icp_df.to_csv(out_path, index=False)
    print(f"ICP scores: {len(icp_df)} entries  →  {out_path}")

    # DAG edges
    edges_df = build_dag_edges(icp_df, flight, top_n=50)
    edges_df.to_csv(OUT_DIR / "dag_edges.csv", index=False)

    # Top causal stressor per tissue
    summary = {"n_genes": int(icp_df["gene"].nunique()),
                "n_environments": int(env_matrix["environment"].nunique()),
                "top_causal_edges": []}
    for _, row in edges_df.head(15).iterrows():
        summary["top_causal_edges"].append({
            "from": row["node_from"],
            "to": row["node_to"],
            "icp_score": round(row["icp_score"], 4),
            "n_genes": int(row["n_genes"]),
        })

    # Per-stressor aggregate ICP
    per_stressor = icp_df.groupby("stressor")["icp_score"].agg(["mean", "max", "count"])
    summary["stressor_aggregate_icp"] = per_stressor.round(4).to_dict("index")

    # Top 20 globally invariant genes
    top_genes = icp_df.groupby("gene")["icp_score"].mean().nlargest(20)
    summary["top20_invariant_genes"] = top_genes.round(4).to_dict()

    (OUT_DIR / "icp_dag_summary.json").write_text(json.dumps(summary, indent=2))
    print(f"Wrote {OUT_DIR}/icp_dag_summary.json")

    # Print readable summary
    print("\n=== Top causal stressor → tissue edges ===")
    for edge in summary["top_causal_edges"][:10]:
        print(f"  {edge['from']:35s} → {edge['to']:25s}  ICP={edge['icp_score']:.4f}  n={edge['n_genes']}")

    print("\n=== Stressor ICP aggregate ===")
    for s, vals in sorted(summary["stressor_aggregate_icp"].items(),
                           key=lambda x: -x[1]["mean"]):
        print(f"  {s:12s}  mean_ICP={vals['mean']:.4f}  max={vals['max']:.4f}  n={vals['count']}")

    print("\n=== Top 10 environment-invariant genes ===")
    for gene, score in list(summary["top20_invariant_genes"].items())[:10]:
        print(f"  {gene:12s}  mean_ICP={score:.4f}")


if __name__ == "__main__":
    main()
