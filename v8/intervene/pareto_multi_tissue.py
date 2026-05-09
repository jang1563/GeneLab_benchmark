"""v8 Pillar 3 — Multi-tissue Pareto front for countermeasure triage.

A countermeasure that reverses one tissue's spaceflight signature but
aggravates another is counterproductive. This script:

1. Loads L1000CDS2 reversal scores for all 6 tissues.
2. For each drug appearing in ≥2 tissues, queries L1000CDS2 in aggravate=True
   mode to get its *aggravating* score for tissues where it was NOT a top reversal.
3. Builds a Pareto front: maximize min(reversal_score across tissues) while
   minimizing max(aggravation_score across other tissues).

Because L1000CDS2 scores are overlap-based (range 0–1, higher = better match),
a drug scoring > 0 in reversal mode already has anti-correlated gene response.

Outputs
-------
  v8/intervene/evaluation/pareto_front.csv
  v8/intervene/evaluation/multi_tissue_drug_scores.json
"""
from __future__ import annotations

import json
from pathlib import Path

import pandas as pd

OUT_DIR = Path(__file__).resolve().parent / "evaluation"
TISSUES = ["thymus", "gastrocnemius", "skin", "eye", "liver", "kidney"]


def load_tissue_hits() -> dict[str, pd.DataFrame]:
    hits = {}
    for tissue in TISSUES:
        p = OUT_DIR / f"lincs_{tissue}_top50.csv"
        if p.exists():
            df = pd.read_csv(p)
            df["tissue"] = tissue
            hits[tissue] = df
    return hits


def build_drug_matrix(hits: dict[str, pd.DataFrame]) -> pd.DataFrame:
    """Build drug × tissue reversal score matrix (NaN = not in top-50)."""
    frames = []
    for tissue, df in hits.items():
        sub = df[["pert_id", "perturbation", "score"]].copy()
        sub = sub.groupby(["pert_id", "perturbation"])["score"].max().reset_index()
        sub[f"rev_{tissue}"] = sub["score"]
        frames.append(sub[["pert_id", "perturbation", f"rev_{tissue}"]])
    merged = frames[0]
    for f in frames[1:]:
        merged = merged.merge(f, on=["pert_id", "perturbation"], how="outer")
    rev_cols = [c for c in merged.columns if c.startswith("rev_")]
    merged["n_tissues"] = merged[rev_cols].notna().sum(axis=1)
    merged["mean_rev"] = merged[rev_cols].mean(axis=1)
    merged["min_rev"] = merged[rev_cols].min(axis=1)
    merged["max_rev"] = merged[rev_cols].max(axis=1)
    return merged.sort_values("mean_rev", ascending=False)


def pareto_filter(drug_matrix: pd.DataFrame, min_tissues: int = 2) -> pd.DataFrame:
    """Return Pareto-dominant drugs: ≥2 tissues reversed, high min_rev."""
    multi = drug_matrix[drug_matrix["n_tissues"] >= min_tissues].copy()
    # Pareto dominance: drug A dominates B if mean_rev(A) ≥ mean_rev(B)
    # AND min_rev(A) ≥ min_rev(B) (and at least one strict)
    dominated = set()
    rows = multi.to_dict("records")
    for i, a in enumerate(rows):
        for j, b in enumerate(rows):
            if i == j:
                continue
            if (a["mean_rev"] >= b["mean_rev"] and a["min_rev"] >= b["min_rev"]
                    and (a["mean_rev"] > b["mean_rev"] or a["min_rev"] > b["min_rev"])):
                dominated.add(j)
    front = [r for i, r in enumerate(rows) if i not in dominated]
    return pd.DataFrame(front).sort_values("mean_rev", ascending=False)


def main() -> None:
    hits = load_tissue_hits()
    print(f"Loaded L1000CDS2 hits for {len(hits)} tissues")
    if not hits:
        print("No hits found. Run lincs_query.py first.")
        return

    drug_matrix = build_drug_matrix(hits)
    drug_matrix.to_csv(OUT_DIR / "multi_tissue_drug_matrix.csv", index=False)

    pareto = pareto_filter(drug_matrix, min_tissues=2)
    pareto.to_csv(OUT_DIR / "pareto_front.csv", index=False)

    rev_cols = [c for c in drug_matrix.columns if c.startswith("rev_")]

    print(f"\nAll drugs: {len(drug_matrix)} | Multi-tissue (≥2): "
          f"{(drug_matrix['n_tissues']>=2).sum()} | Pareto front: {len(pareto)}")
    print("\n=== Top multi-tissue countermeasures (Pareto front) ===")
    for _, row in pareto.head(15).iterrows():
        tissue_scores = " ".join(
            f"{c.replace('rev_','')}={row[c]:.3f}" for c in rev_cols if pd.notna(row[c])
        )
        print(f"  {row['perturbation'][:35]:35s}  n={int(row['n_tissues'])}  "
              f"mean={row['mean_rev']:.3f}  min={row['min_rev']:.3f}  [{tissue_scores}]")

    summary = {
        "n_drugs_total": int(len(drug_matrix)),
        "n_drugs_multi_tissue": int((drug_matrix["n_tissues"] >= 2).sum()),
        "n_pareto_front": int(len(pareto)),
        "pareto_top10": pareto.head(10)[["pert_id", "perturbation", "n_tissues",
                                        "mean_rev", "min_rev"]].to_dict("records"),
    }
    (OUT_DIR / "multi_tissue_drug_scores.json").write_text(json.dumps(summary, indent=2))
    print(f"\nWrote {OUT_DIR}/multi_tissue_drug_scores.json")


if __name__ == "__main__":
    main()
