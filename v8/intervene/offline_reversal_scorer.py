"""v8 Pillar 3 (INTERVENE) — offline drug reversal scoring.

A LINCS-independent first pass. For each tissue, rank drugs by how strongly
their DGIdb target set overlaps the tissue's spaceflight signature. Serves as:
  - baseline for the LINCS query (network-dependent; to be added next)
  - orthogonal sanity check for the forthcoming LINCS rankings

Method
------
  For each drug D with DGIdb target gene set T_D:
    For each tissue t:
      up_t, dn_t = top-K up- and down-regulated genes in tissue t
      fg_hits = |T_D ∩ up_t| + |T_D ∩ dn_t|      (any spaceflight-responsive)
      bg_hits = |T_D|                             (targets overall)
      coherence = (|T_D ∩ up_t| - |T_D ∩ dn_t|) / (|T_D ∩ (up_t ∪ dn_t)| + 1)
      enrichment_p = hypergeometric(k=fg_hits, M=n_ranked, n=|up_t∪dn_t|, N=|T_D|)

  Drugs targeting many flight-responsive genes in coherent direction rank highest.
  A true *reversal* score requires knowing drug-gene interaction *direction*
  (inhibitor vs agonist) and gene up/down direction — flagged for LINCS follow-up.

Inputs
------
  v5/evaluation/consensus_biomarker_panel.json     — 20-gene panel with drugs
  v5/evaluation/drug_targets.json                  — 271 druggable genes, 3154 interactions
  v8/intervene/signatures/{tissue}_ranked.csv      — per-tissue signed-z rank

Outputs
-------
  v8/intervene/evaluation/offline_reversal_{tissue}.csv
  v8/intervene/evaluation/offline_reversal_summary.json
"""
from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import hypergeom

REPO_ROOT = Path(__file__).resolve().parents[2]
V5_CONSENSUS = REPO_ROOT / "v5" / "evaluation" / "consensus_biomarker_panel.json"
V5_DRUGS = REPO_ROOT / "v5" / "evaluation" / "drug_targets.json"
SIG_DIR = Path(__file__).resolve().parent / "signatures"
OUT_DIR = Path(__file__).resolve().parent / "evaluation"

TOP_K = 150
TISSUES = ["thymus", "gastrocnemius", "skin", "eye", "liver", "kidney"]


def build_drug_target_map() -> dict[str, set[str]]:
    """Invert the gene→drugs mapping from v5 consensus panel into drug→genes.

    Falls back to empty if sources missing.
    """
    drug_to_targets: dict[str, set[str]] = {}
    if V5_CONSENSUS.exists():
        panel = json.loads(V5_CONSENSUS.read_text())
        for entry in panel.get("panel", []):
            gene = (entry.get("human_symbol") or entry.get("gene") or "").upper()
            for d in entry.get("drugs", []):
                name = str(d.get("drug", "")).upper().strip()
                if name:
                    drug_to_targets.setdefault(name, set()).add(gene)
    return drug_to_targets


def load_signature(tissue: str, top_k: int = TOP_K) -> tuple[set[str], set[str], int]:
    p = SIG_DIR / f"{tissue}_ranked.csv"
    if not p.exists():
        return set(), set(), 0
    df = pd.read_csv(p)
    df = df.sort_values("mean_z", ascending=False)
    up = set(df.head(top_k)["gene"].astype(str).str.upper())
    dn = set(df.tail(top_k)["gene"].astype(str).str.upper())
    return up, dn, int(len(df))


def score_drug(T_D: set[str], up: set[str], dn: set[str], n_ranked: int) -> dict:
    up_hits = len(T_D & up)
    dn_hits = len(T_D & dn)
    fg_union = up | dn
    fg_hits = len(T_D & fg_union)
    bg = len(T_D)
    if bg == 0 or n_ranked == 0 or not fg_union:
        return {"up_hits": up_hits, "dn_hits": dn_hits, "bg": bg,
                "coherence": np.nan, "pval": np.nan}
    coherence = (up_hits - dn_hits) / (fg_hits + 1)
    # hypergeom p = P(X >= fg_hits)
    pval = float(hypergeom.sf(fg_hits - 1, n_ranked, len(fg_union), bg))
    return {"up_hits": up_hits, "dn_hits": dn_hits, "bg": bg,
            "coherence": coherence, "pval": pval}


def main() -> None:
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    drug_targets = build_drug_target_map()
    print(f"Loaded {len(drug_targets)} drugs with DGIdb targets from v5 consensus panel.")

    summary: dict[str, dict] = {"drugs_scored": len(drug_targets), "tissues": {}}

    for tissue in TISSUES:
        up, dn, n_ranked = load_signature(tissue)
        if not up:
            print(f"[{tissue}] SKIP — no signature available yet")
            continue
        rows = []
        for drug, T_D in drug_targets.items():
            s = score_drug(T_D, up, dn, n_ranked)
            s.update({"drug": drug, "targets": "|".join(sorted(T_D))})
            rows.append(s)
        df = pd.DataFrame(rows)
        df = df[df["bg"] >= 1]
        df = df.sort_values(["pval", "coherence"], ascending=[True, False])
        df.to_csv(OUT_DIR / f"offline_reversal_{tissue}.csv", index=False)

        top = df.head(10)
        print(f"\n[{tissue}] top-10 drugs by hypergeom p (n_ranked={n_ranked}, up∪dn={len(up|dn)})")
        for _, r in top.iterrows():
            print(f"  {r['drug']:<30} p={r['pval']:.2e}  up={r['up_hits']} dn={r['dn_hits']}  "
                  f"coh={r['coherence']:+.2f}  bg={r['bg']}  targets={r['targets'][:50]}")

        summary["tissues"][tissue] = {
            "n_up": len(up), "n_dn": len(dn), "n_ranked": n_ranked,
            "top5": top.head(5)[["drug", "pval", "up_hits", "dn_hits", "coherence", "bg"]].to_dict(orient="records"),
        }

    (OUT_DIR / "offline_reversal_summary.json").write_text(json.dumps(summary, indent=2, default=str))
    print(f"\nWrote {OUT_DIR}/offline_reversal_summary.json")


if __name__ == "__main__":
    main()
