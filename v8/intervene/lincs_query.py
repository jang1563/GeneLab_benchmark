"""v8 Pillar 3 (INTERVENE) — LINCS L2S2 / SigCom signature reversal query.

Queries the SigCom LINCS API (https://maayanlab.cloud/sigcom-lincs) with each
tissue's top-150 up/down gene signature to rank small-molecule perturbations
by connectivity score (negative = reversal).

API: POST /api/v1/data/enrich/overlap
  body: {"up": [genes], "down": [genes], "filter": {"limit": 100}}

Outputs
-------
  v8/intervene/evaluation/lincs_{tissue}_top50.csv
  v8/intervene/evaluation/lincs_summary.json
"""
from __future__ import annotations

import json
import time
from pathlib import Path

import pandas as pd
import requests

SIG_DIR = Path(__file__).resolve().parent / "signatures"
OUT_DIR = Path(__file__).resolve().parent / "evaluation"

L1000CDS2_URL = "https://maayanlab.cloud/L1000CDS2/query"

TISSUES = ["thymus", "gastrocnemius", "skin", "eye", "liver", "kidney"]


def _load_genes(tissue: str) -> tuple[list[str], list[str]]:
    up = (SIG_DIR / f"{tissue}_up.txt").read_text().strip().splitlines()
    dn = (SIG_DIR / f"{tissue}_dn.txt").read_text().strip().splitlines()
    return [g.strip() for g in up if g.strip()], [g.strip() for g in dn if g.strip()]


def query_l1000cds2(up: list[str], dn: list[str], n_results: int = 50) -> list[dict]:
    """Query L1000CDS2 for signature reversal.

    aggravate=False → find perturbations that REVERSE the query signature.
    Returns list of top hits ranked by score (most reversing first).
    """
    payload = {
        "data": {"upGenes": up[:50], "dnGenes": dn[:50]},
        "config": {
            "aggravate": False,
            "searchMethod": "geneSet",
            "share": False,
            "combination": False,
            "db-version": "latest",
        },
        "metadata": [],
    }
    headers = {"Content-Type": "application/json", "Accept": "application/json"}
    try:
        r = requests.post(L1000CDS2_URL, json=payload, headers=headers, timeout=60)
        r.raise_for_status()
        data = r.json()
        return data.get("topMeta", [])[:n_results]
    except requests.exceptions.Timeout:
        print("  TIMEOUT on L1000CDS2 API")
        return []
    except requests.exceptions.HTTPError as e:
        print(f"  HTTP error: {e} — {r.text[:300]}")
        return []
    except Exception as e:
        print(f"  Error: {e}")
        return []


def _parse_results(raw: list[dict]) -> pd.DataFrame:
    rows = []
    for item in raw:
        rows.append({
            "perturbation": item.get("pert_desc", ""),
            "pert_id": item.get("pert_id", ""),
            "cell_line": item.get("cell_id", ""),
            "dose_um": item.get("pert_dose", ""),
            "time_h": item.get("pert_time", ""),
            "score": item.get("score", None),
            "pubchem_id": item.get("pubchem_id", ""),
            "sig_id": item.get("sig_id", ""),
            "overlap_up_dn": "|".join(item.get("overlap", {}).get("up/dn", [])),
            "overlap_dn_up": "|".join(item.get("overlap", {}).get("dn/up", [])),
        })
    df = pd.DataFrame(rows)
    if df.empty:
        return df
    df["score"] = pd.to_numeric(df["score"], errors="coerce")
    return df.sort_values("score").reset_index(drop=True)


def main() -> None:
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    summary: dict = {"api": L1000CDS2_URL, "tissues": {}}

    for tissue in TISSUES:
        up, dn = _load_genes(tissue)
        print(f"\n[{tissue}]  up={len(up)}  dn={len(dn)}")
        raw = query_l1000cds2(up, dn, n_results=50)
        if not raw:
            print("  No results returned.")
            summary["tissues"][tissue] = {"n_results": 0, "note": "no API results"}
            continue

        df = _parse_results(raw)
        if df.empty:
            summary["tissues"][tissue] = {"n_results": 0, "note": "parse failed"}
            continue

        out_path = OUT_DIR / f"lincs_{tissue}_top50.csv"
        df.head(50).to_csv(out_path, index=False)

        top5_rev = df.head(5)[["perturbation", "pert_id", "cell_line", "score"]].to_dict("records")
        summary["tissues"][tissue] = {
            "n_results": len(df),
            "top5_reversal": top5_rev,
        }
        print(f"  Top reversal: " + ", ".join(r["perturbation"] for r in top5_rev[:3]))
        print(f"  Wrote {out_path}")
        time.sleep(2)  # rate-limit courtesy

    (OUT_DIR / "lincs_summary.json").write_text(json.dumps(summary, indent=2))
    print(f"\nWrote {OUT_DIR}/lincs_summary.json")


if __name__ == "__main__":
    main()
