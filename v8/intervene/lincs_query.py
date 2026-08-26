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

import argparse
import json
import time
from pathlib import Path
from typing import Any

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


def _safe_name(value: str) -> str:
    return "".join(c if c.isalnum() or c in {"-", "_"} else "_" for c in value)


def _write_raw_record(raw_dir: Path | None, name: str, record: dict[str, Any]) -> None:
    if raw_dir is None:
        return
    raw_dir.mkdir(parents=True, exist_ok=True)
    (raw_dir / f"{_safe_name(name)}.json").write_text(json.dumps(record, indent=2, sort_keys=True))


def _l1000_payload(up: list[str], dn: list[str]) -> dict[str, Any]:
    return {
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


def query_l1000cds2(
    up: list[str],
    dn: list[str],
    n_results: int = 50,
    raw_dir: Path | None = None,
    raw_name: str | None = None,
) -> list[dict]:
    """Query L1000CDS2 for signature reversal.

    aggravate=False → find perturbations that REVERSE the query signature.
    Returns list of top hits ranked by score (most reversing first).
    """
    payload = _l1000_payload(up, dn)
    headers = {"Content-Type": "application/json", "Accept": "application/json"}
    try:
        r = requests.post(L1000CDS2_URL, json=payload, headers=headers, timeout=60)
        r.raise_for_status()
        data = r.json()
        _write_raw_record(
            raw_dir,
            raw_name or "l1000cds2_response",
            {
                "endpoint": L1000CDS2_URL,
                "request_payload": payload,
                "response_status_code": r.status_code,
                "response_headers": {
                    key: value
                    for key, value in r.headers.items()
                    if key.lower() in {"content-type", "date", "server"}
                },
                "response_json": data,
            },
        )
        return data.get("topMeta", [])[:n_results]
    except requests.exceptions.Timeout:
        print("  TIMEOUT on L1000CDS2 API")
        _write_raw_record(
            raw_dir,
            raw_name or "l1000cds2_response",
            {"endpoint": L1000CDS2_URL, "request_payload": payload, "error": "timeout"},
        )
        return []
    except requests.exceptions.HTTPError as e:
        print(f"  HTTP error: {e} — {r.text[:300]}")
        _write_raw_record(
            raw_dir,
            raw_name or "l1000cds2_response",
            {
                "endpoint": L1000CDS2_URL,
                "request_payload": payload,
                "response_status_code": getattr(r, "status_code", None),
                "response_text": getattr(r, "text", "")[:10000],
                "error": str(e),
            },
        )
        return []
    except Exception as e:
        print(f"  Error: {e}")
        _write_raw_record(
            raw_dir,
            raw_name or "l1000cds2_response",
            {"endpoint": L1000CDS2_URL, "request_payload": payload, "error": str(e)},
        )
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
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--out-dir", default=str(OUT_DIR), help="Directory for parsed CSV/JSON outputs.")
    parser.add_argument("--raw-dir", default=None, help="Optional directory for raw API response archives.")
    args = parser.parse_args()

    out_dir = Path(args.out_dir)
    raw_dir = Path(args.raw_dir) if args.raw_dir else None
    out_dir.mkdir(parents=True, exist_ok=True)
    summary: dict = {"api": L1000CDS2_URL, "tissues": {}}

    for tissue in TISSUES:
        up, dn = _load_genes(tissue)
        print(f"\n[{tissue}]  up={len(up)}  dn={len(dn)}")
        raw = query_l1000cds2(
            up,
            dn,
            n_results=50,
            raw_dir=raw_dir,
            raw_name=f"l1000cds2_{tissue}_response",
        )
        if not raw:
            print("  No results returned.")
            summary["tissues"][tissue] = {"n_results": 0, "note": "no API results"}
            continue

        df = _parse_results(raw)
        if df.empty:
            summary["tissues"][tissue] = {"n_results": 0, "note": "parse failed"}
            continue

        out_path = out_dir / f"lincs_{tissue}_top50.csv"
        df.head(50).to_csv(out_path, index=False)

        top5_rev = df.head(5)[["perturbation", "pert_id", "cell_line", "score"]].to_dict("records")
        summary["tissues"][tissue] = {
            "n_results": len(df),
            "top5_reversal": top5_rev,
        }
        print(f"  Top reversal: " + ", ".join(r["perturbation"] for r in top5_rev[:3]))
        print(f"  Wrote {out_path}")
        time.sleep(2)  # rate-limit courtesy

    (out_dir / "lincs_summary.json").write_text(json.dumps(summary, indent=2))
    print(f"\nWrote {out_dir}/lincs_summary.json")


if __name__ == "__main__":
    main()
