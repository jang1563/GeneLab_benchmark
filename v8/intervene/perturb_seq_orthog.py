"""v8 Pillar 3 — Perturb-seq/CRISPR orthogonal validation.

Queries the Enrichr LINCS_L1000_CRISPR_KO_Consensus_Sigs library with each
tissue's up/down signatures. For each Pillar 3b top drug hit, checks whether
its primary target (or related pathway members) also appear as top CRISPR
KO reversers — orthogonal validation.

Primary drug targets (from ChEMBL/DrugBank):
  AZD-5438, CGP-60474, AT-7519 → CDK1, CDK2, CDK5, CDK9
  Geldanamycin, NVP-AUY922     → HSP90AA1, HSP90AB1
  Dorsomorphin                 → PRKAA1, PRKAA2 (AMPK α1/α2), BMPR1A, BMPR1B, ALK2
  KU-0063794                   → MTOR, RPTOR, RICTOR
  GDC-0980                     → PIK3CA, PIK3CB, MTOR
  PD-184352                    → MAP2K1, MAP2K2 (MEK1/2)
  Mitoxantrone                 → TOP2A, TOP2B

A hit in the CRISPR KO library for any of these targets — in reversal direction —
orthogonally supports the chemical-reversal finding.

Enrichr API: https://maayanlab.cloud/Enrichr/
  POST /addList  →  userListId
  GET /enrich?userListId=X&backgroundType=LINCS_L1000_CRISPR_KO_Consensus_Sigs

Outputs
-------
  v8/intervene/evaluation/crispr_orthog_{tissue}_up.csv     (KO Down hits = reversal of UP)
  v8/intervene/evaluation/crispr_orthog_{tissue}_dn.csv     (KO Up hits = reversal of DN)
  v8/intervene/evaluation/crispr_orthog_summary.json
"""
from __future__ import annotations

import argparse
import io
import json
import time
from pathlib import Path
from typing import Any

import pandas as pd
import requests

SIG_DIR = Path(__file__).resolve().parent / "signatures"
OUT_DIR = Path(__file__).resolve().parent / "evaluation"

ENRICHR_ADD = "https://maayanlab.cloud/Enrichr/addList"
ENRICHR_ENRICH = "https://maayanlab.cloud/Enrichr/enrich"
LIBRARY = "LINCS_L1000_CRISPR_KO_Consensus_Sigs"

TISSUES = ["thymus", "gastrocnemius", "skin", "eye", "liver", "kidney"]

# Drug → primary targets (from ChEMBL/MoA)
DRUG_TARGETS = {
    "AZD-5438": {"CDK1", "CDK2", "CDK9"},
    "CGP-60474": {"CDK1", "CDK2", "CDK5"},
    "AT-7519": {"CDK1", "CDK2", "CDK4", "CDK6", "CDK9"},
    "Geldanamycin": {"HSP90AA1", "HSP90AB1"},
    "NVP-AUY922": {"HSP90AA1", "HSP90AB1"},
    "Dorsomorphin": {"PRKAA1", "PRKAA2", "BMPR1A", "BMPR1B", "ACVR1"},
    "KU-0063794": {"MTOR", "RPTOR", "RICTOR"},
    "GDC-0980": {"PIK3CA", "PIK3CB", "MTOR"},
    "PD-184352": {"MAP2K1", "MAP2K2"},
    "Mitoxantrone": {"TOP2A", "TOP2B"},
}


def _safe_name(value: str) -> str:
    return "".join(c if c.isalnum() or c in {"-", "_"} else "_" for c in value)


def _write_raw_record(raw_dir: Path | None, name: str, record: dict[str, Any]) -> None:
    if raw_dir is None:
        return
    raw_dir.mkdir(parents=True, exist_ok=True)
    (raw_dir / f"{_safe_name(name)}.json").write_text(json.dumps(record, indent=2, sort_keys=True))


def upload_list(
    genes: list[str],
    description: str,
    raw_dir: Path | None = None,
    raw_name: str | None = None,
) -> int | None:
    payload = {"list": (None, "\n".join(genes)), "description": (None, description)}
    try:
        r = requests.post(ENRICHR_ADD, files=payload, timeout=30)
        r.raise_for_status()
        data = r.json()
        _write_raw_record(
            raw_dir,
            raw_name or f"enrichr_{description}_addList",
            {
                "endpoint": ENRICHR_ADD,
                "request": {"description": description, "genes": genes},
                "response_status_code": r.status_code,
                "response_headers": {
                    key: value
                    for key, value in r.headers.items()
                    if key.lower() in {"content-type", "date", "server"}
                },
                "response_json": data,
            },
        )
        return data["userListId"]
    except Exception as e:
        print(f"  Enrichr addList failed: {e}")
        _write_raw_record(
            raw_dir,
            raw_name or f"enrichr_{description}_addList",
            {"endpoint": ENRICHR_ADD, "request": {"description": description, "genes": genes}, "error": str(e)},
        )
        return None


def enrich(
    user_list_id: int,
    library: str,
    raw_dir: Path | None = None,
    raw_name: str | None = None,
) -> list[list]:
    try:
        r = requests.get(ENRICHR_ENRICH,
                         params={"userListId": user_list_id, "backgroundType": library},
                         timeout=60)
        r.raise_for_status()
        data = r.json()
        _write_raw_record(
            raw_dir,
            raw_name or f"enrichr_{user_list_id}_{library}",
            {
                "endpoint": ENRICHR_ENRICH,
                "request_params": {"userListId": user_list_id, "backgroundType": library},
                "response_status_code": r.status_code,
                "response_headers": {
                    key: value
                    for key, value in r.headers.items()
                    if key.lower() in {"content-type", "date", "server"}
                },
                "response_json": data,
            },
        )
        return data.get(library, [])
    except Exception as e:
        print(f"  Enrichr enrich failed: {e}")
        _write_raw_record(
            raw_dir,
            raw_name or f"enrichr_{user_list_id}_{library}",
            {
                "endpoint": ENRICHR_ENRICH,
                "request_params": {"userListId": user_list_id, "backgroundType": library},
                "error": str(e),
            },
        )
        return []


def parse_hits(raw: list[list]) -> pd.DataFrame:
    """Enrichr result row: [rank, term, pvalue, odds, combined_score, overlap_genes, adj_p, old_p, old_adj_p]."""
    rows = []
    for r in raw:
        if len(r) < 7:
            continue
        term = r[1]
        parts = term.rsplit(" ", 1)
        target = parts[0] if len(parts) == 2 else term
        direction = parts[1] if len(parts) == 2 else ""
        rows.append({
            "term": term,
            "target": target,
            "ko_direction": direction,  # "Up" or "Down"
            "pvalue": r[2],
            "odds": r[3],
            "combined_score": r[4],
            "overlap_genes": "|".join(r[5]) if isinstance(r[5], list) else str(r[5]),
            "adj_pvalue": r[6],
        })
    return pd.DataFrame(rows).sort_values("combined_score", ascending=False)


def check_drug_validation(crispr_hits: dict[str, pd.DataFrame], tissue: str) -> dict:
    """For each drug in DRUG_TARGETS, check if its target appears as top CRISPR hit."""
    validated = {}
    # For spaceflight UP: CRISPR KO Down = reversal
    # For spaceflight DN: CRISPR KO Up = reversal
    up_rev = crispr_hits.get(f"{tissue}_up", pd.DataFrame())
    dn_rev = crispr_hits.get(f"{tissue}_dn", pd.DataFrame())
    for drug, targets in DRUG_TARGETS.items():
        drug_validated = []
        for _, row in up_rev.iterrows():
            if row["ko_direction"] == "Down" and row["target"] in targets:
                drug_validated.append({
                    "target": row["target"], "ko_dir": "Down",
                    "reverses": "UP genes", "adj_p": row["adj_pvalue"],
                    "overlap": row["overlap_genes"],
                })
        for _, row in dn_rev.iterrows():
            if row["ko_direction"] == "Up" and row["target"] in targets:
                drug_validated.append({
                    "target": row["target"], "ko_dir": "Up",
                    "reverses": "DN genes", "adj_p": row["adj_pvalue"],
                    "overlap": row["overlap_genes"],
                })
        if drug_validated:
            validated[drug] = drug_validated
    return validated


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--out-dir", default=str(OUT_DIR), help="Directory for parsed CSV/JSON outputs.")
    parser.add_argument("--raw-dir", default=None, help="Optional directory for raw API response archives.")
    args = parser.parse_args()

    out_dir = Path(args.out_dir)
    raw_dir = Path(args.raw_dir) if args.raw_dir else None
    out_dir.mkdir(parents=True, exist_ok=True)
    summary: dict = {"library": LIBRARY, "tissues": {}}
    all_hits: dict[str, pd.DataFrame] = {}

    for tissue in TISSUES:
        up_genes = (SIG_DIR / f"{tissue}_up.txt").read_text().strip().splitlines()
        dn_genes = (SIG_DIR / f"{tissue}_dn.txt").read_text().strip().splitlines()
        up_genes = [g.strip() for g in up_genes if g.strip()][:100]
        dn_genes = [g.strip() for g in dn_genes if g.strip()][:100]

        print(f"\n[{tissue}]  up={len(up_genes)}  dn={len(dn_genes)}")

        # Upload UP list → find KO signatures with "Down" direction (= reversers)
        ulid_up = upload_list(
            up_genes,
            f"{tissue}_up_v8",
            raw_dir=raw_dir,
            raw_name=f"enrichr_{tissue}_up_addList",
        )
        if ulid_up:
            time.sleep(1)
            raw_up = enrich(
                ulid_up,
                LIBRARY,
                raw_dir=raw_dir,
                raw_name=f"enrichr_{tissue}_up_enrich",
            )
            df_up = parse_hits(raw_up)
            # Keep only "Down" direction for UP genes (reversal)
            df_up_rev = df_up[df_up["ko_direction"] == "Down"].head(50)
            df_up_rev.to_csv(out_dir / f"crispr_orthog_{tissue}_up.csv", index=False)
            all_hits[f"{tissue}_up"] = df_up_rev
            top3 = df_up_rev.head(3)["target"].tolist()
            print(f"  UP reversers (KO→Down): {top3}")

        time.sleep(2)

        # Upload DN list → find KO signatures with "Up" direction (= reversers)
        ulid_dn = upload_list(
            dn_genes,
            f"{tissue}_dn_v8",
            raw_dir=raw_dir,
            raw_name=f"enrichr_{tissue}_dn_addList",
        )
        if ulid_dn:
            time.sleep(1)
            raw_dn = enrich(
                ulid_dn,
                LIBRARY,
                raw_dir=raw_dir,
                raw_name=f"enrichr_{tissue}_dn_enrich",
            )
            df_dn = parse_hits(raw_dn)
            df_dn_rev = df_dn[df_dn["ko_direction"] == "Up"].head(50)
            df_dn_rev.to_csv(out_dir / f"crispr_orthog_{tissue}_dn.csv", index=False)
            all_hits[f"{tissue}_dn"] = df_dn_rev
            top3 = df_dn_rev.head(3)["target"].tolist()
            print(f"  DN reversers (KO→Up):   {top3}")

        time.sleep(2)

        # Drug validation: which Pillar-3b drug targets surface as CRISPR reversers?
        validated = check_drug_validation(all_hits, tissue)
        summary["tissues"][tissue] = {
            "n_up_reversers": int(len(all_hits.get(f"{tissue}_up", pd.DataFrame()))),
            "n_dn_reversers": int(len(all_hits.get(f"{tissue}_dn", pd.DataFrame()))),
            "validated_drugs": validated,
        }
        if validated:
            print(f"  ✓ Validated drugs via CRISPR: {list(validated.keys())}")

    (out_dir / "crispr_orthog_summary.json").write_text(json.dumps(summary, indent=2, default=str))
    print(f"\nWrote {out_dir}/crispr_orthog_summary.json")


if __name__ == "__main__":
    main()
