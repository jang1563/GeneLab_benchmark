"""Build compact hashes for the current INTERVENE API-result snapshot.

This script does not call L1000CDS2 or Enrichr. It records deterministic hashes
for the query payloads that would be sent, plus hashes for the tracked parsed
outputs. Raw API response dumps remain outside Git unless explicitly archived
for a beta/frozen release.
"""
from __future__ import annotations

import hashlib
import json
from pathlib import Path
from typing import Any

SIG_DIR = Path(__file__).resolve().parent / "signatures"
OUT_DIR = Path(__file__).resolve().parent / "evaluation"

L1000CDS2_URL = "https://maayanlab.cloud/L1000CDS2/query"
ENRICHR_ADD = "https://maayanlab.cloud/Enrichr/addList"
ENRICHR_ENRICH = "https://maayanlab.cloud/Enrichr/enrich"
LIBRARY = "LINCS_L1000_CRISPR_KO_Consensus_Sigs"
TISSUES = ["thymus", "gastrocnemius", "skin", "eye", "liver", "kidney"]


def _sha256_bytes(data: bytes) -> str:
    return hashlib.sha256(data).hexdigest()


def _sha256_file(path: Path) -> str | None:
    if not path.exists():
        return None
    h = hashlib.sha256()
    with path.open("rb") as f:
        for chunk in iter(lambda: f.read(1024 * 1024), b""):
            h.update(chunk)
    return h.hexdigest()


def _canonical_json_hash(payload: dict[str, Any]) -> str:
    blob = json.dumps(payload, sort_keys=True, separators=(",", ":")).encode("utf-8")
    return _sha256_bytes(blob)


def _load_genes(path: Path, limit: int | None = None) -> list[str]:
    genes = [g.strip() for g in path.read_text().splitlines() if g.strip()]
    return genes[:limit] if limit is not None else genes


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


def _enrichr_payload_hash(genes: list[str], description: str) -> str:
    payload = {"description": description, "genes": genes}
    return _canonical_json_hash(payload)


def main() -> None:
    OUT_DIR.mkdir(parents=True, exist_ok=True)

    l1000_tissues = {}
    for tissue in TISSUES:
        up = _load_genes(SIG_DIR / f"{tissue}_up.txt")
        dn = _load_genes(SIG_DIR / f"{tissue}_dn.txt")
        payload = _l1000_payload(up, dn)
        out_csv = OUT_DIR / f"lincs_{tissue}_top50.csv"
        l1000_tissues[tissue] = {
            "up_signature_sha256": _canonical_json_hash({"genes": up}),
            "dn_signature_sha256": _canonical_json_hash({"genes": dn}),
            "payload_sha256": _canonical_json_hash(payload),
            "payload_top_n_per_direction": 50,
            "parsed_output": str(out_csv.relative_to(OUT_DIR.parent.parent)),
            "parsed_output_sha256": _sha256_file(out_csv),
        }

    crispr_tissues = {}
    for tissue in TISSUES:
        up = _load_genes(SIG_DIR / f"{tissue}_up.txt", limit=100)
        dn = _load_genes(SIG_DIR / f"{tissue}_dn.txt", limit=100)
        up_csv = OUT_DIR / f"crispr_orthog_{tissue}_up.csv"
        dn_csv = OUT_DIR / f"crispr_orthog_{tissue}_dn.csv"
        crispr_tissues[tissue] = {
            "up_query": {
                "description": f"{tissue}_up_v8",
                "n_genes": len(up),
                "payload_sha256": _enrichr_payload_hash(up, f"{tissue}_up_v8"),
                "parsed_output": str(up_csv.relative_to(OUT_DIR.parent.parent)),
                "parsed_output_sha256": _sha256_file(up_csv),
            },
            "dn_query": {
                "description": f"{tissue}_dn_v8",
                "n_genes": len(dn),
                "payload_sha256": _enrichr_payload_hash(dn, f"{tissue}_dn_v8"),
                "parsed_output": str(dn_csv.relative_to(OUT_DIR.parent.parent)),
                "parsed_output_sha256": _sha256_file(dn_csv),
            },
        }

    manifest = {
        "schema_version": "0.1.0",
        "snapshot_id": "intervene.api_snapshot.current",
        "api_refresh": "not_recalled_by_this_script",
        "l1000cds2": {
            "endpoint": L1000CDS2_URL,
            "config": {
                "aggravate": False,
                "searchMethod": "geneSet",
                "db-version": "latest",
            },
            "summary_sha256": _sha256_file(OUT_DIR / "lincs_summary.json"),
            "tissues": l1000_tissues,
        },
        "enrichr_crispr": {
            "add_list_endpoint": ENRICHR_ADD,
            "enrich_endpoint": ENRICHR_ENRICH,
            "library": LIBRARY,
            "summary_sha256": _sha256_file(OUT_DIR / "crispr_orthog_summary.json"),
            "tissues": crispr_tissues,
        },
        "multi_tissue_outputs": {
            "pareto_front_sha256": _sha256_file(OUT_DIR / "pareto_front.csv"),
            "multi_tissue_drug_matrix_sha256": _sha256_file(OUT_DIR / "multi_tissue_drug_matrix.csv"),
            "multi_tissue_drug_scores_sha256": _sha256_file(OUT_DIR / "multi_tissue_drug_scores.json"),
        },
        "limitations": [
            "Hashes cover deterministic query payloads and tracked parsed outputs, not raw API response dumps.",
            "L1000CDS2 db-version remains recorded as latest in the historical snapshot; a beta freeze should pin a concrete upstream version if available.",
            "This manifest intentionally does not re-call external APIs to avoid response drift.",
        ],
    }

    out_path = OUT_DIR / "api_snapshot_manifest.json"
    out_path.write_text(json.dumps(manifest, indent=2))
    print(json.dumps({"snapshot": str(out_path), "l1000_tissues": len(l1000_tissues), "crispr_tissues": len(crispr_tissues)}, indent=2))


if __name__ == "__main__":
    main()
