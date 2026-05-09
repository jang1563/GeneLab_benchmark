"""Build a toxicity-aware triage table for INTERVENE candidates.

This script does not recommend clinical use. It joins the current multi-tissue
L1000CDS2 reversal matrix with a small, explicit class-level curation table so
countermeasure hypotheses are ranked with safety liabilities visible.
"""
from __future__ import annotations

import json
from pathlib import Path
from typing import Any

import pandas as pd

OUT_DIR = Path(__file__).resolve().parent / "evaluation"
REV_COLS = ["rev_thymus", "rev_gastrocnemius", "rev_skin", "rev_eye", "rev_liver", "rev_kidney"]
TISSUES = [c.replace("rev_", "") for c in REV_COLS]

CANONICAL = {
    "BRD-K79090631": "CGP-60474",
    "BRD-A45889380": "QUINACRINE HYDROCHLORIDE",
    "BRD-K72414522": "AZD-5438",
    "BRD-A19500257": "geldanamycin",
    "BRD-K13390322": "AT-7519",
    "BRD-K21680192": "mitoxantrone",
    "BRD-K54233340": "Dorsomorphin dihydrochloride",
    "BRD-K79131256": "ALBENDAZOLE",
}

CURATION = {
    "BRD-K79090631": {
        "target_class": "CDK1/2/5 inhibitor",
        "broad_mechanism": "Cell-cycle kinase inhibition",
        "known_toxicity_class": "anti-proliferative kinase liability",
        "safety_gate": "high: mechanism axis only; cell-cycle targets require marrow/gut/gonad safety work",
    },
    "BRD-A45889380": {
        "target_class": "Acridine antimalarial / DNA-interacting compound",
        "broad_mechanism": "Immune-modulating lysosomotropic/DNA-interacting stress response",
        "known_toxicity_class": "systemic antiprotozoal/acridine liability",
        "safety_gate": "moderate-high: approved-use history does not transfer to spaceflight countermeasure use",
    },
    "BRD-K72414522": {
        "target_class": "CDK1/2/9 inhibitor",
        "broad_mechanism": "Cell-cycle and transcriptional kinase inhibition",
        "known_toxicity_class": "anti-proliferative oncology-like kinase liability",
        "safety_gate": "high: broad CDK inhibition is a target-axis flag, not a deployable candidate",
    },
    "BRD-A19500257": {
        "target_class": "HSP90 inhibitor",
        "broad_mechanism": "Heat-shock client protein destabilization",
        "known_toxicity_class": "proteostasis stress / oncology-like liability",
        "safety_gate": "high: HSP90 perturbation needs independent toxicity-aware validation",
    },
    "BRD-K13390322": {
        "target_class": "CDK1/2/4/6/9 inhibitor",
        "broad_mechanism": "Multi-CDK cell-cycle inhibition",
        "known_toxicity_class": "anti-proliferative kinase liability",
        "safety_gate": "high: CRISPR support increases target plausibility, not safety",
    },
    "BRD-K21680192": {
        "target_class": "TOP2 inhibitor / anthracenedione",
        "broad_mechanism": "DNA topoisomerase II poisoning",
        "known_toxicity_class": "cytotoxic DNA-damage liability",
        "safety_gate": "very high: cytotoxic mechanism is incompatible with direct countermeasure language",
    },
    "BRD-K54233340": {
        "target_class": "AMPK/BMP pathway tool inhibitor",
        "broad_mechanism": "Kinase-pathway modulation around AMPK/BMP signaling",
        "known_toxicity_class": "tool-compound/off-target kinase uncertainty",
        "safety_gate": "high uncertainty: use as pathway pointer only",
    },
    "BRD-K79131256": {
        "target_class": "Microtubule/anthelmintic",
        "broad_mechanism": "Tubulin and proliferation-axis perturbation",
        "known_toxicity_class": "anthelmintic systemic-use liability",
        "safety_gate": "moderate-high: single-tissue signature hit; requires dose, liver, and marrow checks",
    },
}


def _load_scores() -> pd.DataFrame:
    df = pd.read_csv(OUT_DIR / "multi_tissue_drug_matrix.csv")
    df["candidate"] = df["pert_id"].map(CANONICAL).fillna(df["perturbation"])
    # Merge synonym rows with the same perturbagen id by taking the best observed
    # reversal score per tissue from the tracked L1000CDS2 snapshot.
    grouped = df.groupby(["pert_id", "candidate"], as_index=False)[REV_COLS].max(numeric_only=True)
    grouped["n_tissues"] = grouped[REV_COLS].notna().sum(axis=1)
    grouped["mean_rev"] = grouped[REV_COLS].mean(axis=1)
    grouped["min_rev"] = grouped[REV_COLS].min(axis=1)
    grouped["max_rev"] = grouped[REV_COLS].max(axis=1)
    grouped["tissues"] = grouped.apply(
        lambda row: ";".join(t for t, c in zip(TISSUES, REV_COLS) if pd.notna(row[c])),
        axis=1,
    )
    return grouped


def _load_pareto_ids() -> set[str]:
    p = OUT_DIR / "pareto_front.csv"
    if not p.exists():
        return set()
    return set(pd.read_csv(p)["pert_id"].astype(str))


def _load_crispr_support() -> dict[str, str]:
    p = OUT_DIR / "crispr_orthog_summary.json"
    if not p.exists():
        return {}
    summary = json.loads(p.read_text())
    support: dict[str, list[str]] = {}
    for tissue, payload in summary.get("tissues", {}).items():
        for drug, hits in payload.get("validated_drugs", {}).items():
            for hit in hits:
                target = hit.get("target", "target")
                adj_p = hit.get("adj_p")
                detail = f"{tissue}:{target} KO {hit.get('reverses', 'reversal')}"
                if isinstance(adj_p, (int, float)):
                    detail += f" adj_p={adj_p:.3g}"
                support.setdefault(drug.upper(), []).append(detail)
    return {k: "; ".join(v) for k, v in support.items()}


def main() -> None:
    scores = _load_scores()
    pareto_ids = _load_pareto_ids()
    crispr = _load_crispr_support()

    rows: list[dict[str, Any]] = []
    for pert_id, curation in CURATION.items():
        match = scores[scores["pert_id"] == pert_id]
        if match.empty:
            continue
        record = match.iloc[0].to_dict()
        candidate = str(record["candidate"])
        row = {
            "pert_id": pert_id,
            "candidate": candidate,
            "target_class": curation["target_class"],
            "known_toxicity_class": curation["known_toxicity_class"],
            "broad_mechanism": curation["broad_mechanism"],
            "n_tissues": int(record["n_tissues"]),
            "tissues": record["tissues"],
            "mean_rev": round(float(record["mean_rev"]), 6),
            "min_rev": round(float(record["min_rev"]), 6),
            "max_rev": round(float(record["max_rev"]), 6),
            "pareto_front": pert_id in pareto_ids,
            "orthogonal_support": crispr.get(candidate.upper(), "none in current CRISPR-KO screen"),
            "safety_gate": curation["safety_gate"],
            "recommendation_scope": "hypothesis-generating target/pathway triage only",
        }
        rows.append(row)

    out = pd.DataFrame(rows).sort_values(
        ["pareto_front", "n_tissues", "mean_rev"],
        ascending=[False, False, False],
    )
    out.to_csv(OUT_DIR / "safety_triage.csv", index=False)

    summary = {
        "n_candidates": int(len(out)),
        "n_pareto_front_candidates": int(out["pareto_front"].sum()),
        "n_with_crispr_orthogonal_support": int((out["orthogonal_support"] != "none in current CRISPR-KO screen").sum()),
        "high_level_interpretation": (
            "Current multi-tissue reversers are dominated by CDK, HSP90, TOP2, "
            "and tool-compound mechanisms; these remain pathway hypotheses, not "
            "candidate countermeasures."
        ),
    }
    (OUT_DIR / "safety_triage_summary.json").write_text(json.dumps(summary, indent=2))
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
