#!/usr/bin/env python3
"""Audit regenerated premium visual values against source-of-truth files."""

from __future__ import annotations

import csv
import importlib.util
import json
import re
from pathlib import Path
from typing import Any


ROOT = Path(__file__).resolve().parents[1]
OUT_DIR = ROOT / "output" / "premium_audit"
SCRIPT_PATH = ROOT / "scripts" / "build_premium_visuals.py"


def load_visual_module() -> Any:
    spec = importlib.util.spec_from_file_location("build_premium_visuals", SCRIPT_PATH)
    if spec is None or spec.loader is None:
        raise RuntimeError(f"Could not load {SCRIPT_PATH}")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


bp = load_visual_module()


def clean_cell(value: str) -> str:
    value = value.strip()
    value = value.replace("**", "")
    value = value.replace("`", "")
    value = value.replace("<br>", " ")
    return value.strip()


def first_float(value: str) -> float:
    match = re.search(r"[-+]?\d+(?:\.\d+)?", value)
    if not match:
        raise ValueError(f"No numeric value in {value!r}")
    return float(match.group(0))


def ci_pair(value: str) -> tuple[float, float]:
    nums = re.findall(r"[-+]?\d+(?:\.\d+)?", value)
    if len(nums) < 2:
        raise ValueError(f"No CI pair in {value!r}")
    return float(nums[0]), float(nums[1])


def md_table_after(path: Path, heading: str) -> list[dict[str, str]]:
    lines = path.read_text(encoding="utf-8").splitlines()
    start = None
    for idx, line in enumerate(lines):
        if heading in line:
            start = idx + 1
            break
    if start is None:
        raise ValueError(f"Heading not found: {heading}")
    table_lines: list[str] = []
    in_table = False
    for line in lines[start:]:
        if line.strip().startswith("|"):
            in_table = True
            table_lines.append(line)
        elif in_table:
            break
    if len(table_lines) < 2:
        raise ValueError(f"No markdown table after heading: {heading}")
    headers = [clean_cell(cell) for cell in table_lines[0].strip("|").split("|")]
    rows: list[dict[str, str]] = []
    for line in table_lines[2:]:
        cells = [clean_cell(cell) for cell in line.strip("|").split("|")]
        if len(cells) != len(headers):
            continue
        rows.append(dict(zip(headers, cells)))
    return rows


def csv_rows(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as fh:
        return list(csv.DictReader(fh))


def json_data(path: Path) -> Any:
    with path.open(encoding="utf-8") as fh:
        return json.load(fh)


def add_check(
    checks: list[dict[str, Any]],
    figure: str,
    item: str,
    field: str,
    expected: Any,
    observed: Any,
    source: str,
    tolerance: float = 5e-4,
) -> None:
    numeric = isinstance(expected, (int, float)) and isinstance(observed, (int, float))
    if numeric:
        diff = abs(float(expected) - float(observed))
        status = "PASS" if diff <= tolerance else "FAIL"
    else:
        diff = "" if str(expected) == str(observed) else "mismatch"
        status = "PASS" if str(expected) == str(observed) else "FAIL"
    checks.append(
        {
            "figure": figure,
            "item": item,
            "field": field,
            "expected": expected,
            "observed": observed,
            "difference": diff,
            "tolerance": tolerance if numeric else "",
            "status": status,
            "source": source,
        }
    )


def audit_mouse_core(checks: list[dict[str, Any]]) -> None:
    results_path = ROOT / "evaluation" / "RESULTS_SUMMARY.md"
    category_b = md_table_after(results_path, "## Category B: Cross-Mission Transfer")
    category_b_by_tissue = {row["Tissue"].lower(): row for row in category_b}
    for row in bp.CATEGORY_B_TRANSFER:
        tissue = row["tissue"]
        observed = category_b_by_tissue[tissue]
        ci_low, ci_high = ci_pair(observed["95% CI"])
        for field, obs in [
            ("auroc", first_float(observed["Mean AUROC"])),
            ("ci_low", ci_low),
            ("ci_high", ci_high),
            ("tier", int(first_float(observed["Tier"]))),
        ]:
            add_check(checks, "fig1_tissue_transfer", tissue, field, row[field], obs, row["source"])


def audit_pathway(checks: list[dict[str, Any]]) -> None:
    results_path = ROOT / "evaluation" / "RESULTS_SUMMARY.md"
    category_d = md_table_after(results_path, "## Category D: Condition/Confounder Prediction")
    category_d_lookup = {}
    for row in category_d:
        task = row["Task"].split()[0]
        tissue = row["Tissue"].lower()
        category_d_lookup[(task, tissue)] = row
    for row in bp.CONFOUNDER_FEATURE_COMPARISON:
        task, tissue = str(row["task"]).split("\n")
        observed = category_d_lookup[(str(row["note"]), tissue)]
        add_check(
            checks,
            "fig2_pathway",
            f"{row['note']}:{tissue}:gene",
            "macro_f1",
            row["gene"],
            first_float(observed["Gene F1"]),
            row["source"],
        )
        add_check(
            checks,
            "fig2_pathway",
            f"{row['note']}:{tissue}:pathway",
            "macro_f1",
            row["pathway"],
            first_float(observed["Pathway F1"]),
            row["source"],
        )

    j5 = md_table_after(results_path, "### Category A — Spaceflight Detection")
    j5_by_tissue = {row["Tissue"].lower(): row for row in j5}
    for row in bp.GENE_PATHWAY_DETECTION:
        tissue = row["tissue"]
        observed = j5_by_tissue[tissue]
        add_check(checks, "fig2_pathway", tissue, "gene_auroc", row["gene"], first_float(observed["Gene"]), row["source"])
        add_check(
            checks,
            "fig2_pathway",
            tissue,
            "pathway_auroc",
            row["pathway"],
            first_float(observed["Pathway"]),
            row["source"],
        )

    nes = md_table_after(results_path, "## NES Conservation vs Cross-Mission Transfer")
    nes_by_tissue = {row["Tissue"].lower(): row for row in nes}
    for row in bp.NES_TRANSFER:
        tissue = row["tissue"]
        observed = nes_by_tissue[tissue]
        add_check(checks, "fig2_pathway", tissue, "nes_r", row["nes_r"], first_float(observed["NES Mean r"]), row["source"])
        add_check(
            checks,
            "fig2_pathway",
            tissue,
            "transfer_auroc",
            row["transfer_auroc"],
            first_float(observed["Transfer AUROC"]),
            row["source"],
        )


def audit_models(checks: list[dict[str, Any]]) -> None:
    canonical = md_table_after(ROOT / "docs" / "CANONICAL_RESULTS_V7_1.md", "## Canonical FM / LLM Snapshot")
    canonical_by_model = {row["Model"].replace(" baseline", ""): row for row in canonical}
    tier3 = md_table_after(ROOT / "evaluation" / "RESULTS_SUMMARY.md", "## Tier 3: LLM Zero-Shot Classification")
    tier3_by_model = {row["Model"].replace("(ref)", "").strip(): row for row in tier3}
    for row in bp.MODEL_TIER_MEANS:
        model = row["model"]
        if model in canonical_by_model:
            observed = first_float(canonical_by_model[model]["AUROC / score"])
            source = row["source"]
        else:
            observed = first_float(tier3_by_model[model]["Mean"])
            source = "evaluation/RESULTS_SUMMARY.md:Tier 3 LLM Zero-Shot Classification"
        add_check(checks, "fig3_model", model, "score", row["score"], observed, source)

    tier2 = md_table_after(ROOT / "evaluation" / "RESULTS_SUMMARY.md", "## Tier 2: scGPT")
    tier2_rows = [row for row in tier2 if row["Tissue"].lower() not in {"6 tissues"}]
    tier2_by_tissue = {row["Tissue"].lower(): row for row in tier2_rows}
    for row in bp.MODEL_TIER_DELTAS:
        tissue = row["tissue"]
        observed_row = tier2_by_tissue[tissue]
        if row["model"] == "scGPT":
            observed = first_float(observed_row["Δ vs Baseline"])
        else:
            observed = first_float(observed_row["Geneformer AUROC"]) - first_float(observed_row["Baseline AUROC"])
        add_check(checks, "fig3_model", f"{row['model']}:{tissue}", "delta_vs_baseline", row["delta"], observed, "evaluation/RESULTS_SUMMARY.md:Tier 2")


def audit_v9_tables(checks: list[dict[str, Any]]) -> None:
    gap_summary = csv_rows(ROOT / "v9" / "reports" / "public_bulk_alpha_gap_matrix" / "public_bulk_alpha_gap_summary.csv")[0]
    source_inventory = csv_rows(ROOT / "v9" / "source_inventory.csv")
    organoid_source_inventory = csv_rows(ROOT / "v9" / "human_organoid" / "source_inventory.draft.csv")
    organoid_samples = csv_rows(ROOT / "v9" / "human_organoid" / "sample_factors.draft.csv")
    de_manifest = json_data(ROOT / "v9" / "human_organoid" / "de_references" / "human_organoid_de_reference_manifest.draft.json")
    multispecies_samples = csv_rows(ROOT / "v9" / "multispecies" / "sample_factors.draft.csv")
    sc_assets = csv_rows(ROOT / "v9" / "sc_spaceflight" / "asset_inventory_summary.csv")[0]
    sc_audit = csv_rows(ROOT / "v9" / "sc_spaceflight" / "obs_var_audit_summary.csv")[0]

    observed_values = {
        ("Bulk", "public source rows"): len(source_inventory),
        ("Bulk", "held-out mission folds"): int(gap_summary["bulk_fold_count"]),
        ("Bulk", "baseline rows"): int(gap_summary["baseline_row_count"]),
        ("Organoid", "public samples"): len(organoid_samples),
        ("Organoid", "gene contrasts"): int(de_manifest["totals"]["n_contrasts"]),
        ("Organoid", "gene reference rows"): int(de_manifest["totals"]["n_rows"]),
        ("Multispecies", "public samples"): len(multispecies_samples),
        ("Single-cell", "legacy assets indexed"): int(sc_assets["total_asset_count"]),
        ("Single-cell", "data blockers"): int(sc_audit["blocker_count"]),
    }
    for row in bp.V9_FOOTPRINT_ROWS:
        key = (row["group"], row["metric"])
        add_check(checks, "fig4_v9_platform", ":".join(key), "value", row["value"], observed_values[key], row["source"])

    decision = csv_rows(ROOT / "v9" / "reports" / "public_bulk_alpha_snapshot_decision" / "snapshot_decision_summary.csv")[0]
    summary_fields = {
        "bulk_task_count": 8,
        "bulk_fold_count": 33,
        "bulk_source_count": 22,
        "checksum_parsed_source_count": 22,
        "freeze_ready_source_count": 0,
        "baseline_row_count": 24,
    }
    for field, expected in summary_fields.items():
        add_check(checks, "fig5_public_bulk_boundary", "snapshot_decision_summary", field, expected, int(decision[field]), "v9/reports/public_bulk_alpha_snapshot_decision/snapshot_decision_summary.csv")


def audit_organoid(checks: list[dict[str, Any]]) -> None:
    source_inventory = csv_rows(ROOT / "v9" / "human_organoid" / "source_inventory.draft.csv")
    samples = csv_rows(ROOT / "v9" / "human_organoid" / "sample_factors.draft.csv")
    de_manifest = json_data(ROOT / "v9" / "human_organoid" / "de_references" / "human_organoid_de_reference_manifest.draft.json")
    nearest = csv_rows(ROOT / "v9" / "human_organoid" / "reports" / "nearest_centroid" / "human_organoid_baseline_summary.csv")[0]
    source_transfer = json_data(ROOT / "v9" / "human_organoid" / "reports" / "source_transfer_signature" / "metrics.json")["metrics"]
    feature_effect = json_data(ROOT / "v9" / "human_organoid" / "reports" / "logistic_feature_effect" / "metrics.json")["metrics"]

    footprint_values = {
        "public sources": len(source_inventory),
        "samples": len(samples),
        "flight-ground gene contrasts": int(de_manifest["totals"]["n_contrasts"]),
        "gene reference rows": int(de_manifest["totals"]["n_rows"]),
        "significant gene rows": int(de_manifest["totals"]["n_significant_fdr_0_05"]),
        "pilot status": "draft",
    }
    for row in bp.ORGANOID_FOOTPRINT_ROWS:
        add_check(checks, "fig6_organoid", row["metric"], "value", row["display"] if row["metric"] == "pilot status" else row["value"], footprint_values[row["metric"]], row["source"])

    metric_values = {
        ("Primary task", "AUROC"): float(nearest["auroc"]),
        ("Primary task", "Macro-F1"): float(nearest["macro_f1"]),
        ("Flight-response pattern", "Direction match"): source_transfer["de_direction_match"]["value"],
        ("Flight-response pattern", "Rank correlation"): source_transfer["signature_rank_correlation"]["value"],
        ("Model gene effects", "Direction match"): feature_effect["feature_effect_direction_match"]["value"],
        ("Model gene effects", "Rank correlation"): feature_effect["feature_effect_rank_correlation"]["value"],
    }
    for row in bp.ORGANOID_METRIC_ROWS:
        key = (row["family"], row["metric"])
        add_check(checks, "fig6_organoid", ":".join(key), "value", row["value"], metric_values[key], row["source"])

    topk = feature_effect["feature_effect_top_k_de_overlap"]["value"]
    for row in bp.ORGANOID_TOPK_ROWS:
        observed = topk[str(row["k"])]
        add_check(checks, "fig6_organoid", f"top_{row['k']}", "observed", row["observed"], observed["n_overlap"], row["source"])
        add_check(checks, "fig6_organoid", f"top_{row['k']}", "expected", row["expected"], observed["expected_overlap"], row["source"])
        add_check(checks, "fig6_organoid", f"top_{row['k']}", "enrichment", row["enrichment"], observed["enrichment"], row["source"])
        add_check(checks, "fig6_organoid", f"top_{row['k']}", "p_value", row["p_value"], observed["hypergeometric_p_value_greater_equal"], row["source"])


def write_outputs(checks: list[dict[str, Any]]) -> None:
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    csv_path = OUT_DIR / "premium_visual_numeric_audit.csv"
    json_path = OUT_DIR / "premium_visual_numeric_audit.json"
    md_path = ROOT / "docs" / "VISUAL_NUMERIC_SOURCE_AUDIT_2026_05_31.md"
    keys = ["figure", "item", "field", "expected", "observed", "difference", "tolerance", "status", "source"]
    with csv_path.open("w", newline="", encoding="utf-8") as fh:
        writer = csv.DictWriter(fh, fieldnames=keys)
        writer.writeheader()
        writer.writerows(checks)
    with json_path.open("w", encoding="utf-8") as fh:
        json.dump(checks, fh, indent=2)
        fh.write("\n")
    total = len(checks)
    failed = [row for row in checks if row["status"] != "PASS"]
    by_figure: dict[str, list[dict[str, Any]]] = {}
    for row in checks:
        by_figure.setdefault(row["figure"], []).append(row)
    lines = [
        "# Visual Numeric Source Audit",
        "",
        "Date: 2026-05-31",
        "",
        "Purpose: verify that regenerated premium figure/table values match their cited source files.",
        "",
        f"Total checks: {total}",
        f"Passing checks: {total - len(failed)}",
        f"Failing checks: {len(failed)}",
        "",
        "Generated audit files:",
        "",
        "- `output/premium_audit/premium_visual_numeric_audit.csv`",
        "- `output/premium_audit/premium_visual_numeric_audit.json`",
        "",
        "## Figure Summary",
        "",
        "| Figure | Checks | Fails | Status |",
        "|---|---:|---:|---|",
    ]
    for figure, rows in sorted(by_figure.items()):
        fails = sum(row["status"] != "PASS" for row in rows)
        lines.append(f"| {figure} | {len(rows)} | {fails} | {'PASS' if fails == 0 else 'FAIL'} |")
    lines.extend(["", "## Failing Checks", ""])
    if failed:
        lines.append("| Figure | Item | Field | Expected | Observed | Source |")
        lines.append("|---|---|---|---:|---:|---|")
        for row in failed:
            lines.append(f"| {row['figure']} | {row['item']} | {row['field']} | {row['expected']} | {row['observed']} | {row['source']} |")
    else:
        lines.append("No failing numeric checks.")
    lines.extend(
        [
            "",
            "## Interpretation",
            "",
            "- A PASS means the scripted figure value matches the cited source within the audit tolerance.",
            "- This audit does not validate whether a source file is itself final or publication-frozen.",
            "- Source freeze and caption wording remain separate L4 requirements.",
            "",
        ]
    )
    md_path.write_text("\n".join(lines), encoding="utf-8")


def main() -> None:
    checks: list[dict[str, Any]] = []
    audit_mouse_core(checks)
    audit_pathway(checks)
    audit_models(checks)
    audit_v9_tables(checks)
    audit_organoid(checks)
    write_outputs(checks)
    fail_count = sum(row["status"] != "PASS" for row in checks)
    print(json.dumps({"checks": len(checks), "failures": fail_count, "output": "output/premium_audit/premium_visual_numeric_audit.csv"}, indent=2))
    if fail_count:
        raise SystemExit(1)


if __name__ == "__main__":
    main()
