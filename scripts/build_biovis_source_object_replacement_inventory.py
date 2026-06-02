#!/usr/bin/env python3
"""Build a source-object replacement inventory for v0.4 proof modules.

This is not a final slide generator. It records which real source objects can
replace v0.4 placeholders before a production deck is assembled.
"""

from __future__ import annotations

import csv
import json
from pathlib import Path

from PIL import Image, ImageDraw, ImageFont


ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "output" / "biovis_source_object_replacement_inventory_v0_1"
CREATED = "2026-06-02"

COLORS = {
    "deep": "#0D1720",
    "deep2": "#111D27",
    "deep3": "#172636",
    "white": "#F2F6F8",
    "soft": "#B8C4CF",
    "muted": "#6D7D8D",
    "blue": "#2D6F9F",
    "sky": "#6BAED6",
    "green": "#178B63",
    "teal": "#1AA090",
    "amber": "#B4832F",
    "rose": "#B45A7E",
    "red": "#B23A3A",
}


def rgb(token: str) -> tuple[int, int, int]:
    value = COLORS[token].lstrip("#")
    return tuple(int(value[i : i + 2], 16) for i in (0, 2, 4))


def rgba(token: str, alpha: int) -> tuple[int, int, int, int]:
    return rgb(token) + (alpha,)


def rel(path: Path) -> str:
    return str(path.relative_to(ROOT))


def load_json(path: str) -> object:
    with (ROOT / path).open() as handle:
        return json.load(handle)


def font(size: int, *, bold: bool = False):
    candidates = [
        "/System/Library/Fonts/Supplemental/Arial Bold.ttf" if bold else "/System/Library/Fonts/Supplemental/Arial.ttf",
        "/System/Library/Fonts/Helvetica.ttc",
    ]
    for candidate in candidates:
        try:
            return ImageFont.truetype(candidate, size=size)
        except OSError:
            continue
    return ImageFont.load_default()


F = {
    "title": font(48, bold=True),
    "subtitle": font(23),
    "h": font(25, bold=True),
    "body": font(19),
    "small": font(15),
    "tiny": font(12),
}


def compact_sha(value: str) -> str:
    if not value or value == "NA":
        return value
    return f"{value[:10]}...{value[-8:]}"


def write_wrapped(draw: ImageDraw.ImageDraw, xy: tuple[int, int], text: str, font_obj, fill, width: int, line_height: int, max_lines: int) -> int:
    words = text.split()
    lines: list[str] = []
    current = ""
    for word in words:
        candidate = f"{current} {word}".strip()
        if draw.textlength(candidate, font=font_obj) <= width:
            current = candidate
        else:
            if current:
                lines.append(current)
            current = word
    if current:
        lines.append(current)
    if len(lines) > max_lines:
        lines = lines[:max_lines]
        lines[-1] = lines[-1].rstrip(".") + "..."
    x, y = xy
    for idx, line in enumerate(lines):
        draw.text((x, y + idx * line_height), line, font=font_obj, fill=fill)
    return y + len(lines) * line_height


def build_records() -> list[dict[str, str]]:
    human = load_json("v9/human_organoid/task_manifests/draft_human_organoid_spaceflight.json")
    human_matrix = human["split"]["expression_matrix_sources"]
    human_metrics = load_json("v9/human_organoid/reports/source_transfer_signature/metrics.json")
    multi37 = load_json("v9/multispecies/task_manifests/draft_osd37_arabidopsis_seedling_spaceflight.json")
    multi207 = load_json("v9/multispecies/task_manifests/draft_osd207_drosophila_whole_body_spaceflight.json")
    multi120 = load_json("v9/multispecies/interaction_task_manifests/draft_osd120_arabidopsis_root_light_interaction_spaceflight.json")
    sc = load_json("v9/sc_spaceflight/task_manifests/draft_rrrm1_blood_single_cell_spaceflight.json")
    sc_payload = load_json("v9/sc_spaceflight/obs_var_audit_summary.json")[0]

    human_matrix_line = "; ".join(
        [
            f"{sid}: {entry['n_feature_rows']} genes x {entry['n_sample_columns']} samples, sha256 {compact_sha(entry['sha256'])}"
            for sid, entry in human_matrix.items()
        ]
    )
    organoid_folds = human["split"]["candidate_folds"]
    organoid_fold_line = "; ".join(
        [f"{fold['heldout_factor']}={fold['heldout_value']} train {fold['n_train']} test {fold['n_test']}" for fold in organoid_folds[:4]]
    )
    organoid_result_line = (
        f"source-transfer signature: AUROC {human_metrics['metrics']['auroc']['value']:.3f}, "
        f"balanced accuracy {human_metrics['metrics']['balanced_accuracy']['value']:.3f}, "
        f"DE direction match {human_metrics['metrics']['de_direction_match']['value']:.3f}; diagnostic only"
    )

    multi_matrix_line = "; ".join(
        [
            f"OSD-37: {multi37['split']['expression_matrix_sources']['OSD-37']['n_feature_rows']} genes x {multi37['split']['expression_matrix_sources']['OSD-37']['n_sample_columns']} samples",
            f"OSD-207: {multi207['split']['expression_matrix_sources']['OSD-207']['n_feature_rows']} genes x {multi207['split']['expression_matrix_sources']['OSD-207']['n_sample_columns']} samples",
            f"OSD-120: {multi120['split']['expression_matrix_sources']['OSD-120']['n_feature_rows']} genes x {multi120['split']['expression_matrix_sources']['OSD-120']['n_sample_columns']} samples",
        ]
    )
    osd120_primary = "; ".join(
        [
            f"{fold['heldout_value']} train {fold['n_train']} test {fold['n_test']}"
            for fold in multi120["split"]["primary_candidate_folds"]
        ]
    )
    sc_record = sc["source_records"][0]
    sc_line = (
        f"{sc['sample_qc']['n_cells']} cells, {sc['sample_qc']['n_genes']} genes, "
        f"{sc['split']['sample_count']} animal samples; payload status {sc_payload['payload_path_status']}"
    )

    return [
        {
            "proof_id": "organoid_osdr_source_records",
            "priority": "P0",
            "lane": "human_organoid",
            "source_ids": "OSD-863;OSD-871",
            "module": "source_dataset_record_plate",
            "replacement_object": "official OSDR record screenshots for OSD-863 and OSD-871",
            "local_evidence": "v9/human_organoid/source_inventory.draft.json; v9/human_organoid/task_manifests/draft_human_organoid_spaceflight.json",
            "source_url": "https://osdr.nasa.gov/bio/repo/data/studies/OSD-863; https://osdr.nasa.gov/bio/repo/data/studies/OSD-871",
            "evidence_summary": "OSD-863/871 are public human iPSC neural organoid RNA-seq sources; OSD-863 lists GLDS-716, DOI, GSE259421, and related OSD-871.",
            "claim_boundary": "source record proof only; does not by itself validate benchmark performance",
            "next_capture_action": "capture OSDR header and metadata tabs; include source date and URL footer",
            "recommended_slide_role": "main bridge source proof or appendix source inventory",
        },
        {
            "proof_id": "organoid_expression_matrix_qc",
            "priority": "P0",
            "lane": "human_organoid",
            "source_ids": "OSD-863;OSD-871",
            "module": "expression_matrix_proof_plate",
            "replacement_object": "matrix audit table plus compact matrix/QC crop",
            "local_evidence": "v9/human_organoid/expression_matrix_audit.draft.json",
            "source_url": "https://osdr.nasa.gov/geode-py/ws/studies/OSD-863/download?source=datamanager&file=GLDS-716_rna_seq_Normalized_Counts_GLbulkRNAseq.csv; https://osdr.nasa.gov/geode-py/ws/studies/OSD-871/download?source=datamanager&file=GLDS-720_rna_seq_Normalized_Counts_GLbulkRNAseq.csv",
            "evidence_summary": human_matrix_line,
            "claim_boundary": "matrix downloaded and sample-aligned; source payload hash/freeze remains draft",
            "next_capture_action": "render matrix audit proof crop with sample counts, sha256, and matching sample columns",
            "recommended_slide_role": "appendix proof or methods bridge",
        },
        {
            "proof_id": "organoid_blocked_split",
            "priority": "P0",
            "lane": "human_organoid",
            "source_ids": "OSD-863;OSD-871",
            "module": "held_out_task_proof_plate",
            "replacement_object": "blocked fold manifest visual for organoid type and microglia condition",
            "local_evidence": "v9/human_organoid/task_manifests/draft_human_organoid_spaceflight.json",
            "source_url": "https://osdr.nasa.gov/bio/repo/data/studies/OSD-863",
            "evidence_summary": organoid_fold_line,
            "claim_boundary": "single ISS mission pilot; not a mission-held-out benchmark",
            "next_capture_action": "draw fold lanes directly from manifest and label draft/sample-count-backed status",
            "recommended_slide_role": "methods explainer",
        },
        {
            "proof_id": "organoid_diagnostic_result_claim",
            "priority": "P1",
            "lane": "human_organoid",
            "source_ids": "OSD-863;OSD-871",
            "module": "result_claim_plate",
            "replacement_object": "diagnostic result surface with explicit non-leaderboard caveat",
            "local_evidence": "v9/human_organoid/reports/source_transfer_signature/metrics.json",
            "source_url": "v9/human_organoid/reports/source_transfer_signature/metrics.json",
            "evidence_summary": organoid_result_line,
            "claim_boundary": "diagnostic-only; no frozen leaderboard or broad organoid benchmark claim",
            "next_capture_action": "show metrics only if the caveat footer is visible",
            "recommended_slide_role": "appendix result proof",
        },
        {
            "proof_id": "multispecies_source_matrix_inventory",
            "priority": "P0",
            "lane": "multispecies",
            "source_ids": "OSD-37;OSD-207;OSD-120",
            "module": "source_dataset_record_plate;expression_matrix_proof_plate",
            "replacement_object": "OSDR record screenshot plus matrix audit rows for plant/fly extension",
            "local_evidence": "v9/multispecies/source_inventory.draft.json; v9/multispecies/expression_matrix_audit.draft.json",
            "source_url": "https://osdr.nasa.gov/bio/repo/data/studies/OSD-37; https://osdr.nasa.gov/bio/repo/data/studies/OSD-207; https://osdr.nasa.gov/bio/repo/data/studies/OSD-120",
            "evidence_summary": multi_matrix_line,
            "claim_boundary": "species-native within-source pilots; do not present as raw-gene cross-species benchmark",
            "next_capture_action": "capture OSDR source headers and local matrix audit rows",
            "recommended_slide_role": "appendix inventory or species-extension bridge",
        },
        {
            "proof_id": "osd120_interaction_split",
            "priority": "P0",
            "lane": "multispecies",
            "source_ids": "OSD-120",
            "module": "held_out_task_proof_plate",
            "replacement_object": "genotype-by-light split manifest visual",
            "local_evidence": "v9/multispecies/interaction_task_manifests/draft_osd120_arabidopsis_root_light_interaction_spaceflight.json",
            "source_url": "https://osdr.nasa.gov/bio/repo/data/studies/OSD-120",
            "evidence_summary": osd120_primary,
            "claim_boundary": "within-source interaction diagnostic; not leave-one-mission-out",
            "next_capture_action": "draw primary genotype holdouts and secondary light-treatment holdouts as separate lanes",
            "recommended_slide_role": "main methods explainer candidate",
        },
        {
            "proof_id": "rrrm1_blood_sc_source_contract",
            "priority": "P1",
            "lane": "single_cell",
            "source_ids": "OSD-918",
            "module": "source_dataset_record_plate;single_cell_embedding_proof_plate",
            "replacement_object": "OSDR source record plus AnnData contract blocker panel",
            "local_evidence": "v9/sc_spaceflight/task_manifests/draft_rrrm1_blood_single_cell_spaceflight.json",
            "source_url": sc_record["source_url"],
            "evidence_summary": sc_line,
            "claim_boundary": sc["claim_boundary"],
            "next_capture_action": "capture OSDR-918 source record; do not render embedding until local h5ad exists",
            "recommended_slide_role": "future single-cell lane proof with blocker badge",
        },
    ]


def write_records(records: list[dict[str, str]]) -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    (OUT / "source_object_replacement_inventory.json").write_text(json.dumps({"created": CREATED, "records": records}, indent=2) + "\n")
    with (OUT / "source_object_replacement_inventory.csv").open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(records[0].keys()))
        writer.writeheader()
        writer.writerows(records)


def draw_badge(draw: ImageDraw.ImageDraw, xy: tuple[int, int], text: str, tone: str) -> None:
    x, y = xy
    width = max(88, int(draw.textlength(text, font=F["small"])) + 38)
    draw.rounded_rectangle((x, y, x + width, y + 30), radius=12, fill=rgba("deep2", 245), outline=rgb(tone), width=2)
    draw.ellipse((x + 12, y + 11, x + 20, y + 19), fill=rgb(tone))
    draw.text((x + 28, y + 8), text, font=F["small"], fill=rgb("white"))


def render_contact_sheet(records: list[dict[str, str]]) -> None:
    width, height = 3200, 1900
    canvas = Image.new("RGBA", (width, height), rgb("deep") + (255,))
    draw = ImageDraw.Draw(canvas)
    draw.text((80, 58), "Source-object replacement inventory v0.1", font=F["title"], fill=rgb("white"))
    draw.text(
        (80, 120),
        "Real source objects to replace v0.4 placeholders before production slide assembly.",
        font=F["subtitle"],
        fill=rgb("soft"),
    )
    draw.line((80, 176, width - 80, 176), fill=rgba("deep3", 255), width=2)
    for x in [560, 1180, 1980, 2600]:
        draw.line((x, 200, x, height - 90), fill=rgba("deep3", 150), width=1)

    y = 230
    row_h = 210
    lane_tones = {"human_organoid": "sky", "multispecies": "green", "single_cell": "rose"}
    for idx, record in enumerate(records):
        tone = lane_tones.get(record["lane"], "amber")
        band_alpha = 48 if idx % 2 else 28
        draw.rounded_rectangle((70, y - 16, width - 70, y + row_h - 28), radius=18, fill=rgba("deep2", band_alpha), outline=rgba("deep3", 190), width=1)
        draw.rectangle((88, y - 3, 98, y + row_h - 45), fill=rgb(tone))
        write_wrapped(draw, (120, y), record["proof_id"], F["h"], rgb("white"), 400, 29, 2)
        draw_badge(draw, (120, y + 38), record["priority"], tone)
        draw_badge(draw, (220, y + 38), record["lane"], tone)
        draw.text((120, y + 86), record["source_ids"], font=F["body"], fill=rgb("soft"))

        write_wrapped(draw, (640, y), record["replacement_object"], F["body"], rgb("white"), 470, 25, 4)
        write_wrapped(draw, (1260, y), record["evidence_summary"], F["small"], rgb("soft"), 650, 21, 5)
        write_wrapped(draw, (2050, y), record["claim_boundary"], F["small"], rgb("soft"), 500, 21, 4)
        write_wrapped(draw, (2670, y), record["module"], F["small"], rgb("white"), 380, 21, 4)
        y += row_h

    draw.text((80, height - 64), "Use rule: source screenshots and local proof crops go first; interpretation text is added only after status/caveat remains visible.", font=F["small"], fill=rgb("muted"))
    canvas.convert("RGB").save(OUT / "source_object_replacement_contact_sheet.png")


def write_qa(records: list[dict[str, str]]) -> None:
    qa = {
        "created": CREATED,
        "record_count": len(records),
        "automatic_checks": {
            "json_exists": (OUT / "source_object_replacement_inventory.json").exists(),
            "csv_exists": (OUT / "source_object_replacement_inventory.csv").exists(),
            "contact_sheet_exists": (OUT / "source_object_replacement_contact_sheet.png").exists(),
            "all_records_have_claim_boundary": all(bool(record["claim_boundary"]) for record in records),
            "all_records_have_next_capture_action": all(bool(record["next_capture_action"]) for record in records),
        },
        "manual_review": {
            "initial_verdict": "inventory ready for source-object replacement sprint",
            "main_slide_candidates": ["organoid_osdr_source_records", "osd120_interaction_split"],
            "appendix_candidates": [
                "organoid_expression_matrix_qc",
                "organoid_diagnostic_result_claim",
                "multispecies_source_matrix_inventory",
            ],
            "blocked_candidate": "rrrm1_blood_sc_source_contract remains source-record only until local AnnData payload exists",
        },
    }
    (OUT / "source_object_replacement_inventory_qa.json").write_text(json.dumps(qa, indent=2) + "\n")


def main() -> None:
    records = build_records()
    write_records(records)
    render_contact_sheet(records)
    write_qa(records)
    print(json.dumps({"output": rel(OUT), "records": len(records), "contact_sheet": True}, indent=2))


if __name__ == "__main__":
    main()
