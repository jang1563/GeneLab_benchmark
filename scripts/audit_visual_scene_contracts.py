#!/usr/bin/env python3
"""Audit premium visual scene contracts before or after rendering."""

from __future__ import annotations

import argparse
import csv
import json
import re
from pathlib import Path
from typing import Any


ROOT = Path(__file__).resolve().parents[1]
DEFAULT_ROOT = ROOT / "output" / "premium_bridge_scenes"
SAFE_AREA = {"x_min": 0.055, "x_max": 0.945, "y_min": 0.075, "y_max": 0.900}
REQUIRED_MANIFEST_FIELDS = [
    "slide_id",
    "created",
    "stage",
    "content_brief",
    "evidence_sources",
    "claim_boundary",
    "generator",
    "outputs",
]
REQUIRED_OVERLAY_FIELDS = [
    "slide_id",
    "canvas",
    "coordinate_system",
    "text",
    "status_labels",
    "focus_marks",
    "forbidden_visible_terms",
]


def load_json(path: Path) -> Any:
    with path.open(encoding="utf-8") as fh:
        return json.load(fh)


def rel(path: Path) -> str:
    return str(path.relative_to(ROOT))


def resolve_repo_path(value: str) -> Path:
    path = Path(value)
    if path.is_absolute():
        return path
    return ROOT / path


def add_check(
    rows: list[dict[str, Any]],
    slide_id: str,
    check: str,
    status: bool,
    detail: str,
    severity: str = "error",
) -> None:
    rows.append(
        {
            "slide_id": slide_id,
            "check": check,
            "status": "PASS" if status else "FAIL",
            "severity": severity,
            "detail": detail,
        }
    )


def visible_text(overlay: dict[str, Any]) -> str:
    items = list(overlay.get("text", [])) + list(overlay.get("status_labels", []))
    return " ".join(str(item.get("content", "")) for item in items)


def count_words(text: str) -> int:
    return len(re.findall(r"[A-Za-z0-9]+(?:[-'][A-Za-z0-9]+)?", text))


def check_canvas(rows: list[dict[str, Any]], slide_id: str, overlay: dict[str, Any]) -> None:
    canvas = overlay.get("canvas", {})
    add_check(rows, slide_id, "canvas_width", canvas.get("width_px") == 3840, str(canvas.get("width_px")))
    add_check(rows, slide_id, "canvas_height", canvas.get("height_px") == 2160, str(canvas.get("height_px")))
    add_check(rows, slide_id, "canvas_aspect", canvas.get("aspect_ratio") == "16:9", str(canvas.get("aspect_ratio")))
    add_check(
        rows,
        slide_id,
        "coordinate_system",
        overlay.get("coordinate_system") == "normalized_0_1",
        str(overlay.get("coordinate_system")),
    )


def check_required_fields(rows: list[dict[str, Any]], slide_id: str, data: dict[str, Any], fields: list[str], label: str) -> None:
    for field in fields:
        add_check(rows, slide_id, f"{label}_field_{field}", field in data and data[field] not in (None, "", []), repr(data.get(field)))


def check_paths(rows: list[dict[str, Any]], slide_id: str, manifest: dict[str, Any], qa: dict[str, Any]) -> None:
    for field in ["content_brief", "technical_preflight", "generator"]:
        if field in manifest:
            path = resolve_repo_path(str(manifest[field]))
            add_check(rows, slide_id, f"path_exists_{field}", path.exists(), str(path))

    for idx, source in enumerate(manifest.get("evidence_sources", [])):
        path_value = source.get("path")
        if not path_value:
            add_check(rows, slide_id, f"evidence_source_{idx}_path_declared", False, repr(source))
            continue
        path = resolve_repo_path(str(path_value))
        add_check(rows, slide_id, f"evidence_source_{idx}_exists", path.exists(), str(path))

    outputs = manifest.get("outputs", {})
    for key in ["scene_plate", "rendered_preview_png", "rendered_preview_pdf", "overlay_spec", "manifest", "qa"]:
        add_check(rows, slide_id, f"output_path_declared_{key}", bool(outputs.get(key)), repr(outputs.get(key)))

    stage = str(qa.get("stage", manifest.get("stage", "")))
    if stage == "post_render":
        for key in ["scene_plate", "rendered_preview_png", "rendered_preview_pdf"]:
            path_value = outputs.get(key)
            path = resolve_repo_path(str(path_value)) if path_value else ROOT / "__missing__"
            add_check(rows, slide_id, f"render_output_exists_{key}", path.exists(), str(path))
    else:
        add_check(rows, slide_id, "render_outputs_allowed_pending", True, f"stage={stage}", severity="info")


def check_text(rows: list[dict[str, Any]], slide_id: str, overlay: dict[str, Any]) -> None:
    text = visible_text(overlay)
    words = count_words(text)
    budget = int(overlay.get("visible_word_budget", 45))
    add_check(rows, slide_id, "visible_word_budget", words <= budget, f"{words}/{budget}")

    for term in overlay.get("forbidden_visible_terms", []):
        term_text = str(term)
        if not term_text:
            continue
        present = term_text.lower() in text.lower()
        add_check(rows, slide_id, f"forbidden_term_absent_{term_text}", not present, term_text)

    add_check(rows, slide_id, "no_absolute_local_path_in_visible_text", "/Users/" not in text, text)


def check_overlay_items(rows: list[dict[str, Any]], slide_id: str, overlay: dict[str, Any]) -> None:
    for group in ["text", "status_labels"]:
        for item in overlay.get(group, []):
            item_id = item.get("id", "<missing>")
            x = item.get("x")
            y = item.get("y")
            font = item.get("font_pt")
            in_bounds = isinstance(x, (int, float)) and isinstance(y, (int, float)) and 0 <= x <= 1 and 0 <= y <= 1
            add_check(rows, slide_id, f"{group}_{item_id}_coord_normalized", in_bounds, f"x={x}, y={y}")
            safe = (
                isinstance(x, (int, float))
                and isinstance(y, (int, float))
                and SAFE_AREA["x_min"] <= x <= SAFE_AREA["x_max"]
                and SAFE_AREA["y_min"] <= y <= SAFE_AREA["y_max"]
            )
            add_check(rows, slide_id, f"{group}_{item_id}_coord_safe_area", safe, f"x={x}, y={y}", severity="warn")
            font_ok = isinstance(font, (int, float)) and font >= 6.5
            add_check(rows, slide_id, f"{group}_{item_id}_font_min", font_ok, f"font_pt={font}")
            line_len = max((len(line) for line in str(item.get("content", "")).splitlines()), default=0)
            add_check(rows, slide_id, f"{group}_{item_id}_line_length", line_len <= 58, f"max_line_chars={line_len}", severity="warn")

    for item in overlay.get("focus_marks", []):
        item_id = item.get("id", "<missing>")
        coords = []
        for key in ["x", "y", "x0", "x1", "y0", "y1"]:
            if key in item:
                coords.append((key, item[key]))
        normalized = all(isinstance(value, (int, float)) and 0 <= value <= 1 for _, value in coords)
        add_check(rows, slide_id, f"focus_{item_id}_coords_normalized", normalized, repr(coords))


def audit_slide_dir(slide_dir: Path) -> list[dict[str, Any]]:
    manifest_path = slide_dir / "manifest.json"
    overlay_path = slide_dir / "overlay_spec.json"
    qa_path = slide_dir / "qa.json"
    slide_id = slide_dir.name
    rows: list[dict[str, Any]] = []

    for label, path in [("manifest", manifest_path), ("overlay_spec", overlay_path), ("qa", qa_path)]:
        add_check(rows, slide_id, f"{label}_exists", path.exists(), str(path))
    if not manifest_path.exists() or not overlay_path.exists() or not qa_path.exists():
        return rows

    manifest = load_json(manifest_path)
    overlay = load_json(overlay_path)
    qa = load_json(qa_path)
    slide_id = str(manifest.get("slide_id", slide_id))

    check_required_fields(rows, slide_id, manifest, REQUIRED_MANIFEST_FIELDS, "manifest")
    check_required_fields(rows, slide_id, overlay, REQUIRED_OVERLAY_FIELDS, "overlay")

    add_check(rows, slide_id, "slide_id_match_overlay", manifest.get("slide_id") == overlay.get("slide_id"), f"{manifest.get('slide_id')} vs {overlay.get('slide_id')}")
    add_check(rows, slide_id, "slide_id_match_qa", manifest.get("slide_id") == qa.get("slide_id"), f"{manifest.get('slide_id')} vs {qa.get('slide_id')}")
    add_check(rows, slide_id, "stage_match_qa", manifest.get("stage") == qa.get("stage"), f"{manifest.get('stage')} vs {qa.get('stage')}")

    check_canvas(rows, slide_id, overlay)
    check_paths(rows, slide_id, manifest, qa)
    check_text(rows, slide_id, overlay)
    check_overlay_items(rows, slide_id, overlay)

    pre_gate = qa.get("pre_render_gate", {})
    for key, value in pre_gate.items():
        if isinstance(value, bool):
            add_check(rows, slide_id, f"pre_render_gate_{key}", value, str(value))
    return rows


def discover_slide_dirs(root: Path) -> list[Path]:
    if (root / "manifest.json").exists():
        return [root]
    return sorted(path for path in root.iterdir() if path.is_dir() and (path / "manifest.json").exists())


def write_reports(rows: list[dict[str, Any]], output_root: Path) -> dict[str, Path]:
    output_root.mkdir(parents=True, exist_ok=True)
    json_path = output_root / "visual_scene_contract_audit.json"
    csv_path = output_root / "visual_scene_contract_audit.csv"
    json_path.write_text(json.dumps(rows, indent=2) + "\n", encoding="utf-8")
    with csv_path.open("w", encoding="utf-8", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=["slide_id", "check", "status", "severity", "detail"])
        writer.writeheader()
        writer.writerows(rows)
    return {"json": json_path, "csv": csv_path}


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", default=str(DEFAULT_ROOT), help="Scene contract root or one slide directory.")
    parser.add_argument("--output-root", default=str(DEFAULT_ROOT), help="Directory for audit reports.")
    args = parser.parse_args()

    root = resolve_repo_path(args.root)
    rows: list[dict[str, Any]] = []
    for slide_dir in discover_slide_dirs(root):
        rows.extend(audit_slide_dir(slide_dir))

    output_root = resolve_repo_path(args.output_root)
    reports = write_reports(rows, output_root)
    failures = [row for row in rows if row["status"] == "FAIL" and row["severity"] == "error"]
    warnings = [row for row in rows if row["status"] == "FAIL" and row["severity"] == "warn"]
    print(
        json.dumps(
            {
                "checks": len(rows),
                "error_failures": len(failures),
                "warnings": len(warnings),
                "json": rel(reports["json"]),
                "csv": rel(reports["csv"]),
            },
            indent=2,
        )
    )
    if failures:
        raise SystemExit(1)


if __name__ == "__main__":
    main()
