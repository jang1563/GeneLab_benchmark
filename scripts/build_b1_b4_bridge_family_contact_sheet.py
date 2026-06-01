#!/usr/bin/env python3
"""Build a family QA contact sheet for premium B1-B4 bridge candidates."""

from __future__ import annotations

import json
import os
from pathlib import Path
from typing import Any


ROOT = Path(__file__).resolve().parents[1]
os.environ.setdefault("MPLCONFIGDIR", str(ROOT / "output" / ".matplotlib"))

import matplotlib

matplotlib.use("Agg")
import matplotlib.image as mpimg
import matplotlib.pyplot as plt

try:
    from PIL import Image
except ImportError:  # pragma: no cover - local QA fallback only
    Image = None  # type: ignore[assignment]


OUT_DIR = ROOT / "output" / "premium_bridge_family_review"
CREATED = "2026-06-01"
CONTACT_SHEET = OUT_DIR / "b1_b4_premium_family_contact_sheet.png"
MANIFEST_PATH = OUT_DIR / "b1_b4_premium_family_contact_sheet_manifest.json"
QA_PATH = OUT_DIR / "b1_b4_premium_family_contact_sheet_qa.json"

SLIDES = [
    {
        "slide_id": "b1_evaluation_layer",
        "panel": "B1",
        "title": "Evaluation layer",
        "claim": "Public studies become evaluable.",
        "source": "output/premium_bridge_scenes/b1_evaluation_layer/rendered_preview.png",
        "qa": "output/premium_bridge_scenes/b1_evaluation_layer/qa.json",
        "overlay": "output/premium_bridge_scenes/b1_evaluation_layer/overlay_spec.json",
    },
    {
        "slide_id": "b2_study_to_task_premium",
        "panel": "B2",
        "title": "Benchmark task",
        "claim": "Source context stays attached.",
        "source": "output/premium_bridge_rebuild_scenes/b2_study_to_task_premium/rendered_preview.png",
        "qa": "output/premium_bridge_rebuild_scenes/b2_study_to_task_premium/qa.json",
        "overlay": "output/premium_bridge_rebuild_scenes/b2_study_to_task_premium/overlay_spec.json",
    },
    {
        "slide_id": "b3_mission_held_out_premium",
        "panel": "B3",
        "title": "Held-out mission",
        "claim": "The test set is an entire mission.",
        "source": "output/premium_bridge_rebuild_scenes/b3_mission_held_out_premium/rendered_preview.png",
        "qa": "output/premium_bridge_rebuild_scenes/b3_mission_held_out_premium/qa.json",
        "overlay": "output/premium_bridge_rebuild_scenes/b3_mission_held_out_premium/overlay_spec.json",
    },
    {
        "slide_id": "b4_train_only_guard_premium",
        "panel": "B4",
        "title": "Train-only guard",
        "claim": "Training choices stop before the hidden mission.",
        "source": "output/premium_bridge_rebuild_scenes/b4_train_only_guard_premium/rendered_preview.png",
        "qa": "output/premium_bridge_rebuild_scenes/b4_train_only_guard_premium/qa.json",
        "overlay": "output/premium_bridge_rebuild_scenes/b4_train_only_guard_premium/overlay_spec.json",
    },
]


def rel(path: Path) -> str:
    return str(path.relative_to(ROOT))


def write_json(path: Path, data: Any) -> None:
    path.write_text(json.dumps(data, indent=2) + "\n", encoding="utf-8")


def read_json(path: Path) -> dict[str, Any]:
    return json.loads(path.read_text(encoding="utf-8"))


def image_dimensions(path: Path) -> list[int] | None:
    if Image is None:
        return None
    with Image.open(path) as img:
        return [int(img.width), int(img.height)]


def slide_record(slide: dict[str, str]) -> dict[str, Any]:
    source = ROOT / slide["source"]
    qa_path = ROOT / slide["qa"]
    overlay_path = ROOT / slide["overlay"]
    qa = read_json(qa_path) if qa_path.exists() else {}
    overlay = read_json(overlay_path) if overlay_path.exists() else {}
    return {
        **slide,
        "source_exists": source.exists(),
        "qa_exists": qa_path.exists(),
        "overlay_exists": overlay_path.exists(),
        "image_dimensions": image_dimensions(source) if source.exists() else None,
        "visible_word_count": overlay.get("visible_word_count"),
        "visual_verdict": qa.get("post_render_gate", {}).get("visual_verdict"),
    }


def build_contact_sheet(records: list[dict[str, Any]]) -> None:
    fig = plt.figure(figsize=(18.5, 11.0), dpi=220)
    grid = fig.add_gridspec(2, 2, left=0.028, right=0.985, bottom=0.060, top=0.835, wspace=0.034, hspace=0.175)
    fig.suptitle(
        "Premium bridge family QA: B1-B4",
        x=0.028,
        y=0.955,
        ha="left",
        fontsize=13.5,
        fontweight="bold",
        color="#17212B",
    )
    fig.text(
        0.028,
        0.905,
        "Review target: one evidence surface per slide, consistent rail grammar, readable claim, quiet source/status layer.",
        fontsize=8.2,
        color="#5D6978",
        ha="left",
    )
    fig.text(
        0.028,
        0.872,
        "This is a QA sheet, not deck art.",
        fontsize=7.5,
        color="#5D6978",
        ha="left",
    )
    for idx, record in enumerate(records):
        ax = fig.add_subplot(grid[idx // 2, idx % 2])
        ax.imshow(mpimg.imread(ROOT / record["source"]))
        title = f"{record['panel']} | {record['title']} | {record['claim']}"
        ax.set_title(title, loc="left", fontsize=7.6, pad=4, color="#17212B")
        ax.axis("off")
    fig.savefig(CONTACT_SHEET, dpi=220, facecolor="white")
    plt.close(fig)


def build() -> dict[str, Any]:
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    records = [slide_record(slide) for slide in SLIDES]
    build_contact_sheet(records)
    manifest = {
        "created": CREATED,
        "contact_sheet": rel(CONTACT_SHEET),
        "purpose": "family QA for premium B1-B4 bridge candidates",
        "slides": records,
        "generator": "scripts/build_b1_b4_bridge_family_contact_sheet.py",
    }
    qa = {
        "created": CREATED,
        "contact_sheet": rel(CONTACT_SHEET),
        "all_sources_exist": all(record["source_exists"] for record in records),
        "all_overlay_specs_exist": all(record["overlay_exists"] for record in records),
        "all_qa_files_exist": all(record["qa_exists"] for record in records),
        "source_dimensions": {record["slide_id"]: record["image_dimensions"] for record in records},
        "visible_word_counts": {record["slide_id"]: record["visible_word_count"] for record in records},
        "manual_family_review": {
            "claim_readability": "pending",
            "rail_grammar_consistency": "pending",
            "evidence_surface_scale": "pending",
            "status_source_layer_consistency": "pending",
            "promotion_verdict": "pending",
        },
    }
    write_json(MANIFEST_PATH, manifest)
    write_json(QA_PATH, qa)
    return {"manifest": manifest, "qa": qa}


def main() -> None:
    print(json.dumps(build(), indent=2))


if __name__ == "__main__":
    main()
