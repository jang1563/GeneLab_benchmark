#!/usr/bin/env python3
"""Assemble the SpaceBio-Bench detailed deck into a PPTX.

The detailed-deck slides are already designed and QA'd as 3840x2160 PNG
surfaces. This builder preserves those audited surfaces in an export-ready PPTX
by placing each image full-bleed on a 16:9 artifact-tool slide, then rendering
previews and a contact sheet through the standard Presentations skill builder.
"""

from __future__ import annotations

import csv
import json
import os
import shutil
import subprocess
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
CREATED = "2026-06-13"
THREAD_ID = os.environ.get("CODEX_THREAD_ID", "019e8bd1-3b42-7821-9e1c-beebfa8e2ece")
TASK_SLUG = "spacebiobench-detailed-deck"
WORKSPACE = ROOT / "outputs" / THREAD_ID / "presentations" / TASK_SLUG
SLIDES_DIR = WORKSPACE / "slides"
PREVIEW_DIR = WORKSPACE / "preview" / "pptx-final"
LAYOUT_DIR = WORKSPACE / "layout" / "pptx-final"
QA_DIR = WORKSPACE / "qa" / "pptx-final"
WORKSPACE_OUTPUT_DIR = WORKSPACE / "output"
PROJECT_OUTPUT_DIR = ROOT / "output" / "spacebiobench_detailed_deck_v0_1"

SKILL_DIR = Path(
    "<PATH>"
)
NODE = Path(
    os.environ.get(
        "PRESENTATIONS_NODE",
        "<PATH>",
    )
)
BUILD_ARTIFACT_DECK = SKILL_DIR / "scripts" / "build_artifact_deck.mjs"

ASSEMBLY_PLAN = WORKSPACE / "spacebiobench_detailed_deck_assembly_plan_v2.tsv"
FINAL_PPTX = WORKSPACE_OUTPUT_DIR / "spacebiobench_detailed_deck_v0_1.pptx"
CONTACT_SHEET = PREVIEW_DIR / "spacebiobench_detailed_deck_v0_1_pptx_contact_sheet.png"
BUILD_MANIFEST = WORKSPACE_OUTPUT_DIR / "spacebiobench_detailed_deck_v0_1_artifact_manifest.json"


SHARED_MODULE = r'''
export async function addImageSlide(presentation, ctx, spec) {
  const slide = presentation.slides.add();
  ctx.addShape(slide, {
    x: 0,
    y: 0,
    w: ctx.W,
    h: ctx.H,
    geometry: "rect",
    fill: "#070A0E",
    line: ctx.line("#00000000", 0),
  });
  await ctx.addImage(slide, {
    path: spec.asset,
    x: 0,
    y: 0,
    w: ctx.W,
    h: ctx.H,
    fit: "contain",
    alt: `${spec.slide}: ${spec.title}`,
    name: `spacebiobench_detailed_slide_${String(spec.slide).padStart(2, "0")}`,
  });
  if (slide.speakerNotes?.setText) {
    slide.speakerNotes.setText(spec.notes || "");
  }
  return slide;
}
'''


def rel(path: Path) -> str:
    return str(path.relative_to(ROOT))


def js(value: object) -> str:
    return json.dumps(value, ensure_ascii=True)


def load_deck_spec() -> list[dict[str, object]]:
    if not ASSEMBLY_PLAN.exists():
        raise FileNotFoundError(f"Missing assembly plan: {ASSEMBLY_PLAN}")
    with ASSEMBLY_PLAN.open(newline="", encoding="utf-8") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))

    deck_spec: list[dict[str, object]] = []
    for row in rows:
        slide = int(row["slide"])
        asset = ROOT / row["asset_or_action"]
        if row["status"] != "ready":
            raise ValueError(f"Slide {slide} is not ready: {row['status']}")
        if not asset.exists():
            raise FileNotFoundError(f"Slide {slide} asset missing: {asset}")
        deck_spec.append(
            {
                "slide": slide,
                "act": row["act"],
                "title": row["title"],
                "main_question": row["main_question"],
                "visual_object": row["proof_object"],
                "asset": str(asset),
                "notes": (
                    f"{row['title']}\n"
                    f"Act: {row['act']}\n"
                    f"Question: {row['main_question']}\n"
                    f"Visual: {row['proof_object']}"
                ),
            }
        )

    expected = list(range(1, len(deck_spec) + 1))
    observed = [int(item["slide"]) for item in deck_spec]
    if observed != expected:
        raise ValueError(f"Slide numbers are not contiguous: {observed[:5]} ... {observed[-5:]}")
    return deck_spec


def write_workspace_notes(deck_spec: list[dict[str, object]]) -> None:
    WORKSPACE.mkdir(parents=True, exist_ok=True)
    profile_lines = [
        "Task mode: create / final PPTX assembly from audited slide surfaces",
        "Primary deck-profile: appendix-heavy scientific methods/results deck",
        "Secondary gates: engineering-platform clarity; scientific pedagogy for ML/foundation-model methods",
        "Output rule: preserve the audited 3840x2160 slide PNGs as full-bleed PPTX slides.",
        "QA gate: artifact-tool preview contact sheet plus detailed-deck full walkthrough OCR/visual QA.",
    ]
    (WORKSPACE / "profile-plan.txt").write_text("\n".join(profile_lines) + "\n", encoding="utf-8")

    claim_lines = ["# SpaceBio-Bench Detailed Deck PPTX Assembly Spine", "", f"Date: {CREATED}", ""]
    for item in deck_spec:
        claim_lines.append(
            f"{item['slide']:02d}. {item['title']} | {item['act']} | {item['main_question']}"
        )
    (WORKSPACE / "claim-spine.txt").write_text("\n".join(claim_lines) + "\n", encoding="utf-8")

    contact_lines = [
        "Contact-sheet plan",
        "",
        "1-12 methods setup and audience teaching scaffolds",
        "13-21 core result readout",
        "22-29 model-family and foundation-model comparison",
        "30-33 robustness and core benchmark synthesis",
        "34-42 biology and translation interpretation",
        "43-49 v8 incubator modules",
        "50-58 platform and extension readiness",
        "59-60 close and release synthesis",
    ]
    (WORKSPACE / "contact-sheet-plan.txt").write_text("\n".join(contact_lines) + "\n", encoding="utf-8")


def write_slide_modules(deck_spec: list[dict[str, object]]) -> None:
    if SLIDES_DIR.exists():
        shutil.rmtree(SLIDES_DIR)
    SLIDES_DIR.mkdir(parents=True, exist_ok=True)
    (SLIDES_DIR / "shared.mjs").write_text(SHARED_MODULE.strip() + "\n", encoding="utf-8")
    for item in deck_spec:
        module = f"""import {{ addImageSlide }} from "./shared.mjs";

const spec = {js(item)};

export async function slide{int(item['slide']):02d}(presentation, ctx) {{
  return addImageSlide(presentation, ctx, spec);
}}
"""
        (SLIDES_DIR / f"slide-{int(item['slide']):02d}.mjs").write_text(module, encoding="utf-8")


def run_build(deck_spec: list[dict[str, object]]) -> dict[str, object]:
    WORKSPACE_OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    PREVIEW_DIR.mkdir(parents=True, exist_ok=True)
    LAYOUT_DIR.mkdir(parents=True, exist_ok=True)
    QA_DIR.mkdir(parents=True, exist_ok=True)
    cmd = [
        str(NODE),
        str(BUILD_ARTIFACT_DECK),
        "--workspace",
        str(WORKSPACE),
        "--slides-dir",
        str(SLIDES_DIR),
        "--out",
        str(FINAL_PPTX),
        "--preview-dir",
        str(PREVIEW_DIR),
        "--layout-dir",
        str(LAYOUT_DIR),
        "--contact-sheet",
        str(CONTACT_SHEET),
        "--manifest",
        str(BUILD_MANIFEST),
        "--slide-count",
        str(len(deck_spec)),
        "--slide-size",
        "1280x720",
        "--scale",
        "1",
    ]
    result = subprocess.run(cmd, cwd=ROOT, text=True, capture_output=True, check=True)
    return json.loads(result.stdout)


def copy_final_deliverables(build_manifest: dict[str, object], deck_spec: list[dict[str, object]]) -> dict[str, object]:
    PROJECT_OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    copied = {
        "pptx": PROJECT_OUTPUT_DIR / "spacebiobench_detailed_deck_v0_1.pptx",
        "contact_sheet": PROJECT_OUTPUT_DIR / "spacebiobench_detailed_deck_v0_1_pptx_contact_sheet.png",
        "artifact_manifest": PROJECT_OUTPUT_DIR / "spacebiobench_detailed_deck_v0_1_artifact_manifest.json",
        "speaker_notes": PROJECT_OUTPUT_DIR / "spacebiobench_detailed_deck_v0_1_speaker_notes.md",
    }
    shutil.copy2(FINAL_PPTX, copied["pptx"])
    shutil.copy2(CONTACT_SHEET, copied["contact_sheet"])
    shutil.copy2(BUILD_MANIFEST, copied["artifact_manifest"])

    note_lines = ["# SpaceBio-Bench Detailed Deck Speaker Notes", "", f"Date: {CREATED}", ""]
    for item in deck_spec:
        note_lines.extend(
            [
                f"## Slide {item['slide']}: {item['title']}",
                "",
                f"Act: {item['act']}",
                f"Question: {item['main_question']}",
                f"Visual: {item['visual_object']}",
                "",
            ]
        )
    copied["speaker_notes"].write_text("\n".join(note_lines), encoding="utf-8")

    production_manifest = {
        "created": CREATED,
        "artifact_role": "60-slide SpaceBio-Bench detailed deck assembled from audited 3840x2160 slide surfaces",
        "workspace": rel(WORKSPACE),
        "presentation_skill": "artifact-tool presentation build",
        "outputs": {key: rel(path) for key, path in copied.items()},
        "workspace_outputs": {
            "pptx": rel(FINAL_PPTX),
            "contact_sheet": rel(CONTACT_SHEET),
            "artifact_manifest": rel(BUILD_MANIFEST),
            "preview_dir": rel(PREVIEW_DIR),
            "layout_dir": rel(LAYOUT_DIR),
        },
        "qa_summary": {
            "slide_count": build_manifest.get("slideCount"),
            "pptx_bytes": build_manifest.get("outputBytes"),
            "assembly_plan": rel(ASSEMBLY_PLAN),
            "source_surfaces": "Each slide is inserted full-bleed from its ready PNG in the assembly plan.",
        },
    }
    production_path = PROJECT_OUTPUT_DIR / "spacebiobench_detailed_deck_v0_1_manifest.json"
    production_path.write_text(json.dumps(production_manifest, indent=2) + "\n", encoding="utf-8")
    copied["production_manifest"] = production_path
    return production_manifest


def main() -> None:
    deck_spec = load_deck_spec()
    write_workspace_notes(deck_spec)
    write_slide_modules(deck_spec)
    build_manifest = run_build(deck_spec)
    production_manifest = copy_final_deliverables(build_manifest, deck_spec)
    print(json.dumps(production_manifest, indent=2))


if __name__ == "__main__":
    main()
