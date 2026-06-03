#!/usr/bin/env python3
"""Build the slides 1-14 editable PPTX skeleton with artifact-tool.

The script writes temporary slide modules inside a Presentations skill workspace,
then uses the bundled artifact-tool builder to export PPTX, rendered previews,
layout JSON, and a contact sheet. Final deliverables are copied into the
project's output directory for the SpaceBio-Bench deck lineage.
"""

from __future__ import annotations

import json
import os
import re
import shutil
import subprocess
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
CREATED = "2026-06-03"
THREAD_ID = os.environ.get("CODEX_THREAD_ID", "manual-20260603-slides-1-14-pptx-skeleton-v0-1")
TASK_SLUG = "spacebiobench-slides-1-14-pptx-skeleton"
WORKSPACE = ROOT / "outputs" / THREAD_ID / "presentations" / TASK_SLUG
SLIDES_DIR = WORKSPACE / "slides"
PREVIEW_DIR = WORKSPACE / "preview"
LAYOUT_DIR = WORKSPACE / "layout" / "final"
QA_DIR = WORKSPACE / "qa"
WORKSPACE_OUTPUT_DIR = WORKSPACE / "output"
PROJECT_OUTPUT_DIR = ROOT / "output" / "slides_1_14_pptx_skeleton_v0_1"

SKILL_DIR = Path(
    "/Users/jak4013/.codex/plugins/cache/openai-primary-runtime/presentations/26.601.10930/skills/presentations"
)
NODE = Path(
    os.environ.get(
        "PRESENTATIONS_NODE",
        "/Users/jak4013/.cache/codex-runtimes/codex-primary-runtime/dependencies/node/bin/node",
    )
)
BUILD_ARTIFACT_DECK = SKILL_DIR / "scripts" / "build_artifact_deck.mjs"

AUDIT_JSON = ROOT / "output" / "slides_1_14_pptx_readiness_audit_v0_1" / "slides_1_14_pptx_readiness_audit.json"
FINAL_PPTX = WORKSPACE_OUTPUT_DIR / "spacebiobench_slides_1_14_pptx_skeleton_v0_1.pptx"
CONTACT_SHEET = PREVIEW_DIR / "slides_1_14_pptx_skeleton_contact_sheet.png"
BUILD_MANIFEST = WORKSPACE_OUTPUT_DIR / "slides_1_14_pptx_skeleton_artifact_manifest.json"


COLORS = {
    "void": "#070A0E",
    "panel": "#0B1117",
    "panel2": "#101923",
    "ink": "#F4F7F8",
    "soft": "#B9C7D2",
    "muted": "#788898",
    "rule": "#33465A",
    "teal": "#1AA090",
    "sky": "#6BAED6",
    "amber": "#B4832F",
    "rose": "#B45A7E",
    "green": "#178B63",
    "violet": "#7B68A8",
}

SECTION_ACCENT = {
    "Opening": "teal",
    "Core methods": "sky",
    "v1-v7 results": "amber",
    "core result spine": "green",
}


def rel(path: Path) -> str:
    return str(path.relative_to(ROOT))


def clean_text(text: str) -> str:
    return re.sub(r"\s+", " ", text or "").strip()


def wrap_source(text: str, max_len: int = 92) -> str:
    text = clean_text(text)
    if len(text) <= max_len:
        return text
    return text[: max_len - 1].rstrip() + "."


def js(value) -> str:
    return json.dumps(value, ensure_ascii=True)


def choose_overlay(row: dict) -> dict:
    visible = [clean_text(item) for item in row.get("visible_text_sample", []) if clean_text(item)]
    section = row["section"]
    slide_number = row["slide"]

    if slide_number in {1, 2, 3, 4, 5}:
        eyebrow = visible[0] if visible else section.upper()
        title = visible[1] if len(visible) > 1 else row["headline"]
        subtitle = visible[2] if len(visible) > 2 else ""
        bridge = visible[3] if len(visible) > 3 else ""
        caveat = visible[4] if len(visible) > 4 else row.get("claim_boundary", "")
        source = visible[5] if len(visible) > 5 else row.get("editable_overlay_source", "")
        if slide_number == 1:
            title = "SpaceBio-Bench"
            subtitle = "Testing biological AI under spaceflight domain shift."
        if slide_number == 3:
            title = "From benchmark results to platform"
            subtitle = "v1-v9 separates completed evidence, hypothesis-only translation, and platformization."
        return {
            "eyebrow": eyebrow,
            "title": title,
            "subtitle": subtitle,
            "bridge": bridge,
            "caveat": caveat,
            "source": source,
            "layout": "opening" if slide_number <= 3 else "method",
        }

    if slide_number == 6:
        return {
            "eyebrow": "FEATURE LAYERS",
            "title": row["headline"],
            "subtitle": visible[1] if len(visible) > 1 else "Same samples, different feature views.",
            "bridge": "Define the feature view before reading score surfaces.",
            "caveat": row.get("claim_boundary", ""),
            "source": "feature-layer bridge overlay spec + v1-v7 benchmark core",
            "layout": "method",
        }

    if slide_number in {7, 8, 9}:
        return {
            "eyebrow": "EARLY RESULT FAMILY",
            "title": visible[0] if visible else row["headline"],
            "subtitle": visible[1] if len(visible) > 1 else row["headline"],
            "bridge": row.get("headline", ""),
            "caveat": row.get("claim_boundary", ""),
            "source": visible[3] if len(visible) > 3 else row.get("editable_overlay_source", ""),
            "layout": "result",
        }

    return {
        "eyebrow": "CORE RESULT SPINE",
        "title": visible[0] if visible else row["headline"],
        "subtitle": visible[1] if len(visible) > 1 else row["headline"],
        "bridge": "",
        "caveat": visible[2] if len(visible) > 2 else row.get("claim_boundary", ""),
        "source": "core result caption pack + slide-specific v4-v6 evaluation artifacts",
        "layout": "core",
    }


def speaker_notes(row: dict, overlay: dict) -> str:
    lines = [
        f"Slide {row['slide']}: {row['title']}",
        f"Claim: {overlay['title']}",
        f"PPTX instruction: {row.get('pptx_instruction', '')}",
        f"Claim boundary: {row.get('claim_boundary', '')}",
    ]
    if row.get("viewer_language_watch_hits"):
        lines.append("Viewer-language watch: " + ", ".join(row["viewer_language_watch_hits"]))
    if row.get("speaker_notes_required"):
        lines.append("Required presenter bridge:")
        lines.extend(f"- {item}" for item in row["speaker_notes_required"])
    return "\n".join(line for line in lines if clean_text(line))


def build_deck_spec() -> list[dict]:
    audit = json.loads(AUDIT_JSON.read_text(encoding="utf-8"))
    deck_spec: list[dict] = []
    for row in audit["slides"]:
        overlay = choose_overlay(row)
        accent = COLORS[SECTION_ACCENT.get(row["section"], "teal")]
        deck_spec.append(
            {
                "slide": row["slide"],
                "section": row["section"],
                "title": row["title"],
                "role": row["role"],
                "asset": str((ROOT / row["scene_plate"]).resolve()),
                "accent": accent,
                "overlay": overlay,
                "notes": speaker_notes(row, overlay),
                "watch": row.get("viewer_language_watch_hits", []),
            }
        )
    return deck_spec


def write_workspace_notes(deck_spec: list[dict]) -> None:
    WORKSPACE.mkdir(parents=True, exist_ok=True)
    QA_DIR.mkdir(parents=True, exist_ok=True)
    profile = [
        "task mode: create",
        "primary deck-profile: engineering-platform",
        "secondary profile gates: appendix-heavy source discipline; scientific figure claim-boundary discipline",
        "required proof objects: workflow/method bridges, feature-layer bridge, metric result proof scenes, translation-boundary slide",
        "source/asset requirements: use audited scene plates from local output; rebuild interpretation text as editable PPTX objects",
        "brand authenticity constraints: no external logos or invented identity marks",
        "profile-specific QA gates: preserve technical precision, keep claim attached to proof object, render contact sheet, inspect overlaps",
        "known missing inputs: no final institutional template supplied; this is a skeleton production draft",
    ]
    (WORKSPACE / "profile-plan.txt").write_text("\n".join(profile) + "\n", encoding="utf-8")

    source_lines = [
        "Source notes",
        f"Date: {CREATED}",
        "",
        "Primary local source: output/slides_1_14_pptx_readiness_audit_v0_1/slides_1_14_pptx_readiness_audit.json",
        "All slide proof surfaces use audited local scene plates. No external logo or identity asset is embedded.",
        "",
    ]
    for slide in deck_spec:
        source_lines.append(f"Slide {slide['slide']}: {rel(Path(slide['asset']))}")
    (WORKSPACE / "source-notes.txt").write_text("\n".join(source_lines) + "\n", encoding="utf-8")

    claim_lines = ["Claim spine", f"Date: {CREATED}", ""]
    for slide in deck_spec:
        claim_lines.append(f"{slide['slide']}. {slide['overlay']['title']} | {slide['overlay']['caveat']}")
    (WORKSPACE / "claim-spine.txt").write_text("\n".join(claim_lines) + "\n", encoding="utf-8")

    contact_lines = [
        "Contact-sheet plan",
        f"Date: {CREATED}",
        "",
        "Macro-layout rhythm:",
        "1-3 opening scenes with atmospheric proof surfaces",
        "4-6 methods bridge scenes with workflow proof surfaces",
        "7-9 early result proof scenes",
        "10-14 hardened core result proof scenes",
        "",
        "Hard rule: no rendered preview should become the final all-in-one slide; text shell remains editable.",
    ]
    (WORKSPACE / "contact-sheet-plan.txt").write_text("\n".join(contact_lines) + "\n", encoding="utf-8")


SHARED_MODULE = r'''
const COLORS = {
  void: "#070A0E",
  panel: "#0B1117",
  panel2: "#101923",
  ink: "#F4F7F8",
  soft: "#B9C7D2",
  muted: "#788898",
  rule: "#33465A",
};

function text(slide, ctx, opts) {
  return ctx.addText(slide, {
    ...opts,
    typeface: opts.typeface || "Aptos",
    insets: opts.insets || { left: 0, right: 0, top: 0, bottom: 0 },
  });
}

function rule(slide, ctx, x, y, w, color) {
  ctx.addShape(slide, {
    x,
    y,
    w,
    h: 2,
    geometry: "rect",
    fill: color,
    line: ctx.line("#00000000", 0),
  });
}

function matte(slide, ctx, y, h) {
  ctx.addShape(slide, {
    x: 0,
    y,
    w: ctx.W,
    h,
    geometry: "rect",
    fill: COLORS.void,
    line: ctx.line("#00000000", 0),
  });
}

function fontForTitle(title, layout) {
  if (layout === "opening") return title.length > 54 ? 35 : 42;
  if (layout === "core") return title.length > 42 ? 27 : 32;
  return title.length > 48 ? 29 : 34;
}

export async function addShellSlide(presentation, ctx, spec) {
  const slide = presentation.slides.add();
  await ctx.addImage(slide, {
    path: spec.asset,
    x: 0,
    y: 0,
    w: ctx.W,
    h: ctx.H,
    fit: "cover",
    alt: `${spec.title} audited scene plate`,
    name: `scene_plate_${spec.slide}`,
  });

  matte(slide, ctx, 0, spec.overlay.layout === "opening" ? 178 : spec.overlay.layout === "core" ? 154 : 170);
  matte(slide, ctx, 638, 82);

  rule(slide, ctx, 58, spec.overlay.layout === "opening" ? 158 : spec.overlay.layout === "core" ? 136 : 150, 388, spec.accent);
  text(slide, ctx, {
    text: `${String(spec.slide).padStart(2, "0")} / ${spec.overlay.eyebrow}`,
    x: 58,
    y: 30,
    w: 440,
    h: 22,
    fontSize: 11,
    bold: true,
    color: spec.accent,
    name: `eyebrow_${spec.slide}`,
  });

  text(slide, ctx, {
    text: spec.overlay.title,
    x: 58,
    y: spec.overlay.layout === "opening" ? 54 : spec.overlay.layout === "core" ? 52 : 58,
    w: spec.overlay.layout === "opening" ? 790 : 720,
    h: spec.overlay.layout === "opening" ? 56 : spec.overlay.layout === "core" ? 38 : 42,
    fontSize: fontForTitle(spec.overlay.title, spec.overlay.layout),
    bold: true,
    color: COLORS.ink,
    typeface: "Aptos Display",
    name: `title_${spec.slide}`,
  });

  text(slide, ctx, {
    text: spec.overlay.subtitle,
    x: 60,
    y: spec.overlay.layout === "opening" ? 121 : spec.overlay.layout === "core" ? 96 : 108,
    w: spec.overlay.layout === "opening" ? 740 : 690,
    h: spec.overlay.layout === "core" ? 38 : 40,
    fontSize: spec.overlay.layout === "core" ? 13 : 15,
    color: COLORS.soft,
    name: `subtitle_${spec.slide}`,
  });

  if (spec.overlay.bridge) {
    text(slide, ctx, {
      text: spec.overlay.bridge,
      x: 858,
      y: spec.overlay.layout === "opening" ? 52 : 45,
      w: 318,
      h: 58,
      fontSize: spec.overlay.layout === "opening" ? 17 : 15,
      color: COLORS.ink,
      bold: spec.overlay.layout === "opening",
      name: `bridge_${spec.slide}`,
    });
    rule(slide, ctx, 858, spec.overlay.layout === "opening" ? 114 : 104, 230, spec.accent);
  }

  text(slide, ctx, {
    text: spec.overlay.caveat,
    x: 58,
    y: 653,
    w: 792,
    h: 28,
    fontSize: 10.5,
    color: COLORS.soft,
    name: `claim_boundary_${spec.slide}`,
  });

  text(slide, ctx, {
    text: spec.overlay.source,
    x: 58,
    y: 681,
    w: 1048,
    h: 18,
    fontSize: 8.5,
    color: COLORS.muted,
    name: `source_${spec.slide}`,
  });

  if (slide.speakerNotes?.setText) {
    slide.speakerNotes.setText(spec.notes || "");
  }

  return slide;
}
'''


def write_slide_modules(deck_spec: list[dict]) -> None:
    if SLIDES_DIR.exists():
        shutil.rmtree(SLIDES_DIR)
    SLIDES_DIR.mkdir(parents=True, exist_ok=True)
    (SLIDES_DIR / "shared.mjs").write_text(SHARED_MODULE.strip() + "\n", encoding="utf-8")
    for slide in deck_spec:
        module = f"""import {{ addShellSlide }} from "./shared.mjs";

const spec = {js(slide)};

export async function slide{slide['slide']:02d}(presentation, ctx) {{
  return addShellSlide(presentation, ctx, spec);
}}
"""
        (SLIDES_DIR / f"slide-{slide['slide']:02d}.mjs").write_text(module, encoding="utf-8")


def run_build() -> dict:
    WORKSPACE_OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    PREVIEW_DIR.mkdir(parents=True, exist_ok=True)
    LAYOUT_DIR.mkdir(parents=True, exist_ok=True)
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
        "14",
        "--slide-size",
        "1280x720",
        "--scale",
        "1",
    ]
    result = subprocess.run(cmd, cwd=ROOT, text=True, capture_output=True, check=True)
    return json.loads(result.stdout)


def copy_final_deliverables(build_manifest: dict, deck_spec: list[dict]) -> dict:
    PROJECT_OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    copied = {
        "pptx": PROJECT_OUTPUT_DIR / "spacebiobench_slides_1_14_pptx_skeleton_v0_1.pptx",
        "contact_sheet": PROJECT_OUTPUT_DIR / "slides_1_14_pptx_skeleton_contact_sheet.png",
        "artifact_manifest": PROJECT_OUTPUT_DIR / "slides_1_14_pptx_skeleton_artifact_manifest.json",
        "speaker_notes": PROJECT_OUTPUT_DIR / "slides_1_14_pptx_skeleton_speaker_notes.md",
    }
    shutil.copy2(FINAL_PPTX, copied["pptx"])
    shutil.copy2(CONTACT_SHEET, copied["contact_sheet"])
    shutil.copy2(BUILD_MANIFEST, copied["artifact_manifest"])

    note_lines = ["# Slides 1-14 PPTX Skeleton Speaker Notes", "", f"Date: {CREATED}", ""]
    for slide in deck_spec:
        note_lines.extend([f"## Slide {slide['slide']}: {slide['title']}", "", slide["notes"], ""])
    copied["speaker_notes"].write_text("\n".join(note_lines), encoding="utf-8")

    production_manifest = {
        "created": CREATED,
        "artifact_role": "slides 1-14 editable PPTX skeleton",
        "workspace": rel(WORKSPACE),
        "presentation_skill": "artifact-tool presentation build",
        "decision_source": "output/slides_1_14_pptx_readiness_audit_v0_1/slides_1_14_pptx_readiness_audit.json",
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
            "editable_shell": "title, subtitle, bridge, source, claim boundary, and speaker notes are built as artifact-tool text objects",
            "scene_plate_rule": "audited scene plates are used instead of rendered previews",
            "known_status": "production skeleton, not final polished presentation export",
        },
    }
    production_path = PROJECT_OUTPUT_DIR / "slides_1_14_pptx_skeleton_manifest.json"
    production_path.write_text(json.dumps(production_manifest, indent=2), encoding="utf-8")
    copied["production_manifest"] = production_path
    return production_manifest


def main() -> None:
    deck_spec = build_deck_spec()
    write_workspace_notes(deck_spec)
    write_slide_modules(deck_spec)
    build_manifest = run_build()
    production_manifest = copy_final_deliverables(build_manifest, deck_spec)
    print(json.dumps(production_manifest, indent=2))


if __name__ == "__main__":
    main()
