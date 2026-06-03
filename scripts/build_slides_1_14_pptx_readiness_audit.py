#!/usr/bin/env python3
"""Build a PPTX production-readiness audit for slides 1-14."""

from __future__ import annotations

import csv
import json
import re
from pathlib import Path
from typing import Any

try:
    from PIL import Image, ImageDraw, ImageFilter, ImageFont
except ImportError:  # pragma: no cover - text outputs remain useful without Pillow.
    Image = None
    ImageDraw = None
    ImageFilter = None
    ImageFont = None


ROOT = Path(__file__).resolve().parents[1]
CREATED = "2026-06-03"
OUT = ROOT / "output" / "slides_1_14_pptx_readiness_audit_v0_1"
AUDIT_QA = OUT / "qa"

ASSEMBLY_MANIFEST = "output/slides_1_14_deck_assembly_bridge_v0_1/slides_1_14_deck_assembly_manifest.json"
ASSEMBLY_QA = "output/slides_1_14_deck_assembly_bridge_v0_1/slides_1_14_deck_assembly_qa.json"
ASSEMBLY_OVERLAY_RULES = "output/slides_1_14_deck_assembly_bridge_v0_1/slides_1_14_deck_overlay_rules.json"

OPENING_OVERLAY = "output/premium_opening_slides_1_3_v0_1/opening_slides_1_3_overlay_spec.json"
OPENING_QA = "output/premium_opening_slides_1_3_v0_1/opening_slides_1_3_qa.json"
METHODS_OVERLAY = "output/premium_methods_dark_variants_slides_4_5_v0_1/slides_4_5_dark_methods_overlay_spec.json"
METHODS_QA = "output/premium_methods_dark_variants_slides_4_5_v0_1/slides_4_5_dark_methods_qa.json"
FEATURE_OVERLAY = "output/premium_feature_layer_bridge_v0_4/feature_layer_bridge_v0_4_overlay_spec.json"
FEATURE_QA = "output/premium_feature_layer_bridge_v0_4/feature_layer_bridge_v0_4_qa.json"
CORE_CAPTION_PACK = "output/premium_core_result_slides/core_result_set_qa_v0_1/core_result_slides_10_14_caption_pack.json"
CORE_INTEGRATED_QA = "output/premium_core_result_slides/core_result_set_qa_v0_1/core_result_slides_10_14_integrated_qa.json"

EARLY_RESULT_OVERLAYS = {
    7: "output/premium_slide_scenes/fig1_tissue_transfer_editable_overlay_spec.json",
    8: "output/premium_slide_scenes/fig2_pathway_editable_overlay_spec.json",
    9: "output/premium_slide_scenes/fig3_model_tier_editable_overlay_spec.json",
}

EARLY_RESULT_SCENE_PLATES = {
    7: "output/premium_slide_scenes/fig1_tissue_transfer_scene_plate.png",
    8: "output/premium_slide_scenes/fig2_pathway_scene_plate.png",
    9: "output/premium_slide_scenes/fig3_model_tier_scene_plate.png",
}

CORE_MANIFESTS = {
    10: "output/premium_core_result_slides/slide10_v4_hardening_v0_1/slide10_v4_hardening_manifest.json",
    11: "output/premium_core_result_slides/slide11_temporal_rrrm_guardrails_v0_1/slide11_temporal_rrrm_guardrails_manifest.json",
    12: "output/premium_core_result_slides/slide12_negative_boundary_v0_1/slide12_negative_boundary_manifest.json",
    13: "output/premium_core_result_slides/slide13_biological_interpretation_v0_1/slide13_biological_interpretation_manifest.json",
    14: "output/premium_core_result_slides/slide14_human_translation_ladder_v0_1/slide14_human_translation_ladder_manifest.json",
}

SPEAKER_NOTES_REQUIRED = {
    4: [
        "Define the task unit: one source study becomes auditable benchmark tasks with source context retained.",
        "Explain that B2 is intentionally kept in notes/backup to avoid a dense implementation diagram.",
    ],
    5: [
        "Define mission-held-out: the hidden test set is an entire mission, not random sample validation.",
        "State that train-only choices must be fixed before the held-out mission is scored.",
    ],
    6: [
        "Bridge gene-level views to pathway-summary views before showing tissue and model comparisons.",
    ],
    10: [
        "Clarify that hardening means coverage and robustness, not a new leaderboard.",
    ],
    11: [
        "Clarify preservation, descriptive recovery projection, and underpowered RRRM composition testing.",
    ],
    12: [
        "Clarify that negative results define current task limits, not absence of biology.",
    ],
    13: [
        "Clarify that hits become follow-up hypotheses, not mechanisms, treatments, or countermeasures.",
    ],
    14: [
        "Clarify partial pathway/target-tier alignment and weak direct gene transfer.",
    ],
}

PPTX_INSTRUCTIONS = {
    1: "Use the scene plate as atmospheric Z0-Z2; rebuild title, subtitle, bridge, caveat, and source as editable text.",
    2: "Keep the external map as positioning only; avoid ranking or firstness language.",
    3: "Preserve evidence-level separation between completed results and v8/v9 extensions.",
    4: "Use dark methods scene as main slide; put the task-unit definition in notes or backup.",
    5: "Use dark methods scene as main slide; put train-only guardrail mechanics in notes or backup.",
    6: "Use the feature bridge to define what the model sees before any score surfaces.",
    7: "Keep tissue-specific transfer as the first result; do not generalize across all tissues.",
    8: "Frame pathway summaries as selected-task rescue, not universal superiority.",
    9: "Frame model comparison as tested-setting evidence; do not universalize foundation-model failure.",
    10: "Keep the hardening claim as coverage/robustness evidence, not leaderboard rank.",
    11: "Keep temporal/RRRM statements as guardrails; do not over-read recovery or composition.",
    12: "Keep the negative boundary as task-limit evidence, not universal biological absence.",
    13: "Keep biological interpretation as hypotheses and target triage only.",
    14: "Keep human translation as partial alignment only; no clinical or countermeasure implication.",
}

JARGON_WATCH_TERMS = [
    "AUROC",
    "RRRM",
    "PCA-LR",
    "KEGG",
    "ISS-T",
    "LAR",
    "GeneLab",
    "OSDR",
    "v1-v9",
    "v8/v9",
    "foundation-model",
    "task unit",
    "held-out",
]

FORBIDDEN_VISIBLE_PATTERNS = {
    "placeholder": re.compile(r"placeholder", re.IGNORECASE),
    "production brief": re.compile(r"production brief", re.IGNORECASE),
    "wireframe": re.compile(r"wireframe", re.IGNORECASE),
    "micro-plan": re.compile(r"micro-plan", re.IGNORECASE),
    "local absolute path": re.compile(r"/Users/", re.IGNORECASE),
    "sklearn": re.compile(r"sklearn", re.IGNORECASE),
    "function": re.compile(r"\bfunction\b", re.IGNORECASE),
    "class": re.compile(r"\bclass\b", re.IGNORECASE),
    "payload": re.compile(r"\bpayload\b", re.IGNORECASE),
    "LOMO": re.compile(r"\bLOMO\b"),
    "random CV": re.compile(r"random CV", re.IGNORECASE),
    "cross-validation": re.compile(r"cross-validation", re.IGNORECASE),
    "internal decision": re.compile(r"internal decision", re.IGNORECASE),
    "draft": re.compile(r"\bdraft\b", re.IGNORECASE),
    "generated schematic texture": re.compile(r"generated schematic texture", re.IGNORECASE),
    "not source evidence": re.compile(r"not source evidence", re.IGNORECASE),
    "replacement test": re.compile(r"replacement test", re.IGNORECASE),
    "same split": re.compile(r"same split", re.IGNORECASE),
    "M1": re.compile(r"\bM1\b"),
    "M2": re.compile(r"\bM2\b"),
    "M3": re.compile(r"\bM3\b"),
    "M4": re.compile(r"\bM4\b"),
}

PRE_SLIDE_20_EXTENSION_TERMS = {
    "organoid": re.compile(r"organoid", re.IGNORECASE),
    "OSD-120": re.compile(r"OSD-120", re.IGNORECASE),
}


def rel(path: Path) -> str:
    return str(path.relative_to(ROOT))


def p(path: str) -> Path:
    return ROOT / path


def load_json(path: str) -> dict[str, Any]:
    with p(path).open("r", encoding="utf-8") as handle:
        return json.load(handle)


def maybe_load_json(path: str) -> dict[str, Any]:
    candidate = p(path)
    if not candidate.exists():
        return {}
    return load_json(path)


def dedupe(items: list[str]) -> list[str]:
    seen = set()
    out: list[str] = []
    for item in items:
        if not item:
            continue
        normalized = " ".join(str(item).split())
        if normalized and normalized not in seen:
            out.append(normalized)
            seen.add(normalized)
    return out


def image_size(path: str | None) -> list[int] | None:
    if not path or Image is None:
        return None
    candidate = p(path)
    if not candidate.exists():
        return None
    with Image.open(candidate) as image:
        return [image.width, image.height]


def ensure_audit_grayscale(slide_number: int, title: str, asset: str | None) -> str | None:
    if Image is None or not asset:
        return None
    candidate = p(asset)
    if not candidate.exists():
        return None
    AUDIT_QA.mkdir(parents=True, exist_ok=True)
    safe_title = re.sub(r"[^a-z0-9]+", "_", title.lower()).strip("_")
    output = AUDIT_QA / f"slide{slide_number:02d}_{safe_title}_grayscale.png"
    if not output.exists():
        with Image.open(candidate).convert("RGB") as image:
            image.convert("L").convert("RGB").save(output, quality=95)
    return rel(output)


def find_slide(spec: dict[str, Any], slide_number: int) -> dict[str, Any]:
    for slide in spec.get("slides", []):
        if int(slide.get("slide", -1)) == slide_number:
            return slide
    return {}


def visible_text_from_overlay(slide_number: int, manifest_slide: dict[str, Any], loaded: dict[str, dict[str, Any]]) -> tuple[list[str], str, str | None]:
    if slide_number in {1, 2, 3}:
        item = find_slide(loaded["opening_overlay"], slide_number)
        editable = item.get("editable_text", {})
        return dedupe(list(editable.values())), OPENING_OVERLAY, item.get("scene_plate")

    if slide_number in {4, 5}:
        item = find_slide(loaded["methods_overlay"], slide_number)
        editable = item.get("editable_text", {})
        return dedupe(list(editable.values())), METHODS_OVERLAY, item.get("scene_plate")

    if slide_number == 6:
        spec = loaded["feature_overlay"]
        texts = [item.get("content", "") for item in spec.get("text", [])]
        scene_plate = "output/premium_feature_layer_bridge_v0_4/feature_layer_bridge_v0_4_scene_plate.png"
        return dedupe(texts), FEATURE_OVERLAY, scene_plate

    if slide_number in EARLY_RESULT_OVERLAYS:
        spec = loaded[f"early_overlay_{slide_number}"]
        texts = [item.get("text", "") for item in spec.get("editable_overlay_elements", [])]
        texts.append(manifest_slide.get("claim_boundary", ""))
        return dedupe(texts), EARLY_RESULT_OVERLAYS[slide_number], EARLY_RESULT_SCENE_PLATES[slide_number]

    if slide_number in CORE_MANIFESTS:
        core_slide = find_slide(loaded["core_caption_pack"], slide_number)
        texts = [
            core_slide.get("headline", ""),
            core_slide.get("plain_caption", ""),
            core_slide.get("claim_boundary", ""),
        ]
        core_manifest = loaded[f"core_manifest_{slide_number}"]
        scene_plate = core_manifest.get("outputs", {}).get("scene_plate")
        return dedupe(texts), CORE_CAPTION_PACK, scene_plate

    fallback = [
        manifest_slide.get("headline", ""),
        manifest_slide.get("caption", ""),
        manifest_slide.get("claim_boundary", ""),
    ]
    return dedupe(fallback), ASSEMBLY_MANIFEST, None


def term_hits(texts: list[str], patterns: dict[str, re.Pattern[str]]) -> list[str]:
    joined = "\n".join(texts)
    return sorted([name for name, pattern in patterns.items() if pattern.search(joined)])


def jargon_hits(texts: list[str]) -> list[str]:
    joined_lower = "\n".join(texts).lower()
    return sorted([term for term in JARGON_WATCH_TERMS if term.lower() in joined_lower])


def word_count(texts: list[str]) -> int:
    joined = " ".join(texts)
    return len(re.findall(r"[A-Za-z0-9][A-Za-z0-9_/.-]*", joined))


def source_or_scope_present(slide_number: int, texts: list[str], loaded: dict[str, dict[str, Any]]) -> bool:
    joined_lower = "\n".join(texts).lower()
    source_scope_tokens = [
        "source",
        "gene",
        "mouse",
        "diagnostic check",
        "shared",
        "public",
        "study records",
        "master outline",
        "review",
        "narration pack",
        "split rule",
        "mission-held-out",
        "v1-v7",
    ]
    if any(token in joined_lower for token in source_scope_tokens):
        return True
    if slide_number in CORE_MANIFESTS:
        return bool(loaded[f"core_manifest_{slide_number}"].get("source_documents"))
    return False


def grayscale_path_for(slide_number: int, loaded: dict[str, dict[str, Any]]) -> str | None:
    if slide_number in {1, 2, 3}:
        item = find_slide(loaded["opening_overlay"], slide_number)
        preview = item.get("rendered_preview")
        if preview:
            stem = Path(preview).stem
            return f"output/premium_opening_slides_1_3_v0_1/qa/{stem}_grayscale.png"
    if slide_number in {4, 5}:
        item = find_slide(loaded["methods_overlay"], slide_number)
        preview = item.get("rendered_preview")
        if preview:
            stem = Path(preview).stem
            return f"output/premium_methods_dark_variants_slides_4_5_v0_1/qa/{stem}_grayscale.png"
    if slide_number == 6:
        return "output/premium_feature_layer_bridge_v0_4/qa/feature_layer_bridge_v0_4_rendered_preview_grayscale.png"
    if slide_number in CORE_MANIFESTS:
        return loaded[f"core_manifest_{slide_number}"].get("outputs", {}).get("grayscale_qa")
    return None


def build_audit() -> dict[str, Any]:
    loaded: dict[str, dict[str, Any]] = {
        "assembly_manifest": load_json(ASSEMBLY_MANIFEST),
        "assembly_qa": load_json(ASSEMBLY_QA),
        "assembly_overlay_rules": load_json(ASSEMBLY_OVERLAY_RULES),
        "opening_overlay": load_json(OPENING_OVERLAY),
        "opening_qa": load_json(OPENING_QA),
        "methods_overlay": load_json(METHODS_OVERLAY),
        "methods_qa": load_json(METHODS_QA),
        "feature_overlay": load_json(FEATURE_OVERLAY),
        "feature_qa": load_json(FEATURE_QA),
        "core_caption_pack": load_json(CORE_CAPTION_PACK),
        "core_integrated_qa": load_json(CORE_INTEGRATED_QA),
    }
    for slide_number, overlay in EARLY_RESULT_OVERLAYS.items():
        loaded[f"early_overlay_{slide_number}"] = load_json(overlay)
    for slide_number, manifest in CORE_MANIFESTS.items():
        loaded[f"core_manifest_{slide_number}"] = load_json(manifest)

    slide_rows: list[dict[str, Any]] = []
    for manifest_slide in loaded["assembly_manifest"]["slides"]:
        slide_number = int(manifest_slide["slide"])
        asset = manifest_slide.get("asset")
        backup_assets = manifest_slide.get("backup_assets", [])
        visible_texts, overlay_source, scene_plate = visible_text_from_overlay(slide_number, manifest_slide, loaded)
        forbidden_hits = term_hits(visible_texts, FORBIDDEN_VISIBLE_PATTERNS)
        extension_hits = term_hits(visible_texts, PRE_SLIDE_20_EXTENSION_TERMS)
        watch_hits = jargon_hits(visible_texts)
        grayscale = grayscale_path_for(slide_number, loaded)

        asset_exists = bool(asset and p(asset).exists())
        scene_plate_exists = bool(scene_plate and p(scene_plate).exists())
        overlay_source_exists = bool(overlay_source and p(overlay_source).exists())
        claim_boundary_present = bool(manifest_slide.get("claim_boundary")) or any(
            item for item in visible_texts if any(token in item.lower() for token in ["caveat", "not ", "only", "no "])
        )
        source_scope = source_or_scope_present(slide_number, visible_texts, loaded)
        if not grayscale or not p(grayscale).exists():
            grayscale = ensure_audit_grayscale(slide_number, manifest_slide.get("title", f"slide_{slide_number}"), asset)
        grayscale_exists = bool(grayscale and p(grayscale).exists())
        backup_existing = [item for item in backup_assets if p(item).exists()]
        backup_missing = [item for item in backup_assets if not p(item).exists()]

        blockers: list[str] = []
        warnings: list[str] = []
        if not asset_exists:
            blockers.append("missing rendered preview asset")
        if not scene_plate_exists:
            blockers.append("missing scene plate")
        if not overlay_source_exists:
            blockers.append("missing editable overlay/caption source")
        if not claim_boundary_present:
            blockers.append("missing claim boundary")
        if forbidden_hits:
            blockers.append("forbidden visible terms present")
        if extension_hits:
            blockers.append("pre-slide-20 extension term visible")
        if not grayscale_exists:
            warnings.append("missing grayscale QA output")
        if not source_scope:
            warnings.append("source or scope line needs explicit PPTX text")
        if backup_missing:
            warnings.append("backup assets missing")
        if slide_number in SPEAKER_NOTES_REQUIRED:
            warnings.append("speaker notes required")
        if watch_hits:
            warnings.append("viewer-language bridge needed for jargon")

        slide_rows.append(
            {
                "slide": slide_number,
                "section": manifest_slide.get("section"),
                "title": manifest_slide.get("title"),
                "role": manifest_slide.get("role"),
                "headline": manifest_slide.get("headline"),
                "asset": asset,
                "asset_exists": asset_exists,
                "asset_size_px": image_size(asset),
                "scene_plate": scene_plate,
                "scene_plate_exists": scene_plate_exists,
                "grayscale_qa": grayscale,
                "grayscale_qa_exists": grayscale_exists,
                "editable_overlay_source": overlay_source,
                "editable_overlay_ready": overlay_source_exists,
                "claim_boundary": manifest_slide.get("claim_boundary"),
                "claim_boundary_present": claim_boundary_present,
                "source_or_scope_present": source_scope,
                "backup_assets_existing": backup_existing,
                "backup_assets_missing": backup_missing,
                "visible_text_word_count": word_count(visible_texts),
                "visible_text_sample": visible_texts,
                "forbidden_visible_hits": forbidden_hits,
                "pre_slide_20_extension_hits": extension_hits,
                "viewer_language_watch_hits": watch_hits,
                "speaker_notes_required": SPEAKER_NOTES_REQUIRED.get(slide_number, []),
                "pptx_instruction": PPTX_INSTRUCTIONS.get(slide_number, ""),
                "blockers": blockers,
                "warnings": warnings,
                "pptx_status": "blocked" if blockers else "ready_with_required_overlay_rebuild",
            }
        )

    blockers = [
        {"slide": row["slide"], "title": row["title"], "blockers": row["blockers"]}
        for row in slide_rows
        if row["blockers"]
    ]
    warnings = [
        {"slide": row["slide"], "title": row["title"], "warnings": row["warnings"]}
        for row in slide_rows
        if row["warnings"]
    ]
    return {
        "created": CREATED,
        "scope": "slides 1-14 PPTX production-readiness audit",
        "decision": "ready for PPTX assembly as a production draft if every slide is rebuilt as layered PNG scene plus editable scientific interpretation",
        "inputs": {
            "assembly_manifest": ASSEMBLY_MANIFEST,
            "assembly_qa": ASSEMBLY_QA,
            "assembly_overlay_rules": ASSEMBLY_OVERLAY_RULES,
            "opening_overlay": OPENING_OVERLAY,
            "methods_overlay": METHODS_OVERLAY,
            "feature_overlay": FEATURE_OVERLAY,
            "early_result_overlays": EARLY_RESULT_OVERLAYS,
            "core_caption_pack": CORE_CAPTION_PACK,
            "core_manifests": CORE_MANIFESTS,
        },
        "automatic_summary": {
            "slides_audited": len(slide_rows),
            "slides_ready": sum(1 for row in slide_rows if not row["blockers"]),
            "slides_blocked": len(blockers),
            "slides_with_speaker_note_requirements": sum(1 for row in slide_rows if row["speaker_notes_required"]),
            "slides_with_viewer_language_watch_hits": sum(1 for row in slide_rows if row["viewer_language_watch_hits"]),
            "all_assets_exist": all(row["asset_exists"] for row in slide_rows),
            "all_scene_plates_exist": all(row["scene_plate_exists"] for row in slide_rows),
            "all_overlay_sources_exist": all(row["editable_overlay_ready"] for row in slide_rows),
            "all_grayscale_qa_exist": all(row["grayscale_qa_exists"] for row in slide_rows),
            "forbidden_visible_terms_pass": all(not row["forbidden_visible_hits"] for row in slide_rows),
            "pre_slide_20_extension_terms_pass": all(not row["pre_slide_20_extension_hits"] for row in slide_rows),
        },
        "required_actions_before_final_export": [
            "Do not paste rendered previews as final all-in-one slides; use scene plates as Z0-Z2 proof surfaces and rebuild text as editable PPTX objects.",
            "Insert speaker notes for slides 4, 5, 6, and 10-14 before presenter rehearsal.",
            "Keep B2 and B4 mechanics in notes or backup unless the methods section is expanded.",
            "Keep organoid and OSD-120 extension visuals out of slides 1-14; introduce them only after the completed v1-v7 result spine.",
            "Move table-like evidence to backup/manuscript tables, not main-figure cards.",
        ],
        "blockers": blockers,
        "warnings": warnings,
        "slides": slide_rows,
    }


def write_json(audit: dict[str, Any]) -> Path:
    path = OUT / "slides_1_14_pptx_readiness_audit.json"
    path.write_text(json.dumps(audit, indent=2), encoding="utf-8")
    return path


def write_csv(audit: dict[str, Any]) -> Path:
    path = OUT / "slides_1_14_pptx_readiness_matrix.csv"
    fields = [
        "slide",
        "section",
        "title",
        "role",
        "asset_exists",
        "scene_plate_exists",
        "grayscale_qa_exists",
        "editable_overlay_ready",
        "claim_boundary_present",
        "source_or_scope_present",
        "visible_text_word_count",
        "forbidden_visible_hits",
        "pre_slide_20_extension_hits",
        "viewer_language_watch_hits",
        "speaker_notes_required",
        "pptx_status",
    ]
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        for row in audit["slides"]:
            writer.writerow(
                {
                    key: "; ".join(row[key]) if isinstance(row.get(key), list) else row.get(key)
                    for key in fields
                }
            )
    return path


def markdown_bool(value: bool) -> str:
    return "yes" if value else "no"


def write_markdown(audit: dict[str, Any]) -> Path:
    path = OUT / "slides_1_14_pptx_readiness_summary.md"
    summary = audit["automatic_summary"]
    lines = [
        "# Slides 1-14 PPTX Readiness Summary",
        "",
        f"Date: {CREATED}",
        "",
        "## Decision",
        "",
        audit["decision"],
        "",
        "## Automatic Gate",
        "",
        f"- Slides audited: {summary['slides_audited']}",
        f"- Slides ready: {summary['slides_ready']}",
        f"- Slides blocked: {summary['slides_blocked']}",
        f"- All rendered assets exist: {markdown_bool(summary['all_assets_exist'])}",
        f"- All scene plates exist: {markdown_bool(summary['all_scene_plates_exist'])}",
        f"- All editable overlay/caption sources exist: {markdown_bool(summary['all_overlay_sources_exist'])}",
        f"- All grayscale QA outputs exist: {markdown_bool(summary['all_grayscale_qa_exist'])}",
        f"- Forbidden visible terms pass: {markdown_bool(summary['forbidden_visible_terms_pass'])}",
        f"- Pre-slide-20 extension terms pass: {markdown_bool(summary['pre_slide_20_extension_terms_pass'])}",
        "",
        "## Required Actions Before Final Export",
        "",
    ]
    lines.extend([f"- {item}" for item in audit["required_actions_before_final_export"]])
    lines.extend(
        [
            "",
            "## Slide Matrix",
            "",
            "| Slide | Title | Asset | Scene | Overlay | Boundary | Source/scope | Watch | Status |",
            "|---:|---|---|---|---|---|---|---|---|",
        ]
    )
    for row in audit["slides"]:
        watch = []
        if row["speaker_notes_required"]:
            watch.append("notes")
        if row["viewer_language_watch_hits"]:
            watch.append("jargon")
        if row["warnings"] and not watch:
            watch.append("warn")
        watch_text = ", ".join(watch) if watch else "clear"
        lines.append(
            "| {slide} | {title} | {asset} | {scene} | {overlay} | {boundary} | {source} | {watch} | {status} |".format(
                slide=row["slide"],
                title=row["title"],
                asset=markdown_bool(row["asset_exists"]),
                scene=markdown_bool(row["scene_plate_exists"]),
                overlay=markdown_bool(row["editable_overlay_ready"]),
                boundary=markdown_bool(row["claim_boundary_present"]),
                source=markdown_bool(row["source_or_scope_present"]),
                watch=watch_text,
                status=row["pptx_status"],
            )
        )
    lines.extend(["", "## Viewer-Language Watchpoints", ""])
    for row in audit["slides"]:
        if row["viewer_language_watch_hits"]:
            lines.append(
                f"- Slide {row['slide']} ({row['title']}): {', '.join(row['viewer_language_watch_hits'])}"
            )
    lines.extend(["", "## Speaker Notes Required", ""])
    for row in audit["slides"]:
        if row["speaker_notes_required"]:
            lines.append(f"- Slide {row['slide']} ({row['title']}): {'; '.join(row['speaker_notes_required'])}")
    if not audit["blockers"]:
        lines.extend(["", "## Blockers", "", "None."])
    path.write_text("\n".join(lines), encoding="utf-8")
    return path


COLORS = {
    "bg": "#070A0E",
    "panel": "#101923",
    "ink": "#F4F7F8",
    "soft": "#B9C7D2",
    "muted": "#788898",
    "line": "#33465A",
    "teal": "#1AA090",
    "amber": "#B4832F",
    "rose": "#B45A7E",
}


def rgb(color: str) -> tuple[int, int, int]:
    value = COLORS[color].lstrip("#")
    return tuple(int(value[i : i + 2], 16) for i in (0, 2, 4))


def rgba(color: str, alpha: int) -> tuple[int, int, int, int]:
    return rgb(color) + (alpha,)


def font(size: int, *, bold: bool = False):
    if ImageFont is None:
        return None
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


def draw_status_board(audit: dict[str, Any]) -> Path | None:
    if Image is None or ImageDraw is None or ImageFilter is None:
        return None
    board = Image.new("RGBA", (3840, 2160), rgb("bg") + (255,))
    draw = ImageDraw.Draw(board)
    f_eyebrow = font(22, bold=True)
    f_title = font(62, bold=True)
    f_body = font(24)
    f_small = font(17)
    f_tiny = font(13)

    draw.text((132, 92), "PPTX READINESS AUDIT", font=f_eyebrow, fill=rgb("teal"))
    draw.text((132, 142), "Slides 1-14 are ready for layered assembly", font=f_title, fill=rgb("ink"))
    draw.text(
        (136, 224),
        "Gate: use scene plates as proof surfaces, then rebuild interpretation, caveat, source, and notes as editable slide objects.",
        font=f_body,
        fill=rgb("soft"),
    )

    rows = audit["slides"]
    card_w, card_h = 430, 242
    x0, y0 = 130, 350
    xgap, ygap = 54, 230
    for index, row in enumerate(rows):
        col = index % 7
        grid_row = index // 7
        x = x0 + col * (card_w + xgap)
        y = y0 + grid_row * (card_h + ygap)
        shadow = Image.new("RGBA", board.size, (0, 0, 0, 0))
        sd = ImageDraw.Draw(shadow)
        sd.rounded_rectangle((x + 10, y + 14, x + card_w + 10, y + card_h + 14), radius=12, fill=(0, 0, 0, 120))
        board.alpha_composite(shadow.filter(ImageFilter.GaussianBlur(10)))
        if row["asset_exists"]:
            thumb = Image.open(p(row["asset"])).convert("RGB")
            resample = getattr(getattr(Image, "Resampling", Image), "LANCZOS")
            thumb.thumbnail((card_w, card_h), resample)
            plate = Image.new("RGB", (card_w, card_h), rgb("panel"))
            plate.paste(thumb, ((card_w - thumb.width) // 2, (card_h - thumb.height) // 2))
            board.alpha_composite(plate.convert("RGBA"), (x, y))
        else:
            draw.rounded_rectangle((x, y, x + card_w, y + card_h), radius=10, fill=rgb("panel"))
        outline = "teal" if not row["blockers"] else "rose"
        draw.rounded_rectangle((x - 1, y - 1, x + card_w + 1, y + card_h + 1), radius=11, outline=rgba(outline, 190), width=3)
        status = "READY" if not row["blockers"] else "BLOCKED"
        draw.rounded_rectangle((x + 18, y + 18, x + 128, y + 48), radius=6, fill=rgba(outline, 215))
        draw.text((x + 30, y + 25), status, font=f_tiny, fill=rgb("ink"))
        draw.text((x, y + card_h + 22), f"{row['slide']}. {row['title']}", font=f_small, fill=rgb("ink"))
        watch = []
        if row["speaker_notes_required"]:
            watch.append("notes")
        if row["viewer_language_watch_hits"]:
            watch.append("jargon")
        watch_text = " + ".join(watch) if watch else "clear"
        draw.text((x, y + card_h + 48), f"watch: {watch_text}", font=f_tiny, fill=rgb("muted"))

    y = 1510
    panels = [
        ("assets", audit["automatic_summary"]["all_assets_exist"]),
        ("scene plates", audit["automatic_summary"]["all_scene_plates_exist"]),
        ("editable overlays", audit["automatic_summary"]["all_overlay_sources_exist"]),
        ("grayscale QA", audit["automatic_summary"]["all_grayscale_qa_exist"]),
        ("visible terms", audit["automatic_summary"]["forbidden_visible_terms_pass"]),
        ("no early extensions", audit["automatic_summary"]["pre_slide_20_extension_terms_pass"]),
    ]
    for idx, (label, passed) in enumerate(panels):
        x = 160 + idx * 590
        color = "teal" if passed else "rose"
        draw.rounded_rectangle((x, y, x + 510, y + 106), radius=10, fill=rgba("panel", 225), outline=rgba(color, 160), width=2)
        draw.text((x + 28, y + 28), label, font=f_body, fill=rgb("ink"))
        draw.text((x + 28, y + 66), "pass" if passed else "fail", font=f_small, fill=rgb(color))

    draw.rounded_rectangle((130, 1875, 3680, 2044), radius=12, fill=rgba("panel", 232), outline=rgba("line", 150), width=2)
    draw.text((166, 1914), "Residual production work", font=f_body, fill=rgb("ink"))
    draw.text(
        (166, 1964),
        "Rebuild all titles, caveats, sources, and bridge lines as editable PPTX objects; add notes for methods and boundary-heavy result slides.",
        font=f_small,
        fill=rgb("soft"),
    )
    output = OUT / "slides_1_14_pptx_readiness_status_board.png"
    board.convert("RGB").save(output, quality=95)
    return output


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    audit = build_audit()
    paths = {
        "json": rel(write_json(audit)),
        "csv": rel(write_csv(audit)),
        "markdown": rel(write_markdown(audit)),
    }
    board = draw_status_board(audit)
    if board:
        paths["status_board"] = rel(board)
    manifest = {
        "created": CREATED,
        "artifact_role": "slides 1-14 PPTX production-readiness audit",
        "outputs": paths,
        "summary": audit["automatic_summary"],
    }
    manifest_path = OUT / "slides_1_14_pptx_readiness_manifest.json"
    manifest_path.write_text(json.dumps(manifest, indent=2), encoding="utf-8")
    print(json.dumps(manifest, indent=2))


if __name__ == "__main__":
    main()
