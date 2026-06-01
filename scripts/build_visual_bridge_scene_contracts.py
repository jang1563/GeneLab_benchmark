#!/usr/bin/env python3
"""Build pre-render visual production contracts for bridge slides."""

from __future__ import annotations

import json
from pathlib import Path
from typing import Any


ROOT = Path(__file__).resolve().parents[1]
OUT_ROOT = ROOT / "output" / "premium_bridge_scenes"
CANVAS = {"width_px": 3840, "height_px": 2160, "aspect_ratio": "16:9"}
CREATED = "2026-06-01"


def rel(path: Path) -> str:
    return str(path.relative_to(ROOT))


def word_count(items: list[dict[str, Any]]) -> int:
    words = 0
    for item in items:
        content = str(item.get("content", ""))
        words += len([token for token in content.replace("/", " ").split() if token])
    return words


SLIDES: list[dict[str, Any]] = [
    {
        "slide_id": "b3_mission_held_out",
        "stage": "pre_render",
        "decision_headline": "The test set is a mission, not a random sample",
        "audience_question": "Why is mission-held-out evaluation central to this benchmark?",
        "claim_boundary": "held-out mission evaluation; not random sample cross-validation",
        "content_brief": "docs/VISUAL_BRIDGE_CONTENT_BRIEFS_B1_B4_2026_06_01.md",
        "technical_preflight": "docs/VISUAL_BRIDGE_TECHNICAL_PREFLIGHT_B3_B4_2026_06_01.md",
        "visual_move": "one mission moves behind a clean boundary while remaining missions feed training",
        "evidence_sources": [
            {
                "path": "docs/VISUAL_BRIDGE_CONTENT_BRIEFS_B1_B4_2026_06_01.md",
                "role": "B3 content brief and production gate",
            },
            {
                "path": "docs/VISUAL_METHODS_EXPLANATION_GAP_MAP_2026_06_01.md",
                "role": "viewer confusion map for held-out mission evaluation",
            },
            {
                "path": "docs/VISUAL_METHODS_STORYBOARD_2026_06_01.md",
                "role": "methods section sequence and bridge role",
            },
            {
                "path": "docs/V1_V9_PRESENTATION_AND_MANUSCRIPT_MASTER_OUTLINE_2026_05_31.md",
                "role": "overall manuscript/deck narrative spine",
            },
            {
                "path": "docs/PROJECT_SLIDE_CONTENT_INVENTORY_V1_TO_V9_2026_05_31.md",
                "role": "project-level mission-held-out benchmark framing",
            },
        ],
        "forbidden_visible_terms": [
            "LOMO",
            "random CV",
            "cross-validation",
            "payload",
            "RRRM",
            "alpha",
            "NES",
            "macro-F1",
            "/Users/",
        ],
        "overlay": {
            "text": [
                {
                    "id": "headline",
                    "role": "decision_headline",
                    "content": "The test set is a mission, not a random sample",
                    "x": 0.065,
                    "y": 0.135,
                    "font_pt": 27,
                    "color": "ink",
                    "max_lines": 2,
                    "z": "Z3",
                },
                {
                    "id": "train_label",
                    "role": "primary_callout",
                    "content": "train on prior missions",
                    "x": 0.205,
                    "y": 0.545,
                    "font_pt": 11,
                    "color": "muted",
                    "max_lines": 1,
                    "z": "Z3",
                },
                {
                    "id": "hidden_label",
                    "role": "primary_callout",
                    "content": "hide one mission",
                    "x": 0.705,
                    "y": 0.545,
                    "font_pt": 12,
                    "color": "red",
                    "max_lines": 1,
                    "z": "Z3",
                },
                {
                    "id": "score_label",
                    "role": "secondary_label",
                    "content": "score after training",
                    "x": 0.755,
                    "y": 0.710,
                    "font_pt": 8.5,
                    "color": "muted",
                    "max_lines": 1,
                    "z": "Z3",
                },
            ],
            "status_labels": [
                {
                    "id": "scope",
                    "role": "claim_boundary",
                    "content": "held-out mission evaluation",
                    "x": 0.065,
                    "y": 0.865,
                    "font_pt": 7.8,
                    "color": "muted",
                    "z": "Z4",
                },
                {
                    "id": "caveat",
                    "role": "trust_caveat",
                    "content": "No hidden-mission samples are used for training.",
                    "x": 0.065,
                    "y": 0.900,
                    "font_pt": 7.4,
                    "color": "muted",
                    "z": "Z4",
                },
            ],
            "focus_marks": [
                {
                    "id": "test_boundary",
                    "role": "boundary",
                    "shape": "vertical_rule",
                    "x": 0.640,
                    "y0": 0.320,
                    "y1": 0.725,
                    "color": "red",
                    "z": "Z5",
                }
            ],
        },
    },
    {
        "slide_id": "b4_train_only_guard",
        "stage": "pre_render",
        "decision_headline": "Feature processing stays on the training side",
        "audience_question": "How does the benchmark avoid learning from the mission it is supposed to test?",
        "claim_boundary": "feature selection, scaling, and model fitting are learned from training missions only",
        "content_brief": "docs/VISUAL_BRIDGE_CONTENT_BRIEFS_B1_B4_2026_06_01.md",
        "technical_preflight": "docs/VISUAL_BRIDGE_TECHNICAL_PREFLIGHT_B3_B4_2026_06_01.md",
        "visual_move": "two lanes move in parallel until the hidden mission joins only at scoring",
        "evidence_sources": [
            {
                "path": "docs/VISUAL_BRIDGE_CONTENT_BRIEFS_B1_B4_2026_06_01.md",
                "role": "B4 content brief and production gate",
            },
            {
                "path": "docs/VISUAL_METHODS_EXPLANATION_GAP_MAP_2026_06_01.md",
                "role": "viewer confusion map for leakage prevention",
            },
            {
                "path": "docs/VISUAL_METHODS_STORYBOARD_2026_06_01.md",
                "role": "methods sequence and bridge role",
            },
            {
                "path": "docs/VISUAL_METHODS_EXPLAINER_PILOT_2026_06_01.md",
                "role": "first methods explainer review and text simplification record",
            },
            {
                "path": "output/premium_methods_scenes/methods_data_to_evaluation_overview_manifest.json",
                "role": "current methods overview artifact manifest",
            },
        ],
        "forbidden_visible_terms": [
            "LOMO",
            "payload",
            "artifact",
            "macro-F1",
            "/Users/",
            "sklearn",
            "StandardScaler",
            "fit_transform",
        ],
        "overlay": {
            "text": [
                {
                    "id": "headline",
                    "role": "decision_headline",
                    "content": "Feature processing stays on the training side",
                    "x": 0.065,
                    "y": 0.135,
                    "font_pt": 27,
                    "color": "ink",
                    "max_lines": 2,
                    "z": "Z3",
                },
                {
                    "id": "choose_features",
                    "role": "process_verb",
                    "content": "choose features",
                    "x": 0.355,
                    "y": 0.470,
                    "font_pt": 9.0,
                    "color": "blue",
                    "max_lines": 1,
                    "z": "Z3",
                },
                {
                    "id": "scale",
                    "role": "process_verb",
                    "content": "scale",
                    "x": 0.505,
                    "y": 0.470,
                    "font_pt": 9.0,
                    "color": "teal",
                    "max_lines": 1,
                    "z": "Z3",
                },
                {
                    "id": "fit_model",
                    "role": "process_verb",
                    "content": "fit model",
                    "x": 0.650,
                    "y": 0.470,
                    "font_pt": 9.0,
                    "color": "purple",
                    "max_lines": 1,
                    "z": "Z3",
                },
                {
                    "id": "score_hidden",
                    "role": "process_verb",
                    "content": "score hidden mission",
                    "x": 0.805,
                    "y": 0.610,
                    "font_pt": 9.0,
                    "color": "red",
                    "max_lines": 1,
                    "z": "Z3",
                },
            ],
            "status_labels": [
                {
                    "id": "scope",
                    "role": "claim_boundary",
                    "content": "train-only processing",
                    "x": 0.065,
                    "y": 0.865,
                    "font_pt": 7.8,
                    "color": "muted",
                    "z": "Z4",
                },
                {
                    "id": "caveat",
                    "role": "trust_caveat",
                    "content": "Feature choices are learned before scoring.",
                    "x": 0.065,
                    "y": 0.900,
                    "font_pt": 7.4,
                    "color": "muted",
                    "z": "Z4",
                },
            ],
            "focus_marks": [
                {
                    "id": "train_test_guard",
                    "role": "guard",
                    "shape": "vertical_rule",
                    "x": 0.730,
                    "y0": 0.335,
                    "y1": 0.705,
                    "color": "red",
                    "z": "Z5",
                }
            ],
        },
    },
]


def build_slide_contract(slide: dict[str, Any]) -> dict[str, Path]:
    slide_id = slide["slide_id"]
    out_dir = OUT_ROOT / slide_id
    out_dir.mkdir(parents=True, exist_ok=True)

    outputs = {
        "scene_plate": rel(out_dir / "scene_plate.png"),
        "rendered_preview_png": rel(out_dir / "rendered_preview.png"),
        "rendered_preview_pdf": rel(out_dir / "rendered_preview.pdf"),
        "overlay_spec": rel(out_dir / "overlay_spec.json"),
        "manifest": rel(out_dir / "manifest.json"),
        "qa": rel(out_dir / "qa.json"),
    }

    text_items = slide["overlay"]["text"] + slide["overlay"]["status_labels"]
    visible_words = word_count(text_items)

    overlay_spec = {
        "slide_id": slide_id,
        "stage": slide["stage"],
        "canvas": CANVAS,
        "coordinate_system": "normalized_0_1",
        "text": slide["overlay"]["text"],
        "status_labels": slide["overlay"]["status_labels"],
        "focus_marks": slide["overlay"]["focus_marks"],
        "visible_word_count": visible_words,
        "visible_word_budget": 45,
        "forbidden_visible_terms": slide["forbidden_visible_terms"],
    }

    manifest = {
        "slide_id": slide_id,
        "created": CREATED,
        "stage": slide["stage"],
        "content_brief": slide["content_brief"],
        "technical_preflight": slide["technical_preflight"],
        "audience_question": slide["audience_question"],
        "decision_headline": slide["decision_headline"],
        "visual_move": slide["visual_move"],
        "claim_boundary": slide["claim_boundary"],
        "evidence_sources": slide["evidence_sources"],
        "generator": "scripts/build_visual_bridge_scene_contracts.py",
        "planned_renderer": "scripts/build_bridge_method_scenes.py",
        "outputs": outputs,
        "qa": {
            "stage": slide["stage"],
            "pre_render_gate_expected": True,
            "post_render_required_before_use": True,
        },
    }

    qa = {
        "slide_id": slide_id,
        "stage": slide["stage"],
        "created": CREATED,
        "pre_render_gate": {
            "content_brief_declared": bool(slide["content_brief"]),
            "technical_preflight_declared": bool(slide["technical_preflight"]),
            "evidence_sources_declared": bool(slide["evidence_sources"]),
            "claim_boundary_declared": bool(slide["claim_boundary"]),
            "visible_text_word_count": visible_words,
            "visible_text_budget": 45,
            "output_paths_declared": True,
            "overlay_spec_declared": True,
            "manifest_declared": True,
        },
        "manual_review_pending": [
            "full_size_render_inspection",
            "thumbnail_contact_sheet_inspection",
            "text_overlap_check",
            "caveat_visibility_check",
        ],
    }

    paths = {
        "overlay_spec": out_dir / "overlay_spec.json",
        "manifest": out_dir / "manifest.json",
        "qa": out_dir / "qa.json",
    }
    paths["overlay_spec"].write_text(json.dumps(overlay_spec, indent=2) + "\n", encoding="utf-8")
    paths["manifest"].write_text(json.dumps(manifest, indent=2) + "\n", encoding="utf-8")
    paths["qa"].write_text(json.dumps(qa, indent=2) + "\n", encoding="utf-8")
    return paths


def main() -> None:
    built = {}
    for slide in SLIDES:
        paths = build_slide_contract(slide)
        built[slide["slide_id"]] = {key: rel(path) for key, path in paths.items()}
    print(json.dumps({"built": built}, indent=2))


if __name__ == "__main__":
    main()
