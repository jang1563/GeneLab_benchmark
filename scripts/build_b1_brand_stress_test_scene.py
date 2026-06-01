#!/usr/bin/env python3
"""Render the B1 brand stress-test bridge scene for SpaceBio-Bench."""

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
import numpy as np
from matplotlib.patches import Circle, FancyArrowPatch, Rectangle

try:
    from PIL import Image
except ImportError:  # pragma: no cover - local QA fallback only
    Image = None  # type: ignore[assignment]


TOKEN_PATH = ROOT / "config" / "visual_identity" / "spacebio_bench_visual_identity_v0_1.json"
OUT_ROOT = ROOT / "output" / "premium_bridge_scenes"
SLIDE_ID = "b1_evaluation_layer"
SLIDE_DIR = OUT_ROOT / SLIDE_ID
CREATED = "2026-06-01"
SLIDE_DPI = 240
CANVAS = {"width_px": 3840, "height_px": 2160, "aspect_ratio": "16:9"}

FORBIDDEN_VISIBLE_TERMS = [
    "payload",
    "RRRM",
    "alpha",
    "LOMO",
    "macro-F1",
    "artifact",
    "wireframe",
    "micro-plan",
    "/Users/",
    "sklearn",
    "function",
    "class",
]

EVIDENCE_SOURCES = [
    {
        "path": "docs/PROJECT_SLIDE_CONTENT_INVENTORY_V1_TO_V9_2026_05_31.md",
        "role": "project version and slide-content inventory",
    },
    {
        "path": "docs/PROJECT_RESULTS_LOCATION_INVENTORY_2026_05_31.md",
        "role": "result and asset location inventory",
    },
    {
        "path": "docs/SIMILAR_PROJECTS_AND_POSITIONING_SCAN_2026_05_31.md",
        "role": "positioning scan for benchmark distinction",
    },
    {"path": "v9/source_inventory.csv", "role": "public source table"},
    {"path": "v9/task_manifest_index.csv", "role": "task record index"},
    {
        "path": "docs/VISUAL_BRIDGE_CONTENT_BRIEFS_B1_B4_2026_06_01.md",
        "role": "B1 bridge content brief",
    },
    {
        "path": "docs/VISUAL_BRIDGE_TECHNICAL_PREFLIGHT_B1_B2_2026_06_01.md",
        "role": "B1 technical preflight",
    },
    {
        "path": "docs/VISUAL_PRIOR_PREMIUM_DECKS_AND_AGENTIC_WORKFLOW_REVIEW_2026_06_01.md",
        "role": "reference-deck calibration review",
    },
]

OVERLAY_TEXT = [
    {
        "id": "headline",
        "role": "decision_headline",
        "content": "Public space omics needs an evaluation layer.",
        "x": 0.065,
        "y": 0.118,
        "font_pt": 27,
        "color": "ink",
        "max_lines": 1,
        "z": "Z4",
    },
    {
        "id": "support",
        "role": "supporting_claim",
        "content": "SpaceBio-Bench organizes public studies into traceable tasks and audited scores.",
        "x": 0.066,
        "y": 0.190,
        "font_pt": 11.2,
        "color": "muted",
        "max_lines": 1,
        "z": "Z4",
    },
    {
        "id": "public_studies",
        "role": "primary_callout",
        "content": "Public studies",
        "x": 0.104,
        "y": 0.337,
        "font_pt": 10.5,
        "color": "source_blue",
        "z": "Z4",
    },
    {
        "id": "alignment_layer",
        "role": "primary_callout",
        "content": "Missions, tissues, samples, labels",
        "x": 0.366,
        "y": 0.337,
        "font_pt": 10.0,
        "color": "bio_green",
        "z": "Z4",
    },
    {
        "id": "benchmark_tasks",
        "role": "primary_callout",
        "content": "Benchmark tasks",
        "x": 0.716,
        "y": 0.337,
        "font_pt": 10.5,
        "color": "model_purple",
        "z": "Z4",
    },
    {
        "id": "audited_scores",
        "role": "primary_callout",
        "content": "Audited scores",
        "x": 0.858,
        "y": 0.337,
        "font_pt": 10.5,
        "color": "assay_teal",
        "z": "Z4",
    },
]

STATUS_LABELS = [
    {
        "id": "trace_status",
        "role": "claim_boundary",
        "content": "public-data resource; source-traceable",
        "x": 0.707,
        "y": 0.866,
        "font_pt": 7.8,
        "color": "muted",
        "z": "Z6",
    },
    {
        "id": "source_note",
        "role": "source_note",
        "content": "Source: project inventories and source records.",
        "x": 0.065,
        "y": 0.900,
        "font_pt": 7.2,
        "color": "muted",
        "z": "Z6",
    },
]

FOCUS_MARKS = [
    {
        "id": "source_to_benchmark_path",
        "role": "dominant_motion",
        "shape": "converging_transfer_path",
        "x0": 0.145,
        "x1": 0.900,
        "y": 0.510,
        "color": "source_blue",
        "z": "Z7",
    }
]


def rel(path: Path) -> str:
    return str(path.relative_to(ROOT))


def write_json(path: Path, data: Any) -> None:
    path.write_text(json.dumps(data, indent=2) + "\n", encoding="utf-8")


def load_tokens() -> dict[str, Any]:
    return json.loads(TOKEN_PATH.read_text(encoding="utf-8"))


def count_words(items: list[dict[str, Any]]) -> int:
    words = 0
    for item in items:
        content = str(item.get("content", ""))
        normalized = content.replace("/", " ").replace(";", " ").replace("+", " ")
        words += len([token for token in normalized.split() if token])
    return words


def y_from_slide(value: float) -> float:
    return 1.0 - value


def make_figure() -> tuple[plt.Figure, plt.Axes]:
    fig = plt.figure(figsize=(16, 9), dpi=SLIDE_DPI)
    ax = fig.add_axes([0, 0, 1, 1])
    ax.set_axis_off()
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)
    return fig, ax


def color(tokens: dict[str, Any], name: str) -> str:
    return tokens["colors"].get(name, name)


def rect(
    ax: plt.Axes,
    tokens: dict[str, Any],
    x: float,
    y: float,
    w: float,
    h: float,
    *,
    face: str,
    edge: str | None = None,
    alpha: float = 1.0,
    lw: float = 0.8,
    z: float = 3,
) -> None:
    ax.add_patch(
        Rectangle(
            (x, y),
            w,
            h,
            transform=ax.transAxes,
            facecolor=color(tokens, face),
            edgecolor=color(tokens, edge) if edge else "none",
            linewidth=lw,
            alpha=alpha,
            zorder=z,
        )
    )


def circle(
    ax: plt.Axes,
    tokens: dict[str, Any],
    x: float,
    y: float,
    r: float,
    *,
    face: str,
    edge: str | None = None,
    alpha: float = 1.0,
    lw: float = 1.0,
    z: float = 5,
) -> None:
    ax.add_patch(
        Circle(
            (x, y),
            r,
            transform=ax.transAxes,
            facecolor=color(tokens, face),
            edgecolor=color(tokens, edge) if edge else "none",
            linewidth=lw,
            alpha=alpha,
            zorder=z,
        )
    )


def shadow_rect(
    ax: plt.Axes,
    tokens: dict[str, Any],
    x: float,
    y: float,
    w: float,
    h: float,
    *,
    alpha: float = 0.10,
    z: float = 2,
) -> None:
    rect(ax, tokens, x + 0.009, y - 0.010, w, h, face="shadow", edge=None, alpha=alpha, z=z)


def arrow(
    ax: plt.Axes,
    tokens: dict[str, Any],
    start: tuple[float, float],
    end: tuple[float, float],
    *,
    arrowstyle: str = "-|>",
    color_name: str = "muted",
    alpha: float = 0.50,
    lw: float = 1.05,
    z: float = 8,
    mutation_scale: float = 12,
) -> None:
    ax.add_patch(
        FancyArrowPatch(
            start,
            end,
            arrowstyle=arrowstyle,
            mutation_scale=mutation_scale,
            linewidth=lw,
            color=color(tokens, color_name),
            alpha=alpha,
            transform=ax.transAxes,
            zorder=z,
            connectionstyle="arc3,rad=0.0",
        )
    )


def draw_canvas(ax: plt.Axes, tokens: dict[str, Any]) -> None:
    h, w = 720, 1280
    yy, xx = np.mgrid[0:h, 0:w]
    rng = np.random.default_rng(20260601)
    warm = np.array([0.955, 0.940, 0.910])
    cool = np.array([0.925, 0.955, 0.965])
    gradient = (xx / w)[..., None]
    base = warm * (1.0 - gradient * 0.32) + cool * (gradient * 0.32)
    grain = rng.normal(0, 0.0042, size=(h, w, 1))
    vignette = ((xx / w - 0.52) ** 2 + (yy / h - 0.46) ** 2)[..., None] * np.array([0.043, 0.035, 0.024])
    texture = np.clip(base + grain - vignette, 0, 1)
    ax.imshow(texture, extent=(0, 1, 0, 1), origin="lower", zorder=0)
    ax.set_aspect("auto")

    for y in [0.190, 0.328, 0.500, 0.672, 0.835]:
        ax.plot([0.055, 0.945], [y, y], color=color(tokens, "rule"), alpha=0.18, linewidth=0.75, transform=ax.transAxes, zorder=1)
    for x in [0.075, 0.275, 0.500, 0.725, 0.925]:
        ax.plot([x, x], [0.095, 0.870], color=color(tokens, "rule"), alpha=0.10, linewidth=0.70, transform=ax.transAxes, zorder=1)


def draw_source_mark(ax: plt.Axes, tokens: dict[str, Any], x: float, y: float, size: float, accent: str, *, z: float) -> None:
    shadow_rect(ax, tokens, x - size * 0.48, y - size * 0.55, size * 0.95, size * 1.22, alpha=0.045, z=z - 1)
    rect(
        ax,
        tokens,
        x - size * 0.52,
        y - size * 0.50,
        size * 0.95,
        size * 1.18,
        face="#FDFDFB",
        edge="rule",
        alpha=0.74,
        lw=0.45,
        z=z,
    )
    ax.plot(
        [x - size * 0.33, x + size * 0.24],
        [y + size * 0.26, y + size * 0.26],
        color=color(tokens, accent),
        alpha=0.55,
        linewidth=0.65,
        transform=ax.transAxes,
        zorder=z + 1,
    )
    for offset in [0.08, -0.09, -0.26]:
        ax.plot(
            [x - size * 0.33, x + size * 0.28],
            [y + size * offset, y + size * offset],
            color=color(tokens, "rule"),
            alpha=0.42,
            linewidth=0.50,
            transform=ax.transAxes,
            zorder=z + 1,
        )


def draw_fragmented_sources(ax: plt.Axes, tokens: dict[str, Any]) -> list[tuple[float, float]]:
    rng = np.random.default_rng(20260602)
    accents = ["source_blue", "bio_green", "assay_teal", "label_amber"]
    points: list[tuple[float, float]] = []
    for idx in range(22):
        x = float(rng.uniform(0.095, 0.282))
        y = float(rng.uniform(0.385, 0.642))
        size = float(rng.uniform(0.016, 0.024))
        accent = accents[idx % len(accents)]
        draw_source_mark(ax, tokens, x, y, size, accent, z=5 + idx * 0.01)
        points.append((x, y))
    for x, y in points[::3]:
        circle(ax, tokens, x + 0.018, y - 0.021, 0.0048, face="source_blue", edge=None, alpha=0.35, z=6)
    return points


def draw_alignment_rail(ax: plt.Axes, tokens: dict[str, Any]) -> None:
    rail_y = 0.505
    ax.plot([0.330, 0.680], [rail_y, rail_y], color=color(tokens, "rule"), alpha=0.75, linewidth=1.25, transform=ax.transAxes, zorder=5)
    tick_specs = [
        (0.365, "bio_green", 0.022),
        (0.455, "assay_teal", 0.019),
        (0.545, "source_blue", 0.019),
        (0.635, "label_amber", 0.022),
    ]
    for idx, (x, accent, radius) in enumerate(tick_specs):
        ax.plot([x, x], [rail_y - 0.058, rail_y + 0.058], color=color(tokens, "rule"), alpha=0.20, linewidth=0.85, transform=ax.transAxes, zorder=4)
        circle(ax, tokens, x, rail_y, radius, face="#FDFDFB", edge=accent, alpha=0.96, lw=1.05, z=7)
        circle(ax, tokens, x, rail_y, radius * 0.38, face=accent, edge=None, alpha=0.70, z=8)
        if idx in {0, 3}:
            circle(ax, tokens, x, rail_y, radius * 1.85, face="#FFFFFF", edge=accent, alpha=0.16, lw=0.90, z=6)


def draw_task_surface(ax: plt.Axes, tokens: dict[str, Any]) -> None:
    x, y, w, h = 0.714, 0.396, 0.125, 0.235
    shadow_rect(ax, tokens, x, y, w, h, alpha=0.13, z=5)
    rect(ax, tokens, x, y, w, h, face="#FDFDFB", edge="rule", alpha=0.98, lw=0.70, z=7)
    ax.plot([x + 0.018, x + w - 0.018], [y + h - 0.047, y + h - 0.047], color=color(tokens, "model_purple"), alpha=0.48, linewidth=1.0, transform=ax.transAxes, zorder=8)
    for line_y in [y + 0.158, y + 0.123, y + 0.088, y + 0.053]:
        ax.plot([x + 0.020, x + w - 0.020], [line_y, line_y], color=color(tokens, "rule"), alpha=0.55, linewidth=0.68, transform=ax.transAxes, zorder=8)
    for idx, accent in enumerate(["bio_green", "source_blue", "label_amber", "assay_teal"]):
        circle(ax, tokens, x + 0.028 + idx * 0.024, y + 0.026, 0.006, face=accent, edge=None, alpha=0.82, z=8)


def draw_score_ledger(ax: plt.Axes, tokens: dict[str, Any]) -> None:
    x, y, w, h = 0.858, 0.418, 0.074, 0.190
    shadow_rect(ax, tokens, x, y, w, h, alpha=0.12, z=5)
    rect(ax, tokens, x, y, w, h, face="#FCFCFA", edge="rule", alpha=0.97, lw=0.68, z=7)
    bar_xs = [x + 0.022, x + 0.039, x + 0.056]
    bar_heights = [0.072, 0.103, 0.127]
    for bar_x, bar_h, accent in zip(bar_xs, bar_heights, ["source_blue", "bio_green", "assay_teal"]):
        ax.plot([bar_x, bar_x], [y + 0.035, y + 0.035 + bar_h], color=color(tokens, accent), alpha=0.70, linewidth=3.0, transform=ax.transAxes, zorder=8)
    circle(ax, tokens, x + 0.038, y + h - 0.038, 0.020, face="#FDFDFB", edge="assay_teal", alpha=0.96, lw=1.0, z=8)
    ax.plot([x + 0.029, x + 0.036, x + 0.050], [y + h - 0.038, y + h - 0.047, y + h - 0.026], color=color(tokens, "assay_teal"), alpha=0.78, linewidth=1.15, transform=ax.transAxes, zorder=9)
    ax.plot([x + 0.017, x + w - 0.017], [y + 0.026, y + 0.026], color=color(tokens, "rule"), alpha=0.35, linewidth=0.6, transform=ax.transAxes, zorder=8)


def draw_motion(ax: plt.Axes, tokens: dict[str, Any], source_points: list[tuple[float, float]]) -> None:
    rail_targets = [(0.365, 0.505), (0.455, 0.505), (0.545, 0.505), (0.635, 0.505)]
    for idx, point in enumerate(source_points[0:12:2]):
        target = rail_targets[idx % len(rail_targets)]
        arrow(ax, tokens, point, target, arrowstyle="-", color_name="source_blue", alpha=0.13, lw=0.70, z=4, mutation_scale=1)
    arrow(ax, tokens, (0.286, 0.515), (0.342, 0.507), color_name="source_blue", alpha=0.42, lw=1.05, z=8)
    arrow(ax, tokens, (0.665, 0.505), (0.706, 0.505), color_name="muted", alpha=0.52, lw=1.05, z=8)
    arrow(ax, tokens, (0.840, 0.505), (0.858, 0.505), color_name="assay_teal", alpha=0.40, lw=1.0, z=8)
    ax.plot(
        [0.128, 0.245, 0.365, 0.545, 0.705, 0.898],
        [0.616, 0.585, 0.535, 0.508, 0.506, 0.513],
        color=color(tokens, "source_blue"),
        alpha=0.18,
        linewidth=6.0,
        solid_capstyle="round",
        transform=ax.transAxes,
        zorder=3,
    )


def draw_scene(ax: plt.Axes, tokens: dict[str, Any]) -> None:
    draw_canvas(ax, tokens)
    rect(ax, tokens, 0.058, 0.305, 0.884, 0.395, face="#FBFAF6", edge=None, alpha=0.47, z=2)
    rect(ax, tokens, 0.078, 0.345, 0.220, 0.300, face="#F2F5F4", edge=None, alpha=0.48, z=2.2)
    rect(ax, tokens, 0.350, 0.415, 0.300, 0.180, face="#EEF5F3", edge=None, alpha=0.34, z=2.2)
    rect(ax, tokens, 0.694, 0.350, 0.246, 0.310, face="#F6F2FA", edge=None, alpha=0.29, z=2.2)
    ax.plot([0.080, 0.920], [0.505, 0.505], color=color(tokens, "rule"), alpha=0.23, linewidth=1.0, transform=ax.transAxes, zorder=3)

    source_points = draw_fragmented_sources(ax, tokens)
    draw_alignment_rail(ax, tokens)
    draw_task_surface(ax, tokens)
    draw_score_ledger(ax, tokens)
    draw_motion(ax, tokens, source_points)

    for x in [0.328, 0.688, 0.848]:
        ax.plot([x, x], [0.352, 0.648], color=color(tokens, "rule"), alpha=0.15, linewidth=0.75, transform=ax.transAxes, zorder=3)


def build_overlay() -> dict[str, Any]:
    text_items = OVERLAY_TEXT + STATUS_LABELS
    return {
        "slide_id": SLIDE_ID,
        "stage": "pre_render",
        "canvas": CANVAS,
        "coordinate_system": "normalized_0_1_from_top_left",
        "text": OVERLAY_TEXT,
        "status_labels": STATUS_LABELS,
        "focus_marks": FOCUS_MARKS,
        "visible_word_count": count_words(text_items),
        "visible_word_budget": 45,
        "forbidden_visible_terms": FORBIDDEN_VISIBLE_TERMS,
    }


def render_overlay(ax: plt.Axes, tokens: dict[str, Any], overlay: dict[str, Any]) -> None:
    for item in list(overlay["text"]) + list(overlay["status_labels"]):
        role = str(item.get("role", ""))
        weight = "bold" if role in {"decision_headline", "primary_callout", "claim_boundary"} else "normal"
        va = "top" if role in {"decision_headline", "supporting_claim"} else "center"
        ax.text(
            float(item["x"]),
            y_from_slide(float(item["y"])),
            str(item["content"]),
            color=color(tokens, str(item.get("color", "ink"))),
            fontsize=float(item.get("font_pt", 8.0)),
            fontweight=weight,
            ha="left",
            va=va,
            transform=ax.transAxes,
            zorder=20,
        )


def output_paths() -> dict[str, str]:
    return {
        "scene_plate": rel(SLIDE_DIR / "scene_plate.png"),
        "rendered_preview_png": rel(SLIDE_DIR / "rendered_preview.png"),
        "rendered_preview_pdf": rel(SLIDE_DIR / "rendered_preview.pdf"),
        "overlay_spec": rel(SLIDE_DIR / "overlay_spec.json"),
        "manifest": rel(SLIDE_DIR / "manifest.json"),
        "qa": rel(SLIDE_DIR / "qa.json"),
        "qa_zoom_evaluation_surface": rel(SLIDE_DIR / "qa_zoom_evaluation_surface.png"),
        "qa_zoom_status_scores": rel(SLIDE_DIR / "qa_zoom_status_scores.png"),
        "inspection_sheet": rel(OUT_ROOT / "b1_evaluation_layer_inspection_sheet.png"),
        "inspection_sheet_manifest": rel(OUT_ROOT / "b1_evaluation_layer_inspection_sheet_manifest.json"),
    }


def scan_forbidden_terms(overlay: dict[str, Any]) -> list[str]:
    visible = "\n".join(
        [item["content"] for item in overlay["text"]]
        + [item["content"] for item in overlay["status_labels"]]
    ).lower()
    return [term for term in FORBIDDEN_VISIBLE_TERMS if term.lower() in visible]


def evidence_presence() -> list[dict[str, Any]]:
    records = []
    for item in EVIDENCE_SOURCES:
        path = ROOT / item["path"]
        records.append({**item, "exists": path.exists()})
    return records


def image_dimensions(path: Path) -> list[int] | None:
    if Image is None:
        return None
    with Image.open(path) as img:
        return [int(img.width), int(img.height)]


def crop_image(source: Path, target: Path, box_norm_top_left: tuple[float, float, float, float]) -> None:
    if Image is None:
        return
    with Image.open(source) as img:
        w, h = img.size
        x0, y0, x1, y1 = box_norm_top_left
        crop = img.crop((int(x0 * w), int(y0 * h), int(x1 * w), int(y1 * h)))
        crop.save(target)


def build_inspection_sheet(tokens: dict[str, Any], paths: dict[str, str]) -> None:
    preview = ROOT / paths["rendered_preview_png"]
    zoom_surface = ROOT / paths["qa_zoom_evaluation_surface"]
    zoom_status = ROOT / paths["qa_zoom_status_scores"]
    sheet = ROOT / paths["inspection_sheet"]

    fig = plt.figure(figsize=(16.8, 6.6), dpi=220)
    grid = fig.add_gridspec(1, 3, left=0.025, right=0.985, bottom=0.075, top=0.805, wspace=0.035)
    fig.suptitle(
        "B1 brand stress-test inspection: fragmentation to evaluation layer",
        x=0.025,
        y=0.948,
        ha="left",
        fontsize=12.5,
        fontweight="bold",
        color=color(tokens, "ink"),
    )
    fig.text(
        0.025,
        0.880,
        "Gate: one movement, source-traceable status, no internal terms, no card-box layout.",
        fontsize=8.0,
        color=color(tokens, "muted"),
        ha="left",
    )
    panels = [
        (preview, "full slide"),
        (zoom_surface, "evaluation surface zoom"),
        (zoom_status, "status and score zoom"),
    ]
    for idx, (source, title) in enumerate(panels):
        ax = fig.add_subplot(grid[0, idx])
        ax.imshow(mpimg.imread(source))
        ax.set_title(title, loc="left", fontsize=7.5, pad=4, color=color(tokens, "ink"))
        ax.axis("off")
    fig.savefig(sheet, dpi=220, facecolor="white")
    plt.close(fig)

    write_json(
        ROOT / paths["inspection_sheet_manifest"],
        {
            "slide_id": SLIDE_ID,
            "created": CREATED,
            "inspection_sheet": paths["inspection_sheet"],
            "panels": [
                {"role": "full slide", "path": paths["rendered_preview_png"]},
                {"role": "evaluation surface zoom", "path": paths["qa_zoom_evaluation_surface"]},
                {"role": "status and score zoom", "path": paths["qa_zoom_status_scores"]},
            ],
        },
    )


def build_manifest(overlay: dict[str, Any], paths: dict[str, str]) -> dict[str, Any]:
    return {
        "slide_id": SLIDE_ID,
        "created": CREATED,
        "stage": "pre_render",
        "slide_type": "thesis",
        "visual_identity_token_version": "0.1",
        "style_direction": "hybrid",
        "reference_calibration_role": {
            "primary": "kmu_proof_stage",
            "secondary": ["pbs_benchmark_concept", "isgp_live_rhythm"],
            "rationale": "B1 must explain an operating model to a first-time viewer while opening the live deck spine.",
        },
        "decision_headline": "Public space omics is valuable but fragmented; SpaceBio-Bench adds an auditable evaluation layer.",
        "audience_question": "Why is this project more than another reanalysis of public GeneLab/OSDR data?",
        "visual_move": "scattered public-study marks compress into one horizontal evaluation surface",
        "claim_boundary": "conceptual bridge slide; not a quantitative result claim",
        "content_brief": "docs/VISUAL_BRIDGE_CONTENT_BRIEFS_B1_B4_2026_06_01.md",
        "technical_preflight": "docs/VISUAL_BRIDGE_TECHNICAL_PREFLIGHT_B1_B2_2026_06_01.md",
        "technical_protocol": "docs/VISUAL_TECHNICAL_PRODUCTION_PROTOCOL_2026_06_01.md",
        "evidence_sources": evidence_presence(),
        "semantic_color_roles_used": [
            "source_blue: public source marks and transfer path",
            "bio_green: alignment/context",
            "assay_teal: audited score check",
            "label_amber: label tick",
            "model_purple: benchmark task surface",
            "muted/rule: infrastructure",
        ],
        "z_layers_used": [
            "Z0 canvas_atmosphere",
            "Z1 measurement_infrastructure",
            "Z2 evidence_surface",
            "Z3 proof_objects",
            "Z4 interpretation_overlay",
            "Z6 status_caveat",
            "Z7 motion_focus",
        ],
        "depth_tokens_used": ["recessed_band", "proof_object", "dominant_evidence", "flat"],
        "anti_pattern_checks": {
            "no_card_box_layout": True,
            "no_dashboard_tiles": True,
            "no_nested_cards": True,
            "no_decorative_space_wallpaper": True,
            "no_raw_local_paths_visible": True,
            "one_dominant_motion_path": True,
        },
        "visible_terms": [item["content"] for item in overlay["text"] + overlay["status_labels"]],
        "terms_pushed_to_speaker_notes": [
            "GeneLab/OSDR",
            "mission-held-out benchmark contract",
            "source records",
            "task definitions",
            "metrics",
            "audit trails",
        ],
        "outputs": paths,
        "generator": "scripts/build_b1_brand_stress_test_scene.py",
        "random_seeds": {"canvas": 20260601, "source_marks": 20260602},
        "qa": {"stage": "pre_render", "pre_render_gate_expected": True, "post_render_required_before_use": True},
    }


def build_qa(overlay: dict[str, Any], paths: dict[str, str], preview_path: Path) -> dict[str, Any]:
    forbidden_hits = scan_forbidden_terms(overlay)
    dims = image_dimensions(preview_path)
    expected_dims = [CANVAS["width_px"], CANVAS["height_px"]]
    evidence = evidence_presence()
    return {
        "slide_id": SLIDE_ID,
        "stage": "post_render",
        "created": CREATED,
        "pre_render_gate": {
            "content_brief_declared": True,
            "technical_preflight_declared": True,
            "technical_protocol_declared": True,
            "evidence_sources_declared": True,
            "all_evidence_sources_exist": all(item["exists"] for item in evidence),
            "claim_boundary_declared": True,
            "visible_text_word_count": overlay["visible_word_count"],
            "visible_text_budget": overlay["visible_word_budget"],
            "forbidden_visible_terms_absent": len(forbidden_hits) == 0,
            "forbidden_visible_term_hits": forbidden_hits,
            "output_paths_declared": True,
            "overlay_spec_declared": True,
            "manifest_declared": True,
            "dominant_visual_move_declared": True,
        },
        "post_render_gate": {
            "rendered_outputs": {
                "scene_plate": paths["scene_plate"],
                "rendered_preview_png": paths["rendered_preview_png"],
                "rendered_preview_pdf": paths["rendered_preview_pdf"],
                "inspection_sheet": paths["inspection_sheet"],
            },
            "image_dimensions": {"rendered_preview_png": dims, "expected": expected_dims},
            "expected_dimensions_match": dims == expected_dims,
            "manual_full_size_inspection": "pending",
            "manual_zoom_inspection": "pending",
            "manual_thumbnail_inspection": "pending",
            "text_overlap_check": "pending",
            "jargon_check": "pending",
            "internal_language_check": "pending",
            "premium_quality_gate": "pending",
            "visual_verdict": "brand stress-test render; awaiting manual inspection",
        },
    }


def render() -> dict[str, Any]:
    tokens = load_tokens()
    SLIDE_DIR.mkdir(parents=True, exist_ok=True)
    paths = output_paths()
    overlay = build_overlay()
    manifest = build_manifest(overlay, paths)

    scene_plate = ROOT / paths["scene_plate"]
    preview_png = ROOT / paths["rendered_preview_png"]
    preview_pdf = ROOT / paths["rendered_preview_pdf"]

    fig, ax = make_figure()
    draw_scene(ax, tokens)
    fig.savefig(scene_plate, dpi=SLIDE_DPI, facecolor=color(tokens, "canvas_warm"))
    plt.close(fig)

    fig, ax = make_figure()
    draw_scene(ax, tokens)
    render_overlay(ax, tokens, overlay)
    fig.savefig(preview_png, dpi=SLIDE_DPI, facecolor=color(tokens, "canvas_warm"))
    fig.savefig(preview_pdf, facecolor=color(tokens, "canvas_warm"))
    plt.close(fig)

    crop_image(preview_png, ROOT / paths["qa_zoom_evaluation_surface"], (0.055, 0.285, 0.945, 0.720))
    crop_image(preview_png, ROOT / paths["qa_zoom_status_scores"], (0.680, 0.315, 0.950, 0.910))
    build_inspection_sheet(tokens, paths)

    overlay["stage"] = "post_render"
    manifest["stage"] = "post_render"
    manifest["renderer"] = "scripts/build_b1_brand_stress_test_scene.py"
    manifest["qa"]["stage"] = "post_render"
    manifest["qa"]["render_outputs_exist"] = True
    qa = build_qa(overlay, paths, preview_png)

    write_json(ROOT / paths["overlay_spec"], overlay)
    write_json(ROOT / paths["manifest"], manifest)
    write_json(ROOT / paths["qa"], qa)
    return {"rendered": paths, "qa": qa}


def main() -> None:
    print(json.dumps(render(), indent=2))


if __name__ == "__main__":
    main()
