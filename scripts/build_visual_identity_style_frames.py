#!/usr/bin/env python3
"""Render visual identity style-frame candidates for SpaceBio-Bench."""

from __future__ import annotations

import argparse
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


TOKEN_PATH = ROOT / "config" / "visual_identity" / "spacebio_bench_visual_identity_v0_1.json"
OUT_ROOT = ROOT / "output" / "visual_identity_style_frames"
CREATED = "2026-06-01"


def load_tokens() -> dict[str, Any]:
    return json.loads(TOKEN_PATH.read_text(encoding="utf-8"))


def rel(path: Path) -> str:
    return str(path.relative_to(ROOT))


def write_json(path: Path, data: Any) -> None:
    path.write_text(json.dumps(data, indent=2) + "\n", encoding="utf-8")


def y_from_slide(value: float) -> float:
    return 1.0 - value


def make_figure(tokens: dict[str, Any]) -> tuple[plt.Figure, plt.Axes]:
    canvas = tokens["layout"]["canvas"]
    fig = plt.figure(figsize=(16, 9), dpi=int(canvas["dpi"]))
    ax = fig.add_axes([0, 0, 1, 1])
    ax.set_axis_off()
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)
    return fig, ax


def count_words(items: list[dict[str, Any]]) -> int:
    return sum(len(str(item.get("content", "")).replace("/", " ").split()) for item in items)


def rect(
    ax: plt.Axes,
    colors: dict[str, str],
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
            facecolor=colors.get(face, face),
            edgecolor=colors.get(edge, edge) if edge else "none",
            alpha=alpha,
            linewidth=lw,
            zorder=z,
        )
    )


def circle(
    ax: plt.Axes,
    colors: dict[str, str],
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
            facecolor=colors.get(face, face),
            edgecolor=colors.get(edge, edge) if edge else "none",
            alpha=alpha,
            linewidth=lw,
            zorder=z,
        )
    )


def arrow(
    ax: plt.Axes,
    colors: dict[str, str],
    start: tuple[float, float],
    end: tuple[float, float],
    *,
    color: str,
    alpha: float = 0.7,
    lw: float = 1.2,
    z: float = 8,
) -> None:
    ax.add_patch(
        FancyArrowPatch(
            start,
            end,
            arrowstyle="-|>",
            mutation_scale=13,
            linewidth=lw,
            color=colors[color],
            alpha=alpha,
            transform=ax.transAxes,
            zorder=z,
        )
    )


def draw_texture(ax: plt.Axes, colors: dict[str, str], *, seed: int, warm: bool = True, dark: bool = False) -> None:
    h, w = 720, 1280
    yy, xx = np.mgrid[0:h, 0:w]
    rng = np.random.default_rng(seed)
    if dark:
        base = np.zeros((h, w, 3), dtype=float)
        base[..., 0] = 0.055
        base[..., 1] = 0.135
        base[..., 2] = 0.190
        grain = rng.normal(0, 0.006, size=(h, w, 1))
        glow = np.exp(-(((xx / w - 0.72) ** 2) / 0.090 + ((yy / h - 0.48) ** 2) / 0.180))[..., None]
        texture = np.clip(base + grain + glow * np.array([0.018, 0.052, 0.074]), 0, 1)
    else:
        base = np.zeros((h, w, 3), dtype=float)
        if warm:
            base[..., 0] = 0.956
            base[..., 1] = 0.942
            base[..., 2] = 0.910
        else:
            base[..., 0] = 0.945
            base[..., 1] = 0.958
            base[..., 2] = 0.962
        grain = rng.normal(0, 0.004, size=(h, w, 1))
        vignette = ((xx / w - 0.55) ** 2 + (yy / h - 0.46) ** 2)[..., None] * np.array([0.040, 0.034, 0.024])
        texture = np.clip(base + grain - vignette, 0, 1)
    ax.imshow(texture, extent=(0, 1, 0, 1), origin="lower", zorder=0)
    ax.set_aspect("auto")


def draw_grid(ax: plt.Axes, colors: dict[str, str], *, dark: bool = False, alpha: float = 0.14) -> None:
    rule_color = "#B9D4E5" if dark else colors["rule"]
    for y in [0.190, 0.335, 0.500, 0.665, 0.835]:
        ax.plot([0.055, 0.945], [y, y], color=rule_color, alpha=alpha, linewidth=0.75, transform=ax.transAxes, zorder=1)
    for x in [0.075, 0.275, 0.500, 0.725, 0.925]:
        ax.plot([x, x], [0.095, 0.870], color=rule_color, alpha=alpha * 0.75, linewidth=0.70, transform=ax.transAxes, zorder=1)


def draw_source_stack(ax: plt.Axes, colors: dict[str, str], x: float, y: float, *, dark: bool = False, z: float = 5) -> None:
    paper = "#FDFDFB" if not dark else "#DDEAF1"
    edge = colors["rule"] if not dark else "#8DB5CA"
    shadow = colors["shadow"] if not dark else "#020A10"
    for idx, (dx, dy, alpha) in enumerate([(0.030, 0.026, 0.22), (0.015, 0.013, 0.44), (0.0, 0.0, 0.95)]):
        rect(ax, {**colors, "paper_local": paper, "edge_local": edge, "shadow_local": shadow}, x + dx + 0.008, y + dy - 0.010, 0.142, 0.190, face="shadow_local", alpha=0.10 * alpha, z=z - 1)
        rect(ax, {**colors, "paper_local": paper, "edge_local": edge}, x + dx, y + dy, 0.142, 0.190, face="paper_local", edge="edge_local", alpha=alpha, lw=0.65, z=z + idx * 0.03)
    line = "#A8B3BF" if not dark else "#719BAF"
    for yy in [y + 0.145, y + 0.118, y + 0.091, y + 0.064]:
        ax.plot([x + 0.023, x + 0.115], [yy, yy], color=line, alpha=0.60, linewidth=0.8, transform=ax.transAxes, zorder=z + 1)


def draw_task_surface(ax: plt.Axes, colors: dict[str, str], x: float, y: float, *, dark: bool = False, z: float = 6) -> None:
    paper = "#FFFFFF" if not dark else "#DDEAF1"
    edge = colors["rule"] if not dark else "#8DB5CA"
    shadow = colors["shadow"] if not dark else "#020A10"
    rect(ax, {**colors, "shadow_local": shadow}, x + 0.011, y - 0.012, 0.150, 0.245, face="shadow_local", alpha=0.14, z=z - 1)
    rect(ax, {**colors, "paper_local": paper, "edge_local": edge}, x, y, 0.150, 0.245, face="paper_local", edge="edge_local", alpha=0.98, lw=0.75, z=z)
    for yy in [y + 0.195, y + 0.158, y + 0.121, y + 0.084, y + 0.047]:
        ax.plot([x + 0.025, x + 0.125], [yy, yy], color=edge, alpha=0.55, linewidth=0.75, transform=ax.transAxes, zorder=z + 1)
    for xx, color in [(x + 0.035, "source_blue"), (x + 0.063, "bio_green"), (x + 0.091, "label_amber"), (x + 0.119, "assay_teal")]:
        circle(ax, colors, xx, y + 0.025, 0.0065, face=color, alpha=0.86, z=z + 1)


def draw_mission_gate(ax: plt.Axes, colors: dict[str, str], x: float, *, dark: bool = False) -> None:
    red = colors["test_red"]
    ax.plot([x, x], [0.335, 0.690], color=red, alpha=0.78, linewidth=1.35, transform=ax.transAxes, zorder=9)
    ax.plot([x + 0.010, x + 0.010], [0.360, 0.665], color=red, alpha=0.20, linewidth=0.85, transform=ax.transAxes, zorder=8)
    label_color = "#F9CACA" if dark else red
    ax.text(x + 0.025, 0.665, "held-out mission", color=label_color, fontsize=8.3, fontweight="bold", ha="left", va="center", transform=ax.transAxes, zorder=10)
    circle(ax, colors, x + 0.102, 0.520, 0.033, face="paper" if not dark else "#102E40", edge="test_red", alpha=0.98, lw=1.2, z=10)
    ax.text(x + 0.102, 0.520, "M4", color=label_color, fontsize=6.8, fontweight="bold", ha="center", va="center", transform=ax.transAxes, zorder=11)


def draw_overlay(ax: plt.Axes, tokens: dict[str, Any], style: dict[str, Any]) -> None:
    colors = tokens["colors"]
    overlay = style["overlay"]
    for item in overlay["text"] + overlay["status_labels"]:
        color = item.get("color", "ink")
        value = style.get("dark_text", {}).get(color, colors.get(color, color)) if style.get("dark") else colors.get(color, color)
        role = str(item.get("role", ""))
        weight = "bold" if role in {"headline", "callout", "status_strong"} else "normal"
        va = "top" if role == "headline" else "center"
        ax.text(
            float(item["x"]),
            y_from_slide(float(item["y"])),
            str(item["content"]),
            color=value,
            fontsize=float(item["font_pt"]),
            fontweight=weight,
            ha="left",
            va=va,
            transform=ax.transAxes,
            zorder=20,
        )


def draw_editorial(ax: plt.Axes, tokens: dict[str, Any], *, overlay: bool) -> None:
    colors = tokens["colors"]
    ax.set_facecolor("#FFFFFF")
    rect(ax, colors, 0, 0, 1, 1, face="#FFFFFF", z=0)
    draw_grid(ax, colors, alpha=0.10)
    ax.plot([0.065, 0.925], [0.305, 0.305], color=colors["ink"], alpha=0.18, linewidth=0.8, transform=ax.transAxes, zorder=2)
    ax.plot([0.065, 0.925], [0.745, 0.745], color=colors["ink"], alpha=0.12, linewidth=0.8, transform=ax.transAxes, zorder=2)
    draw_source_stack(ax, colors, 0.105, 0.435, z=5)
    ax.plot([0.270, 0.755], [0.530, 0.530], color=colors["rule"], alpha=0.80, linewidth=1.1, transform=ax.transAxes, zorder=4)
    for x, label, color in [(0.335, "mission", "bio_green"), (0.460, "samples", "source_blue"), (0.585, "labels", "label_amber")]:
        circle(ax, colors, x, 0.530, 0.038, face="paper", edge=color, alpha=0.98, lw=1.20, z=6)
        ax.text(x, 0.530, label, color=colors[color], fontsize=6.7, ha="center", va="center", fontweight="bold", transform=ax.transAxes, zorder=7)
    draw_task_surface(ax, colors, 0.745, 0.407, z=7)
    arrow(ax, colors, (0.255, 0.530), (0.296, 0.530), color="muted", alpha=0.55, lw=1.0)
    arrow(ax, colors, (0.625, 0.530), (0.735, 0.530), color="muted", alpha=0.55, lw=1.0)
    draw_mission_gate(ax, colors, 0.855)
    if overlay:
        draw_overlay(ax, tokens, STYLE_FRAMES["editorial"])


def draw_mission_system(ax: plt.Axes, tokens: dict[str, Any], *, overlay: bool) -> None:
    colors = tokens["colors"]
    draw_texture(ax, colors, seed=20260631, dark=True)
    draw_grid(ax, colors, dark=True, alpha=0.20)
    rect(ax, {**colors, "mission_surface": "#0C2635"}, 0.062, 0.320, 0.876, 0.385, face="mission_surface", edge=None, alpha=0.72, z=2)
    ax.plot([0.100, 0.840], [0.515, 0.515], color="#8AD4E8", alpha=0.33, linewidth=1.1, transform=ax.transAxes, zorder=3)
    theta = np.linspace(0.0, np.pi, 160)
    ax.plot(0.465 + np.cos(theta) * 0.330, 0.515 + np.sin(theta) * 0.155, color="#8AD4E8", alpha=0.24, linewidth=1.1, transform=ax.transAxes, zorder=3)
    draw_source_stack(ax, colors, 0.105, 0.430, dark=True, z=5)
    for x, label in [(0.340, "M1"), (0.435, "M2"), (0.530, "M3")]:
        circle(ax, colors, x, 0.515, 0.034, face="#102E40", edge="bio_green", alpha=0.98, lw=1.2, z=6)
        ax.text(x, 0.515, label, color="#BCEEDB", fontsize=7.0, fontweight="bold", ha="center", va="center", transform=ax.transAxes, zorder=7)
    for x, color in [(0.640, "source_blue"), (0.695, "assay_teal"), (0.750, "model_purple")]:
        rect(ax, colors, x - 0.020, 0.445, 0.040, 0.140, face="#0E3144", edge=color, alpha=0.94, lw=0.95, z=6)
    draw_task_surface(ax, colors, 0.812, 0.400, dark=True, z=7)
    draw_mission_gate(ax, colors, 0.705, dark=True)
    arrow(ax, colors, (0.250, 0.515), (0.300, 0.515), color="source_blue", alpha=0.62, lw=1.0)
    arrow(ax, colors, (0.560, 0.515), (0.635, 0.515), color="source_blue", alpha=0.62, lw=1.0)
    if overlay:
        draw_overlay(ax, tokens, STYLE_FRAMES["mission_system"])


def draw_hybrid(ax: plt.Axes, tokens: dict[str, Any], *, overlay: bool) -> None:
    colors = tokens["colors"]
    draw_texture(ax, colors, seed=20260632, warm=True)
    draw_grid(ax, colors, alpha=0.16)
    rect(ax, colors, 0.055, 0.330, 0.890, 0.350, face="paper", edge=None, alpha=0.48, z=2)
    rect(ax, colors, 0.090, 0.370, 0.220, 0.245, face="canvas_cool", edge=None, alpha=0.32, z=3)
    rect(ax, colors, 0.705, 0.370, 0.190, 0.245, face="#F7F4FD", edge=None, alpha=0.34, z=3)
    draw_source_stack(ax, colors, 0.125, 0.420, z=5)
    ax.plot([0.285, 0.750], [0.515, 0.515], color=colors["rule"], alpha=0.55, linewidth=1.05, transform=ax.transAxes, zorder=4)
    for x, label, color in [(0.350, "mission", "bio_green"), (0.470, "sample", "source_blue"), (0.590, "label", "label_amber")]:
        circle(ax, colors, x, 0.515, 0.043, face="paper", edge=color, alpha=0.98, lw=1.15, z=6)
        ax.text(x, 0.515, label, color=colors[color], fontsize=7.0, fontweight="bold", ha="center", va="center", transform=ax.transAxes, zorder=7)
    draw_task_surface(ax, colors, 0.735, 0.395, z=7)
    draw_mission_gate(ax, colors, 0.855)
    arrow(ax, colors, (0.270, 0.515), (0.312, 0.515), color="muted", alpha=0.48, lw=1.0)
    arrow(ax, colors, (0.630, 0.515), (0.722, 0.515), color="muted", alpha=0.48, lw=1.0)
    if overlay:
        draw_overlay(ax, tokens, STYLE_FRAMES["hybrid"])


STYLE_FRAMES: dict[str, dict[str, Any]] = {
    "editorial": {
        "style_frame_id": "editorial",
        "title": "Scientific editorial",
        "dark": False,
        "renderer": draw_editorial,
        "used_layers": ["Z1", "Z2", "Z3", "Z4", "Z5", "Z6", "Z7"],
        "depth_tokens": ["flat", "proof_object", "focus_boundary"],
        "overlay": {
            "text": [
                {"id": "headline", "role": "headline", "content": "Public studies become auditable benchmark tasks", "x": 0.065, "y": 0.135, "font_pt": 26, "color": "ink"},
                {"id": "source", "role": "callout", "content": "source records", "x": 0.108, "y": 0.350, "font_pt": 9.2, "color": "source_blue"},
                {"id": "task", "role": "callout", "content": "task contract", "x": 0.745, "y": 0.350, "font_pt": 9.2, "color": "model_purple"},
            ],
            "status_labels": [
                {"id": "scope", "role": "status_strong", "content": "editorial: manuscript-portable, lowest atmosphere", "x": 0.065, "y": 0.860, "font_pt": 7.8, "color": "muted"},
                {"id": "caveat", "role": "status", "content": "Best for paper figures; may feel underpowered for opening deck slides.", "x": 0.065, "y": 0.895, "font_pt": 7.4, "color": "muted"},
            ],
        },
    },
    "mission_system": {
        "style_frame_id": "mission_system",
        "title": "Mission system",
        "dark": True,
        "dark_text": {
            "ink": "#EAF4F8",
            "muted": "#A9C3D2",
            "source_blue": "#8AD4E8",
            "bio_green": "#BCEEDB",
            "model_purple": "#D7C6FF",
            "test_red": "#F9CACA"
        },
        "renderer": draw_mission_system,
        "used_layers": ["Z0", "Z1", "Z2", "Z3", "Z4", "Z5", "Z6", "Z7"],
        "depth_tokens": ["recessed_band", "proof_object", "dominant_evidence", "focus_boundary"],
        "overlay": {
            "text": [
                {"id": "headline", "role": "headline", "content": "Mission-held-out evaluation becomes the operating system", "x": 0.065, "y": 0.135, "font_pt": 25.5, "color": "ink"},
                {"id": "source", "role": "callout", "content": "public source", "x": 0.108, "y": 0.350, "font_pt": 9.2, "color": "source_blue"},
                {"id": "task", "role": "callout", "content": "benchmark surface", "x": 0.812, "y": 0.350, "font_pt": 9.2, "color": "model_purple"},
            ],
            "status_labels": [
                {"id": "scope", "role": "status_strong", "content": "mission-system: strongest identity, highest atmosphere", "x": 0.065, "y": 0.860, "font_pt": 7.8, "color": "muted"},
                {"id": "caveat", "role": "status", "content": "Use carefully: too much darkness can become sci-fi rather than scientific.", "x": 0.065, "y": 0.895, "font_pt": 7.4, "color": "muted"},
            ],
        },
    },
    "hybrid": {
        "style_frame_id": "hybrid",
        "title": "Mission-grade editorial hybrid",
        "dark": False,
        "renderer": draw_hybrid,
        "used_layers": ["Z0", "Z1", "Z2", "Z3", "Z4", "Z5", "Z6", "Z7"],
        "depth_tokens": ["recessed_band", "proof_object", "dominant_evidence", "focus_boundary"],
        "overlay": {
            "text": [
                {"id": "headline", "role": "headline", "content": "Space omics becomes a mission-grade benchmark system", "x": 0.065, "y": 0.135, "font_pt": 26, "color": "ink"},
                {"id": "source", "role": "callout", "content": "source evidence", "x": 0.125, "y": 0.350, "font_pt": 9.2, "color": "source_blue"},
                {"id": "task", "role": "callout", "content": "auditable task", "x": 0.735, "y": 0.350, "font_pt": 9.2, "color": "model_purple"},
            ],
            "status_labels": [
                {"id": "scope", "role": "status_strong", "content": "hybrid: recommended deck direction", "x": 0.065, "y": 0.860, "font_pt": 7.8, "color": "muted"},
                {"id": "caveat", "role": "status", "content": "Keeps scientific restraint while adding mission-grade depth and hierarchy.", "x": 0.065, "y": 0.895, "font_pt": 7.4, "color": "muted"},
            ],
        },
    },
}


def frame_contract(tokens: dict[str, Any], style: dict[str, Any]) -> dict[str, Any]:
    frame_id = style["style_frame_id"]
    frame_dir = OUT_ROOT / frame_id
    outputs = {
        "scene_plate": rel(frame_dir / "scene_plate.png"),
        "rendered_preview_png": rel(frame_dir / "rendered_preview.png"),
        "rendered_preview_pdf": rel(frame_dir / "rendered_preview.pdf"),
        "overlay_spec": rel(frame_dir / "overlay_spec.json"),
        "manifest": rel(frame_dir / "manifest.json"),
        "qa": rel(frame_dir / "qa.json"),
    }
    text_items = style["overlay"]["text"] + style["overlay"]["status_labels"]
    overlay_spec = {
        "slide_id": frame_id,
        "stage": "post_render",
        "canvas": tokens["layout"]["canvas"],
        "coordinate_system": "normalized_0_1",
        "text": style["overlay"]["text"],
        "status_labels": style["overlay"]["status_labels"],
        "focus_marks": [{"id": "mission_grade_path", "role": "motion_focus", "shape": "source_to_task_to_boundary", "z": "Z7"}],
        "visible_word_count": count_words(text_items),
        "visible_word_budget": 55,
        "forbidden_visible_terms": ["NASA logo", "patch", "dashboard", "card box", "payload", "RRRM", "alpha", "/Users/"],
    }
    manifest = {
        "slide_id": frame_id,
        "created": CREATED,
        "stage": "post_render",
        "visual_identity_token": rel(TOKEN_PATH),
        "style_direction": style["title"],
        "slide_type": "thesis",
        "content_brief": "docs/VISUAL_IDENTITY_RESEARCH_AND_DEPTH_STRATEGY_2026_06_01.md",
        "technical_preflight": "docs/VISUAL_TECHNICAL_PRODUCTION_PROTOCOL_2026_06_01.md",
        "evidence_sources": [
            {"path": "docs/VISUAL_IDENTITY_RESEARCH_AND_DEPTH_STRATEGY_2026_06_01.md", "role": "identity and z-depth strategy"},
            {"path": "config/visual_identity/spacebio_bench_visual_identity_v0_1.json", "role": "visual identity tokens"},
            {"path": "docs/VISUAL_METHODS_BRIDGE_AND_CONSULTING_BRIEF_2026_06_01.md", "role": "bridge narrative spine"},
        ],
        "claim_boundary": "style-frame only; not a final science slide",
        "used_z_layers": style["used_layers"],
        "used_depth_tokens": style["depth_tokens"],
        "anti_pattern_checks": {
            "no_card_box_layout": True,
            "no_nasa_logo_imitation": True,
            "no_decorative_space_wallpaper": True,
            "no_dense_table_inside_figure": True,
            "red_only_for_trust_boundary": True,
        },
        "generator": "scripts/build_visual_identity_style_frames.py",
        "outputs": outputs,
    }
    qa = {
        "slide_id": frame_id,
        "stage": "post_render",
        "created": CREATED,
        "depth_audit": {
            "used_z_layers": style["used_layers"],
            "used_depth_tokens": style["depth_tokens"],
            "depth_semantics": "depth separates canvas, measurement, evidence, proof objects, interpretation, trust boundary, and caveat",
            "dominant_depth_event": "one central source-to-task evidence scene",
            "card_box_risk": "low" if frame_id != "mission_system" else "medium: dark surface can read as interface if overused",
        },
        "style_scores": style_score(frame_id),
        "manual_review_pending": [
            "full_size_render_inspection",
            "thumbnail_contact_sheet_inspection",
            "direction_selection",
        ],
    }
    return {"overlay_spec": overlay_spec, "manifest": manifest, "qa": qa}


def style_score(frame_id: str) -> dict[str, int]:
    if frame_id == "editorial":
        return {
            "scientific_trust": 5,
            "premium_impression": 3,
            "thumbnail_recognizability": 3,
            "manuscript_portability": 5,
            "technical_reproducibility": 5,
        }
    if frame_id == "mission_system":
        return {
            "scientific_trust": 3,
            "premium_impression": 5,
            "thumbnail_recognizability": 5,
            "manuscript_portability": 2,
            "technical_reproducibility": 4,
        }
    return {
        "scientific_trust": 5,
        "premium_impression": 5,
        "thumbnail_recognizability": 4,
        "manuscript_portability": 4,
        "technical_reproducibility": 5,
    }


def render_frame(tokens: dict[str, Any], frame_id: str) -> dict[str, str]:
    style = STYLE_FRAMES[frame_id]
    frame_dir = OUT_ROOT / frame_id
    frame_dir.mkdir(parents=True, exist_ok=True)
    contract = frame_contract(tokens, style)
    for name, data in contract.items():
        write_json(frame_dir / f"{name}.json", data)

    fig, ax = make_figure(tokens)
    style["renderer"](ax, tokens, overlay=False)
    scene_plate = frame_dir / "scene_plate.png"
    fig.savefig(scene_plate, dpi=tokens["layout"]["canvas"]["dpi"], facecolor=fig.get_facecolor())
    plt.close(fig)

    fig, ax = make_figure(tokens)
    style["renderer"](ax, tokens, overlay=True)
    rendered_preview = frame_dir / "rendered_preview.png"
    rendered_pdf = frame_dir / "rendered_preview.pdf"
    fig.savefig(rendered_preview, dpi=tokens["layout"]["canvas"]["dpi"], facecolor=fig.get_facecolor())
    fig.savefig(rendered_pdf, facecolor=fig.get_facecolor())
    plt.close(fig)
    return {
        "scene_plate": rel(scene_plate),
        "rendered_preview_png": rel(rendered_preview),
        "rendered_preview_pdf": rel(rendered_pdf),
    }


def render_contact_sheet(tokens: dict[str, Any], frame_ids: list[str]) -> dict[str, str]:
    OUT_ROOT.mkdir(parents=True, exist_ok=True)
    output = OUT_ROOT / "visual_identity_style_frame_contact_sheet.png"
    fig = plt.figure(figsize=(19.0, 5.8), dpi=220)
    grid = fig.add_gridspec(1, len(frame_ids), left=0.018, right=0.988, top=0.790, bottom=0.060, wspace=0.026)
    fig.suptitle("Visual identity style-frame candidates", x=0.018, y=0.946, ha="left", fontsize=12.8, fontweight="bold")
    fig.text(
        0.018,
        0.882,
        "Same message, three directions: editorial portability, mission-system identity, and the recommended hybrid.",
        fontsize=7.8,
        color=tokens["colors"]["muted"],
        ha="left",
    )
    for idx, frame_id in enumerate(frame_ids):
        ax = fig.add_subplot(grid[0, idx])
        source = OUT_ROOT / frame_id / "rendered_preview.png"
        ax.imshow(mpimg.imread(source))
        style = STYLE_FRAMES[frame_id]
        scores = style_score(frame_id)
        score_text = f"trust {scores['scientific_trust']} | premium {scores['premium_impression']} | thumb {scores['thumbnail_recognizability']} | paper {scores['manuscript_portability']}"
        ax.set_title(f"{idx + 1}. {style['title']} | {score_text}", loc="left", fontsize=6.7, pad=3)
        ax.axis("off")
    fig.savefig(output, dpi=220, facecolor="white")
    plt.close(fig)
    manifest = {
        "contact_sheet": rel(output),
        "visual_identity_token": rel(TOKEN_PATH),
        "frames": [
            {
                "style_frame_id": frame_id,
                "title": STYLE_FRAMES[frame_id]["title"],
                "source": rel(OUT_ROOT / frame_id / "rendered_preview.png"),
                "scores": style_score(frame_id),
            }
            for frame_id in frame_ids
        ],
        "recommended_direction": "hybrid",
        "recommendation_reason": "best balance of scientific trust, premium impression, thumbnail recognizability, manuscript portability, and reproducibility",
    }
    manifest_path = OUT_ROOT / "visual_identity_style_frame_contact_sheet_manifest.json"
    write_json(manifest_path, manifest)
    return {"contact_sheet": rel(output), "manifest": rel(manifest_path)}


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--style", default="all", choices=["all", *STYLE_FRAMES.keys()])
    args = parser.parse_args()
    tokens = load_tokens()
    OUT_ROOT.mkdir(parents=True, exist_ok=True)
    frame_ids = list(STYLE_FRAMES) if args.style == "all" else [args.style]
    rendered = {frame_id: render_frame(tokens, frame_id) for frame_id in frame_ids}
    contact_sheet = render_contact_sheet(tokens, frame_ids) if len(frame_ids) > 1 else None
    print(json.dumps({"rendered": rendered, "contact_sheet": contact_sheet}, indent=2))


if __name__ == "__main__":
    main()
