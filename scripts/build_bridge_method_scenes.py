#!/usr/bin/env python3
"""Render bridge-method slide prototypes from scene contracts."""

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
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.patches import Circle, FancyArrowPatch, Rectangle


SCENE_ROOT = ROOT / "output" / "premium_bridge_scenes"
SLIDE_DPI = 240

COLORS = {
    "bg": "#F7F4EF",
    "paper": "#FCFCFB",
    "ink": "#1F2933",
    "muted": "#687385",
    "rule": "#B9C1CC",
    "blue": "#2F6C9F",
    "green": "#1F8F5F",
    "amber": "#D99A22",
    "red": "#C81E1E",
    "teal": "#1F9D8A",
    "purple": "#7C3AED",
    "shadow": "#B7B0A6",
}


def load_json(path: Path) -> Any:
    with path.open(encoding="utf-8") as fh:
        return json.load(fh)


def write_json(path: Path, data: Any) -> None:
    path.write_text(json.dumps(data, indent=2) + "\n", encoding="utf-8")


def rel(path: Path) -> str:
    return str(path.relative_to(ROOT))


def resolve_repo_path(value: str) -> Path:
    path = Path(value)
    if path.is_absolute():
        return path
    return ROOT / path


def y_from_slide(value: float) -> float:
    """Convert slide-style top-origin normalized y to Matplotlib axes y."""
    return 1.0 - value


def y_range_from_slide(y0: float, y1: float) -> tuple[float, float]:
    converted = [y_from_slide(y0), y_from_slide(y1)]
    return min(converted), max(converted)


def draw_canvas(ax: plt.Axes) -> None:
    ax.add_patch(Rectangle((0, 0), 1, 1, color=COLORS["bg"], transform=ax.transAxes, zorder=0))
    h, w = 540, 960
    yy, xx = np.mgrid[0:h, 0:w]
    base = np.zeros((h, w, 3), dtype=float)
    base[..., 0] = 0.968
    base[..., 1] = 0.952
    base[..., 2] = 0.923
    grain = np.random.default_rng(20260610).normal(0, 0.0048, size=(h, w, 1))
    diagonal = ((xx + yy) / (w + h))[..., None] * np.array([0.010, 0.008, 0.004])
    texture = np.clip(base + grain - diagonal, 0, 1)
    ax.imshow(texture, extent=(0, 1, 0, 1), origin="lower", zorder=0)
    ax.set_aspect("auto")

    for y in [0.180, 0.355, 0.530, 0.705, 0.880]:
        ax.plot([0.055, 0.945], [y, y], color=COLORS["rule"], alpha=0.20, linewidth=0.75, zorder=1)
    for x in [0.055, 0.300, 0.640, 0.945]:
        ax.plot([x, x], [0.075, 0.900], color=COLORS["rule"], alpha=0.14, linewidth=0.75, zorder=1)


def draw_arrow(
    ax: plt.Axes,
    start: tuple[float, float],
    end: tuple[float, float],
    color: str = "muted",
    *,
    alpha: float = 0.70,
    width: float = 1.0,
    zorder: int = 4,
) -> None:
    ax.add_patch(
        FancyArrowPatch(
            start,
            end,
            arrowstyle="-|>",
            mutation_scale=10,
            linewidth=width,
            color=COLORS[color],
            alpha=alpha,
            transform=ax.transAxes,
            zorder=zorder,
        )
    )


def draw_mission_node(
    ax: plt.Axes,
    x: float,
    y: float,
    label: str,
    color: str,
    *,
    radius: float = 0.025,
    alpha: float = 1.0,
    zorder: int = 5,
) -> None:
    ax.add_patch(
        Circle(
            (x, y),
            radius,
            transform=ax.transAxes,
            facecolor=COLORS["paper"],
            edgecolor=COLORS[color],
            linewidth=1.35,
            alpha=alpha,
            zorder=zorder,
        )
    )
    ax.text(x, y, label, color=COLORS[color], fontsize=7.3, fontweight="bold", ha="center", va="center", zorder=zorder + 1)


def draw_operation_tick(ax: plt.Axes, x: float, y: float, color: str, *, zorder: int = 5) -> None:
    ax.plot([x, x], [y - 0.045, y + 0.045], color=COLORS[color], linewidth=1.35, alpha=0.86, zorder=zorder, transform=ax.transAxes)


def draw_scene_base(slide_id: str, ax: plt.Axes, focus_marks: list[dict[str, Any]]) -> None:
    draw_canvas(ax)
    if slide_id == "b3_mission_held_out":
        draw_b3_base(ax, focus_marks)
    elif slide_id == "b4_train_only_guard":
        draw_b4_base(ax, focus_marks)
    else:
        raise ValueError(f"Unsupported slide_id: {slide_id}")


def draw_b3_base(ax: plt.Axes, focus_marks: list[dict[str, Any]]) -> None:
    ax.add_patch(
        Rectangle((0.092, 0.330), 0.472, 0.300, transform=ax.transAxes, facecolor="#FCFCFB", edgecolor="none", alpha=0.50, zorder=2)
    )
    ax.add_patch(
        Rectangle((0.690, 0.330), 0.170, 0.300, transform=ax.transAxes, facecolor="#FCFCFB", edgecolor="none", alpha=0.56, zorder=2)
    )
    y = 0.465
    train_nodes = [(0.170, "M1"), (0.270, "M2"), (0.370, "M3")]
    for x, label in train_nodes:
        draw_mission_node(ax, x, y, label, "green")
    ax.plot([0.145, 0.420], [y, y], color=COLORS["green"], alpha=0.24, linewidth=4.0, zorder=3, transform=ax.transAxes)
    draw_arrow(ax, (0.420, y), (0.565, y), "muted", alpha=0.72, width=1.0)
    ax.add_patch(
        Circle((0.565, y), 0.037, transform=ax.transAxes, facecolor="#F6F8FA", edgecolor=COLORS["rule"], linewidth=0.9, zorder=4)
    )
    for y0 in [0.452, 0.465, 0.478]:
        ax.plot([0.548, 0.582], [y0, y0], color=COLORS["rule"], alpha=0.85, linewidth=1.0, zorder=5, transform=ax.transAxes)
    boundary_x = focus_marks[0].get("x", 0.640) if focus_marks else 0.640
    if focus_marks:
        boundary_y0, boundary_y1 = y_range_from_slide(float(focus_marks[0].get("y0", 0.320)), float(focus_marks[0].get("y1", 0.725)))
    else:
        boundary_y0, boundary_y1 = 0.305, 0.725
    ax.plot([boundary_x, boundary_x], [boundary_y0, boundary_y1], color=COLORS["red"], linewidth=1.25, alpha=0.70, zorder=4, transform=ax.transAxes)
    draw_mission_node(ax, 0.745, y, "M4", "red", radius=0.028)
    draw_arrow(ax, (0.745, 0.505), (0.785, 0.635), "red", alpha=0.50, width=0.95)
    ax.add_patch(
        Circle((0.805, 0.665), 0.030, transform=ax.transAxes, facecolor="#FCFCFB", edgecolor=COLORS["red"], linewidth=1.0, alpha=0.90, zorder=5)
    )
    ax.plot([0.792, 0.818], [0.665, 0.665], color=COLORS["red"], alpha=0.75, linewidth=1.0, zorder=6, transform=ax.transAxes)
    ax.plot([0.805, 0.805], [0.652, 0.678], color=COLORS["red"], alpha=0.75, linewidth=1.0, zorder=6, transform=ax.transAxes)
    for x in [0.170, 0.270, 0.370, 0.745]:
        ax.plot([x, x], [0.250, 0.278], color=COLORS["rule"], alpha=0.45, linewidth=0.9, zorder=2, transform=ax.transAxes)
    ax.plot([0.170, 0.745], [0.264, 0.264], color=COLORS["rule"], alpha=0.35, linewidth=0.9, zorder=2, transform=ax.transAxes)


def draw_b4_base(ax: plt.Axes, focus_marks: list[dict[str, Any]]) -> None:
    train_y = 0.435
    hidden_y = 0.610
    ax.plot([0.140, 0.705], [train_y, train_y], color=COLORS["rule"], linewidth=1.1, alpha=0.78, zorder=3, transform=ax.transAxes)
    ax.plot([0.140, 0.765], [hidden_y, hidden_y], color=COLORS["rule"], linewidth=1.0, alpha=0.52, zorder=3, transform=ax.transAxes)
    for x, label in [(0.160, "M1"), (0.225, "M2"), (0.290, "M3")]:
        draw_mission_node(ax, x, train_y, label, "green", radius=0.020)
    draw_mission_node(ax, 0.160, hidden_y, "M4", "red", radius=0.022)
    for x, color in [(0.355, "blue"), (0.505, "teal"), (0.650, "purple")]:
        draw_operation_tick(ax, x, train_y, color)
    guard_x = focus_marks[0].get("x", 0.730) if focus_marks else 0.730
    if focus_marks:
        guard_y0, guard_y1 = y_range_from_slide(float(focus_marks[0].get("y0", 0.335)), float(focus_marks[0].get("y1", 0.705)))
    else:
        guard_y0, guard_y1 = 0.320, 0.710
    ax.plot([guard_x, guard_x], [guard_y0, guard_y1], color=COLORS["red"], linewidth=1.25, alpha=0.70, zorder=4, transform=ax.transAxes)
    draw_arrow(ax, (0.705, train_y), (0.785, 0.545), "muted", alpha=0.55, width=0.9)
    draw_arrow(ax, (0.205, hidden_y), (0.785, hidden_y), "red", alpha=0.36, width=0.95)
    draw_operation_tick(ax, 0.805, hidden_y, "red")
    ax.add_patch(
        Rectangle((0.315, 0.372), 0.372, 0.125, transform=ax.transAxes, facecolor="#FCFCFB", edgecolor="none", alpha=0.46, zorder=2)
    )
    ax.add_patch(
        Rectangle((0.112, 0.560), 0.148, 0.098, transform=ax.transAxes, facecolor="#FCFCFB", edgecolor="none", alpha=0.40, zorder=2)
    )


def render_overlay(ax: plt.Axes, overlay: dict[str, Any]) -> None:
    for item in list(overlay.get("text", [])) + list(overlay.get("status_labels", [])):
        color = COLORS.get(str(item.get("color", "ink")), COLORS["ink"])
        role = str(item.get("role", ""))
        weight = "bold" if role in {"decision_headline", "primary_callout", "process_verb", "claim_boundary"} else "normal"
        va = "top" if role == "decision_headline" else "center"
        ax.text(
            float(item["x"]),
            y_from_slide(float(item["y"])),
            str(item["content"]),
            color=color,
            fontsize=float(item.get("font_pt", 8.0)),
            fontweight=weight,
            ha="left",
            va=va,
            linespacing=1.05,
            transform=ax.transAxes,
            zorder=8,
        )


def make_figure() -> tuple[plt.Figure, plt.Axes]:
    fig = plt.figure(figsize=(16, 9), dpi=SLIDE_DPI)
    ax = fig.add_axes([0, 0, 1, 1])
    ax.set_axis_off()
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)
    return fig, ax


def save_scene(slide_dir: Path, manifest: dict[str, Any], overlay: dict[str, Any]) -> dict[str, Path]:
    outputs = manifest["outputs"]
    scene_plate = resolve_repo_path(outputs["scene_plate"])
    preview_png = resolve_repo_path(outputs["rendered_preview_png"])
    preview_pdf = resolve_repo_path(outputs["rendered_preview_pdf"])
    slide_id = str(manifest["slide_id"])
    focus_marks = list(overlay.get("focus_marks", []))

    fig, ax = make_figure()
    draw_scene_base(slide_id, ax, focus_marks)
    fig.savefig(scene_plate, dpi=SLIDE_DPI, facecolor=COLORS["bg"])
    plt.close(fig)

    fig, ax = make_figure()
    draw_scene_base(slide_id, ax, focus_marks)
    render_overlay(ax, overlay)
    fig.savefig(preview_png, dpi=SLIDE_DPI, facecolor=COLORS["bg"])
    fig.savefig(preview_pdf, facecolor=COLORS["bg"])
    plt.close(fig)

    return {"scene_plate": scene_plate, "rendered_preview_png": preview_png, "rendered_preview_pdf": preview_pdf}


def update_contract_after_render(slide_dir: Path, manifest: dict[str, Any], overlay: dict[str, Any], outputs: dict[str, Path]) -> None:
    manifest["stage"] = "post_render"
    manifest["renderer"] = "scripts/build_bridge_method_scenes.py"
    manifest["qa"]["stage"] = "post_render"
    manifest["qa"]["render_outputs_exist"] = True
    write_json(slide_dir / "manifest.json", manifest)

    overlay["stage"] = "post_render"
    write_json(slide_dir / "overlay_spec.json", overlay)

    qa = load_json(slide_dir / "qa.json")
    qa["stage"] = "post_render"
    qa["post_render_gate"] = {
        "rendered_outputs": {key: rel(path) for key, path in outputs.items()},
        "manual_full_size_inspection": "pending",
        "manual_thumbnail_inspection": "pending",
        "text_overlap_check": "pending",
        "caveat_visibility_check": "pending",
        "visual_verdict": "draft render; awaiting manual inspection",
    }
    write_json(slide_dir / "qa.json", qa)


def render_slide(slide_id: str) -> dict[str, str]:
    slide_dir = SCENE_ROOT / slide_id
    manifest = load_json(slide_dir / "manifest.json")
    overlay = load_json(slide_dir / "overlay_spec.json")
    outputs = save_scene(slide_dir, manifest, overlay)
    update_contract_after_render(slide_dir, manifest, overlay, outputs)
    return {key: rel(path) for key, path in outputs.items()}


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--slide", default="b3_mission_held_out", help="Slide ID to render, or 'all'.")
    args = parser.parse_args()

    if args.slide == "all":
        slide_ids = ["b3_mission_held_out", "b4_train_only_guard"]
    else:
        slide_ids = [args.slide]
    rendered = {slide_id: render_slide(slide_id) for slide_id in slide_ids}
    print(json.dumps({"rendered": rendered}, indent=2))


if __name__ == "__main__":
    main()
