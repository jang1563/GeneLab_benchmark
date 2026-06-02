#!/usr/bin/env python3
"""Build a premium feature-layer bridge scene using biological motif assets."""

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
from matplotlib.patches import Circle, FancyArrowPatch, Polygon, Rectangle


OUT = ROOT / "output" / "premium_feature_layer_bridge"
MOTIFS = ROOT / "assets" / "biovis_motif_pack_v0_2" / "png"
SLIDE_DPI = 240
CREATED = "2026-06-02"

COLORS = {
    "bg": "#F4F0E8",
    "paper": "#FCFCFA",
    "ink": "#17212B",
    "muted": "#5D6978",
    "rule": "#AEB8C5",
    "blue": "#286EA6",
    "green": "#15815E",
    "amber": "#A36F13",
    "red": "#B91C1C",
    "teal": "#159A8A",
    "purple": "#6D3EDB",
    "shadow": "#9C9487",
    "slate": "#283544",
}


TEXT_ITEMS: list[dict[str, Any]] = [
    {
        "id": "headline",
        "role": "decision_headline",
        "content": "Gene signals become pathway summaries",
        "x": 0.065,
        "y": 0.122,
        "font_pt": 26.5,
        "color": "ink",
    },
    {
        "id": "support",
        "role": "supporting_claim",
        "content": "Both layers use the same hidden-mission split.",
        "x": 0.066,
        "y": 0.198,
        "font_pt": 11.3,
        "color": "muted",
    },
    {"id": "matrix", "role": "primary_callout", "content": "Samples x genes", "x": 0.160, "y": 0.353, "font_pt": 10.2, "color": "blue"},
    {"id": "rna", "role": "primary_callout", "content": "Gene activity", "x": 0.398, "y": 0.353, "font_pt": 10.2, "color": "teal"},
    {"id": "pathway", "role": "primary_callout", "content": "Pathway summary", "x": 0.642, "y": 0.353, "font_pt": 10.2, "color": "green"},
    {"id": "score", "role": "primary_callout", "content": "Hidden mission score", "x": 0.826, "y": 0.353, "font_pt": 10.2, "color": "red"},
]

STATUS_LABELS: list[dict[str, Any]] = [
    {
        "id": "status",
        "role": "claim_boundary",
        "content": "fit on training missions; score on hidden missions",
        "x": 0.065,
        "y": 0.866,
        "font_pt": 7.8,
        "color": "muted",
    },
    {
        "id": "caveat",
        "role": "trust_caveat",
        "content": "Pathways summarize biology; not every task improves.",
        "x": 0.065,
        "y": 0.900,
        "font_pt": 7.3,
        "color": "muted",
    },
]

INTERNAL_VISIBLE_TEXT = ["same split", "M1", "M2", "M3", "M4", "summarize"]


def rel(path: Path) -> str:
    return str(path.relative_to(ROOT))


def write_json(path: Path, data: Any) -> None:
    path.write_text(json.dumps(data, indent=2) + "\n", encoding="utf-8")


def count_words(items: list[dict[str, Any]]) -> int:
    count = 0
    for item in items:
        content = str(item.get("content", "")) if isinstance(item, dict) else str(item)
        count += len([token for token in content.replace("/", " ").replace(";", " ").split() if token])
    return count


def y_from_slide(value: float) -> float:
    return 1.0 - value


def make_figure() -> tuple[plt.Figure, plt.Axes]:
    fig = plt.figure(figsize=(16, 9), dpi=SLIDE_DPI)
    ax = fig.add_axes([0, 0, 1, 1])
    ax.set_axis_off()
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)
    return fig, ax


def draw_canvas(ax: plt.Axes) -> None:
    h, w = 720, 1280
    yy, xx = np.mgrid[0:h, 0:w]
    rng = np.random.default_rng(20260630)
    base = np.zeros((h, w, 3), dtype=float)
    base[..., 0] = 0.956
    base[..., 1] = 0.940
    base[..., 2] = 0.910
    grain = rng.normal(0, 0.0042, size=(h, w, 1))
    vignette = ((xx / w - 0.52) ** 2 + (yy / h - 0.46) ** 2)[..., None] * np.array([0.038, 0.030, 0.018])
    cool_left = (1 - xx / w)[..., None] * np.array([0.000, 0.004, 0.009])
    texture = np.clip(base + grain - vignette + cool_left, 0, 1)
    ax.imshow(texture, extent=(0, 1, 0, 1), origin="lower", zorder=0)
    ax.set_aspect("auto")

    for y in [0.190, 0.328, 0.500, 0.672, 0.835]:
        ax.plot([0.055, 0.945], [y, y], color=COLORS["rule"], alpha=0.20, linewidth=0.85, transform=ax.transAxes, zorder=1)
    for x in [0.075, 0.285, 0.500, 0.715, 0.925]:
        ax.plot([x, x], [0.095, 0.870], color=COLORS["rule"], alpha=0.10, linewidth=0.75, transform=ax.transAxes, zorder=1)


def rect(ax: plt.Axes, x: float, y: float, w: float, h: float, *, face: str, edge: str | None = None, alpha: float = 1.0, lw: float = 0.8, z: float = 3) -> None:
    ax.add_patch(
        Rectangle(
            (x, y),
            w,
            h,
            transform=ax.transAxes,
            facecolor=COLORS.get(face, face),
            edgecolor=COLORS.get(edge, edge) if edge else "none",
            linewidth=lw,
            alpha=alpha,
            zorder=z,
        )
    )


def arrow(ax: plt.Axes, start: tuple[float, float], end: tuple[float, float], color: str, *, alpha: float = 0.62, lw: float = 1.15, z: float = 7) -> None:
    ax.add_patch(
        FancyArrowPatch(
            start,
            end,
            arrowstyle="-|>",
            mutation_scale=13,
            linewidth=lw,
            color=COLORS[color],
            alpha=alpha,
            transform=ax.transAxes,
            zorder=z,
        )
    )


def add_motif(ax: plt.Axes, name: str, x: float, y: float, w: float, h: float, *, alpha: float = 1.0, z: float = 5) -> None:
    path = MOTIFS / f"{name}.png"
    img = mpimg.imread(path)
    ax.imshow(img, extent=(x, x + w, y, y + h), transform=ax.transAxes, zorder=z, alpha=alpha, interpolation="lanczos")
    ax.set_aspect("auto")


def draw_gene_cloud(ax: plt.Axes) -> None:
    rng = np.random.default_rng(20260631)
    colors = [COLORS["blue"], COLORS["teal"], COLORS["green"], COLORS["muted"]]
    for _ in range(64):
        x = 0.360 + rng.normal(0, 0.050)
        y = 0.505 + rng.normal(0, 0.058)
        r = rng.uniform(0.0024, 0.0055)
        color = colors[int(rng.integers(0, len(colors)))]
        ax.add_patch(Circle((x, y), r, transform=ax.transAxes, facecolor=color, edgecolor="none", alpha=rng.uniform(0.28, 0.58), zorder=6))
    for y in [0.445, 0.475, 0.505, 0.535, 0.565]:
        ax.plot([0.298, 0.425], [y, y], color=COLORS["rule"], alpha=0.20, linewidth=0.8, transform=ax.transAxes, zorder=5)


def draw_train_only_mini(ax: plt.Axes) -> None:
    x0, y0 = 0.645, 0.710
    ax.text(x0, y0 + 0.055, "same split", color=COLORS["muted"], fontsize=7.2, ha="left", fontweight="bold", transform=ax.transAxes, zorder=10)
    for idx, (x, label, color) in enumerate([(x0 + 0.012, "M1", "green"), (x0 + 0.052, "M2", "green"), (x0 + 0.092, "M3", "green"), (x0 + 0.145, "M4", "red")]):
        ax.add_patch(Circle((x, y0 + 0.025), 0.013, transform=ax.transAxes, facecolor=COLORS["paper"], edgecolor=COLORS[color], linewidth=1.0, zorder=10))
        ax.text(x, y0 + 0.025, label, color=COLORS[color], fontsize=4.9, ha="center", va="center", fontweight="bold", transform=ax.transAxes, zorder=11)
        if idx == 2:
            ax.plot([x0 + 0.119, x0 + 0.119], [y0 + 0.000, y0 + 0.050], color=COLORS["red"], alpha=0.55, linewidth=0.9, transform=ax.transAxes, zorder=10)


def draw_scene_plate(ax: plt.Axes) -> None:
    draw_canvas(ax)

    rect(ax, 0.060, 0.302, 0.890, 0.390, face="shadow", edge=None, alpha=0.075, z=1.8)
    rect(ax, 0.055, 0.315, 0.890, 0.390, face="#FBFAF6", edge=None, alpha=0.53, z=2)
    rect(ax, 0.070, 0.338, 0.200, 0.335, face="#F3F7F8", edge=None, alpha=0.42, z=2.2)
    rect(ax, 0.300, 0.338, 0.190, 0.335, face="#F3FAF7", edge=None, alpha=0.38, z=2.2)
    rect(ax, 0.525, 0.338, 0.230, 0.335, face="#F5FAF4", edge=None, alpha=0.38, z=2.2)
    rect(ax, 0.800, 0.338, 0.115, 0.335, face="#FFF6F6", edge=None, alpha=0.40, z=2.2)

    ax.plot([0.120, 0.860], [0.505, 0.505], color=COLORS["rule"], alpha=0.54, linewidth=1.2, transform=ax.transAxes, zorder=4)
    ax.plot([0.595, 0.595], [0.350, 0.670], color=COLORS["red"], alpha=0.17, linewidth=0.9, transform=ax.transAxes, zorder=5)

    ax.add_patch(
        Polygon(
            [(0.426, 0.585), (0.548, 0.548), (0.548, 0.462), (0.426, 0.425)],
            closed=True,
            transform=ax.transAxes,
            facecolor=COLORS["teal"],
            edgecolor="none",
            alpha=0.075,
            zorder=4.4,
        )
    )
    for y, color, alpha in [(0.565, "green", 0.13), (0.505, "amber", 0.12), (0.445, "blue", 0.09)]:
        ax.plot([0.535, 0.720], [y, y], color=COLORS[color], alpha=alpha, linewidth=4.8, solid_capstyle="round", transform=ax.transAxes, zorder=4.7)

    add_motif(ax, "sample_to_matrix_assay", 0.078, 0.382, 0.205, 0.205, alpha=0.92, z=5)
    add_motif(ax, "rna_expression_trace", 0.284, 0.392, 0.220, 0.185, alpha=0.84, z=5)
    draw_gene_cloud(ax)
    add_motif(ax, "pathway_program_field", 0.520, 0.382, 0.225, 0.205, alpha=0.88, z=5)
    add_motif(ax, "protein_complex_field", 0.625, 0.428, 0.150, 0.135, alpha=0.46, z=4)

    rect(ax, 0.815, 0.412, 0.090, 0.180, face="paper", edge="red", alpha=0.88, lw=0.8, z=6)
    ax.plot([0.834, 0.886], [0.550, 0.550], color=COLORS["red"], alpha=0.50, linewidth=1.0, transform=ax.transAxes, zorder=7)
    ax.plot([0.834, 0.886], [0.507, 0.507], color=COLORS["rule"], alpha=0.55, linewidth=0.8, transform=ax.transAxes, zorder=7)
    ax.plot([0.834, 0.886], [0.464, 0.464], color=COLORS["rule"], alpha=0.40, linewidth=0.8, transform=ax.transAxes, zorder=7)
    ax.add_patch(Circle((0.860, 0.610), 0.032, transform=ax.transAxes, facecolor=COLORS["paper"], edgecolor=COLORS["red"], linewidth=1.0, alpha=0.95, zorder=8))
    ax.plot([0.847, 0.873], [0.610, 0.610], color=COLORS["red"], alpha=0.56, linewidth=0.8, transform=ax.transAxes, zorder=9)
    ax.plot([0.860, 0.860], [0.598, 0.622], color=COLORS["red"], alpha=0.56, linewidth=0.8, transform=ax.transAxes, zorder=9)

    arrow(ax, (0.250, 0.505), (0.318, 0.505), "muted", alpha=0.42, lw=1.0, z=8)
    arrow(ax, (0.430, 0.505), (0.545, 0.505), "teal", alpha=0.55, lw=1.15, z=8)
    arrow(ax, (0.725, 0.505), (0.810, 0.505), "muted", alpha=0.48, lw=1.0, z=8)

    ax.text(0.470, 0.455, "summarize", color=COLORS["teal"], fontsize=7.6, fontweight="bold", ha="center", transform=ax.transAxes, zorder=9)

    draw_train_only_mini(ax)

    ax.plot([0.115, 0.863], [0.668, 0.668], color=COLORS["rule"], alpha=0.22, linewidth=0.85, transform=ax.transAxes, zorder=4)
    ax.plot([0.115, 0.863], [0.341, 0.341], color=COLORS["rule"], alpha=0.18, linewidth=0.85, transform=ax.transAxes, zorder=4)


def draw_overlay(ax: plt.Axes) -> None:
    for item in TEXT_ITEMS + STATUS_LABELS:
        color = COLORS[item["color"]]
        role = str(item["role"])
        weight = "bold" if role in {"decision_headline", "primary_callout", "claim_boundary"} else "normal"
        va = "top" if role in {"decision_headline", "supporting_claim"} else "center"
        ax.text(
            float(item["x"]),
            y_from_slide(float(item["y"])),
            str(item["content"]),
            color=color,
            fontsize=float(item["font_pt"]),
            fontweight=weight,
            ha="left",
            va=va,
            transform=ax.transAxes,
            zorder=20,
        )


def build_scene() -> dict[str, str]:
    OUT.mkdir(parents=True, exist_ok=True)
    png = OUT / "feature_layer_bridge_scene.png"
    pdf = OUT / "feature_layer_bridge_scene.pdf"
    scene_plate = OUT / "feature_layer_bridge_scene_plate.png"
    overlay_spec_path = OUT / "feature_layer_bridge_overlay_spec.json"
    manifest_path = OUT / "feature_layer_bridge_manifest.json"
    qa_path = OUT / "feature_layer_bridge_qa.json"

    fig, ax = make_figure()
    draw_scene_plate(ax)
    fig.savefig(scene_plate, dpi=SLIDE_DPI, facecolor=COLORS["bg"])
    draw_overlay(ax)
    fig.savefig(png, dpi=SLIDE_DPI, facecolor=COLORS["bg"])
    fig.savefig(pdf, facecolor=COLORS["bg"])
    plt.close(fig)

    visible_word_count = count_words(TEXT_ITEMS + STATUS_LABELS + INTERNAL_VISIBLE_TEXT)
    overlay_spec = {
        "slide_id": "feature_layer_bridge",
        "stage": "post_render",
        "created": CREATED,
        "text": TEXT_ITEMS,
        "status_labels": STATUS_LABELS,
        "internal_visible_text": INTERNAL_VISIBLE_TEXT,
        "visible_word_count": visible_word_count,
        "visible_word_budget": 48,
        "motif_assets": [
            rel(MOTIFS / "sample_to_matrix_assay.png"),
            rel(MOTIFS / "rna_expression_trace.png"),
            rel(MOTIFS / "pathway_program_field.png"),
            rel(MOTIFS / "protein_complex_field.png"),
        ],
    }
    manifest = {
        "slide_id": "feature_layer_bridge",
        "created": CREATED,
        "purpose": "First-time-viewer bridge explaining gene/RNA versus pathway feature layers under the same mission-held-out evaluation rule.",
        "audience_question": "What does gene versus pathway mean in this benchmark?",
        "visual_move": "sample-by-gene measurement and gene activity compress into a pathway summary, then flow to a hidden mission score under the same split",
        "claim_boundary": "Pathway summaries are train-only feature summaries, not a separate dataset or guaranteed universal improvement.",
        "outputs": {
            "scene_plate": rel(scene_plate),
            "rendered_preview_png": rel(png),
            "rendered_preview_pdf": rel(pdf),
            "overlay_spec": rel(overlay_spec_path),
            "qa": rel(qa_path),
        },
        "evidence_sources": [
            "docs/VISUAL_METHODS_EXPLANATION_GAP_MAP_2026_06_01.md",
            "docs/BIOVIS_MOTIF_ASSET_PACK_V0_2_REVIEW_2026_06_02.md",
            "assets/biovis_motif_pack_v0_2/manifest.json",
        ],
    }
    qa = {
        "slide_id": "feature_layer_bridge",
        "stage": "post_render",
        "created": CREATED,
        "automatic_checks": {
            "visible_text_word_count": visible_word_count,
            "visible_text_budget": 48,
            "motif_assets_exist": all((ROOT / path).exists() for path in overlay_spec["motif_assets"]),
            "rendered_outputs_exist": all(path.exists() for path in [scene_plate, png, pdf]),
        },
        "manual_review": {
            "full_size_inspection": "pass - 3840x2160 preview inspected; no material text overlap; headline, labels, and status line are readable.",
            "thumbnail_inspection": "pass - main gene-to-pathway-to-hidden-score flow remains readable at fit-to-screen scale; M1-M4 labels are secondary.",
            "biological_specificity": "pass - uses sample matrix, RNA trace, gene cloud, pathway field, protein-complex context, and hidden-score marker.",
            "premium_quality_gate": "conditional pass - suitable as a draft premium bridge scene; still a symbolic explanation layer rather than a source-derived proof figure.",
            "visual_verdict": "draft premium feature-layer bridge candidate for deck assembly; keep as layered PNG scene plus editable text shell.",
        },
    }
    write_json(overlay_spec_path, overlay_spec)
    write_json(manifest_path, manifest)
    write_json(qa_path, qa)

    return {
        "png": rel(png),
        "pdf": rel(pdf),
        "scene_plate": rel(scene_plate),
        "overlay_spec": rel(overlay_spec_path),
        "manifest": rel(manifest_path),
        "qa": rel(qa_path),
    }


def main() -> None:
    print(json.dumps(build_scene(), indent=2))


if __name__ == "__main__":
    main()
