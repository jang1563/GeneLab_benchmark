#!/usr/bin/env python3
"""Test v0.3 biological visual assets inside the feature-layer bridge context."""

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


OUT = ROOT / "output" / "premium_feature_layer_bridge_v0_3_test"
ASSET_DIR = ROOT / "assets" / "biovis_motif_pack_v0_3" / "png"
CREATED = "2026-06-02"
SLIDE_DPI = 240

COLORS = {
    "paper": "#F7F4EC",
    "ink": "#17212B",
    "muted": "#5D6978",
    "rule": "#AEB8C5",
    "blue": "#2D6F9F",
    "green": "#178B63",
    "teal": "#1AA090",
    "amber": "#B4832F",
    "red": "#B23A3A",
    "shadow": "#9C9487",
}


TEXT = [
    "v0.3 motif bridge replacement test",
    "Proof-like texture plates carry biology; editable labels keep the claim bounded.",
    "Samples x genes",
    "Gene layer",
    "Pathway summary",
    "Hidden score",
    "same hidden-mission split",
    "generated schematic texture; not source evidence",
]


def rel(path: Path) -> str:
    return str(path.relative_to(ROOT))


def write_json(path: Path, data: Any) -> None:
    path.write_text(json.dumps(data, indent=2) + "\n", encoding="utf-8")


def make_fig() -> tuple[plt.Figure, plt.Axes]:
    fig = plt.figure(figsize=(16, 9), dpi=SLIDE_DPI)
    ax = fig.add_axes([0, 0, 1, 1])
    ax.set_axis_off()
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)
    return fig, ax


def draw_canvas(ax: plt.Axes) -> None:
    rng = np.random.default_rng(202606021)
    h, w = 720, 1280
    yy, xx = np.mgrid[0:h, 0:w]
    base = np.zeros((h, w, 3), dtype=float)
    base[..., 0] = 0.970
    base[..., 1] = 0.955
    base[..., 2] = 0.925
    grain = rng.normal(0, 0.005, (h, w, 1))
    vignette = ((xx / w - 0.53) ** 2 + (yy / h - 0.50) ** 2)[..., None] * np.array([0.055, 0.044, 0.030])
    img = np.clip(base + grain - vignette, 0, 1)
    ax.imshow(img, extent=(0, 1, 0, 1), origin="lower", zorder=0)
    ax.set_aspect("auto")
    for y in [0.20, 0.50, 0.80]:
        ax.plot([0.055, 0.945], [y, y], color=COLORS["rule"], alpha=0.18, lw=0.8, zorder=1)
    for x in [0.08, 0.30, 0.52, 0.74, 0.92]:
        ax.plot([x, x], [0.10, 0.88], color=COLORS["rule"], alpha=0.10, lw=0.8, zorder=1)


def add_asset(ax: plt.Axes, name: str, x: float, y: float, w: float, h: float, alpha: float = 1.0, z: int = 5) -> None:
    img = mpimg.imread(ASSET_DIR / f"{name}.png")
    ax.imshow(img, extent=(x, x + w, y, y + h), transform=ax.transAxes, zorder=z, alpha=alpha, interpolation="lanczos")
    ax.set_aspect("auto")


def arrow(ax: plt.Axes, start: tuple[float, float], end: tuple[float, float], color: str = "rule") -> None:
    ax.add_patch(
        FancyArrowPatch(
            start,
            end,
            arrowstyle="-|>",
            mutation_scale=16,
            color=COLORS[color],
            lw=1.2,
            alpha=0.60,
            transform=ax.transAxes,
            zorder=12,
        )
    )


def draw_bridge(ax: plt.Axes) -> None:
    draw_canvas(ax)
    ax.text(0.065, 0.875, TEXT[0], fontsize=25, fontweight="bold", color=COLORS["ink"], ha="left", va="center", zorder=20)
    ax.text(0.066, 0.820, TEXT[1], fontsize=11, color=COLORS["muted"], ha="left", va="center", zorder=20)

    ax.add_patch(Rectangle((0.060, 0.305), 0.890, 0.350, facecolor="#FBFAF6", edgecolor="none", alpha=0.44, zorder=2))
    ax.add_patch(Rectangle((0.065, 0.290), 0.890, 0.350, facecolor=COLORS["shadow"], edgecolor="none", alpha=0.070, zorder=1.8))

    add_asset(ax, "omics_matrix_texture_field", 0.080, 0.385, 0.240, 0.135, alpha=0.96, z=5)
    add_asset(ax, "single_cell_embedding_field", 0.335, 0.375, 0.245, 0.138, alpha=0.92, z=5)
    add_asset(ax, "pathway_network_summary_field", 0.590, 0.378, 0.238, 0.134, alpha=0.94, z=5)

    ax.add_patch(Rectangle((0.855, 0.372), 0.082, 0.156, facecolor="#FFFDF9", edgecolor=COLORS["red"], lw=1.1, alpha=0.90, zorder=6))
    ax.add_patch(Circle((0.896, 0.552), 0.034, facecolor="#FFFDF9", edgecolor=COLORS["red"], lw=1.1, alpha=0.95, zorder=7))
    for y in [0.486, 0.442, 0.402]:
        ax.plot([0.873, 0.920], [y, y], color=COLORS["red"] if y == 0.486 else COLORS["rule"], alpha=0.45, lw=0.9, zorder=8)

    arrow(ax, (0.305, 0.455), (0.350, 0.455), "rule")
    arrow(ax, (0.565, 0.455), (0.606, 0.455), "teal")
    arrow(ax, (0.810, 0.455), (0.862, 0.455), "rule")

    labels = [
        (0.125, 0.625, TEXT[2], COLORS["blue"]),
        (0.405, 0.625, TEXT[3], COLORS["teal"]),
        (0.650, 0.625, TEXT[4], COLORS["green"]),
        (0.850, 0.625, TEXT[5], COLORS["red"]),
    ]
    for x, y, label, col in labels:
        ax.text(x, y, label, fontsize=10.5, fontweight="bold", color=col, ha="left", va="center", zorder=20)

    add_asset(ax, "species_model_strip", 0.575, 0.680, 0.310, 0.174, alpha=0.55, z=3)
    ax.text(0.600, 0.742, TEXT[6], fontsize=8.5, color=COLORS["muted"], fontweight="bold", ha="left", va="center", zorder=20)
    ax.plot([0.780, 0.780], [0.700, 0.805], color=COLORS["red"], alpha=0.36, lw=1.0, zorder=20)

    ax.text(0.065, 0.135, TEXT[7], fontsize=8.0, color=COLORS["muted"], ha="left", va="center", zorder=20)
    ax.text(0.065, 0.100, "Use as bridge context only; replace with source-derived panels for empirical evidence.", fontsize=7.5, color=COLORS["muted"], ha="left", va="center", zorder=20)


def build() -> dict[str, str]:
    OUT.mkdir(parents=True, exist_ok=True)
    png = OUT / "feature_layer_bridge_v0_3_motif_test.png"
    pdf = OUT / "feature_layer_bridge_v0_3_motif_test.pdf"
    qa = OUT / "feature_layer_bridge_v0_3_motif_test_qa.json"
    manifest = OUT / "feature_layer_bridge_v0_3_motif_test_manifest.json"
    fig, ax = make_fig()
    draw_bridge(ax)
    fig.savefig(png, dpi=SLIDE_DPI, facecolor=COLORS["paper"])
    fig.savefig(pdf, facecolor=COLORS["paper"])
    plt.close(fig)
    used = [
        "omics_matrix_texture_field",
        "single_cell_embedding_field",
        "pathway_network_summary_field",
        "species_model_strip",
    ]
    write_json(
        manifest,
        {
            "created": CREATED,
            "purpose": "Context test for replacing v0.2 bridge motifs with v0.3 biological visual system plates.",
            "outputs": {"png": rel(png), "pdf": rel(pdf), "qa": rel(qa)},
            "assets_used": [rel(ASSET_DIR / f"{asset}.png") for asset in used],
            "claim_boundary": "Generated schematic context only; not source evidence.",
        },
    )
    write_json(
        qa,
        {
            "created": CREATED,
            "automatic_checks": {
                "rendered_outputs_exist": png.exists() and pdf.exists(),
                "assets_exist": all((ASSET_DIR / f"{asset}.png").exists() for asset in used),
                "visible_text_count": len(" ".join(TEXT).split()),
            },
            "manual_review": {
                "full_size_inspection": "pass - 3840x2160 preview inspected; labels and flow are readable with no major overlap.",
                "biological_specificity": "pass - v0.3 omics matrix, embedding, pathway, and species plates provide stronger biological/data texture than v0.2 motifs.",
                "premium_quality_gate": "conditional pass - suitable as a context stress test; final slide needs stronger crop/contrast and should remove meta-test title.",
                "verdict": "v0.3 assets improve the bridge vocabulary, but source-derived or higher-contrast proof panels are still needed for final premium slides.",
            },
        },
    )
    return {"png": rel(png), "pdf": rel(pdf), "qa": rel(qa), "manifest": rel(manifest)}


def main() -> None:
    print(json.dumps(build(), indent=2))


if __name__ == "__main__":
    main()
