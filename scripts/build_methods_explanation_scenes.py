#!/usr/bin/env python3
"""Build first-time-viewer methods explainer slide scenes."""

from __future__ import annotations

import json
import os
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
os.environ.setdefault("MPLCONFIGDIR", str(ROOT / "output" / ".matplotlib"))

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.patches import Arc, Circle, FancyArrowPatch, Rectangle


OUT_DIR = ROOT / "output" / "premium_methods_scenes"
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
}


def draw_canvas(ax: plt.Axes) -> None:
    ax.add_patch(Rectangle((0, 0), 1, 1, color=COLORS["bg"], transform=ax.transAxes, zorder=0))
    h, w = 540, 960
    yy, xx = np.mgrid[0:h, 0:w]
    base = np.zeros((h, w, 3), dtype=float)
    base[..., 0] = 0.968
    base[..., 1] = 0.952
    base[..., 2] = 0.923
    grain = np.random.default_rng(20260608).normal(0, 0.0055, size=(h, w, 1))
    vertical = (yy / h)[..., None] * np.array([0.010, 0.008, 0.004])
    texture = np.clip(base + grain - vertical, 0, 1)
    ax.imshow(texture, extent=(0, 1, 0, 1), origin="lower", zorder=0)
    ax.set_aspect("auto")

    for y in [0.18, 0.36, 0.56, 0.76]:
        ax.plot([0.055, 0.945], [y, y], color=COLORS["rule"], alpha=0.23, linewidth=0.8, zorder=1)
    for x in [0.055, 0.945]:
        ax.plot([x, x], [0.080, 0.900], color=COLORS["rule"], alpha=0.18, linewidth=0.8, zorder=1)


def draw_step(
    ax: plt.Axes,
    x: float,
    y: float,
    number: str,
    title: str,
    detail: str,
    color: str,
    *,
    emphasize: bool = False,
) -> None:
    radius = 0.023 if not emphasize else 0.027
    ax.add_patch(Circle((x, y), radius, transform=ax.transAxes, facecolor=COLORS["paper"], edgecolor=color, linewidth=1.7, zorder=5))
    ax.text(x, y, number, color=color, fontsize=8.6, fontweight="bold", ha="center", va="center", zorder=6)
    ax.text(x, y - 0.057, title, color=COLORS["ink"], fontsize=9.8, fontweight="bold", ha="center", va="top", zorder=6)
    ax.text(x, y - 0.092, detail, color=COLORS["muted"], fontsize=7.1, ha="center", va="top", linespacing=1.18, zorder=6)


def draw_arrow(ax: plt.Axes, start: tuple[float, float], end: tuple[float, float], color: str = "#808A98") -> None:
    ax.add_patch(
        FancyArrowPatch(
            start,
            end,
            arrowstyle="-|>",
            mutation_scale=9.5,
            linewidth=1.0,
            color=color,
            alpha=0.82,
            transform=ax.transAxes,
            zorder=4,
        )
    )


def draw_hidden_mission_inset(ax: plt.Axes) -> None:
    x0, y0 = 0.075, 0.223
    ax.text(x0, y0 + 0.130, "One mission stays hidden", color=COLORS["ink"], fontsize=10.4, fontweight="bold", ha="left", zorder=6)
    ax.text(
        x0,
        y0 + 0.101,
        "Train on other missions; score only on the hidden mission.",
        color=COLORS["muted"],
        fontsize=7.5,
        ha="left",
        zorder=6,
    )
    missions = [
        ("M1", COLORS["green"], False),
        ("M2", COLORS["green"], False),
        ("M3", COLORS["green"], False),
        ("M4", COLORS["red"], True),
    ]
    for idx, (label, color, hidden) in enumerate(missions):
        x = x0 + idx * 0.058
        ax.add_patch(Circle((x, y0 + 0.050), 0.018, transform=ax.transAxes, facecolor=COLORS["paper"], edgecolor=color, linewidth=1.2, zorder=6))
        ax.text(x, y0 + 0.050, label, color=color, fontsize=6.7, fontweight="bold", ha="center", va="center", zorder=7)
        ax.text(
            x,
            y0 + 0.010,
            "test" if hidden else "train",
            color=COLORS["muted"],
            fontsize=5.9,
            ha="center",
            va="top",
            zorder=7,
        )
    ax.plot([x0 + 0.200, x0 + 0.200], [y0 + 0.005, y0 + 0.085], color=COLORS["red"], linewidth=1.1, alpha=0.65, zorder=6)
    ax.text(x0 + 0.217, y0 + 0.066, "no test leakage", color=COLORS["red"], fontsize=7.1, fontweight="bold", ha="left", zorder=6)


def draw_train_only_guard(ax: plt.Axes) -> None:
    x0, y0 = 0.650, 0.220
    ax.text(x0, y0 + 0.128, "Train-only guard", color=COLORS["ink"], fontsize=10.4, fontweight="bold", ha="left", zorder=6)
    ax.text(x0, y0 + 0.100, "Feature selection and scaling are learned inside each fold.", color=COLORS["muted"], fontsize=7.5, ha="left", zorder=6)
    ax.plot([x0 + 0.000, x0 + 0.250], [y0 + 0.055, y0 + 0.055], color=COLORS["rule"], linewidth=1.0, alpha=0.82, zorder=5)
    for x, label, color in [
        (x0 + 0.030, "select", COLORS["blue"]),
        (x0 + 0.108, "scale", COLORS["teal"]),
        (x0 + 0.186, "fit", COLORS["purple"]),
        (x0 + 0.250, "score", COLORS["red"]),
    ]:
        ax.plot([x, x], [y0 + 0.031, y0 + 0.079], color=color, linewidth=1.2, zorder=6)
        ax.text(x, y0 + 0.006, label, color=COLORS["muted"], fontsize=6.2, ha="center", va="top", zorder=6)
    ax.plot([x0 + 0.220, x0 + 0.220], [y0 + 0.020, y0 + 0.090], color=COLORS["red"], linewidth=1.0, alpha=0.70, zorder=6)
    ax.text(x0 + 0.232, y0 + 0.075, "hidden mission", color=COLORS["red"], fontsize=6.6, ha="left", zorder=6)


def draw_gene_to_pathway(ax: plt.Axes) -> None:
    x0, y0 = 0.380, 0.214
    ax.text(x0, y0 + 0.137, "Genes become pathway views", color=COLORS["ink"], fontsize=10.4, fontweight="bold", ha="left", zorder=6)
    ax.text(x0, y0 + 0.108, "Many gene signals can be summarized into biological programs.", color=COLORS["muted"], fontsize=7.5, ha="left", zorder=6)
    rng = np.random.default_rng(20260609)
    for _ in range(28):
        x = x0 + rng.uniform(0.002, 0.095)
        y = y0 + rng.uniform(0.025, 0.083)
        color = rng.choice([COLORS["blue"], COLORS["teal"], COLORS["muted"]])
        ax.add_patch(Circle((x, y), rng.uniform(0.0022, 0.0042), transform=ax.transAxes, facecolor=color, edgecolor="none", alpha=0.45, zorder=6))
    draw_arrow(ax, (x0 + 0.115, y0 + 0.055), (x0 + 0.165, y0 + 0.055), COLORS["teal"])
    for idx, color in enumerate([COLORS["green"], COLORS["amber"], COLORS["blue"]]):
        y = y0 + 0.026 + idx * 0.025
        ax.plot([x0 + 0.185, x0 + 0.275], [y, y], color=color, linewidth=3.2, alpha=0.82, zorder=6)
    ax.text(x0 + 0.185, y0 + 0.000, "pathway summaries", color=COLORS["muted"], fontsize=6.2, ha="left", zorder=6)


def build_methods_overview() -> dict[str, Path]:
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    fig = plt.figure(figsize=(16, 9), dpi=SLIDE_DPI)
    ax = fig.add_axes([0, 0, 1, 1])
    ax.set_axis_off()
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)
    draw_canvas(ax)

    ax.text(
        0.066,
        0.842,
        "How public studies\nbecome a benchmark",
        color=COLORS["ink"],
        fontsize=23.0,
        fontweight="bold",
        ha="left",
        va="top",
        linespacing=1.02,
        zorder=8,
    )
    ax.text(
        0.068,
        0.675,
        "The key move is simple: hide one mission,\ntrain only on the rest, then audit every output.",
        color=COLORS["muted"],
        fontsize=10.0,
        ha="left",
        va="top",
        linespacing=1.22,
        zorder=8,
    )

    y = 0.585
    xs = [0.145, 0.275, 0.405, 0.535, 0.665, 0.795]
    steps = [
        ("1", "Public studies", "public sample\nand data tables", COLORS["green"]),
        ("2", "Clean sample table", "source records\nand labels", COLORS["blue"]),
        ("3", "Hide a mission", "other missions train\nhidden mission tests", COLORS["red"]),
        ("4", "Build features", "genes and\npathways", COLORS["teal"]),
        ("5", "Train models", "baseline and\nfoundation models", COLORS["purple"]),
        ("6", "Score and audit", "held-out scores\nand audit trail", COLORS["amber"]),
    ]
    for i, (number, title, detail, color) in enumerate(steps):
        draw_step(ax, xs[i], y, number, title, detail, color, emphasize=(number == "3"))
        if i < len(xs) - 1:
            draw_arrow(ax, (xs[i] + 0.033, y), (xs[i + 1] - 0.033, y))

    ax.add_patch(Arc((0.405, 0.595), 0.132, 0.142, theta1=205, theta2=335, color=COLORS["red"], linewidth=1.2, alpha=0.45, zorder=4))

    draw_hidden_mission_inset(ax)
    draw_gene_to_pathway(ax)
    draw_train_only_guard(ax)

    ax.text(
        0.066,
        0.078,
        "Methods explainer | public studies -> hidden mission tasks -> train-only features -> held-out scores -> audited release boundary",
        color=COLORS["muted"],
        fontsize=7.8,
        ha="left",
        alpha=0.90,
        zorder=8,
    )

    png = OUT_DIR / "methods_data_to_evaluation_overview.png"
    pdf = OUT_DIR / "methods_data_to_evaluation_overview.pdf"
    fig.savefig(png, dpi=SLIDE_DPI, facecolor=COLORS["bg"])
    fig.savefig(pdf, facecolor=COLORS["bg"])
    plt.close(fig)

    manifest = {
        "slide_id": "methods_data_to_evaluation_overview",
        "created": "2026-06-01",
        "purpose": "first-time-viewer explanation of data collection, processing, held-out evaluation, and audit flow",
        "visible_claim": "How public studies become a benchmark",
        "must_explain": [
            "public OSDR/GeneLab studies",
            "sample metadata cleanup",
            "one mission held out",
            "train-only feature processing",
            "gene and pathway feature views",
            "model scoring and audit trail",
        ],
        "outputs": {
            "png": str(png.relative_to(ROOT)),
            "pdf": str(pdf.relative_to(ROOT)),
        },
    }
    manifest_path = OUT_DIR / "methods_data_to_evaluation_overview_manifest.json"
    manifest_path.write_text(json.dumps(manifest, indent=2) + "\n")
    return {"png": png, "pdf": pdf, "manifest": manifest_path}


def main() -> None:
    outputs = build_methods_overview()
    print(json.dumps({key: str(path.relative_to(ROOT)) for key, path in outputs.items()}, indent=2))


if __name__ == "__main__":
    main()
