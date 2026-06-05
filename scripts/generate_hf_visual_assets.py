#!/usr/bin/env python3
"""Generate static visual assets for the Hugging Face dataset card."""

from __future__ import annotations

import os
import tempfile
from pathlib import Path

os.environ.setdefault("MPLCONFIGDIR", str(Path(tempfile.gettempdir()) / "matplotlib"))

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch


BASE_DIR = Path(__file__).resolve().parent.parent
OUT_DIR = BASE_DIR / "docs" / "assets"
OUT_PATH = OUT_DIR / "hf_benchmark_summary.png"

TISSUE_RESULTS = [
    ("Thymus", 0.948, "PCA-LR / KEGG", True),
    ("Colon", 0.921, "PCA-LR / KEGG", True),
    ("Lung", 0.901, "PCA-LR / Gene", True),
    ("Kidney", 0.829, "ElasticNet-LR / Hallmark", True),
    ("Eye", 0.823, "PCA-LR / Hallmark", False),
    ("Skin", 0.819, "PCA-LR / Gene", False),
    ("Gastrocnemius", 0.776, "PCA-LR / Gene", False),
    ("Liver", 0.670, "PCA-LR / Gene", False),
]

METRICS = [
    ("4", "public GO tasks"),
    ("15", "fold packages"),
    ("8", "tissues"),
    ("600+", "processed samples"),
    ("256", "ML evaluations"),
]


def add_card(ax, xy, width, height, radius=0.025, fc="#ffffff", ec="#d6dde8", lw=1.2):
    patch = FancyBboxPatch(
        xy,
        width,
        height,
        boxstyle=f"round,pad=0.012,rounding_size={radius}",
        linewidth=lw,
        edgecolor=ec,
        facecolor=fc,
        transform=ax.transAxes,
        clip_on=False,
    )
    ax.add_patch(patch)
    return patch


def draw_metric_cards(ax):
    left = 0.055
    gap = 0.018
    width = (0.89 - gap * (len(METRICS) - 1)) / len(METRICS)
    for idx, (value, label) in enumerate(METRICS):
        x = left + idx * (width + gap)
        add_card(ax, (x, 0.625), width, 0.13, radius=0.018, fc="#ffffff")
        ax.text(
            x + 0.025,
            0.70,
            value,
            transform=ax.transAxes,
            ha="left",
            va="center",
            fontsize=28,
            fontweight="bold",
            color="#14202e",
        )
        ax.text(
            x + 0.025,
            0.65,
            label,
            transform=ax.transAxes,
            ha="left",
            va="center",
            fontsize=11.5,
            color="#465465",
        )


def draw_right_panel(ax):
    add_card(ax, (0.635, 0.055), 0.31, 0.53, radius=0.02, fc="#ffffff")
    ax.text(
        0.665,
        0.535,
        "Review-ready package boundary",
        transform=ax.transAxes,
        fontsize=17,
        fontweight="bold",
        color="#14202e",
    )

    bullets = [
        ("Direct HF folds", "features, labels, metadata, fold_info, selected_genes"),
        ("Leakage guard", "gene selection is recomputed on training missions only"),
        ("Held-out validation", "RR-23 thymus AUROC 0.905; RR-7 skin AUROC 0.885"),
        ("Viewer disabled", "high-dimensional matrices are accessed via downloads"),
        ("Stable A6 path", "fold_OSD-397_test replaces stale fold_TBD_test"),
    ]
    y = 0.48
    for title, body in bullets:
        ax.text(
            0.665,
            y,
            title,
            transform=ax.transAxes,
            fontsize=12.8,
            fontweight="bold",
            color="#245b73",
        )
        ax.text(
            0.665,
            y - 0.033,
            body,
            transform=ax.transAxes,
            fontsize=11.3,
            color="#465465",
            wrap=True,
        )
        y -= 0.081


def draw_figure():
    fig = plt.figure(figsize=(16, 9), dpi=150)
    fig.patch.set_facecolor("#f4f7fb")
    ax = fig.add_axes([0, 0, 1, 1])
    ax.set_axis_off()

    ax.text(
        0.055,
        0.93,
        "GeneLab Spaceflight Transcriptomics Benchmark",
        transform=ax.transAxes,
        fontsize=29,
        fontweight="bold",
        color="#101820",
    )
    ax.text(
        0.055,
        0.89,
        "Cross-mission detection of spaceflight transcriptomic signatures from public NASA OSDR mouse RNA-seq.",
        transform=ax.transAxes,
        fontsize=14.5,
        color="#465465",
    )
    ax.text(
        0.055,
        0.842,
        "Dataset freeze: 2026-03-01 | Maintainer/citation: Jihoon Kim, Weill Cornell Medicine",
        transform=ax.transAxes,
        fontsize=12,
        color="#687789",
    )
    ax.text(
        0.055,
        0.812,
        "Public fold package: self-contained features, labels, metadata, and provenance",
        transform=ax.transAxes,
        fontsize=12,
        color="#687789",
    )

    draw_metric_cards(ax)

    # Bar chart card.
    add_card(ax, (0.055, 0.055), 0.54, 0.53, radius=0.02, fc="#ffffff")
    ax.text(
        0.085,
        0.535,
        "Best AUROC by tissue",
        transform=ax.transAxes,
        fontsize=17,
        fontweight="bold",
        color="#14202e",
    )
    ax.text(
        0.085,
        0.505,
        "Best method-feature row by tissue; detailed methods are tabulated below.",
        transform=ax.transAxes,
        fontsize=11.3,
        color="#687789",
    )

    bx = fig.add_axes([0.115, 0.125, 0.405, 0.335])
    tissues = [row[0] for row in reversed(TISSUE_RESULTS)]
    aurocs = [row[1] for row in reversed(TISSUE_RESULTS)]
    significant = [row[3] for row in reversed(TISSUE_RESULTS)]
    colors = ["#1b7f8a" if sig else "#8ea0b5" for sig in significant]
    bx.barh(range(len(tissues)), aurocs, color=colors, height=0.62)
    bx.axvline(0.776, color="#d97b34", linewidth=1.7, linestyle="--")
    bx.text(0.781, 7.48, "PCA-LR mean 0.776", color="#a65f24", fontsize=10.8)
    bx.set_xlim(0.55, 1.0)
    bx.set_yticks(range(len(tissues)))
    bx.set_yticklabels(tissues, fontsize=11.7)
    bx.set_xlabel("AUROC", fontsize=11.6, color="#465465")
    bx.tick_params(axis="x", labelsize=10.8, colors="#465465")
    bx.tick_params(axis="y", length=0, colors="#14202e")
    bx.grid(axis="x", color="#e3e8ef", linewidth=1)
    bx.set_axisbelow(True)
    for spine in bx.spines.values():
        spine.set_visible(False)
    for idx, value in enumerate(aurocs):
        bx.text(value + 0.008, idx, f"{value:.3f}", va="center", fontsize=11.2, color="#14202e")

    draw_right_panel(ax)

    OUT_DIR.mkdir(parents=True, exist_ok=True)
    fig.savefig(OUT_PATH, facecolor=fig.get_facecolor(), bbox_inches="tight", pad_inches=0.08)
    plt.close(fig)


if __name__ == "__main__":
    draw_figure()
    print(OUT_PATH)
