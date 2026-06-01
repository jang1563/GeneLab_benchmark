#!/usr/bin/env python3
"""Render grayscale and approximate color-vision QA sheets for premium figures."""

from __future__ import annotations

import csv
import json
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.image as mpimg
import matplotlib.pyplot as plt
import numpy as np


ROOT = Path(__file__).resolve().parents[1]
FIGURE_DIR = ROOT / "output" / "premium_figures"
OUT_DIR = FIGURE_DIR / "color_qa"


FIGURES = [
    {
        "figure_id": "fig1_tissue_transfer",
        "source": FIGURE_DIR / "premium_fig1_tissue_transfer_hierarchy.png",
    },
    {
        "figure_id": "fig2_pathway_manuscript",
        "source": FIGURE_DIR / "manuscript_variants" / "premium_fig2_pathway_rescue_manuscript.png",
    },
    {
        "figure_id": "fig3_model_manuscript",
        "source": FIGURE_DIR / "manuscript_variants" / "premium_fig3_model_tier_comparison_manuscript.png",
    },
    {
        "figure_id": "fig6_organoid_manuscript",
        "source": FIGURE_DIR / "manuscript_variants" / "premium_fig6_organoid_biology_check_manuscript.png",
    },
]


COLOR_TRANSFORMS = {
    "original": np.eye(3),
    "deuteranopia_approx": np.array(
        [
            [0.625, 0.375, 0.000],
            [0.700, 0.300, 0.000],
            [0.000, 0.300, 0.700],
        ]
    ),
    "protanopia_approx": np.array(
        [
            [0.567, 0.433, 0.000],
            [0.558, 0.442, 0.000],
            [0.000, 0.242, 0.758],
        ]
    ),
}


def load_rgb(path: Path) -> np.ndarray:
    image = mpimg.imread(path).astype(float)
    if image.max() > 1.0:
        image = image / 255.0
    if image.shape[-1] == 4:
        alpha = image[..., 3:4]
        image = image[..., :3] * alpha + (1.0 - alpha)
    return image[..., :3]


def to_grayscale(image: np.ndarray) -> np.ndarray:
    gray = image @ np.array([0.2126, 0.7152, 0.0722])
    return np.repeat(gray[..., None], 3, axis=2)


def apply_matrix(image: np.ndarray, matrix: np.ndarray) -> np.ndarray:
    return np.clip(image @ matrix.T, 0.0, 1.0)


def save_transforms() -> list[dict[str, object]]:
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    rows: list[dict[str, object]] = []
    for item in FIGURES:
        source = Path(item["source"])
        image = load_rgb(source)
        variants = {"grayscale": to_grayscale(image)}
        for name, matrix in COLOR_TRANSFORMS.items():
            variants[name] = apply_matrix(image, matrix)
        for variant_name, variant_image in variants.items():
            output = OUT_DIR / f"{item['figure_id']}_{variant_name}.png"
            plt.imsave(output, variant_image)
            rows.append(
                {
                    "figure_id": item["figure_id"],
                    "variant": variant_name,
                    "source": str(source.relative_to(ROOT)),
                    "output": str(output.relative_to(ROOT)),
                    "width_px": int(variant_image.shape[1]),
                    "height_px": int(variant_image.shape[0]),
                }
            )
    return rows


def write_contact_sheet(rows: list[dict[str, object]]) -> Path:
    by_key = {(row["figure_id"], row["variant"]): Path(row["output"]) for row in rows}
    figure_ids = [item["figure_id"] for item in FIGURES]
    variants = ["original", "grayscale", "deuteranopia_approx", "protanopia_approx"]
    fig = plt.figure(figsize=(13.5, 10.5), dpi=210)
    grid = fig.add_gridspec(
        len(figure_ids),
        len(variants),
        left=0.035,
        right=0.99,
        top=0.93,
        bottom=0.035,
        wspace=0.04,
        hspace=0.18,
    )
    fig.suptitle(
        "Premium figure grayscale and color-vision QA",
        x=0.035,
        y=0.985,
        ha="left",
        fontsize=14,
        fontweight="bold",
    )
    for col, variant in enumerate(variants):
        fig.text(0.035 + col * 0.238, 0.95, variant.replace("_", " "), fontsize=8.5, fontweight="bold")
    for row_idx, figure_id in enumerate(figure_ids):
        for col_idx, variant in enumerate(variants):
            ax = fig.add_subplot(grid[row_idx, col_idx])
            image = mpimg.imread(ROOT / by_key[(figure_id, variant)])
            ax.imshow(image)
            ax.axis("off")
            if col_idx == 0:
                ax.set_title(figure_id, loc="left", fontsize=8, pad=3)
    output = OUT_DIR / "premium_color_qa_contact_sheet.png"
    fig.savefig(output, dpi=210)
    plt.close(fig)
    return output


def main() -> None:
    rows = save_transforms()
    contact_sheet = write_contact_sheet(rows)
    csv_path = OUT_DIR / "premium_color_qa_manifest.csv"
    json_path = OUT_DIR / "premium_color_qa_manifest.json"
    keys = list(rows[0].keys())
    with csv_path.open("w", newline="", encoding="utf-8") as fh:
        writer = csv.DictWriter(fh, fieldnames=keys)
        writer.writeheader()
        writer.writerows(rows)
    with json_path.open("w", encoding="utf-8") as fh:
        json.dump({"transforms": rows, "contact_sheet": str(contact_sheet.relative_to(ROOT))}, fh, indent=2)
        fh.write("\n")
    print(json.dumps({"transforms": len(rows), "contact_sheet": str(contact_sheet.relative_to(ROOT))}, indent=2))


if __name__ == "__main__":
    main()
