#!/usr/bin/env python3
"""Render manuscript-width previews from premium figure PNG exports.

These previews do not replace source figures. They answer one narrow QA
question: how does the current 16:9 figure art behave when shown at a
two-column manuscript width?
"""

from __future__ import annotations

import csv
import json
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.image as mpimg
import matplotlib.pyplot as plt


ROOT = Path(__file__).resolve().parents[1]
FIGURE_DIR = ROOT / "output" / "premium_figures"
OUT_DIR = ROOT / "output" / "premium_figures" / "manuscript_previews"


FIGURES = [
    {
        "figure_id": "fig1_tissue_transfer",
        "source": FIGURE_DIR / "premium_fig1_tissue_transfer_hierarchy.png",
        "role": "main_candidate",
        "intended_width": "two_column",
    },
    {
        "figure_id": "fig2_pathway",
        "source": FIGURE_DIR / "premium_fig2_pathway_artifact_rescue.png",
        "role": "main_candidate",
        "intended_width": "two_column",
    },
    {
        "figure_id": "fig3_model",
        "source": FIGURE_DIR / "premium_fig3_model_tier_comparison.png",
        "role": "main_candidate_deck_first",
        "intended_width": "two_column",
    },
    {
        "figure_id": "fig4_platform_architecture",
        "source": FIGURE_DIR / "premium_fig4_v9_platform_architecture.png",
        "role": "resource_overview",
        "intended_width": "two_column",
    },
    {
        "figure_id": "fig5_release_boundary",
        "source": FIGURE_DIR / "premium_fig5_public_bulk_release_boundary_schematic.png",
        "role": "status_or_supplement",
        "intended_width": "two_column_or_slide",
    },
    {
        "figure_id": "fig6_organoid",
        "source": FIGURE_DIR / "premium_fig6_organoid_diagnostic_surface.png",
        "role": "main_candidate_or_extension",
        "intended_width": "two_column",
    },
]


def render_preview(source: Path, output: Path, width_in: float = 7.2, dpi: int = 300) -> dict[str, object]:
    image = mpimg.imread(source)
    height_px, width_px = image.shape[:2]
    height_in = width_in * height_px / width_px
    fig = plt.figure(figsize=(width_in, height_in), dpi=dpi)
    ax = fig.add_axes([0, 0, 1, 1])
    ax.imshow(image)
    ax.axis("off")
    fig.savefig(output, dpi=dpi)
    plt.close(fig)
    return {
        "source_width_px": width_px,
        "source_height_px": height_px,
        "preview_width_in": width_in,
        "preview_height_in": round(height_in, 3),
        "preview_dpi": dpi,
        "preview_width_px": int(round(width_in * dpi)),
        "preview_height_px": int(round(height_in * dpi)),
    }


def write_contact_sheet(rows: list[dict[str, object]]) -> Path:
    images = [(row, mpimg.imread(row["preview_path"])) for row in rows]
    fig = plt.figure(figsize=(11, 8.5), dpi=220)
    grid = fig.add_gridspec(3, 2, left=0.04, right=0.98, top=0.92, bottom=0.05, wspace=0.1, hspace=0.22)
    fig.suptitle("Premium figure manuscript-width previews", x=0.04, y=0.975, ha="left", fontsize=14, fontweight="bold")
    for idx, (row, image) in enumerate(images):
        ax = fig.add_subplot(grid[idx // 2, idx % 2])
        ax.imshow(image)
        ax.axis("off")
        ax.set_title(f"{row['figure_id']} ({row['role']})", loc="left", fontsize=8, pad=4)
    output = OUT_DIR / "manuscript_preview_contact_sheet.png"
    fig.savefig(output, dpi=220)
    plt.close(fig)
    return output


def main() -> None:
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    rows: list[dict[str, object]] = []
    for item in FIGURES:
        source = Path(item["source"])
        output = OUT_DIR / f"{item['figure_id']}_two_column_7p2in.png"
        metrics = render_preview(source, output)
        rows.append(
            {
                "figure_id": item["figure_id"],
                "role": item["role"],
                "intended_width": item["intended_width"],
                "source": str(source.relative_to(ROOT)),
                "preview_path": str(output.relative_to(ROOT)),
                **metrics,
            }
        )
    contact_sheet = write_contact_sheet(rows)
    csv_path = OUT_DIR / "manuscript_preview_manifest.csv"
    json_path = OUT_DIR / "manuscript_preview_manifest.json"
    keys = list(rows[0].keys())
    with csv_path.open("w", newline="", encoding="utf-8") as fh:
        writer = csv.DictWriter(fh, fieldnames=keys)
        writer.writeheader()
        writer.writerows(rows)
    with json_path.open("w", encoding="utf-8") as fh:
        json.dump({"previews": rows, "contact_sheet": str(contact_sheet.relative_to(ROOT))}, fh, indent=2)
        fh.write("\n")
    print(json.dumps({"rendered": len(rows), "contact_sheet": str(contact_sheet.relative_to(ROOT))}, indent=2))


if __name__ == "__main__":
    main()
