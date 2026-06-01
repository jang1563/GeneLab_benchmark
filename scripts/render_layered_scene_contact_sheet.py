#!/usr/bin/env python3
"""Render a contact sheet for layered slide-scene family review."""

from __future__ import annotations

import json
import os
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
os.environ.setdefault("MPLCONFIGDIR", str(ROOT / "output" / ".matplotlib"))

import matplotlib

matplotlib.use("Agg")
import matplotlib.image as mpimg
import matplotlib.pyplot as plt


SCENE_DIR = ROOT / "output" / "premium_slide_scenes"


SCENES = [
    {
        "scene_id": "fig1_tissue_transfer",
        "claim": "Some tissues transfer",
        "scene_plate": SCENE_DIR / "fig1_tissue_transfer_scene_plate.png",
        "final": SCENE_DIR / "fig1_tissue_transfer_layered_scene.png",
    },
    {
        "scene_id": "fig2_pathway",
        "claim": "Pathways suppress unwanted labels",
        "scene_plate": SCENE_DIR / "fig2_pathway_scene_plate.png",
        "final": SCENE_DIR / "fig2_pathway_layered_scene.png",
    },
    {
        "scene_id": "fig3_model_tier",
        "claim": "Scale alone does not transfer",
        "scene_plate": SCENE_DIR / "fig3_model_tier_scene_plate.png",
        "final": SCENE_DIR / "fig3_model_tier_layered_scene.png",
    },
    {
        "scene_id": "fig6_organoid",
        "claim": "Organoids add biology checks",
        "scene_plate": SCENE_DIR / "fig6_organoid_scene_plate.png",
        "final": SCENE_DIR / "fig6_organoid_layered_scene.png",
    },
]


def write_contact_sheet() -> Path:
    fig = plt.figure(figsize=(19.2, 9.5), dpi=210)
    grid = fig.add_gridspec(
        2,
        len(SCENES),
        left=0.025,
        right=0.99,
        top=0.90,
        bottom=0.045,
        wspace=0.035,
        hspace=0.105,
    )
    fig.suptitle(
        "Layered slide-scene family review",
        x=0.025,
        y=0.965,
        ha="left",
        fontsize=14,
        fontweight="bold",
    )
    fig.text(0.025, 0.915, "Final composite", fontsize=8.5, fontweight="bold")
    fig.text(0.025, 0.472, "Scene plate only", fontsize=8.5, fontweight="bold")

    for col, item in enumerate(SCENES):
        for row, key in enumerate(["final", "scene_plate"]):
            ax = fig.add_subplot(grid[row, col])
            image = mpimg.imread(item[key])
            ax.imshow(image)
            ax.axis("off")
            if row == 0:
                ax.set_title(f"{item['scene_id']} | {item['claim']}", loc="left", fontsize=7.6, pad=3)

    output = SCENE_DIR / "layered_scene_family_contact_sheet.png"
    fig.savefig(output, dpi=210, facecolor="white")
    plt.close(fig)
    return output


def main() -> None:
    SCENE_DIR.mkdir(parents=True, exist_ok=True)
    output = write_contact_sheet()
    manifest = {
        "contact_sheet": str(output.relative_to(ROOT)),
        "scenes": [
            {
                "scene_id": item["scene_id"],
                "claim": item["claim"],
                "scene_plate": str(item["scene_plate"].relative_to(ROOT)),
                "final": str(item["final"].relative_to(ROOT)),
            }
            for item in SCENES
        ],
    }
    manifest_path = SCENE_DIR / "layered_scene_family_contact_sheet_manifest.json"
    manifest_path.write_text(json.dumps(manifest, indent=2) + "\n")
    print(json.dumps({"contact_sheet": str(output.relative_to(ROOT))}, indent=2))


if __name__ == "__main__":
    main()
