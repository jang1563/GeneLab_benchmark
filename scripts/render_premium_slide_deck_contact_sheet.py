#!/usr/bin/env python3
"""Render a deck-order contact sheet for premium slide-scene candidates."""

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


RESULT_DIR = ROOT / "output" / "premium_slide_scenes"
V9_DIR = ROOT / "output" / "premium_v9_document_scenes"
OUT_DIR = ROOT / "output" / "premium_deck_review"


SLIDES = [
    {
        "slide_id": "fig1_tissue_transfer",
        "section": "result",
        "claim": "Some tissues transfer",
        "source": RESULT_DIR / "fig1_tissue_transfer_layered_scene.png",
    },
    {
        "slide_id": "fig2_pathway",
        "section": "result",
        "claim": "Pathways suppress unwanted labels",
        "source": RESULT_DIR / "fig2_pathway_layered_scene.png",
    },
    {
        "slide_id": "fig3_model_tier",
        "section": "result",
        "claim": "Scale alone does not transfer",
        "source": RESULT_DIR / "fig3_model_tier_layered_scene.png",
    },
    {
        "slide_id": "fig6_organoid",
        "section": "extension",
        "claim": "Organoids add biology checks",
        "source": RESULT_DIR / "fig6_organoid_layered_scene.png",
    },
    {
        "slide_id": "v9_platform",
        "section": "resource",
        "claim": "V9 is a staged evidence resource",
        "source": V9_DIR / "v9_platform_provenance_document_scene.png",
    },
    {
        "slide_id": "v9_public_bulk_boundary",
        "section": "resource",
        "claim": "Public bulk is metadata-ready",
        "source": V9_DIR / "v9_public_bulk_boundary_document_scene.png",
    },
]


def main() -> None:
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    fig = plt.figure(figsize=(18.5, 10.0), dpi=210)
    grid = fig.add_gridspec(2, 3, left=0.025, right=0.985, top=0.865, bottom=0.055, wspace=0.030, hspace=0.120)
    fig.suptitle(
        "Premium slide deck-order contact sheet",
        x=0.025,
        y=0.955,
        ha="left",
        fontsize=14,
        fontweight="bold",
    )
    fig.text(
        0.025,
        0.902,
        "Dark result grammar transitions to light provenance-document grammar for v9 resource/release-boundary slides.",
        fontsize=8.4,
        color="#45515F",
        ha="left",
    )
    for idx, slide in enumerate(SLIDES):
        ax = fig.add_subplot(grid[idx // 3, idx % 3])
        ax.imshow(mpimg.imread(slide["source"]))
        ax.set_title(f"{idx + 1}. {slide['slide_id']} | {slide['section']} | {slide['claim']}", loc="left", fontsize=7.0, pad=3)
        ax.axis("off")
    output = OUT_DIR / "premium_slide_deck_order_contact_sheet.png"
    fig.savefig(output, dpi=210, facecolor="white")
    plt.close(fig)

    manifest = {
        "contact_sheet": str(output.relative_to(ROOT)),
        "slides": [
            {
                "slide_id": slide["slide_id"],
                "section": slide["section"],
                "claim": slide["claim"],
                "source": str(slide["source"].relative_to(ROOT)),
            }
            for slide in SLIDES
        ],
    }
    manifest_path = OUT_DIR / "premium_slide_deck_order_contact_sheet_manifest.json"
    manifest_path.write_text(json.dumps(manifest, indent=2) + "\n")
    print(json.dumps({"contact_sheet": str(output.relative_to(ROOT))}, indent=2))


if __name__ == "__main__":
    main()
