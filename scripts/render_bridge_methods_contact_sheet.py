#!/usr/bin/env python3
"""Render a contact sheet for bridge-method slide prototypes."""

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


OUT_DIR = ROOT / "output" / "premium_bridge_scenes"

SLIDES = [
    {
        "slide_id": "b3_mission_held_out",
        "claim": "The test set is a mission",
        "source": OUT_DIR / "b3_mission_held_out" / "rendered_preview.png",
    },
    {
        "slide_id": "b4_train_only_guard",
        "claim": "Feature processing stays on train side",
        "source": OUT_DIR / "b4_train_only_guard" / "rendered_preview.png",
    },
]


def main() -> None:
    fig = plt.figure(figsize=(17.8, 6.2), dpi=220)
    grid = fig.add_gridspec(1, 2, left=0.022, right=0.986, top=0.810, bottom=0.055, wspace=0.030)
    fig.suptitle(
        "Bridge-method slide prototypes: B3-B4",
        x=0.022,
        y=0.945,
        ha="left",
        fontsize=13.2,
        fontweight="bold",
    )
    fig.text(
        0.022,
        0.880,
        "Both slides use the light methods/provenance grammar: paper field, thin rules, one trust boundary, and editable interpretation overlays.",
        fontsize=7.8,
        color="#45515F",
        ha="left",
    )
    for idx, slide in enumerate(SLIDES):
        ax = fig.add_subplot(grid[0, idx])
        ax.imshow(mpimg.imread(slide["source"]))
        ax.set_title(f"{idx + 1}. {slide['slide_id']} | {slide['claim']}", loc="left", fontsize=7.0, pad=3)
        ax.axis("off")
    output = OUT_DIR / "bridge_methods_b3_b4_contact_sheet.png"
    fig.savefig(output, dpi=220, facecolor="white")
    plt.close(fig)

    manifest = {
        "contact_sheet": str(output.relative_to(ROOT)),
        "slides": [
            {
                "slide_id": slide["slide_id"],
                "claim": slide["claim"],
                "source": str(slide["source"].relative_to(ROOT)),
            }
            for slide in SLIDES
        ],
    }
    manifest_path = OUT_DIR / "bridge_methods_b3_b4_contact_sheet_manifest.json"
    manifest_path.write_text(json.dumps(manifest, indent=2) + "\n", encoding="utf-8")
    print(json.dumps({"contact_sheet": str(output.relative_to(ROOT))}, indent=2))


if __name__ == "__main__":
    main()
