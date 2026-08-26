#!/usr/bin/env python3
"""Render the detailed-deck slides 7-21 assembly contact sheet."""

from __future__ import annotations

import json
from pathlib import Path

from PIL import Image, ImageDraw, ImageFont


ROOT = Path(__file__).resolve().parent.parent
WORKSPACE = (
    ROOT
    / "outputs"
    / "019e8bd1-3b42-7821-9e1c-beebfa8e2ece"
    / "presentations"
    / "spacebiobench-detailed-deck"
)
ASSET_DIR = WORKSPACE / "assets"
OUT_DIR = ASSET_DIR / "act1_2_assembly"
OUT_DIR.mkdir(parents=True, exist_ok=True)

SLIDES = [
    {
        "slide": 7,
        "act": "Method",
        "title": "The Test Set Is An Entire Mission",
        "role": "hidden-unit rule",
        "path": ASSET_DIR / "mission_holdout" / "mission_heldout_protocol_premium.png",
    },
    {
        "slide": 8,
        "act": "Method",
        "title": "Training Choices Stop Before The Hidden Mission",
        "role": "leakage-control rule",
        "path": ASSET_DIR / "train_only_guard" / "train_only_feature_guard_premium.png",
    },
    {
        "slide": 9,
        "act": "Method",
        "title": "How To Read AUROC",
        "role": "metric reader",
        "path": ASSET_DIR / "metric_primer" / "metric_primer_auroc_uncertainty_premium.png",
    },
    {
        "slide": 10,
        "act": "Method",
        "title": "One Task Contract, Three Model Families",
        "role": "family bridge",
        "path": ASSET_DIR / "model_family_bridge" / "model_family_bridge_premium.png",
    },
    {
        "slide": 11,
        "act": "Method",
        "title": "Same Benchmark, Three Input Surfaces",
        "role": "input surface primer",
        "path": ASSET_DIR / "model_family" / "model_family_input_surface_premium.png",
    },
    {
        "slide": 12,
        "act": "Method",
        "title": "What Counts As Evidence",
        "role": "scope/readout legend",
        "path": ASSET_DIR / "evidence_scope_ladder" / "evidence_scope_ladder_premium.png",
    },
    {
        "slide": 13,
        "act": "Core Result",
        "title": "Read One Tissue Score",
        "role": "worked result row",
        "path": ASSET_DIR / "worked_tissue_score" / "worked_tissue_score_thymus_transfer_premium.png",
    },
    {
        "slide": 14,
        "act": "Core Result",
        "title": "Thymus Leads The Transfer Hierarchy",
        "role": "tissue ranking",
        "path": ASSET_DIR / "tissue_hierarchy" / "tissue_transfer_hierarchy_premium.png",
    },
    {
        "slide": 15,
        "act": "Core Result",
        "title": "More Liver Missions Expose Heterogeneity",
        "role": "counter-example explainer",
        "path": ASSET_DIR / "liver_heterogeneity" / "liver_mission_heterogeneity_premium.png",
    },
    {
        "slide": 16,
        "act": "Core Result",
        "title": "The Transfer Matrix Behind The Ranking",
        "role": "pair-level proof",
        "path": ASSET_DIR / "transfer_matrix" / "transfer_matrix_behind_ranking_premium.png",
    },
    {
        "slide": 17,
        "act": "Core Result",
        "title": "Pathway Features Reduce Selected Nuisance Signals",
        "role": "feature-view diagnostic",
        "path": ASSET_DIR / "pathway_nuisance" / "pathway_features_reduce_selected_nuisance_premium.png",
    },
    {
        "slide": 18,
        "act": "Core Result",
        "title": "Conserved Pathways Predict Transfer",
        "role": "pathway conservation proof",
        "path": ASSET_DIR / "nes_conservation" / "nes_conservation_predicts_transfer_premium.png",
    },
    {
        "slide": 19,
        "act": "Core Result",
        "title": "Screen Transfer Feasibility Before Training",
        "role": "preflight triage workflow",
        "path": ASSET_DIR / "transfer_feasibility" / "screen_transfer_feasibility_before_training_premium.png",
    },
    {
        "slide": 20,
        "act": "Core Result",
        "title": "Held-Out Missions Confirm The Signal",
        "role": "independent mission check",
        "path": ASSET_DIR / "heldout_validation" / "heldout_missions_confirm_signal_premium.png",
    },
    {
        "slide": 21,
        "act": "Core Result",
        "title": "Negative Controls Anchor The Readout",
        "role": "control gate",
        "path": ASSET_DIR / "negative_controls" / "negative_controls_anchor_readout_premium.png",
    },
]

COLORS = {
    "bg": "#0C111A",
    "panel": "#101823",
    "line": "#2A394D",
    "text": "#F3F7FC",
    "muted": "#A8B4C4",
    "teal": "#5FD3C4",
    "sky": "#73A7FF",
    "amber": "#F4C26B",
    "green": "#84D278",
}


def font(size: int, *, bold: bool = False) -> ImageFont.ImageFont:
    candidates = [
        "/System/Library/Fonts/Supplemental/Arial Bold.ttf" if bold else "/System/Library/Fonts/Supplemental/Arial.ttf",
        "/System/Library/Fonts/Supplemental/Helvetica Bold.ttf" if bold else "/System/Library/Fonts/Supplemental/Helvetica.ttf",
    ]
    for candidate in candidates:
        try:
            return ImageFont.truetype(candidate, size)
        except OSError:
            continue
    return ImageFont.load_default()


F = {
    "title": font(54, bold=True),
    "subtitle": font(25),
    "label": font(28, bold=True),
    "small": font(22),
    "tiny": font(18),
    "num": font(34, bold=True),
}


def fit(image: Image.Image, size: tuple[int, int]) -> Image.Image:
    image = image.convert("RGB").copy()
    resample = getattr(getattr(Image, "Resampling", Image), "LANCZOS")
    image.thumbnail(size, resample)
    result = Image.new("RGB", size, COLORS["bg"])
    result.paste(image, ((size[0] - image.width) // 2, (size[1] - image.height) // 2))
    return result


def main() -> None:
    thumb_w, thumb_h = 760, 428
    label_h = 98
    cols = 3
    rows = (len(SLIDES) + cols - 1) // cols
    pad_x = 72
    gap_x = 52
    gap_y = 58
    header_h = 190
    footer_h = 130
    width = pad_x * 2 + cols * thumb_w + (cols - 1) * gap_x
    height = header_h + rows * (thumb_h + label_h) + (rows - 1) * gap_y + footer_h

    canvas = Image.new("RGB", (width, height), COLORS["bg"])
    draw = ImageDraw.Draw(canvas)
    draw.text((pad_x, 36), "Detailed deck Act 1-2 assembly preview", font=F["title"], fill=COLORS["text"])
    draw.text(
        (pad_x, 104),
        "Slides 7-21: hidden mission -> leakage guard -> reader primers -> first result and control sequence.",
        font=F["subtitle"],
        fill=COLORS["muted"],
    )
    draw.rounded_rectangle((pad_x, 142, pad_x + 615, 180), radius=16, fill="#121B28", outline=COLORS["line"], width=1)
    draw.text((pad_x + 22, 151), "METHOD: slides 7-12", font=F["tiny"], fill=COLORS["teal"])
    draw.rounded_rectangle((pad_x + 650, 142, pad_x + 1260, 180), radius=16, fill="#121B28", outline=COLORS["line"], width=1)
    draw.text((pad_x + 672, 151), "CORE RESULT ENTRY: slides 13-21", font=F["tiny"], fill=COLORS["amber"])

    manifest_entries = []
    for i, slide in enumerate(SLIDES):
        row, col = divmod(i, cols)
        x = pad_x + col * (thumb_w + gap_x)
        y = header_h + row * (thumb_h + label_h + gap_y)
        color = COLORS["teal"] if slide["act"] == "Method" else COLORS["amber"]
        image = fit(Image.open(slide["path"]), (thumb_w, thumb_h))
        draw.rounded_rectangle((x - 4, y - 4, x + thumb_w + 4, y + thumb_h + 4), radius=18, outline=COLORS["line"], width=2)
        canvas.paste(image, (x, y))
        badge_w = 78
        draw.rounded_rectangle((x + 16, y + 16, x + 16 + badge_w, y + 62), radius=18, fill=color)
        draw.text((x + 38, y + 22), f"{slide['slide']:02d}", font=F["num"], fill="#081018")
        draw.text((x, y + thumb_h + 18), f"{slide['slide']}. {slide['title']}", font=F["label"], fill=COLORS["text"])
        draw.text((x, y + thumb_h + 53), f"{slide['act']} / {slide['role']}", font=F["small"], fill=color)
        manifest_entries.append(
            {
                "slide": slide["slide"],
                "act": slide["act"],
                "title": slide["title"],
                "role": slide["role"],
                "path": str(slide["path"].relative_to(ROOT)),
            }
        )

    footer_y = height - footer_h + 26
    draw.text((pad_x, footer_y), "QA focus", font=F["label"], fill=COLORS["text"])
    draw.text(
        (pad_x + 135, footer_y + 5),
        "Check rhythm, visual density, and whether held-out validation flows into the control gate.",
        font=F["small"],
        fill=COLORS["muted"],
    )

    out = OUT_DIR / "act1_2_ready_sequence_contact_sheet.png"
    canvas.save(out, quality=95)
    manifest = OUT_DIR / "act1_2_ready_sequence_contact_sheet_manifest.json"
    manifest.write_text(
        json.dumps({"contact_sheet": str(out.relative_to(ROOT)), "slides": manifest_entries}, indent=2) + "\n"
    )
    print(json.dumps({"contact_sheet": str(out.relative_to(ROOT)), "manifest": str(manifest.relative_to(ROOT))}, indent=2))


if __name__ == "__main__":
    main()
