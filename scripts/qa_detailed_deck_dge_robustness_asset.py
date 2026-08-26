#!/usr/bin/env python3
"""Build strict visual QA crops for the DGE robustness asset."""

from __future__ import annotations

import json
from pathlib import Path

from PIL import Image, ImageDraw, ImageFont


ROOT = Path(__file__).resolve().parent.parent
ASSET_DIR = (
    ROOT
    / "outputs"
    / "019e8bd1-3b42-7821-9e1c-beebfa8e2ece"
    / "presentations"
    / "spacebiobench-detailed-deck"
    / "assets"
    / "dge_robustness"
)
SOURCE = ASSET_DIR / "dge_pipeline_robustness_premium.png"
QA_DIR = ASSET_DIR / "qa_crops_v2"
STRICT_DIR = ASSET_DIR / "strict_qa"


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


F = {"label": font(28, bold=True)}


def fit(image: Image.Image, size: tuple[int, int]) -> Image.Image:
    image = image.convert("RGB").copy()
    resample = getattr(getattr(Image, "Resampling", Image), "LANCZOS")
    image.thumbnail(size, resample)
    result = Image.new("RGB", size, "#0C111A")
    result.paste(image, ((size[0] - image.width) // 2, (size[1] - image.height) // 2))
    return result


def write_contact_sheet(crops: list[dict[str, object]], directory: Path, sheet_name: str) -> Path:
    width = 1700
    thumb_w = 1560
    y = 42
    blocks: list[tuple[str, Image.Image, int]] = []
    for crop in crops:
        name = str(crop["name"])
        image = Image.open(directory / f"dge_robustness_{name}.png")
        target_h = 420 if image.width > image.height * 1.8 else 620
        fitted = fit(image, (thumb_w, target_h))
        blocks.append((name, fitted, y))
        y += fitted.height + 98

    canvas = Image.new("RGB", (width, y + 20), "#0C111A")
    draw = ImageDraw.Draw(canvas)
    for name, image, top in blocks:
        draw.text((70, top - 32), name.replace("_", " "), font=F["label"], fill="#F3F7FC")
        draw.rounded_rectangle((66, top - 2, 1634, top + image.height + 2), radius=18, outline="#2A394D", width=2)
        canvas.paste(image, (70, top))
    out = directory / sheet_name
    canvas.save(out, quality=95)
    return out


def crop_set(target_dir: Path, specs: list[dict[str, object]], sheet_name: str) -> Path:
    target_dir.mkdir(parents=True, exist_ok=True)
    source = Image.open(SOURCE).convert("RGB")
    manifest = []
    for spec in specs:
        name = str(spec["name"])
        box = tuple(int(v) for v in spec["box"])
        out = target_dir / f"dge_robustness_{name}.png"
        source.crop(box).save(out, quality=95)
        manifest.append({"name": name, "box": box, "path": str(out.relative_to(ROOT))})
    sheet = write_contact_sheet(specs, target_dir, sheet_name)
    (target_dir / "dge_robustness_crop_manifest.json").write_text(
        json.dumps({"source": str(SOURCE.relative_to(ROOT)), "crops": manifest, "contact_sheet": str(sheet.relative_to(ROOT))}, indent=2) + "\n"
    )
    return sheet


def main() -> None:
    qa_specs = [
        {"name": "01_header_badges", "box": (120, 40, 3720, 270)},
        {"name": "02_workflow_left", "box": (130, 320, 1150, 1830)},
        {"name": "03_rank_panel", "box": (1140, 320, 2515, 1210)},
        {"name": "04_rank_rows_axis", "box": (1530, 490, 2515, 1125)},
        {"name": "05_list_panel", "box": (1140, 1200, 2515, 1830)},
        {"name": "06_list_rows_axis", "box": (1530, 1340, 2515, 1765)},
        {"name": "07_readout_right", "box": (2505, 320, 3720, 1830)},
        {"name": "08_footer", "box": (130, 1860, 3710, 2070)},
        {"name": "09_full_downscaled", "box": (0, 0, 3840, 2160)},
    ]
    strict_specs = [
        {"name": "01_full_thumbnail", "box": (0, 0, 3840, 2160)},
        {"name": "02_header_title_badges", "box": (120, 40, 3720, 270)},
        {"name": "03_workflow_contrast_callers", "box": (180, 500, 1090, 1135)},
        {"name": "04_workflow_outputs_readout", "box": (180, 1270, 1090, 1775)},
        {"name": "05_rank_stat", "box": (1210, 460, 1585, 720)},
        {"name": "06_rank_rows_axis", "box": (1530, 490, 2515, 1125)},
        {"name": "07_list_stat", "box": (1210, 1340, 1585, 1605)},
        {"name": "08_list_rows_axis", "box": (1530, 1340, 2515, 1765)},
        {"name": "09_right_layer_cards", "box": (2585, 500, 3640, 1085)},
        {"name": "10_right_scope_box", "box": (2585, 1090, 3640, 1415)},
        {"name": "11_footer", "box": (150, 1880, 3690, 2045)},
    ]
    qa_sheet = crop_set(QA_DIR, qa_specs, "dge_robustness_qa_crop_contact_sheet.png")
    strict_sheet = crop_set(STRICT_DIR, strict_specs, "dge_robustness_STRICT_QA_contact_sheet.png")
    print(json.dumps({"qa_contact_sheet": str(qa_sheet.relative_to(ROOT)), "strict_contact_sheet": str(strict_sheet.relative_to(ROOT))}, indent=2))


if __name__ == "__main__":
    main()
