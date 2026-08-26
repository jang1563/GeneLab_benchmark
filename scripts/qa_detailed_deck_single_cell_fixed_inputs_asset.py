#!/usr/bin/env python3
"""Build strict visual QA crops for the single-cell fixed-inputs asset."""

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
    / "single_cell_fixed_inputs"
)
SLIDE = ASSET_DIR / "single_cell_scoring_needs_fixed_inputs_premium.png"
QA_DIR = ASSET_DIR / "qa_crops"
STRICT_DIR = ASSET_DIR / "strict_qa"
QA_DIR.mkdir(parents=True, exist_ok=True)
STRICT_DIR.mkdir(parents=True, exist_ok=True)


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


F_TITLE = font(44, bold=True)


def fit(image: Image.Image, size: tuple[int, int]) -> Image.Image:
    image = image.convert("RGB").copy()
    resample = getattr(getattr(Image, "Resampling", Image), "LANCZOS")
    image.thumbnail(size, resample)
    canvas = Image.new("RGB", size, "#0C111A")
    canvas.paste(image, ((size[0] - image.width) // 2, (size[1] - image.height) // 2))
    return canvas


def crop_set(target_dir: Path, specs: list[tuple[str, tuple[int, int, int, int]]], filename: str) -> Path:
    slide = Image.open(SLIDE).convert("RGB")
    entries = []
    for name, box in specs:
        crop = slide.crop(box)
        out = target_dir / f"single_cell_fixed_inputs_{name}.png"
        crop.save(out, quality=95)
        entries.append({"name": name, "box": list(box), "path": str(out.relative_to(ROOT))})

    thumb_w, thumb_h = 1480, 620
    pad_x, pad_y = 70, 88
    width = thumb_w + pad_x * 2
    height = 56 + len(specs) * (thumb_h + pad_y) + 40
    canvas = Image.new("RGB", (width, height), "#0C111A")
    draw = ImageDraw.Draw(canvas)
    y = 28
    for name, _box in specs:
        crop = Image.open(target_dir / f"single_cell_fixed_inputs_{name}.png")
        draw.text((pad_x, y), name.replace("_", " "), font=F_TITLE, fill="#F3F7FC")
        y += 54
        draw.rounded_rectangle((pad_x - 6, y - 6, pad_x + thumb_w + 6, y + thumb_h + 6), radius=18, outline="#2A394D", width=2)
        canvas.paste(fit(crop, (thumb_w, thumb_h)), (pad_x, y))
        y += thumb_h + pad_y
    sheet = target_dir / filename
    canvas.save(sheet, quality=95)
    (target_dir / "single_cell_fixed_inputs_crop_manifest.json").write_text(json.dumps(entries, indent=2) + "\n")
    return sheet


def main() -> None:
    qa_specs = [
        ("full_thumbnail", (0, 0, 3840, 2160)),
        ("title_badges", (90, 42, 3745, 520)),
        ("metric_contract_panel", (90, 580, 1355, 1495)),
        ("fixed_input_ladder", (1335, 580, 2565, 1495)),
        ("evaluator_panel", (2545, 580, 3745, 1495)),
        ("run_path", (90, 1500, 3745, 1875)),
        ("footer", (90, 1880, 3745, 2105)),
    ]
    strict_specs = [
        ("01_full_thumbnail", (0, 0, 3840, 2160)),
        ("02_header_title_badges", (90, 42, 3745, 520)),
        ("03_metric_footprint", (555, 820, 1220, 1175)),
        ("04_metric_profile_cards", (160, 1185, 1260, 1410)),
        ("05_input_ladder_rows", (1470, 810, 2470, 1358)),
        ("06_evaluator_table", (2625, 820, 3655, 1085)),
        ("07_input_field_rows", (2625, 1100, 3655, 1318)),
        ("08_run_path_arrows", (230, 1608, 3605, 1748)),
        ("09_run_path_chips", (230, 1740, 2795, 1835)),
        ("10_footer_takeaway", (120, 1910, 3720, 2076)),
    ]
    qa_sheet = crop_set(QA_DIR, qa_specs, "single_cell_fixed_inputs_qa_crop_contact_sheet.png")
    strict_sheet = crop_set(STRICT_DIR, strict_specs, "single_cell_fixed_inputs_STRICT_QA_contact_sheet.png")
    print(json.dumps({"qa_contact_sheet": str(qa_sheet.relative_to(ROOT)), "strict_contact_sheet": str(strict_sheet.relative_to(ROOT))}, indent=2))


if __name__ == "__main__":
    main()
