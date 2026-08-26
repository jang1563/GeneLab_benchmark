#!/usr/bin/env python3
"""Build strict visual QA crops for the mission-held-out methodology asset."""

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
    / "mission_holdout"
)
SOURCE = ASSET_DIR / "mission_heldout_protocol_premium.png"

QA_DIR = ASSET_DIR / "qa_crops"
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


F = {
    "label": font(28, bold=True),
    "small": font(20),
}


def fit(image: Image.Image, size: tuple[int, int]) -> Image.Image:
    image = image.convert("RGB").copy()
    resample = getattr(getattr(Image, "Resampling", Image), "LANCZOS")
    image.thumbnail(size, resample)
    result = Image.new("RGB", size, "#0C111A")
    x = (size[0] - image.width) // 2
    y = (size[1] - image.height) // 2
    result.paste(image, (x, y))
    return result


def write_contact_sheet(crops: list[dict[str, object]], directory: Path, sheet_name: str) -> Path:
    thumbs: list[tuple[str, Image.Image]] = []
    for crop in crops:
        name = str(crop["name"])
        path = directory / f"mission_holdout_{name}.png"
        thumbs.append((name, Image.open(path)))

    width = 1600
    thumb_w = 1460
    y = 42
    blocks: list[tuple[str, Image.Image, int]] = []
    for name, image in thumbs:
        fitted = fit(image, (thumb_w, 420 if image.width > image.height * 1.8 else 560))
        blocks.append((name, fitted, y))
        y += fitted.height + 96

    canvas = Image.new("RGB", (width, y + 20), "#0C111A")
    draw = ImageDraw.Draw(canvas)
    for name, image, top in blocks:
        draw.text((70, top - 32), name.replace("_", " "), font=F["label"], fill="#F3F7FC")
        draw.rounded_rectangle((66, top - 2, 1534, top + image.height + 2), radius=18, outline="#2A394D", width=2)
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
        out = target_dir / f"mission_holdout_{name}.png"
        source.crop(box).save(out, quality=95)
        manifest.append({"name": name, "box": box, "path": str(out.relative_to(ROOT))})
    sheet = write_contact_sheet(specs, target_dir, sheet_name)
    (target_dir / "mission_holdout_crop_manifest.json").write_text(
        json.dumps({"source": str(SOURCE.relative_to(ROOT)), "crops": manifest, "contact_sheet": str(sheet.relative_to(ROOT))}, indent=2) + "\n"
    )
    return sheet


def main() -> None:
    qa_specs = [
        {"name": "01_header_badges", "box": (120, 40, 3720, 260)},
        {"name": "02_random_panel", "box": (130, 330, 930, 1820)},
        {"name": "03_protocol_panel", "box": (935, 330, 2620, 1820)},
        {"name": "04_inventory_panel", "box": (2625, 330, 3710, 1820)},
        {"name": "05_footer", "box": (130, 1860, 3710, 2065)},
        {"name": "06_full_downscaled", "box": (0, 0, 3840, 2160)},
    ]
    strict_specs = [
        {"name": "01_full_thumbnail", "box": (0, 0, 3840, 2160)},
        {"name": "02_header_title_badges", "box": (120, 40, 3720, 260)},
        {"name": "03_random_top_split", "box": (150, 350, 910, 790)},
        {"name": "04_random_problem_note", "box": (150, 1495, 910, 1800)},
        {"name": "05_protocol_header_rail", "box": (955, 350, 2600, 850)},
        {"name": "06_protocol_cards", "box": (955, 940, 2600, 1285)},
        {"name": "07_protocol_reader_rule", "box": (955, 1650, 2600, 1800)},
        {"name": "08_inventory_stats", "box": (2645, 500, 3690, 705)},
        {"name": "09_inventory_rows", "box": (2645, 675, 3690, 1490)},
        {"name": "10_inventory_boundary", "box": (2645, 1510, 3690, 1800)},
        {"name": "11_footer", "box": (150, 1880, 3690, 2045)},
    ]
    qa_sheet = crop_set(QA_DIR, qa_specs, "mission_holdout_qa_crop_contact_sheet.png")
    strict_sheet = crop_set(STRICT_DIR, strict_specs, "mission_holdout_STRICT_QA_contact_sheet.png")
    print(
        json.dumps(
            {
                "qa_contact_sheet": str(qa_sheet.relative_to(ROOT)),
                "strict_contact_sheet": str(strict_sheet.relative_to(ROOT)),
            },
            indent=2,
        )
    )


if __name__ == "__main__":
    main()
