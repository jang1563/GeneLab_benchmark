#!/usr/bin/env python3
"""Build section crops for the first-33 walkthrough contact-sheet QA."""

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
ASSEMBLY_DIR = WORKSPACE / "assets" / "first33_walkthrough_assembly"
SOURCE = ASSEMBLY_DIR / "slides_1_33_walkthrough_contact_sheet.png"
GRAY = ASSEMBLY_DIR / "slides_1_33_walkthrough_contact_sheet_grayscale.png"
QA_DIR = ASSEMBLY_DIR / "section_qa"
QA_DIR.mkdir(parents=True, exist_ok=True)


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


F_TITLE = font(38, bold=True)


def fit(image: Image.Image, size: tuple[int, int]) -> Image.Image:
    image = image.convert("RGB").copy()
    resample = getattr(getattr(Image, "Resampling", Image), "LANCZOS")
    image.thumbnail(size, resample)
    canvas = Image.new("RGB", size, "#0C111A")
    canvas.paste(image, ((size[0] - image.width) // 2, (size[1] - image.height) // 2))
    return canvas


def write_review_sheet(crops: list[dict[str, object]]) -> Path:
    width = 1800
    y = 62
    blocks = []
    for crop in crops:
        image = Image.open(QA_DIR / str(crop["file"]))
        fitted = fit(image, (1650, int(crop["target_h"])))
        blocks.append((str(crop["label"]), fitted, y))
        y += fitted.height + 108

    canvas = Image.new("RGB", (width, y + 20), "#0C111A")
    draw = ImageDraw.Draw(canvas)
    for label, image, top in blocks:
        draw.text((75, top - 40), label, font=F_TITLE, fill="#F3F7FC")
        draw.rounded_rectangle((70, top - 6, 1730, top + image.height + 6), radius=20, outline="#2A394D", width=2)
        canvas.paste(image, (75, top))
    out = QA_DIR / "slides_1_33_walkthrough_section_QA_sheet.png"
    canvas.save(out, quality=95)
    return out


def main() -> None:
    source = Image.open(SOURCE).convert("RGB")
    gray = Image.open(GRAY).convert("RGB")
    specs = [
        {"name": "01_open_method_setup", "label": "Slides 1-12: opening and method setup", "box": (0, 0, 2022, 2475), "target_h": 2020},
        {"name": "02_core_result", "label": "Slides 13-21: core result chain", "box": (0, 2230, 2022, 4000), "target_h": 1450},
        {"name": "03_models_robustness", "label": "Slides 22-33: models, robustness, takeaway", "box": (0, 3740, 2022, 6040), "target_h": 1880},
        {"name": "04_grayscale_full", "label": "Grayscale full sequence", "box": (0, 0, 2022, 6040), "target_h": 2800, "gray": True},
    ]
    manifest = []
    for spec in specs:
        image = gray if spec.get("gray") else source
        out = QA_DIR / f"first33_{spec['name']}.png"
        image.crop(tuple(int(v) for v in spec["box"])).save(out, quality=95)
        manifest.append({"name": spec["name"], "label": spec["label"], "box": spec["box"], "path": str(out.relative_to(ROOT))})
    sheet = write_review_sheet([{**spec, "file": f"first33_{spec['name']}.png"} for spec in specs])
    manifest_path = QA_DIR / "slides_1_33_walkthrough_section_QA_manifest.json"
    manifest_path.write_text(
        json.dumps({"source": str(SOURCE.relative_to(ROOT)), "review_sheet": str(sheet.relative_to(ROOT)), "crops": manifest}, indent=2) + "\n",
        encoding="utf-8",
    )
    print(json.dumps({"review_sheet": str(sheet.relative_to(ROOT)), "manifest": str(manifest_path.relative_to(ROOT))}, indent=2))


if __name__ == "__main__":
    main()
