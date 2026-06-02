#!/usr/bin/env python3
"""Build a methods-bridge shortlist contact sheet.

This is a review artifact, not a final manuscript figure. It decides which
source-proof scenes can enter the slide deck and which should remain appendix
or reviewer backup material.
"""

from __future__ import annotations

import csv
import json
from pathlib import Path

from PIL import Image, ImageDraw, ImageFilter, ImageFont, ImageOps


ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "output" / "biovis_methods_bridge_shortlist_v0_1"
CONTACT = OUT / "contact_sheet"
QA = OUT / "qa"
CREATED = "2026-06-02"

COLORS = {
    "void": "#070A0E",
    "deep": "#0B1117",
    "deep2": "#101923",
    "panel": "#142232",
    "panel2": "#172B3A",
    "ink": "#F4F7F8",
    "soft": "#B7C5D1",
    "muted": "#758695",
    "rule": "#33465A",
    "blue": "#2D6F9F",
    "sky": "#6BAED6",
    "teal": "#1AA090",
    "green": "#178B63",
    "amber": "#B4832F",
    "rose": "#B45A7E",
    "red": "#B23A3A",
}


def rgb(token: str) -> tuple[int, int, int]:
    value = COLORS[token].lstrip("#")
    return tuple(int(value[i : i + 2], 16) for i in (0, 2, 4))


def rgba(token: str, alpha: int) -> tuple[int, int, int, int]:
    return rgb(token) + (alpha,)


def rel(path: Path) -> str:
    return str(path.relative_to(ROOT))


def font(size: int, *, bold: bool = False):
    candidates = [
        "/System/Library/Fonts/Supplemental/Arial Bold.ttf" if bold else "/System/Library/Fonts/Supplemental/Arial.ttf",
        "/System/Library/Fonts/Helvetica.ttc",
    ]
    for candidate in candidates:
        try:
            return ImageFont.truetype(candidate, size=size)
        except OSError:
            continue
    return ImageFont.load_default()


F = {
    "eyebrow": font(22, bold=True),
    "title": font(64, bold=True),
    "subtitle": font(28),
    "section": font(25, bold=True),
    "h": font(25, bold=True),
    "body": font(20),
    "small": font(16),
    "tiny": font(13),
}


CANDIDATES = [
    {
        "id": "organoid_source_to_matrix",
        "title": "Human organoid source-to-matrix",
        "image": ROOT / "output" / "biovis_organoid_audience_matrix_proof_v0_2" / "panels" / "01_dark_organoid_clean_source_to_matrix.png",
        "proposed_use": "Main deck",
        "slot": "Methods bridge: extension readiness",
        "rationale": "Shows official OSDR pages becoming downloaded expression matrices. Good first-view explanation for how organoid data enters v9.",
        "guardrail": "Draft pilot only; not mission-held-out and not a final benchmark result.",
        "status": "keep-main-candidate",
        "tone": "teal",
        "rank": 1,
    },
    {
        "id": "osd120_source_to_task",
        "title": "OSD-120 source-to-task",
        "image": ROOT / "output" / "biovis_osd120_audience_split_proof_v0_2" / "panels" / "01_dark_osd120_clean_source_to_task.png",
        "proposed_use": "Main deck",
        "slot": "Methods bridge: same-study task check",
        "rationale": "Explains how one OSDR study becomes a constrained task check with held-out groups and visible sample balance.",
        "guardrail": "Same-study only; not mission-held-out and not cross-species generalization.",
        "status": "keep-main-candidate",
        "tone": "green",
        "rank": 2,
    },
    {
        "id": "osdr_source_contact_sheet",
        "title": "OSDR source-record contact sheet",
        "image": ROOT / "output" / "biovis_osdr_source_record_proof_panel_v0_1" / "panels" / "01_osdr_source_record_contact_sheet.png",
        "proposed_use": "Appendix",
        "slot": "Reviewer provenance receipt",
        "rationale": "Efficient proof that OSD-863, OSD-871, OSD-120, and OSD-918 official source records were captured.",
        "guardrail": "Use as receipt material; do not turn it into a main-story methods slide.",
        "status": "appendix-reviewer-backup",
        "tone": "amber",
        "rank": 3,
    },
    {
        "id": "extension_source_lanes",
        "title": "Extension source lanes",
        "image": ROOT / "output" / "biovis_osdr_source_record_proof_panel_v0_1" / "panels" / "04_dark_extension_source_lanes.png",
        "proposed_use": "Appendix",
        "slot": "Extension map backup",
        "rationale": "Useful for showing organoid, plant, and single-cell lanes as possible extensions around the v9 provenance layer.",
        "guardrail": "Keep behind the core deck unless extension scope is explicitly discussed.",
        "status": "appendix-reviewer-backup",
        "tone": "amber",
        "rank": 4,
    },
    {
        "id": "organoid_source_record_proof_old",
        "title": "Earlier organoid source proof",
        "image": ROOT / "output" / "biovis_osdr_source_record_proof_panel_v0_1" / "panels" / "02_dark_organoid_source_record_proof.png",
        "proposed_use": "Superseded",
        "slot": "Replace with v0.2 matrix proof",
        "rationale": "The newer organoid panel explains the local matrix step more clearly and uses cleaner audience-facing text.",
        "guardrail": "Do not use in main deck unless the v0.2 panel fails slide-level QA.",
        "status": "superseded",
        "tone": "red",
        "rank": 5,
    },
    {
        "id": "osd120_source_record_proof_old",
        "title": "Earlier OSD-120 source proof",
        "image": ROOT / "output" / "biovis_osdr_source_record_proof_panel_v0_1" / "panels" / "03_dark_osd120_source_record_proof.png",
        "proposed_use": "Superseded",
        "slot": "Replace with v0.2 task proof",
        "rationale": "The newer OSD-120 panel better explains the held-out task and avoids overly internal wording.",
        "guardrail": "Keep only as provenance backup, not as a slide candidate.",
        "status": "superseded",
        "tone": "red",
        "rank": 6,
    },
]


VISIBLE_TEXT = [
    "Methods-bridge shortlist",
    "Main deck",
    "Appendix",
    "Superseded",
    "Evidence layers only",
    "Keep v1-v7 results first",
    "Draft pilot only",
    "Same-study only",
]

BLOCKED_VISIBLE_TERMS = [
    "diagnostic scorer",
    "hidden",
    "within-source",
    "payloads not hashed",
    "single-mission",
]


def ensure() -> None:
    for directory in [OUT, CONTACT, QA]:
        directory.mkdir(parents=True, exist_ok=True)


def draw_wrapped(
    draw: ImageDraw.ImageDraw,
    xy: tuple[int, int],
    text: str,
    font_obj,
    fill,
    *,
    width: int,
    line_height: int,
    max_lines: int,
) -> int:
    words = text.split()
    lines: list[str] = []
    current = ""
    for word in words:
        candidate = f"{current} {word}".strip()
        if draw.textlength(candidate, font=font_obj) <= width:
            current = candidate
        else:
            if current:
                lines.append(current)
            current = word
    if current:
        lines.append(current)
    if len(lines) > max_lines:
        lines = lines[:max_lines]
        lines[-1] = lines[-1].rstrip(".") + "..."
    x, y = xy
    for idx, line in enumerate(lines):
        draw.text((x, y + idx * line_height), line, font=font_obj, fill=fill)
    return y + len(lines) * line_height


def fit_cover(image: Image.Image, size: tuple[int, int]) -> Image.Image:
    return ImageOps.fit(image.convert("RGB"), size, method=Image.Resampling.LANCZOS, centering=(0.5, 0.48))


def shadowed_plate(
    canvas: Image.Image,
    image: Image.Image,
    box: tuple[int, int, int, int],
    *,
    accent: str,
) -> None:
    x, y, w, h = box
    shadow = Image.new("RGBA", canvas.size, (0, 0, 0, 0))
    sd = ImageDraw.Draw(shadow)
    sd.rounded_rectangle((x + 18, y + 20, x + w + 18, y + h + 20), radius=13, fill=(0, 0, 0, 136))
    shadow = shadow.filter(ImageFilter.GaussianBlur(18))
    canvas.alpha_composite(shadow)

    draw = ImageDraw.Draw(canvas)
    draw.rounded_rectangle((x - 4, y - 4, x + w + 4, y + h + 4), radius=15, fill=rgba("panel", 255), outline=rgba(accent, 210), width=3)
    thumb = fit_cover(image, (w, h)).convert("RGBA")
    mask = Image.new("L", (w, h), 0)
    md = ImageDraw.Draw(mask)
    md.rounded_rectangle((0, 0, w, h), radius=12, fill=255)
    canvas.paste(thumb, (x, y), mask)
    draw.rounded_rectangle((x, y, x + w, y + h), radius=12, outline=rgba("ink", 34), width=1)


def draw_background(canvas: Image.Image) -> None:
    draw = ImageDraw.Draw(canvas)
    w, h = canvas.size
    draw.rectangle((0, 0, w, h), fill=rgb("void"))

    top = rgb("void")
    bottom = rgb("deep2")
    for y in range(0, h, 2):
        t = y / max(1, h - 1)
        color = tuple(int(top[i] * (1 - t) + bottom[i] * t) for i in range(3))
        draw.line((0, y, w, y), fill=color + (255,), width=2)

    for x in range(170, w, 170):
        draw.line((x, 0, x, h), fill=rgba("rule", 28), width=1)
    for y in range(160, h, 160):
        draw.line((0, y, w, y), fill=rgba("rule", 24), width=1)

    center = (int(w * 0.74), int(h * 0.27))
    for idx, radius in enumerate([820, 1050, 1280]):
        bbox = (
            center[0] - radius,
            center[1] - int(radius * 0.32),
            center[0] + radius,
            center[1] + int(radius * 0.32),
        )
        draw.arc(bbox, 200, 350, fill=rgba("sky", 50 - idx * 11), width=3)

    for idx in range(42):
        x = int((idx * 263) % w)
        y = int((idx * 419) % h)
        r = 2 + idx % 3
        draw.ellipse((x - r, y - r, x + r, y + r), fill=rgba("soft", 34))


def draw_pill(draw: ImageDraw.ImageDraw, xy: tuple[int, int], text: str, tone: str) -> None:
    x, y = xy
    width = int(draw.textlength(text, font=F["small"])) + 42
    draw.rounded_rectangle((x, y, x + width, y + 34), radius=12, fill=rgba("deep2", 236), outline=rgb(tone), width=2)
    draw.ellipse((x + 14, y + 12, x + 24, y + 22), fill=rgb(tone))
    draw.text((x + 34, y + 9), text, font=F["small"], fill=rgb("ink"))


def draw_candidate(canvas: Image.Image, candidate: dict, box: tuple[int, int, int, int]) -> None:
    draw = ImageDraw.Draw(canvas)
    x, y, w, h = box
    tone = candidate["tone"]

    draw.rounded_rectangle((x, y, x + w, y + h), radius=12, fill=rgba("deep", 222), outline=rgba("rule", 190), width=2)
    draw.rectangle((x, y, x + 8, y + h), fill=rgba(tone, 238))

    image = Image.open(candidate["image"])
    thumb_box = (x + 34, y + 46, 560, 315)
    shadowed_plate(canvas, image, thumb_box, accent=tone)

    tx = x + 640
    draw_pill(draw, (tx, y + 42), candidate["proposed_use"], tone)
    draw.text((tx, y + 92), candidate["title"], font=F["h"], fill=rgb("ink"))
    draw.text((tx, y + 128), candidate["slot"], font=F["small"], fill=rgb("sky" if tone != "red" else "soft"))

    y2 = draw_wrapped(
        draw,
        (tx, y + 172),
        candidate["rationale"],
        F["body"],
        rgb("soft"),
        width=w - 680,
        line_height=28,
        max_lines=3,
    )
    draw_wrapped(
        draw,
        (tx, y2 + 14),
        f"Guardrail: {candidate['guardrail']}",
        F["small"],
        rgb("muted") if tone != "red" else rgb("rose"),
        width=w - 680,
        line_height=23,
        max_lines=2,
    )

    rank = f"{candidate['rank']:02d}"
    draw.text((x + w - 64, y + 32), rank, font=F["section"], fill=rgba(tone, 225))


def render_contact_sheet() -> Path:
    canvas = Image.new("RGBA", (3840, 2160), (0, 0, 0, 255))
    draw_background(canvas)
    draw = ImageDraw.Draw(canvas)

    draw.text((130, 92), "Methods-bridge shortlist", font=F["title"], fill=rgb("ink"))
    draw_wrapped(
        draw,
        (134, 176),
        "Selection board for source-proof scenes. Use them as evidence layers; keep v1-v7 benchmark results as the deck core.",
        F["subtitle"],
        rgb("soft"),
        width=1900,
        line_height=38,
        max_lines=2,
    )
    draw.text((3170, 102), "2026-06-02", font=F["small"], fill=rgb("muted"))
    draw.text((3018, 132), "review artifact, not final slide", font=F["small"], fill=rgb("muted"))

    rule_x = 2470
    draw.rounded_rectangle((rule_x, 98, 3708, 230), radius=16, fill=rgba("panel", 215), outline=rgba("rule", 190), width=2)
    draw.text((rule_x + 30, 122), "Deck rule", font=F["section"], fill=rgb("ink"))
    draw_wrapped(
        draw,
        (rule_x + 30, 158),
        "At most one organoid bridge and one OSD-120 bridge. Appendix absorbs receipts. Core result slides stay first.",
        F["small"],
        rgb("soft"),
        width=1130,
        line_height=23,
        max_lines=2,
    )

    left_x = 130
    right_x = 1995
    row_y = [340, 895, 1450]
    box_w = 1715
    box_h = 456

    section_labels = [
        (left_x, 292, "Main-deck candidates", "teal"),
        (right_x, 292, "Appendix / reviewer backup", "amber"),
        (left_x, 1402, "Superseded for audience use", "red"),
    ]
    for sx, sy, text, tone in section_labels:
        draw.text((sx, sy), text, font=F["section"], fill=rgb("ink"))
        draw.line((sx, sy + 42, sx + 500, sy + 42), fill=rgba(tone, 220), width=3)

    positions = [
        (left_x, row_y[0], box_w, box_h),
        (left_x, row_y[1], box_w, box_h),
        (right_x, row_y[0], box_w, box_h),
        (right_x, row_y[1], box_w, box_h),
        (left_x, row_y[2], box_w, box_h),
        (right_x, row_y[2], box_w, box_h),
    ]
    for candidate, box in zip(CANDIDATES, positions):
        draw_candidate(canvas, candidate, box)

    footer_y = 2038
    draw.line((130, footer_y - 34, 3710, footer_y - 34), fill=rgba("rule", 170), width=1)
    draw.text((130, footer_y), "Next: place only selected bridges into the v1-v9 deck outline; stop open-ended proof generation.", font=F["body"], fill=rgb("soft"))
    draw.text((2735, footer_y), "Claim boundary: provenance/methods only", font=F["body"], fill=rgb("muted"))

    out = CONTACT / "methods_bridge_shortlist_contact_sheet.png"
    canvas.convert("RGB").save(out, quality=95)
    return out


def write_manifest(contact_sheet: Path, grayscale: Path) -> tuple[Path, Path]:
    candidate_records = []
    for candidate in CANDIDATES:
        candidate_records.append(
            {
                "id": candidate["id"],
                "title": candidate["title"],
                "source_image": rel(candidate["image"]),
                "proposed_use": candidate["proposed_use"],
                "deck_slot": candidate["slot"],
                "status": candidate["status"],
                "rationale": candidate["rationale"],
                "guardrail": candidate["guardrail"],
                "rank": candidate["rank"],
            }
        )

    classification = {
        "created": CREATED,
        "verdict": "Use only two main-deck bridge candidates; keep appendix and superseded assets out of the core result spine.",
        "deck_rule": "At most one human organoid bridge and at most one OSD-120/multispecies bridge; v1-v7 benchmark result slides remain the core.",
        "candidates": candidate_records,
    }

    classification_path = OUT / "methods_bridge_slide_classification.json"
    with classification_path.open("w") as handle:
        json.dump(classification, handle, indent=2)
        handle.write("\n")

    csv_path = OUT / "methods_bridge_slide_classification.csv"
    with csv_path.open("w", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=["rank", "id", "title", "proposed_use", "deck_slot", "status", "source_image", "guardrail"],
        )
        writer.writeheader()
        for row in candidate_records:
            writer.writerow({field: row[field] for field in writer.fieldnames})

    visible_blob = "\n".join(
        [
            *VISIBLE_TEXT,
            *[candidate["title"] for candidate in CANDIDATES],
            *[candidate["proposed_use"] for candidate in CANDIDATES],
            *[candidate["slot"] for candidate in CANDIDATES],
            *[candidate["rationale"] for candidate in CANDIDATES],
            *[candidate["guardrail"] for candidate in CANDIDATES],
        ]
    ).lower()
    blocked_results = {term: term not in visible_blob for term in BLOCKED_VISIBLE_TERMS}

    manifest = {
        "created": CREATED,
        "purpose": "Review contact sheet for classifying methods-bridge visual assets before slide-deck assembly.",
        "outputs": {
            "contact_sheet": rel(contact_sheet),
            "grayscale_qa": rel(grayscale),
            "classification_json": rel(classification_path),
            "classification_csv": rel(csv_path),
        },
        "candidate_count": len(CANDIDATES),
        "main_deck_candidate_count": sum(1 for c in CANDIDATES if c["proposed_use"] == "Main deck"),
        "appendix_or_reviewer_candidate_count": sum(1 for c in CANDIDATES if c["proposed_use"] == "Appendix"),
        "superseded_candidate_count": sum(1 for c in CANDIDATES if c["proposed_use"] == "Superseded"),
        "visible_text_policy": {
            "preferred_terms": VISIBLE_TEXT,
            "blocked_internal_terms_absent_from_new_visible_text": blocked_results,
        },
        "claim_boundary": "This contact sheet supports deck placement decisions only; it is not a benchmark result figure.",
    }

    manifest_path = OUT / "methods_bridge_shortlist_manifest.json"
    with manifest_path.open("w") as handle:
        json.dump(manifest, handle, indent=2)
        handle.write("\n")
    return manifest_path, classification_path


def render_grayscale(contact_sheet: Path) -> Path:
    image = Image.open(contact_sheet).convert("L").convert("RGB")
    out = QA / "methods_bridge_shortlist_contact_sheet_grayscale.png"
    image.save(out, quality=95)
    return out


def main() -> None:
    ensure()
    missing = [candidate["image"] for candidate in CANDIDATES if not candidate["image"].exists()]
    if missing:
        raise FileNotFoundError("Missing source images: " + ", ".join(str(path) for path in missing))
    contact_sheet = render_contact_sheet()
    grayscale = render_grayscale(contact_sheet)
    manifest_path, classification_path = write_manifest(contact_sheet, grayscale)
    print(json.dumps({"contact_sheet": rel(contact_sheet), "manifest": rel(manifest_path), "classification": rel(classification_path)}, indent=2))


if __name__ == "__main__":
    main()
