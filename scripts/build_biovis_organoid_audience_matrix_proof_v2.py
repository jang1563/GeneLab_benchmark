#!/usr/bin/env python3
"""Build audience-facing human organoid matrix proof visuals.

v0.2 keeps the real matrix-derived heatmap evidence, but lowers internal audit
language for presentation use.
"""

from __future__ import annotations

import csv
import json
import math
from pathlib import Path

from PIL import Image, ImageDraw, ImageFilter, ImageFont, ImageOps


ROOT = Path(__file__).resolve().parents[1]
MANIFEST_PATH = ROOT / "v9" / "human_organoid" / "task_manifests" / "draft_human_organoid_spaceflight.json"
CAPTURE_DIR = ROOT / "output" / "biovis_osdr_source_record_captures_v0_1"
OUT = ROOT / "output" / "biovis_organoid_audience_matrix_proof_v0_2"
PROOF = OUT / "proof"
PANELS = OUT / "panels"
QA = OUT / "qa"
CREATED = "2026-06-02"

COLORS = {
    "void": "#080B0F",
    "deep": "#0C1218",
    "deep2": "#101922",
    "deep3": "#172231",
    "panel": "#13202B",
    "ink": "#F4F7F8",
    "soft": "#B7C5D1",
    "muted": "#718293",
    "rule": "#33465A",
    "blue": "#2D6F9F",
    "sky": "#6BAED6",
    "teal": "#1AA090",
    "green": "#178B63",
    "amber": "#B4832F",
    "rose": "#B45A7E",
    "red": "#B23A3A",
    "violet": "#6750A4",
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
    "eyebrow": font(20, bold=True),
    "title": font(48, bold=True),
    "subtitle": font(25),
    "h1": font(34, bold=True),
    "h2": font(26, bold=True),
    "body": font(20),
    "small": font(16),
    "tiny": font(13),
    "slide_title": font(66, bold=True),
    "slide_subtitle": font(32),
}


VISIBLE_TEXT = [
    "Human organoid matrix proof",
    "Official OSDR pages",
    "Downloaded expression matrices",
    "Cortical organoids",
    "Dopaminergic organoids",
    "42 samples",
    "Draft pilot",
    "Not mission-held-out",
]


def ensure() -> None:
    for directory in [OUT, PROOF, PANELS, QA]:
        directory.mkdir(parents=True, exist_ok=True)


def load_manifest() -> dict:
    with MANIFEST_PATH.open() as handle:
        return json.load(handle)


def fit_cover_top(image: Image.Image, size: tuple[int, int]) -> Image.Image:
    return ImageOps.fit(image.convert("RGB"), size, method=Image.Resampling.LANCZOS, centering=(0.5, 0.0))


def fit_contain(image: Image.Image, size: tuple[int, int]) -> Image.Image:
    copy = image.convert("RGB")
    copy.thumbnail(size, Image.Resampling.LANCZOS)
    return copy


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


def card(draw: ImageDraw.ImageDraw, xywh: tuple[int, int, int, int], *, alpha: int = 226, outline: str = "rule", radius: int = 24) -> None:
    x, y, w, h = xywh
    draw.rounded_rectangle((x, y, x + w, y + h), radius=radius, fill=rgba("deep2", alpha), outline=rgba(outline, 230), width=2)


def pill(draw: ImageDraw.ImageDraw, xy: tuple[int, int], text: str, tone: str, *, fill_token: str = "deep2") -> None:
    x, y = xy
    width = max(112, int(draw.textlength(text, font=F["small"])) + 42)
    draw.rounded_rectangle((x, y, x + width, y + 34), radius=14, fill=rgba(fill_token, 236), outline=rgb(tone), width=2)
    draw.ellipse((x + 14, y + 12, x + 24, y + 22), fill=rgb(tone))
    draw.text((x + 34, y + 9), text, font=F["small"], fill=rgb("ink"))


def draw_background(canvas: Image.Image, *, accent: str = "teal") -> None:
    draw = ImageDraw.Draw(canvas)
    w, h = canvas.size
    draw.rectangle((0, 0, w, h), fill=rgb("void"))
    for y in range(0, h, 2):
        t = y / max(1, h - 1)
        color = tuple(int(rgb("void")[i] * (1 - t) + rgb("deep2")[i] * t) for i in range(3))
        draw.line((0, y, w, y), fill=color + (255,), width=2)
    for x in range(160, w, 160):
        draw.line((x, 0, x, h), fill=rgba("rule", 36), width=1)
    for y in range(140, h, 140):
        draw.line((0, y, w, y), fill=rgba("rule", 30), width=1)
    center = (int(w * 0.68), int(h * 0.34))
    for idx, radius in enumerate([720, 940, 1160]):
        bbox = (center[0] - radius, center[1] - int(radius * 0.36), center[0] + radius, center[1] + int(radius * 0.36))
        draw.arc(bbox, 198, 348, fill=rgba(accent, 58 - idx * 12), width=3)
    for idx in range(38):
        x = int((idx * 251) % w)
        y = int((idx * 431) % h)
        r = 2 + idx % 3
        draw.ellipse((x - r, y - r, x + r, y + r), fill=rgba("soft", 36))


def draw_cell_motif(draw: ImageDraw.ImageDraw, xy: tuple[int, int], *, scale: float = 1.0, tone: str = "teal") -> None:
    x, y = xy
    nodes = [(0, 0, 32), (54, -10, 24), (94, 26, 30), (34, 54, 22), (86, 82, 18), (132, 66, 20)]
    for i, (dx, dy, r) in enumerate(nodes):
        cx = x + int(dx * scale)
        cy = y + int(dy * scale)
        rr = int(r * scale)
        draw.ellipse((cx - rr, cy - rr, cx + rr, cy + rr), outline=rgba(tone, 125), width=max(2, int(3 * scale)))
        draw.ellipse((cx - rr // 3, cy - rr // 3, cx + rr // 3, cy + rr // 3), fill=rgba(tone, 72))
        if i > 0:
            px = x + int(nodes[i - 1][0] * scale)
            py = y + int(nodes[i - 1][1] * scale)
            draw.line((px, py, cx, cy), fill=rgba(tone, 70), width=max(1, int(2 * scale)))


def compact_sha(value: str) -> str:
    return f"{value[:10]}...{value[-8:]}" if value else "NA"


def read_matrix_crop(path: Path, *, rows: int = 42) -> tuple[list[str], list[str], list[list[float]]]:
    with path.open(newline="") as handle:
        reader = csv.reader(handle)
        header = next(reader)
        samples = header[1:]
        genes: list[str] = []
        values: list[list[float]] = []
        for row in reader:
            if len(row) < 2:
                continue
            try:
                numeric = [math.log2(float(value) + 1.0) for value in row[1:]]
            except ValueError:
                continue
            genes.append(row[0].strip('"'))
            values.append(numeric)
            if len(values) >= rows:
                break
    return genes, samples, values


def value_color(value: float, low: float, high: float) -> tuple[int, int, int]:
    if high <= low:
        t = 0.5
    else:
        t = max(0.0, min(1.0, (value - low) / (high - low)))
    stops = [
        (0.0, (20, 35, 51)),
        (0.32, rgb("blue")),
        (0.62, rgb("teal")),
        (0.82, rgb("amber")),
        (1.0, rgb("rose")),
    ]
    for idx in range(len(stops) - 1):
        t0, c0 = stops[idx]
        t1, c1 = stops[idx + 1]
        if t0 <= t <= t1:
            local = (t - t0) / (t1 - t0)
            return tuple(int(c0[channel] + (c1[channel] - c0[channel]) * local) for channel in range(3))
    return stops[-1][1]


def draw_heatmap(
    draw: ImageDraw.ImageDraw,
    xy: tuple[int, int],
    values: list[list[float]],
    *,
    cell_w: int,
    cell_h: int,
    gap: int = 2,
) -> None:
    flat = [value for row in values for value in row]
    flat_sorted = sorted(flat)
    lo = flat_sorted[int(len(flat_sorted) * 0.04)]
    hi = flat_sorted[int(len(flat_sorted) * 0.96)]
    x0, y0 = xy
    for r, row in enumerate(values):
        row_mean = sum(row) / len(row)
        centered = [value - row_mean for value in row]
        for col, value in enumerate(centered):
            color = value_color(value, lo - row_mean, hi - row_mean)
            x = x0 + col * (cell_w + gap)
            y = y0 + r * (cell_h + gap)
            draw.rounded_rectangle((x, y, x + cell_w, y + cell_h), radius=2, fill=color)


def organoid_label(source_id: str) -> tuple[str, str, str]:
    if source_id == "OSD-863":
        return "Cortical organoids", "GLDS-716", "sky"
    return "Dopaminergic organoids", "GLDS-720", "teal"


def render_clean_matrix_proof(manifest: dict) -> Path:
    split = manifest["split"]
    matrices = split["expression_matrix_sources"]
    sources = ["OSD-863", "OSD-871"]
    canvas = Image.new("RGBA", (2200, 1250), rgba("deep", 255))
    draw = ImageDraw.Draw(canvas)
    draw_cell_motif(draw, (1900, 96), scale=1.05, tone="teal")
    draw.text((70, 54), "Human organoid matrix proof", font=F["title"], fill=rgb("ink"))
    draw_wrapped(
        draw,
        (70, 110),
        "Two official OSDR organoid studies are represented by downloaded RNA-seq expression matrices.",
        F["subtitle"],
        rgb("soft"),
        width=1580,
        line_height=32,
        max_lines=2,
    )
    draw.line((70, 172, 2130, 172), fill=rgba("rule", 255), width=2)
    pill(draw, (70, 202), "OSD-863 / OSD-871", "teal")
    pill(draw, (268, 202), "42 samples", "sky")
    pill(draw, (416, 202), "GSE259421", "violet")
    pill(draw, (566, 202), "draft pilot", "amber")
    pill(draw, (710, 202), "not mission-held-out", "red")

    positions = [(70, 292), (1130, 292)]
    for source_id, (x, y) in zip(sources, positions):
        entry = matrices[source_id]
        label, glds_id, tone = organoid_label(source_id)
        _, samples, values = read_matrix_crop(ROOT / entry["local_matrix_path"])
        card(draw, (x, y, 1000, 770), alpha=216, outline="rule")
        draw.text((x + 36, y + 38), source_id, font=F["h1"], fill=rgb("ink"))
        draw.text((x + 36, y + 78), label, font=F["body"], fill=rgb("soft"))
        pill(draw, (x + 36, y + 122), glds_id, tone)
        pill(draw, (x + 170, y + 122), "downloaded matrix", "green")
        heat_x, heat_y = x + 46, y + 192
        cell_w, cell_h, gap = 32, 8, 2
        draw_heatmap(draw, (heat_x, heat_y), values, cell_w=cell_w, cell_h=cell_h, gap=gap)
        heat_w = len(samples) * (cell_w + gap) - gap
        heat_h = len(values) * (cell_h + gap) - gap
        draw.rectangle((heat_x, heat_y, heat_x + heat_w, heat_y + heat_h), outline=rgba("soft", 120), width=1)
        draw.text((x + 48, y + 632), f"{entry['n_feature_rows']:,} genes x {entry['n_sample_columns']} samples", font=F["h2"], fill=rgb("ink"))
        draw.text((x + 48, y + 674), f"sample columns aligned; sha256 {compact_sha(entry['sha256'])}", font=F["small"], fill=rgb("soft"))
        draw.text((x + 48, y + 706), entry["matrix_file"], font=F["tiny"], fill=rgb("muted"))

    card(draw, (70, 1084, 2060, 96), alpha=190, outline="rule", radius=18)
    boundary_label = "Claim boundary"
    draw.text((104, 1118), boundary_label, font=F["h2"], fill=rgb("ink"))
    body_x = 104 + int(draw.textlength(boundary_label, font=F["h2"])) + 34
    draw.text((body_x, 1122), "This proves local matrix presence and sample alignment. It is not a final benchmark result or mission-held-out claim.", font=F["body"], fill=rgb("soft"))
    out = PROOF / "human_organoid_downloaded_matrix_proof.png"
    canvas.convert("RGB").save(out)
    return out


def source_crop(source_id: str) -> Image.Image:
    path = CAPTURE_DIR / f"{source_id.lower().replace('-', '_')}_viewport.png"
    source = Image.open(path).convert("RGB")
    return source.crop((292, 145, 1580, 1120))


def shadowed_surface(
    canvas: Image.Image,
    image: Image.Image,
    xy: tuple[int, int],
    *,
    width: int,
    height: int,
    radius: int,
    outline: tuple[int, int, int, int],
    contain: bool = False,
) -> None:
    x, y = xy
    surface = (fit_contain(image, (width - 70, height - 70)) if contain else fit_cover_top(image, (width, height))).convert("RGBA")
    shadow = Image.new("RGBA", canvas.size, (0, 0, 0, 0))
    sd = ImageDraw.Draw(shadow)
    sd.rounded_rectangle((x + 18, y + 26, x + width + 18, y + height + 26), radius=radius, fill=(0, 0, 0, 155))
    shadow = shadow.filter(ImageFilter.GaussianBlur(34))
    canvas.alpha_composite(shadow)
    draw = ImageDraw.Draw(canvas)
    draw.rounded_rectangle((x, y, x + width, y + height), radius=radius, fill=rgba("deep2", 238), outline=outline, width=3)
    if contain:
        px = x + (width - surface.width) // 2
        py = y + (height - surface.height) // 2
        canvas.alpha_composite(surface, (px, py))
    else:
        mask = Image.new("L", (width, height), 0)
        md = ImageDraw.Draw(mask)
        md.rounded_rectangle((0, 0, width, height), radius=radius, fill=255)
        canvas.paste(surface, (x, y), mask)
        draw.rounded_rectangle((x, y, x + width, y + height), radius=radius, outline=outline, width=3)


def arrow(draw: ImageDraw.ImageDraw, start: tuple[int, int], end: tuple[int, int], *, tone: str = "teal") -> None:
    x1, y1 = start
    x2, y2 = end
    mid = ((x1 + x2) // 2, min(y1, y2) - 80)
    points = []
    for i in range(44):
        t = i / 43
        x = int((1 - t) ** 2 * x1 + 2 * (1 - t) * t * mid[0] + t**2 * x2)
        y = int((1 - t) ** 2 * y1 + 2 * (1 - t) * t * mid[1] + t**2 * y2)
        points.append((x, y))
    draw.line(points, fill=rgba(tone, 205), width=5)
    angle = math.atan2(y2 - points[-3][1], x2 - points[-3][0])
    size = 23
    p1 = (x2, y2)
    p2 = (int(x2 - size * math.cos(angle - 0.45)), int(y2 - size * math.sin(angle - 0.45)))
    p3 = (int(x2 - size * math.cos(angle + 0.45)), int(y2 - size * math.sin(angle + 0.45)))
    draw.polygon([p1, p2, p3], fill=rgba(tone, 220))


def caption(draw: ImageDraw.ImageDraw, xywh: tuple[int, int, int, int], title: str, body: str, *, tone: str) -> None:
    x, y, w, h = xywh
    draw.rounded_rectangle((x, y, x + w, y + h), radius=24, fill=rgba("deep2", 228), outline=rgba(tone, 160), width=2)
    draw.text((x + 30, y + 28), title, font=F["h2"], fill=rgb("ink"))
    draw_wrapped(draw, (x + 30, y + 76), body, F["small"], rgb("soft"), width=w - 60, line_height=27, max_lines=4)


def render_source_to_matrix_panel(proof_path: Path) -> Path:
    canvas = Image.new("RGBA", (3840, 2160), rgba("void", 255))
    draw_background(canvas, accent="teal")
    draw = ImageDraw.Draw(canvas)
    draw_cell_motif(draw, (3220, 230), scale=2.2, tone="teal")
    draw.text((120, 92), "Human organoid source to matrix proof", font=F["slide_title"], fill=rgb("ink"))
    draw_wrapped(
        draw,
        (122, 176),
        "Official OSDR pages on the left; downloaded expression matrices on the right.",
        F["slide_subtitle"],
        rgb("soft"),
        width=2150,
        line_height=42,
        max_lines=2,
    )
    pill(draw, (122, 260), "official source pages", "teal")
    pill(draw, (408, 260), "downloaded matrices", "green")
    pill(draw, (688, 260), "42 samples", "sky")
    pill(draw, (848, 260), "draft pilot", "amber")
    pill(draw, (1002, 260), "claim boundary visible", "red")

    source_863 = source_crop("OSD-863")
    source_871 = source_crop("OSD-871")
    proof = Image.open(proof_path).convert("RGB")
    draw.text((120, 350), "Official OSDR source pages", font=F["h2"], fill=rgb("ink"))
    draw.text((1780, 350), "Audience-facing local proof", font=F["h2"], fill=rgb("ink"))
    shadowed_surface(canvas, source_863, (120, 405), width=1480, height=500, radius=34, outline=rgba("teal", 230))
    shadowed_surface(canvas, source_871, (245, 925), width=1480, height=500, radius=34, outline=rgba("sky", 230))
    shadowed_surface(canvas, proof, (1780, 405), width=1900, height=1080, radius=38, outline=rgba("teal", 230), contain=True)
    arrow(draw, (1600, 650), (1780, 795), tone="teal")
    arrow(draw, (1725, 1175), (1780, 1085), tone="sky")

    caption(
        draw,
        (120, 1518, 1480, 244),
        "Viewer path",
        "Start from the official organoid study records, then move to the local expression matrices used by the draft task.",
        tone="teal",
    )
    caption(
        draw,
        (1780, 1518, 1420, 244),
        "What this proves",
        "The downloaded matrices exist locally and their sample columns are aligned. This does not validate final benchmark performance.",
        tone="amber",
    )
    draw.rounded_rectangle((120, 1948, 3720, 2040), radius=22, fill=rgba("deep2", 230), outline=rgba("rule", 220), width=2)
    draw.text((158, 1981), "Claim boundary", font=F["h2"], fill=rgb("ink"))
    draw_wrapped(
        draw,
        (410, 1984),
        "Use this as a methods bridge for the organoid extension. Keep source URLs and draft status in the deck footer.",
        F["small"],
        rgb("soft"),
        width=3140,
        line_height=27,
        max_lines=2,
    )
    out = PANELS / "01_dark_organoid_clean_source_to_matrix.png"
    canvas.convert("RGB").save(out)
    return out


def grayscale_exports(paths: list[Path]) -> list[Path]:
    outs = []
    for path in paths:
        out = QA / f"{path.stem}_grayscale.png"
        Image.open(path).convert("L").save(out)
        outs.append(out)
    return outs


def luminance(color: tuple[int, int, int]) -> float:
    values = [v / 255 for v in color]
    linear = [v / 12.92 if v <= 0.03928 else ((v + 0.055) / 1.055) ** 2.4 for v in values]
    return 0.2126 * linear[0] + 0.7152 * linear[1] + 0.0722 * linear[2]


def contrast(c1: tuple[int, int, int], c2: tuple[int, int, int]) -> float:
    l1 = luminance(c1)
    l2 = luminance(c2)
    return (max(l1, l2) + 0.05) / (min(l1, l2) + 0.05)


def write_metadata(proof_path: Path, panel_path: Path, grayscale: list[Path], manifest: dict) -> None:
    split = manifest["split"]
    blocked_terms = ["payloads not hashed", "diagnostic scorer", "single-mission"]
    visible = " ".join(VISIBLE_TEXT).lower()
    data = {
        "created": CREATED,
        "source_manifest": rel(MANIFEST_PATH),
        "source_records": {
            "OSD-863": "https://osdr.nasa.gov/bio/repo/data/studies/OSD-863",
            "OSD-871": "https://osdr.nasa.gov/bio/repo/data/studies/OSD-871",
        },
        "outputs": {
            "proof": rel(proof_path),
            "panel": rel(panel_path),
            "grayscale_qa": [rel(path) for path in grayscale],
        },
        "sample_counts": {
            "n_samples": split["n_samples"],
            "label_distribution": split["label_distribution"],
            "organoid_type_distribution": split["organoid_type_distribution"],
            "microglia_condition_distribution": split["microglia_condition_distribution"],
            "source_sample_rows": split["source_sample_rows"],
        },
        "matrix_sources": {
            source_id: {
                "matrix_file": entry["matrix_file"],
                "n_feature_rows": entry["n_feature_rows"],
                "n_sample_columns": entry["n_sample_columns"],
                "sha256": entry["sha256"],
            }
            for source_id, entry in split["expression_matrix_sources"].items()
        },
        "visible_text_policy": {
            "preferred_terms": VISIBLE_TEXT,
            "blocked_internal_terms_absent_from_new_visible_text": {term: term.lower() not in visible for term in blocked_terms},
        },
        "claim_boundary": "Downloaded matrix and sample-alignment proof only; not final benchmark performance or mission-held-out evidence.",
    }
    (OUT / "organoid_audience_matrix_proof_manifest.json").write_text(json.dumps(data, indent=2) + "\n")
    qa = {
        "created": CREATED,
        "contrast_checks": {
            "ink_on_deep": round(contrast(rgb("ink"), rgb("deep")), 2),
            "soft_on_deep": round(contrast(rgb("soft"), rgb("deep")), 2),
            "ink_on_void": round(contrast(rgb("ink"), rgb("void")), 2),
        },
        "visual_checks": [
            "Audience-facing proof foregrounds source pages, downloaded matrices, dimensions, and claim boundary.",
            "Heatmaps are derived from real normalized count matrices.",
            "Internal audit language is moved into metadata rather than visible labels.",
            "Color is reinforced by text labels and grayscale QA outputs.",
        ],
        "review_required": [
            "Final deck should use editable text for callouts, with this PNG as the evidence surface.",
            "Do not use as final result or mission-held-out evidence.",
        ],
    }
    (OUT / "organoid_audience_matrix_proof_qa.json").write_text(json.dumps(qa, indent=2) + "\n")


def main() -> None:
    ensure()
    manifest = load_manifest()
    proof_path = render_clean_matrix_proof(manifest)
    panel_path = render_source_to_matrix_panel(proof_path)
    grayscale = grayscale_exports([proof_path, panel_path])
    write_metadata(proof_path, panel_path, grayscale, manifest)
    print(
        json.dumps(
            {
                "output": rel(OUT),
                "proof": rel(proof_path),
                "panel": rel(panel_path),
                "grayscale_qa": [rel(path) for path in grayscale],
            },
            indent=2,
        )
    )


if __name__ == "__main__":
    main()
