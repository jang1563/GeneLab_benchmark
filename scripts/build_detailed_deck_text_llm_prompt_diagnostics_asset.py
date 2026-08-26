#!/usr/bin/env python3
"""Build the detailed-deck text-LLM prompt diagnostic asset."""

from __future__ import annotations

import glob
import json
from pathlib import Path

from PIL import Image, ImageDraw, ImageFont


ROOT = Path(__file__).resolve().parent.parent
PROMPT_SPEC = ROOT / "docs" / "text_llm_format.md"
PROMPT_SUMMARY = ROOT / "v2" / "processed" / "llm_prompts" / "generation_summary.json"
PARSE_GLOB = str(ROOT / "v2" / "evaluation" / "parse_metrics_*_A*.json")
EVAL_GLOB = str(ROOT / "evaluation" / "submission_*_zeroshot_A*_eval.json")
OUT_DIR = (
    ROOT
    / "outputs"
    / "019e8bd1-3b42-7821-9e1c-beebfa8e2ece"
    / "presentations"
    / "spacebiobench-detailed-deck"
    / "assets"
    / "text_llm_prompt_diagnostics"
)
OUT_DIR.mkdir(parents=True, exist_ok=True)

W, H = 3840, 2160
OUT_PATH = OUT_DIR / "text_llm_prompt_diagnostics_premium.png"
GRAY_PATH = OUT_DIR / "text_llm_prompt_diagnostics_premium_grayscale.png"
MANIFEST_PATH = OUT_DIR / "text_llm_prompt_diagnostics_manifest.json"

COLORS = {
    "bg": "#0C111A",
    "panel": "#101823",
    "panel2": "#151F2D",
    "panel3": "#121B28",
    "grid": "#2A3546",
    "axis": "#98A7BA",
    "text": "#F3F7FC",
    "muted": "#A8B4C4",
    "dim": "#667487",
    "blue": "#66A6E8",
    "sky": "#73A7FF",
    "teal": "#5FD3C4",
    "green": "#8BD17C",
    "amber": "#F4C26B",
    "rose": "#E17882",
    "violet": "#B39DFF",
    "white": "#FFFFFF",
}


def hex_to_rgb(value: str) -> tuple[int, int, int]:
    value = value.lstrip("#")
    return tuple(int(value[i : i + 2], 16) for i in (0, 2, 4))


def rgba(value: str, alpha: int) -> tuple[int, int, int, int]:
    return (*hex_to_rgb(value), alpha)


def load_font(size: int, bold: bool = False) -> ImageFont.FreeTypeFont:
    candidates = [
        "/System/Library/Fonts/Supplemental/Arial Bold.ttf" if bold else "/System/Library/Fonts/Supplemental/Arial.ttf",
        "/System/Library/Fonts/Supplemental/Helvetica Bold.ttf" if bold else "/System/Library/Fonts/Supplemental/Helvetica.ttf",
        "/Library/Fonts/Arial Bold.ttf" if bold else "/Library/Fonts/Arial.ttf",
    ]
    for candidate in candidates:
        try:
            return ImageFont.truetype(candidate, size)
        except OSError:
            continue
    return ImageFont.load_default()


F = {
    "kicker": load_font(34, True),
    "title": load_font(76, True),
    "subtitle": load_font(35, False),
    "h2": load_font(42, True),
    "h3": load_font(32, True),
    "body": load_font(29, False),
    "body_bold": load_font(29, True),
    "small": load_font(25, False),
    "small_bold": load_font(25, True),
    "tiny": load_font(21, False),
    "tiny_bold": load_font(21, True),
    "micro": load_font(18, False),
    "micro_bold": load_font(18, True),
    "mono": load_font(24, False),
    "mono_bold": load_font(24, True),
    "stat": load_font(66, True),
}


def rounded(
    draw: ImageDraw.ImageDraw,
    box: tuple[int, int, int, int],
    radius: int,
    fill: str | tuple[int, int, int, int],
    outline: str | tuple[int, int, int, int] | None = None,
    width: int = 1,
) -> None:
    draw.rounded_rectangle(box, radius=radius, fill=fill, outline=outline, width=width)


def text(
    draw: ImageDraw.ImageDraw,
    xy: tuple[int | float, int | float],
    value: str,
    font: ImageFont.ImageFont,
    fill: str | tuple[int, int, int] | tuple[int, int, int, int] = COLORS["text"],
    anchor: str | None = None,
) -> None:
    draw.text(xy, value, font=font, fill=fill, anchor=anchor)


def wrap_lines(draw: ImageDraw.ImageDraw, label: str, font: ImageFont.ImageFont, max_width: int) -> list[str]:
    words = label.split()
    lines: list[str] = []
    current: list[str] = []
    for word in words:
        trial = " ".join(current + [word])
        if draw.textlength(trial, font=font) <= max_width:
            current.append(word)
        else:
            if current:
                lines.append(" ".join(current))
            current = [word]
    if current:
        lines.append(" ".join(current))
    return lines


def paragraph(
    draw: ImageDraw.ImageDraw,
    xy: tuple[int, int],
    body: str,
    font: ImageFont.ImageFont,
    max_width: int,
    fill: str | tuple[int, int, int, int],
    leading: int = 7,
) -> int:
    x, y = xy
    for line in wrap_lines(draw, body, font, max_width):
        text(draw, (x, y), line, font, fill)
        y += font.size + leading
    return y


def arrow(draw: ImageDraw.ImageDraw, start: tuple[int, int], end: tuple[int, int], color: str, width: int = 4) -> None:
    x0, y0 = start
    x1, y1 = end
    draw.line((x0, y0, x1, y1), fill=rgba(color, 180), width=width)
    if abs(x1 - x0) >= abs(y1 - y0):
        direction = 1 if x1 >= x0 else -1
        points = [(x1, y1), (x1 - direction * 18, y1 - 11), (x1 - direction * 18, y1 + 11)]
    else:
        direction = 1 if y1 >= y0 else -1
        points = [(x1, y1), (x1 - 11, y1 - direction * 18), (x1 + 11, y1 - direction * 18)]
    draw.polygon(points, fill=rgba(color, 210))


def draw_background(canvas: Image.Image) -> None:
    overlay = Image.new("RGBA", (W, H), (0, 0, 0, 0))
    draw = ImageDraw.Draw(overlay)
    for x in range(0, W, 160):
        draw.line((x, 0, x, H), fill=rgba("#1E2A3A", 28), width=1)
    for y in range(0, H, 160):
        draw.line((0, y, W, y), fill=rgba("#1E2A3A", 22), width=1)
    draw.rectangle((0, 0, W, 310), fill=rgba("#121B28", 170))
    draw.rectangle((0, 1855, W, H), fill=rgba("#090E15", 145))
    canvas.alpha_composite(overlay)


def provider_from_parse_name(provider: str) -> str:
    aliases = {
        "deepseek-chat": "DeepSeek-V3",
        "gemini-2.5-flash": "Gemini-2.5-Flash",
        "groq-llama-3.3-70b": "Llama-3.3-70B",
    }
    return aliases.get(provider, provider)


def load_data() -> dict[str, object]:
    prompt = json.loads(PROMPT_SUMMARY.read_text())
    parse: dict[str, dict[str, int]] = {}
    for path in glob.glob(PARSE_GLOB):
        data = json.loads(Path(path).read_text())
        provider = provider_from_parse_name(str(data["provider"]))
        summary = data["parse_summary"]
        if provider not in parse:
            parse[provider] = {"total": 0, "parsed": 0, "failed": 0}
        parse[provider]["total"] += int(summary["total"])
        parse[provider]["parsed"] += int(summary["parsed"])
        parse[provider]["failed"] += int(summary["failed"])

    task_means: dict[str, list[float]] = {}
    pooled: dict[str, list[float]] = {}
    folds: dict[str, int] = {}
    for path in glob.glob(EVAL_GLOB):
        eval_data = json.loads(Path(path).read_text())
        stem = Path(path).stem
        provider = stem[len("submission_") :].split("_zeroshot_")[0]
        task_means.setdefault(provider, []).append(float(eval_data["summary"]["mean_auroc"]))
        pooled.setdefault(provider, []).append(float(eval_data["summary"]["overall_auroc_pooled"]))
        folds[provider] = folds.get(provider, 0) + int(eval_data["n_folds"])

    providers: list[dict[str, object]] = []
    for provider in ["DeepSeek-V3", "Llama-3.3-70B", "Gemini-2.5-Flash"]:
        p = parse[provider]
        providers.append(
            {
                "provider": provider,
                "task_mean_auroc": sum(task_means[provider]) / len(task_means[provider]),
                "pooled_auroc": sum(pooled[provider]) / len(pooled[provider]),
                "parse_rate": p["parsed"] / p["total"],
                "parsed": p["parsed"],
                "total": p["total"],
                "failed": p["failed"],
                "tasks": len(task_means[provider]),
                "folds": folds[provider],
            }
        )

    return {
        "generated_at": prompt["generated_at"],
        "n_genes": int(prompt["n_genes_per_prompt"]),
        "total_prompts": int(prompt["total_prompts"]),
        "tasks": prompt["tasks"],
        "n_tasks": len(prompt["tasks"]),
        "providers": providers,
        "n_eval_files": len(glob.glob(EVAL_GLOB)),
    }


def stat_badge(draw: ImageDraw.ImageDraw, x: int, y: int, w: int, label: str, value: str, accent: str) -> None:
    rounded(draw, (x, y, x + w, y + 104), 24, "#121B28", "#2A394D", 2)
    text(draw, (x + 26, y + 23), label, F["tiny_bold"], accent)
    text(draw, (x + 26, y + 57), value, F["small_bold"], COLORS["text"])


def draw_header(draw: ImageDraw.ImageDraw, data: dict[str, object]) -> None:
    text(draw, (136, 76), "MODELS / TEXT LLM DIAGNOSTIC", F["kicker"], COLORS["sky"])
    text(draw, (136, 122), "Text LLM checks are prompt diagnostics", F["title"], COLORS["text"])
    text(
        draw,
        (138, 222),
        "The model sees ranked gene text; a parser converts A/B plus confidence into a flight-probability score.",
        F["subtitle"],
        COLORS["muted"],
    )
    badges = [
        ("PROMPTS", f"{data['total_prompts']} samples", 250, COLORS["teal"]),
        ("GENES", f"top-{data['n_genes']}", 178, COLORS["violet"]),
        ("PROVIDERS", "3 LLMs", 208, COLORS["sky"]),
        ("EVALS", f"{data['n_eval_files']} files", 196, COLORS["amber"]),
        ("OUTPUT", "A/B + confidence", 282, COLORS["green"]),
    ]
    bx = 2468
    for label, value, width, accent in badges:
        stat_badge(draw, bx, 72, width, label, value, accent)
        bx += width + 18


def draw_prompt_section(
    draw: ImageDraw.ImageDraw,
    x: int,
    y: int,
    w: int,
    h: int,
    label: str,
    title: str,
    body: list[str],
    accent: str,
) -> int:
    rounded(draw, (x, y, x + w, y + h), 24, COLORS["panel3"], "#2A394D", 2)
    rounded(draw, (x + 24, y + 24, x + 160, y + 62), 17, rgba(accent, 42), rgba(accent, 150), 1)
    text(draw, (x + 92, y + 44), label, F["micro_bold"], accent, "mm")
    text(draw, (x + 188, y + 28), title, F["small_bold"], COLORS["text"])
    yy = y + 78
    for line in body:
        if line == "":
            yy += 18
            continue
        font = F["mono_bold"] if line.startswith(">") else F["mono"]
        fill = COLORS["text"] if line.startswith(">") else COLORS["muted"]
        clean = line[1:].strip() if line.startswith(">") else line
        yy = paragraph(draw, (x + 30, yy), clean, font, w - 60, fill, leading=6)
    return y + h


def draw_prompt_panel(draw: ImageDraw.ImageDraw, box: tuple[int, int, int, int], data: dict[str, object]) -> None:
    x0, y0, x1, y1 = box
    rounded(draw, box, 34, COLORS["panel"], "#29374A", 2)
    text(draw, (x0 + 48, y0 + 46), "A. What the LLM sees", F["h2"], COLORS["text"])
    paragraph(
        draw,
        (x0 + 48, y0 + 98),
        "Each test sample becomes a compact text prompt built from train-selected genes.",
        F["small"],
        x1 - x0 - 96,
        COLORS["muted"],
        leading=6,
    )

    card_x = x0 + 54
    card_y = y0 + 190
    card_w = x1 - x0 - 108
    rounded(draw, (card_x, card_y, card_x + card_w, y1 - 282), 30, "#0E1621", "#2A394D", 2)
    draw_prompt_section(
        draw,
        card_x + 34,
        card_y + 36,
        card_w - 68,
        182,
        "SYSTEM",
        "role instruction",
        ["bioinformatics expert specializing in", "spaceflight transcriptomics and", "mouse gene expression"],
        COLORS["sky"],
    )
    draw_prompt_section(
        draw,
        card_x + 34,
        card_y + 252,
        card_w - 68,
        408,
        "USER",
        "sample prompt surface",
        [
            f"A mouse tissue RNA-seq sample shows the highest expression variability in {data['n_genes']} genes.",
            "",
            "Gene list is ordered by absolute z-score:",
            "> Fbxo5 (ENSMUSG...), Mki67 (ENSMUSG...), Top2a (ENSMUSG...), ...",
            "",
            "Predict whether this sample is:",
            "(A) Flight    (B) Ground",
        ],
        COLORS["teal"],
    )
    draw_prompt_section(
        draw,
        card_x + 34,
        card_y + 694,
        card_w - 68,
        180,
        "FORMAT",
        "machine-readable answer",
        ["> A 0.82", "or", "> B 0.65"],
        COLORS["green"],
    )
    rounded(draw, (card_x + 34, card_y + 912, card_x + card_w - 34, card_y + 1016), 24, "#121B28", "#2A394D", 1)
    text(draw, (card_x + 66, card_y + 950), "Input surface", F["small_bold"], COLORS["violet"])
    paragraph(
        draw,
        (card_x + 66, card_y + 984),
        "Ranked gene symbols and Ensembl IDs are the full test-time observation passed to the language model.",
        F["tiny_bold"],
        card_w - 132,
        COLORS["text"],
        leading=4,
    )

    stat_y = y1 - 228
    rounded(draw, (x0 + 58, stat_y, x1 - 58, y1 - 58), 30, COLORS["panel2"], COLORS["teal"], 2)
    text(draw, (x0 + 96, stat_y + 38), f"{data['total_prompts']} prompts across {data['n_tasks']} tasks", F["h3"], COLORS["teal"])
    paragraph(
        draw,
        (x0 + 96, stat_y + 86),
        "Prompt generation uses the documented gene filter and top-50 ordering protocol.",
        F["small_bold"],
        x1 - x0 - 192,
        COLORS["text"],
        leading=5,
    )


def draw_pipeline_step(
    draw: ImageDraw.ImageDraw,
    x: int,
    y: int,
    w: int,
    title: str,
    body: str,
    accent: str,
    index: int,
) -> None:
    rounded(draw, (x, y, x + w, y + 136), 28, COLORS["panel2"], "#2A394D", 2)
    rounded(draw, (x + 28, y + 34, x + 96, y + 102), 18, rgba(accent, 50), rgba(accent, 180), 2)
    text(draw, (x + 62, y + 58), str(index), F["small_bold"], accent, "mm")
    text(draw, (x + 126, y + 26), title, F["h3"], COLORS["text"])
    paragraph(draw, (x + 126, y + 74), body, F["tiny"], w - 166, COLORS["muted"], leading=4)


def draw_pipeline_panel(draw: ImageDraw.ImageDraw, box: tuple[int, int, int, int]) -> None:
    x0, y0, x1, y1 = box
    rounded(draw, box, 34, COLORS["panel"], "#29374A", 2)
    text(draw, (x0 + 48, y0 + 46), "B. Prompt diagnostic pipeline", F["h2"], COLORS["text"])
    paragraph(
        draw,
        (x0 + 48, y0 + 98),
        "The text response is converted into the same held-out scoring grammar used elsewhere in the deck.",
        F["small"],
        x1 - x0 - 96,
        COLORS["muted"],
        leading=6,
    )
    steps = [
        ("Train-fold gene selection", "Select high-variance protein-coding genes inside the training fold.", COLORS["violet"]),
        ("Top-50 gene prompt", "Serialize ranked symbols and mouse Ensembl IDs into natural language.", COLORS["teal"]),
        ("Text LLM response", "Collect the provider answer as an A/B label plus confidence.", COLORS["sky"]),
        ("Parser", "A c maps to P(flight)=c; B c maps to P(flight)=1-c.", COLORS["amber"]),
        ("AUROC input", "Evaluate flight probabilities against hidden mission labels.", COLORS["green"]),
    ]
    step_x = x0 + 74
    step_w = x1 - x0 - 148
    y = y0 + 226
    for i, (title, body, accent) in enumerate(steps, start=1):
        draw_pipeline_step(draw, step_x, y, step_w, title, body, accent, i)
        if i < len(steps):
            arrow(draw, (step_x + step_w // 2, y + 146), (step_x + step_w // 2, y + 184), accent, width=4)
        y += 196

    bridge_y = y1 - 260
    rounded(draw, (x0 + 74, bridge_y, x1 - 74, y1 - 60), 32, "#121B28", rgba("#73A7FF", 145), 2)
    text(draw, (x0 + 112, bridge_y + 42), "Takeaway", F["h3"], COLORS["sky"])
    paragraph(
        draw,
        (x0 + 112, bridge_y + 92),
        "Treat the row as a prompt/input-surface diagnostic: prompt design, parser behavior, and task signal are coupled.",
        F["body_bold"],
        x1 - x0 - 224,
        COLORS["text"],
        leading=6,
    )


def chart_x(x0: int, width: int, value: float) -> int:
    lo, hi = 0.35, 0.65
    return int(x0 + (max(lo, min(hi, value)) - lo) / (hi - lo) * width)


def draw_diagnostic_chart(draw: ImageDraw.ImageDraw, box: tuple[int, int, int, int], providers: list[dict[str, object]]) -> None:
    x0, y0, x1, y1 = box
    rounded(draw, box, 30, "#121B28", "#2A394D", 2)
    text(draw, (x0 + 34, y0 + 38), "Diagnostic readout", F["h3"], COLORS["text"])
    text(draw, (x0 + 34, y0 + 82), "task-mean AUROC with parse rate", F["tiny_bold"], COLORS["muted"])

    axis_x = x0 + 300
    axis_w = x1 - axis_x - 182
    axis_y = y0 + 150
    draw.line((axis_x, axis_y, axis_x + axis_w, axis_y), fill=rgba("#98A7BA", 90), width=2)
    for val, label in [(0.35, "0.35"), (0.50, "chance"), (0.65, "0.65")]:
        px = chart_x(axis_x, axis_w, val)
        draw.line((px, axis_y - 12, px, axis_y + 210), fill=rgba("#98A7BA", 80 if val != 0.5 else 150), width=2 if val == 0.5 else 1)
        text(draw, (px, axis_y - 42), label, F["micro_bold"], COLORS["axis"], "mm")

    colors = [COLORS["teal"], COLORS["violet"], COLORS["amber"]]
    for idx, provider in enumerate(providers):
        y = y0 + 160 + idx * 76
        name = str(provider["provider"]).replace("Llama-3.3-70B", "Llama-3.3").replace("Gemini-2.5-Flash", "Gemini-2.5")
        auroc = float(provider["task_mean_auroc"])
        parse_rate = float(provider["parse_rate"])
        color = colors[idx]
        text(draw, (x0 + 34, y + 10), name, F["tiny_bold"], COLORS["text"])
        px = chart_x(axis_x, axis_w, auroc)
        draw.line((axis_x, y + 20, axis_x + axis_w, y + 20), fill=rgba("#2A3546", 130), width=10)
        draw.ellipse((px - 18, y + 2, px + 18, y + 38), fill=rgba(color, 220), outline=COLORS["white"], width=2)
        text(draw, (axis_x + axis_w + 30, y + 4), f"{auroc:.3f}", F["tiny_bold"], color)
        text(draw, (axis_x + axis_w + 30, y + 34), f"parse {parse_rate:.0%}", F["micro_bold"], COLORS["muted"])


def draw_readout_card(
    draw: ImageDraw.ImageDraw,
    x: int,
    y: int,
    w: int,
    title: str,
    body: str,
    accent: str,
) -> None:
    rounded(draw, (x, y, x + w, y + 152), 28, COLORS["panel2"], "#2A394D", 2)
    draw.line((x + 24, y + 31, x + 24, y + 120), fill=rgba(accent, 190), width=6)
    text(draw, (x + 54, y + 30), title, F["h3"], COLORS["text"])
    paragraph(draw, (x + 54, y + 78), body, F["tiny"], w - 92, COLORS["muted"], leading=4)


def draw_readout_panel(draw: ImageDraw.ImageDraw, box: tuple[int, int, int, int], data: dict[str, object]) -> None:
    x0, y0, x1, y1 = box
    rounded(draw, box, 34, COLORS["panel"], "#29374A", 2)
    text(draw, (x0 + 50, y0 + 46), "C. Diagnostic decoder", F["h2"], COLORS["text"])
    paragraph(
        draw,
        (x0 + 50, y0 + 98),
        "This slide separates prompt behavior from expression-encoder model behavior.",
        F["small"],
        x1 - x0 - 100,
        COLORS["muted"],
        leading=6,
    )

    draw_diagnostic_chart(draw, (x0 + 58, y0 + 195, x1 - 58, y0 + 645), data["providers"])

    cards = [
        ("Input surface", "Ranked gene text plus the model's language-form biology prior.", COLORS["teal"]),
        ("Observed behavior", "Scores cluster near chance while parse rates differ by provider.", COLORS["amber"]),
        ("Diagnostic use", "Tests whether text-form domain prior carries task signal under the same split.", COLORS["violet"]),
        ("Parser frame", "Prompt design and parsing remain part of the measured result surface.", COLORS["green"]),
    ]
    y = y0 + 692
    for title, body, accent in cards:
        draw_readout_card(draw, x0 + 58, y, x1 - x0 - 116, title, body, accent)
        y += 176


def draw_footer(draw: ImageDraw.ImageDraw) -> None:
    y = 1872
    rounded(draw, (136, y, 3704, 2042), 32, "#101823", "#2A394D", 2)
    text(draw, (180, y + 38), "Slide 25 readout", F["small_bold"], COLORS["sky"])
    paragraph(
        draw,
        (180, y + 76),
        "Prompt-based LLM rows are best read as diagnostics of the text input surface, parser, and language prior.",
        F["body_bold"],
        2260,
        COLORS["text"],
        leading=6,
    )
    paragraph(
        draw,
        (2610, y + 46),
        "Next: compare the tested model tiers after the audience understands each input surface.",
        F["tiny_bold"],
        920,
        COLORS["muted"],
        leading=5,
    )
    text(
        draw,
        (140, 2102),
        "Takeaway: prompt rows test the text input surface, parser, and language prior as one diagnostic package.",
        F["micro"],
        COLORS["dim"],
    )
    text(draw, (3704, 2102), "SPACEBIO-BENCH DETAILED DECK / MODELS", F["micro_bold"], COLORS["dim"], "ra")


def main() -> None:
    data = load_data()
    canvas = Image.new("RGBA", (W, H), COLORS["bg"])
    draw_background(canvas)
    draw = ImageDraw.Draw(canvas)

    draw_header(draw, data)
    draw_prompt_panel(draw, (136, 345, 1196, 1815), data)
    draw_pipeline_panel(draw, (1236, 345, 2546, 1815))
    draw_readout_panel(draw, (2586, 345, 3704, 1815), data)
    draw_footer(draw)

    rgb = canvas.convert("RGB")
    rgb.save(OUT_PATH, quality=95)
    rgb.convert("L").convert("RGB").save(GRAY_PATH, quality=95)
    MANIFEST_PATH.write_text(
        json.dumps(
            {
                "asset": str(OUT_PATH.relative_to(ROOT)),
                "grayscale": str(GRAY_PATH.relative_to(ROOT)),
                "source_files": [
                    str(PROMPT_SPEC.relative_to(ROOT)),
                    str(PROMPT_SUMMARY.relative_to(ROOT)),
                    "v2/evaluation/parse_metrics_*_A*.json",
                    "evaluation/submission_*_zeroshot_A*_eval.json",
                ],
                "metrics": data,
            },
            indent=2,
        )
        + "\n"
    )
    print(json.dumps({"asset": str(OUT_PATH.relative_to(ROOT)), "manifest": str(MANIFEST_PATH.relative_to(ROOT))}, indent=2))


if __name__ == "__main__":
    main()
