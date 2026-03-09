#!/usr/bin/env python3
"""
llm_parse_responses.py — GeneLab_benchmark v2.0: Tier 3 Response Parsing

Parses LLM responses and generates submission JSONs for evaluation.

Output:
  v2/evaluation/
    submission_GPT4o_zeroshot_{task_id}.json

Usage:
  python v2/scripts/llm_parse_responses.py --task A4
  python v2/scripts/llm_parse_responses.py --all
"""

import json
import re
import sys
import argparse
from pathlib import Path
from datetime import datetime

V2_DIR = Path(__file__).resolve().parent.parent
PROJECT_DIR = V2_DIR.parent

RESPONSE_BASE = V2_DIR / "processed" / "llm_responses"
EVAL_DIR = V2_DIR / "evaluation"

# Provider → model name mapping for submission
PROVIDER_MODEL_NAMES = {
    "gemini-2.5-flash": "Gemini-2.5-Flash_zeroshot",
    "groq-llama-3.3-70b": "Llama-3.3-70B_zeroshot",
    "deepseek-chat": "DeepSeek-V3_zeroshot",
    "together-llama-3.3-70b": "Llama-3.3-70B-Together_zeroshot",
    "gpt-4o": "GPT-4o_zeroshot",
    "claude-sonnet": "Claude-Sonnet_zeroshot",
}


def parse_llm_response(response: str) -> tuple:
    """
    Parse LLM response for prediction letter and confidence.
    Returns (label, flight_probability).
    label: 1 = Flight (A), 0 = Ground (B)
    flight_probability: confidence if Flight, 1-confidence if Ground

    Handles multiple formats:
    1. "A 0.82" or "B 0.65" (exact format)
    2. "A, 0.82" or "A: 0.82" (alternate separators)
    3. "**(A) Flight** with confidence 0.85" (verbose with bold)
    4. "classify this sample as (A) Flight" (no explicit confidence → 0.75 default)
    5. "Flight" or "Ground" keywords in last sentences (truncated CoT fallback)
    """
    if response is None:
        return None, 0.5

    resp_upper = response.upper()

    # 1. Exact format: "A 0.82" or "B 0.35"
    pattern = r'\b([AB])\s+(0\.\d+|1\.0|0|1)\b'
    matches = re.findall(pattern, resp_upper)

    # 2. Alternate separators: "A, 0.82" or "A: 0.82"
    if not matches:
        pattern2 = r'\b([AB])[,:\s]+(0\.\d+|1\.0|0|1)\b'
        matches = re.findall(pattern2, resp_upper)

    # 3. Verbose: "(A) Flight" with separate confidence
    if not matches:
        letter_match = re.search(r'\(([AB])\)\s*(FLIGHT|GROUND)', resp_upper)
        conf_match = re.search(r'(?:CONFIDENCE|CONF)[:\s]*(0\.\d+|1\.0)', resp_upper)
        if letter_match:
            letter = letter_match.group(1)
            conf = float(conf_match.group(1)) if conf_match else 0.75
            matches = [(letter, str(conf))]

    # 4. Keyword fallback for truncated CoT responses
    if not matches:
        # Check last 500 chars for classification keywords
        tail = resp_upper[-500:]
        flight_keywords = len(re.findall(r'\bFLIGHT\b', tail))
        ground_keywords = len(re.findall(r'\bGROUND\b', tail))
        # Also check for "classify.*as.*flight/ground"
        if re.search(r'CLASSIF\w*.*\b(FLIGHT|OPTION A|ANSWER.*A)\b', tail):
            flight_keywords += 2
        if re.search(r'CLASSIF\w*.*\b(GROUND|OPTION B|ANSWER.*B)\b', tail):
            ground_keywords += 2

        if flight_keywords > ground_keywords:
            matches = [("A", "0.70")]
        elif ground_keywords > flight_keywords:
            matches = [("B", "0.70")]

    if not matches:
        return None, 0.5

    letter, conf_str = matches[-1]  # take last match (final answer in CoT)
    conf = float(conf_str)
    label = 1 if letter == "A" else 0
    flight_prob = conf if label == 1 else (1.0 - conf)
    return label, flight_prob


def process_task(task_id, provider_dir):
    """Parse all responses for a task and create submission JSON."""
    task_dir = RESPONSE_BASE / provider_dir / task_id
    if not task_dir.exists():
        print(f"  [SKIP] No responses for {task_id} in {provider_dir}")
        return None

    fold_files = sorted(task_dir.glob("fold_*_responses.json"))
    if not fold_files:
        print(f"  [SKIP] No response files in {task_dir}")
        return None

    # Detect model from first response file
    first_data = json.loads(fold_files[0].read_text())
    model_used = first_data.get("model", provider_dir)

    predictions = {}
    parse_stats = {"total": 0, "parsed": 0, "failed": 0}

    for fold_file in fold_files:
        fold_data = json.loads(fold_file.read_text())
        fold_name = fold_data["fold"]

        fold_preds = {}
        for resp in fold_data["responses"]:
            sample_id = resp["sample_id"]
            content = resp.get("content")
            parse_stats["total"] += 1

            label, flight_prob = parse_llm_response(content)
            if label is not None:
                parse_stats["parsed"] += 1
            else:
                parse_stats["failed"] += 1
                print(f"    PARSE FAIL: {sample_id} → {content[:80] if content else 'None'}")

            fold_preds[sample_id] = flight_prob

        predictions[f"{fold_name}"] = fold_preds

    model_name = PROVIDER_MODEL_NAMES.get(provider_dir, f"{provider_dir}_zeroshot")

    # Create submission JSON
    submission = {
        "task_id": task_id,
        "model_name": model_name,
        "model_description": f"{model_used}, zero-shot, top-50 genes, per-sample z-score",
        "tier": "3",
        "llm_variant": "zero_shot",
        "n_genes_in_prompt": 50,
        "submission_date": datetime.now().strftime("%Y-%m-%d"),
        "predictions": predictions,
    }

    EVAL_DIR.mkdir(parents=True, exist_ok=True)
    safe_name = provider_dir.replace("/", "_").replace(" ", "_")
    outpath = EVAL_DIR / f"submission_{safe_name}_zeroshot_{task_id}.json"
    outpath.write_text(json.dumps(submission, indent=2))

    print(f"  {task_id}: {parse_stats['parsed']}/{parse_stats['total']} parsed "
          f"({parse_stats['failed']} failed)")
    print(f"    → {outpath}")

    return submission


def parse_args():
    parser = argparse.ArgumentParser(description="Parse LLM responses → submission JSON")
    parser.add_argument("--task", type=str, help="Task ID")
    parser.add_argument("--all", action="store_true")
    parser.add_argument("--provider", type=str, default=None,
                        help="Provider directory name (e.g., gemini-2.0-flash, gpt-4o). "
                             "Auto-detects if only one provider exists.")
    return parser.parse_args()


def main():
    args = parse_args()

    if not args.task and not args.all:
        print("Specify --task or --all")
        return

    print("=" * 70)
    print("Tier 3: LLM Response Parsing → Submission JSON")
    print("=" * 70)

    # Find available providers
    if args.provider:
        providers = [args.provider]
    else:
        providers = sorted([d.name for d in RESPONSE_BASE.iterdir()
                           if d.is_dir()]) if RESPONSE_BASE.exists() else []

    if not providers:
        print("  No response directories found. Run llm_call_api.py first.")
        return

    for provider_dir in providers:
        print(f"\n  Provider: {provider_dir}")
        response_dir = RESPONSE_BASE / provider_dir

        if args.all:
            tasks = sorted([d.name for d in response_dir.iterdir()
                           if d.is_dir() and d.name.startswith("A")])
        else:
            tasks = [args.task]

        for task_id in tasks:
            process_task(task_id, provider_dir)

    print(f"\n  Submissions: {EVAL_DIR}")
    print(f"  Evaluate with: python scripts/evaluate_submission.py --submission <file>")


if __name__ == "__main__":
    main()
