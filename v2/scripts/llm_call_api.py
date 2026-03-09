#!/usr/bin/env python3
"""
llm_call_api.py — GeneLab_benchmark v2.0: Tier 3 API Calls

Calls LLM API with generated prompts and saves raw responses.
Supports multiple providers (all free-tier compatible):
  - Gemini (default, free): gemini-2.0-flash
  - Groq (free): llama-3.3-70b-versatile
  - DeepSeek (free trial): deepseek-chat
  - Together ($100 free credits): Llama-3.3-70B-Instruct-Turbo
  - OpenAI (paid): gpt-4o
  - Claude (paid): claude-sonnet

API keys loaded from ~/.api_keys if env vars not set.

Output:
  v2/processed/llm_responses/{model_name}/{task_id}/
    fold_{mission}_responses.json

Usage:
  python v2/scripts/llm_call_api.py --task A4                     # Gemini (default)
  python v2/scripts/llm_call_api.py --task A4 --provider groq     # Groq (free)
  python v2/scripts/llm_call_api.py --task A4 --provider deepseek # DeepSeek (free trial)
  python v2/scripts/llm_call_api.py --all
  python v2/scripts/llm_call_api.py --task A4 --dry-run
"""

import json
import os
import sys
import time
import argparse
from pathlib import Path
from datetime import datetime

V2_DIR = Path(__file__).resolve().parent.parent
PROJECT_DIR = V2_DIR.parent

PROMPT_DIR = V2_DIR / "processed" / "llm_prompts"

# ── API Key Loading ──────────────────────────────────────────────────────────

def load_api_keys():
    """Load API keys from ~/.api_keys if not already in environment."""
    keys_file = Path.home() / ".api_keys"
    if not keys_file.exists():
        return
    for line in keys_file.read_text().splitlines():
        line = line.strip()
        if not line or line.startswith("#"):
            continue
        if line.startswith("export "):
            line = line[7:]
        if "=" in line:
            key, val = line.split("=", 1)
            val = val.strip().strip('"').strip("'")
            if key and val and key not in os.environ:
                os.environ[key] = val


# ── Provider Config ───────────────────────────────────────────────────────────
PROVIDERS = {
    "gemini": {
        "model": "gemini-2.5-flash",
        "env_keys": ["GOOGLE_API_KEY", "GEMINI_API_KEY"],
        "dir_name": "gemini-2.5-flash",
        "temperature": 0,
        "max_tokens": 1000,
        "rate_limit_delay": 0.5,
    },
    "groq": {
        "model": "llama-3.3-70b-versatile",
        "env_keys": ["GROQ_API_KEY"],
        "dir_name": "groq-llama-3.3-70b",
        "base_url": "https://api.groq.com/openai/v1",
        "temperature": 0,
        "seed": 42,
        "max_tokens": 1000,
        "rate_limit_delay": 0.2,
    },
    "deepseek": {
        "model": "deepseek-chat",
        "env_keys": ["DEEPSEEK_API_KEY"],
        "dir_name": "deepseek-chat",
        "base_url": "https://api.deepseek.com",
        "temperature": 0,
        "max_tokens": 1000,
        "rate_limit_delay": 0.1,
    },
    "together": {
        "model": "meta-llama/Llama-3.3-70B-Instruct-Turbo",
        "env_keys": ["TOGETHER_API_KEY"],
        "dir_name": "together-llama-3.3-70b",
        "base_url": "https://api.together.xyz/v1",
        "temperature": 0,
        "max_tokens": 1000,
        "rate_limit_delay": 0.1,
    },
    "openai": {
        "model": "gpt-4o-2024-11-20",
        "env_keys": ["OPENAI_API_KEY"],
        "dir_name": "gpt-4o",
        "temperature": 0,
        "seed": 42,
        "max_tokens": 1000,
        "rate_limit_delay": 0.0,
    },
    "claude": {
        "model": "claude-sonnet-4-20250514",
        "env_keys": ["ANTHROPIC_API_KEY"],
        "dir_name": "claude-sonnet",
        "temperature": 0,
        "max_tokens": 1000,
        "rate_limit_delay": 0.0,
    },
}

DEFAULT_PROVIDER = "gemini"
RETRY_DELAY = 2
MAX_RETRIES = 3


# ── API Callers ───────────────────────────────────────────────────────────────

def call_gemini(client, model, system_prompt, user_prompt, config):
    """Call Google Gemini API."""
    for attempt in range(MAX_RETRIES):
        try:
            response = client.models.generate_content(
                model=model,
                contents=user_prompt,
                config={
                    "system_instruction": system_prompt,
                    "temperature": config.get("temperature", 0),
                    "max_output_tokens": config.get("max_tokens", 100),
                },
            )
            return {
                "content": response.text,
                "finish_reason": str(getattr(response.candidates[0], "finish_reason", "STOP")),
                "model": model,
                "usage": {
                    "prompt_tokens": getattr(response.usage_metadata, "prompt_token_count", 0),
                    "completion_tokens": getattr(response.usage_metadata, "candidates_token_count", 0),
                    "total_tokens": getattr(response.usage_metadata, "total_token_count", 0),
                },
            }
        except Exception as e:
            if attempt < MAX_RETRIES - 1:
                print(f"      Retry {attempt + 1}/{MAX_RETRIES}: {e}")
                time.sleep(RETRY_DELAY * (attempt + 1))
            else:
                return {"error": str(e), "content": None}


def call_openai(client, model, system_prompt, user_prompt, config):
    """Call OpenAI API."""
    for attempt in range(MAX_RETRIES):
        try:
            response = client.chat.completions.create(
                model=model,
                messages=[
                    {"role": "system", "content": system_prompt},
                    {"role": "user", "content": user_prompt},
                ],
                temperature=config.get("temperature", 0),
                seed=config.get("seed", 42),
                max_tokens=config.get("max_tokens", 100),
            )
            return {
                "content": response.choices[0].message.content,
                "finish_reason": response.choices[0].finish_reason,
                "model": response.model,
                "system_fingerprint": getattr(response, "system_fingerprint", None),
                "usage": {
                    "prompt_tokens": response.usage.prompt_tokens,
                    "completion_tokens": response.usage.completion_tokens,
                    "total_tokens": response.usage.total_tokens,
                },
            }
        except Exception as e:
            if attempt < MAX_RETRIES - 1:
                print(f"      Retry {attempt + 1}/{MAX_RETRIES}: {e}")
                time.sleep(RETRY_DELAY * (attempt + 1))
            else:
                return {"error": str(e), "content": None}


def call_claude(client, model, system_prompt, user_prompt, config):
    """Call Anthropic Claude API."""
    for attempt in range(MAX_RETRIES):
        try:
            response = client.messages.create(
                model=model,
                system=system_prompt,
                messages=[{"role": "user", "content": user_prompt}],
                temperature=config.get("temperature", 0),
                max_tokens=config.get("max_tokens", 100),
            )
            return {
                "content": response.content[0].text,
                "finish_reason": response.stop_reason,
                "model": response.model,
                "usage": {
                    "prompt_tokens": response.usage.input_tokens,
                    "completion_tokens": response.usage.output_tokens,
                    "total_tokens": response.usage.input_tokens + response.usage.output_tokens,
                },
            }
        except Exception as e:
            if attempt < MAX_RETRIES - 1:
                print(f"      Retry {attempt + 1}/{MAX_RETRIES}: {e}")
                time.sleep(RETRY_DELAY * (attempt + 1))
            else:
                return {"error": str(e), "content": None}


# ── Core ──────────────────────────────────────────────────────────────────────

def setup_client(provider_name):
    """Initialize API client for the given provider."""
    config = PROVIDERS[provider_name]

    # Find API key
    api_key = None
    for key_name in config["env_keys"]:
        api_key = os.environ.get(key_name)
        if api_key:
            break

    if not api_key:
        keys = " or ".join(config["env_keys"])
        raise ValueError(f"Set {keys} environment variable (or add to ~/.api_keys)")

    if provider_name == "gemini":
        from google import genai
        client = genai.Client(api_key=api_key)
        caller = call_gemini
    elif provider_name in ("groq", "deepseek", "together"):
        # OpenAI-compatible providers
        from openai import OpenAI
        client = OpenAI(api_key=api_key, base_url=config["base_url"])
        caller = call_openai
    elif provider_name == "openai":
        from openai import OpenAI
        client = OpenAI(api_key=api_key)
        caller = call_openai
    elif provider_name == "claude":
        import anthropic
        client = anthropic.Anthropic(api_key=api_key)
        caller = call_claude
    else:
        raise ValueError(f"Unknown provider: {provider_name}")

    return client, caller, config


def process_task(task_id, client, caller, config, provider_name, dry_run=False):
    """Process all folds for a task."""
    prompt_task_dir = PROMPT_DIR / task_id
    if not prompt_task_dir.exists():
        print(f"  [SKIP] No prompts for {task_id}. Run llm_generate_prompts.py first.")
        return None

    response_dir = V2_DIR / "processed" / "llm_responses" / config["dir_name"]
    output_dir = response_dir / task_id
    output_dir.mkdir(parents=True, exist_ok=True)

    fold_files = sorted(prompt_task_dir.glob("fold_*_prompts.json"))
    total_calls = 0
    total_tokens = 0

    for fold_file in fold_files:
        fold_data = json.loads(fold_file.read_text())
        fold_name = fold_data["fold"]
        prompts = fold_data["prompts"]
        n = len(prompts)

        # Check if already processed
        response_path = output_dir / f"{fold_name}_responses.json"
        if response_path.exists():
            existing = json.loads(response_path.read_text())
            if len(existing.get("responses", [])) == n:
                print(f"    {fold_name}: already complete ({n} responses)")
                total_calls += n
                continue

        if dry_run:
            print(f"    {fold_name}: {n} prompts (dry-run)")
            total_calls += n
            continue

        print(f"    {fold_name}: {n} prompts...")
        responses = []
        for i, prompt in enumerate(prompts):
            result = caller(
                client, config["model"],
                prompt["system_prompt"], prompt["user_prompt"],
                config
            )
            result["sample_id"] = prompt["sample_id"]
            result["true_label"] = prompt["true_label"]
            responses.append(result)

            if result.get("usage"):
                total_tokens += result["usage"].get("total_tokens", 0)

            if (i + 1) % 10 == 0:
                print(f"      {i + 1}/{n} done")

            # Rate limiting
            delay = config.get("rate_limit_delay", 0)
            if delay > 0:
                time.sleep(delay)

        # Save responses
        fold_output = {
            "fold": fold_name,
            "task_id": task_id,
            "provider": provider_name,
            "model": config["model"],
            "temperature": config.get("temperature", 0),
            "generated_at": datetime.now().isoformat(),
            "n_responses": len(responses),
            "responses": responses,
        }
        response_path.write_text(json.dumps(fold_output, indent=2))
        total_calls += len(responses)
        print(f"      Saved: {response_path}")

    return {"task_id": task_id, "total_calls": total_calls, "total_tokens": total_tokens}


def parse_args():
    parser = argparse.ArgumentParser(description="Call LLM API for Tier 3")
    parser.add_argument("--task", type=str, help="Task ID (e.g., A4)")
    parser.add_argument("--all", action="store_true", help="Process all tasks")
    parser.add_argument("--provider", type=str, default=DEFAULT_PROVIDER,
                        choices=list(PROVIDERS.keys()),
                        help=f"LLM provider (default: {DEFAULT_PROVIDER})")
    parser.add_argument("--dry-run", action="store_true",
                        help="Show what would be called without making API requests")
    return parser.parse_args()


def main():
    args = parse_args()

    if not args.task and not args.all:
        print("Specify --task A1/.../A6 or --all")
        return

    # Load API keys from ~/.api_keys
    load_api_keys()

    provider = args.provider
    config = PROVIDERS[provider]

    print("=" * 70)
    print(f"Tier 3: LLM API Calls ({'DRY RUN' if args.dry_run else 'LIVE'})")
    print(f"  Provider: {provider}")
    print(f"  Model: {config['model']}")
    print(f"  Temperature: {config.get('temperature', 0)}")
    print("=" * 70)

    client, caller = None, None
    if not args.dry_run:
        try:
            client, caller, config = setup_client(provider)
        except (ImportError, ValueError) as e:
            print(f"ERROR: {e}")
            install_hints = {
                "gemini": "pip install google-genai",
                "groq": "pip install openai>=1.0",
                "deepseek": "pip install openai>=1.0",
                "together": "pip install openai>=1.0",
                "openai": "pip install openai>=1.0",
                "claude": "pip install anthropic",
            }
            if provider in install_hints:
                print(f"  Install: {install_hints[provider]}")
            return

    # Find tasks with prompts
    if args.all:
        tasks = sorted([d.name for d in PROMPT_DIR.iterdir()
                        if d.is_dir() and d.name.startswith("A")])
    else:
        tasks = [args.task]

    for task_id in tasks:
        print(f"\n  Task {task_id}:")
        result = process_task(task_id, client, caller, config, provider, dry_run=args.dry_run)

    response_dir = V2_DIR / "processed" / "llm_responses" / config["dir_name"]
    print(f"\n  Done. Responses: {response_dir}")
    if args.dry_run:
        print("  (Dry run — no API calls made)")


if __name__ == "__main__":
    main()
