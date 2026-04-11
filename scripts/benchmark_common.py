#!/usr/bin/env python3
"""
benchmark_common.py — shared benchmark constants and output helpers.

Keeps benchmark-wide standards in one place so runners, submission tooling,
and evaluators cannot silently drift apart.
"""

from __future__ import annotations

import re
from pathlib import Path


PHASE1_AUROC_THRESHOLD = 0.700
PHASE1_CI_LOWER_THRESHOLD = 0.600
PHASE1_PERM_P_THRESHOLD = 0.050


def safe_slug(value: str) -> str:
    """Convert a string into a filesystem-safe identifier."""
    slug = re.sub(r"[^A-Za-z0-9._-]+", "_", value.strip())
    slug = re.sub(r"_+", "_", slug).strip("_")
    return slug or "default"


def task_variant_suffix(task_id: str, task_dir: Path, tasks_dir: Path) -> str:
    """
    Derive a stable output suffix for task-directory variants.

    Canonical task directories return an empty suffix. Variants that extend the
    canonical directory name return only the extra suffix, preserving existing
    names like `_combat` and `_iss_only`.
    """
    candidates = sorted(
        (
            d.name
            for d in tasks_dir.iterdir()
            if d.is_dir() and d.name.startswith(f"{task_id}_")
        ),
        key=lambda name: (len(name), name),
    )
    canonical_name = candidates[0] if candidates else task_dir.name

    if task_dir.name == canonical_name:
        return ""
    if task_dir.name.startswith(canonical_name):
        return task_dir.name[len(canonical_name):]
    return f"_{safe_slug(task_dir.name)}"
