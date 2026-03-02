#!/usr/bin/env python3
"""
generate_tasks.py — GeneLab_benchmark: Task Split Generation

Generates LOMO (Leave-One-Mission-Out) train/test splits for each Category A task.
Applies LOMO-aware variance filtering INSIDE the loop (DD-03).

Design decisions implemented:
  DD-01: Feature = log2(normalized counts), NOT LFC
  DD-03: Variance filter inside LOMO loop (train missions only)
  DD-04: Independence unit = Mission for Category A (LOMO)

Output structure:
  tasks/A1_liver_lomo/
    fold_RR-1_test/
      train_X.csv         # (n_train_samples, n_selected_genes)
      train_y.csv         # (n_train_samples,) — 0/1 binary Flight label
      test_X.csv          # (n_test_samples, n_selected_genes)
      test_y.csv          # (n_test_samples,) — 0/1 binary Flight label
      train_meta.csv      # sample metadata (mission, osd_id, label_raw)
      test_meta.csv
      selected_genes.txt  # genes selected by train-only variance filter
      fold_info.json      # missions, n_train/test, label distribution
    fold_RR-3_test/
      ...
    task_info.json        # task-level summary

Usage:
  python generate_tasks.py --task A1          # generate Task A1 (liver LOMO)
  python generate_tasks.py --tissue liver     # all tasks for liver
  python generate_tasks.py --all              # all tasks
  python generate_tasks.py --list             # list available tasks
  python generate_tasks.py --task A1 --dry-run  # show what would be generated
"""

import json
import argparse
import numpy as np
import pandas as pd
from pathlib import Path
from datetime import datetime

# ── Paths ──────────────────────────────────────────────────────────────────────
BASE_DIR = Path(__file__).resolve().parent.parent
PROCESSED_DIR = BASE_DIR / "processed" / "A_detection"
TASKS_DIR = BASE_DIR / "tasks"
SPLITS_DIR = BASE_DIR / "splits"

# ── Task definitions ───────────────────────────────────────────────────────────
# Each task: tissue + missions to include + split type
TASK_DEFINITIONS = {
    "A1": {
        "name": "Spaceflight detection — Liver",
        "tissue": "liver",
        "split": "LOMO",
        "metric": "AUROC",
        "task_type": "binary",
        "track": "2a",    # C57BL/6J only
        "description": "Binary classification: Flight vs. Ground Control (per sample)",
        "notes": "Primary task. 6 missions. LOMO independence unit = mission.",
    },
    "A2": {
        "name": "Spaceflight detection — Gastrocnemius",
        "tissue": "gastrocnemius",
        "split": "LOMO",
        "metric": "AUROC",
        "task_type": "binary",
        "track": "2a",
        "description": "Binary classification: Flight vs. Ground Control (per sample)",
    },
    "A3": {
        "name": "Spaceflight detection — Kidney",
        "tissue": "kidney",
        "split": "LOMO",
        "metric": "AUROC",
        "task_type": "binary",
        "track": "2a",
        "description": "Binary classification: Flight vs. Ground Control (per sample)",
    },
    "A4": {
        "name": "Spaceflight detection — Thymus",
        "tissue": "thymus",
        "split": "LOMO",
        "metric": "AUROC",
        "task_type": "binary",
        "track": "2a",
        "description": "Binary classification: Flight vs. Ground Control (per sample)",
    },
    "A5": {
        "name": "Spaceflight detection — Skin",
        "tissue": "skin",
        "split": "LOMO",
        "metric": "AUROC",
        "task_type": "binary",
        "track": "2a",
        "description": "3 missions (MHU-2, RR-6, RR-7). MHU-2=dorsal+femoral merged. RR-7=C57BL/6J subset of OSD-254.",
    },
    "A6": {
        "name": "Spaceflight detection — Eye/Retina",
        "tissue": "eye",
        "split": "LOMO",
        "metric": "AUROC",
        "task_type": "binary",
        "track": "2a",
        "description": "Binary classification: Flight vs. Ground Control (per sample)",
    },
}

# ── QC thresholds ──────────────────────────────────────────────────────────────
VARIANCE_PERCENTILE = 0.25   # keep top 75% by variance (train-only)
MIN_SAMPLES_PER_FOLD = 4     # minimum test samples per fold
FLIGHT_LABEL = "Flight"      # canonical Flight label
GROUND_LABELS = {"GC", "VC"} # GC or VC → label 0 (Ground)
# Note: BC (Basal), AG (Artificial Gravity) are excluded from binary classification


# ── Core functions ─────────────────────────────────────────────────────────────

def load_tissue_data(tissue: str, batch_corrected: bool = False) -> tuple[pd.DataFrame, pd.DataFrame] | None:
    """
    Load tissue-wide log2 normalized counts + metadata.
    Returns (counts_df, meta_df) or None if file not found.

    counts_df: rows=samples, columns=genes
    meta_df: rows=samples, with 'mission', 'label' columns

    batch_corrected: if True, loads ComBat-seq corrected counts instead of raw normalized.
    """
    tissue_dir = PROCESSED_DIR / tissue
    meta_path = tissue_dir / f"{tissue}_all_missions_metadata.csv"

    if batch_corrected:
        counts_path = tissue_dir / f"{tissue}_combat_seq_log2_norm.csv"
        if not counts_path.exists():
            print(f"  [ERROR] ComBat-seq corrected counts not found: {counts_path}")
            print(f"          Run: Rscript scripts/batch_correct.R --tissue {tissue}")
            return None
        print(f"  Using ComBat-seq batch-corrected counts: {counts_path.name}")
    else:
        counts_path = tissue_dir / f"{tissue}_all_missions_log2_norm.csv"
        if not counts_path.exists():
            print(f"  [ERROR] Tissue-wide counts not found: {counts_path}")
            print(f"          Run: python scripts/quality_filter.py --tissue {tissue}")
            return None

    if not meta_path.exists():
        print(f"  [ERROR] Tissue-wide metadata not found: {meta_path}")
        return None

    counts = pd.read_csv(counts_path, index_col=0)
    meta = pd.read_csv(meta_path, index_col=0)

    # Drop non-gene columns that quality_filter.py appended (mission, osd_id, etc.)
    non_gene_cols = [c for c in counts.columns if c in {"mission", "osd_id", "label"}]
    if non_gene_cols:
        counts = counts.drop(columns=non_gene_cols)

    # Keep only numeric gene columns
    counts = counts.select_dtypes(include=[np.number])

    # ComBat-seq output: R's cbind() prepends list name (OSD-XX.) to avoid col conflicts.
    # Strip "OSD-XXXXX." prefix from sample index if present.
    if batch_corrected:
        import re
        counts.index = [re.sub(r'^OSD-\d+\.', '', s) for s in counts.index]

    # Align with metadata index (only QC-passing samples)
    common_samples = counts.index.intersection(meta.index)
    if len(common_samples) < len(meta):
        print(f"  [WARN] {len(meta) - len(common_samples)} metadata samples not in counts — dropping")
    if len(common_samples) == 0:
        print(f"  [ERROR] No overlapping samples between counts and metadata.")
        print(f"          Check sample name formats.")
        return None
    counts = counts.loc[common_samples]
    meta   = meta.loc[common_samples]

    print(f"  Loaded: {counts.shape[0]} samples × {counts.shape[1]} genes")
    return counts, meta


def get_binary_labels(meta: pd.DataFrame, flight_label: str = FLIGHT_LABEL,
                      ground_labels: set = GROUND_LABELS) -> pd.Series:
    """
    Extract binary labels: 1=Flight, 0=Ground (GC or VC).
    Returns pd.Series with same index as meta, NaN for excluded samples (BC, AG, Unknown).
    """
    label_col = "label" if "label" in meta.columns else None
    if label_col is None:
        # Try to find label column
        for col in ["label_raw", "condition", "group"]:
            if col in meta.columns:
                label_col = col
                break

    if label_col is None:
        raise ValueError(f"No label column found in metadata. Columns: {list(meta.columns)}")

    labels_raw = meta[label_col]
    binary = pd.Series(np.nan, index=meta.index)
    binary[labels_raw == flight_label] = 1
    binary[labels_raw.isin(ground_labels)] = 0
    return binary


def variance_filter_train(train_X: pd.DataFrame,
                          var_percentile: float = VARIANCE_PERCENTILE) -> list[str]:
    """
    Compute variance filter on TRAIN data only. Returns list of selected gene names.
    This is the DD-03 implementation: variance filter inside LOMO loop.

    var_percentile: filter out bottom var_percentile fraction.
    VARIANCE_PERCENTILE=0.25 → keep top 75% by variance.
    """
    gene_var = train_X.var(axis=0)
    threshold = gene_var.quantile(var_percentile)
    selected = gene_var[gene_var >= threshold].index.tolist()
    return selected


def lomo_split(counts: pd.DataFrame, meta: pd.DataFrame,
               task_def: dict, output_dir: Path,
               dry_run: bool = False, verbose: bool = True) -> dict:
    """
    Generate LOMO folds for a task.
    For each mission M: train = all missions except M, test = M.
    Applies variance filter on train only (DD-03).
    Binary labels: Flight=1, Ground=0. BC/AG excluded.

    Returns summary dict.
    """
    # Get binary labels
    binary_labels = get_binary_labels(meta)
    n_excluded = binary_labels.isna().sum()

    # Get mission list
    if "mission" not in meta.columns:
        raise ValueError("No 'mission' column in metadata.")

    missions = sorted(meta["mission"].unique())

    if verbose:
        print(f"  Missions ({len(missions)}): {missions}")
        if n_excluded > 0:
            excluded_labels = meta.loc[binary_labels.isna(), "label"].value_counts()
            print(f"  Excluded from binary task: {n_excluded} samples {dict(excluded_labels)}")

    folds_summary = []

    for test_mission in missions:
        fold_name = f"fold_{test_mission.replace(' ', '_')}_test"
        fold_dir = output_dir / fold_name

        # Split: samples in test_mission vs. all others
        test_mask = meta["mission"] == test_mission
        train_mask = ~test_mask

        # Remove BC/AG (non-binary) samples from train
        valid_binary = ~binary_labels.isna()
        train_idx = meta.index[train_mask & valid_binary]
        test_idx = meta.index[test_mask & valid_binary]

        if len(test_idx) < MIN_SAMPLES_PER_FOLD:
            if verbose:
                print(f"  [SKIP] {test_mission} test fold: only {len(test_idx)} binary-labeled samples")
            continue

        # Check label balance in test
        test_labels = binary_labels[test_idx]
        n_flight_test = (test_labels == 1).sum()
        n_ground_test = (test_labels == 0).sum()

        if n_flight_test == 0 or n_ground_test == 0:
            if verbose:
                print(f"  [SKIP] {test_mission} test fold: missing Flight or Ground samples "
                      f"(Flight={n_flight_test}, Ground={n_ground_test})")
            continue

        # Extract feature matrices
        train_X = counts.loc[train_idx]
        test_X = counts.loc[test_idx]
        train_y = binary_labels[train_idx]
        test_y = binary_labels[test_idx]

        # ── DD-03: Variance filter on TRAIN only ──────────────────────────────
        selected_genes = variance_filter_train(train_X, VARIANCE_PERCENTILE)
        train_X_filtered = train_X[selected_genes]
        test_X_filtered = test_X[selected_genes]

        # ── Cross-mission z-score (train-only statistics, no leakage) ─────────
        # Per-gene z-score removes cross-study library size and scale differences.
        # Mean and std computed from train samples only; applied to both train+test.
        train_mean = train_X_filtered.mean(axis=0)
        train_std  = train_X_filtered.std(axis=0).replace(0, 1)  # avoid div by 0
        train_X_filtered = (train_X_filtered - train_mean) / train_std
        test_X_filtered  = (test_X_filtered  - train_mean) / train_std

        train_label_dist = {
            "Flight": int((train_y == 1).sum()),
            "Ground": int((train_y == 0).sum()),
        }
        test_label_dist = {
            "Flight": int(n_flight_test),
            "Ground": int(n_ground_test),
        }

        if verbose:
            print(f"\n  Fold: test={test_mission}")
            print(f"    Train: {len(train_idx)} samples, {len(selected_genes)} genes (after var filter)")
            print(f"    Test:  {len(test_idx)} samples")
            print(f"    Train labels: {train_label_dist}")
            print(f"    Test labels:  {test_label_dist}")

        fold_info = {
            "test_mission": test_mission,
            "train_missions": [m for m in missions if m != test_mission],
            "n_train": len(train_idx),
            "n_test": len(test_idx),
            "n_genes_before_var_filter": train_X.shape[1],
            "n_genes_after_var_filter": len(selected_genes),
            "var_percentile_cutoff": VARIANCE_PERCENTILE,
            "train_label_distribution": train_label_dist,
            "test_label_distribution": test_label_dist,
            "excluded_bc_ag_samples": int(n_excluded),
        }

        if not dry_run:
            fold_dir.mkdir(parents=True, exist_ok=True)

            train_X_filtered.to_csv(fold_dir / "train_X.csv")
            test_X_filtered.to_csv(fold_dir / "test_X.csv")
            train_y.to_csv(fold_dir / "train_y.csv", header=True)
            test_y.to_csv(fold_dir / "test_y.csv", header=True)

            # Save metadata
            train_meta = meta.loc[train_idx]
            test_meta = meta.loc[test_idx]
            train_meta.to_csv(fold_dir / "train_meta.csv")
            test_meta.to_csv(fold_dir / "test_meta.csv")

            # Save selected genes
            (fold_dir / "selected_genes.txt").write_text("\n".join(selected_genes))

            # Save fold info
            (fold_dir / "fold_info.json").write_text(
                json.dumps(fold_info, indent=2)
            )

        folds_summary.append(fold_info)

    return {
        "task": task_def["name"],
        "tissue": task_def["tissue"],
        "split": "LOMO",
        "n_missions": len(missions),
        "n_folds_generated": len(folds_summary),
        "folds": folds_summary,
    }


def generate_task(task_id: str, dry_run: bool = False, verbose: bool = True,
                  batch_corrected: bool = False) -> dict | None:
    """Generate all LOMO folds for a single task."""
    if task_id not in TASK_DEFINITIONS:
        print(f"[ERROR] Unknown task: {task_id}. Available: {list(TASK_DEFINITIONS)}")
        return None

    task_def = TASK_DEFINITIONS[task_id]
    tissue = task_def["tissue"]

    print(f"\n{'='*60}")
    print(f"Task {task_id}: {task_def['name']}")
    if batch_corrected:
        print(f"  [ComBat-seq batch-corrected]")
    print(f"{'='*60}")

    # Load data
    result = load_tissue_data(tissue, batch_corrected=batch_corrected)
    if result is None:
        return None
    counts, meta = result

    # Output directory — add suffix for batch-corrected variant
    suffix = "_combat" if batch_corrected else ""
    output_dir = TASKS_DIR / f"{task_id}_{tissue}_lomo{suffix}"
    if not dry_run:
        output_dir.mkdir(parents=True, exist_ok=True)

    # Generate LOMO folds
    summary = lomo_split(counts, meta, task_def, output_dir,
                         dry_run=dry_run, verbose=verbose)
    summary["task_id"] = task_id
    summary["generated_at"] = datetime.now().isoformat()
    summary["dry_run"] = dry_run

    if not dry_run:
        task_info_path = output_dir / "task_info.json"
        task_info_path.write_text(json.dumps(summary, indent=2))
        print(f"\n  ✓ Task {task_id} saved to: {output_dir}")
        print(f"    {summary['n_folds_generated']} folds generated")

    return summary


# ── CLI ────────────────────────────────────────────────────────────────────────

def parse_args():
    parser = argparse.ArgumentParser(
        description="Generate LOMO task splits for GeneLab_benchmark"
    )
    parser.add_argument("--task", type=str,
                        help="Task ID to generate (e.g., A1)")
    parser.add_argument("--tissue", type=str,
                        help="Generate all tasks for a tissue (e.g., liver)")
    parser.add_argument("--all", action="store_true",
                        help="Generate all tasks")
    parser.add_argument("--list", action="store_true",
                        help="List available tasks")
    parser.add_argument("--dry-run", action="store_true",
                        help="Show what would be generated without writing files")
    parser.add_argument("--quiet", action="store_true",
                        help="Suppress per-fold verbose output")
    parser.add_argument("--batch-corrected", action="store_true",
                        help="Use ComBat-seq batch-corrected counts (requires batch_correct.R first)")
    return parser.parse_args()


def main():
    args = parse_args()
    verbose = not args.quiet

    if args.list:
        print("\nAvailable tasks:")
        for tid, tdef in TASK_DEFINITIONS.items():
            print(f"  {tid}: {tdef['name']}")
            print(f"       Tissue: {tdef['tissue']}, Split: {tdef['split']}")
        return

    tasks_to_run = []

    if args.task:
        tasks_to_run = [args.task.upper()]
    elif args.tissue:
        tasks_to_run = [
            tid for tid, tdef in TASK_DEFINITIONS.items()
            if tdef["tissue"] == args.tissue
        ]
        if not tasks_to_run:
            print(f"[ERROR] No tasks found for tissue: {args.tissue}")
            return
    elif args.all:
        tasks_to_run = list(TASK_DEFINITIONS.keys())
    else:
        # Default: run A1
        tasks_to_run = ["A1"]
        print("No task specified. Running A1 (liver) by default.")
        print("Use --task, --tissue, --all to specify tasks.")

    if args.dry_run:
        print("\n[DRY RUN] No files will be written.\n")

    batch_corrected = getattr(args, "batch_corrected", False)

    all_summaries = {}
    for task_id in tasks_to_run:
        summary = generate_task(task_id, dry_run=args.dry_run, verbose=verbose,
                                batch_corrected=batch_corrected)
        if summary:
            all_summaries[task_id] = summary

    # Final summary
    print(f"\n{'='*60}")
    print(f"Task Generation Complete")
    print(f"{'='*60}")
    for task_id, summary in all_summaries.items():
        n_folds = summary["n_folds_generated"]
        n_missions = summary["n_missions"]
        print(f"  {task_id}: {n_folds}/{n_missions} folds generated")
        if not args.dry_run:
            outdir = TASKS_DIR / f"{task_id}_{summary['tissue']}_lomo"
            print(f"    → {outdir}")

    print(f"\nNext: python scripts/run_baselines.py --task A1")


if __name__ == "__main__":
    main()
