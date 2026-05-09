"""v8 Pillar 1 (BRIDGE) — supervised species-transfer experiment.

Task: SpaceOmicsBench I2 (pathway conservation I4→Twins, 452 pathways, binary).
Baseline: 8 I4 PBMC aggregated features.
v8 extension: add 6 mouse-tissue NES features (thymus, gastrocnemius, skin,
eye, liver, kidney; mean across missions). Does the mouse prior lift AUROC?

Evaluation: stratified 5-fold CV, Logistic Regression + Random Forest.
Report baseline AUROC, +mouse AUROC, delta with bootstrap CI.

Secondary: feature importance — which mouse tissue carries signal?
"""
from __future__ import annotations

import json
import os
from pathlib import Path

import numpy as np
import pandas as pd
from sklearn.ensemble import RandomForestClassifier
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import roc_auc_score
from sklearn.model_selection import StratifiedKFold
from sklearn.pipeline import Pipeline
from sklearn.preprocessing import StandardScaler

REPO_ROOT = Path(__file__).resolve().parents[2]
SOB_ROOT = Path(os.environ.get("SPACEOMICS_ROOT", REPO_ROOT.parent / "SpaceOmicsBench")).expanduser()
I2_CSV = SOB_ROOT / "v2_public/data/processed/cross_mission_pathway_features.csv"
OUT_DIR = Path(__file__).resolve().parent / "evaluation"

BASE_FEATURES = ["mean_NES", "std_NES", "mean_ES", "mean_padj", "min_padj",
                 "n_celltypes", "mean_size", "direction_consistency"]
LABEL_COL = "label"
MOUSE_TISSUES = ["thymus", "gastrocnemius", "skin", "eye", "liver", "kidney"]
RNG = 42
N_BOOT = 1000


def load_mouse_nes_matrix() -> pd.DataFrame:
    """Pathway × mouse-tissue mean-NES matrix (reuse tissue_nes_bridge loader)."""
    from tissue_nes_bridge import load_mouse_tissue_nes  # local import avoids top-level side effect
    long = load_mouse_tissue_nes()
    wide = long.pivot_table(index="pathway", columns="tissue", values="mouse_NES", aggfunc="mean")
    wide.columns = [f"mouse_{c}_NES" for c in wide.columns]
    return wide.reset_index()


def _cv_auroc(X: np.ndarray, y: np.ndarray, model) -> tuple[list[float], np.ndarray]:
    skf = StratifiedKFold(n_splits=5, shuffle=True, random_state=RNG)
    aurocs = []
    oof = np.zeros(len(y), dtype=float)
    for tr, te in skf.split(X, y):
        clf = Pipeline([("sc", StandardScaler()), ("clf", model)])
        clf.fit(X[tr], y[tr])
        p = clf.predict_proba(X[te])[:, 1]
        oof[te] = p
        aurocs.append(roc_auc_score(y[te], p))
    return aurocs, oof


def _bootstrap_auroc(y: np.ndarray, p: np.ndarray, n: int = N_BOOT, rng: int = RNG) -> tuple[float, float, float]:
    rs = np.random.default_rng(rng)
    aucs = []
    idx = np.arange(len(y))
    for _ in range(n):
        b = rs.choice(idx, size=len(idx), replace=True)
        if len(np.unique(y[b])) < 2:
            continue
        aucs.append(roc_auc_score(y[b], p[b]))
    aucs = np.asarray(aucs)
    return float(aucs.mean()), float(np.percentile(aucs, 2.5)), float(np.percentile(aucs, 97.5))


def _delta_ci(y: np.ndarray, p_base: np.ndarray, p_aug: np.ndarray, n: int = N_BOOT, rng: int = RNG) -> tuple[float, float, float]:
    rs = np.random.default_rng(rng)
    idx = np.arange(len(y))
    diffs = []
    for _ in range(n):
        b = rs.choice(idx, size=len(idx), replace=True)
        if len(np.unique(y[b])) < 2:
            continue
        diffs.append(roc_auc_score(y[b], p_aug[b]) - roc_auc_score(y[b], p_base[b]))
    diffs = np.asarray(diffs)
    return float(diffs.mean()), float(np.percentile(diffs, 2.5)), float(np.percentile(diffs, 97.5))


def main() -> None:
    OUT_DIR.mkdir(parents=True, exist_ok=True)

    i2 = pd.read_csv(I2_CSV)
    mouse_wide = load_mouse_nes_matrix()
    merged = i2.merge(mouse_wide, on="pathway", how="left")

    coverage = merged[[c for c in merged.columns if c.startswith("mouse_")]].notna().mean().round(3).to_dict()
    print("Per-tissue mouse coverage on I2 pathways:", coverage)

    # Impute missing mouse NES with 0 (no evidence of enrichment)
    for c in [c for c in merged.columns if c.startswith("mouse_")]:
        merged[c] = merged[c].fillna(0.0)

    y = merged["label"].values.astype(int)
    X_base = merged[BASE_FEATURES].values.astype(float)
    mouse_cols = [f"mouse_{t}_NES" for t in MOUSE_TISSUES if f"mouse_{t}_NES" in merged.columns]
    X_aug = merged[BASE_FEATURES + mouse_cols].values.astype(float)

    results = {"n_pathways": int(len(y)), "n_positive": int(y.sum()),
               "mouse_features": mouse_cols, "mouse_coverage": coverage}

    for name, model in [("logreg", LogisticRegression(max_iter=1000, class_weight="balanced")),
                        ("rf", RandomForestClassifier(n_estimators=500, random_state=RNG, n_jobs=-1, class_weight="balanced"))]:
        aurocs_b, oof_b = _cv_auroc(X_base, y, model)
        aurocs_a, oof_a = _cv_auroc(X_aug, y, model)
        boot_b = _bootstrap_auroc(y, oof_b)
        boot_a = _bootstrap_auroc(y, oof_a)
        delta = _delta_ci(y, oof_b, oof_a)

        # Per-tissue ablation: drop one mouse feature at a time
        ablation = {}
        if name == "logreg":
            for drop in mouse_cols:
                cols = BASE_FEATURES + [c for c in mouse_cols if c != drop]
                Xd = merged[cols].values.astype(float)
                aurocs_d, _ = _cv_auroc(Xd, y, LogisticRegression(max_iter=1000, class_weight="balanced"))
                ablation[drop] = {"cv_mean_auroc": float(np.mean(aurocs_d)),
                                  "delta_vs_full": float(np.mean(aurocs_d) - np.mean(aurocs_a))}

        results[name] = {
            "cv_mean_auroc_base": float(np.mean(aurocs_b)),
            "cv_std_auroc_base": float(np.std(aurocs_b)),
            "cv_mean_auroc_aug": float(np.mean(aurocs_a)),
            "cv_std_auroc_aug": float(np.std(aurocs_a)),
            "bootstrap_auroc_base": {"mean": boot_b[0], "ci_low": boot_b[1], "ci_high": boot_b[2]},
            "bootstrap_auroc_aug": {"mean": boot_a[0], "ci_low": boot_a[1], "ci_high": boot_a[2]},
            "delta_aug_minus_base": {"mean": delta[0], "ci_low": delta[1], "ci_high": delta[2]},
            "ablation_leave_one_mouse_out": ablation,
        }
        print(f"\n[{name}] base={np.mean(aurocs_b):.4f}±{np.std(aurocs_b):.4f}  "
              f"aug={np.mean(aurocs_a):.4f}±{np.std(aurocs_a):.4f}  "
              f"delta_boot={delta[0]:+.4f} [{delta[1]:+.4f}, {delta[2]:+.4f}]")
        if ablation:
            print("  Ablation (drop one mouse tissue, delta vs full-augmented):")
            for k, v in sorted(ablation.items(), key=lambda kv: kv[1]["delta_vs_full"]):
                print(f"    {k:<28} delta={v['delta_vs_full']:+.4f}")

    out_path = OUT_DIR / "supervised_conservation.json"
    with open(out_path, "w") as f:
        json.dump(results, f, indent=2)
    print(f"\nWrote {out_path}")


if __name__ == "__main__":
    main()
