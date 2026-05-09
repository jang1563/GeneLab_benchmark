"""BRIDGE supervised-conservation leakage audit.

This audit is intentionally conservative. It does not prove that the upstream
SpaceOmicsBench feature builder never observed labels, but it records the local
contract used by the v8 supervised model:

- the model feature columns are explicit and exclude the target label,
- mouse NES features are merged by pathway and are constructed without the
  supervised label column,
- fold assignments are deterministic and stratified by the binary label,
- single-feature screens do not show an obvious perfect label proxy.
"""
from __future__ import annotations

import hashlib
import json
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd
from sklearn.metrics import roc_auc_score
from sklearn.model_selection import StratifiedKFold

from supervised_conservation import (
    BASE_FEATURES,
    I2_CSV,
    LABEL_COL,
    MOUSE_TISSUES,
    N_BOOT,
    OUT_DIR,
    RNG,
    load_mouse_nes_matrix,
)


def _sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as f:
        for chunk in iter(lambda: f.read(1024 * 1024), b""):
            h.update(chunk)
    return h.hexdigest()


def _hash_values(values: list[str]) -> str:
    h = hashlib.sha256()
    for value in sorted(values):
        h.update(value.encode("utf-8"))
        h.update(b"\n")
    return h.hexdigest()


def _safe_auc(y: np.ndarray, x: pd.Series) -> float | None:
    vals = pd.to_numeric(x, errors="coerce")
    mask = vals.notna()
    if mask.sum() < 2 or len(np.unique(y[mask.to_numpy()])) < 2 or vals[mask].nunique() < 2:
        return None
    auc = roc_auc_score(y[mask.to_numpy()], vals[mask].to_numpy())
    return float(max(auc, 1.0 - auc))


def main() -> None:
    OUT_DIR.mkdir(parents=True, exist_ok=True)

    i2 = pd.read_csv(I2_CSV)
    mouse_wide = load_mouse_nes_matrix()
    merged = i2.merge(mouse_wide, on="pathway", how="left", validate="one_to_one")

    mouse_cols = [f"mouse_{t}_NES" for t in MOUSE_TISSUES if f"mouse_{t}_NES" in merged.columns]
    model_feature_cols = BASE_FEATURES + mouse_cols
    required_cols = ["pathway", LABEL_COL] + BASE_FEATURES

    missing_required = [c for c in required_cols if c not in i2.columns]
    label_like_feature_tokens = ("label", "target", "outcome", "y_true", "conserved")
    label_like_features = [
        c for c in model_feature_cols if any(token in c.lower() for token in label_like_feature_tokens)
    ]

    y = i2[LABEL_COL].astype(int).to_numpy() if LABEL_COL in i2.columns else np.array([], dtype=int)
    feature_auc_screen: dict[str, float | None] = {}
    if len(y):
        for col in model_feature_cols:
            feature_auc_screen[col] = _safe_auc(y, merged[col]) if col in merged.columns else None

    suspicious_auc_features = [
        {"feature": col, "abs_oriented_auc": auc}
        for col, auc in feature_auc_screen.items()
        if auc is not None and auc >= 0.995
    ]

    fold_manifest: list[dict[str, Any]] = []
    if len(y):
        skf = StratifiedKFold(n_splits=5, shuffle=True, random_state=RNG)
        pathways = i2["pathway"].astype(str).to_numpy()
        for fold, (train_idx, test_idx) in enumerate(skf.split(i2[BASE_FEATURES], y), start=1):
            fold_manifest.append({
                "fold": fold,
                "train_n": int(len(train_idx)),
                "test_n": int(len(test_idx)),
                "train_positive": int(y[train_idx].sum()),
                "test_positive": int(y[test_idx].sum()),
                "train_pathway_hash": _hash_values(pathways[train_idx].tolist()),
                "test_pathway_hash": _hash_values(pathways[test_idx].tolist()),
            })

    warnings = []
    if missing_required:
        warnings.append(f"Missing required I2 columns: {missing_required}")
    if label_like_features:
        warnings.append(f"Label-like model feature names found: {label_like_features}")
    if suspicious_auc_features:
        warnings.append("One or more single features nearly perfectly reproduce the label.")
    if i2["pathway"].duplicated().any():
        warnings.append("I2 pathway identifiers are not unique.")
    if mouse_wide["pathway"].duplicated().any():
        warnings.append("Mouse NES pathway identifiers are not unique after aggregation.")

    audit = {
        "schema_version": "0.1.0",
        "audit_id": "bridge.supervised_conservation.leakage_audit",
        "status": "pass" if not warnings else "review",
        "inputs": {
            "i2_features": {
                "path": str(I2_CSV),
                "sha256": _sha256(I2_CSV),
                "n_rows": int(i2.shape[0]),
                "columns": list(i2.columns),
            },
            "mouse_nes_matrix": {
                "n_rows": int(mouse_wide.shape[0]),
                "columns": list(mouse_wide.columns),
            },
        },
        "model_contract": {
            "label_column": LABEL_COL,
            "base_features": BASE_FEATURES,
            "mouse_features": mouse_cols,
            "n_model_features": len(model_feature_cols),
            "label_excluded_from_features": LABEL_COL not in model_feature_cols,
            "missing_required_columns": missing_required,
            "bootstrap": N_BOOT,
            "cv": {
                "class": "StratifiedKFold",
                "n_splits": 5,
                "shuffle": True,
                "random_state": RNG,
            },
        },
        "merge_checks": {
            "merge_key": "pathway",
            "merged_rows": int(merged.shape[0]),
            "i2_pathway_unique": bool(not i2["pathway"].duplicated().any()),
            "mouse_pathway_unique_after_aggregation": bool(not mouse_wide["pathway"].duplicated().any()),
            "mouse_feature_coverage": merged[mouse_cols].notna().mean().round(3).to_dict(),
            "mouse_missing_imputation": "0.0 in supervised_conservation.py",
        },
        "fold_manifest": fold_manifest,
        "single_feature_label_auc_screen": feature_auc_screen,
        "suspicious_auc_features": suspicious_auc_features,
        "warnings": warnings,
        "limitations": [
            "This audit verifies the v8 local modeling contract, not the complete upstream SpaceOmicsBench feature-generation pipeline.",
            "The target label is loaded from the same CSV as baseline features but is not selected into X_base or X_aug.",
            "A full beta freeze should archive the upstream SpaceOmicsBench commit or release tag and feature-builder command.",
        ],
    }

    out_path = OUT_DIR / "bridge_leakage_audit.json"
    with out_path.open("w") as f:
        json.dump(audit, f, indent=2)
    print(json.dumps({"audit": str(out_path), "status": audit["status"], "warnings": warnings}, indent=2))


if __name__ == "__main__":
    main()
