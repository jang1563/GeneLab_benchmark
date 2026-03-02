#!/usr/bin/env python3
"""
shap_analysis.py — GeneLab_benchmark: SHAP Feature Importance (DD-11 Condition 3)

Trains the best baseline model on each LOMO fold and computes SHAP values to
identify biologically relevant features driving spaceflight classification.

Phase 1 Condition 3 (DD-11):
  - Verify that the model learns known spaceflight-responsive genes
  - Report top-N genes by mean |SHAP| across all test samples
  - Check whether tissue-specific target genes appear in top-50

Known target genes by tissue:
  liver:   ANGPTL4, PCK1, G6pc, Fasn, Ppara, Acox1, Cyp7a1
  thymus:  Lef1, Ikzf2, Rag1, Rag2, Il7r, Ccr7, Foxn1, Myc, Trp53

Usage:
  python shap_analysis.py --task A4           # SHAP for A4 thymus
  python shap_analysis.py --task A1           # SHAP for A1 liver
  python shap_analysis.py --task A4 --model rf  # RF SHAP only
  python shap_analysis.py --task A4 --top-n 100  # show top 100 genes
"""

import json
import argparse
import warnings
import numpy as np
import pandas as pd
from pathlib import Path

warnings.filterwarnings("ignore")

BASE_DIR = Path(__file__).resolve().parent.parent
TASKS_DIR = BASE_DIR / "tasks"
RESULTS_DIR = BASE_DIR / "evaluation"
PROCESSED_DIR = BASE_DIR / "processed"


def load_symbol_map() -> dict:
    """Load ENSEMBL → gene symbol map built from DE files."""
    map_file = PROCESSED_DIR / "ensembl_symbol_map.csv"
    if not map_file.exists():
        return {}
    df = pd.read_csv(map_file)
    return dict(zip(df["ENSEMBL"], df["SYMBOL"]))

# ── Tissue-specific target genes for condition 3 check ────────────────────────

TARGET_GENES = {
    "liver": [
        "Angptl4", "ANGPTL4", "angptl4",
        "Pck1", "PCK1", "pck1",
        "G6pc", "G6PC",
        "Fasn", "FASN",
        "Ppara", "PPARA",
        "Acox1", "ACOX1",
        "Cyp7a1", "CYP7A1",
        "Hmgcr", "HMGCR",
        "Scd1", "SCD1",
    ],
    "thymus": [
        # T-cell development & thymic function
        "Lef1",      # Wnt/TCF target, T-cell development
        "Ikzf2",     # Helios, Treg/thymic T-cell TF
        "Rag1",      # V(D)J recombination
        "Rag2",      # V(D)J recombination
        "Il7r",      # IL-7 receptor, thymocyte survival
        "Ccr7",      # mature T-cell homing
        "Foxn1",     # thymic epithelium master TF
        "Dntt",      # TdT, immature thymocyte marker
        "Ptcra",     # pre-TCR alpha
        "Cd4",       # CD4+ T-cell marker
        "Cd8a",      # CD8+ T-cell marker
        "Bcl2",      # anti-apoptotic, thymic selection
        "Maf",       # c-Maf, T-cell TF (Treg, Tfh, γδT)
        "Btla",      # T-cell inhibitory receptor
        "Gzma",      # Granzyme A, cytotoxic T-cell effector
        # DNA damage response (spaceflight radiation)
        "Brca1",     # DNA DSB repair
        "Msh6",      # DNA mismatch repair
        # TCR genes
        "Trbv12-1",  # T-cell receptor beta variable
        # Chromatin / epigenetic
        "H2ac12",    # Histone H2A variant
        "H2ac22",    # Histone H2A variant
        # Known spaceflight genes
        "Apoe",      # Apolipoprotein E, spaceflight lipid metabolism
        "Apod",      # Apolipoprotein D
        "Trp63",     # p63, thymic epithelial marker
    ],
    "gastrocnemius": [
        # Muscle atrophy (spaceflight-specific)
        "Fbxo32",    # Atrogin-1, E3 ubiquitin ligase, muscle atrophy marker
        "Trim63",    # MuRF-1, E3 ubiquitin ligase, muscle atrophy marker
        "Mstn",      # Myostatin, muscle growth inhibitor
        "Foxo3",     # FOXO3, atrophy/autophagy TF
        "Igf1",      # IGF-1, muscle growth/atrophy signaling
        "Klf15",     # KLF15, circadian muscle TF, atrophy
        "Myog",      # Myogenin, myogenesis TF
        "Myh7",      # Myosin heavy chain beta (slow-twitch)
        "Myh4",      # Myosin heavy chain IIb (fast-twitch)
        # Circadian / metabolism (consistently top in spaceflight muscle)
        "Bmal1",     # Master circadian TF (Arntl)
        "Per2",      # Period 2, circadian repressor
        "Dbp",       # D-box binding protein, circadian output
        "Ciart",     # Circadian-associated repressor
        # Oxidative stress / mitochondria
        "Hmox1",     # Heme oxygenase-1, oxidative stress
        "Sod2",      # MnSOD, mitochondrial antioxidant
        "Ppargc1a",  # PGC-1α, mitochondrial biogenesis
    ],
    "kidney": [
        "Havcr1", "HAVCR1",     # KIM-1, injury marker
        "Lcn2", "LCN2",         # NGAL
        "Vim", "VIM",
        "Fn1", "FN1",
        "Tgfb1", "TGFB1",
        "Acta2", "ACTA2",
    ],
    "eye": [
        "Rho", "RHO",
        "Pde6b", "PDE6B",
        "Rdh12", "RDH12",
        "Vegfa", "VEGFA",
        "Hif1a", "HIF1A",
        "Sod2", "SOD2",
        "Txnip", "TXNIP",
    ],
}

# Task → tissue mapping
TASK_TISSUE = {
    "A1": "liver",
    "A2": "gastrocnemius",
    "A3": "kidney",
    "A4": "thymus",
    "A5": "skin",
    "A6": "eye",
}


# ── Model builders ─────────────────────────────────────────────────────────────

def build_rf(seed: int = 42):
    from sklearn.ensemble import RandomForestClassifier
    return RandomForestClassifier(
        n_estimators=200,
        max_features="sqrt",
        class_weight="balanced",
        n_jobs=-1,
        random_state=seed,
    )


def build_lr(seed: int = 42):
    from sklearn.linear_model import LogisticRegression
    from sklearn.preprocessing import StandardScaler
    from sklearn.pipeline import Pipeline
    return Pipeline([
        ("scaler", StandardScaler()),
        ("clf", LogisticRegression(
            penalty="elasticnet",
            solver="saga",
            l1_ratio=0.5,
            C=1.0,
            class_weight="balanced",
            max_iter=2000,
            random_state=seed,
        ))
    ])


# ── SHAP computation ───────────────────────────────────────────────────────────

def compute_shap_rf(model, X_train: np.ndarray, X_test: np.ndarray,
                    gene_names: list) -> pd.DataFrame:
    """
    SHAP TreeExplainer for Random Forest.
    Returns DataFrame: genes × folds, values = mean |SHAP| on test set.
    """
    import shap
    explainer = shap.TreeExplainer(model)
    shap_values = explainer.shap_values(X_test)

    # For binary classification, shap_values is [class0, class1] or (n_samples, n_features, 2)
    # We want the class-1 (Flight) SHAP values
    if isinstance(shap_values, list):
        sv = shap_values[1]           # class 1 = Flight
    elif shap_values.ndim == 3:
        sv = shap_values[:, :, 1]
    else:
        sv = shap_values

    mean_abs_shap = np.mean(np.abs(sv), axis=0)
    return pd.Series(mean_abs_shap, index=gene_names, name="mean_abs_shap")


def compute_shap_lr(model, X_train: np.ndarray, X_test: np.ndarray,
                    gene_names: list) -> pd.Series:
    """
    SHAP LinearExplainer for Logistic Regression (post-StandardScaler).
    Approximates feature importance via coefficient × std(feature).
    Returns Series: genes, values = |coef| × std.
    """
    from sklearn.pipeline import Pipeline
    from sklearn.linear_model import LogisticRegression as LR

    if isinstance(model, Pipeline):
        scaler = model.named_steps.get("scaler")
        clf = model.named_steps.get("clf")
        X_scaled = scaler.transform(X_test)
    else:
        clf = model
        X_scaled = X_test

    coef = clf.coef_[0]  # (n_features,) for binary
    mean_abs_shap = np.abs(coef) * np.std(X_scaled, axis=0)
    return pd.Series(mean_abs_shap, index=gene_names, name="mean_abs_shap")


# ── Per-fold SHAP ──────────────────────────────────────────────────────────────

def shap_fold(fold_dir: Path, model_name: str = "rf") -> pd.Series | None:
    """Train model on fold, compute SHAP, return mean |SHAP| Series."""
    from sklearn.metrics import roc_auc_score

    train_X = pd.read_csv(fold_dir / "train_X.csv", index_col=0)
    test_X  = pd.read_csv(fold_dir / "test_X.csv",  index_col=0)
    train_y = pd.read_csv(fold_dir / "train_y.csv", index_col=0).squeeze()
    test_y  = pd.read_csv(fold_dir / "test_y.csv",  index_col=0).squeeze()

    fold_info = json.loads((fold_dir / "fold_info.json").read_text()) if (fold_dir / "fold_info.json").exists() else {}
    test_mission = fold_info.get("test_mission", fold_dir.name)

    common_genes = train_X.columns.intersection(test_X.columns)
    train_X = train_X[common_genes]
    test_X  = test_X[common_genes]

    X_train = train_X.values.astype(np.float32)
    X_test  = test_X.values.astype(np.float32)
    y_train = train_y.values.astype(int)
    y_test  = test_y.values.astype(int)

    if len(np.unique(y_test)) < 2:
        print(f"    [SKIP] {test_mission}: only one class in test")
        return None

    if model_name == "rf":
        model = build_rf()
    elif model_name == "lr":
        model = build_lr()
    else:
        raise ValueError(f"Unknown model: {model_name}")

    print(f"  Training {model_name} on fold {test_mission} "
          f"(n_train={len(y_train)}, n_test={len(y_test)})...")

    try:
        model.fit(X_train, y_train)
    except Exception as e:
        print(f"    [ERROR] Training failed: {e}")
        return None

    # AUROC sanity check
    if hasattr(model, "predict_proba"):
        y_score = model.predict_proba(X_test)[:, 1]
    else:
        y_score = model.decision_function(X_test)
    auroc = roc_auc_score(y_test, y_score)
    print(f"    AUROC = {auroc:.3f}")

    # Compute SHAP
    if model_name == "rf":
        shap_s = compute_shap_rf(model, X_train, X_test, list(common_genes))
    else:
        shap_s = compute_shap_lr(model, X_train, X_test, list(common_genes))

    shap_s.name = test_mission
    return shap_s


# ── Task-level SHAP aggregation ────────────────────────────────────────────────

def shap_task(task_id: str, model_name: str = "rf", top_n: int = 50) -> dict:
    """Compute SHAP for all folds, aggregate, return top-N genes."""

    task_dirs = list(TASKS_DIR.glob(f"{task_id}_*_lomo"))
    if not task_dirs:
        print(f"[ERROR] No task directory found for {task_id}")
        return {}

    task_dir = task_dirs[0]
    tissue = TASK_TISSUE.get(task_id, "unknown")
    target_genes = TARGET_GENES.get(tissue, [])

    fold_dirs = sorted([d for d in task_dir.iterdir()
                        if d.is_dir() and d.name.startswith("fold_")])

    print(f"\n{'='*60}")
    print(f"SHAP Analysis — Task {task_id} ({tissue})")
    print(f"Model: {model_name.upper()}, Top-{top_n} genes")
    print(f"Folds: {len(fold_dirs)}")
    print(f"{'='*60}")

    fold_shaps = []
    for fold_dir in fold_dirs:
        s = shap_fold(fold_dir, model_name=model_name)
        if s is not None:
            fold_shaps.append(s)

    if not fold_shaps:
        print("[ERROR] No SHAP results computed.")
        return {}

    # Aggregate: mean |SHAP| across folds (on common genes)
    shap_df = pd.concat(fold_shaps, axis=1)   # genes × folds
    shap_df = shap_df.fillna(0.0)
    mean_shap = shap_df.mean(axis=1).sort_values(ascending=False)

    top_genes = mean_shap.head(top_n)

    # Load Ensembl → Symbol map (may be empty for non-Ensembl feature IDs)
    sym_map = load_symbol_map()

    def resolve_symbol(gene_id: str) -> str:
        """Return symbol if Ensembl ID, else return gene_id as-is."""
        return sym_map.get(gene_id, gene_id)

    print(f"\n{'='*60}")
    print(f"Top-{top_n} genes by mean |SHAP| (across {len(fold_shaps)} folds)")
    print(f"{'='*60}")
    print(f"{'Rank':>4}  {'Symbol':20}  {'ID':25}  {'Mean|SHAP|':>10}")
    print("-" * 68)
    for rank, (gene, val) in enumerate(top_genes.items(), 1):
        sym = resolve_symbol(gene)
        marker = ""
        sym_lower = sym.lower()
        for tg in target_genes:
            if sym_lower == tg.lower():
                marker = "  ← TARGET"
                break
        print(f"{rank:>4}  {sym:20}  {gene:25}  {val:>10.6f}{marker}")

    # Build symbol-indexed mean_shap for target gene lookup
    mean_shap_sym = mean_shap.copy()
    mean_shap_sym.index = [resolve_symbol(g) for g in mean_shap.index]

    # Check target genes in top-N (by symbol)
    print(f"\n{'='*60}")
    print(f"Target Gene Check (tissue: {tissue})")
    print(f"{'='*60}")

    # Deduplicate target genes (remove same-symbol duplicates like Lef1/LEF1)
    seen_symbols = set()
    unique_targets = []
    for tg in target_genes:
        if tg.lower() not in seen_symbols:
            seen_symbols.add(tg.lower())
            unique_targets.append(tg)

    found = []
    not_found = []
    for tg in unique_targets:
        # Find all matching ranks by symbol (case-insensitive)
        ranks = [i + 1 for i, sym in enumerate(mean_shap_sym.index)
                 if sym.lower() == tg.lower()]
        if ranks:
            r = ranks[0]
            in_top = r <= top_n
            found.append((tg, r, in_top))
            status = f"rank {r}" + (f" ✓ (top-{top_n})" if in_top else f" (outside top-{top_n})")
        else:
            not_found.append(tg)
            status = "not in features"
        print(f"  {tg:20}: {status}")

    n_in_top = sum(1 for _, _, in_top in found if in_top)
    print(f"\n  In top-{top_n}: {n_in_top}/{len(unique_targets)} target genes")
    cond3_pass = n_in_top >= 1
    print(f"  Condition 3: {'✓ PASS' if cond3_pass else '✗ FAIL'} "
          f"(≥1 target gene in top-{top_n})")

    # Save results (use symbols in output)
    out_data = {
        "task_id": task_id,
        "tissue": tissue,
        "model": model_name,
        "top_n": top_n,
        "n_folds": len(fold_shaps),
        "top_genes": {
            resolve_symbol(gene): {"ensembl": gene, "mean_abs_shap": float(val)}
            for gene, val in top_genes.items()
        },
        "target_gene_ranks": {
            tg: {"rank": r, "in_top_n": in_top}
            for tg, r, in_top in found
        },
        "condition3_pass": cond3_pass,
    }

    out_path = RESULTS_DIR / f"{task_id}_shap_{model_name}.json"
    RESULTS_DIR.mkdir(parents=True, exist_ok=True)
    out_path.write_text(json.dumps(out_data, indent=2))
    print(f"\n  ✓ SHAP results saved: {out_path}")

    return out_data


# ── CLI ────────────────────────────────────────────────────────────────────────

def main():
    parser = argparse.ArgumentParser(
        description="SHAP feature importance analysis for GeneLab_benchmark tasks"
    )
    parser.add_argument("--task", type=str, default="A4")
    parser.add_argument("--model", type=str, default="rf", choices=["rf", "lr"])
    parser.add_argument("--top-n", type=int, default=50)
    args = parser.parse_args()

    shap_task(
        task_id=args.task.upper(),
        model_name=args.model,
        top_n=args.top_n,
    )


if __name__ == "__main__":
    main()
