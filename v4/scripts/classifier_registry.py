#!/usr/bin/env python3
"""
classifier_registry.py — GeneLabBench v4: Standardized classifier builder API

8 classifiers spanning 5 ML families (DD-31):
  Linear:         ElasticNet-LR (existing), PCA-LR (existing)
  Tree-based:     Random Forest (existing), XGBoost (existing)
  Kernel:         SVM-RBF (new)
  Instance-based: kNN (new)
  Neural:         MLP (new), TabNet (new, GPU)

All classifiers expose fit(X, y) and predict_proba(X) interfaces.
No hyperparameter tuning (DD-22): established defaults only.
"""

import numpy as np


# ── Existing classifiers (from scripts/run_baselines.py) ──────────────────────

def build_elasticnet_lr(seed=42):
    """ElasticNet Logistic Regression (L1+L2, SAGA solver).
    Existing v1 'lr' model. L1 component induces sparsity on gene space.
    """
    from sklearn.linear_model import LogisticRegression
    from sklearn.preprocessing import StandardScaler
    from sklearn.pipeline import Pipeline
    import sklearn
    lr_kwargs = dict(
        solver="saga",
        l1_ratio=0.5,
        C=1.0,
        class_weight="balanced",
        max_iter=10000,
        random_state=seed,
    )
    # sklearn ≥1.8 deprecated 'penalty' param; use l1_ratio directly
    if tuple(int(x) for x in sklearn.__version__.split(".")[:2]) < (1, 8):
        lr_kwargs["penalty"] = "elasticnet"
    return Pipeline([
        ("scaler", StandardScaler()),
        ("clf", LogisticRegression(**lr_kwargs))
    ])


def build_pca_lr(n_components=50, seed=42):
    """PCA + L2 Logistic Regression.
    Existing v1 baseline. n_components adapted to min(50, n_train-1) at fit time.
    """
    from sklearn.decomposition import PCA
    from sklearn.linear_model import LogisticRegression
    from sklearn.preprocessing import StandardScaler
    from sklearn.pipeline import Pipeline
    return Pipeline([
        ("scaler", StandardScaler()),
        ("pca", PCA(n_components=n_components, random_state=seed)),
        ("clf", LogisticRegression(
            C=1.0,
            class_weight="balanced",
            max_iter=5000,
            random_state=seed,
        ))
    ])


def build_rf(seed=42):
    """Random Forest. Existing v1 model. TreeSHAP-compatible."""
    from sklearn.ensemble import RandomForestClassifier
    return RandomForestClassifier(
        n_estimators=200,
        max_features="sqrt",
        class_weight="balanced",
        n_jobs=-1,
        random_state=seed,
    )


def build_xgb(seed=42):
    """XGBoost. Existing v1 model. TreeSHAP-compatible."""
    from xgboost import XGBClassifier
    return XGBClassifier(
        n_estimators=200,
        max_depth=3,
        subsample=0.7,
        learning_rate=0.1,
        scale_pos_weight=1,
        eval_metric="logloss",
        verbosity=0,
        random_state=seed,
    )


# ── New classifiers (v4) ─────────────────────────────────────────────────────

def build_svm_rbf(seed=42):
    """SVM with RBF kernel. probability=True enables predict_proba via Platt scaling.
    KernelSHAP required for interpretability.
    """
    from sklearn.svm import SVC
    from sklearn.preprocessing import StandardScaler
    from sklearn.pipeline import Pipeline
    return Pipeline([
        ("scaler", StandardScaler()),
        ("clf", SVC(
            kernel="rbf",
            C=1.0,
            gamma="scale",
            probability=True,
            class_weight="balanced",
            random_state=seed,
        ))
    ])


def build_knn(seed=42):
    """k-Nearest Neighbors (k=5, distance-weighted).
    Instance-based: no model fitting, prediction based on training data proximity.
    seed unused but accepted for API consistency.
    """
    from sklearn.neighbors import KNeighborsClassifier
    from sklearn.preprocessing import StandardScaler
    from sklearn.pipeline import Pipeline
    return Pipeline([
        ("scaler", StandardScaler()),
        ("clf", KNeighborsClassifier(
            n_neighbors=5,
            weights="distance",
            metric="euclidean",
            n_jobs=-1,
        ))
    ])


def build_mlp(seed=42):
    """Multi-Layer Perceptron (sklearn). Two hidden layers (256, 64).
    early_stopping=True with 10% validation split prevents overfitting on small n.
    """
    from sklearn.neural_network import MLPClassifier
    from sklearn.preprocessing import StandardScaler
    from sklearn.pipeline import Pipeline
    return Pipeline([
        ("scaler", StandardScaler()),
        ("clf", MLPClassifier(
            hidden_layer_sizes=(256, 64),
            activation="relu",
            solver="adam",
            alpha=1e-4,
            learning_rate="adaptive",
            learning_rate_init=1e-3,
            early_stopping=True,
            validation_fraction=0.1,
            n_iter_no_change=10,
            max_iter=500,
            random_state=seed,
        ))
    ])


def build_tabnet(seed=42):
    """TabNet (pytorch-tabnet). Attention-based deep learning for tabular data.
    Requires GPU for reasonable performance. Falls back to CPU if unavailable.
    """
    from pytorch_tabnet.tab_model import TabNetClassifier
    import torch

    device = "cuda" if torch.cuda.is_available() else "cpu"
    return TabNetClassifier(
        n_d=8,
        n_a=8,
        n_steps=3,
        gamma=1.3,
        lambda_sparse=1e-3,
        optimizer_fn=torch.optim.Adam,
        optimizer_params=dict(lr=2e-2),
        scheduler_fn=torch.optim.lr_scheduler.StepLR,
        scheduler_params=dict(step_size=10, gamma=0.9),
        mask_type="entmax",
        seed=seed,
        verbose=0,
        device_name=device,
    )


# ── Registry ──────────────────────────────────────────────────────────────────

CLASSIFIERS = {
    # Existing (from run_baselines.py)
    "elasticnet_lr": ("ElasticNet-LR",  build_elasticnet_lr, "cpu"),
    "pca_lr":        ("PCA-LR",         build_pca_lr,        "cpu"),
    "rf":            ("Random Forest",   build_rf,            "cpu"),
    "xgb":           ("XGBoost",         build_xgb,           "cpu"),
    # New (v4)
    "svm_rbf":       ("SVM-RBF",         build_svm_rbf,       "cpu"),
    "knn":           ("kNN",             build_knn,           "cpu"),
    "mlp":           ("MLP",             build_mlp,           "cpu"),
    "tabnet":        ("TabNet",          build_tabnet,        "gpu"),
}

# Method family mapping (for SHAP strategy, DD-23)
SHAP_STRATEGY = {
    "elasticnet_lr": "linear",     # |coef| × std
    "pca_lr":        "pca_linear", # |coef| × PCA_components back to gene space
    "rf":            "tree",       # TreeExplainer
    "xgb":           "tree",       # TreeExplainer
    "svm_rbf":       "kernel",     # KernelExplainer(kmeans=100)
    "knn":           "kernel",     # KernelExplainer(kmeans=100)
    "mlp":           "kernel",     # KernelExplainer(kmeans=100)
    "tabnet":        "kernel",     # KernelExplainer(kmeans=100)
}

CPU_METHODS = [k for k, v in CLASSIFIERS.items() if v[2] == "cpu"]
GPU_METHODS = [k for k, v in CLASSIFIERS.items() if v[2] == "gpu"]


def get_classifier(method_key, seed=42):
    """Build and return a fresh classifier instance.

    Args:
        method_key: Key from CLASSIFIERS dict
        seed: Random state for reproducibility

    Returns:
        (label, model) tuple
    """
    if method_key not in CLASSIFIERS:
        raise ValueError(f"Unknown method: {method_key}. Available: {list(CLASSIFIERS.keys())}")

    label, builder, _ = CLASSIFIERS[method_key]
    model = builder(seed=seed)
    return label, model


def adapt_pca_components(model, n_train, n_features):
    """Adapt PCA n_components for small training sets (avoids n_components > n_train-1).
    Modifies model in-place.
    """
    from sklearn.pipeline import Pipeline
    from sklearn.decomposition import PCA

    if not isinstance(model, Pipeline):
        return

    for _, step in model.steps:
        if isinstance(step, PCA):
            max_comps = min(step.n_components, n_train - 1, n_features)
            if max_comps < step.n_components:
                step.n_components = max_comps


def fit_tabnet_with_eval(model, X_train, y_train, X_val=None, y_val=None,
                         max_epochs=200, patience=20, batch_size=256):
    """Fit TabNet with optional validation set for early stopping.
    TabNet API differs from sklearn — requires special handling.
    """
    fit_params = {
        "X_train": X_train,
        "y_train": y_train,
        "max_epochs": max_epochs,
        "patience": patience,
        "batch_size": batch_size,
    }
    if X_val is not None and y_val is not None:
        fit_params["eval_set"] = [(X_val, y_val)]
        fit_params["eval_name"] = ["val"]
        fit_params["eval_metric"] = ["auc"]

    model.fit(**fit_params)
    return model


# ── Verification ──────────────────────────────────────────────────────────────

if __name__ == "__main__":
    """Quick verification: all 8 classifiers run on dummy data."""
    from sklearn.datasets import make_classification

    X, y = make_classification(n_samples=100, n_features=100, n_informative=20,
                               random_state=42)

    print("Classifier Registry Verification")
    print("=" * 60)

    for key in CLASSIFIERS:
        try:
            label, model = get_classifier(key, seed=42)
            if key == "tabnet":
                # TabNet needs float32 and special fit
                X_f = X.astype(np.float32)
                y_i = y.astype(np.int64)
                fit_tabnet_with_eval(model, X_f[:80], y_i[:80],
                                     X_f[80:], y_i[80:],
                                     max_epochs=10, patience=5)
                proba = model.predict_proba(X_f[80:])
            else:
                model.fit(X[:80], y[:80])
                proba = model.predict_proba(X[80:])

            auroc_approx = np.mean(proba[:, 1][y[80:] == 1]) > np.mean(proba[:, 1][y[80:] == 0])
            print(f"  [{key:15s}] {label:20s} — OK (predict_proba shape: {proba.shape})")
        except Exception as e:
            print(f"  [{key:15s}] {label:20s} — FAILED: {e}")

    print(f"\nTotal: {len(CLASSIFIERS)} classifiers")
    print(f"CPU: {CPU_METHODS}")
    print(f"GPU: {GPU_METHODS}")
    print(f"SHAP strategies: {set(SHAP_STRATEGY.values())}")
