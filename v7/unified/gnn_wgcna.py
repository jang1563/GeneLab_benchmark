"""
GNN on WGCNA Topology — v7 GeneLabBench
========================================
Graph Neural Network classifier using WGCNA topological overlap matrix (TOM)
as the gene-gene graph structure for spaceflight classification.

Each bulk RNA-seq sample -> one graph:
  - Nodes  = top-N MAD genes (matching WGCNA gene count)
  - Edges  = WGCNA TOM, random graph, or no edges
  - Node feature = z-scored expression value of gene in that sample (1-dim)
  - Graph label = 1 (Flight) or 0 (Ground)

Three graph conditions compared:
  1. WGCNA-TOM guided  (biological co-expression topology)
  2. Random graph      (same edge budget, random connectivity)
  3. No-edges MLP      (global mean pooling baseline)

Evaluation: LOMO-CV (same folds as v4 Phase 1)
Comparison: PCA-LR baseline from v4/evaluation/M1_summary.json

Method note:
  To avoid cross-fold leakage, v7 defaults to train-fold feature selection,
  train-fold topology construction, and train-fold scaling.

Usage:
  python gnn_wgcna.py --tissue liver --graph-type wgcna
  python gnn_wgcna.py --tissue all --graph-type all

Output:
  v7/evaluation/GNN_{tissue}_{graph_type}.json
  v7/evaluation/GNN_summary.json
"""

import argparse
import json
import random
import time
import warnings
from pathlib import Path

import numpy as np
import pandas as pd
from sklearn.metrics import roc_auc_score
from sklearn.preprocessing import StandardScaler

warnings.filterwarnings("ignore")

BASE_DIR = Path(__file__).resolve().parent.parent.parent
V4_EVAL = BASE_DIR / "v4" / "evaluation"
PROC_DIR = BASE_DIR / "processed" / "A_detection"
M1_SUMMARY = V4_EVAL / "M1_summary.json"
OUT_DIR = BASE_DIR / "v7" / "evaluation"
OUT_DIR.mkdir(parents=True, exist_ok=True)

TISSUE_MISSIONS = {
    "liver": ["RR-1", "RR-3", "RR-6", "RR-8", "RR-9", "MHU-2"],
    "gastrocnemius": ["RR-1", "RR-5", "RR-9"],
    "kidney": ["RR-1", "RR-3", "RR-7"],
    "thymus": ["MHU-1", "MHU-2", "RR-6", "RR-9"],
    "eye": ["RR-1", "RR-3", "TBD"],
    "skin": ["MHU-2", "RR-6", "RR-7"],
}

EXPR_FILE = {
    "liver": "liver_all_missions_log2_norm_limma_rbe.csv",
    "gastrocnemius": "gastrocnemius_all_missions_log2_norm_limma_rbe.csv",
    "kidney": "kidney_all_missions_log2_norm_limma_rbe.csv",
    "thymus": "thymus_all_missions_log2_norm_limma_rbe.csv",
    "eye": "eye_all_missions_log2_norm_limma_rbe.csv",
    "skin": "skin_all_missions_log2_norm.csv",
}

META_FILE = {t: f"{t}_all_missions_metadata.csv" for t in TISSUE_MISSIONS}

LABEL_MAP = {
    "FLT": 1,
    "Flight": 1,
    "GC": 0,
    "Ground Control": 0,
    "Ground": 0,
    "Basal": 0,
    "BC": 0,
    "VC": 0,
    "Vivarium": 0,
}
EXCLUDE_LABELS = {"AG", "Artificial Gravity"}


def load_tissue_data(tissue):
    """Load expression + metadata, align samples, and encode labels."""
    tdir = PROC_DIR / tissue
    expr = pd.read_csv(tdir / EXPR_FILE[tissue], index_col=0)
    meta = pd.read_csv(tdir / META_FILE[tissue])

    # Align metadata index with expression index
    if "sample_name" in meta.columns:
        meta = meta.set_index("sample_name")
    elif "Unnamed: 0" in meta.columns:
        meta = meta.set_index("Unnamed: 0")
    common = expr.index.intersection(meta.index)
    expr = expr.loc[common]
    meta = meta.loc[common]

    # Find label column
    for candidate in ["label", "condition", "Label", "Condition"]:
        if candidate in meta.columns:
            label_col = candidate
            break
    else:
        label_col = meta.columns[0]
    y_raw = meta[label_col].astype(str).str.strip()
    valid = y_raw.apply(lambda value: value not in EXCLUDE_LABELS and value in LABEL_MAP)
    expr = expr.loc[valid]
    meta = meta.loc[valid]
    y = y_raw[valid].map(LABEL_MAP).values.astype(np.int32)

    if "mission" not in meta.columns:
        meta = meta.copy()
        meta["mission"] = "unknown"

    return expr, y, meta


def filter_top_mad_genes(expr_df, n_genes):
    """Select the top-N MAD genes from the provided dataframe."""
    mad = np.median(np.abs(expr_df.values - np.median(expr_df.values, axis=0)), axis=0)
    top_idx = np.argsort(mad)[-n_genes:]
    return expr_df.iloc[:, top_idx], top_idx


def compute_tom_from_expr(x_genes_samples, beta):
    """Compute a TOM matrix from a genes x samples expression matrix."""
    print(f"  Correlation matrix ({x_genes_samples.shape[0]}x{x_genes_samples.shape[0]})...")
    x = x_genes_samples.astype(np.float32)
    mu = x.mean(axis=1, keepdims=True)
    sd = x.std(axis=1, keepdims=True)
    sd[sd < 1e-8] = 1e-8
    x_std = (x - mu) / sd

    n_samples = x.shape[1]
    cor = (x_std @ x_std.T) / max(n_samples - 1, 1)
    np.fill_diagonal(cor, 1.0)
    np.clip(cor, -1.0, 1.0, out=cor)

    print(f"  Adjacency (beta={beta})...")
    adj = np.power(np.clip(0.5 + 0.5 * cor, 0.0, 1.0), beta, dtype=np.float32)
    np.fill_diagonal(adj, 0.0)
    del cor

    print("  TOM (matrix multiply)...")
    k = adj.sum(axis=0)
    aa = (adj @ adj).astype(np.float32)
    numerator = aa + adj
    min_k = np.minimum(k[:, None], k[None, :]).astype(np.float32)
    denominator = np.maximum(min_k + 1.0 - adj, 1e-8)
    tom = numerator / denominator
    np.fill_diagonal(tom, 0.0)
    np.clip(tom, 0.0, 1.0, out=tom)
    return tom


def tom_to_edge_list(tom, n_edges_per_gene, n_genes, rng=None, graph_type="wgcna"):
    """Convert TOM (or a null graph spec) into a bidirectional edge list."""
    if graph_type == "no_edges":
        return np.array([], dtype=np.int32), np.array([], dtype=np.int32)

    if rng is None:
        rng = np.random.default_rng(42)

    if graph_type == "random":
        n_edges_directed = n_genes * n_edges_per_gene
        src = rng.integers(0, n_genes, size=n_edges_directed)
        dst = rng.integers(0, n_genes, size=n_edges_directed)
        mask = src != dst
        src, dst = src[mask], dst[mask]
        return (
            np.concatenate([src, dst]).astype(np.int32),
            np.concatenate([dst, src]).astype(np.int32),
        )

    src_list, dst_list = [], []
    top_k = min(n_edges_per_gene, n_genes - 1)
    for gene_idx in range(n_genes):
        row = tom[gene_idx].copy()
        row[gene_idx] = 0.0
        top_j = np.argpartition(row, -top_k)[-top_k:]
        src_list.append(np.full(top_k, gene_idx, dtype=np.int32))
        dst_list.append(top_j.astype(np.int32))

    src = np.concatenate(src_list)
    dst = np.concatenate(dst_list)
    return (
        np.concatenate([src, dst]).astype(np.int32),
        np.concatenate([dst, src]).astype(np.int32),
    )


def build_pyg_dataset(x_norm, y, edge_src, edge_dst):
    """Build a PyG dataset with one graph per sample."""
    import torch
    from torch_geometric.data import Data

    dataset = []
    edge_index = torch.tensor(np.stack([edge_src, edge_dst], axis=0), dtype=torch.long)
    for sample_idx in range(x_norm.shape[0]):
        node_x = torch.tensor(x_norm[sample_idx, :, None], dtype=torch.float32)
        label = torch.tensor([y[sample_idx]], dtype=torch.float32)
        dataset.append(Data(x=node_x, edge_index=edge_index, y=label))
    return dataset


def seed_torch(seed):
    """Make model init and shuffled batches reproducible per fold."""
    import torch

    random.seed(seed)
    np.random.seed(seed)
    torch.manual_seed(seed)
    if torch.cuda.is_available():
        torch.cuda.manual_seed_all(seed)


def build_gnn_model(in_channels, hidden_channels, use_attention=True):
    """Build a graph-level classifier."""
    import torch.nn as nn
    import torch.nn.functional as F
    from torch_geometric.nn import GATConv, global_mean_pool

    class GATClassifier(nn.Module):
        def __init__(self):
            super().__init__()
            self.conv1 = GATConv(in_channels, hidden_channels, heads=4, concat=False, dropout=0.2)
            self.conv2 = GATConv(hidden_channels, hidden_channels, heads=4, concat=False, dropout=0.2)
            self.lin = nn.Linear(hidden_channels, 1)
            self.dropout = nn.Dropout(0.3)

        def forward(self, x, edge_index, batch):
            x = F.elu(self.conv1(x, edge_index))
            x = self.dropout(x)
            x = F.elu(self.conv2(x, edge_index))
            x = global_mean_pool(x, batch)
            return self.lin(x).squeeze(-1)

    class NoEdgeClassifier(nn.Module):
        def __init__(self):
            super().__init__()
            self.lin1 = nn.Linear(in_channels, hidden_channels)
            self.lin2 = nn.Linear(hidden_channels, hidden_channels)
            self.lin3 = nn.Linear(hidden_channels, 1)
            self.dropout = nn.Dropout(0.3)

        def forward(self, x, edge_index, batch):
            from torch_geometric.nn import global_mean_pool

            pooled = global_mean_pool(x, batch)
            h = F.elu(self.lin1(pooled))
            h = self.dropout(h)
            h = F.elu(self.lin2(h))
            return self.lin3(h).squeeze(-1)

    return GATClassifier() if use_attention else NoEdgeClassifier()


def prepare_fold_inputs(
    expr,
    train_mask,
    test_mask,
    n_genes,
    beta,
    graph_type,
    n_edges_per_gene,
    rng,
):
    """Prepare train/test matrices plus graph edges for one LOMO fold."""
    expr_train = expr.loc[train_mask]
    expr_test = expr.loc[test_mask]
    expr_train_filt, _ = filter_top_mad_genes(expr_train, n_genes)
    gene_names = expr_train_filt.columns.tolist()
    x_train_raw = expr_train_filt.values.astype(np.float32)
    x_test_raw = expr_test.loc[:, gene_names].values.astype(np.float32)

    if graph_type == "wgcna":
        tom = compute_tom_from_expr(x_train_raw.T, beta)
        edge_src, edge_dst = tom_to_edge_list(tom, n_edges_per_gene, n_genes, rng, graph_type)
    elif graph_type == "random":
        edge_src, edge_dst = tom_to_edge_list(None, n_edges_per_gene, n_genes, rng, graph_type)
    else:
        edge_src, edge_dst = tom_to_edge_list(None, n_edges_per_gene, n_genes, rng, "no_edges")

    scaler = StandardScaler()
    x_train = scaler.fit_transform(x_train_raw).astype(np.float32)
    x_test = scaler.transform(x_test_raw).astype(np.float32)
    return x_train, x_test, edge_src, edge_dst


def bootstrap_mean_ci_from_folds(fold_results, seed):
    """Estimate a CI on mean AUROC by resampling folds."""
    if not fold_results:
        return None, None

    values = np.array([f["auroc"] for f in fold_results], dtype=np.float32)
    if len(values) == 1:
        return float(values[0]), float(values[0])

    rng = np.random.default_rng(seed)
    boot_means = []
    for _ in range(5000):
        idx = rng.integers(0, len(values), size=len(values))
        boot_means.append(float(np.mean(values[idx])))

    return float(np.percentile(boot_means, 2.5)), float(np.percentile(boot_means, 97.5))


def train_gnn_lomo(
    tissue,
    graph_type,
    n_edges_per_gene,
    n_bootstrap,
    n_perm,
    seed,
    topology_scope="train_fold",
):
    """Run the full LOMO pipeline for one tissue x graph type."""
    import torch
    import torch.nn as nn
    from torch_geometric.loader import DataLoader

    print(f"\n{'=' * 60}")
    print(f"Tissue: {tissue}, Graph: {graph_type}")
    print(f"{'=' * 60}")

    expr, y, meta = load_tissue_data(tissue)
    missions = TISSUE_MISSIONS[tissue]

    with open(V4_EVAL / f"WGCNA_{tissue}_modules.json") as handle:
        wgcna = json.load(handle)
    n_genes = int(wgcna["n_genes_input"])
    beta = int(wgcna["soft_threshold_beta"])

    print(f"  Samples: {expr.shape[0]}, Full genes: {expr.shape[1]}, WGCNA target: {n_genes}")
    print(f"  Topology scope: {topology_scope}")

    use_attention = graph_type != "no_edges"
    device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
    print(f"  Device: {device}")

    fold_results = []
    fold_payloads = []
    global_expr_filt = None
    global_edge_src = None
    global_edge_dst = None

    if topology_scope == "full_dataset":
        global_expr_filt, _ = filter_top_mad_genes(expr, n_genes)
        global_rng = np.random.default_rng(seed)
        if graph_type == "wgcna":
            print(f"  Computing full-dataset TOM on {global_expr_filt.shape[0]} samples...")
            tom = compute_tom_from_expr(global_expr_filt.values.astype(np.float32).T, beta)
            global_edge_src, global_edge_dst = tom_to_edge_list(
                tom, n_edges_per_gene, n_genes, global_rng, graph_type
            )
        elif graph_type == "random":
            global_edge_src, global_edge_dst = tom_to_edge_list(
                None, n_edges_per_gene, n_genes, global_rng, graph_type
            )
        else:
            global_edge_src, global_edge_dst = tom_to_edge_list(
                None, n_edges_per_gene, n_genes, global_rng, "no_edges"
            )
        print(f"  Full-dataset edges (undirected): {len(global_edge_src) // 2:,}")

    for fold_idx, test_mission in enumerate(missions):
        test_mask = meta["mission"].values == test_mission
        train_mask = ~test_mask

        if test_mask.sum() == 0:
            print(f"  [SKIP] {test_mission}: no test samples")
            continue
        if len(np.unique(y[test_mask])) < 2:
            print(f"  [SKIP] {test_mission}: single class in test")
            continue
        if train_mask.sum() < 4:
            print(f"  [SKIP] {test_mission}: too few train samples")
            continue

        y_train = y[train_mask]
        y_test = y[test_mask]
        n_tr, n_te = int(train_mask.sum()), int(test_mask.sum())

        if topology_scope == "train_fold":
            rng = np.random.default_rng(seed + fold_idx)
            if graph_type == "wgcna":
                print(f"  Fold {test_mission}: computing train-fold TOM on {n_tr} samples...")
            x_train, x_test, edge_src, edge_dst = prepare_fold_inputs(
                expr,
                train_mask,
                test_mask,
                n_genes,
                beta,
                graph_type,
                n_edges_per_gene,
                rng,
            )
        else:
            gene_names = global_expr_filt.columns.tolist()
            x_train_raw = expr.loc[train_mask, gene_names].values.astype(np.float32)
            x_test_raw = expr.loc[test_mask, gene_names].values.astype(np.float32)
            scaler = StandardScaler()
            x_train = scaler.fit_transform(x_train_raw).astype(np.float32)
            x_test = scaler.transform(x_test_raw).astype(np.float32)
            edge_src, edge_dst = global_edge_src, global_edge_dst

        n_edges = len(edge_src) // 2
        dataset_tr = build_pyg_dataset(x_train, y_train, edge_src, edge_dst)
        dataset_te = build_pyg_dataset(x_test, y_test, edge_src, edge_dst)

        fold_seed = seed + fold_idx
        seed_torch(fold_seed)
        tr_generator = torch.Generator()
        tr_generator.manual_seed(fold_seed)

        loader_tr = DataLoader(dataset_tr, batch_size=32, shuffle=True, generator=tr_generator)
        loader_te = DataLoader(dataset_te, batch_size=64, shuffle=False)

        model = build_gnn_model(1, 64, use_attention=use_attention).to(device)
        optimizer = torch.optim.Adam(model.parameters(), lr=1e-3, weight_decay=1e-4)
        pos_weight = torch.tensor(
            [np.sum(y_train == 0) / max(np.sum(y_train == 1), 1)],
            dtype=torch.float32,
        ).to(device)
        criterion = nn.BCEWithLogitsLoss(pos_weight=pos_weight)

        best_train_loss = np.inf
        patience_cnt = 0
        patience = 15

        print(
            f"  Fold: {test_mission} | train={n_tr}, test={n_te}, edges={n_edges:,}",
            end="",
            flush=True,
        )
        t0 = time.time()
        for epoch in range(150):
            model.train()
            epoch_loss = 0.0
            for batch in loader_tr:
                batch = batch.to(device)
                logits = model(batch.x, batch.edge_index, batch.batch)
                loss = criterion(logits, batch.y.squeeze())
                optimizer.zero_grad()
                loss.backward()
                optimizer.step()
                epoch_loss += loss.item()

            epoch_loss /= max(len(loader_tr), 1)
            if epoch_loss < best_train_loss - 1e-4:
                best_train_loss = epoch_loss
                patience_cnt = 0
            else:
                patience_cnt += 1
            if patience_cnt >= patience:
                break

        model.eval()
        y_scores = []
        with torch.no_grad():
            for batch in loader_te:
                batch = batch.to(device)
                logits = model(batch.x, batch.edge_index, batch.batch)
                y_scores.extend(torch.sigmoid(logits).cpu().numpy().tolist())

        y_scores = np.asarray(y_scores, dtype=np.float32)
        auroc = roc_auc_score(y_test, y_scores)
        duration = time.time() - t0
        print(f" | AUROC={auroc:.3f} | {epoch + 1}ep | {duration:.0f}s")

        rng_boot = np.random.default_rng(seed + 999 + fold_idx)
        boot_aurocs = []
        for _ in range(n_bootstrap):
            idx = rng_boot.integers(0, len(y_test), size=len(y_test))
            if len(np.unique(y_test[idx])) < 2:
                continue
            try:
                boot_aurocs.append(roc_auc_score(y_test[idx], y_scores[idx]))
            except Exception:
                continue
        ci_lo = np.percentile(boot_aurocs, 2.5) if boot_aurocs else np.nan
        ci_hi = np.percentile(boot_aurocs, 97.5) if boot_aurocs else np.nan

        fold_results.append(
            {
                "test_mission": test_mission,
                "auroc": round(float(auroc), 4),
                "ci_lower": round(float(ci_lo), 4),
                "ci_upper": round(float(ci_hi), 4),
                "n_train": n_tr,
                "n_test": n_te,
                "epochs": epoch + 1,
                "n_edges": int(n_edges),
            }
        )
        fold_payloads.append(
            {
                "train_idx": np.where(train_mask)[0],
                "test_idx": np.where(test_mask)[0],
                "x_train": x_train,
                "x_test": x_test,
                "edge_src": edge_src.copy(),
                "edge_dst": edge_dst.copy(),
            }
        )

    if fold_results:
        mean_auroc = float(np.mean([fold["auroc"] for fold in fold_results]))
        mean_edges = float(np.mean([fold["n_edges"] for fold in fold_results]))

        perm_means = []
        print(f"\n  Permutation test ({n_perm} permutations)...")
        perm_rng = np.random.default_rng(seed + 1234)
        for perm_i in range(n_perm):
            y_perm = perm_rng.permutation(y)
            fold_aurocs_perm = []

            for fold_j, payload in enumerate(fold_payloads):
                y_train_p = y_perm[payload["train_idx"]]
                y_test_p = y_perm[payload["test_idx"]]
                if len(np.unique(y_train_p)) < 2 or len(np.unique(y_test_p)) < 2:
                    continue

                dataset_tr_p = build_pyg_dataset(
                    payload["x_train"], y_train_p, payload["edge_src"], payload["edge_dst"]
                )
                dataset_te_p = build_pyg_dataset(
                    payload["x_test"], y_test_p, payload["edge_src"], payload["edge_dst"]
                )

                perm_seed = seed + 10000 + perm_i * 100 + fold_j
                seed_torch(perm_seed)
                tr_generator_p = torch.Generator()
                tr_generator_p.manual_seed(perm_seed)

                loader_tr_p = DataLoader(
                    dataset_tr_p, batch_size=32, shuffle=True, generator=tr_generator_p
                )
                loader_te_p = DataLoader(dataset_te_p, batch_size=64, shuffle=False)

                model_p = build_gnn_model(1, 64, use_attention=use_attention).to(device)
                opt_p = torch.optim.Adam(model_p.parameters(), lr=1e-3, weight_decay=1e-4)
                pw = torch.tensor(
                    [np.sum(y_train_p == 0) / max(np.sum(y_train_p == 1), 1)],
                    dtype=torch.float32,
                ).to(device)
                crit_p = nn.BCEWithLogitsLoss(pos_weight=pw)

                for _ in range(50):
                    model_p.train()
                    for batch in loader_tr_p:
                        batch = batch.to(device)
                        logits = model_p(batch.x, batch.edge_index, batch.batch)
                        loss = crit_p(logits, batch.y.squeeze())
                        opt_p.zero_grad()
                        loss.backward()
                        opt_p.step()

                model_p.eval()
                ys_perm = []
                with torch.no_grad():
                    for batch in loader_te_p:
                        batch = batch.to(device)
                        logits = model_p(batch.x, batch.edge_index, batch.batch)
                        ys_perm.extend(torch.sigmoid(logits).cpu().numpy().tolist())
                try:
                    fold_aurocs_perm.append(roc_auc_score(y_test_p, np.asarray(ys_perm)))
                except Exception:
                    pass

            if fold_aurocs_perm:
                perm_means.append(float(np.mean(fold_aurocs_perm)))

        perm_p = (np.sum(np.asarray(perm_means) >= mean_auroc) + 1) / (len(perm_means) + 1)
        print(f"  perm_p = {perm_p:.3f} (n_perm={n_perm})")
    else:
        mean_auroc = np.nan
        mean_edges = np.nan
        perm_p = np.nan

    mean_ci_lo, mean_ci_hi = bootstrap_mean_ci_from_folds(
        fold_results, seed + 2026
    )

    return {
        "tissue": tissue,
        "graph_type": graph_type,
        "topology_scope": topology_scope,
        "n_edges_per_gene": n_edges_per_gene,
        "n_edges_total": int(round(mean_edges)) if not np.isnan(mean_edges) else None,
        "n_edges_per_fold": [int(fold["n_edges"]) for fold in fold_results],
        "n_genes": n_genes,
        "beta": beta,
        "mean_auroc": round(float(mean_auroc), 4) if not np.isnan(mean_auroc) else None,
        "mean_ci_lower": round(float(mean_ci_lo), 4) if mean_ci_lo is not None else None,
        "mean_ci_upper": round(float(mean_ci_hi), 4) if mean_ci_hi is not None else None,
        "perm_pvalue": round(float(perm_p), 4) if not np.isnan(perm_p) else None,
        "n_folds": len(fold_results),
        "folds": fold_results,
        "seed": seed,
        "n_bootstrap": n_bootstrap,
        "n_perm": n_perm,
    }


def load_pcalr_baseline(tissue):
    """Load PCA-LR AUROC from the v4 summary."""
    try:
        with open(M1_SUMMARY) as handle:
            data = json.load(handle)
        return data[tissue]["gene"]["pca_lr"]["auroc"]
    except Exception:
        return None


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--tissue", default="liver", help="Tissue name or 'all'")
    parser.add_argument(
        "--graph-type",
        default="wgcna",
        choices=["wgcna", "random", "no_edges", "all"],
    )
    parser.add_argument(
        "--topology-scope",
        default="train_fold",
        choices=["train_fold", "full_dataset"],
    )
    parser.add_argument("--n-edges-per-gene", type=int, default=10)
    parser.add_argument("--n-bootstrap", type=int, default=1000)
    parser.add_argument("--n-perm", type=int, default=100)
    parser.add_argument("--seed", type=int, default=42)
    args = parser.parse_args()

    tissues = list(TISSUE_MISSIONS.keys()) if args.tissue == "all" else [args.tissue]
    graph_types = ["wgcna", "random", "no_edges"] if args.graph_type == "all" else [args.graph_type]

    all_results = []
    for tissue in tissues:
        for graph_type in graph_types:
            out_file = OUT_DIR / f"GNN_{tissue}_{graph_type}.json"
            if out_file.exists():
                print(f"[SKIP] {out_file.name} already exists")
                with open(out_file) as handle:
                    all_results.append(json.load(handle))
                continue

            result = train_gnn_lomo(
                tissue=tissue,
                graph_type=graph_type,
                n_edges_per_gene=args.n_edges_per_gene,
                n_bootstrap=args.n_bootstrap,
                n_perm=args.n_perm,
                seed=args.seed,
                topology_scope=args.topology_scope,
            )
            result["pca_lr_auroc"] = load_pcalr_baseline(tissue)
            result["delta_vs_pca_lr"] = round(
                float(result["mean_auroc"] - result["pca_lr_auroc"]), 4
            ) if (result.get("mean_auroc") is not None and result.get("pca_lr_auroc") is not None) else None

            with open(out_file, "w") as handle:
                json.dump(result, handle, indent=2)
            print(f"  Saved: {out_file.name}")
            all_results.append(result)

    if len(all_results) > 1:
        summary = {
            "results": all_results,
            "tissues": tissues,
            "graph_types": graph_types,
            "topology_scope": args.topology_scope,
        }
        summary_path = OUT_DIR / "GNN_summary.json"
        with open(summary_path, "w") as handle:
            json.dump(summary, handle, indent=2)
        print(f"\nSummary saved: {summary_path}")

        print("\n" + "=" * 70)
        print(f"{'Tissue':<14} {'Graph':<12} {'AUROC':<8} {'perm_p':<8} {'PCA-LR':<8}")
        print("-" * 70)
        for result in all_results:
            auroc = result.get("mean_auroc")
            perm_p = result.get("perm_pvalue")
            pcalr = result.get("pca_lr_auroc")
            sig = "*" if perm_p and perm_p < 0.05 else ""
            print(
                f"{result['tissue']:<14} {result['graph_type']:<12} "
                f"{(auroc if auroc is not None else 'N/A'):<8} "
                f"{(str(perm_p) + sig) if perm_p is not None else 'N/A':<8} "
                f"{(pcalr if pcalr is not None else 'N/A'):<8}"
            )


if __name__ == "__main__":
    main()
