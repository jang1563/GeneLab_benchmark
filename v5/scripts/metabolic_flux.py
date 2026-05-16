#!/usr/bin/env python
"""GeneLabBench v5 Phase 3: Metabolic Flux Modeling via E-Flux + COBRApy.

Constrains mouse GEM (iMM1865) reaction bounds by gene expression,
runs FBA for FLT vs GC conditions, compares subsystem-level fluxes.

Usage:
    python metabolic_flux.py --tissue liver
    python metabolic_flux.py  # all 6 LOMO tissues
"""
import argparse
import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd

# ── paths ──
BASE_DIR = Path(__file__).resolve().parent.parent.parent
V5_DATA_DIR = BASE_DIR / "v5" / "data"
V5_EVAL_DIR = BASE_DIR / "v5" / "evaluation"
SYMBOL_MAP_PATH = BASE_DIR / "processed" / "ensembl_symbol_map.csv"
PATHWAY_DIR = BASE_DIR / "processed" / "pathway_scores"

sys.path.insert(0, str(BASE_DIR / "v4" / "scripts"))
from v4_utils import (
    TISSUE_MISSIONS, LOMO_TISSUES,
    load_metadata, load_gene_features, align_features_with_meta, encode_labels,
)

# Only 6 LOMO tissues for metabolic modeling (lung/colon too few samples)
FLUX_TISSUES = LOMO_TISSUES


def load_model():
    """Load iMM1865 mouse GEM via COBRApy."""
    import cobra

    model_path = V5_DATA_DIR / "iMM1865.xml.gz"
    if not model_path.exists():
        # Try uncompressed
        model_path = V5_DATA_DIR / "iMM1865.xml"
        if not model_path.exists():
            raise FileNotFoundError(
                f"iMM1865 model not found at {V5_DATA_DIR}. "
                "Download from https://bigg.ucsd.edu/static/models/iMM1865.xml.gz"
            )

    model = cobra.io.read_sbml_model(str(model_path))
    print(f"  Model: {model.id} — {len(model.genes)} genes, "
          f"{len(model.reactions)} reactions, {len(model.metabolites)} metabolites")
    return model


def map_genes_to_model(model, expr_symbols):
    """Map expression gene symbols to model gene IDs.

    iMM1865 uses Entrez IDs as gene.id (e.g. '11754') and
    mouse symbols as gene.name (e.g. 'Aoc3').
    """
    model_gene_ids = {g.id for g in model.genes}
    # Primary: name → id (symbol → entrez)
    model_gene_names = {g.name: g.id for g in model.genes if g.name}
    # Case-insensitive fallback
    model_gene_names_lower = {g.name.lower(): g.id for g in model.genes if g.name}

    mapped = {}
    for sym in expr_symbols:
        if sym in model_gene_names:
            mapped[sym] = model_gene_names[sym]
        elif sym in model_gene_ids:
            mapped[sym] = sym
        elif sym.lower() in model_gene_names_lower:
            mapped[sym] = model_gene_names_lower[sym.lower()]

    return mapped


def eflux_constrain(model, gene_expression, gene_mapping):
    """Apply E-Flux constraints: set reaction bounds by gene expression.

    For GPR rules:
      AND → min(expression of genes)
      OR  → sum(expression of genes)

    Expression values are normalized to [0, 1] relative to the maximum
    observed expression, then used to scale default reaction bounds.
    This ensures constraints are actually binding (not always > default ub).
    """
    import cobra

    model_copy = model.copy()

    # Build expression lookup for model gene IDs
    expr_lookup = {}
    for sym, model_id in gene_mapping.items():
        if sym in gene_expression:
            expr_lookup[model_id] = max(0, gene_expression[sym])

    # Default expression for unmapped genes (median of all mapped)
    if expr_lookup:
        default_expr = float(np.median(list(expr_lookup.values())))
    else:
        default_expr = 1.0

    # First pass: evaluate all GPR rules to get raw expression values
    rxn_expr = {}
    for rxn in model_copy.reactions:
        if not rxn.gene_reaction_rule:
            continue
        try:
            rxn_expr[rxn.id] = _evaluate_gpr(rxn.gene_reaction_rule, expr_lookup, default_expr)
        except Exception:
            continue

    # Normalize expression values to [0, 1] relative to max
    if rxn_expr:
        max_expr = max(rxn_expr.values())
        if max_expr > 0:
            for rid in rxn_expr:
                rxn_expr[rid] = rxn_expr[rid] / max_expr
        else:
            max_expr = 1.0
    else:
        max_expr = 1.0

    # Second pass: apply normalized constraints
    n_constrained = 0
    n_actually_binding = 0
    for rxn in model_copy.reactions:
        if rxn.id not in rxn_expr:
            continue

        norm_expr = rxn_expr[rxn.id]  # in [0, 1]

        # Scale bounds by normalized expression
        # A gene with max expression keeps original bounds;
        # a gene with 0 expression gets bounds = 0 (reaction shut off)
        orig_ub = rxn.upper_bound
        new_ub = orig_ub * norm_expr
        rxn.upper_bound = new_ub

        if rxn.lower_bound < 0:
            orig_lb = rxn.lower_bound
            rxn.lower_bound = orig_lb * norm_expr

        n_constrained += 1
        if new_ub < orig_ub - 1e-6:
            n_actually_binding += 1

    print(f"    Normalized to [0,1] (max raw expr={max_expr:.2f}), "
          f"{n_actually_binding}/{n_constrained} constraints binding")

    return model_copy, n_constrained


def _evaluate_gpr(gpr_string, expr_lookup, default):
    """Evaluate Gene-Protein-Reaction rule using expression values.

    Simple parser: handles 'and', 'or', parentheses.
    AND → min, OR → sum.
    """
    # Replace gene IDs with expression values
    tokens = gpr_string.replace("(", " ( ").replace(")", " ) ").split()
    values = []
    current_op = None

    def get_expr(gene_id):
        return expr_lookup.get(gene_id, default)

    # Simple recursive evaluation
    def evaluate(tokens, pos=0):
        vals = []
        ops = []
        i = pos
        while i < len(tokens):
            token = tokens[i]
            if token == "(":
                val, i = evaluate(tokens, i + 1)
                vals.append(val)
            elif token == ")":
                break
            elif token.lower() == "and":
                ops.append("and")
            elif token.lower() == "or":
                ops.append("or")
            else:
                # Gene ID
                vals.append(get_expr(token))
            i += 1

        if not vals:
            return default, i

        # Process AND first (higher precedence), then OR
        # Simple approach: if any 'and', use min; if any 'or', use sum
        if "and" in ops:
            result = min(vals)
        elif "or" in ops:
            result = sum(vals)
        else:
            result = vals[0]

        return result, i

    result, _ = evaluate(tokens)
    return result


def run_fba(model, condition_name):
    """Run pFBA (parsimonious FBA) to get unique flux solution.

    pFBA first maximizes biomass, then minimizes total flux subject to
    the optimal objective. This resolves LP degeneracy where standard FBA
    returns arbitrary flux distributions.
    """
    import cobra

    try:
        # Use pFBA to get unique, biologically meaningful flux distribution
        solution = cobra.flux_analysis.pfba(model)
        if solution.status == "optimal":
            fluxes = {rxn.id: solution.fluxes[rxn.id] for rxn in model.reactions}
            obj_value = solution.objective_value
            print(f"    {condition_name}: pFBA optimal, objective={obj_value:.4f}")
            return fluxes, obj_value
        else:
            print(f"    {condition_name}: pFBA status={solution.status}")
            return None, 0.0
    except Exception as e:
        # Fallback to standard FBA if pFBA fails
        print(f"    {condition_name}: pFBA failed ({e}), trying standard FBA...")
        try:
            solution = model.optimize()
            if solution.status == "optimal":
                fluxes = {rxn.id: solution.fluxes[rxn.id] for rxn in model.reactions}
                obj_value = solution.objective_value
                print(f"    {condition_name}: FBA optimal, objective={obj_value:.4f}")
                return fluxes, obj_value
        except Exception as e2:
            print(f"    {condition_name}: FBA also failed — {e2}")
        return None, 0.0


def compare_fluxes(fluxes_flt, fluxes_gc, model):
    """Compare FLT vs GC fluxes, grouped by subsystem."""
    subsystem_changes = {}

    for rxn in model.reactions:
        rid = rxn.id
        subsystem = rxn.subsystem or "Unknown"

        flux_flt = fluxes_flt.get(rid, 0)
        flux_gc = fluxes_gc.get(rid, 0)
        diff = flux_flt - flux_gc

        if subsystem not in subsystem_changes:
            subsystem_changes[subsystem] = {
                "reactions": [],
                "flux_diffs": [],
                "n_reactions": 0,
            }
        subsystem_changes[subsystem]["reactions"].append(rid)
        subsystem_changes[subsystem]["flux_diffs"].append(diff)
        subsystem_changes[subsystem]["n_reactions"] += 1

    # Summarize per subsystem
    subsystem_summary = {}
    for ss, data in subsystem_changes.items():
        diffs = np.array(data["flux_diffs"])
        subsystem_summary[ss] = {
            "n_reactions": data["n_reactions"],
            "mean_flux_diff": round(float(np.mean(diffs)), 6),
            "median_flux_diff": round(float(np.median(diffs)), 6),
            "abs_mean_diff": round(float(np.mean(np.abs(diffs))), 6),
            "n_increased": int(np.sum(diffs > 1e-6)),
            "n_decreased": int(np.sum(diffs < -1e-6)),
        }

    return subsystem_summary


def load_ssgsea_for_validation(tissue):
    """Load ssGSEA pathway scores for validation scatter plot."""
    # Try hallmark
    hallmark_dir = PATHWAY_DIR / tissue
    if not hallmark_dir.exists():
        return None

    # Find hallmark GSVA/ssGSEA files
    for pattern in ["*hallmark*", "*HALLMARK*"]:
        files = list(hallmark_dir.glob(pattern))
        if files:
            df = pd.read_csv(files[0], index_col=0)
            return df
    return None


def run_tissue(tissue, model):
    """Run E-Flux + FBA for one tissue."""
    print(f"\n{'='*60}")
    print(f"Tissue: {tissue}")
    print(f"{'='*60}")

    # Load data
    meta = load_metadata(tissue)
    genes = load_gene_features(tissue)
    genes, meta = align_features_with_meta(genes, meta)
    y, valid = encode_labels(meta)

    genes_valid = genes.loc[valid]
    y_valid = y[valid].astype(int)

    # Map Ensembl → Symbol
    sym_map = pd.read_csv(SYMBOL_MAP_PATH)
    ens_to_sym = dict(zip(sym_map["ENSEMBL"], sym_map["SYMBOL"]))

    expr = genes_valid.rename(columns=ens_to_sym)
    expr = expr[[c for c in expr.columns if not str(c).startswith("ENSMUSG")]]
    expr = expr.loc[:, ~expr.columns.duplicated()]

    print(f"  Samples: {len(expr)} (FLT={sum(y_valid==1)}, GC={sum(y_valid==0)})")
    print(f"  Genes (symbols): {expr.shape[1]}")

    # Map to model genes
    gene_mapping = map_genes_to_model(model, expr.columns.tolist())
    coverage = len(gene_mapping) / len(model.genes) * 100
    print(f"  Gene mapping: {len(gene_mapping)}/{len(model.genes)} "
          f"({coverage:.1f}% coverage)")

    if coverage < 30:
        print(f"  [WARN] Very low gene coverage ({coverage:.1f}%). Results unreliable.")

    # Compute mean expression per condition
    flt_mask = y_valid.values == 1
    gc_mask = y_valid.values == 0

    mean_flt = expr.loc[flt_mask].mean(axis=0).to_dict()
    mean_gc = expr.loc[gc_mask].mean(axis=0).to_dict()

    # E-Flux: constrain model with FLT expression
    print("  Running E-Flux (FLT)...")
    model_flt, n_flt_constrained = eflux_constrain(model, mean_flt, gene_mapping)
    print(f"    Constrained {n_flt_constrained} reactions")

    # E-Flux: constrain model with GC expression
    print("  Running E-Flux (GC)...")
    model_gc, n_gc_constrained = eflux_constrain(model, mean_gc, gene_mapping)
    print(f"    Constrained {n_gc_constrained} reactions")

    # FBA
    fluxes_flt, obj_flt = run_fba(model_flt, "FLT")
    fluxes_gc, obj_gc = run_fba(model_gc, "GC")

    if fluxes_flt is None or fluxes_gc is None:
        return {
            "tissue": tissue,
            "status": "fba_infeasible",
            "gene_coverage_pct": round(coverage, 1),
        }

    # Compare fluxes by subsystem
    subsystem_summary = compare_fluxes(fluxes_flt, fluxes_gc, model)

    # Sort by absolute mean difference
    top_subsystems = sorted(subsystem_summary.items(),
                           key=lambda x: abs(x[1]["abs_mean_diff"]),
                           reverse=True)[:30]

    return {
        "tissue": tissue,
        "status": "success",
        "n_samples_flt": int(sum(flt_mask)),
        "n_samples_gc": int(sum(gc_mask)),
        "gene_coverage_pct": round(coverage, 1),
        "n_genes_mapped": len(gene_mapping),
        "n_model_genes": len(model.genes),
        "fba_objective_flt": round(obj_flt, 6),
        "fba_objective_gc": round(obj_gc, 6),
        "fba_objective_diff": round(obj_flt - obj_gc, 6),
        "n_subsystems": len(subsystem_summary),
        "top_subsystems": {k: v for k, v in top_subsystems},
        "all_subsystems": subsystem_summary,
    }


def main():
    parser = argparse.ArgumentParser(description="v5 Phase 3: Metabolic Flux Modeling")
    parser.add_argument("--tissue", type=str, default=None)
    args = parser.parse_args()

    V5_EVAL_DIR.mkdir(parents=True, exist_ok=True)

    print("Loading iMM1865 model...")
    model = load_model()

    tissues = [args.tissue] if args.tissue else FLUX_TISSUES

    for tissue in tissues:
        result = run_tissue(tissue, model)

        out_path = V5_EVAL_DIR / f"metabolic_flux_{tissue}.json"
        with open(out_path, "w") as f:
            json.dump(result, f, indent=2)
        print(f"  Saved: {out_path.name}")

    print("\nDone!")


if __name__ == "__main__":
    main()
