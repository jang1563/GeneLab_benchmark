#!/usr/bin/env Rscript
# wgcna_preservation.R — WGCNA module preservation across tissue pairs
#
# For each (reference, test) pair from 6 LOMO tissues:
#   modulePreservation() → Zsummary.pres, medianRank.pres per module
#   Zsummary > 10 = high preservation; > 2 = moderate
#
# Input:   v4/wgcna_inputs/{tissue}_expr.csv      (genes × samples)
#          v4/wgcna_outputs/{tissue}/module_assignments.csv
# Output:  v4/evaluation/WGCNA_preservation.json
#
# Usage:
#   Rscript wgcna_preservation.R

suppressPackageStartupMessages({
  library(WGCNA)
  library(jsonlite)
})

allowWGCNAThreads(8)

# ── Paths ─────────────────────────────────────────────────────────────────────
script_dir  <- tryCatch(dirname(sys.frame(1)$ofile), error = function(e) ".")
project_dir <- normalizePath(file.path(script_dir, "..", ".."))
inputs_dir  <- file.path(project_dir, "v4", "wgcna_inputs")
outputs_dir <- file.path(project_dir, "v4", "wgcna_outputs")
eval_dir    <- file.path(project_dir, "v4", "evaluation")
dir.create(eval_dir, recursive = TRUE, showWarnings = FALSE)

LOMO_TISSUES <- c("liver", "gastrocnemius", "kidney", "thymus", "eye", "skin")

# ── Load tissue data ──────────────────────────────────────────────────────────
cat("Loading tissue expression data and module assignments...\n")

tissue_data <- list()

for (tissue in LOMO_TISSUES) {
  expr_file <- file.path(inputs_dir, paste0(tissue, "_expr.csv"))
  mod_file  <- file.path(outputs_dir, tissue, "module_assignments.csv")

  if (!file.exists(expr_file)) {
    cat(sprintf("  SKIP [%s]: expression file not found\n", tissue))
    next
  }
  if (!file.exists(mod_file)) {
    cat(sprintf("  SKIP [%s]: module assignments not found (run wgcna_analysis.R first)\n", tissue))
    next
  }

  # Load expression (genes × samples → transpose to samples × genes)
  expr_mat <- read.csv(expr_file, row.names = 1, check.names = FALSE)
  datExpr  <- t(expr_mat)

  # Load module assignments
  mod_df <- read.csv(mod_file, stringsAsFactors = FALSE)
  # gene names in mod_df should match colnames(datExpr)
  gene_order <- match(colnames(datExpr), mod_df$gene)
  colors <- rep("grey", ncol(datExpr))
  valid  <- !is.na(gene_order)
  colors[valid] <- mod_df$module_color[gene_order[valid]]

  tissue_data[[tissue]] <- list(
    datExpr = datExpr,
    colors  = colors
  )
  cat(sprintf("  Loaded [%s]: %d samples × %d genes, %d modules\n",
              tissue, nrow(datExpr), ncol(datExpr),
              length(setdiff(unique(colors), "grey"))))
}

available_tissues <- names(tissue_data)
cat(sprintf("\nAvailable tissues (%d): %s\n", length(available_tissues),
            paste(available_tissues, collapse = ", ")))

if (length(available_tissues) < 2) {
  stop("Need at least 2 tissues for preservation analysis.")
}

# ── Run module preservation for all pairs ─────────────────────────────────────
pairs <- combn(available_tissues, 2, simplify = FALSE)
cat(sprintf("Running preservation for %d tissue pairs...\n", length(pairs)))

all_results <- list()

for (pair in pairs) {
  ref_tissue  <- pair[1]
  test_tissue <- pair[2]
  pair_key    <- paste0(ref_tissue, "_vs_", test_tissue)

  cat(sprintf("\n── %s ──\n", pair_key))

  ref_data  <- tissue_data[[ref_tissue]]
  test_data <- tissue_data[[test_tissue]]

  # Find common genes
  common_genes <- intersect(colnames(ref_data$datExpr), colnames(test_data$datExpr))
  if (length(common_genes) < 100) {
    cat(sprintf("  SKIP: only %d common genes (need ≥100)\n", length(common_genes)))
    all_results[[pair_key]] <- list(
      ref_tissue   = ref_tissue,
      test_tissue  = test_tissue,
      n_common_genes = length(common_genes),
      error = "Too few common genes"
    )
    next
  }

  cat(sprintf("  Common genes: %d\n", length(common_genes)))

  # Subset to common genes
  ref_expr_sub  <- ref_data$datExpr[, common_genes, drop = FALSE]
  test_expr_sub <- test_data$datExpr[, common_genes, drop = FALSE]
  ref_colors_sub  <- ref_data$colors[match(common_genes, colnames(ref_data$datExpr))]
  test_colors_sub <- test_data$colors[match(common_genes, colnames(test_data$datExpr))]

  # Build multiExpr list for modulePreservation
  multiExpr <- list(
    reference = list(data = ref_expr_sub),
    test      = list(data = test_expr_sub)
  )
  multiColors <- list(
    reference = ref_colors_sub,
    test      = test_colors_sub
  )

  # Run preservation (reference = first, test = second)
  # nPermutations=100 is standard for publication-grade results
  mp <- tryCatch({
    modulePreservation(
      multiData          = multiExpr,
      multiColor         = multiColors,
      dataIsExpr         = TRUE,
      referenceNetworks  = 1,
      nPermutations      = 100,
      randomSeed         = 42,
      networkType        = "signed",
      maxModuleSize      = 400,
      maxGoldModuleSize  = 400,
      quickCor           = 0,
      verbose            = 2
    )
  }, error = function(e) {
    cat(sprintf("  ERROR: %s\n", conditionMessage(e)))
    NULL
  })

  if (is.null(mp)) {
    all_results[[pair_key]] <- list(
      ref_tissue   = ref_tissue,
      test_tissue  = test_tissue,
      n_common_genes = length(common_genes),
      error = "modulePreservation failed"
    )
    next
  }

  # Extract preservation statistics
  # mp$preservation$Z$ref.ref → Z-scores for each module
  Z_stats     <- mp$preservation$Z$ref.ref[["inColumnsAlsoPresentIn.test"]]
  obs_stats   <- mp$preservation$observed$ref.ref[["inColumnsAlsoPresentIn.test"]]
  z_summary   <- Z_stats[, "Zsummary.pres", drop = FALSE]
  median_rank <- Z_stats[, "medianRank.pres", drop = FALSE]

  # Build per-module result
  modules <- rownames(z_summary)
  modules <- modules[!modules %in% c("gold", "grey")]

  module_results <- list()
  for (mod in modules) {
    z_val    <- z_summary[mod, "Zsummary.pres"]
    rank_val <- median_rank[mod, "medianRank.pres"]
    preservation_level <- if (!is.na(z_val) && z_val > 10) "high" else
                          if (!is.na(z_val) && z_val > 2)  "moderate" else "low"

    module_results[[mod]] <- list(
      Zsummary       = round(z_val, 4),
      medianRank     = round(rank_val, 4),
      preservation   = preservation_level
    )
  }

  n_high <- sum(sapply(module_results, function(x) x$preservation == "high"))
  n_mod  <- sum(sapply(module_results, function(x) x$preservation == "moderate"))
  cat(sprintf("  Modules: %d high preservation, %d moderate\n", n_high, n_mod))

  all_results[[pair_key]] <- list(
    ref_tissue       = ref_tissue,
    test_tissue      = test_tissue,
    n_common_genes   = length(common_genes),
    n_ref_modules    = length(setdiff(unique(ref_colors_sub), "grey")),
    n_test_modules   = length(setdiff(unique(test_colors_sub), "grey")),
    n_tested_modules = length(modules),
    modules          = module_results,
    n_high_preserved = n_high,
    n_moderate_preserved = n_mod
  )
}

# ── Save results ──────────────────────────────────────────────────────────────
out_path <- file.path(eval_dir, "WGCNA_preservation.json")

output <- list(
  n_tissues     = length(available_tissues),
  tissues       = available_tissues,
  n_pairs       = length(pairs),
  n_perm        = 100,
  z_high_threshold = 10,
  z_moderate_threshold = 2,
  pairs         = all_results,
  timestamp     = format(Sys.time(), "%Y-%m-%dT%H:%M:%S")
)

write_json(output, out_path, pretty = TRUE, auto_unbox = TRUE)
cat(sprintf("\nPreservation results saved to: %s\n", out_path))

# ── Summary ───────────────────────────────────────────────────────────────────
cat("\n═══ Preservation Summary ═══\n")
for (key in names(all_results)) {
  res <- all_results[[key]]
  if (!is.null(res$error)) {
    cat(sprintf("  %-35s ERROR: %s\n", key, res$error))
  } else {
    cat(sprintf("  %-35s high=%d, moderate=%d (of %d modules)\n",
                key, res$n_high_preserved, res$n_moderate_preserved, res$n_tested_modules))
  }
}
