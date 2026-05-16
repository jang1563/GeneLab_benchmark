#!/usr/bin/env Rscript
# wgcna_analysis.R — Per-tissue WGCNA co-expression network analysis
#
# Input:  v4/wgcna_inputs/{tissue}_expr.csv   (genes × samples)
#         v4/wgcna_inputs/{tissue}_traits.csv  (samples × traits)
# Output: v4/wgcna_outputs/{tissue}/
#           module_assignments.csv   — gene, module_color, kME
#           eigengenes.csv           — samples × module eigengenes
#           trait_correlations.csv   — modules × traits (r, p)
#           soft_threshold.csv       — beta, R2, model.fit
#           summary.json             — metadata: n_modules, beta, etc.
#
# Usage:
#   Rscript wgcna_analysis.R --tissue liver
#   Rscript wgcna_analysis.R --tissue gastrocnemius

suppressPackageStartupMessages({
  library(WGCNA)
  library(argparse)
  library(jsonlite)
})

# ── CLI args ──────────────────────────────────────────────────────────────────
parser <- ArgumentParser(description = "WGCNA per-tissue co-expression analysis")
parser$add_argument("--tissue", required = TRUE,
                    help = "Tissue name (e.g., liver, gastrocnemius)")
parser$add_argument("--n-threads", type = "integer", default = 8,
                    help = "Threads for WGCNA (default: 8)")
args <- parser$parse_args()

tissue    <- args$tissue
n_threads <- args$n_threads

# ── Paths ─────────────────────────────────────────────────────────────────────
# Robust path resolution: works both with Rscript and source()
script_dir <- tryCatch({
  # When invoked via Rscript, commandArgs() contains the script path
  args_all <- commandArgs(trailingOnly = FALSE)
  script_arg <- grep("--file=", args_all, value = TRUE)
  if (length(script_arg) > 0) {
    dirname(normalizePath(sub("--file=", "", script_arg[1])))
  } else {
    dirname(normalizePath(sys.frame(1)$ofile))
  }
}, error = function(e) ".")

project_dir  <- normalizePath(file.path(script_dir, "..", ".."))
inputs_dir   <- file.path(project_dir, "v4", "wgcna_inputs")
outputs_dir  <- file.path(project_dir, "v4", "wgcna_outputs", tissue)
dir.create(outputs_dir, recursive = TRUE, showWarnings = FALSE)

expr_file   <- file.path(inputs_dir, paste0(tissue, "_expr.csv"))
traits_file <- file.path(inputs_dir, paste0(tissue, "_traits.csv"))

if (!file.exists(expr_file))   stop(paste("Expression file not found:", expr_file))
if (!file.exists(traits_file)) stop(paste("Traits file not found:", traits_file))

cat(sprintf("[%s] Loading data...\n", tissue))

# ── Load data ─────────────────────────────────────────────────────────────────
# expr: genes × samples → WGCNA expects samples × genes → transpose
expr_mat <- read.csv(expr_file, row.names = 1, check.names = FALSE)
cat(sprintf("  Expression: %d genes × %d samples\n", nrow(expr_mat), ncol(expr_mat)))

# datExpr: samples × genes (WGCNA convention)
datExpr <- t(expr_mat)
cat(sprintf("  datExpr: %d samples × %d genes\n", nrow(datExpr), ncol(datExpr)))

# traits
datTraits <- read.csv(traits_file, row.names = 1, check.names = FALSE)
# Align samples (expr may have different order)
common_samples <- intersect(rownames(datExpr), rownames(datTraits))
if (length(common_samples) < 5) {
  stop(sprintf("Too few common samples between expr and traits: %d", length(common_samples)))
}
datExpr   <- datExpr[common_samples, , drop = FALSE]
datTraits <- datTraits[common_samples, , drop = FALSE]
cat(sprintf("  Aligned: %d samples\n", length(common_samples)))

# ── WGCNA thread setup ────────────────────────────────────────────────────────
allowWGCNAThreads(n_threads)

# ── Soft threshold selection ──────────────────────────────────────────────────
cat(sprintf("[%s] Selecting soft threshold...\n", tissue))

powers <- c(1:20, seq(22, 30, by = 2))

# For small-n tissues (gastrocnemius n=32, eye n=41): fewer powers
is_small_n <- nrow(datExpr) < 50

sft <- pickSoftThreshold(
  datExpr,
  powerVector = powers,
  networkType  = "signed hybrid",
  verbose      = 2,
  corFnc       = "bicor"  # biweight midcorrelation (more robust than Pearson)
)

# Find first power with R² ≥ 0.80
sft_df <- sft$fitIndices
r2_threshold <- 0.80
good_idx <- which(sft_df$SFT.R.sq >= r2_threshold & sft_df$slope < 0)

if (length(good_idx) > 0) {
  chosen_beta <- sft_df$Power[min(good_idx)]
  cat(sprintf("  Chosen beta=%d (R²=%.3f >= %.2f)\n",
              chosen_beta, sft_df$SFT.R.sq[min(good_idx)], r2_threshold))
} else {
  # Fallback: use beta=6 for signed hybrid (recommended for n<50)
  chosen_beta <- ifelse(is_small_n, 6, 8)
  best_r2 <- max(sft_df$SFT.R.sq, na.rm = TRUE)
  cat(sprintf("  WARNING: R² never reached %.2f (max=%.3f). Using fallback beta=%d\n",
              r2_threshold, best_r2, chosen_beta))
}

# Save soft threshold data
write.csv(sft_df, file.path(outputs_dir, "soft_threshold.csv"), row.names = FALSE)

# ── Adjacency and TOM ─────────────────────────────────────────────────────────
cat(sprintf("[%s] Computing adjacency (beta=%d)...\n", tissue, chosen_beta))

adjacency_mat <- adjacency(
  datExpr,
  power        = chosen_beta,
  type         = "signed hybrid",
  corFnc       = "bicor"
)

cat(sprintf("[%s] Computing TOM dissimilarity...\n", tissue))
dissTOM <- TOMdist(adjacency_mat, TOMType = "signed")
geneTree <- hclust(as.dist(dissTOM), method = "average")

# ── Module detection ──────────────────────────────────────────────────────────
cat(sprintf("[%s] Detecting modules...\n", tissue))

# Small-n tissues: lower minModuleSize threshold
min_mod_size <- ifelse(is_small_n, 20, 30)
cat(sprintf("  minModuleSize=%d (small_n=%s)\n", min_mod_size, is_small_n))

dynamicMods <- cutreeDynamic(
  dendro       = geneTree,
  distM        = dissTOM,
  deepSplit    = 2,
  pamRespectsDendro = FALSE,
  minClusterSize = min_mod_size
)

dynamicColors <- labels2colors(dynamicMods)
n_mod_raw <- length(unique(dynamicColors[dynamicColors != "grey"]))
cat(sprintf("  Before merging: %d modules (%d grey/unassigned)\n",
            n_mod_raw, sum(dynamicColors == "grey")))

# ── Merge similar modules ─────────────────────────────────────────────────────
MEList   <- moduleEigengenes(datExpr, colors = dynamicColors)
MEs      <- MEList$eigengenes
mergedMods <- mergeCloseModules(
  datExpr,
  dynamicColors,
  cutHeight = 0.25,
  verbose   = 2
)
mergedColors <- mergedMods$colors
mergedMEs    <- mergedMods$newMEs

n_mod_final <- length(unique(mergedColors[mergedColors != "grey"]))
cat(sprintf("  After merging: %d modules\n", n_mod_final))

# ── Module eigengenes ─────────────────────────────────────────────────────────
MEs_ordered <- orderMEs(mergedMEs)
module_names <- colnames(MEs_ordered)  # "ME<color>"
cat(sprintf("  Modules: %s\n", paste(sub("^ME", "", module_names), collapse = ", ")))

# ── Module-trait correlations ─────────────────────────────────────────────────
cat(sprintf("[%s] Computing module-trait correlations...\n", tissue))

# Only use numeric trait columns (drop all-NA columns)
trait_cols_numeric <- sapply(datTraits, is.numeric) | sapply(datTraits, is.logical)
datTraits_num <- datTraits[, trait_cols_numeric, drop = FALSE]

# Remove columns with near-zero variance (constant)
trait_var <- apply(datTraits_num, 2, var, na.rm = TRUE)
datTraits_num <- datTraits_num[, !is.na(trait_var) & trait_var > 0, drop = FALSE]

if (ncol(datTraits_num) == 0) {
  cat("  WARNING: No valid trait columns. Skipping module-trait correlation.\n")
  moduleTraitCor <- matrix(NA, nrow = ncol(MEs_ordered), ncol = 0)
  moduleTraitPvalue <- moduleTraitCor
} else {
  # Use bicor for robustness; handle missing trait values (condition=NA for AG)
  n_samples_notNA <- colSums(!is.na(datTraits_num))
  cat(sprintf("  Traits (%d): %s\n", ncol(datTraits_num), paste(colnames(datTraits_num), collapse = ", ")))
  cat(sprintf("  N non-NA per trait: %s\n", paste(n_samples_notNA, collapse = ", ")))

  # Build correlation matrix, treating NA as pairwise complete
  moduleTraitCor <- cor(MEs_ordered, datTraits_num, use = "pairwise.complete.obs", method = "pearson")
  nUsed <- matrix(NA, nrow = nrow(moduleTraitCor), ncol = ncol(moduleTraitCor))
  for (j in seq_len(ncol(datTraits_num))) {
    valid <- !is.na(datTraits_num[, j])
    nUsed[, j] <- sum(valid)
  }
  moduleTraitPvalue <- corPvalueStudent(moduleTraitCor, nUsed)
}

# ── Module membership (kME) ───────────────────────────────────────────────────
cat(sprintf("[%s] Computing module membership (kME)...\n", tissue))

kME <- signedKME(datExpr, MEs_ordered, outputColumnName = "kME")
# kME columns: kME{color} for each module

# ── Save results ──────────────────────────────────────────────────────────────
cat(sprintf("[%s] Saving results to %s\n", tissue, outputs_dir))

# 1. Module assignments
mod_df <- data.frame(
  gene         = colnames(datExpr),
  module_color = mergedColors,
  stringsAsFactors = FALSE
)
# Attach kME for the assigned module
for (i in seq_len(nrow(mod_df))) {
  gene_mod <- mod_df$module_color[i]
  kme_col  <- paste0("kME", gene_mod)
  if (kme_col %in% colnames(kME)) {
    mod_df$kME_assigned[i] <- kME[mod_df$gene[i], kme_col]
  } else {
    mod_df$kME_assigned[i] <- NA
  }
}
# Also attach kME for all modules (full kME matrix)
mod_df_full <- cbind(mod_df, kME[rownames(kME) %in% mod_df$gene, , drop = FALSE])
write.csv(mod_df_full, file.path(outputs_dir, "module_assignments.csv"), row.names = FALSE)

# 2. Module eigengenes (samples × modules)
write.csv(MEs_ordered, file.path(outputs_dir, "eigengenes.csv"), row.names = TRUE)

# 3. Trait correlations
if (nrow(moduleTraitCor) > 0 && ncol(moduleTraitCor) > 0) {
  trait_cor_df <- as.data.frame(moduleTraitCor)
  trait_cor_df$module <- rownames(moduleTraitCor)
  write.csv(trait_cor_df, file.path(outputs_dir, "trait_correlations_r.csv"), row.names = FALSE)

  trait_p_df <- as.data.frame(moduleTraitPvalue)
  trait_p_df$module <- rownames(moduleTraitPvalue)
  write.csv(trait_p_df, file.path(outputs_dir, "trait_correlations_p.csv"), row.names = FALSE)
}

# 4. Module sizes
mod_sizes <- table(mergedColors)
mod_size_df <- data.frame(
  module_color = names(mod_sizes),
  n_genes = as.integer(mod_sizes)
)
write.csv(mod_size_df, file.path(outputs_dir, "module_sizes.csv"), row.names = FALSE)

# 5. Summary JSON
n_sig_traits <- 0
if (nrow(moduleTraitCor) > 0 && ncol(moduleTraitCor) > 0) {
  n_sig_traits <- sum(moduleTraitPvalue < 0.05, na.rm = TRUE)
}

summary_list <- list(
  tissue           = tissue,
  n_samples        = nrow(datExpr),
  n_genes_input    = ncol(datExpr),
  soft_threshold_beta = chosen_beta,
  r2_achieved      = ifelse(length(good_idx) > 0,
                             round(sft_df$SFT.R.sq[min(good_idx)], 4), NA),
  r2_max           = round(max(sft_df$SFT.R.sq, na.rm = TRUE), 4),
  r2_threshold_met = length(good_idx) > 0,
  min_module_size  = min_mod_size,
  n_modules_raw    = n_mod_raw,
  n_modules_final  = n_mod_final,
  n_grey_genes     = sum(mergedColors == "grey"),
  modules          = setdiff(unique(mergedColors), "grey"),
  n_genes_per_module = as.list(setNames(
    as.integer(mod_sizes[names(mod_sizes) != "grey"]),
    names(mod_sizes)[names(mod_sizes) != "grey"]
  )),
  n_traits         = ncol(datTraits_num),
  n_significant_module_trait = n_sig_traits,
  network_type     = "signed hybrid",
  cor_function     = "bicor",
  merge_cut_height = 0.25,
  deep_split       = 2
)

write_json(summary_list, file.path(outputs_dir, "summary.json"), pretty = TRUE, auto_unbox = TRUE)

cat(sprintf("\n[%s] WGCNA complete:\n", tissue))
cat(sprintf("  Modules: %d (+ grey)\n", n_mod_final))
cat(sprintf("  Beta: %d, R²=%.3f\n", chosen_beta, ifelse(length(good_idx) > 0,
                                                           sft_df$SFT.R.sq[min(good_idx)], NA)))
cat(sprintf("  Significant module-trait: %d associations (p<0.05)\n", n_sig_traits))
cat(sprintf("  Outputs: %s\n", outputs_dir))
