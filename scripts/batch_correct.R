#!/usr/bin/env Rscript
# batch_correct.R — GeneLab_benchmark: Batch Correction (DD-10)
#
# Applies ComBat-seq batch correction to cross-mission liver RNA-seq data.
# Used to investigate whether batch effects explain the low A1 LOMO AUROC.
#
# Method: ComBat-seq (sva package) — designed for RNA-seq count data (negative binomial)
# Batch variable: mission (OSD study) — accounts for technical variation between missions
#
# Output: corrected_counts.csv for use in generate_tasks.py
#
# Usage:
#   Rscript batch_correct.R --tissue liver
#   Rscript batch_correct.R --tissue liver --compare  # compare methods

suppressPackageStartupMessages({
  library(sva)      # ComBat-seq
  library(limma)    # removeBatchEffect (for comparison)
  library(argparse)
})

# ── CLI ─────────────────────────────────────────────────────────────────────────
parser <- ArgumentParser(description = "Apply batch correction to liver RNA-seq")
parser$add_argument("--tissue", default = "liver", help = "Tissue to process")
parser$add_argument("--method", default = "combat_seq",
                    choices = c("combat_seq", "limma", "limma_rbe", "none"),
                    help = "Batch correction method: combat_seq | limma | limma_rbe (all tissues)")
parser$add_argument("--compare", action = "store_true",
                    help = "Compare all methods side-by-side")
args <- parser$parse_args()

tissue  <- args$tissue
method  <- args$method

# Get script directory — works with Rscript --file= invocation
args_all  <- commandArgs(trailingOnly = FALSE)
file_args <- args_all[startsWith(args_all, "--file=")]
if (length(file_args) > 0) {
  script_file <- sub("^--file=", "", file_args[1])
  base_dir    <- dirname(dirname(normalizePath(script_file)))
} else {
  base_dir <- normalizePath(".")
}
processed_dir <- file.path(base_dir, "processed", "A_detection", tissue)

# ── Load data ────────────────────────────────────────────────────────────────────
cat(sprintf("\n=== Batch Correction: %s ===\n", toupper(tissue)))

# For limma_rbe, skip raw counts entirely — use existing log2_norm
if (method == "limma_rbe") {
  cat("[limma_rbe] Using existing log2_norm (no raw counts needed)\n")

  log2_norm_f <- file.path(processed_dir,
                           sprintf("%s_all_missions_log2_norm.csv", tissue))
  meta_f      <- file.path(processed_dir,
                           sprintf("%s_all_missions_metadata.csv", tissue))

  if (!file.exists(log2_norm_f)) stop(sprintf("Not found: %s", log2_norm_f))
  if (!file.exists(meta_f))      stop(sprintf("Not found: %s", meta_f))

  log2_expr <- read.csv(log2_norm_f, row.names = 1, check.names = FALSE)
  meta_df   <- read.csv(meta_f, row.names = 1, stringsAsFactors = FALSE)

  cat(sprintf("  Loaded log2_norm: %d samples x %d cols\n",
      nrow(log2_expr), ncol(log2_expr)))

  # Keep ENSMUSG gene columns only
  gene_cols <- grepl("^ENSMUSG", colnames(log2_expr))
  if (sum(gene_cols) > 0) {
    log2_expr <- log2_expr[, gene_cols, drop = FALSE]
  }
  cat(sprintf("  Gene columns: %d\n", ncol(log2_expr)))

  # Filter REMOVE samples
  if ("REMOVE" %in% colnames(meta_df)) {
    meta_df <- meta_df[meta_df$REMOVE != TRUE, , drop = FALSE]
  }

  # Align
  common_samples <- intersect(rownames(log2_expr), rownames(meta_df))
  cat(sprintf("  Common samples: %d\n", length(common_samples)))
  if (length(common_samples) < 10) stop("Too few common samples")

  log2_aligned <- log2_expr[common_samples, ]
  meta_aligned <- meta_df[common_samples, ]

  batch_rbe <- factor(meta_aligned$mission)
  cat(sprintf("  Batches: %s\n", paste(levels(batch_rbe), collapse = ", ")))
  cat(sprintf("  Samples per batch: %s\n",
      paste(table(batch_rbe), collapse = ", ")))

  # Transpose to genes x samples
  expr_mat <- t(as.matrix(log2_aligned))
  expr_mat[is.na(expr_mat)] <- 0

  # Apply limma::removeBatchEffect
  corrected_rbe <- removeBatchEffect(expr_mat, batch = batch_rbe)
  corrected_df  <- as.data.frame(t(corrected_rbe))

  cat(sprintf("  Done: %d x %d\n", nrow(corrected_df), ncol(corrected_df)))

  out_f <- file.path(processed_dir,
                     sprintf("%s_all_missions_log2_norm_limma_rbe.csv", tissue))
  write.csv(corrected_df, out_f)
  cat(sprintf("  Saved: %s\n", out_f))

  cat("\nDone (limma_rbe).\n")
  quit(save = "no", status = 0)
}

# Load unnormalized counts (ComBat-seq works on raw counts)
# For each mission, find the unnormalized count file
tissue_osd_map <- list(
  liver = list(
    "OSD-48"  = list(mission = "RR-1",  dir = "RR-1"),
    "OSD-137" = list(mission = "RR-3",  dir = "RR-3"),
    "OSD-245" = list(mission = "RR-6",  dir = "RR-6"),
    "OSD-379" = list(mission = "RR-8",  dir = "RR-8"),
    "OSD-242" = list(mission = "RR-9",  dir = "RR-9"),
    "OSD-686" = list(mission = "MHU-2", dir = "MHU-2", glds_prefix = "GLDS-617")
  )
)

osd_map <- tissue_osd_map[[tissue]]
if (is.null(osd_map)) stop(sprintf("Tissue '%s' not configured.", tissue))

data_dir <- file.path(base_dir, "data", "mouse", tissue)

# Load each mission's unnormalized counts
cat("Loading raw counts per mission...\n")
count_list   <- list()
batch_vector <- character()
metadata_list <- list()

for (osd_id in names(osd_map)) {
  info    <- osd_map[[osd_id]]
  mission <- info$mission
  mission_dir <- file.path(data_dir, info$dir)
  glds_prefix <- if (!is.null(info$glds_prefix)) info$glds_prefix else gsub("OSD", "GLDS", osd_id)

  # Find unnormalized counts file
  f_candidates <- c(
    file.path(mission_dir, sprintf("%s_rna_seq_RSEM_Unnormalized_Counts_GLbulkRNAseq.csv", glds_prefix)),
    file.path(mission_dir, sprintf("%s_rna_seq_Unnormalized_Counts.csv", glds_prefix))
  )
  f_counts <- f_candidates[file.exists(f_candidates)][1]

  if (is.na(f_counts)) {
    cat(sprintf("  [SKIP] %s: no unnormalized counts found\n", osd_id))
    next
  }

  cat(sprintf("  %s (%s): %s\n", osd_id, mission, basename(f_counts)))
  counts_raw <- read.csv(f_counts, row.names = 1, check.names = FALSE)
  counts_raw <- counts_raw[, sapply(counts_raw, is.numeric), drop = FALSE]
  counts_raw <- round(counts_raw)  # ensure integer (RSEM gives fractional counts)
  counts_raw[counts_raw < 0] <- 0

  count_list[[osd_id]] <- counts_raw
  batch_vector <- c(batch_vector, rep(mission, ncol(counts_raw)))
  metadata_list[[osd_id]] <- data.frame(
    sample  = colnames(counts_raw),
    mission = mission,
    osd_id  = osd_id,
    stringsAsFactors = FALSE
  )
}

if (length(count_list) < 2) {
  stop("Need at least 2 missions for batch correction.")
}

# ── Find common genes across all missions ────────────────────────────────────────
cat("\nFinding common genes...\n")
gene_sets <- lapply(count_list, rownames)
common_genes <- Reduce(intersect, gene_sets)
cat(sprintf("  Common genes: %d\n", length(common_genes)))

# Build combined count matrix (genes × samples)
combined_counts <- do.call(cbind, lapply(count_list, function(m) m[common_genes, , drop = FALSE]))
cat(sprintf("  Combined matrix: %d genes × %d samples\n",
    nrow(combined_counts), ncol(combined_counts)))

batch <- factor(batch_vector)
cat(sprintf("  Batch (mission) levels: %s\n", paste(levels(batch), collapse = ", ")))

# ── Apply batch correction ────────────────────────────────────────────────────────

if (method == "combat_seq" || args$compare) {
  cat("\n[1] ComBat-seq (sva, negative binomial model)...\n")
  tryCatch({
    corrected_combat <- ComBat_seq(
      counts = as.matrix(combined_counts),
      batch  = batch,
      group  = NULL  # no biological variable known a priori for correction
    )
    # Log2 transform
    log2_combat <- log2(corrected_combat + 1)
    cat(sprintf("  ✓ ComBat-seq complete. Shape: %d × %d\n",
        nrow(log2_combat), ncol(log2_combat)))

    # Save
    out_f <- file.path(processed_dir, sprintf("%s_combat_seq_log2_norm.csv", tissue))
    write.csv(t(log2_combat), out_f)  # transpose: samples × genes
    cat(sprintf("  ✓ Saved: %s\n", out_f))
  }, error = function(e) {
    cat(sprintf("  [ERROR] ComBat-seq failed: %s\n", conditionMessage(e)))
  })
}

if (method == "limma" || args$compare) {
  cat("\n[2] limma::removeBatchEffect (log-CPM)...\n")
  tryCatch({
    # Convert to log-CPM first
    library_sizes <- colSums(combined_counts)
    cpm_matrix    <- t(t(combined_counts) / library_sizes) * 1e6
    log2_cpm      <- log2(cpm_matrix + 1)

    # Remove batch effect
    corrected_limma <- removeBatchEffect(log2_cpm, batch = batch)
    cat(sprintf("  ✓ limma correction complete.\n"))

    # Save
    out_f <- file.path(processed_dir, sprintf("%s_limma_batch_corrected.csv", tissue))
    write.csv(t(corrected_limma), out_f)
    cat(sprintf("  ✓ Saved: %s\n", out_f))
  }, error = function(e) {
    cat(sprintf("  [ERROR] limma failed: %s\n", conditionMessage(e)))
  })
}

cat("\n✓ Batch correction complete.\n")
cat("Next: python scripts/generate_tasks.py --task A1 --batch-corrected\n")
