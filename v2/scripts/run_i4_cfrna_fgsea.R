#!/usr/bin/env Rscript
# run_i4_cfrna_fgsea.R — GeneLab_benchmark v2, Category E2
#
# Runs fGSEA on Inspiration4 (I4) human cfRNA DE data against MSigDB Hallmark gene sets.
# I4 = SpaceX Inspiration4 mission (3-day, 4 civilian crew).
# Data source: SpaceOmicsBench 2025_01_08_v2.cfRNA (cfrna_crossmission_r1.csv)
#
# Input:  cfrna_crossmission_r1.csv — 5,346 common genes (I4 ∩ Polaris Dawn)
#         Column: log2FoldChange_I4 (FP1 vs pre-flight, DESeq2)
# Output: v2/processed/E2_human/i4_cfrna_hallmark_fgsea.csv
#
# Usage:
#   Rscript v2/scripts/run_i4_cfrna_fgsea.R

suppressPackageStartupMessages({
  library(fgsea)
  library(data.table)
})

# ── Paths ──────────────────────────────────────────────────────────────────────
REPO_ROOT <- normalizePath(
  file.path(dirname(commandArgs(trailingOnly = FALSE)[grep("--file=", commandArgs(trailingOnly = FALSE))]), "../.."),
  mustWork = FALSE
)
if (!dir.exists(REPO_ROOT)) {
  REPO_ROOT <- normalizePath(".", mustWork = TRUE)
}

I4_PATH <- file.path(
  "/Users/jak4013/Dropbox/Bioinformatics/Claude/SpaceOmicsBench",
  "2025_01_08_v2.cfRNA/SpaceOmicsBench_v2.0/data/cfrna_crossmission_r1.csv"
)
OUT_DIR  <- file.path(REPO_ROOT, "v2/processed/E2_human")
OUT_FILE <- file.path(OUT_DIR, "i4_cfrna_hallmark_fgsea.csv")

dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

cat("=== I4 cfRNA fGSEA (Hallmark, hsapiens) ===\n")
cat(sprintf("Input: %s\n", I4_PATH))

# ── Load I4 cfRNA DE ───────────────────────────────────────────────────────────
cat("Loading I4 cfRNA DE...\n")
de <- fread(I4_PATH)
cat(sprintf("  %d genes loaded\n", nrow(de)))

# Ranking column: log2FoldChange_I4 (positive = up in spaceflight)
rank_col <- "log2FoldChange_I4"
gene_col <- "gene_name"
if (!rank_col %in% names(de)) {
  stop(sprintf("Column '%s' not found. Available: %s", rank_col, paste(names(de), collapse=", ")))
}

# Create named ranking vector
ranks <- de[[rank_col]]
names(ranks) <- de[[gene_col]]

# Remove NA and zero values, deduplicate gene names (keep highest |LFC|)
ranks <- ranks[!is.na(ranks) & ranks != 0]
# Handle duplicates: keep max |LFC|
if (any(duplicated(names(ranks)))) {
  df_tmp <- data.frame(gene=names(ranks), lfc=ranks)
  df_tmp <- df_tmp[order(-abs(df_tmp$lfc)), ]
  df_tmp <- df_tmp[!duplicated(df_tmp$gene), ]
  ranks <- setNames(df_tmp$lfc, df_tmp$gene)
}

# Sort descending
ranks <- sort(ranks, decreasing = TRUE)
cat(sprintf("  Ranked genes: %d\n", length(ranks)))
cat(sprintf("  Score range: [%.4f, %.4f]\n", min(ranks), max(ranks)))

# ── Load Hallmark gene sets (Homo sapiens) from GMT ────────────────────────────
GMT_PATH <- "/tmp/h.all.human.symbols.gmt"
if (!file.exists(GMT_PATH)) {
  cat("Downloading Hallmark GMT from MSigDB...\n")
  download.file(
    "https://data.broadinstitute.org/gsea-msigdb/msigdb/release/7.5.1/h.all.v7.5.1.symbols.gmt",
    destfile = GMT_PATH, quiet = TRUE, timeout = 30
  )
}
cat(sprintf("Loading Hallmark GMT from: %s\n", GMT_PATH))
hallmark_list <- gmtPathways(GMT_PATH)
cat(sprintf("  %d gene sets loaded\n", length(hallmark_list)))

# Check overlap: I4 has only 5,346 genes, so some pathways may be small
overlaps <- sapply(hallmark_list, function(gs) sum(names(ranks) %in% gs))
cat(sprintf("  Pathway gene overlap: min=%d, mean=%.1f, max=%d\n",
            min(overlaps), mean(overlaps), max(overlaps)))
cat(sprintf("  Pathways with overlap <10: %d\n", sum(overlaps < 10)))
small_pw <- names(overlaps[overlaps < 10])
if (length(small_pw) > 0) {
  cat(sprintf("    Excluded (overlap<10): %s\n", paste(small_pw, collapse=", ")))
}

# ── Run fGSEA ─────────────────────────────────────────────────────────────────
# Use minSize=10 (not 15) because I4 has only 5,346 genes
# This allows 49/50 pathways to be computed (only ANGIOGENESIS excluded with 8 genes)
cat("Running fGSEA (minSize=10)...\n")
set.seed(42)
fgsea_res <- fgsea(
  pathways  = hallmark_list,
  stats     = ranks,
  minSize   = 10,
  maxSize   = 500,
  nPermSimple = 10000
)

fgsea_res <- as.data.table(fgsea_res)
fgsea_res$leadingEdge_str <- sapply(fgsea_res$leadingEdge, paste, collapse="; ")
fgsea_res$leadingEdge <- NULL
fgsea_res$tissue <- "blood_cfrna"
fgsea_res$mission <- "I4"
fgsea_res$species <- "human"
fgsea_res$ranking_metric <- rank_col
fgsea_res$n_ranked_genes <- length(ranks)

cat(sprintf("  fGSEA complete: %d pathways\n", nrow(fgsea_res)))

# Print top results
cat("\nTop pathways (by |NES|):\n")
top <- fgsea_res[order(-abs(NES))][1:min(10, nrow(fgsea_res))]
for (i in seq_len(nrow(top))) {
  cat(sprintf("  %-50s NES=%6.3f padj=%.4f size=%d\n",
              top$pathway[i], top$NES[i], top$padj[i], top$size[i]))
}

cat(sprintf("\nSig pathways (padj<0.05): %d\n", sum(fgsea_res$padj < 0.05, na.rm=TRUE)))

# ── Save ──────────────────────────────────────────────────────────────────────
fwrite(fgsea_res, OUT_FILE)
cat(sprintf("\nSaved: %s\n", OUT_FILE))
cat("=== Done ===\n")
