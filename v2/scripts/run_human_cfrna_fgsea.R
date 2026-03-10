#!/usr/bin/env Rscript
# run_human_cfrna_fgsea.R — GeneLab_benchmark v2, Category E1
#
# Runs fGSEA on JAXA CFE human cfRNA DE data against MSigDB Hallmark gene sets
# (Homo sapiens).
#
# Input:  cfrna_3group_de_noleak.csv (SpaceOmicsBench)
# Output: v2/processed/E1_human_cfrna_hallmark_fgsea.csv
#
# Usage:
#   Rscript v2/scripts/run_human_cfrna_fgsea.R

suppressPackageStartupMessages({
  library(fgsea)
  library(data.table)
})

# ── Paths ──────────────────────────────────────────────────────────────────────
REPO_ROOT <- normalizePath(
  file.path(dirname(commandArgs(trailingOnly = FALSE)[grep("--file=", commandArgs(trailingOnly = FALSE))]), "../.."),
  mustWork = FALSE
)
# Fallback: use current working directory
if (!dir.exists(REPO_ROOT)) {
  REPO_ROOT <- normalizePath(".", mustWork = TRUE)
}

CFRNA_PATH <- file.path(
  "/Users/jak4013/Dropbox/Bioinformatics/Claude/SpaceOmicsBench",
  "v2_public/data/processed/cfrna_3group_de_noleak.csv"
)
OUT_DIR <- file.path(REPO_ROOT, "v2/processed/E1_human")
OUT_FILE <- file.path(OUT_DIR, "human_cfrna_hallmark_fgsea.csv")

dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

cat("=== Human cfRNA fGSEA (Hallmark, hsapiens) ===\n")
cat(sprintf("Input: %s\n", CFRNA_PATH))

# ── Load cfRNA DE ──────────────────────────────────────────────────────────────
cat("Loading human cfRNA DE...\n")
de <- fread(CFRNA_PATH)
cat(sprintf("  %d genes loaded\n", nrow(de)))

# Use edge_pre_vs_flight_diff as ranking: positive = up in flight
rank_col <- "edge_pre_vs_flight_diff"
if (!rank_col %in% names(de)) {
  stop(sprintf("Column '%s' not found. Available: %s", rank_col, paste(names(de), collapse=", ")))
}

# Create named ranking vector (gene symbol → diff score)
ranks <- de[[rank_col]]
names(ranks) <- de[["gene"]]

# Remove NA and zero values
ranks <- ranks[!is.na(ranks) & ranks != 0]

# Sort descending
ranks <- sort(ranks, decreasing = TRUE)
cat(sprintf("  Ranked genes: %d\n", length(ranks)))
cat(sprintf("  Score range: [%.4f, %.4f]\n", min(ranks), max(ranks)))

# ── Load Hallmark gene sets (Homo sapiens) from GMT ────────────────────────────
# Uses downloaded GMT file: MSigDB v7.5.1 Hallmark gene sets for Homo sapiens
GMT_PATH <- "/tmp/h.all.human.symbols.gmt"
if (!file.exists(GMT_PATH)) {
  # Try to download if not present
  download.file(
    "https://data.broadinstitute.org/gsea-msigdb/msigdb/release/7.5.1/h.all.v7.5.1.symbols.gmt",
    destfile = GMT_PATH, quiet = TRUE, timeout = 30
  )
}
cat(sprintf("Loading Hallmark GMT from: %s\n", GMT_PATH))
hallmark_list <- gmtPathways(GMT_PATH)
cat(sprintf("  %d gene sets loaded\n", length(hallmark_list)))

# ── Run fGSEA ─────────────────────────────────────────────────────────────────
cat("Running fGSEA...\n")
set.seed(42)
fgsea_res <- fgsea(
  pathways = hallmark_list,
  stats = ranks,
  minSize = 15,
  maxSize = 500,
  nPermSimple = 10000
)

fgsea_res <- as.data.table(fgsea_res)
fgsea_res$leadingEdge_str <- sapply(fgsea_res$leadingEdge, paste, collapse="; ")
fgsea_res$leadingEdge <- NULL  # drop list column for CSV
fgsea_res$tissue <- "blood_cfrna"
fgsea_res$mission <- "JAXA_CFE"
fgsea_res$species <- "human"
fgsea_res$ranking_metric <- rank_col

cat(sprintf("  fGSEA complete: %d pathways\n", nrow(fgsea_res)))

# Print top results
cat("\nTop pathways (by |NES|):\n")
top <- fgsea_res[order(-abs(NES))][1:10]
for (i in seq_len(nrow(top))) {
  cat(sprintf("  %-50s NES=%6.3f padj=%.4f\n",
              top$pathway[i], top$NES[i], top$padj[i]))
}

# ── Save ──────────────────────────────────────────────────────────────────────
fwrite(fgsea_res[, !"leadingEdge", with = FALSE], OUT_FILE)
cat(sprintf("\nSaved: %s\n", OUT_FILE))
cat("=== Done ===\n")
