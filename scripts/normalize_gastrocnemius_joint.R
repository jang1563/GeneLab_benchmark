#!/usr/bin/env Rscript
# normalize_gastrocnemius_joint.R
#
# Joint DESeq2 normalization of all 3 gastrocnemius missions (RR-1, RR-5, RR-9)
# from raw (RSEM unnormalized) counts, using a shared gene universe.
#
# Rationale: Previous per-study normalization created systematic inter-mission
# differences (RR-9 zero rate 30% vs RR-5 ~1%). Joint normalization on a
# common gene set removes this artifact.
#
# Method: DESeq2 estimateSizeFactors on all 32 samples together.
# No ComBat-seq batch correction — just normalize, then LOMO uses
# variance filtering per fold for feature selection.
#
# Usage: Rscript scripts/normalize_gastrocnemius_joint.R

suppressPackageStartupMessages({
  library(DESeq2)
})

# ── Paths ────────────────────────────────────────────────────────────────────────
args_all  <- commandArgs(trailingOnly = FALSE)
file_args <- args_all[startsWith(args_all, "--file=")]
if (length(file_args) > 0) {
  script_file <- sub("^--file=", "", file_args[1])
  base_dir    <- dirname(dirname(normalizePath(script_file)))
} else {
  base_dir <- normalizePath(".")
}

data_dir  <- file.path(base_dir, "data", "mouse", "gastrocnemius")
out_dir   <- file.path(base_dir, "processed", "A_detection", "gastrocnemius")

missions <- list(
  "RR-1" = list(
    counts = file.path(data_dir, "RR-1",
      "GLDS-101_rna_seq_RSEM_Unnormalized_Counts_GLbulkRNAseq.csv"),
    runsheet = file.path(data_dir, "RR-1",
      "GLDS-101_rna_seq_bulkRNASeq_v1_runsheet.csv")
  ),
  "RR-5" = list(
    counts = file.path(data_dir, "RR-5",
      "GLDS-401_rna_seq_RSEM_Unnormalized_Counts_GLbulkRNAseq.csv"),
    runsheet = file.path(data_dir, "RR-5",
      "GLDS-401_rna_seq_bulkRNASeq_v1_runsheet.csv")
  ),
  "RR-9" = list(
    counts = file.path(data_dir, "RR-9",
      "GLDS-326_rna_seq_RSEM_Unnormalized_Counts_GLbulkRNAseq.csv"),
    runsheet = file.path(data_dir, "RR-9",
      "GLDS-326_rna_seq_bulkRNASeq_v2_runsheet.csv")
  )
)

cat("=== Joint DESeq2 Normalization: Gastrocnemius (RR-1 + RR-5 + RR-9) ===\n\n")

# ── Load per-mission counts & labels ─────────────────────────────────────────────
count_list <- list()
coldata_list <- list()

for (mission in names(missions)) {
  info <- missions[[mission]]
  cat(sprintf("Loading %s...\n", mission))

  counts_raw <- read.csv(info$counts, row.names = 1, check.names = FALSE)
  rs <- read.csv(info$runsheet, check.names = FALSE)

  # Build coldata
  cond <- ifelse(rs[["Factor Value[Spaceflight]"]] == "Space Flight", "Flight", "GC")
  cd <- data.frame(
    sample  = rs[["Sample Name"]],
    mission = mission,
    condition = cond,
    row.names = rs[["Sample Name"]],
    stringsAsFactors = FALSE
  )

  # Align samples
  common_samples <- intersect(colnames(counts_raw), rownames(cd))
  counts_raw <- counts_raw[, common_samples]
  cd <- cd[common_samples, , drop = FALSE]

  cat(sprintf("  %s: %d genes × %d samples (Flight=%d, GC=%d)\n",
    mission, nrow(counts_raw), ncol(counts_raw),
    sum(cd$condition == "Flight"), sum(cd$condition == "GC")))

  count_list[[mission]] <- counts_raw
  coldata_list[[mission]] <- cd
}

# ── Find shared gene universe ─────────────────────────────────────────────────────
gene_sets <- lapply(count_list, rownames)
common_genes <- Reduce(intersect, gene_sets)
cat(sprintf("\nCommon genes across all 3 missions: %d\n", length(common_genes)))
cat(sprintf("  RR-1: %d, RR-5: %d, RR-9: %d\n",
  nrow(count_list[["RR-1"]]), nrow(count_list[["RR-5"]]),
  nrow(count_list[["RR-9"]])))

# Restrict to common genes
for (m in names(count_list)) {
  count_list[[m]] <- count_list[[m]][common_genes, ]
}

# ── Merge into single count matrix ───────────────────────────────────────────────
counts_merged <- do.call(cbind, count_list)
coldata_merged <- do.call(rbind, coldata_list)
cat(sprintf("\nMerged matrix: %d genes × %d samples\n",
  nrow(counts_merged), ncol(counts_merged)))
cat("Sample breakdown:\n")
print(table(coldata_merged$mission, coldata_merged$condition))

# ── DESeq2 joint normalization ────────────────────────────────────────────────────
cat("\nRunning DESeq2 on merged counts...\n")

counts_int <- round(as.matrix(counts_merged))
storage.mode(counts_int) <- "integer"

# Use mission + condition to prevent confounding
# But estimate size factors on all samples jointly
dds <- DESeqDataSetFromMatrix(
  countData = counts_int,
  colData   = coldata_merged,
  design    = ~ mission + condition
)

# Low-count filter: ≥10 counts in ≥2 samples (applied globally)
keep <- rowSums(counts(dds) >= 10) >= 2
dds  <- dds[keep, ]
cat(sprintf("  Genes after low-count filter: %d / %d\n", sum(keep), length(keep)))

# Size factor estimation
dds <- estimateSizeFactors(dds)
cat("  Size factors per mission:\n")
for (m in names(missions)) {
  sf <- sizeFactors(dds)[coldata_merged$mission == m]
  cat(sprintf("    %s: min=%.3f, max=%.3f, mean=%.3f\n",
    m, min(sf), max(sf), mean(sf)))
}

# Normalized counts
norm_counts <- counts(dds, normalized = TRUE)
log2_norm   <- log2(norm_counts + 1)

cat(sprintf("\nNormalized: %d genes × %d samples\n", nrow(log2_norm), ncol(log2_norm)))

# ── Zero rate check ───────────────────────────────────────────────────────────────
cat("\nZero rate check (after joint normalization):\n")
for (m in names(missions)) {
  idx <- coldata_merged$mission == m
  vals <- log2_norm[, idx]
  cat(sprintf("  %s: zero rate = %.1f%%, median = %.2f\n",
    m, 100 * mean(vals == 0), median(vals)))
}

# ── Save outputs ──────────────────────────────────────────────────────────────────
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

# Transpose to samples × genes
log2_norm_t <- as.data.frame(t(log2_norm))

out_norm <- file.path(out_dir, "gastrocnemius_joint_log2_norm.csv")
write.csv(log2_norm_t, out_norm)
cat(sprintf("\n✓ Saved joint log2_norm: %s\n", out_norm))

# Save metadata
meta_out <- data.frame(
  label     = coldata_merged$condition,
  mission   = coldata_merged$mission,
  REMOVE    = FALSE,
  tissue    = "gastrocnemius",
  row.names = rownames(coldata_merged),
  stringsAsFactors = FALSE
)
out_meta <- file.path(out_dir, "gastrocnemius_joint_metadata.csv")
write.csv(meta_out, out_meta)
cat(sprintf("✓ Saved joint metadata: %s\n", out_meta))

cat(sprintf("\nSummary: %d samples × %d genes\n", nrow(log2_norm_t), ncol(log2_norm_t)))
cat(sprintf("  RR-1: %d, RR-5: %d, RR-9: %d\n",
  sum(meta_out$mission == "RR-1"),
  sum(meta_out$mission == "RR-5"),
  sum(meta_out$mission == "RR-9")))
