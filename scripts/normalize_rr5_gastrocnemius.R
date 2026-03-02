#!/usr/bin/env Rscript
# normalize_rr5_gastrocnemius.R
#
# DESeq2 normalization for RR-5 (GLDS-401) gastrocnemius data.
# Produces log2(DESeq2_normalized + 1) matrix in sample × gene format
# matching the existing gastrocnemius_all_missions_log2_norm.csv.
#
# Usage: Rscript scripts/normalize_rr5_gastrocnemius.R

suppressPackageStartupMessages({
  library(DESeq2)
})

# ── Paths ───────────────────────────────────────────────────────────────────────
args_all  <- commandArgs(trailingOnly = FALSE)
file_args <- args_all[startsWith(args_all, "--file=")]
if (length(file_args) > 0) {
  script_file <- sub("^--file=", "", file_args[1])
  base_dir    <- dirname(dirname(normalizePath(script_file)))
} else {
  base_dir <- normalizePath(".")
}

rr5_dir   <- file.path(base_dir, "data", "mouse", "gastrocnemius", "RR-5")
out_dir   <- file.path(base_dir, "processed", "A_detection", "gastrocnemius")

counts_file   <- file.path(rr5_dir, "GLDS-401_rna_seq_RSEM_Unnormalized_Counts_GLbulkRNAseq.csv")
runsheet_file <- file.path(rr5_dir, "GLDS-401_rna_seq_bulkRNASeq_v1_runsheet.csv")
out_file      <- file.path(out_dir, "gastrocnemius_RR-5_log2_norm.csv")
meta_file     <- file.path(out_dir, "gastrocnemius_RR-5_metadata.csv")

cat("=== DESeq2 Normalization: Gastrocnemius RR-5 (GLDS-401) ===\n")
cat(sprintf("Input counts: %s\n", counts_file))
cat(sprintf("Output: %s\n", out_file))

# ── Load counts ─────────────────────────────────────────────────────────────────
cat("\nLoading raw counts...\n")
counts_raw <- read.csv(counts_file, row.names = 1, check.names = FALSE)
cat(sprintf("  Loaded: %d genes × %d samples\n", nrow(counts_raw), ncol(counts_raw)))

# ── Load runsheet ────────────────────────────────────────────────────────────────
cat("Loading runsheet...\n")
runsheet <- read.csv(runsheet_file, check.names = FALSE)
cat(sprintf("  Samples: %d\n", nrow(runsheet)))

# Build coldata
coldata <- data.frame(
  sample_name = runsheet[["Sample Name"]],
  condition   = runsheet[["Factor Value[Spaceflight]"]],
  row.names   = runsheet[["Sample Name"]],
  stringsAsFactors = FALSE
)
coldata$condition <- ifelse(coldata$condition == "Space Flight", "Flight", "GC")
cat("\nSample groups:\n")
print(table(coldata$condition))

# Ensure samples match between counts and coldata
common_samples <- intersect(colnames(counts_raw), rownames(coldata))
cat(sprintf("\nMatching samples: %d\n", length(common_samples)))
counts_raw <- counts_raw[, common_samples]
coldata    <- coldata[common_samples, , drop = FALSE]

# ── DESeq2 normalization ─────────────────────────────────────────────────────────
cat("\nRunning DESeq2 size factor estimation...\n")

# Round to integer (RSEM may produce non-integer counts)
counts_int <- round(as.matrix(counts_raw))
storage.mode(counts_int) <- "integer"

dds <- DESeqDataSetFromMatrix(
  countData = counts_int,
  colData   = coldata,
  design    = ~ condition
)

# Filter low-count genes (keep genes with ≥ 10 counts in ≥ 2 samples)
keep <- rowSums(counts(dds) >= 10) >= 2
dds  <- dds[keep, ]
cat(sprintf("  Genes after low-count filter: %d / %d\n", sum(keep), length(keep)))

# Estimate size factors
dds <- estimateSizeFactors(dds)
cat("  Size factors:\n")
sf <- sizeFactors(dds)
for (nm in names(sf)) cat(sprintf("    %s: %.4f\n", nm, sf[nm]))

# Extract normalized counts
norm_counts <- counts(dds, normalized = TRUE)
cat(sprintf("\n  Normalized counts: %d genes × %d samples\n", nrow(norm_counts), ncol(norm_counts)))

# Log2 transform with +1 pseudocount
log2_norm <- log2(norm_counts + 1)

# Transpose to samples × genes (matching existing format)
log2_norm_t <- as.data.frame(t(log2_norm))
cat(sprintf("  Output matrix: %d samples × %d genes\n", nrow(log2_norm_t), ncol(log2_norm_t))  )

# ── Save outputs ─────────────────────────────────────────────────────────────────
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

write.csv(log2_norm_t, out_file)
cat(sprintf("\n✓ Saved log2_norm matrix: %s\n", out_file))

# Save metadata
meta_out <- data.frame(
  sample_id    = rownames(log2_norm_t),
  label        = coldata[rownames(log2_norm_t), "condition"],
  mission      = "RR-5",
  osd_id       = "OSD-401",
  tissue       = "gastrocnemius",
  stringsAsFactors = FALSE
)
write.csv(meta_out, meta_file, row.names = FALSE)
cat(sprintf("✓ Saved metadata: %s\n", meta_file))

cat(sprintf("\nSummary:\n  n_samples = %d (Flight=%d, GC=%d)\n  n_genes = %d\n",
    nrow(log2_norm_t),
    sum(meta_out$label == "Flight"),
    sum(meta_out$label == "GC"),
    ncol(log2_norm_t)))
