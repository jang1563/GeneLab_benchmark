#!/usr/bin/env Rscript
# run_dge_comparison.R — GeneLab_benchmark: J2 DGE Pipeline Comparison
#
# Runs DESeq2, edgeR QLF, and limma-voom on raw counts to compare DGE pipelines.
# Uses the same FLT vs GC contrast across all three methods.
#
# Scope: liver (6 missions) + thymus (3 missions) = 9 × 3 = 27 DGE runs
# Skin excluded: RR-7 has no raw counts available.
#
# Output: processed/J2_dge_comparison/{tissue}/{mission}_{pipeline}_dge.csv
#   Columns: ENSEMBL, log2FC, stat, pvalue, adj_pvalue
#
# Usage:
#   Rscript scripts/run_dge_comparison.R --tissue liver --mission RR-1
#   Rscript scripts/run_dge_comparison.R --tissue liver --all
#   Rscript scripts/run_dge_comparison.R --all-tissues

# ── Package loading ────────────────────────────────────────────────────────────
cat("Loading packages...\n")
suppressPackageStartupMessages({
  if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager", repos = "https://cloud.r-project.org")
  for (pkg in c("DESeq2", "edgeR", "limma")) {
    if (!requireNamespace(pkg, quietly = TRUE))
      BiocManager::install(pkg, ask = FALSE, update = FALSE)
  }
  library(DESeq2)
  library(edgeR)
  library(limma)
  library(argparse)
})

# ── Path setup ─────────────────────────────────────────────────────────────────
args_all  <- commandArgs(trailingOnly = FALSE)
file_args <- args_all[startsWith(args_all, "--file=")]
if (length(file_args) > 0) {
  script_file <- sub("^--file=", "", file_args[1])
  BASE_DIR    <- dirname(dirname(normalizePath(script_file)))
} else {
  BASE_DIR <- normalizePath(".")
}

DATA_DIR <- file.path(BASE_DIR, "data", "mouse")
OUT_DIR  <- file.path(BASE_DIR, "processed", "J2_dge_comparison")

# ── CLI ────────────────────────────────────────────────────────────────────────
parser <- ArgumentParser(description = "J2: DGE Pipeline Comparison (DESeq2/edgeR/limma-voom)")
parser$add_argument("--tissue", default = NULL,
                    help = "Tissue to process (liver or thymus)")
parser$add_argument("--mission", default = NULL,
                    help = "Single mission to process (e.g., RR-1)")
parser$add_argument("--all", action = "store_true",
                    help = "Run all missions for the specified tissue")
parser$add_argument("--all-tissues", action = "store_true",
                    help = "Run all tissues and all missions")
parser$add_argument("--pipeline", default = NULL,
                    choices = c("deseq2", "edger", "limma_voom"),
                    help = "Run only one pipeline (default: all three)")
args <- parser$parse_args()

# ── Tissue-Mission mapping (liver + thymus only) ──────────────────────────────
TISSUE_MISSIONS <- list(
  liver = list(
    list(mission = "RR-1",  dir = "RR-1",  glds = "GLDS-48"),
    list(mission = "RR-3",  dir = "RR-3",  glds = "GLDS-137"),
    list(mission = "RR-6",  dir = "RR-6",  glds = "GLDS-245"),
    list(mission = "RR-8",  dir = "RR-8",  glds = "GLDS-379"),
    list(mission = "RR-9",  dir = "RR-9",  glds = "GLDS-242"),
    list(mission = "MHU-2", dir = "MHU-2", glds = "GLDS-617")
  ),
  thymus = list(
    list(mission = "RR-6",  dir = "RR-6",  glds = "GLDS-244"),
    list(mission = "MHU-2", dir = "MHU-2", glds = "GLDS-289"),
    list(mission = "RR-9",  dir = "RR-9",  glds = "GLDS-421")
  )
)

PIPELINES <- c("deseq2", "edger", "limma_voom")

# ── Helper: find raw counts file ──────────────────────────────────────────────
find_counts_file <- function(tissue, mission_dir, glds) {
  data_path <- file.path(DATA_DIR, tissue, mission_dir)

  # Try new naming convention first
  new_pattern <- sprintf("%s_rna_seq_RSEM_Unnormalized_Counts_GLbulkRNAseq.csv", glds)
  new_file <- file.path(data_path, new_pattern)
  if (file.exists(new_file)) return(new_file)

  # Try old naming convention
  old_pattern <- sprintf("%s_rna_seq_Unnormalized_Counts.csv", glds)
  old_file <- file.path(data_path, old_pattern)
  if (file.exists(old_file)) return(old_file)

  stop(sprintf("No raw counts file found for %s/%s (%s)", tissue, mission_dir, glds))
}

# ── Helper: find SampleTable file ─────────────────────────────────────────────
find_sample_table <- function(tissue, mission_dir, glds) {
  data_path <- file.path(DATA_DIR, tissue, mission_dir)

  new_file <- file.path(data_path,
    sprintf("%s_rna_seq_SampleTable_GLbulkRNAseq.csv", glds))
  if (file.exists(new_file)) return(new_file)

  old_file <- file.path(data_path,
    sprintf("%s_rna_seq_SampleTable.csv", glds))
  if (file.exists(old_file)) return(old_file)

  stop(sprintf("No SampleTable found for %s/%s (%s)", tissue, mission_dir, glds))
}

# ── Helper: classify condition to FLT/GC binary ──────────────────────────────
# Returns "FLT", "GC", or NA (to exclude)
classify_condition <- function(cond) {
  # Exclude centrifuge / artificial gravity
  if (grepl("1G\\.by\\.centrifug|1G.by.centrifug", cond, ignore.case = TRUE)) return(NA)

  # Space Flight (uG or general) → FLT
  # Also handles MHU-prefixed conditions like "MHU.2...Space.Flight...uG"
  if (grepl("Space\\.Flight|Space Flight", cond, ignore.case = TRUE)) return("FLT")

  # Ground Control → GC (anywhere in string for MHU-prefixed conditions)
  if (grepl("Ground\\.Control|Ground Control", cond, ignore.case = TRUE)) return("GC")

  # Vivarium Control → GC (used as control in some missions like MHU-2 thymus)
  if (grepl("Vivarium\\.Control|Vivarium Control", cond, ignore.case = TRUE)) return("GC")

  # Abbreviated forms (RR-9 style)
  if (grepl("^FLT", cond)) return("FLT")
  if (grepl("^GC", cond)) return("GC")
  # CC (cage control) maps to GC in RR-9
  if (grepl("^CC", cond)) return("GC")

  # Everything else (BSL, VIV_C2, etc.) → exclude
  return(NA)
}

# ── Helper: load and filter data ──────────────────────────────────────────────
load_data <- function(tissue, mission_dir, glds) {
  counts_file <- find_counts_file(tissue, mission_dir, glds)
  sample_file <- find_sample_table(tissue, mission_dir, glds)

  cat(sprintf("  Counts: %s\n", basename(counts_file)))
  cat(sprintf("  Samples: %s\n", basename(sample_file)))

  # Read counts (genes × samples, first col = ENSMUSG IDs as row names)
  counts <- read.csv(counts_file, row.names = 1, check.names = FALSE)

  # Read sample table
  st <- read.csv(sample_file, row.names = 1, check.names = FALSE)

  # Classify conditions
  st$group <- sapply(st$condition, classify_condition)

  # Filter to FLT/GC only
  keep_samples <- rownames(st)[!is.na(st$group)]
  st <- st[keep_samples, , drop = FALSE]

  # Align counts with filtered samples
  common_samples <- intersect(colnames(counts), keep_samples)
  if (length(common_samples) < 4) {
    stop(sprintf("Only %d FLT/GC samples found (need >= 4)", length(common_samples)))
  }

  counts <- counts[, common_samples, drop = FALSE]
  st <- st[common_samples, , drop = FALSE]

  # Ensure integer counts (some files have decimal RSEM values)
  counts <- round(counts)

  # Filter low-expression genes: keep genes with >= 10 counts in >= min_group_size samples
  min_group <- min(table(st$group))
  keep_genes <- rowSums(counts >= 10) >= min_group
  counts <- counts[keep_genes, , drop = FALSE]

  # Remove non-ENSMUSG rows (safety check for metadata contamination)
  ensmusg_rows <- grepl("^ENSMUSG", rownames(counts))
  counts <- counts[ensmusg_rows, , drop = FALSE]

  cat(sprintf("  Samples: %d FLT, %d GC\n",
              sum(st$group == "FLT"), sum(st$group == "GC")))
  cat(sprintf("  Genes after filtering: %d\n", nrow(counts)))

  list(counts = counts, meta = st)
}

# ── Pipeline: DESeq2 Wald test ────────────────────────────────────────────────
run_deseq2 <- function(counts, meta) {
  meta$group <- factor(meta$group, levels = c("GC", "FLT"))
  dds <- DESeqDataSetFromMatrix(countData = counts,
                                 colData = meta,
                                 design = ~ group)
  dds <- DESeq(dds, quiet = TRUE)
  res <- results(dds, contrast = c("group", "FLT", "GC"))

  data.frame(
    ENSEMBL = rownames(res),
    log2FC = res$log2FoldChange,
    stat = res$stat,
    pvalue = res$pvalue,
    adj_pvalue = res$padj,
    row.names = NULL,
    stringsAsFactors = FALSE
  )
}

# ── Pipeline: edgeR quasi-likelihood F-test ───────────────────────────────────
run_edger <- function(counts, meta) {
  meta$group <- factor(meta$group, levels = c("GC", "FLT"))
  dge <- DGEList(counts = counts, group = meta$group)
  dge <- calcNormFactors(dge)
  design <- model.matrix(~ group, data = meta)
  dge <- estimateDisp(dge, design)
  fit <- glmQLFit(dge, design)
  qlf <- glmQLFTest(fit, coef = "groupFLT")
  tt <- topTags(qlf, n = Inf, sort.by = "none")$table

  data.frame(
    ENSEMBL = rownames(tt),
    log2FC = tt$logFC,
    stat = tt$F,
    pvalue = tt$PValue,
    adj_pvalue = tt$FDR,
    row.names = NULL,
    stringsAsFactors = FALSE
  )
}

# ── Pipeline: limma-voom with quality weights ─────────────────────────────────
run_limma_voom <- function(counts, meta) {
  meta$group <- factor(meta$group, levels = c("GC", "FLT"))
  dge <- DGEList(counts = counts, group = meta$group)
  dge <- calcNormFactors(dge)
  design <- model.matrix(~ group, data = meta)

  # voomWithQualityWeights for sample quality weighting
  v <- tryCatch(
    voomWithQualityWeights(dge, design, plot = FALSE),
    error = function(e) {
      cat("    voomWithQualityWeights failed, falling back to voom\n")
      voom(dge, design, plot = FALSE)
    }
  )

  fit <- lmFit(v, design)
  fit <- eBayes(fit)
  tt <- topTable(fit, coef = "groupFLT", number = Inf, sort.by = "none")

  data.frame(
    ENSEMBL = rownames(tt),
    log2FC = tt$logFC,
    stat = tt$t,
    pvalue = tt$P.Value,
    adj_pvalue = tt$adj.P.Val,
    row.names = NULL,
    stringsAsFactors = FALSE
  )
}

# ── Main runner ───────────────────────────────────────────────────────────────
run_mission <- function(tissue, mission_info, pipelines) {
  mission <- mission_info$mission
  mission_dir <- mission_info$dir
  glds <- mission_info$glds

  cat(sprintf("\n=== %s / %s (%s) ===\n", tissue, mission, glds))

  # Load and filter data
  data <- tryCatch(
    load_data(tissue, mission_dir, glds),
    error = function(e) {
      cat(sprintf("  ERROR: %s\n", e$message))
      return(NULL)
    }
  )
  if (is.null(data)) return(NULL)

  # Create output directory
  out_path <- file.path(OUT_DIR, tissue)
  dir.create(out_path, recursive = TRUE, showWarnings = FALSE)

  results <- list()
  for (pipeline in pipelines) {
    cat(sprintf("  Running %s...", pipeline))
    t0 <- proc.time()

    res <- tryCatch({
      switch(pipeline,
        deseq2     = run_deseq2(data$counts, data$meta),
        edger      = run_edger(data$counts, data$meta),
        limma_voom = run_limma_voom(data$counts, data$meta)
      )
    }, error = function(e) {
      cat(sprintf(" ERROR: %s\n", e$message))
      return(NULL)
    })

    if (is.null(res)) next

    elapsed <- (proc.time() - t0)["elapsed"]
    n_sig <- sum(res$adj_pvalue < 0.05, na.rm = TRUE)
    cat(sprintf(" done (%.1fs, %d DEGs at FDR<0.05)\n", elapsed, n_sig))

    # Write output
    out_file <- file.path(out_path, sprintf("%s_%s_dge.csv", mission, pipeline))
    write.csv(res, out_file, row.names = FALSE)
    cat(sprintf("    → %s\n", out_file))

    results[[pipeline]] <- list(
      n_genes = nrow(res),
      n_deg_005 = n_sig,
      n_deg_001 = sum(res$adj_pvalue < 0.01, na.rm = TRUE),
      elapsed = elapsed
    )
  }

  # Print comparison summary
  if (length(results) > 1) {
    cat("\n  Pipeline comparison:\n")
    cat(sprintf("  %-12s %7s %10s %10s\n", "Pipeline", "Genes", "DEG(0.05)", "DEG(0.01)"))
    for (p in names(results)) {
      r <- results[[p]]
      cat(sprintf("  %-12s %7d %10d %10d\n", p, r$n_genes, r$n_deg_005, r$n_deg_001))
    }
  }

  return(results)
}

# ── Main ──────────────────────────────────────────────────────────────────────
main <- function() {
  pipelines <- if (!is.null(args$pipeline)) args$pipeline else PIPELINES

  if (args$all_tissues) {
    # Run all tissues, all missions
    all_results <- list()
    for (tissue in names(TISSUE_MISSIONS)) {
      for (mi in TISSUE_MISSIONS[[tissue]]) {
        key <- sprintf("%s/%s", tissue, mi$mission)
        all_results[[key]] <- run_mission(tissue, mi, pipelines)
      }
    }
  } else if (!is.null(args$tissue)) {
    tissue <- args$tissue
    if (!(tissue %in% names(TISSUE_MISSIONS))) {
      stop(sprintf("Unknown tissue: %s. Available: %s",
                   tissue, paste(names(TISSUE_MISSIONS), collapse = ", ")))
    }

    if (args$all) {
      # All missions for this tissue
      for (mi in TISSUE_MISSIONS[[tissue]]) {
        run_mission(tissue, mi, pipelines)
      }
    } else if (!is.null(args$mission)) {
      # Single mission
      mi_match <- Filter(function(x) x$mission == args$mission, TISSUE_MISSIONS[[tissue]])
      if (length(mi_match) == 0) {
        stop(sprintf("Mission %s not found for tissue %s", args$mission, tissue))
      }
      run_mission(tissue, mi_match[[1]], pipelines)
    } else {
      stop("Specify --mission or --all")
    }
  } else {
    stop("Specify --tissue (with --mission or --all) or --all-tissues")
  }

  cat("\nDone.\n")
}

main()
