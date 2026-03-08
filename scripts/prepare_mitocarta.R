#!/usr/bin/env Rscript
# prepare_mitocarta.R — Download MitoCarta3.0 and convert to GMT format
#
# Downloads Mouse.MitoCarta3.0.xls from Broad Institute,
# parses MitoPathways annotations, and saves as GMT for fGSEA/GSVA.
#
# Output: processed/gene_sets/mitocarta3_mouse.gmt
#
# Usage:
#   Rscript scripts/prepare_mitocarta.R
#   Rscript scripts/prepare_mitocarta.R --force   # re-download even if cached

suppressPackageStartupMessages({
  library(argparse)
})

# ── CLI ─────────────────────────────────────────────────────────────────────────
parser <- ArgumentParser(description = "Download MitoCarta3.0 and convert to GMT")
parser$add_argument("--force", action = "store_true",
                    help = "Force re-download even if cached")
args <- parser$parse_args()

# ── Paths ─────────────────────────────────────────────────────────────────────
args_all  <- commandArgs(trailingOnly = FALSE)
file_args <- args_all[startsWith(args_all, "--file=")]
if (length(file_args) > 0) {
  script_file <- sub("^--file=", "", file_args[1])
  BASE_DIR    <- dirname(dirname(normalizePath(script_file)))
} else {
  BASE_DIR <- normalizePath(".")
}

GENE_SETS_DIR <- file.path(BASE_DIR, "processed", "gene_sets")
XLS_CACHE     <- file.path(GENE_SETS_DIR, "Mouse.MitoCarta3.0.xls")
GMT_OUTPUT    <- file.path(GENE_SETS_DIR, "mitocarta3_mouse.gmt")
MITOCARTA_URL <- "https://personal.broadinstitute.org/scalvo/MitoCarta3.0/Mouse.MitoCarta3.0.xls"

dir.create(GENE_SETS_DIR, recursive = TRUE, showWarnings = FALSE)

# ── Download ──────────────────────────────────────────────────────────────────
cat("=== MitoCarta3.0 GMT Preparation ===\n\n")

if (!file.exists(XLS_CACHE) || args$force) {
  cat(sprintf("Downloading MitoCarta3.0 from:\n  %s\n", MITOCARTA_URL))
  download.file(MITOCARTA_URL, destfile = XLS_CACHE, mode = "wb", quiet = FALSE)
  cat(sprintf("  Saved: %s (%.1f MB)\n\n", basename(XLS_CACHE),
              file.size(XLS_CACHE) / 1e6))
} else {
  cat(sprintf("Using cached: %s (%.1f MB)\n\n", XLS_CACHE,
              file.size(XLS_CACHE) / 1e6))
}

# ── Read Excel ────────────────────────────────────────────────────────────────
if (!requireNamespace("readxl", quietly = TRUE)) {
  stop("Package 'readxl' required. Install with: install.packages('readxl')")
}

cat("Reading Sheet 'A Mouse MitoCarta3.0'...\n")
mc <- readxl::read_xls(XLS_CACHE, sheet = "A Mouse MitoCarta3.0")
cat(sprintf("  %d genes × %d columns\n", nrow(mc), ncol(mc)))

# ── Identify columns ─────────────────────────────────────────────────────────
# Find Symbol column
symbol_col <- grep("^Symbol$", colnames(mc), value = TRUE)
if (length(symbol_col) == 0) {
  # Fallback: look for case-insensitive match
  symbol_col <- grep("symbol", colnames(mc), ignore.case = TRUE, value = TRUE)[1]
}
if (is.na(symbol_col) || length(symbol_col) == 0) {
  stop("Cannot find 'Symbol' column in MitoCarta sheet")
}

# Find MitoPathways column
pathway_col <- grep("MitoPathway", colnames(mc), value = TRUE)
if (length(pathway_col) == 0) {
  pathway_col <- grep("mito.*pathway", colnames(mc), ignore.case = TRUE, value = TRUE)
}
if (length(pathway_col) == 0) {
  stop("Cannot find MitoPathways column. Available columns:\n  ",
       paste(head(colnames(mc), 20), collapse = "\n  "))
}
# Use the last match (most specific, e.g., MitoCarta3.0_MitoPathways)
pathway_col <- pathway_col[length(pathway_col)]

cat(sprintf("  Symbol column: '%s'\n", symbol_col))
cat(sprintf("  Pathway column: '%s'\n", pathway_col))

# ── Parse pathways ────────────────────────────────────────────────────────────
cat("\nParsing MitoPathways...\n")

symbols  <- mc[[symbol_col]]
pathways <- mc[[pathway_col]]

# Build pathway → gene list
gs <- list()
n_assigned <- 0

for (i in seq_along(symbols)) {
  sym <- symbols[i]
  pw  <- pathways[i]

  # Skip if no symbol or no pathway
  if (is.na(sym) || sym == "" || is.na(pw) || pw == "" || pw == "0") next

  # Split by "|" (MitoCarta uses pipe separator)
  pws <- trimws(strsplit(as.character(pw), "\\|")[[1]])

  for (p in pws) {
    if (p == "" || is.na(p) || p == "0") next

    # Clean pathway name: preserve hierarchy with underscores
    # e.g., "Metabolism > TCA cycle" → "MITOCARTA_METABOLISM_TCA_CYCLE"
    p_clean <- toupper(p)
    p_clean <- gsub("\\s*>\\s*", "_", p_clean)    # hierarchy separator
    p_clean <- gsub("[^A-Z0-9_]", "_", p_clean)    # special chars
    p_clean <- gsub("_+", "_", p_clean)             # collapse multiple _
    p_clean <- gsub("^_|_$", "", p_clean)           # trim leading/trailing
    p_clean <- paste0("MITOCARTA_", p_clean)

    gs[[p_clean]] <- c(gs[[p_clean]], sym)
    n_assigned <- n_assigned + 1
  }
}

# Deduplicate genes within each pathway
gs <- lapply(gs, unique)

cat(sprintf("  %d pathways extracted\n", length(gs)))
cat(sprintf("  %d unique genes across all pathways\n", length(unique(unlist(gs)))))
cat(sprintf("  %d gene-pathway assignments\n", n_assigned))

# Size distribution
sizes <- sapply(gs, length)
cat(sprintf("  Pathway sizes: min=%d, median=%d, max=%d\n",
            min(sizes), median(sizes), max(sizes)))
cat(sprintf("  Pathways with >= 15 genes: %d / %d\n",
            sum(sizes >= 15), length(gs)))

# ── Write GMT ─────────────────────────────────────────────────────────────────
cat(sprintf("\nWriting GMT: %s\n", GMT_OUTPUT))

con <- file(GMT_OUTPUT, "w")
for (pw_name in sort(names(gs))) {
  genes <- gs[[pw_name]]
  # GMT format: pathway_name<TAB>description<TAB>gene1<TAB>gene2<TAB>...
  line <- paste(c(pw_name, paste0("MitoCarta3.0: ", pw_name), genes),
                collapse = "\t")
  writeLines(line, con)
}
close(con)

cat(sprintf("  %d lines written\n", length(gs)))

# ── Summary ──────────────────────────────────────────────────────────────────
cat("\n=== Summary ===\n")
cat(sprintf("GMT file: %s\n", GMT_OUTPUT))
cat(sprintf("Pathways: %d total (%d with >=15 genes)\n",
            length(gs), sum(sizes >= 15)))
cat(sprintf("Genes: %d unique\n", length(unique(unlist(gs)))))

# Print top 20 pathways by size
cat("\nTop 20 pathways by gene count:\n")
top_idx <- head(order(-sizes), 20)
for (i in top_idx) {
  cat(sprintf("  %-60s %d genes\n", names(gs)[i], sizes[i]))
}

cat("\nDone.\n")
