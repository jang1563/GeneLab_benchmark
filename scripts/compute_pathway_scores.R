#!/usr/bin/env Rscript
# compute_pathway_scores.R — GeneLab_benchmark: Sample-Level Pathway Scores (DD-15)
#
# Computes per-sample pathway activity scores using GSVA.
# Used as ML features across ALL categories (A, B, C, D) for gene vs pathway comparison.
#
# Input:  log2(normalized counts) from processed/A_detection/{tissue}/
# Output: processed/pathway_scores/{tissue}/{mission}_gsva_{db}.csv
#         (samples × pathways matrix)
#
# Usage:
#   Rscript scripts/compute_pathway_scores.R --tissue liver --all
#   Rscript scripts/compute_pathway_scores.R --tissue liver --mission RR-1 --db hallmark
#   Rscript scripts/compute_pathway_scores.R --tissue liver --all --method ssgsea
#   Rscript scripts/compute_pathway_scores.R --all-tissues

suppressPackageStartupMessages({
  library(GSVA)
  library(msigdbr)
  library(argparse)
})

# ── CLI ─────────────────────────────────────────────────────────────────────────
parser <- ArgumentParser(description = "Compute per-sample pathway scores (GSVA/ssGSEA)")
parser$add_argument("--tissue", default = NULL,
                    help = "Tissue to process")
parser$add_argument("--mission", default = NULL,
                    help = "Specific mission (e.g., RR-1)")
parser$add_argument("--all", action = "store_true",
                    help = "Process all missions for the given tissue")
parser$add_argument("--all-tissues", action = "store_true",
                    help = "Process all tissues and all missions")
parser$add_argument("--db", default = "hallmark,kegg,reactome",
                    help = "Comma-separated gene set DBs")
parser$add_argument("--method", default = "gsva",
                    choices = c("gsva", "ssgsea"),
                    help = "Scoring method: gsva (default) or ssgsea")
parser$add_argument("--min-size", type = "integer", default = 15L,
                    help = "Minimum gene set size")
parser$add_argument("--max-size", type = "integer", default = 500L,
                    help = "Maximum gene set size")
args <- parser$parse_args()

# ── Paths ───────────────────────────────────────────────────────────────────────
args_all  <- commandArgs(trailingOnly = FALSE)
file_args <- args_all[startsWith(args_all, "--file=")]
if (length(file_args) > 0) {
  script_file <- sub("^--file=", "", file_args[1])
  BASE_DIR    <- dirname(dirname(normalizePath(script_file)))
} else {
  BASE_DIR <- normalizePath(".")
}

A_DETECTION_DIR  <- file.path(BASE_DIR, "processed", "A_detection")
PATHWAY_DIR      <- file.path(BASE_DIR, "processed", "pathway_scores")
SYMBOL_MAP_FILE  <- file.path(BASE_DIR, "processed", "ensembl_symbol_map.csv")

# ── Tissue-Mission mapping (same as run_fgsea.R) ───────────────────────────────
TISSUE_MISSIONS <- list(
  liver = list(
    list(mission = "RR-1",  dir = "RR-1"),
    list(mission = "RR-3",  dir = "RR-3"),
    list(mission = "RR-6",  dir = "RR-6"),
    list(mission = "RR-8",  dir = "RR-8"),
    list(mission = "RR-9",  dir = "RR-9"),
    list(mission = "MHU-2", dir = "MHU-2")
  ),
  kidney = list(
    list(mission = "RR-1",  dir = "RR-1"),
    list(mission = "RR-3",  dir = "RR-3"),
    list(mission = "RR-7",  dir = "RR-7")
  ),
  thymus = list(
    list(mission = "RR-6",  dir = "RR-6"),
    list(mission = "MHU-2", dir = "MHU-2"),
    list(mission = "RR-9",  dir = "RR-9")
  ),
  gastrocnemius = list(
    list(mission = "RR-1",  dir = "RR-1"),
    list(mission = "RR-5",  dir = "RR-5"),
    list(mission = "RR-9",  dir = "RR-9")
  ),
  eye = list(
    list(mission = "RR-1",  dir = "RR-1"),
    list(mission = "RR-3",  dir = "RR-3"),
    list(mission = "TBD",   dir = "TBD")
  )
)

# ── Gene Set Loading (shared with run_fgsea.R) ─────────────────────────────────

load_gene_sets <- function(db_names = c("hallmark", "kegg", "reactome")) {
  gs_list <- list()

  if ("hallmark" %in% db_names) {
    cat("  Loading MSigDB Hallmark (Mus musculus)...\n")
    h <- msigdbr(species = "Mus musculus", category = "H")
    gs_list[["hallmark"]] <- split(h$gene_symbol, h$gs_name)
    cat(sprintf("    %d gene sets loaded\n", length(gs_list[["hallmark"]])))
  }

  if ("kegg" %in% db_names) {
    cat("  Loading MSigDB KEGG (Mus musculus)...\n")
    k <- tryCatch(
      msigdbr(species = "Mus musculus", category = "C2", subcategory = "CP:KEGG_MEDICUS"),
      error = function(e) {
        cat("    KEGG_MEDICUS not found, trying KEGG_LEGACY...\n")
        msigdbr(species = "Mus musculus", category = "C2", subcategory = "CP:KEGG_LEGACY")
      }
    )
    if (nrow(k) == 0) {
      k <- msigdbr(species = "Mus musculus", category = "C2", subcategory = "CP:KEGG")
    }
    gs_list[["kegg"]] <- split(k$gene_symbol, k$gs_name)
    cat(sprintf("    %d gene sets loaded\n", length(gs_list[["kegg"]])))
  }

  if ("reactome" %in% db_names) {
    cat("  Loading MSigDB Reactome (Mus musculus)...\n")
    r <- msigdbr(species = "Mus musculus", category = "C2", subcategory = "CP:REACTOME")
    gs_list[["reactome"]] <- split(r$gene_symbol, r$gs_name)
    cat(sprintf("    %d gene sets loaded\n", length(gs_list[["reactome"]])))
  }

  return(gs_list)
}

# ── Find Normalized Counts File ─────────────────────────────────────────────────

find_norm_file <- function(tissue, mission_name) {
  tissue_dir <- file.path(A_DETECTION_DIR, tissue)

  # Pattern: {tissue}_{mission}_log2_norm.csv
  patterns <- c(
    sprintf("%s_%s_log2_norm.csv", tissue, mission_name),
    sprintf("%s_%s_log2_norm.csv", tissue, gsub("-", "", mission_name))
  )

  for (p in patterns) {
    f <- file.path(tissue_dir, p)
    if (file.exists(f)) return(f)
  }

  # Fuzzy match
  all_files <- list.files(tissue_dir, pattern = sprintf(".*%s.*log2_norm\\.csv",
                                                        gsub("-", ".", mission_name)),
                          full.names = TRUE, ignore.case = TRUE)
  if (length(all_files) > 0) return(all_files[1])

  return(NULL)
}

# ── Map Ensembl to Symbol ───────────────────────────────────────────────────────

map_to_symbols <- function(expr_matrix, symbol_map_file) {
  # expr_matrix: samples (rows) × genes (columns), colnames = ENSMUSG...
  gene_ids <- colnames(expr_matrix)

  # Check if already symbols (not Ensembl)
  if (!any(grepl("^ENSMUSG", gene_ids))) {
    cat("    Gene IDs appear to be symbols already\n")
    return(expr_matrix)
  }

  if (!file.exists(symbol_map_file)) {
    warning("ensembl_symbol_map.csv not found. Cannot map Ensembl to symbols.")
    return(expr_matrix)
  }

  sym_map <- read.csv(symbol_map_file, stringsAsFactors = FALSE)
  idx <- match(gene_ids, sym_map$ENSEMBL)
  new_names <- ifelse(is.na(idx), gene_ids, sym_map$SYMBOL[idx])

  # Remove unmapped (still Ensembl) and empty
  valid <- !grepl("^ENSMUSG", new_names) & new_names != "" & !is.na(new_names)
  expr_matrix <- expr_matrix[, valid, drop = FALSE]
  new_names <- new_names[valid]

  # Handle duplicates: keep column with highest variance
  if (any(duplicated(new_names))) {
    col_var <- apply(expr_matrix, 2, var, na.rm = TRUE)
    keep <- !duplicated(new_names) | col_var == ave(col_var, new_names, FUN = max)
    # Even after above, ensure unique
    expr_matrix <- expr_matrix[, !duplicated(new_names), drop = FALSE]
    new_names <- new_names[!duplicated(new_names)]
  }

  colnames(expr_matrix) <- new_names
  cat(sprintf("    Mapped to symbols: %d genes\n", ncol(expr_matrix)))
  return(expr_matrix)
}

# ── GSVA Computation ────────────────────────────────────────────────────────────

compute_scores <- function(expr_matrix, gene_sets, method, min_size, max_size) {
  # expr_matrix: samples × genes → need to transpose for GSVA (genes × samples)
  expr_t <- t(expr_matrix)

  cat(sprintf("    Input matrix: %d genes × %d samples\n", nrow(expr_t), ncol(expr_t)))

  # Support both old (1.40.x) and new (1.50+) GSVA API
  has_new_api <- exists("gsvaParam", where = asNamespace("GSVA"), inherits = FALSE)

  if (has_new_api) {
    # New API (GSVA >= 1.50)
    if (method == "gsva") {
      param <- gsvaParam(
        exprData = as.matrix(expr_t),
        geneSets = gene_sets,
        kcdf     = "Gaussian",
        minSize  = min_size,
        maxSize  = max_size
      )
    } else {
      param <- ssgseaParam(
        exprData = as.matrix(expr_t),
        geneSets = gene_sets,
        minSize  = min_size,
        maxSize  = max_size
      )
    }
    scores <- gsva(param)
  } else {
    # Old API (GSVA 1.40.x)
    scores <- gsva(
      expr      = as.matrix(expr_t),
      gset.idx.list = gene_sets,
      method    = method,
      kcdf      = "Gaussian",
      min.sz    = min_size,
      max.sz    = max_size,
      verbose   = FALSE
    )
  }

  return(scores)  # pathways × samples matrix
}

# ── Main Processing ─────────────────────────────────────────────────────────────

process_mission_gsva <- function(tissue, mission_info, gene_sets, method, min_size, max_size) {
  mission <- mission_info$mission
  cat(sprintf("\n  [%s / %s]\n", tissue, mission))

  # Find normalized counts
  norm_file <- find_norm_file(tissue, mission)
  if (is.null(norm_file)) {
    cat(sprintf("    [SKIP] No log2_norm file found for %s/%s\n", tissue, mission))
    return(NULL)
  }
  cat(sprintf("    Norm file: %s\n", basename(norm_file)))

  # Load expression matrix — CSV is genes (rows) × samples (columns)
  expr_raw <- read.csv(norm_file, row.names = 1, check.names = FALSE)

  # Detect orientation: if many more rows than columns, it's genes×samples
  if (nrow(expr_raw) > ncol(expr_raw) * 10) {
    cat(sprintf("    Raw matrix: %d genes × %d samples (transposing)\n",
                nrow(expr_raw), ncol(expr_raw)))
    expr <- as.data.frame(t(expr_raw))  # → samples × genes
  } else {
    expr <- expr_raw
  }
  cat(sprintf("    Expression: %d samples × %d genes\n", nrow(expr), ncol(expr)))

  if (nrow(expr) < 3) {
    cat("    [SKIP] Too few samples (<3)\n")
    return(NULL)
  }

  # Map Ensembl IDs to symbols if needed
  expr <- map_to_symbols(expr, SYMBOL_MAP_FILE)

  # Compute GSVA per DB
  results <- list()
  out_dir <- file.path(PATHWAY_DIR, tissue)
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

  for (db_name in names(gene_sets)) {
    cat(sprintf("    Computing %s [%s] scores...\n", toupper(method), db_name))

    scores <- tryCatch({
      compute_scores(expr, gene_sets[[db_name]], method, min_size, max_size)
    }, error = function(e) {
      cat(sprintf("    [ERROR] %s %s failed: %s\n", method, db_name, conditionMessage(e)))
      return(NULL)
    })

    if (!is.null(scores)) {
      # Transpose: pathways × samples → samples × pathways
      scores_df <- as.data.frame(t(scores))
      cat(sprintf("    Result: %d samples × %d pathways\n", nrow(scores_df), ncol(scores_df)))

      # Save
      method_str <- method  # "gsva" or "ssgsea"
      out_file <- file.path(out_dir, sprintf("%s_%s_%s.csv", mission, method_str, db_name))
      write.csv(scores_df, out_file)
      cat(sprintf("    Saved: %s\n", basename(out_file)))

      results[[db_name]] <- scores_df
    }
  }

  return(results)
}

# ── Entry Point ─────────────────────────────────────────────────────────────────

cat(sprintf("\n=== GeneLab_benchmark: %s Pathway Scores (DD-15) ===\n",
            toupper(args$method)))

# Parse DB list
db_names <- trimws(strsplit(args$db, ",")[[1]])
cat(sprintf("Gene set DBs: %s\n", paste(db_names, collapse = ", ")))
cat(sprintf("Method: %s\n", args$method))

# Load gene sets
cat("\nLoading gene sets...\n")
gene_sets <- load_gene_sets(db_names)

# Determine tissues
if (args$all_tissues) {
  tissues_to_run <- names(TISSUE_MISSIONS)
} else if (!is.null(args$tissue)) {
  tissues_to_run <- args$tissue
} else {
  stop("Specify --tissue, --all, or --all-tissues")
}

# Process
for (tissue in tissues_to_run) {
  cat(sprintf("\n\n========== %s ==========\n", toupper(tissue)))

  if (!(tissue %in% names(TISSUE_MISSIONS))) {
    cat(sprintf("[SKIP] Unknown tissue: %s\n", tissue))
    next
  }

  missions <- TISSUE_MISSIONS[[tissue]]

  if (!is.null(args$mission) && !args$all && !args$all_tissues) {
    missions <- Filter(function(m) m$mission == args$mission, missions)
  }

  for (m_info in missions) {
    process_mission_gsva(tissue, m_info, gene_sets, args$method,
                         args$min_size, args$max_size)
  }
}

cat("\n\nDone.\n")
