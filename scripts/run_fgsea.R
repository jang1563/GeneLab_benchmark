#!/usr/bin/env Rscript
# run_fgsea.R — GeneLab_benchmark: Group-Level Pathway Enrichment (DD-15)
#
# Runs fGSEA on GeneLab DESeq2 DGE output (Wald statistic as ranking metric).
# OSDR does NOT provide fGSEA results — this script fills that gap.
#
# Gene set DBs: MSigDB Hallmark (primary), KEGG, Reactome (secondary)
# Species: Mus musculus (via msigdbr)
#
# Output: processed/fgsea/{tissue}/{mission}_fgsea_{db}.csv
#
# Usage:
#   Rscript scripts/run_fgsea.R --tissue liver --mission RR-1
#   Rscript scripts/run_fgsea.R --tissue liver --all
#   Rscript scripts/run_fgsea.R --all-tissues
#   Rscript scripts/run_fgsea.R --tissue liver --all --db hallmark,kegg

suppressPackageStartupMessages({
  library(fgsea)
  library(msigdbr)
  library(data.table)
  library(argparse)
})

# ── CLI ─────────────────────────────────────────────────────────────────────────
parser <- ArgumentParser(description = "Run fGSEA on GeneLab DGE files (DD-15)")
parser$add_argument("--tissue", default = NULL,
                    help = "Tissue to process (liver, kidney, thymus, etc.)")
parser$add_argument("--mission", default = NULL,
                    help = "Specific mission (e.g., RR-1, RR-3)")
parser$add_argument("--all", action = "store_true",
                    help = "Process all missions for the given tissue")
parser$add_argument("--all-tissues", action = "store_true",
                    help = "Process all tissues and all missions")
parser$add_argument("--db", default = "hallmark,kegg,reactome",
                    help = "Comma-separated gene set DBs (hallmark,kegg,reactome)")
parser$add_argument("--contrast", default = NULL,
                    help = "Override contrast column name (regex pattern for Stat_ column)")
parser$add_argument("--min-size", type = "integer", default = 15L,
                    help = "Minimum gene set size (default: 15)")
parser$add_argument("--max-size", type = "integer", default = 500L,
                    help = "Maximum gene set size (default: 500)")
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

DATA_DIR      <- file.path(BASE_DIR, "data", "mouse")
PROCESSED_DIR <- file.path(BASE_DIR, "processed")
FGSEA_DIR     <- file.path(PROCESSED_DIR, "fgsea")
SYMBOL_MAP    <- file.path(PROCESSED_DIR, "ensembl_symbol_map.csv")

# ── Tissue-Mission-OSD mapping ──────────────────────────────────────────────────
# Maps tissue → list of (mission_name, mission_dir, glds_prefix, osd_id)
TISSUE_MISSIONS <- list(
  liver = list(
    list(mission = "RR-1",  dir = "RR-1",  glds = "GLDS-48",  osd = "OSD-48"),
    list(mission = "RR-3",  dir = "RR-3",  glds = "GLDS-137", osd = "OSD-137"),
    list(mission = "RR-6",  dir = "RR-6",  glds = "GLDS-245", osd = "OSD-245"),
    list(mission = "RR-8",  dir = "RR-8",  glds = "GLDS-379", osd = "OSD-379"),
    list(mission = "RR-9",  dir = "RR-9",  glds = "GLDS-242", osd = "OSD-242"),
    list(mission = "MHU-2", dir = "MHU-2", glds = "GLDS-617", osd = "OSD-686")
  ),
  kidney = list(
    list(mission = "RR-1",  dir = "RR-1",  glds = "GLDS-102", osd = "OSD-102"),
    list(mission = "RR-3",  dir = "RR-3",  glds = "GLDS-163", osd = "OSD-163"),
    list(mission = "RR-7",  dir = "RR-7",  glds = "GLDS-253", osd = "OSD-253")
  ),
  thymus = list(
    list(mission = "RR-6",  dir = "RR-6",  glds = "GLDS-244", osd = "OSD-244"),
    list(mission = "MHU-2", dir = "MHU-2", glds = "GLDS-289", osd = "OSD-289"),
    list(mission = "RR-9",  dir = "RR-9",  glds = "GLDS-421", osd = "OSD-421")
  ),
  gastrocnemius = list(
    list(mission = "RR-1",  dir = "RR-1",  glds = "GLDS-101", osd = "OSD-101"),
    list(mission = "RR-5",  dir = "RR-5",  glds = "GLDS-401", osd = "OSD-401"),
    list(mission = "RR-9",  dir = "RR-9",  glds = "GLDS-326", osd = "OSD-326")
  ),
  eye = list(
    list(mission = "RR-1",  dir = "RR-1",  glds = "GLDS-100", osd = "OSD-100"),
    list(mission = "RR-3",  dir = "RR-3",  glds = "GLDS-194", osd = "OSD-194"),
    list(mission = "TBD",   dir = "TBD",   glds = "GLDS-397", osd = "OSD-397")
  ),
  skin = list(
    list(mission = "RR-6",          dir = "RR-6",            glds = "GLDS-243", osd = "OSD-243"),
    list(mission = "MHU-2_dorsal",  dir = "MHU-2_(dorsal)",  glds = "GLDS-238", osd = "OSD-238"),
    list(mission = "MHU-2_femoral", dir = "MHU-2_(femoral)", glds = "GLDS-239", osd = "OSD-239")
    # RR-7 (GLDS-254): normalized counts only, no DGE file → excluded
  )
)

# ── Gene Set Loading ────────────────────────────────────────────────────────────

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
    # Try KEGG_MEDICUS first, fall back to KEGG_LEGACY
    k <- tryCatch(
      msigdbr(species = "Mus musculus", category = "C2", subcategory = "CP:KEGG_MEDICUS"),
      error = function(e) {
        cat("    KEGG_MEDICUS not found, trying KEGG_LEGACY...\n")
        msigdbr(species = "Mus musculus", category = "C2", subcategory = "CP:KEGG_LEGACY")
      }
    )
    if (nrow(k) == 0) {
      # Final fallback: try without MEDICUS/LEGACY suffix
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

# ── Contrast Auto-Detection ─────────────────────────────────────────────────────
#
# GeneLab convention: Stat_(A)v(B) → positive stat = A > B (numerator higher).
# We normalize to CONTROL-FIRST direction: Stat_(GC)v(FLT) across all missions.
# This ensures consistent NES signs: positive NES = upregulated in GC (down in Flight).
#
# Known issues this handles:
#   RR-9:  Both (FLT_C1)v(GC_C2) and (GC_C2)v(FLT_C1) exist → pick GC-first
#   MHU-2: Three groups (GC, uG, centrifuge) → pick GC vs uG, not GC vs centrifuge

detect_flight_contrast <- function(col_names, override_pattern = NULL) {
  stat_cols <- col_names[grepl("^Stat_", col_names)]

  if (length(stat_cols) == 0) {
    # Fallback: try Log2fc_ columns
    log2fc_cols <- col_names[grepl("^Log2fc_", col_names)]
    if (length(log2fc_cols) == 0) {
      stop("No Stat_ or Log2fc_ columns found in DGE file")
    }
    warning("No Stat_ columns found. Falling back to Log2fc_ (less optimal for fGSEA)")
    stat_cols <- log2fc_cols
  }

  # User override
  if (!is.null(override_pattern)) {
    matched <- stat_cols[grepl(override_pattern, stat_cols, ignore.case = TRUE)]
    if (length(matched) >= 1) return(matched[1])
    warning(sprintf("Override pattern '%s' matched no columns. Using auto-detect.", override_pattern))
  }

  # --- Helper: prefer control-first direction among candidates ---
  # Parses Stat_(A)v(B) and prefers candidates where A is a control group.
  prefer_control_first <- function(candidates) {
    if (length(candidates) <= 1) return(candidates)
    control_first <- sapply(candidates, function(col) {
      m <- regmatches(col, regexec("Stat_\\((.+?)\\)v\\(", col))
      if (length(m[[1]]) < 2) return(FALSE)
      numerator <- m[[1]][2]
      grepl("^GC|^CC|^VIV|^BSL|Ground.?Control", numerator, ignore.case = TRUE)
    })
    if (any(control_first)) return(candidates[control_first][1])
    return(candidates[1])
  }

  # --- Helper: for multi-group experiments (MHU-2), prefer true microgravity ---
  # Excludes centrifuge (artificial gravity) contrasts when uG alternatives exist.
  filter_microgravity <- function(candidates) {
    if (length(candidates) <= 1) return(candidates)
    has_ug <- grepl("uG", candidates)
    has_centrifuge <- grepl("centrifug", candidates, ignore.case = TRUE)
    if (any(has_ug)) {
      ug_only <- candidates[has_ug & !has_centrifuge]
      if (length(ug_only) > 0) return(ug_only)
    }
    non_centrifuge <- candidates[!has_centrifuge]
    if (length(non_centrifuge) > 0) return(non_centrifuge)
    return(candidates)
  }

  # Priority 1: Space Flight vs Ground Control (full names)
  flight_gc <- stat_cols[grepl("Space.?Flight", stat_cols, ignore.case = TRUE) &
                         grepl("Ground.?Control", stat_cols, ignore.case = TRUE)]
  if (length(flight_gc) >= 1) {
    flight_gc <- filter_microgravity(flight_gc)
    return(prefer_control_first(flight_gc))
  }

  # Priority 2: Flight vs Vivarium
  flight_vc <- stat_cols[grepl("Space.?Flight", stat_cols, ignore.case = TRUE) &
                         grepl("Vivarium", stat_cols, ignore.case = TRUE)]
  if (length(flight_vc) >= 1) return(prefer_control_first(flight_vc))

  # Priority 3: FLT vs GC (abbreviated form, e.g., RR-9)
  flt_gc <- stat_cols[grepl("FLT", stat_cols) & grepl("GC", stat_cols)]
  if (length(flt_gc) >= 1) return(prefer_control_first(flt_gc))

  # Priority 4: Any column with "Flight"
  flight_any <- stat_cols[grepl("Flight|FLT", stat_cols, ignore.case = TRUE)]
  if (length(flight_any) >= 1) return(flight_any[1])

  stop(sprintf("No Flight-related contrast found among %d Stat_ columns.\nAvailable: %s",
               length(stat_cols), paste(head(stat_cols, 5), collapse = "\n  ")))
}

# ── Ranking Vector ──────────────────────────────────────────────────────────────

build_rank_vector <- function(dge_df, stat_col, symbol_map = NULL) {
  # Use SYMBOL column; supplement with ensembl_symbol_map if available
  if ("SYMBOL" %in% colnames(dge_df)) {
    symbols <- dge_df$SYMBOL
  } else if ("ENSEMBL" %in% colnames(dge_df) && !is.null(symbol_map)) {
    # Map Ensembl to Symbol
    idx <- match(dge_df$ENSEMBL, symbol_map$ENSEMBL)
    symbols <- symbol_map$SYMBOL[idx]
  } else {
    stop("DGE file has neither SYMBOL nor ENSEMBL column")
  }

  ranks <- as.numeric(dge_df[[stat_col]])
  names(ranks) <- symbols

  # Remove NA, empty, and duplicate symbols
  valid <- !is.na(ranks) & !is.na(names(ranks)) & names(ranks) != ""
  ranks <- ranks[valid]

  # Handle duplicates: keep the one with highest absolute value
  if (any(duplicated(names(ranks)))) {
    dt <- data.table(symbol = names(ranks), stat = ranks)
    dt[, abs_stat := abs(stat)]
    dt <- dt[order(-abs_stat)]
    dt <- dt[!duplicated(symbol)]
    ranks <- dt$stat
    names(ranks) <- dt$symbol
  }

  ranks <- sort(ranks, decreasing = TRUE)

  cat(sprintf("    Ranking vector: %d genes (after dedup + NA removal)\n", length(ranks)))
  return(ranks)
}

# ── fGSEA Execution ────────────────────────────────────────────────────────────

run_fgsea_for_db <- function(ranks, pathways, db_name, min_size, max_size) {
  cat(sprintf("    Running fGSEA [%s]: %d pathways, minSize=%d, maxSize=%d ...\n",
              db_name, length(pathways), min_size, max_size))

  res <- fgsea(
    pathways = pathways,
    stats    = ranks,
    minSize  = min_size,
    maxSize  = max_size,
    eps      = 0  # exact p-value (DD-15)
  )

  res$db <- db_name

  # Convert leadingEdge list to semicolon-separated string for CSV
  res$leadingEdge_str <- sapply(res$leadingEdge, paste, collapse = "; ")
  res$leadingEdge <- NULL

  # Summary
  n_sig <- sum(res$padj < 0.05, na.rm = TRUE)
  cat(sprintf("    Results: %d pathways tested, %d significant (padj < 0.05)\n",
              nrow(res), n_sig))

  return(as.data.frame(res))
}

# ── Find DGE File ───────────────────────────────────────────────────────────────

find_dge_file <- function(tissue, mission_info) {
  tissue_dir <- file.path(DATA_DIR, tissue, mission_info$dir)
  glds <- mission_info$glds

  # GeneLab naming patterns (v2 then v1)
  patterns <- c(
    sprintf("%s_rna_seq_differential_expression_GLbulkRNAseq.csv", glds),
    sprintf("%s_rna_seq_differential_expression.csv", glds),
    sprintf("%s_differential_expression_GLbulkRNAseq.csv", glds),
    sprintf("%s_differential_expression.csv", glds)
  )

  for (p in patterns) {
    f <- file.path(tissue_dir, p)
    if (file.exists(f)) return(f)
  }

  # Last resort: any file with "differential_expression" in name
  all_files <- list.files(tissue_dir, pattern = "differential_expression", full.names = TRUE)
  if (length(all_files) > 0) return(all_files[1])

  return(NULL)
}

# ── Main Processing ─────────────────────────────────────────────────────────────

process_mission <- function(tissue, mission_info, gene_sets, contrast_override, min_size, max_size) {
  mission <- mission_info$mission
  cat(sprintf("\n  [%s / %s] (%s)\n", tissue, mission, mission_info$glds))

  # Find DGE file
  dge_file <- find_dge_file(tissue, mission_info)
  if (is.null(dge_file)) {
    cat(sprintf("    [SKIP] No DGE file found for %s\n", mission_info$glds))
    return(NULL)
  }
  cat(sprintf("    DGE file: %s\n", basename(dge_file)))

  # Load DGE
  dge <- read.csv(dge_file, check.names = FALSE, stringsAsFactors = FALSE)
  cat(sprintf("    DGE: %d genes × %d columns\n", nrow(dge), ncol(dge)))

  # Detect flight contrast
  stat_col <- tryCatch(
    detect_flight_contrast(colnames(dge), override_pattern = contrast_override),
    error = function(e) {
      cat(sprintf("    [ERROR] %s\n", conditionMessage(e)))
      return(NULL)
    }
  )
  if (is.null(stat_col)) return(NULL)

  # Determine if this is Stat_ or Log2fc_
  is_stat <- grepl("^Stat_", stat_col)
  cat(sprintf("    Contrast: %s (%s)\n",
              substr(stat_col, 1, 80),
              ifelse(is_stat, "Wald statistic", "Log2FC fallback")))

  # Load symbol map for Ensembl fallback
  sym_map <- NULL
  if (file.exists(SYMBOL_MAP)) {
    sym_map <- read.csv(SYMBOL_MAP, stringsAsFactors = FALSE)
  }

  # Build ranking vector
  ranks <- tryCatch(
    build_rank_vector(dge, stat_col, symbol_map = sym_map),
    error = function(e) {
      cat(sprintf("    [ERROR] Rank vector failed: %s\n", conditionMessage(e)))
      return(NULL)
    }
  )
  if (is.null(ranks) || length(ranks) < 1000) {
    cat(sprintf("    [WARN] Only %d ranked genes — results may be unreliable\n",
                length(ranks)))
  }

  # Gene symbol matching report
  for (db_name in names(gene_sets)) {
    all_gs_genes <- unique(unlist(gene_sets[[db_name]]))
    n_matched <- sum(all_gs_genes %in% names(ranks))
    pct <- round(100 * n_matched / length(all_gs_genes), 1)
    cat(sprintf("    Gene match [%s]: %d / %d (%.1f%%)\n",
                db_name, n_matched, length(all_gs_genes), pct))
  }

  # Run fGSEA per DB
  results <- list()
  out_dir <- file.path(FGSEA_DIR, tissue)
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

  for (db_name in names(gene_sets)) {
    res <- tryCatch(
      run_fgsea_for_db(ranks, gene_sets[[db_name]], db_name, min_size, max_size),
      error = function(e) {
        cat(sprintf("    [ERROR] fGSEA %s failed: %s\n", db_name, conditionMessage(e)))
        return(NULL)
      }
    )
    if (!is.null(res)) {
      # Add metadata columns
      res$tissue  <- tissue
      res$mission <- mission
      res$glds    <- mission_info$glds

      # Save
      out_file <- file.path(out_dir, sprintf("%s_fgsea_%s.csv", mission, db_name))
      write.csv(res, out_file, row.names = FALSE)
      cat(sprintf("    Saved: %s\n", basename(out_file)))

      results[[db_name]] <- res
    }
  }

  return(results)
}

# ── Entry Point ─────────────────────────────────────────────────────────────────

cat("\n=== GeneLab_benchmark: fGSEA Pathway Enrichment (DD-15) ===\n")

# Parse DB list
db_names <- trimws(strsplit(args$db, ",")[[1]])
cat(sprintf("Gene set DBs: %s\n", paste(db_names, collapse = ", ")))

# Load gene sets (once for all missions)
cat("\nLoading gene sets...\n")
gene_sets <- load_gene_sets(db_names)

# Determine which tissues/missions to process
if (args$all_tissues) {
  tissues_to_run <- names(TISSUE_MISSIONS)
} else if (!is.null(args$tissue)) {
  tissues_to_run <- args$tissue
} else {
  stop("Specify --tissue, --all, or --all-tissues")
}

# Process each tissue
all_results <- list()

for (tissue in tissues_to_run) {
  cat(sprintf("\n\n========== %s ==========\n", toupper(tissue)))

  if (!(tissue %in% names(TISSUE_MISSIONS))) {
    cat(sprintf("[SKIP] Unknown tissue: %s\n", tissue))
    next
  }

  missions <- TISSUE_MISSIONS[[tissue]]

  # Filter to specific mission if requested
  if (!is.null(args$mission) && !args$all && !args$all_tissues) {
    missions <- Filter(function(m) m$mission == args$mission, missions)
    if (length(missions) == 0) {
      cat(sprintf("[SKIP] Mission %s not found for tissue %s\n", args$mission, tissue))
      next
    }
  }

  for (m_info in missions) {
    res <- process_mission(tissue, m_info, gene_sets, args$contrast,
                           args$min_size, args$max_size)
    if (!is.null(res)) {
      key <- sprintf("%s_%s", tissue, m_info$mission)
      all_results[[key]] <- res
    }
  }
}

# ── Summary Report ──────────────────────────────────────────────────────────────

cat("\n\n=== Summary ===\n")
cat(sprintf("Processed: %d tissue-mission combinations\n", length(all_results)))

# Create combined summary per DB
for (db_name in db_names) {
  combined <- do.call(rbind, lapply(all_results, function(r) r[[db_name]]))
  if (!is.null(combined) && nrow(combined) > 0) {
    out_file <- file.path(FGSEA_DIR, "summary", sprintf("all_fgsea_%s.csv", db_name))
    write.csv(combined, out_file, row.names = FALSE)
    cat(sprintf("Combined summary [%s]: %d rows → %s\n",
                db_name, nrow(combined), basename(out_file)))

    # Top enriched pathways across all missions (padj < 0.05)
    sig <- combined[!is.na(combined$padj) & combined$padj < 0.05, ]
    if (nrow(sig) > 0) {
      top_pw <- head(sig[order(sig$padj), c("pathway", "NES", "padj", "mission", "tissue")], 10)
      cat(sprintf("\n  Top 10 significant pathways [%s]:\n", db_name))
      for (i in 1:nrow(top_pw)) {
        cat(sprintf("    %2d. %-50s NES=%+.2f  padj=%.1e  [%s/%s]\n",
                    i, top_pw$pathway[i], top_pw$NES[i], top_pw$padj[i],
                    top_pw$tissue[i], top_pw$mission[i]))
      }
    }
  }
}

cat("\nDone.\n")
