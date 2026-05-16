#!/usr/bin/env Rscript
# F1: Cell-type-specific fGSEA on I4 PBMC snRNA-seq (GLDS-562)
# Data: SpaceOmicsBench v2_public/data/transcriptomics/single_cell/
#       GLDS-562_snRNA-Seq_PBMC_Gene_Expression_snRNA-seq_Processed_Data.xlsx
# Sheet: I4-FP1 → (R+1) vs (L-92, L-44, L-3) per cell type
# Output: v2/processed/F1_scrna/i4_snrnaseq_celltype_fgsea.csv

suppressPackageStartupMessages({
  library(fgsea)
  library(readxl)
  library(dplyr)
  library(tibble)
})

set.seed(42)

# ── Paths ──────────────────────────────────────────────────────────────────
SOB_DIR <- file.path(
  Sys.getenv("SPACEOMICS_ROOT", file.path(Sys.getenv("GENELAB_ROOT", "/path/to/GeneLab_benchmark"), "..", "SpaceOmicsBench")),
  "v2_public/data/transcriptomics/single_cell")
EXCEL_PATH <- file.path(SOB_DIR,
  "GLDS-562_snRNA-Seq_PBMC_Gene_Expression_snRNA-seq_Processed_Data.xlsx")
GMT_PATH   <- "/tmp/h.all.human.symbols.gmt"
OUT_DIR    <- "v2/processed/F1_scrna"
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

# ── Download Hallmark GMT if needed ────────────────────────────────────────
if (!file.exists(GMT_PATH)) {
  message("Downloading Hallmark GMT...")
  download.file(
    "https://data.broadinstitute.org/gsea-msigdb/msigdb/release/2023.2.Hs/h.all.v2023.2.Hs.symbols.gmt",
    GMT_PATH, quiet = TRUE
  )
}
pathways <- gmtPathways(GMT_PATH)
message(sprintf("Loaded %d Hallmark pathways", length(pathways)))

# ── Load snRNA-seq data (I4-FP1: R+1 vs pre-flight) ───────────────────────
message("Reading GLDS-562 snRNA-seq I4-FP1...")
df_raw <- read_excel(EXCEL_PATH, sheet = "I4-FP1", skip = 7, col_names = TRUE)
colnames(df_raw) <- c("Gene", "p_val", "avg_log2FC", "pct1", "pct2", "p_val_adj", "Cell_Type")
df_raw <- df_raw %>% filter(!is.na(Gene), !is.na(avg_log2FC))

cell_types <- unique(df_raw$Cell_Type)
message(sprintf("Cell types (%d): %s", length(cell_types), paste(cell_types, collapse = ", ")))

# ── Run fGSEA per cell type ────────────────────────────────────────────────
run_fgsea_celltype <- function(ct, df, pathways) {
  sub <- df %>% filter(Cell_Type == ct)

  # Deduplicate: keep max |avg_log2FC| per gene
  sub <- sub %>%
    group_by(Gene) %>%
    slice_max(abs(avg_log2FC), n = 1, with_ties = FALSE) %>%
    ungroup()

  ranks <- setNames(sub$avg_log2FC, sub$Gene)
  ranks <- sort(ranks, decreasing = TRUE)

  message(sprintf("  %s: %d genes", ct, length(ranks)))

  res <- fgsea(
    pathways   = pathways,
    stats      = ranks,
    minSize    = 10,
    maxSize    = 500,
    nPermSimple = 10000,
    eps        = 0
  )

  res %>%
    as_tibble() %>%
    mutate(
      Cell_Type  = ct,
      mission    = "I4",
      comparison = "FP1_R1_vs_preflight",
      leadingEdge_str = sapply(leadingEdge, paste, collapse = ";")
    ) %>%
    select(-leadingEdge)
}

message("\nRunning fGSEA per cell type...")
results_list <- lapply(cell_types, function(ct) {
  tryCatch(
    run_fgsea_celltype(ct, df_raw, pathways),
    error = function(e) {
      message(sprintf("  ERROR in %s: %s", ct, e$message))
      NULL
    }
  )
})

results_all <- bind_rows(results_list)

# ── Save ───────────────────────────────────────────────────────────────────
out_path <- file.path(OUT_DIR, "i4_snrnaseq_celltype_fgsea.csv")
write.csv(results_all, out_path, row.names = FALSE)

# ── Summary ────────────────────────────────────────────────────────────────
message(sprintf("\nDone! %d rows saved to %s", nrow(results_all), out_path))
message("\n=== Significant pathways (padj<0.05) per cell type ===")
results_all %>%
  filter(padj < 0.05) %>%
  count(Cell_Type, name = "n_sig") %>%
  arrange(desc(n_sig)) %>%
  print(n = 20)

message("\n=== Top NES pathways across all cell types (padj<0.05) ===")
results_all %>%
  filter(padj < 0.05) %>%
  arrange(desc(abs(NES))) %>%
  select(pathway, Cell_Type, NES, padj) %>%
  head(20) %>%
  print(n = 20)
