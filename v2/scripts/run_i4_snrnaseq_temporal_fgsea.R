#!/usr/bin/env Rscript
# F1-temporal: fGSEA across all I4 snRNA-seq timepoint comparisons
# Sheets: FP1 (R+1 vs pre), LP1 (R+45 vs pre), LP3 (R+45+R+82 vs pre),
#         RP1 (R+45 vs R+1), RP3 (R+45+R+82 vs R+1)
# Output: v2/processed/F1_scrna/i4_snrnaseq_temporal_fgsea.csv

suppressPackageStartupMessages({
  library(fgsea); library(readxl); library(dplyr); library(tibble)
})
set.seed(42)

EXCEL_PATH <- file.path(Sys.getenv("HOME"),
  "Dropbox/Bioinformatics/Claude/SpaceOmicsBench/v2_public/data/transcriptomics/single_cell",
  "GLDS-562_snRNA-Seq_PBMC_Gene_Expression_snRNA-seq_Processed_Data.xlsx")
GMT_PATH   <- "/tmp/h.all.human.symbols.gmt"
OUT_DIR    <- "v2/processed/F1_scrna"
dir.create(OUT_DIR, recursive=TRUE, showWarnings=FALSE)

if (!file.exists(GMT_PATH)) {
  download.file(
    "https://data.broadinstitute.org/gsea-msigdb/msigdb/release/2023.2.Hs/h.all.v2023.2.Hs.symbols.gmt",
    GMT_PATH, quiet=TRUE)
}
pathways <- gmtPathways(GMT_PATH)

# Sheets to process: I4-FP1 (already done), I4-LP1, I4-LP3, I4-RP1, I4-RP3
SHEETS <- list(
  "I4-FP1" = "(R+1) vs pre-flight",
  "I4-LP1" = "(R+45) vs pre-flight",
  "I4-LP3" = "(R+45, R+82) vs pre-flight",
  "I4-RP1" = "(R+45) vs (R+1)  [recovery]",
  "I4-RP3" = "(R+45, R+82) vs (R+1)  [recovery]"
)

run_sheet <- function(sheet_id) {
  message(sprintf("\n=== Sheet: %s ===", sheet_id))
  df <- read_excel(EXCEL_PATH, sheet=sheet_id, skip=7, col_names=TRUE)
  colnames(df) <- c("Gene","p_val","avg_log2FC","pct1","pct2","p_val_adj","Cell_Type")
  df <- filter(df, !is.na(Gene), !is.na(avg_log2FC))

  bind_rows(lapply(unique(df$Cell_Type), function(ct) {
    sub <- df %>% filter(Cell_Type==ct) %>%
      group_by(Gene) %>%
      slice_max(abs(avg_log2FC), n=1, with_ties=FALSE) %>% ungroup()
    ranks <- setNames(sub$avg_log2FC, sub$Gene)
    ranks <- sort(ranks, decreasing=TRUE)
    message(sprintf("  %s: %d genes", ct, length(ranks)))
    res <- tryCatch(
      fgsea(pathways=pathways, stats=ranks,
            minSize=10, maxSize=500, nPermSimple=10000, eps=0),
      error=function(e){ message("  ERROR:", e$message); NULL })
    if (is.null(res)) return(NULL)
    res %>% as_tibble() %>%
      mutate(Cell_Type=ct, sheet=sheet_id,
             leadingEdge_str=sapply(leadingEdge, paste, collapse=";")) %>%
      select(-leadingEdge)
  }))
}

# FP1 already computed; run the others
other_sheets <- setdiff(names(SHEETS), "I4-FP1")
results_list <- lapply(other_sheets, run_sheet)
results_new  <- bind_rows(results_list)

# Load FP1 and combine
fp1 <- read.csv(file.path(OUT_DIR, "i4_snrnaseq_celltype_fgsea.csv")) %>%
  mutate(sheet="I4-FP1")
results_all <- bind_rows(fp1, results_new)

out_path <- file.path(OUT_DIR, "i4_snrnaseq_temporal_fgsea.csv")
write.csv(results_all, out_path, row.names=FALSE)

message(sprintf("\nDone! %d rows saved to %s", nrow(results_all), out_path))
message("\n=== Sig pathways (padj<0.05) per sheet × cell type ===")
results_all %>%
  filter(padj < 0.05) %>%
  count(sheet, Cell_Type, name="n_sig") %>%
  arrange(sheet, desc(n_sig)) %>%
  print(n=40)
