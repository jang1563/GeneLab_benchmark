#!/usr/bin/env Rscript
# Quick inspection of RRRM-2 Seurat objects
# Runs on Cayuga with seurat.v5.R.4.3.3 conda env

library(Seurat)

BASE <- "/athena/masonlab/scratch/projects/GeneLab/ISS_single_cell/ISS_mouse_single_cell/result/GEX/2022_08_06_Seurat_object"

for (glds in c("402", "403", "404", "405")) {
  rds_path <- file.path(BASE, paste0("GLDS-", glds), "UMAP",
                         paste0("GLDS.", glds, ".harmony.integrated.labeled.rds"))

  cat("\n=== GLDS-", glds, " ===\n", sep = "")
  cat("Path:", rds_path, "\n")
  cat("File size:", round(file.info(rds_path)$size / 1e9, 2), "GB\n")

  sobj <- readRDS(rds_path)

  cat("Cells:", ncol(sobj), "\n")
  cat("Genes:", nrow(sobj), "\n")
  cat("Assays:", paste(Assays(sobj), collapse = ", "), "\n")
  cat("Reductions:", paste(Reductions(sobj), collapse = ", "), "\n")

  cat("\nMetadata columns:\n")
  for (col in colnames(sobj@meta.data)) {
    vals <- unique(sobj@meta.data[[col]])
    if (length(vals) <= 20) {
      cat("  ", col, " (n=", length(vals), "): ",
          paste(head(vals, 10), collapse = ", "), "\n", sep = "")
    } else {
      cat("  ", col, " (n=", length(vals), "): [too many unique values]\n", sep = "")
    }
  }

  cat("\nCondition counts:\n")
  if ("condition" %in% colnames(sobj@meta.data)) {
    print(table(sobj@meta.data$condition))
  } else if ("orig.ident" %in% colnames(sobj@meta.data)) {
    print(table(sobj@meta.data$orig.ident))
  }

  cat("\n")
  rm(sobj)
  gc()
}
