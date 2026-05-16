#!/usr/bin/env Rscript
# Convert RRRM-2 Seurat objects → MTX + metadata + embeddings
# Usage: Rscript convert_rrrm2_seurat_to_h5ad.R <GLDS_NUMBER>
#
# Outputs: GLDS_{N}_counts.mtx, _genes.txt, _barcodes.txt, _metadata.csv,
#          _umap.csv, _pca.csv, _harmony.csv (if available)
#
# Requires: seurat.v5.R.4.3.3 conda env on Cayuga

suppressPackageStartupMessages({
  library(Seurat)
  library(Matrix)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 1) {
  stop("Usage: Rscript convert_rrrm2_seurat_to_h5ad.R <GLDS_NUMBER>")
}
glds <- args[1]

BASE_DIR <- Sys.getenv("RRRM2_SEURAT_DIR",
  "/athena/masonlab/scratch/projects/GeneLab/ISS_single_cell/ISS_mouse_single_cell/result/GEX/2022_08_06_Seurat_object")
OUT_DIR <- file.path(Sys.getenv("SCRATCH_DIR", Sys.getenv("HOME")),
  "huggingface/benchmark/GeneLab_benchmark/v3/data/rrrm2")
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

glds_dir <- file.path(BASE_DIR, paste0("GLDS-", glds))

cat("=== GLDS-", glds, " ===\n", sep = "")

# --- Step 1: Load Seurat object ---
# Priority: Cell_annotation labeled file (has predicted.id) > UMAP labeled file
ca_path <- file.path(glds_dir, "Cell_annotation", "temp.droplet.labeled.sobj")
ca_path_alt <- file.path(glds_dir, "Cell_annotation", "2022_10_26",
                          "temp.droplet.labeled.sobj")
umap_path <- file.path(glds_dir, "UMAP",
                         paste0("GLDS.", glds, ".harmony.integrated.labeled.rds"))

if (file.exists(ca_path)) {
  cat("Loading Cell_annotation file...\n")
  sobj <- readRDS(ca_path)
  cat("Source: Cell_annotation/temp.droplet.labeled.sobj\n")
} else if (file.exists(ca_path_alt)) {
  cat("Loading Cell_annotation (2022_10_26) file...\n")
  sobj <- readRDS(ca_path_alt)
  cat("Source: Cell_annotation/2022_10_26/temp.droplet.labeled.sobj\n")
} else {
  cat("No Cell_annotation file. Loading UMAP labeled file...\n")
  sobj <- readRDS(umap_path)
  cat("Source: UMAP/labeled.rds\n")
}

cat("Cells:", ncol(sobj), "\n")
cat("Genes:", nrow(sobj), "\n")
cat("Assays:", paste(Assays(sobj), collapse = ", "), "\n")

# --- Step 2: Metadata summary ---
cat("\nMetadata columns:\n")
for (col in colnames(sobj@meta.data)) {
  vals <- unique(sobj@meta.data[[col]])
  if (length(vals) <= 30) {
    cat("  ", col, " (n=", length(vals), "): ",
        paste(head(vals, 10), collapse = ", "), "\n", sep = "")
  } else {
    cat("  ", col, " (n=", length(vals), "): [many]\n", sep = "")
  }
}

if ("exp" %in% colnames(sobj@meta.data)) {
  cat("\nCondition (exp) counts:\n")
  print(table(sobj@meta.data$exp))
}
if ("predicted.id" %in% colnames(sobj@meta.data)) {
  cat("\nCell type (predicted.id) counts:\n")
  print(sort(table(sobj@meta.data$predicted.id), decreasing = TRUE))
}

# --- Step 3: Extract data ---
DefaultAssay(sobj) <- "RNA"

cat("\nExtracting counts...\n")
counts <- GetAssayData(sobj, assay = "RNA", layer = "counts")
cat("Counts matrix:", nrow(counts), "x", ncol(counts), "\n")

# --- Step 4: Save outputs ---
out_prefix <- file.path(OUT_DIR, paste0("GLDS_", glds))

# Sparse counts (genes × cells)
cat("\nSaving counts matrix (MTX)...\n")
writeMM(counts, paste0(out_prefix, "_counts.mtx"))

# Gene names and cell barcodes
writeLines(rownames(counts), paste0(out_prefix, "_genes.txt"))
writeLines(colnames(counts), paste0(out_prefix, "_barcodes.txt"))

# Full metadata CSV
cat("Saving metadata...\n")
write.csv(sobj@meta.data, paste0(out_prefix, "_metadata.csv"), row.names = TRUE)

# Embeddings
for (red in Reductions(sobj)) {
  tryCatch({
    emb <- Embeddings(sobj, reduction = red)
    emb_path <- paste0(out_prefix, "_", red, ".csv")
    write.csv(emb, emb_path, row.names = TRUE)
    cat("Saved embedding:", red, "- dims:", ncol(emb), "\n")
  }, error = function(e) {
    cat("Skipping reduction:", red, "-", e$message, "\n")
  })
}

# Report
cat("\nOutput files:\n")
for (f in list.files(OUT_DIR, pattern = paste0("GLDS_", glds))) {
  fpath <- file.path(OUT_DIR, f)
  cat("  ", f, " (", round(file.info(fpath)$size / 1e6, 1), " MB)\n", sep = "")
}

cat("\n=== GLDS-", glds, " extraction complete ===\n", sep = "")
rm(sobj, counts)
gc()
