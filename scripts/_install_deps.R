#!/usr/bin/env Rscript
# Install fGSEA dependencies

has_biocm <- requireNamespace("BiocManager", quietly = TRUE)
if (!has_biocm) {
  install.packages("BiocManager", repos = "https://cran.r-project.org")
}

has_fgsea <- requireNamespace("fgsea", quietly = TRUE)
if (!has_fgsea) {
  BiocManager::install("fgsea", ask = FALSE, update = FALSE)
}

has_msigdbr <- requireNamespace("msigdbr", quietly = TRUE)
if (!has_msigdbr) {
  install.packages("msigdbr", repos = "https://cran.r-project.org")
}

has_gsva <- requireNamespace("GSVA", quietly = TRUE)
if (!has_gsva) {
  BiocManager::install("GSVA", ask = FALSE, update = FALSE)
}

has_sva <- requireNamespace("sva", quietly = TRUE)
if (!has_sva) {
  BiocManager::install("sva", ask = FALSE, update = FALSE)
}

has_limma <- requireNamespace("limma", quietly = TRUE)
if (!has_limma) {
  BiocManager::install("limma", ask = FALSE, update = FALSE)
}

has_argparse <- requireNamespace("argparse", quietly = TRUE)
if (!has_argparse) {
  install.packages("argparse", repos = "https://cran.r-project.org")
}

cat("Package check:\n")
cat("  fgsea:    ", requireNamespace("fgsea", quietly = TRUE), "\n")
cat("  msigdbr:  ", requireNamespace("msigdbr", quietly = TRUE), "\n")
cat("  GSVA:     ", requireNamespace("GSVA", quietly = TRUE), "\n")
cat("  sva:      ", requireNamespace("sva", quietly = TRUE), "\n")
cat("  limma:    ", requireNamespace("limma", quietly = TRUE), "\n")
cat("  argparse: ", requireNamespace("argparse", quietly = TRUE), "\n")
cat("Done.\n")
