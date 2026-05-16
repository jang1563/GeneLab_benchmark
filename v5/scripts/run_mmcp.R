#!/usr/bin/env Rscript
# GeneLabBench v5: Run mMCP-counter for mouse immune deconvolution
# Input: genes × samples CSV (gene symbols as row names)
# Output: cell types × samples CSV

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2) {
  stop("Usage: Rscript run_mmcp.R <input_csv> <output_csv>")
}

input_csv <- args[1]
output_csv <- args[2]

cat("Loading mMCPcounter...\n")

# Try mMCPcounter first, fall back to immunedeconv
tryCatch({
  library(mMCPcounter)
  cat("Using mMCPcounter directly\n")

  expr <- read.csv(input_csv, row.names = 1, check.names = FALSE)
  cat(sprintf("Expression matrix: %d genes × %d samples\n", nrow(expr), ncol(expr)))

  result <- mMCPcounter.estimate(expr, features = "Gene.Symbol")
  cat(sprintf("Deconvolution result: %d cell types × %d samples\n", nrow(result), ncol(result)))

  write.csv(as.data.frame(result), output_csv)
  cat(sprintf("Saved to %s\n", output_csv))

}, error = function(e) {
  cat(sprintf("mMCPcounter failed: %s\n", e$message))
  cat("Trying immunedeconv fallback...\n")

  tryCatch({
    library(immunedeconv)

    expr <- read.csv(input_csv, row.names = 1, check.names = FALSE)
    cat(sprintf("Expression matrix: %d genes × %d samples\n", nrow(expr), ncol(expr)))

    result <- deconvolute_mouse(expr, method = "mmcp_counter")

    # immunedeconv returns a tibble with cell_type column + sample columns
    result_df <- as.data.frame(result)
    rownames(result_df) <- result_df$cell_type
    result_df$cell_type <- NULL

    cat(sprintf("Deconvolution result: %d cell types × %d samples\n",
                nrow(result_df), ncol(result_df)))

    write.csv(result_df, output_csv)
    cat(sprintf("Saved to %s\n", output_csv))

  }, error = function(e2) {
    stop(sprintf("Both mMCPcounter and immunedeconv failed:\n  mMCPcounter: %s\n  immunedeconv: %s",
                 e$message, e2$message))
  })
})
