#!/usr/bin/env Rscript
# ==============================================================================
# CLI Utility: Reproduce / Compute Biological Annotation Scoring Matrix
# Author: FDRreg Research Team
# Reference: Supplementary Methods and Materials, Section 1
# ==============================================================================

suppressPackageStartupMessages({
  library(fdrregGenomics)
  library(utils)
})

args <- commandArgs(trailingOnly = TRUE)

if (length(args) == 0 || any(args %in% c("-h", "--help"))) {
  cat("
Usage: Rscript score_david_annotations.r --input=<raw_annotation.csv> --output=<scored_output.csv> [options]

Options:
  --input=<file>       Path to input CSV/TSV file containing raw annotation query text.
  --output=<file>      Path to output CSV file for scored numeric features.
  --id_col=<name>      Gene ID column name (default: 'ID').
  --sep=<sep>          Delimiter for input file (default: ',' or auto-detected).
  -h, --help           Show this help message.

Example:
  Rscript score_david_annotations.r --input=raw_david.csv --output=scored_bio_matrix.csv --id_col=ID
\n")
  quit(save = "no", status = 0)
}

parse_arg <- function(arg_name, default_val = NULL) {
  match_arg <- grep(paste0("^--", arg_name, "="), args, value = TRUE)
  if (length(match_arg) > 0) {
    return(sub(paste0("^--", arg_name, "="), "", match_arg[length(match_arg)]))
  }
  default_val
}

input_file  <- parse_arg("input")
output_file <- parse_arg("output")
id_col      <- parse_arg("id_col", default_val = "ID")

if (is.null(input_file) || !file.exists(input_file)) {
  stop("Input file is required and must exist. Use --input=<file>", call. = FALSE)
}
if (is.null(output_file)) {
  stop("Output file is required. Use --output=<file>", call. = FALSE)
}

cat(sprintf("[score_david_annotations] Reading input: %s\n", input_file))
raw_df <- utils::read.csv(input_file, stringsAsFactors = FALSE, check.names = FALSE)

cat(sprintf("[score_david_annotations] Scoring %d genes on ID column '%s'...\n", nrow(raw_df), id_col))
scored_df <- score_david_annotations(raw_df, id_col = id_col)

cat(sprintf("[score_david_annotations] Saving scored matrix (%d features) to %s\n",
            ncol(scored_df) - 1, output_file))
utils::write.csv(scored_df, file = output_file, row.names = FALSE, quote = FALSE)
cat("[score_david_annotations] Completed successfully.\n")
