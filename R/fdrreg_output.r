# ===========================================================================
# R/fdrreg_output.R — S3 class fdrreg_result with print/summary/significant
# ===========================================================================

#' Constructor for fdrreg_result S3 Class
#' @noRd
new_fdrreg_result <- function(tier, full_results, summary_counts,
                              assessment = NULL, models = list(),
                              varselect_model = NULL,
                              covariates_used = character(),
                              params = list()) {
  structure(
    list(tier = tier, full_results = full_results, summary = summary_counts,
         assessment = assessment, models = models,
         varselect_model = varselect_model,
         covariates_used = covariates_used, params = params),
    class = "fdrreg_result"
  )
}

#' Print an fdrreg_result Object
#'
#' @param x An object of class \code{fdrreg_result}.
#' @param ... Additional arguments (ignored).
#' @return Invisibly returns \code{x}.
#' @export
print.fdrreg_result <- function(x, ...) {
  cat("\n-- FDRreg Genomics Result ---------------------------------\n")
  cat(sprintf("  Tier            : %s\n", x$tier))
  cat(sprintf("  Total features  : %d\n", nrow(x$full_results)))
  cat(sprintf("  Covariates used : %s\n",
              paste(x$covariates_used, collapse = ", ")))
  cat(sprintf("  Var. selection  : %s\n",
              ifelse(is.null(x$varselect_model), "none",
                     class(x$varselect_model)[1])))
  cat("\n-- Discoveries at FDR Thresholds --------------------------\n")
  if (!is.null(x$summary) && nrow(x$summary) > 0) {
    print(x$summary, right = FALSE)
  }
  if (!is.null(x$assessment) && nrow(x$assessment) > 0) {
    cat("\n-- Top Features (Model Assessment) ------------------------\n")
    top <- x$assessment[order(x$assessment$p_value), , drop = FALSE]
    top <- head(top, 10)
    print(top, row.names = FALSE, digits = 4)
  }
  cat("\n-- Reproducibility ----------------------------------------\n")
  cat(sprintf("  Seed      : %s\n", x$params$seed))
  cat(sprintf("  Timestamp : %s\n", x$params$timestamp))
  cat(sprintf("  R version : %s\n", x$params$r_version))
  cat("------------------------------------------------------------\n\n")
  invisible(x)
}

#' Summarize an fdrreg_result Object
#'
#' @param object An object of class \code{fdrreg_result}.
#' @param ... Additional arguments (ignored).
#' @return Invisibly returns \code{object}.
#' @export
summary.fdrreg_result <- function(object, ...) {
  cat(sprintf("Tier: %s\n", object$tier))
  cat(sprintf("Total features: %d\n", nrow(object$full_results)))
  if (!is.null(object$full_results$fdr_theo)) {
    n_sig <- sum(object$full_results$fdr_theo < 0.05, na.rm = TRUE)
    cat(sprintf("Significant at FDR < 0.05 (theoretical): %d\n", n_sig))
  }
  if (!is.null(object$full_results$fdr_emp)) {
    n_sig <- sum(object$full_results$fdr_emp < 0.05, na.rm = TRUE)
    cat(sprintf("Significant at FDR < 0.05 (empirical): %d\n", n_sig))
  }
  invisible(object)
}

#' Extract Significant Findings
#'
#' Returns rows from the full results table where FDR is below
#' the specified threshold.
#'
#' @param x An object of class \code{fdrreg_result}.
#' @param threshold Numeric, FDR threshold (default 0.05).
#' @param type Character, "theoretical" or "empirical".
#' @param ... Additional arguments (ignored).
#' @return A data frame of significant findings.
#' @export
significant <- function(x, threshold = 0.05,
                        type = c("theoretical", "empirical"), ...) {
  UseMethod("significant")
}

#' @export
#' @rdname significant
significant.fdrreg_result <- function(x, threshold = 0.05,
                                      type = c("theoretical", "empirical"),
                                      ...) {
  type <- match.arg(type)
  fdr_col <- if (type == "theoretical") "fdr_theo" else "fdr_emp"
  if (!fdr_col %in% colnames(x$full_results)) {
    stop(sprintf("Column '%s' not found. Was '%s' nulltype used?", fdr_col, type),
         call. = FALSE)
  }
  idx <- which(x$full_results[[fdr_col]] < threshold)
  x$full_results[idx, , drop = FALSE]
}
