# ===========================================================================
# R/features.R — Feature matrix construction and transformation
# ===========================================================================

#' Split Z-Score Matrix into Positive and Negative Components
#'
#' For each column, creates two columns: positive part and negative part.
#' This doubles the number of features, allowing FDRreg to learn
#' asymmetric effects.
#'
#' @param mat Numeric matrix with column names.
#' @return Numeric matrix with twice the number of columns.
#' @noRd
split_z_positive_negative <- function(mat) {
  cnames <- colnames(mat)
  if (is.null(cnames)) cnames <- paste0("V", seq_len(ncol(mat)))
  pos <- pmax(mat, 0)
  neg <- pmax(-mat, 0)
  result <- matrix(0, nrow = nrow(mat), ncol = ncol(mat) * 2)
  for (j in seq_len(ncol(mat))) {
    result[, 2 * j - 1] <- pos[, j]
    result[, 2 * j]     <- neg[, j]
  }
  colnames(result) <- as.vector(rbind(
    paste0(cnames, "_pos"), paste0(cnames, "_neg")
  ))
  result
}

#' Transform Feature Matrix
#' @param feat_matrix Numeric matrix.
#' @param transform Character, one of "signed", "abs", "split".
#' @return Transformed matrix.
#' @noRd
transform_features <- function(feat_matrix, transform = c("signed", "abs", "split")) {
  transform <- match.arg(transform)
  switch(transform,
    signed = feat_matrix,
    abs    = { out <- abs(feat_matrix); colnames(out) <- colnames(feat_matrix); out },
    split  = split_z_positive_negative(feat_matrix)
  )
}

#' Build Combined Feature Matrix
#' @param aux_features Numeric matrix or NULL.
#' @param annotations Numeric matrix or NULL.
#' @return Numeric matrix.
#' @noRd
build_combined_features <- function(aux_features = NULL, annotations = NULL) {
  parts <- list()
  if (!is.null(aux_features) && ncol(aux_features) > 0) parts[["aux"]] <- aux_features
  if (!is.null(annotations) && ncol(annotations) > 0) parts[["annot"]] <- annotations
  if (length(parts) == 0) {
    stop("No features provided. Supply 'aux' or 'annotations'.", call. = FALSE)
  }
  result <- do.call(cbind, parts)
  result[!is.finite(result)] <- 0
  result
}
