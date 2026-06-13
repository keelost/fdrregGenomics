# ===========================================================================
# R/decorrelation.R — LDSC-based covariance matrix and Matpow decorrelation
# ===========================================================================

#' Build Covariance Matrix from LDSC Intercepts
#'
#' @param target_name Character, target trait name.
#' @param overlap_names Character vector, overlapping trait names.
#' @param ldsc_intercepts Data frame with LDSC intercept values.
#' @param col_trait1 Character, column name for trait 1.
#' @param col_trait2 Character, column name for trait 2.
#' @param col_intercept Character, column name for intercept.
#' @return Named numeric matrix.
#' @noRd
build_ldsc_cov_matrix <- function(target_name, overlap_names,
                                  ldsc_intercepts,
                                  col_trait1 = "trait1",
                                  col_trait2 = "trait2",
                                  col_intercept = "intercept") {
  all_traits <- c(target_name, overlap_names)
  n <- length(all_traits)
  validate_columns(ldsc_intercepts, c(col_trait1, col_trait2, col_intercept),
                   "ldsc_intercepts")
  cov_mat <- diag(n)
  rownames(cov_mat) <- all_traits
  colnames(cov_mat) <- all_traits
  for (i in seq_len(n)) {
    for (j in seq_len(n)) {
      if (i == j) next
      trait_i <- all_traits[i]
      trait_j <- all_traits[j]
      row_match <- which(
        (ldsc_intercepts[[col_trait1]] == trait_i &
           ldsc_intercepts[[col_trait2]] == trait_j) |
          (ldsc_intercepts[[col_trait1]] == trait_j &
             ldsc_intercepts[[col_trait2]] == trait_i)
      )
      if (length(row_match) >= 1) {
        cov_mat[i, j] <- ldsc_intercepts[[col_intercept]][row_match[1]]
        cov_mat[j, i] <- cov_mat[i, j]
      } else {
        warning(sprintf("LDSC intercept not found for (%s, %s); set to 0.",
                        trait_i, trait_j), call. = FALSE)
      }
    }
  }
  cov_mat
}

#' Decorrelate Z-Scores via Matpow
#'
#' @param z_matrix Numeric matrix (n_snps x n_traits).
#' @param cov_matrix Numeric matrix (n_traits x n_traits).
#' @return Numeric matrix (n_traits x n_snps).
#' @noRd
decorrelate_z_scores <- function(z_matrix, cov_matrix) {
  if (!requireNamespace("powerplus", quietly = TRUE)) {
    stop("Package 'powerplus' is required for decorrelation.", call. = FALSE)
  }
  Matpow(cov_matrix, -0.5) %*% t(z_matrix)
}
