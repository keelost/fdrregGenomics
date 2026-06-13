# ===========================================================================
# R/fdrreg_core.R — Core FDRreg fitting, assessment, output assembly
# ===========================================================================

#' Fit FDRreg Models
#' @param target_z Numeric vector.
#' @param features Numeric matrix.
#' @param nulltype Character.
#' @param method Character.
#' @return Named list of FDRreg objects.
#' @noRd
fit_fdrreg_models <- function(target_z, features,
                              nulltype = c("both", "theoretical", "empirical"),
                              method = "pr") {
  nulltype <- match.arg(nulltype)
  models <- list()
  if (nulltype %in% c("both", "theoretical")) {
    message("[FDRreg] Fitting with nulltype = 'theoretical' ...")
    models$theoretical <- tryCatch(
      FDRreg(target_z, features, nulltype = "theoretical", method = method),
      error = function(e) {
        warning("FDRreg (theoretical) failed: ", e$message, call. = FALSE)
        NULL
      }
    )
  }
  if (nulltype %in% c("both", "empirical")) {
    message("[FDRreg] Fitting with nulltype = 'empirical' ...")
    models$empirical <- tryCatch(
      FDRreg(target_z, features, nulltype = "empirical", method = method),
      error = function(e) {
        warning("FDRreg (empirical) failed: ", e$message, call. = FALSE)
        NULL
      }
    )
  }
  models
}

#' Extract Feature Assessment
#' @param fdrreg_model Fitted model.
#' @param feature_names Character vector.
#' @return Data frame or NULL.
#' @noRd
extract_model_assessment <- function(fdrreg_model, feature_names = NULL) {
  if (is.null(fdrreg_model) || is.null(fdrreg_model$model)) return(NULL)
  tryCatch({
    coefs <- fdrreg_model$model$coef
    hessian <- fdrreg_model$model$hessian
    if (is.null(coefs) || is.null(hessian)) return(NULL)
    ses <- SEfromHessian(hessian)
    beta <- coefs[-1]
    se <- ses[-1]
    z_score <- beta / se
    p_value <- 2 * pnorm(abs(z_score), lower.tail = FALSE)
    if (is.null(feature_names)) feature_names <- paste0("feature_", seq_along(beta))
    data.frame(feature = feature_names, beta = beta, se = se,
               z_score = z_score, p_value = p_value,
               stringsAsFactors = FALSE)
  }, error = function(e) {
    warning("Feature assessment failed: ", e$message, call. = FALSE)
    NULL
  })
}

#' Assemble Full Results Table
#' @param ids Character vector.
#' @param z_scores Numeric vector.
#' @param p_values Numeric vector.
#' @param models Named list.
#' @return Data frame.
#' @noRd
assemble_full_results <- function(ids, z_scores, p_values, models) {
  result <- data.frame(id = ids, z = z_scores, p = p_values,
                       stringsAsFactors = FALSE)
  if (!is.null(models$theoretical)) {
    result$fdr_theo  <- models$theoretical$FDR
    result$lfdr_theo <- models$theoretical$fdr
    result$pep_theo  <- models$theoretical$pep
  }
  if (!is.null(models$empirical)) {
    result$fdr_emp  <- models$empirical$FDR
    result$lfdr_emp <- models$empirical$fdr
    result$pep_emp  <- models$empirical$pep
  }
  result
}
