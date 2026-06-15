# ===========================================================================
# R/fdrreg_core.R — Core FDRreg fitting, assessment, output assembly
# ===========================================================================

#' Fit FDRreg Models with Diagnostic Messages
#'
#' Wraps FDRreg with clear diagnostic messages for common failure modes:
#' optimization failure, non-finite values, and fallback to no-covariates model.
#'
#' @param target_z Numeric vector of target z-scores.
#' @param features Numeric matrix of covariates.
#' @param nulltype Character, "both", "theoretical", or "empirical".
#' @param method Character, passed to FDRreg (default "pr").
#' @return A list with elements:
#'   \describe{
#'     \item{models}{Named list of FDRreg model objects (or NULL on failure).}
#'     \item{fit_status}{Named character vector of fit outcomes.}
#'   }
#' @noRd
fit_fdrreg_models <- function(target_z, features,
                              nulltype = c("both", "theoretical", "empirical"),
                              method = "pr") {
  nulltype <- match.arg(nulltype)
  models      <- list()
  fit_status  <- list()

  run_one <- function(nt) {
    message(sprintf("[FDRreg] Fitting with nulltype = '%s' ...", nt))

    # ---- Try to fit ----
    result <- tryCatch({
      fit <- FDRreg(target_z, features, nulltype = nt, method = method)
      list(model = fit, error_msg = NULL)
    }, error = function(e) {
      list(model = NULL, error_msg = e$message)
    })

    # ---- Diagnose ----
    if (is.null(result$model)) {
      # FDRreg threw an unrecoverable error
      message(sprintf(
        paste0(
          "\n  [FDRreg ERROR] Fitting failed for nulltype='%s'.\n",
          "    Error: %s\n",
          "    -> FDR columns will be filled with NA.\n",
          "    -> Possible causes: features contain all-NA columns, or\n",
          "       sample size is too small for the number of features.\n",
          "    -> Suggestion: reduce features (var_select='lasso') or\n",
          "       check your input data.\n"
        ),
        nt, result$error_msg
      ))
      return(list(model = NULL, status = "failed"))
    }

    fit <- result$model

    # Check if FDRreg internally fell back to no-covariates model.
    # Indicator: model$coef has only intercept (length 1) or is NULL.
    has_covariates <- !is.null(fit$model) &&
                      !is.null(fit$model$coef) &&
                      length(fit$model$coef) > 1

    if (!has_covariates) {
      message(sprintf(
        paste0(
          "\n  [FDRreg WARNING] nulltype='%s': optimization failed; ",
          "reverted to no-covariates model.\n",
          "    -> FDR values are still returned (computed without feature effects).\n",
          "    -> This usually means features have numerical issues.\n",
          "    -> Suggestions:\n",
          "       1. Check for constant/zero-variance columns in features.\n",
          "       2. Reduce the number of features (var_select='lasso').\n",
          "       3. Ensure features and target have matching sample sizes.\n"
        ),
        nt
      ))
      return(list(model = fit, status = "reverted (no covariates)"))
    }

    # Check that FDR field exists
    if (is.null(fit$FDR)) {
      warning(sprintf("FDRreg (%s): model returned NULL FDR values.", nt),
              call. = FALSE)
      return(list(model = fit, status = "warning (NULL FDR)"))
    }

    return(list(model = fit, status = "success"))
  }

  # ---- Run for each nulltype ----
  if (nulltype %in% c("both", "theoretical")) {
    res <- run_one("theoretical")
    models$theoretical     <- res$model
    fit_status$theoretical <- res$status
  }

  if (nulltype %in% c("both", "empirical")) {
    res <- run_one("empirical")
    models$empirical     <- res$model
    fit_status$empirical <- res$status
  }

  list(models = models, fit_status = fit_status)
}


#' Extract Feature Assessment from FDRreg Model
#'
#' @param fdrreg_model A fitted FDRreg model object.
#' @param feature_names Character vector of feature names.
#' @return Data frame or NULL.
#' @noRd
extract_model_assessment <- function(fdrreg_model, feature_names = NULL) {
  if (is.null(fdrreg_model)) return(NULL)

  tryCatch({
    coefs   <- fdrreg_model$model$coef
    hessian <- fdrreg_model$model$hessian
    if (is.null(coefs) || is.null(hessian)) return(NULL)
    if (length(coefs) <= 1) return(NULL)  # intercept-only model

    ses      <- SEfromHessian(hessian)
    beta     <- coefs[-1]
    se       <- ses[-1]
    z_score  <- beta / se
    p_value  <- 2 * pnorm(abs(z_score), lower.tail = FALSE)

    if (is.null(feature_names)) {
      feature_names <- paste0("feature_", seq_along(beta))
    }

    data.frame(
      feature = feature_names,
      beta    = beta,
      se      = se,
      z_score = z_score,
      p_value = p_value,
      stringsAsFactors = FALSE
    )
  }, error = function(e) {
    message(sprintf("[assessment] Could not extract feature assessment: %s", e$message))
    NULL
  })
}


#' Assemble Full Results Table
#'
#' Always creates all standard columns (filled with NA if the corresponding
#' model is unavailable).
#'
#' @param ids Character vector.
#' @param z_scores Numeric vector.
#' @param p_values Numeric vector.
#' @param models Named list from \code{fit_fdrreg_models}.
#' @return Data frame.
#' @noRd
assemble_full_results <- function(ids, z_scores, p_values, models) {
  n <- length(ids)

  result <- data.frame(
    id = ids, z = z_scores, p = p_values,
    stringsAsFactors = FALSE
  )

  # Always create all FDR columns (NA by default)
  result$fdr_theo  <- rep(NA_real_, n)
  result$lfdr_theo <- rep(NA_real_, n)
  result$pep_theo  <- rep(NA_real_, n)
  result$fdr_emp   <- rep(NA_real_, n)
  result$lfdr_emp  <- rep(NA_real_, n)
  result$pep_emp   <- rep(NA_real_, n)

  # Fill from theoretical model
  if (!is.null(models$theoretical)) {
    mod <- models$theoretical
    if (!is.null(mod$FDR)) result$fdr_theo  <- mod$FDR
    if (!is.null(mod$fdr)) result$lfdr_theo <- mod$fdr
    if (!is.null(mod$pep)) result$pep_theo  <- mod$pep
  }

  # Fill from empirical model
  if (!is.null(models$empirical)) {
    mod <- models$empirical
    if (!is.null(mod$FDR)) result$fdr_emp  <- mod$FDR
    if (!is.null(mod$fdr)) result$lfdr_emp <- mod$fdr
    if (!is.null(mod$pep)) result$pep_emp  <- mod$pep
  }

  result
}