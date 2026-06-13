#' Run FDRreg on TWAS Results (S-PrediXcan / S-MultiXcan)
#'
#' @param target Data frame of TWAS results.
#' @param twas_type Character, "smultixcan" or "spredixcan".
#' @param aux Named list of TWAS data frames.
#' @param annotations Data frame of biological annotations.
#' @param id_col Character (default "gene").
#' @param z_col Character (default "zscore").
#' @param p_col Character (default "pvalue").
#' @param annot_id_col Character, ID column in annotations (default NULL).
#' @param feature_transform Character, "abs" (default), "signed", "split".
#' @param var_select Character, "none" (default), "lasso", "marginal", "elasticnet".
#' @param fdrreg_nulltype Character, "both" (default).
#' @param fdrreg_method Character, "pr" (default).
#' @param seed Integer.
#' @return An object of class \code{fdrreg_result}.
#' @export
run_fdrreg_twas <- function(target, twas_type = c("smultixcan", "spredixcan"),
                            aux = list(), annotations = NULL,
                            id_col = "gene", z_col = "zscore", p_col = "pvalue",
                            annot_id_col = NULL,
                            feature_transform = c("abs", "signed", "split"),
                            var_select = c("none", "lasso", "marginal", "elasticnet"),
                            fdrreg_nulltype = c("both", "theoretical", "empirical"),
                            fdrreg_method = "pr", seed = 42) {
  twas_type         <- match.arg(twas_type)
  feature_transform <- match.arg(feature_transform)
  var_select        <- match.arg(var_select)
  fdrreg_nulltype   <- match.arg(fdrreg_nulltype)
  call <- match.call()

  if (twas_type == "spredixcan") {
    validate_columns(target, c(id_col, z_col, p_col), "target")
  } else {
    validate_columns(target, c(id_col, p_col), "target")
  }
  for (nm in names(aux)) {
    req_cols <- if (twas_type == "spredixcan") c(id_col, z_col) else c(id_col, p_col)
    validate_columns(aux[[nm]], req_cols, sprintf("aux$%s", nm))
  }

  aligned <- match_and_align(target, aux, id_col)
  target  <- aligned$target
  aux     <- aligned$aux
  message(sprintf("[run_fdrreg_twas] %d genes after alignment.", nrow(target)))

  target_z <- extract_z_scores(target, z_col = z_col, p_col = p_col,
                               random_sign = TRUE, seed = seed)

  aux_features <- if (length(aux) > 0) {
    mats <- lapply(aux, function(df) {
      abs(extract_z_scores(df, z_col = z_col, p_col = p_col, random_sign = FALSE))
    })
    out <- do.call(cbind, mats)
    colnames(out) <- names(aux)
    out
  } else NULL

  annot_mat <- NULL
  if (!is.null(annotations)) {
    annot_mat <- align_annotations(target[[id_col]], annotations,
                                   id_col, annot_id_col)
  }

  features     <- build_combined_features(aux_features, annot_mat)
  features     <- transform_features(features, feature_transform)
  covariates_used <- colnames(features)

  vs_result    <- perform_variable_selection(target_z, features,
                                             method = var_select, seed = seed)
  features     <- vs_result$features
  covariates_used <- vs_result$selected_cols

  set.seed(seed)
  fit <- fit_fdrreg_models(target_z, features,
                           nulltype = fdrreg_nulltype, method = fdrreg_method)
  models      <- fit$models
  fit_status  <- fit$fit_status

  full_results <- assemble_full_results(
    ids = target[[id_col]], z_scores = target_z,
    p_values = target[[p_col]], models = models
  )

  thresholds   <- standard_thresholds()
  summary_list <- list(threshold = thresholds)
  if (!is.null(models$theoretical) && !is.null(models$theoretical$FDR))
    summary_list$fdr_theo <- count_discoveries(models$theoretical$FDR, thresholds)
  if (!is.null(models$empirical) && !is.null(models$empirical$FDR))
    summary_list$fdr_emp <- count_discoveries(models$empirical$FDR, thresholds)

  assessment_model <- if (!is.null(models$theoretical)) models$theoretical else models$empirical
  params <- record_params(seed, call)

  new_fdrreg_result(tier = paste0("twas_", twas_type), full_results = full_results,
                    summary_counts = as.data.frame(summary_list, stringsAsFactors = FALSE),
                    assessment = extract_model_assessment(assessment_model, covariates_used),
                    models = models, varselect_model = vs_result$model,
                    covariates_used = covariates_used, params = params,
                    fit_status = fit_status)
}