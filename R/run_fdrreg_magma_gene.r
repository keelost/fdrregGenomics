#' Run FDRreg at the MAGMA Gene Level
#'
#' @param target Data frame of MAGMA gene-level results.
#' @param aux Named list of gene-level data frames from other traits.
#' @param annotations Data frame of biological annotations.
#' @param id_col Character, gene ID column (default "GENE").
#' @param z_col Character, z-score column (default "ZSTAT").
#' @param p_col Character, p-value column (default "P").
#' @param annot_id_col Character, ID column in annotations (default NULL
#'   = same as \code{id_col}). Use when annotations use a different ID
#'   column name (e.g., "gene" vs "GENE").
#' @param feature_transform Character, "abs" (default), "signed", "split".
#' @param var_select Character, "none" (default), "lasso", "marginal", "elasticnet".
#' @param fdrreg_nulltype Character, "both" (default).
#' @param fdrreg_method Character, "pr" (default).
#' @param seed Integer.
#' @return An object of class \code{fdrreg_result}.
#' @export
run_fdrreg_magma_gene <- function(target, aux = list(), annotations = NULL,
                                  id_col = "GENE", z_col = "ZSTAT",
                                  p_col = "P", annot_id_col = NULL,
                                  feature_transform = c("abs", "signed", "split"),
                                  var_select = c("none", "lasso", "marginal", "elasticnet"),
                                  fdrreg_nulltype = c("both", "theoretical", "empirical"),
                                  fdrreg_method = "pr", seed = 42) {
  feature_transform <- match.arg(feature_transform)
  var_select        <- match.arg(var_select)
  fdrreg_nulltype   <- match.arg(fdrreg_nulltype)
  call <- match.call()

  validate_columns(target, c(id_col, z_col, p_col), "target")
  for (nm in names(aux)) {
    validate_columns(aux[[nm]], c(id_col, z_col), sprintf("aux$%s", nm))
  }

  aligned <- match_and_align(target, aux, id_col)
  target  <- aligned$target
  aux     <- aligned$aux
  message(sprintf("[run_fdrreg_magma_gene] %d genes after alignment.", nrow(target)))

  target_z <- target[[z_col]]

  aux_features <- if (length(aux) > 0) {
    out <- do.call(cbind, lapply(aux, `[[`, z_col))
    colnames(out) <- names(aux)
    out
  } else NULL

  annot_mat <- NULL
  if (!is.null(annotations)) {
    annot_mat <- align_annotations(target[[id_col]], annotations,
                                   id_col, annot_id_col)
  }

  features <- build_combined_features(aux_features, annot_mat)
  features <- transform_features(features, feature_transform)
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

  new_fdrreg_result(tier = "magma_gene", full_results = full_results,
                    summary_counts = as.data.frame(summary_list, stringsAsFactors = FALSE),
                    assessment = extract_model_assessment(assessment_model, covariates_used),
                    models = models, varselect_model = vs_result$model,
                    covariates_used = covariates_used, params = params,
                    fit_status = fit_status)
}