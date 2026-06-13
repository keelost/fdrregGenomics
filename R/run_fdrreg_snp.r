#' Run FDRreg at the SNP Level
#'
#' Performs FDRreg on SNP-level GWAS summary statistics with optional
#' sample-overlap decorrelation, biological annotations, and variable selection.
#'
#' @param target Data frame of target GWAS. Must contain columns
#'   specified by \code{id_col}, \code{z_col}, \code{p_col}.
#' @param aux Named list of data frames (auxiliary GWAS).
#' @param overlap_traits Character vector of names in \code{aux} that
#'   have sample overlap with the target.
#' @param cov_matrix Pre-built covariance matrix. If NULL and
#'   \code{ldsc_intercepts} is provided, built from LDSC.
#' @param ldsc_intercepts Data frame with LDSC intercept values.
#' @param target_name Character, target trait name (default "target").
#' @param annotations Data frame of biological annotations.
#' @param id_col Character, SNP ID column (default "snpid").
#' @param z_col Character, z-score column (default "z").
#' @param p_col Character, p-value column (default "pval").
#' @param ldsc_col Named list: trait1, trait2, intercept column names.
#' @param feature_transform Character, "signed" (default), "abs", or "split".
#' @param var_select Character, "none" (default), "lasso", "marginal",
#'   or "elasticnet".
#' @param fdrreg_nulltype Character, "both" (default), "theoretical",
#'   or "empirical".
#' @param fdrreg_method Character, passed to FDRreg (default "pr").
#' @param seed Integer, random seed (default 42).
#'
#' @return An object of class \code{fdrreg_result}.
#'
#' @examples
#' \dontrun{
#' sim <- simulate_example_data(n_snps = 2000, seed = 42)
#' result <- run_fdrreg_snp(
#'   target = sim$snp_target,
#'   aux = sim$snp_aux,
#'   overlap_traits = sim$overlap_traits,
#'   ldsc_intercepts = sim$ldsc_intercepts
#' )
#' print(result)
#' significant(result, threshold = 0.05)
#' }
#'
#' @export
run_fdrreg_snp <- function(target,
                           aux = list(),
                           overlap_traits = NULL,
                           cov_matrix = NULL,
                           ldsc_intercepts = NULL,
                           target_name = "target",
                           annotations = NULL,
                           id_col = "snpid",
                           z_col = "z",
                           p_col = "pval",
                           ldsc_col = list(trait1 = "trait1",
                                           trait2 = "trait2",
                                           intercept = "intercept"),
                           feature_transform = c("signed", "abs", "split"),
                           var_select = c("none", "lasso", "marginal",
                                          "elasticnet"),
                           fdrreg_nulltype = c("both", "theoretical",
                                               "empirical"),
                           fdrreg_method = "pr",
                           seed = 42) {
  feature_transform <- match.arg(feature_transform)
  var_select <- match.arg(var_select)
  fdrreg_nulltype <- match.arg(fdrreg_nulltype)
  call <- match.call()

  # Step 0: Validate
  validate_columns(target, c(id_col, z_col, p_col), "target")
  for (nm in names(aux)) {
    validate_columns(aux[[nm]], c(id_col, z_col), sprintf("aux$%s", nm))
  }

  # Step 1: Align
  aligned <- match_and_align(target, aux, id_col)
  target <- aligned$target
  aux <- aligned$aux
  n_snps <- nrow(target)
  message(sprintf("[run_fdrreg_snp] %d SNPs after alignment.", n_snps))

  # Step 2: Separate overlap / no-overlap
  all_aux_names <- names(aux)
  if (is.null(overlap_traits)) overlap_traits <- character(0)
  no_overlap_names <- setdiff(all_aux_names, overlap_traits)
  overlap_list <- aux[intersect(overlap_traits, all_aux_names)]
  no_overlap_list <- aux[intersect(no_overlap_names, all_aux_names)]

  # Step 3: Decorrelation
  target_z <- target[[z_col]]
  if (length(overlap_list) > 0) {
    if (is.null(cov_matrix) && is.null(ldsc_intercepts)) {
      stop("Either 'cov_matrix' or 'ldsc_intercepts' required when overlap_traits is non-empty.",
           call. = FALSE)
    }
    if (is.null(cov_matrix)) {
      cov_matrix <- build_ldsc_cov_matrix(
        target_name, names(overlap_list), ldsc_intercepts,
        col_trait1 = ldsc_col$trait1, col_trait2 = ldsc_col$trait2,
        col_intercept = ldsc_col$intercept
      )
    }
    z_mat <- cbind(target_z,
                   do.call(cbind, lapply(overlap_list, `[[`, z_col)))
    z_decor <- decorrelate_z_scores(z_mat, cov_matrix)
    target_z <- z_decor[1, ]
    overlap_z_decor <- t(z_decor[-1, , drop = FALSE])
    colnames(overlap_z_decor) <- names(overlap_list)
  } else {
    overlap_z_decor <- NULL
  }

  # Step 4: Construct features
  no_overlap_z <- if (length(no_overlap_list) > 0) {
    out <- do.call(cbind, lapply(no_overlap_list, `[[`, z_col))
    colnames(out) <- names(no_overlap_list)
    out
  } else NULL

  aux_features <- cbind(overlap_z_decor, no_overlap_z)
  annot_mat <- NULL
  if (!is.null(annotations)) {
    annot_mat <- align_annotations(target[[id_col]], annotations, id_col)
  }
  features <- build_combined_features(aux_features, annot_mat)

  # Step 5: Transform
  features <- transform_features(features, feature_transform)
  covariates_used <- colnames(features)

  # Step 6: Variable selection
  vs_result <- perform_variable_selection(target_z, features,
                                          method = var_select, seed = seed)
  features <- vs_result$features
  covariates_used <- vs_result$selected_cols

  # Step 7: FDRreg
  set.seed(seed)
  models <- fit_fdrreg_models(target_z, features,
                              nulltype = fdrreg_nulltype,
                              method = fdrreg_method)

  # Step 8: Assemble
  full_results <- assemble_full_results(
    ids = target[[id_col]], z_scores = target_z,
    p_values = target[[p_col]], models = models
  )

  thresholds <- standard_thresholds()
  summary_list <- list(threshold = thresholds)
  if (!is.null(models$theoretical))
    summary_list$fdr_theo <- count_discoveries(models$theoretical$FDR, thresholds)
  if (!is.null(models$empirical))
    summary_list$fdr_emp <- count_discoveries(models$empirical$FDR, thresholds)
  summary_counts <- as.data.frame(summary_list, stringsAsFactors = FALSE)

  assessment_model <- if (!is.null(models$theoretical)) models$theoretical else models$empirical
  assessment <- extract_model_assessment(assessment_model, covariates_used)
  params <- record_params(seed, call)

  new_fdrreg_result(tier = "snp", full_results = full_results,
                    summary_counts = summary_counts, assessment = assessment,
                    models = models, varselect_model = vs_result$model,
                    covariates_used = covariates_used, params = params)
}
