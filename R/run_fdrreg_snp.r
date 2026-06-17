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
#' @param annot_id_col Character, ID column name in \code{annotations}.
#'   Defaults to \code{id_col}. Useful when annotations use a different
#'   ID column name (e.g., "snpid" vs "SNP").
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
                           annot_id_col = NULL,
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
  var_select        <- match.arg(var_select)
  fdrreg_nulltype   <- match.arg(fdrreg_nulltype)
  call <- match.call()

  # ============================================================
  # NEW: Auto-detect new simulation mode and adjust parameters
  # ============================================================
  # More robust detection based on column names
  is_new_mode <- FALSE
  
  if (is.data.frame(target)) {
    # Check for simulation_mode attribute
    if (!is.null(attr(target, "simulation_mode")) && 
        attr(target, "simulation_mode") == "summary_only") {
      is_new_mode <- TRUE
    }
    # Check for column patterns
    else if ("id" %in% colnames(target) && 
             !"snpid" %in% colnames(target) &&
             "z" %in% colnames(target) &&
             "pval" %in% colnames(target)) {
      # Check for annotation columns
      annot_cols <- grep("^annot\\d+$", colnames(target), value = TRUE)
      if (length(annot_cols) > 0) {
        is_new_mode <- TRUE
      }
    }
  }
  
  if (is_new_mode) {
    # Auto-adjust id_col for new simulation mode
    if (id_col == "snpid" && !"snpid" %in% colnames(target)) {
      id_col <- "id"
      message("[run_fdrreg_snp] Auto-detected new simulation mode. Using 'id' as target ID column.")
    }
    
    # For new mode, aux is typically empty
    if (missing(aux) || is.null(aux)) {
      aux <- list()
      message("[run_fdrreg_snp] No auxiliary traits provided. Using only target and annotations.")
    }
    
    # Auto-set annot_id_col if not provided
    if (is.null(annot_id_col) && !is.null(annotations)) {
      # Check if annotations have ID column
      if (id_col %in% colnames(annotations)) {
        annot_id_col <- id_col
        message(sprintf("[run_fdrreg_snp] Auto-matched annotation ID column: '%s'", id_col))
      } else if ("id" %in% colnames(annotations)) {
        annot_id_col <- "id"
        message("[run_fdrreg_snp] Auto-matched annotation ID column: 'id'")
      }
    }
    
    # If annotations are NULL, extract from target data
    if (is.null(annotations)) {
      # Check if target has annotation columns
      annot_cols <- grep("^annot\\d+$", colnames(target), value = TRUE)
      if (length(annot_cols) > 0) {
        annotations <- target[, c(id_col, annot_cols), drop = FALSE]
        annot_id_col <- id_col
        message("[run_fdrreg_snp] Extracted annotations from target data.")
      }
    }
  }
  # ============================================================
  # END NEW CODE
  # ============================================================

  # Step 0: Validate
  validate_columns(target, c(id_col, z_col, p_col), "target")
  for (nm in names(aux)) {
    validate_columns(aux[[nm]], c(id_col, z_col), sprintf("aux$%s", nm))
  }

  # Step 1: Align
  aligned <- match_and_align(target, aux, id_col)
  target  <- aligned$target
  aux     <- aligned$aux
  message(sprintf("[run_fdrreg_snp] %d SNPs after alignment.", nrow(target)))

  # Step 2: Separate overlap / no-overlap
  all_aux_names   <- names(aux)
  if (is.null(overlap_traits)) overlap_traits <- character(0)
  no_overlap_names <- setdiff(all_aux_names, overlap_traits)
  overlap_list     <- aux[intersect(overlap_traits, all_aux_names)]
  no_overlap_list  <- aux[intersect(no_overlap_names, all_aux_names)]

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
    z_mat  <- cbind(target_z, do.call(cbind, lapply(overlap_list, `[[`, z_col)))
    z_decor <- decorrelate_z_scores(z_mat, cov_matrix)
    target_z         <- z_decor[1, ]
    overlap_z_decor  <- t(z_decor[-1, , drop = FALSE])
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
  annot_mat    <- NULL
  if (!is.null(annotations)) {
    annot_mat <- align_annotations(target[[id_col]], annotations,
                                   id_col, annot_id_col)
  }
  features <- build_combined_features(aux_features, annot_mat)
  features <- transform_features(features, feature_transform)
  covariates_used <- colnames(features)

  # ============================================================
  # NEW: Handle case when no features are provided
  # ============================================================
  if (is.null(features) || ncol(features) == 0) {
    warning("[run_fdrreg_snp] No features available for FDRreg analysis. ",
            "Consider providing 'aux' or 'annotations'.",
            call. = FALSE)
    # Create minimal result with only target z-scores
    full_results <- data.frame(
      id = target[[id_col]],
      z = target_z,
      p = target[[p_col]],
      fdr_theo = NA_real_,
      lfdr_theo = NA_real_,
      pep_theo = NA_real_,
      fdr_emp = NA_real_,
      lfdr_emp = NA_real_,
      pep_emp = NA_real_,
      stringsAsFactors = FALSE
    )
    
    # Return minimal result
    return(new_fdrreg_result(
      tier = "snp",
      full_results = full_results,
      summary_counts = data.frame(threshold = standard_thresholds(),
                                  fdr_theo = 0, fdr_emp = 0),
      assessment = data.frame(beta = NA, se = NA, z_score = NA, p_value = NA),
      models = list(theoretical = NULL, empirical = NULL),
      varselect_model = NULL,
      covariates_used = character(0),
      params = record_params(seed, call),
      fit_status = "no_features"
    ))
  }
  # ============================================================
  # END NEW CODE
  # ============================================================

  # Step 5: Variable selection
  vs_result    <- perform_variable_selection(target_z, features,
                                             method = var_select, seed = seed)
  features     <- vs_result$features
  covariates_used <- vs_result$selected_cols

  # Step 6: FDRreg
  set.seed(seed)
  fit <- fit_fdrreg_models(target_z, features,
                           nulltype = fdrreg_nulltype,
                           method = fdrreg_method)
  models      <- fit$models
  fit_status  <- fit$fit_status

  # Step 7: Assemble output
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
  summary_counts <- as.data.frame(summary_list, stringsAsFactors = FALSE)

  assessment_model <- if (!is.null(models$theoretical)) models$theoretical else models$empirical
  assessment <- extract_model_assessment(assessment_model, covariates_used)
  params     <- record_params(seed, call)

  new_fdrreg_result(tier = "snp", full_results = full_results,
                    summary_counts = summary_counts, assessment = assessment,
                    models = models, varselect_model = vs_result$model,
                    covariates_used = covariates_used, params = params,
                    fit_status = fit_status)
}