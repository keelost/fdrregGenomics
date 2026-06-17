# ===========================================================================
# R/utils.R — Input validation, ID matching, and helper utilities
# ===========================================================================

#' Validate Required Columns
#' @param df A data frame.
#' @param required_cols Character vector of required column names.
#' @param df_name Character, name for error messages.
#' @return Invisible TRUE.
#' @noRd
validate_columns <- function(df, required_cols, df_name = "data") {
  if (!is.data.frame(df)) {
    stop(sprintf("'%s' must be a data frame.", df_name), call. = FALSE)
  }
  missing_cols <- setdiff(required_cols, colnames(df))
  if (length(missing_cols) > 0) {
    stop(sprintf("Missing columns in '%s': %s",
                 df_name, paste(missing_cols, collapse = ", ")),
         call. = FALSE)
  }
  invisible(TRUE)
}

#' Match and Align Datasets by Common ID
#' @param target A data frame.
#' @param aux_list A named list of data frames.
#' @param id_col Character, ID column name.
#' @return List with aligned target and aux.
#' @noRd
match_and_align <- function(target, aux_list, id_col) {
  if (is.null(aux_list) || length(aux_list) == 0) {
    return(list(target = target, aux = list()))
  }
  target_ids <- as.character(target[[id_col]])
  all_ids <- list(target_ids)
  for (i in seq_along(aux_list)) {
    all_ids[[i + 1]] <- as.character(aux_list[[i]][[id_col]])
  }
  common_ids <- Reduce(intersect, all_ids)
  if (length(common_ids) == 0) {
    stop("No common IDs found between target and auxiliary datasets.",
         call. = FALSE)
  }
  n_target <- nrow(target)
  if (length(common_ids) < n_target) {
    message(sprintf("[match_and_align] %d / %d target IDs retained.",
                    length(common_ids), n_target))
  }
  target <- target[match(common_ids, target_ids), , drop = FALSE]
  rownames(target) <- NULL
  aux_aligned <- lapply(aux_list, function(df) {
    ids <- as.character(df[[id_col]])
    sub <- df[match(common_ids, ids), , drop = FALSE]
    rownames(sub) <- NULL
    sub
  })
  list(target = target, aux = aux_aligned)
}

#' Align Annotations to Target IDs
#'
#' Supports different ID column names in target vs annotations via
#' the \code{annot_id_col} parameter. Falls back to case-insensitive
#' matching if the exact column is not found.
#'
#' @param target_ids Character vector of target IDs.
#' @param annotations Data frame of annotations.
#' @param id_col Character, ID column name used in target.
#' @param annot_id_col Character, ID column name in annotations.
#'   Defaults to \code{id_col}. If not found, tries case-insensitive match.
#' @return Numeric matrix or NULL.
#' @noRd
align_annotations <- function(target_ids, annotations, id_col,
                              annot_id_col = NULL) {
  if (is.null(annotations)) return(NULL)

  # Resolve annotation ID column
  if (is.null(annot_id_col)) annot_id_col <- id_col

  if (!annot_id_col %in% colnames(annotations)) {
    # Case-insensitive fallback
    lower_cols   <- tolower(colnames(annotations))
    lower_target <- tolower(annot_id_col)
    idx <- which(lower_cols == lower_target)
    if (length(idx) > 0) {
      annot_id_col <- colnames(annotations)[idx[1]]
      message(sprintf(
        "[align_annotations] Annotation ID column '%s' not found; using '%s' (case-insensitive match).",
        id_col, annot_id_col
      ))
    } else {
      stop(sprintf(
        "Annotation ID column '%s' not found.\nAvailable columns: %s",
        annot_id_col, paste(colnames(annotations), collapse = ", ")
      ), call. = FALSE)
    }
  }

  annot_ids <- as.character(annotations[[annot_id_col]])
  idx <- match(target_ids, annot_ids)
  na_count <- sum(is.na(idx))

  if (na_count > 0) {
    message(sprintf(
      "[align_annotations] %d / %d target IDs have no matching annotation; filled with 0.",
      na_count, length(target_ids)
    ))
  }

  # Keep only numeric columns (excluding the ID column)
  numeric_cols <- names(annotations)[sapply(annotations, is.numeric)]
  numeric_cols <- setdiff(numeric_cols, annot_id_col)

  if (length(numeric_cols) == 0) {
    warning("No numeric annotation columns found.", call. = FALSE)
    return(NULL)
  }

  annot_mat <- as.matrix(annotations[match(target_ids, annot_ids),
                                     numeric_cols, drop = FALSE])
  annot_mat[is.na(annot_mat)] <- 0
  rownames(annot_mat) <- target_ids
  annot_mat
}

#' Record Reproducibility Parameters
#' @param seed Integer.
#' @param call The call object.
#' @return List.
#' @noRd
record_params <- function(seed, call) {
  list(
    call = deparse(call, width.cutoff = 500),
    seed = seed,
    r_version = R.version.string,
    timestamp = format(Sys.time(), "%Y-%m-%d %H:%M:%S"),
    package_versions = tryCatch(
      list(
        fdrregGenomics = as.character(packageVersion("fdrregGenomics")),
        FDRreg = as.character(packageVersion("FDRreg")),
        powerplus = as.character(packageVersion("powerplus")),
        glmnet = as.character(packageVersion("glmnet"))
      ),
      error = function(e) list()
    )
  )
}

#' Extract Z-Scores from Data Frame
#' @param df Data frame.
#' @param z_col Column name for z-score (NULL if not available).
#' @param p_col Column name for p-value.
#' @param random_sign Logical.
#' @param seed Integer.
#' @return Numeric vector.
#' @noRd
extract_z_scores <- function(df, z_col = NULL, p_col = "pvalue",
                             random_sign = TRUE, seed = 42) {
  if (!is.null(z_col) && z_col %in% colnames(df)) {
    return(df[[z_col]])
  }
  if (!p_col %in% colnames(df)) {
    stop("Neither z-score column '", z_col, "' nor p-value column '",
         p_col, "' found.", call. = FALSE)
  }
  pvals <- df[[p_col]]
  pvals[is.na(pvals)] <- 1
  pvals <- pmax(pmin(pvals, 1), .Machine$double.eps)
  z_mag <- abs(qnorm(pvals / 2))
  if (random_sign) {
    set.seed(seed)
    signs <- sample(c(-1, 1), length(z_mag), replace = TRUE)
    return(z_mag * signs)
  }
  -z_mag
}

#' Standard FDR Thresholds
#' @return Numeric vector.
#' @noRd
standard_thresholds <- function() {
  c(0.5, 0.4, 0.3, 0.2, 0.1, 0.05, 0.04, 0.03, 0.02, 0.01,
    0.001, 5e-04, 5e-06, 5e-08)
}

#' Count Discoveries at Thresholds
#' @param fdr_vec Numeric vector.
#' @param thresholds Numeric vector.
#' @return Integer vector.
#' @noRd
count_discoveries <- function(fdr_vec, thresholds) {
  vapply(thresholds, function(th) {
    as.integer(sum(fdr_vec < th, na.rm = TRUE))
  }, integer(1))
}

#' Evaluate FDR Control Performance
#'
#' Calculate False Discovery Proportion (FDP), power, and F1 score
#' for different FDR thresholds.
#'
#' @param fdr_values Numeric vector of FDR values from FDRreg.
#' @param true_signals Logical vector indicating true signals.
#' @param thresholds Numeric vector of FDR thresholds (default c(0.01, 0.05, 0.1, 0.2)).
#'
#' @return Data frame with performance metrics for each threshold.
#'
#' @examples
#' # Example usage with simulated data
#' fdr_vals <- runif(100, 0, 0.3)
#' true_sigs <- rbinom(100, 1, 0.1) == 1
#' perf <- evaluate_fdr_performance(fdr_vals, true_sigs)
#' print(perf)
#'
#' @export
evaluate_fdr_performance <- function(fdr_values, true_signals, 
                                     thresholds = c(0.01, 0.05, 0.1, 0.2)) {
  # Ensure inputs are vectors
  fdr_values <- as.vector(fdr_values)
  true_signals <- as.logical(true_signals)
  
  # Initialize results data frame
  results <- data.frame(
    threshold = thresholds,
    FDP = NA_real_,
    power = NA_real_,
    F1 = NA_real_,
    stringsAsFactors = FALSE
  )
  
  # Calculate metrics for each threshold
  for (i in seq_along(thresholds)) {
    threshold <- thresholds[i]
    findings <- which(fdr_values <= threshold)
    
    if (length(findings) > 0) {
      TP <- sum(true_signals[findings])
      FP <- sum(!true_signals[findings])
      FN <- sum(true_signals) - TP
      
      precision <- TP / (TP + FP)
      recall <- TP / (TP + FN)
      
      results$FDP[i] <- FP / length(findings)
      results$power[i] <- recall
      results$F1[i] <- ifelse(precision + recall > 0, 
                               2 * (precision * recall) / (precision + recall), 
                               0)
    }
  }
  
  return(results)
}

#' Evaluate Variable Selection Performance
#'
#' Calculate precision, recall, and F1 score for variable selection results.
#'
#' @param selected_vars Vector of selected variable indices or logical vector.
#' @param true_vars Vector of true variable indices or logical vector.
#' @param total_vars Total number of variables (optional, required if indices are used).
#'
#' @return List with performance metrics including precision, recall, F1,
#'   and basic counts (TP, FP, FN, TN).
#'
#' @examples
#' # Example usage
#' selected <- sample(1:100, 15)
#' true <- sample(1:100, 20)
#' perf <- evaluate_variable_selection(selected, true, total_vars = 100)
#' cat("Precision:", perf$precision, "\n")
#' cat("Recall:", perf$recall, "\n")
#'
#' @export
evaluate_variable_selection <- function(selected_vars, true_vars, total_vars = NULL) {
  # Convert indices to logical vectors if needed
  if (is.numeric(selected_vars) && !is.null(total_vars)) {
    selected_logical <- rep(FALSE, total_vars)
    selected_logical[selected_vars] <- TRUE
  } else {
    selected_logical <- as.logical(selected_vars)
  }
  
  if (is.numeric(true_vars) && !is.null(total_vars)) {
    true_logical <- rep(FALSE, total_vars)
    true_logical[true_vars] <- TRUE
  } else {
    true_logical <- as.logical(true_vars)
  }
  
  # Calculate basic counts
  TP <- sum(selected_logical & true_logical)
  FP <- sum(selected_logical & !true_logical)
  FN <- sum(!selected_logical & true_logical)
  TN <- sum(!selected_logical & !true_logical)
  
  # Calculate metrics
  precision <- ifelse(TP + FP > 0, TP / (TP + FP), 0)
  recall <- ifelse(TP + FN > 0, TP / (TP + FN), 0)
  F1 <- ifelse(precision + recall > 0, 
               2 * (precision * recall) / (precision + recall), 
               0)
  
  return(list(
    precision = precision,
    recall = recall,
    F1 = F1,
    TP = TP,
    FP = FP,
    FN = FN,
    TN = TN,
    support = sum(selected_logical),  # Number of selected variables
    true_count = sum(true_logical)    # Number of true variables
  ))
}