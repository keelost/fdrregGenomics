# ===========================================================================
# R/varselect.R — Variable selection: LASSO, Marginal, Elastic Net
# ===========================================================================

#' Perform Variable Selection
#' @param target_z Numeric vector.
#' @param features Numeric matrix.
#' @param method Character.
#' @param seed Integer.
#' @return List with features, selected_cols, model.
#' @noRd
perform_variable_selection <- function(target_z, features,
                                       method = c("none", "lasso", "marginal",
                                                  "elasticnet"),
                                       seed = 42) {
  method <- match.arg(method)
  if (method == "none") {
    return(list(features = features, selected_cols = colnames(features),
                model = NULL))
  }
  set.seed(seed)
  result <- switch(method,
    lasso      = vselect_lasso(target_z, features),
    marginal   = vselect_marginal(target_z, features),
    elasticnet = vselect_elasticnet(target_z, features)
  )
  if (length(result$selected_cols) == 0) {
    warning("Variable selection dropped all features. Using all.", call. = FALSE)
    return(list(features = features, selected_cols = colnames(features),
                model = result$model))
  }
  message(sprintf("[var_select] %s: %d / %d features retained.",
                  method, length(result$selected_cols), ncol(features)))
  result
}

#' @noRd
vselect_lasso <- function(target_z, features, s = "lambda.min") {
  cv_fit <- cv.glmnet(x = features, y = abs(target_z),
                      family = "gaussian", alpha = 1, nlambda = 50,
                      standardize = TRUE)
  coefs <- as.matrix(coef(cv_fit, s = s))
  coef_vals <- coefs[-1, , drop = FALSE]
  selected <- rownames(coef_vals)[coef_vals[, 1] != 0]
  list(features = features[, selected, drop = FALSE],
       selected_cols = selected, model = cv_fit)
}

#' @noRd
vselect_elasticnet <- function(target_z, features, s = "lambda.min") {
  cv_fit <- cv.glmnet(x = features, y = abs(target_z),
                      family = "gaussian", alpha = 0.5, nlambda = 50,
                      standardize = TRUE)
  coefs <- as.matrix(coef(cv_fit, s = s))
  coef_vals <- coefs[-1, , drop = FALSE]
  selected <- rownames(coef_vals)[coef_vals[, 1] != 0]
  list(features = features[, selected, drop = FALSE],
       selected_cols = selected, model = cv_fit)
}

#' @noRd
vselect_marginal <- function(target_z, features, top_k = 0.5) {
  abs_z_target <- abs(target_z)
  cors <- apply(features, 2, function(col) {
    cor(abs_z_target, abs(col), use = "complete.obs")
  })
  cors[is.na(cors)] <- 0
  if (top_k < 1) top_k <- max(1, round(top_k * ncol(features)))
  top_k <- min(top_k, ncol(features))
  selected <- names(sort(abs(cors), decreasing = TRUE))[seq_len(top_k)]
  list(features = features[, selected, drop = FALSE],
       selected_cols = selected, model = list(correlations = cors, top_k = top_k))
}
