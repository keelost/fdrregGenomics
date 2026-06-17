# ===========================================================================
# R/simulate_data.R — Simulated example data for fdrregGenomics
# ===========================================================================

#' Simulate Example Data for fdrregGenomics
#'
#' Generates simulated GWAS summary statistics, auxiliary trait data,
#' LDSC intercepts, MAGMA gene-level results, TWAS results, and
#' biological annotations for demonstration and testing.
#'
#' @param n_snps Integer, number of SNPs (default 5000).
#' @param n_genes Integer, number of genes (default 500).
#' @param n_aux Integer, number of auxiliary traits (default 3).
#' @param n_annot Integer, number of biological annotation columns (default 10).
#' @param prop_signal Numeric, proportion of true signals (default 0.02).
#' @param overlap_prop Numeric, proportion of auxiliary traits with
#'   sample overlap (default 0.5).
#' @param seed Integer, random seed (default 42).
#' @param simulation_mode Character, simulation mode:
#'   \code{"full"} (default) generates complete dataset for all tiers,
#'   \code{"summary_only"} generates z-scores and p-values only,
#'   \code{"raw_only"} generates raw data with features and response.
#' @param signal_model Character, signal probability model:
#'   \code{"simple"} (default) uses basic logistic model with random coefficients,
#'   \code{"complex"} uses a user-provided function.
#' @param signal_function Function for complex signal model (default NULL).
#'   Should accept multiple columns as separate arguments (e.g., function(x1, x2, x3)).
#' @param intercept_range Numeric vector of length 2, range for intercept in simple model (default c(-3, -1)).
#' @param coefficient_range Numeric vector of length 2, range for coefficients in simple model (default c(0.1, 1.5)).
#' @param effect_distribution List specifying effect distribution (default NULL).
#'   If NULL, uses standard normal. Otherwise, specify:
#'   - type: "normal" or "mixture"
#'   - weights: vector of mixture weights (required for mixture)
#'   - means: vector of mixture means (required for mixture)
#'   - sds: vector of mixture standard deviations (required for mixture)
#'
#' @return For simulation_mode = "full": a named list with simulated datasets.
#' For simulation_mode = "summary_only": a list with SNP and gene-level z-scores/p-values.
#' For simulation_mode = "raw_only": a list with raw data and features.
#'
#' @examples
#' # Full simulation (original behavior)
#' sim_full <- simulate_example_data(n_snps = 1000, seed = 42)
#' 
#' # Summary statistics only with complex signal model
#' sim_summary <- simulate_example_data(
#'   n_snps = 1000,
#'   simulation_mode = "summary_only",
#'   signal_model = "complex",
#'   signal_function = function(x1, x2) {
#'     -3 + 0.8*x1 + 1.0*x2
#'   }
#' )
#'
#' @export
simulate_example_data <- function(n_snps = 5000, n_genes = 500,
                                  n_aux = 3, n_annot = 10,
                                  prop_signal = 0.02,
                                  overlap_prop = 0.5, seed = 42,
                                  # New parameters
                                  simulation_mode = c("full", "summary_only", "raw_only"),
                                  signal_model = c("simple", "complex"),
                                  signal_function = NULL,
                                  intercept_range = c(-3, -1),
                                  coefficient_range = c(0.1, 1.5),
                                  effect_distribution = NULL) {
  # Match arguments
  simulation_mode <- match.arg(simulation_mode)
  signal_model <- match.arg(signal_model)
  
  # Internal helper functions
  ilogit <- function(x) 1 / (1 + exp(-x))
  rnormix <- function(n, weights, means, sds) {
    k <- length(weights)
    component <- sample(1:k, size = n, replace = TRUE, prob = weights)
    rnorm(n, mean = means[component], sd = sds[component])
  }
  
  # Call appropriate simulation function based on mode
  if (simulation_mode == "full") {
    return(simulate_full_data(n_snps, n_genes, n_aux, n_annot, 
                              prop_signal, overlap_prop, seed))
  } else if (simulation_mode == "summary_only") {
    return(simulate_summary_data(n_snps, n_genes, n_annot, seed,
                                 signal_model, signal_function,
                                 intercept_range, coefficient_range,
                                 effect_distribution))
  } else if (simulation_mode == "raw_only") {
    return(simulate_raw_data(n_snps, n_genes, n_annot, seed,
                             signal_model, signal_function,
                             intercept_range, coefficient_range,
                             effect_distribution))
  }
}

#' @keywords internal
#' @noRd
simulate_full_data <- function(n_snps, n_genes, n_aux, n_annot,
                               prop_signal, overlap_prop, seed) {
  set.seed(seed)

  # =========================================================================
  # SNP-level data
  # =========================================================================

  # Latent signal indicators
  is_signal <- rbinom(n_snps, 1, prop_signal)
  n_signal <- sum(is_signal)
  message(sprintf("[simulate] SNP-level: %d / %d true signals.", n_signal, n_snps))

  # SNP IDs and positions
  snp_ids <- paste0("rs", seq_len(n_snps))
  chr <- sample(1:22, n_snps, replace = TRUE)
  bp <- sample(1:250000000, n_snps, replace = TRUE)

  # Auxiliary trait effect sizes (correlated with signal)
  aux_effects <- matrix(runif(n_aux, 0.5, 2.0), nrow = 1)

  # Auxiliary z-scores: inflated for true signals
  aux_z <- matrix(rnorm(n_snps * n_aux), n_snps, n_aux)
  for (j in seq_len(n_aux)) {
    aux_z[, j] <- aux_z[, j] + is_signal * aux_effects[1, j]
  }

  # Target z-scores
  target_effect <- runif(1, 1.0, 3.0)
  target_z <- rnorm(n_snps) + is_signal * target_effect

  # Add correlation with overlapping auxiliaries
  n_overlap <- max(1, round(n_aux * overlap_prop))
  overlap_idx <- seq_len(n_overlap)
  for (j in overlap_idx) {
    rho <- runif(1, 0.1, 0.4)
    target_z <- target_z + rho * aux_z[, j] * is_signal * 0.3
  }

  target_p <- 2 * pnorm(abs(target_z), lower.tail = FALSE)

  # Build SNP target data frame
  snp_target <- data.frame(
    snpid = snp_ids, chr = chr, bp = bp,
    a1 = sample(c("A", "C", "G", "T"), n_snps, replace = TRUE),
    a2 = sample(c("A", "C", "G", "T"), n_snps, replace = TRUE),
    z = target_z, pval = target_p,
    stringsAsFactors = FALSE
  )

  # Build auxiliary data frames
  aux_names <- paste0("trait", seq_len(n_aux))
  snp_aux <- lapply(seq_len(n_aux), function(j) {
    data.frame(
      snpid = snp_ids, chr = chr, bp = bp,
      a1 = snp_target$a1, a2 = snp_target$a2,
      z = aux_z[, j],
      pval = 2 * pnorm(abs(aux_z[, j]), lower.tail = FALSE),
      stringsAsFactors = FALSE
    )
  })
  names(snp_aux) <- aux_names

  # Overlap traits
  overlap_traits <- aux_names[overlap_idx]

  # LDSC intercepts
  all_traits_ldsc <- c("target", overlap_traits)
  ldsc_rows <- list()
  for (i in seq_along(all_traits_ldsc)) {
    for (j in seq_along(all_traits_ldsc)) {
      if (i < j) {
        ldsc_rows <- c(ldsc_rows, list(data.frame(
          trait1 = all_traits_ldsc[i],
          trait2 = all_traits_ldsc[j],
          intercept = round(runif(1, 0.05, 0.30), 4),
          stringsAsFactors = FALSE
        )))
      }
    }
  }
  ldsc_intercepts <- do.call(rbind, ldsc_rows)

  # Pre-built covariance matrix
  n_cov <- length(all_traits_ldsc)
  cov_matrix <- diag(n_cov)
  rownames(cov_matrix) <- all_traits_ldsc
  colnames(cov_matrix) <- all_traits_ldsc
  for (i in seq_len(nrow(ldsc_intercepts))) {
    r <- match(ldsc_intercepts$trait1[i], all_traits_ldsc)
    cc <- match(ldsc_intercepts$trait2[i], all_traits_ldsc)
    cov_matrix[r, cc] <- ldsc_intercepts$intercept[i]
    cov_matrix[cc, r] <- ldsc_intercepts$intercept[i]
  }

  # SNP-level biological annotations
  annot_cols <- paste0("annot_", seq_len(n_annot))
  annot_matrix <- matrix(
    sample(0:1, n_snps * n_annot, replace = TRUE, prob = c(0.8, 0.2)),
    n_snps, n_annot
  )
  for (j in seq_len(min(3, n_annot))) {
    if (n_signal > 0) {
      annot_matrix[is_signal == 1, j] <- rbinom(n_signal, 1, 0.6)
    }
  }
  colnames(annot_matrix) <- annot_cols
  annotations <- data.frame(snpid = snp_ids, annot_matrix,
                            stringsAsFactors = FALSE)

  # =========================================================================
  # Gene-level data (MAGMA)
  # =========================================================================
  gene_ids <- paste0("ENSG", formatC(seq_len(n_genes), width = 11, flag = "0"))
  gene_names <- paste0("GENE", seq_len(n_genes))

  gene_signal <- rbinom(n_genes, 1, prop_signal * 2)
  n_gene_signal <- sum(gene_signal)
  message(sprintf("[simulate] Gene-level: %d / %d true signals.",
                  n_gene_signal, n_genes))

  gene_target_z <- rnorm(n_genes) + gene_signal * runif(1, 1.5, 3.0)
  gene_target_p <- 2 * pnorm(abs(gene_target_z), lower.tail = FALSE)

  magma_target <- data.frame(
    GENE   = gene_ids,
    CHR    = sample(1:22, n_genes, replace = TRUE),
    START  = sample(1:200000000, n_genes, replace = TRUE),
    STOP   = sample(1:200000000, n_genes, replace = TRUE),
    NSNPS  = sample(10:500, n_genes, replace = TRUE),
    NPARAM = sample(5:100, n_genes, replace = TRUE),
    N      = rep(100000, n_genes),
    ZSTAT  = gene_target_z,
    P      = gene_target_p,
    stringsAsFactors = FALSE
  )

  # Auxiliary gene-level data
  magma_aux <- lapply(seq_len(n_aux), function(j) {
    z <- rnorm(n_genes) + gene_signal * runif(1, 0.5, 1.5)
    data.frame(
      GENE  = gene_ids,
      ZSTAT = z,
      P     = 2 * pnorm(abs(z), lower.tail = FALSE),
      stringsAsFactors = FALSE
    )
  })
  names(magma_aux) <- aux_names

  # Gene-level annotations
  gene_annot <- matrix(
    sample(0:1, n_genes * n_annot, replace = TRUE, prob = c(0.7, 0.3)),
    n_genes, n_annot
  )
  for (j in seq_len(min(2, n_annot))) {
    if (n_gene_signal > 0) {
      gene_annot[gene_signal == 1, j] <- rbinom(n_gene_signal, 1, 0.7)
    }
  }
  colnames(gene_annot) <- annot_cols
  gene_annotations <- data.frame(gene = gene_ids, gene_annot,
                                 stringsAsFactors = FALSE)

  # =========================================================================
  # TWAS data (S-PrediXcan and S-MultiXcan)
  # =========================================================================

  spredixcan_target <- data.frame(
    gene = gene_ids, gene_name = gene_names,
    zscore = gene_target_z,
    effect_size = rep(NA_real_, n_genes),
    pvalue = gene_target_p,
    stringsAsFactors = FALSE
  )

  spredixcan_aux <- lapply(seq_len(n_aux), function(j) {
    z <- magma_aux[[j]]$ZSTAT
    data.frame(
      gene = gene_ids, gene_name = gene_names,
      zscore = z, effect_size = rep(NA_real_, n_genes),
      pvalue = 2 * pnorm(abs(z), lower.tail = FALSE),
      stringsAsFactors = FALSE
    )
  })
  names(spredixcan_aux) <- aux_names

  smultixcan_target <- data.frame(
    gene = gene_ids, gene_name = gene_names,
    pvalue = gene_target_p,
    n = rep(10, n_genes),
    z_mean = gene_target_z,
    stringsAsFactors = FALSE
  )

  smultixcan_aux <- lapply(seq_len(n_aux), function(j) {
    z <- magma_aux[[j]]$ZSTAT
    data.frame(
      gene = gene_ids, gene_name = gene_names,
      pvalue = 2 * pnorm(abs(z), lower.tail = FALSE),
      n = rep(10, n_genes),
      z_mean = z,
      stringsAsFactors = FALSE
    )
  })
  names(smultixcan_aux) <- aux_names

  # =========================================================================
  # Return
  # =========================================================================
  list(
    mode = "full",
    snp_target         = snp_target,
    snp_aux            = snp_aux,
    overlap_traits     = overlap_traits,
    ldsc_intercepts    = ldsc_intercepts,
    cov_matrix         = cov_matrix,
    annotations        = annotations,
    magma_target       = magma_target,
    magma_aux          = magma_aux,
    gene_annotations   = gene_annotations,
    spredixcan_target  = spredixcan_target,
    spredixcan_aux     = spredixcan_aux,
    smultixcan_target  = smultixcan_target,
    smultixcan_aux     = smultixcan_aux
  )
}

#' @keywords internal
#' @noRd
simulate_summary_data <- function(n_snps, n_genes, n_annot, seed,
                                  signal_model, signal_function,
                                  intercept_range, coefficient_range,
                                  effect_distribution) {
  set.seed(seed)
  
  # Internal helper functions
  ilogit <- function(x) 1 / (1 + exp(-x))
  rnormix <- function(n, weights, means, sds) {
    k <- length(weights)
    component <- sample(1:k, size = n, replace = TRUE, prob = weights)
    rnorm(n, mean = means[component], sd = sds[component])
  }
  
  # Generate feature matrix (annotations)
  n_total <- n_snps + n_genes
  if (n_annot > 0) {
    annotations <- matrix(rnorm(n_total * n_annot), ncol = n_annot)
    colnames(annotations) <- paste0("annot", 1:n_annot)
  } else {
    annotations <- matrix(NA, nrow = n_total, ncol = 0)
  }
  
  # Generate signal probabilities
  if (signal_model == "complex" && !is.null(signal_function)) {
  # Use user-provided complex function
  # Split matrix into individual columns
  col_list <- lapply(seq_len(ncol(annotations)), function(i) annotations[, i])
  
  # Check number of expected arguments
  func_args <- formals(signal_function)
  n_expected <- length(func_args)
  
  # Handle special case: if function has ... argument, pass all columns
  if ("..." %in% names(func_args)) {
    log_odds <- do.call(signal_function, col_list)
  } else {
    # If expected args < columns, use first n_expected columns
    if (n_expected < length(col_list)) {
      col_list <- col_list[1:n_expected]
      warning(sprintf(
        "[simulate] Signal function expects %d arguments but %d annotation columns provided. Using first %d columns.",
        n_expected, ncol(annotations), n_expected
      ))
    }
    # If expected args > columns, recycle columns
    if (n_expected > length(col_list)) {
      col_list <- rep(col_list, length.out = n_expected)
      warning(sprintf(
        "[simulate] Signal function expects %d arguments but only %d annotation columns provided. Recycling columns.",
        n_expected, ncol(annotations)
      ))
    }
    log_odds <- do.call(signal_function, col_list)
  }
  signal_prob <- ilogit(log_odds)
}
  
  # Generate true signal indicators
  is_signal <- rbinom(n_total, 1, signal_prob)
  signal_rate <- mean(is_signal)
  message(sprintf("[simulate] Signal rate: %.4f (%d/%d)", signal_rate, sum(is_signal), n_total))
  
  # Generate effect sizes
  effects <- rep(0, n_total)
  if (sum(is_signal) > 0) {
    if (!is.null(effect_distribution) && 
        !is.null(effect_distribution$weights) && 
        effect_distribution$type == "mixture") {
      # Mixture distribution
      effects[is_signal == 1] <- rnormix(
        sum(is_signal),
        weights = effect_distribution$weights,
        means = effect_distribution$means,
        sds = effect_distribution$sds
      )
    } else {
      # Standard normal distribution
      effects[is_signal == 1] <- rnorm(sum(is_signal))
    }
  }
  
  # Generate observed z-scores (effect + standard normal noise)
  z_scores <- effects + rnorm(n_total)
  p_values <- 2 * (1 - pnorm(abs(z_scores)))
  
  # Create data frame
  snp_ids <- paste0("rs", 1:n_snps)
  gene_ids <- paste0("ENSG", sprintf("%011d", 1:n_genes))
  all_ids <- c(snp_ids, gene_ids)
  
  data <- data.frame(
    id = all_ids,
    z = z_scores,
    pval = p_values,
    stringsAsFactors = FALSE
  )
  
  # Add annotations if any
  if (n_annot > 0) {
    data <- cbind(data, annotations)
  }
  
  # Create true info for evaluation
  true_info <- data.frame(
    id = all_ids,
    is_signal = is_signal,
    true_effect = effects,
    stringsAsFactors = FALSE
  )
  
  # Split SNP and gene data
  snp_data <- data[1:n_snps, ]
  gene_data <- data[(n_snps + 1):n_total, ]
  
  true_info_snp <- true_info[1:n_snps, ]
  true_info_gene <- true_info[(n_snps + 1):n_total, ]
  
  return(list(
    mode = "summary_only",
    snp = snp_data,
    gene = gene_data,
    true_info = list(snp = true_info_snp, gene = true_info_gene),
    annotations = if(n_annot > 0) annotations else NULL,
    signal_model = list(
      type = signal_model,
      prob = signal_prob,
      distribution = effect_distribution
    )
  ))
}

#' @keywords internal
#' @noRd
simulate_raw_data <- function(n_snps, n_genes, n_annot, seed,
                              signal_model, signal_function,
                              intercept_range, coefficient_range,
                              effect_distribution) {
  set.seed(seed)
  
  # Internal helper functions
  ilogit <- function(x) 1 / (1 + exp(-x))
  rnormix <- function(n, weights, means, sds) {
    k <- length(weights)
    component <- sample(1:k, size = n, replace = TRUE, prob = weights)
    rnorm(n, mean = means[component], sd = sds[component])
  }
  
  # Generate feature matrix (annotations)
  n_total <- n_snps + n_genes
  if (n_annot > 0) {
    annotations <- matrix(rnorm(n_total * n_annot), ncol = n_annot)
    colnames(annotations) <- paste0("annot", 1:n_annot)
  } else {
    annotations <- matrix(NA, nrow = n_total, ncol = 0)
  }
  
  # Generate signal probabilities
  if (signal_model == "complex" && !is.null(signal_function)) {
    # Use user-provided complex function
    col_list <- lapply(seq_len(ncol(annotations)), function(i) annotations[, i])
    log_odds <- do.call(signal_function, col_list)
    signal_prob <- ilogit(log_odds)
  } else {
    # Default simple logistic model
    beta_intercept <- runif(1, intercept_range[1], intercept_range[2])
    beta_coef <- runif(min(n_annot, 3), coefficient_range[1], coefficient_range[2])
    
    if (n_annot >= 3) {
      log_odds <- beta_intercept + 
        beta_coef[1] * annotations[, 1] + 
        beta_coef[2] * annotations[, 2] +
        beta_coef[3] * annotations[, 3]
    } else if (n_annot >= 2) {
      log_odds <- beta_intercept + 
        beta_coef[1] * annotations[, 1] + 
        beta_coef[2] * annotations[, 2]
    } else if (n_annot >= 1) {
      log_odds <- beta_intercept + 
        beta_coef[1] * annotations[, 1]
    } else {
      log_odds <- rep(beta_intercept, n_total)
    }
    signal_prob <- ilogit(log_odds)
  }
  
  # Generate true signal indicators
  is_signal <- rbinom(n_total, 1, signal_prob)
  signal_rate <- mean(is_signal)
  message(sprintf("[simulate] Signal rate: %.4f (%d/%d)", signal_rate, sum(is_signal), n_total))
  
  # Generate effect sizes
  effects <- rep(0, n_total)
  if (sum(is_signal) > 0) {
    if (!is.null(effect_distribution) && 
        !is.null(effect_distribution$weights) && 
        effect_distribution$type == "mixture") {
      effects[is_signal == 1] <- rnormix(
        sum(is_signal),
        weights = effect_distribution$weights,
        means = effect_distribution$means,
        sds = effect_distribution$sds
      )
    } else {
      effects[is_signal == 1] <- rnorm(sum(is_signal))
    }
  }
  
  # Generate response variable y
  y <- effects + rnorm(n_total)
  
  # Create feature matrix (including relevant and irrelevant features)
  n_relevant <- n_annot
  n_irrelevant <- 50 - n_annot  # Total features set to 50
  if (n_irrelevant < 0) n_irrelevant <- 0
  
  # Generate irrelevant features
  if (n_irrelevant > 0) {
    irrelevant_features <- matrix(rnorm(n_total * n_irrelevant), ncol = n_irrelevant)
    colnames(irrelevant_features) <- paste0("irrel", 1:n_irrelevant)
  } else {
    irrelevant_features <- matrix(NA, nrow = n_total, ncol = 0)
  }
  
  # Combine features
  features <- cbind(annotations, irrelevant_features)
  
  # Create data frame
  snp_ids <- paste0("rs", 1:n_snps)
  gene_ids <- paste0("ENSG", sprintf("%011d", 1:n_genes))
  
  raw_data <- data.frame(
    id = c(snp_ids, gene_ids),
    y = y,
    is_signal = is_signal,
    true_effect = effects,
    features,
    stringsAsFactors = FALSE
  )
  
  return(list(
    mode = "raw_only",
    raw = raw_data,
    features = features,
    true_info = data.frame(
      id = raw_data$id,
      is_signal = is_signal,
      true_effect = effects,
      stringsAsFactors = FALSE
    ),
    signal_model = list(
      type = signal_model,
      prob = signal_prob,
      distribution = effect_distribution
    )
  ))
}