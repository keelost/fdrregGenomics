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
#'
#' @return A named list with simulated datasets:
#' \describe{
#'   \item{snp_target}{Data frame of target GWAS summary statistics.}
#'   \item{snp_aux}{Named list of auxiliary GWAS data frames.}
#'   \item{overlap_traits}{Character vector of overlapping trait names.}
#'   \item{ldsc_intercepts}{Data frame of LDSC intercepts.}
#'   \item{cov_matrix}{Pre-built covariance matrix.}
#'   \item{annotations}{Data frame of biological annotations.}
#'   \item{magma_target}{Data frame of MAGMA gene-level target results.}
#'   \item{magma_aux}{Named list of MAGMA gene-level auxiliary results.}
#'   \item{gene_annotations}{Data frame of gene-level biological annotations.}
#'   \item{spredixcan_target}{Data frame of S-PrediXcan target results.}
#'   \item{spredixcan_aux}{Named list of S-PrediXcan auxiliary results.}
#'   \item{smultixcan_target}{Data frame of S-MultiXcan target results.}
#'   \item{smultixcan_aux}{Named list of S-MultiXcan auxiliary results.}
#' }
#'
#' @examples
#' sim <- simulate_example_data(n_snps = 2000, n_genes = 200, seed = 42)
#' names(sim)
#' nrow(sim$snp_target)
#'
#' @export
simulate_example_data <- function(n_snps = 5000, n_genes = 500,
                                  n_aux = 3, n_annot = 10,
                                  prop_signal = 0.02,
                                  overlap_prop = 0.5, seed = 42) {
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
