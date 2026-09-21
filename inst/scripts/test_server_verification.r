#!/usr/bin/env Rscript
# ==============================================================================
# Comprehensive Smoke & Integration Test for fdrregGenomics v0.3.0
# Tests all three tiers (SNP, MAGMA Gene, TWAS), annotation modes, and scoring.
# ==============================================================================

cat("================================================================\n")
cat("Starting fdrregGenomics Server Verification Test\n")
cat("================================================================\n\n")

# 1. Check Package Loading
cat("[1/6] Testing package loading and version check...\n")
suppressPackageStartupMessages({
  library(fdrregGenomics)
})
pkg_ver <- packageVersion("fdrregGenomics")
cat(sprintf("  -> Successfully loaded fdrregGenomics version %s\n\n", pkg_ver))

# 2. Test Annotation Loading (Option A: Curated Benchmark)
cat("[2/6] Testing Option A: Curated benchmark annotations (Supplementary Table S1.21)...\n")
annot_entrez <- load_builtin_annotations(key_type = "entrez", version = "curated")
annot_ensembl <- load_builtin_annotations(key_type = "ensembl", version = "curated")
annot_symbol  <- load_builtin_annotations(key_type = "symbol", version = "curated")

stopifnot(nrow(annot_entrez) == 19503)
stopifnot(ncol(annot_entrez) == 22) # ID + 21 bio features
stopifnot("dise.brain" %in% colnames(annot_entrez))
cat(sprintf("  -> Curated annotations loaded: %d genes x %d features\n",
            nrow(annot_entrez), ncol(annot_entrez) - 1))
cat("  -> Key types verified: Entrez, Ensembl, Symbol\n\n")

# 3. Test Annotation Scoring Engine (Option B: Automated Fresh Scoring)
cat("[3/6] Testing Option B: Automated deterministic scoring engine (score_david_annotations)...\n")
raw_toy <- data.frame(
  ID = c("1001", "1002", "1003", "1004"),
  tissue = c("Fetal brain, cerebellum", "Liver", NA, "Hippocampus, synapse"),
  disease = c("Schizophrenia, Bipolar", "Hypertension", NA, "Drug abuse"),
  pathway = c("Dopaminergic synapse", "Glycolysis", NA, "GABAergic synapse"),
  tfbs_count = c(5, 1, NA, "siteA, siteB, siteC"),
  process = c("Synaptic transmission", "Cell cycle", NA, "Circadian regulation"),
  interaction = c("Dopamine receptor binding", "Insulin receptor", NA, "Serotonin transporter"),
  stringsAsFactors = FALSE
)

scored_df <- score_david_annotations(raw_toy, id_col = "ID")
stopifnot(identical(scored_df$expre.brain, c(1L, 0L, 0L, 1L)))
stopifnot(identical(scored_df$dise.brain, c(2L, 1L, 0L, 2L)))
stopifnot(identical(scored_df$pathway.brain, c(2L, 1L, 0L, 2L)))
stopifnot(identical(scored_df$tfbs, c(1L, 0L, 0L, 1L)))
stopifnot(identical(scored_df$bio.brain, c(2L, 1L, 0L, 2L)))
stopifnot(identical(scored_df$bio.inter, c(2L, 1L, 0L, 2L)))
cat("  -> Deterministic scoring matched exact expected benchmark scores.\n")

# Also test through load_builtin_annotations dual-mode interface
scored_via_loader <- load_builtin_annotations("entrez", version = "raw_scored", raw_df = raw_toy)
stopifnot(identical(scored_df, scored_via_loader))
cat("  -> Dual-mode loader interface (version = 'raw_scored') verified.\n\n")

# 4. Test Custom Annotation Preparation & Alignment
cat("[4/6] Testing custom user annotation alignment pipeline...\n")
custom_user_annot <- data.frame(
  GENE = c("1001", "1002", "1005"),
  cortex_scRNA = c(4.5, 0.2, 3.1),
  epigenetic_enhancer = c(1, 0, 1)
)
ready_annot <- prepare_custom_annotations(
  user_annotations = custom_user_annot,
  id_col = "GENE",
  base_annotations = annot_entrez[1:10, 1:4],
  base_id_col = "ID",
  fill_na = 0
)
cat(sprintf("  -> Merged & prepared custom annotation matrix: %d rows x %d cols\n\n",
            nrow(ready_annot), ncol(ready_annot)))

# 5. Test Simulation Framework
cat("[5/6] Generating simulated multi-tier genomic dataset...\n")
sim <- simulate_example_data(
  n_snps = 2000,
  n_genes = 300,
  n_aux = 2,
  n_annot = 5,
  prop_signal = 0.03,
  seed = 42
)
cat(sprintf("  -> Simulated %d SNPs, %d Genes, %d Auxiliary traits\n\n",
            nrow(sim$snp_target), nrow(sim$magma_target), length(sim$snp_aux)))

# 6. Test Multi-Tier FDR Regression Pipelines
cat("[6/6] Executing FDRreg across all three genomic tiers...\n")

# Tier 1: SNP Level (with LDSC Sample-Overlap Decorrelation + LASSO)
cat("  6.1 Testing Tier 1: SNP Level (LDSC decorrelation + LASSO)...\n")
res_snp <- run_fdrreg_snp(
  target          = sim$snp_target,
  aux             = sim$snp_aux,
  overlap_traits  = sim$overlap_traits,
  ldsc_intercepts = sim$ldsc_intercepts,
  annotations     = sim$annotations,
  feature_transform = "signed",
  var_select      = "lasso",
  fdrreg_nulltype = "both",
  seed            = 42
)
sig_snp <- significant(res_snp, threshold = 0.05, type = "theoretical")
cat(sprintf("      -> Model fit status: %s (Theoretical discoveries: %d)\n",
            res_snp$fit_status$theoretical, nrow(sig_snp)))

# Tier 2: MAGMA Gene Level (with Built-in Curated Annotations + LASSO)
cat("  6.2 Testing Tier 2: MAGMA Gene Level (Curated bio annotations)...\n")
# Ensure simulated gene IDs match across target and aux
real_gene_ids <- as.character(annot_entrez$ID[1:nrow(sim$magma_target)])
sim$magma_target$GENE <- real_gene_ids
for (nm in names(sim$magma_aux)) {
  sim$magma_aux[[nm]]$GENE <- real_gene_ids
}
res_magma <- run_fdrreg_magma_gene(
  target          = sim$magma_target,
  aux             = sim$magma_aux,
  annotations     = annot_entrez,
  id_col          = "GENE",
  feature_transform = "abs",
  var_select      = "lasso",
  fdrreg_nulltype = "both",
  seed            = 42
)
sig_magma <- significant(res_magma, threshold = 0.05, type = "theoretical")
cat(sprintf("      -> Model fit status: %s (Theoretical discoveries: %d)\n",
            res_magma$fit_status$theoretical, nrow(sig_magma)))

# Tier 3: TWAS Level (S-PrediXcan + Ensembl Curated Annotations)
cat("  6.3 Testing Tier 3: TWAS Level (S-PrediXcan)...\n")
real_ens_ids <- as.character(annot_ensembl$ENSEMBL_GENE_ID[1:nrow(sim$spredixcan_target)])
sim$spredixcan_target$gene <- real_ens_ids
for (nm in names(sim$spredixcan_aux)) {
  sim$spredixcan_aux[[nm]]$gene <- real_ens_ids
}
res_twas <- run_fdrreg_twas(
  target          = sim$spredixcan_target,
  twas_type       = "spredixcan",
  aux             = sim$spredixcan_aux,
  annotations     = annot_ensembl,
  id_col          = "gene",
  feature_transform = "abs",
  var_select      = "none",
  seed            = 42
)
sig_twas <- significant(res_twas, threshold = 0.05, type = "theoretical")
cat(sprintf("      -> Model fit status: %s (Theoretical discoveries: %d)\n",
            res_twas$fit_status$theoretical, nrow(sig_twas)))

# Evaluation Utility Check
cat("  6.4 Testing FDR Evaluation Utilities...\n")
fdr_eval <- evaluate_fdr_performance(
  fdr_values   = res_snp$full_results$fdr_theo,
  true_signals = sim$true_info$snp$is_signal
)
stopifnot(!is.null(fdr_eval))
cat(sprintf("      -> FDR Evaluation passed (True Discoveries count: %d)\n", fdr_eval$true_discoveries[1]))

cat("\n================================================================\n")
cat("All tests PASSED successfully! fdrregGenomics v0.3.0 is verified.\n")
cat("================================================================\n")
