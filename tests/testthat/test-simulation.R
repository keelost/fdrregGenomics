library(testthat)
library(fdrregGenomics)

test_that("simulate_example_data works with all modes", {
  # Define test parameters
  n_snps_test <- 100
  n_genes_test <- 20
  
  # Test 1: Full mode (default)
  sim_full <- simulate_example_data(
    n_snps = n_snps_test, 
    n_genes = n_genes_test, 
    seed = 42
  )
  expect_true(is.list(sim_full))
  expect_true("snp_target" %in% names(sim_full))
  expect_true("annotations" %in% names(sim_full))
  expect_true(nrow(sim_full$snp_target) == n_snps_test)
  expect_true(ncol(sim_full$annotations) == 11)  # snpid + 10 annotations
  
  # Test 2: Summary only mode
  sim_summary <- simulate_example_data(
    n_snps = n_snps_test,
    n_genes = n_genes_test,
    simulation_mode = "summary_only",
    seed = 42
  )
  expect_true(is.list(sim_summary))
  expect_true("snp" %in% names(sim_summary))
  expect_true("annotations" %in% names(sim_summary))
  # In summary_only mode, nrow(snp) should be n_snps (not n_snps + n_genes)
  expect_true(nrow(sim_summary$snp) == n_snps_test)
  expect_true("id" %in% colnames(sim_summary$snp))
  expect_true("id" %in% colnames(sim_summary$annotations))
  
  # Test 3: Raw only mode
  sim_raw <- simulate_example_data(
    n_snps = n_snps_test,
    n_genes = n_genes_test,
    simulation_mode = "raw_only",
    seed = 42
  )
  expect_true(is.list(sim_raw))
  expect_true("raw" %in% names(sim_raw))
  # In raw_only mode, raw data contains both SNPs and genes
  expect_true(nrow(sim_raw$raw) == n_snps_test + n_genes_test)
  expect_true("y" %in% colnames(sim_raw$raw))
})

test_that("complex signal models work", {
  n_snps_test <- 100
  
  # Test complex signal model
  sim_complex <- simulate_example_data(
    n_snps = n_snps_test,
    simulation_mode = "summary_only",
    signal_model = "complex",
    signal_function = function(x1, x2, x3) -3 + 0.8*x1 + 1.0*x2 + 1.2*x3,
    n_annot = 3,
    seed = 42
  )
  
  expect_true(is.list(sim_complex))
  expect_true(nrow(sim_complex$snp) == n_snps_test)
  expect_true(ncol(sim_complex$annotations) == 4)  # id + 3 annotations
  expect_true(mean(sim_complex$true_info$snp$is_signal) > 0)
})

test_that("evaluate functions work", {
  
  # Create test data
  n_test <- 100
  fdr_values <- runif(n_test, 0, 0.3)
  true_signals <- rbinom(n_test, 1, 0.1) == 1
  
  # Test FDR evaluation
  fdr_eval <- evaluate_fdr_performance(fdr_values, true_signals)
  expect_true(is.data.frame(fdr_eval))
  expect_true("threshold" %in% colnames(fdr_eval))
  expect_true("FDP" %in% colnames(fdr_eval))
  expect_true("power" %in% colnames(fdr_eval))
  
  # Test variable selection evaluation
  selected_vars <- sample(1:n_test, 15)
  true_vars <- sample(1:n_test, 20)
  
  var_eval <- evaluate_variable_selection(selected_vars, true_vars, total_vars = n_test)
  expect_true(is.list(var_eval))
  expect_true("precision" %in% names(var_eval))
  expect_true("recall" %in% names(var_eval))
  expect_true("F1" %in% names(var_eval))
})
