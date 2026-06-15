test_that("simulate_example_data works correctly", {
  # 测试默认参数
  sim <- simulate_example_data()
  
  expect_true(is.list(sim))
  expect_true(all(c("snp_target", "magma_target", "spredixcan_target") %in% names(sim)))
  
  # 检查数据结构
  expect_true(nrow(sim$snp_target) > 0)
  expect_true(all(c("snpid", "z", "pval") %in% names(sim$snp_target)))
  
  # 测试自定义参数
  sim_custom <- simulate_example_data(n_snps = 500, n_genes = 50, seed = 123)
  expect_equal(nrow(sim_custom$snp_target), 500)
})

test_that("run_fdrreg_snp works with simulated data", {
  sim <- simulate_example_data(n_snps = 1000, n_genes = 100)
  
  # 基础测试
  result <- run_fdrreg_snp(
    target = sim$snp_target,
    aux = sim$snp_aux,
    feature_transform = "signed",
    seed = 42
  )
  
  expect_s3_class(result, "fdrreg_result")
  expect_true("full_results" %in% names(result))
  expect_true("fdr_theo" %in% names(result$full_results))
  
  # 测试显著性函数
  sig <- significant(result, threshold = 0.05)
  expect_true(is.data.frame(sig))
})

test_that("run_fdrreg_magma_gene works with simulated data", {
  sim <- simulate_example_data(n_snps = 1000, n_genes = 100)
  
  result <- run_fdrreg_magma_gene(
    target = sim$magma_target,
    aux = sim$magma_aux,
    annotations = sim$gene_annotations,
    id_col = "GENE",
    annot_id_col = "gene",
    seed = 42
  )
  
  expect_s3_class(result, "fdrreg_result")
  expect_true("tier" %in% names(result))
  expect_true(result$tier == "magma_gene")
})
