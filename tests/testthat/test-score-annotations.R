test_that("score_david_annotations computes exact deterministic scores", {
  raw_df <- data.frame(
    ID = c("1", "2", "3", "4"),
    tissue = c("Fetal brain, cerebellum", "Liver", NA, "Hippocampus, synapse"),
    disease = c("Schizophrenia, Bipolar", "Hypertension", NA, "Drug abuse"),
    pathway = c("Dopaminergic synapse", "Glycolysis", NA, "GABAergic synapse"),
    tfbs_count = c(5, 1, NA, "siteA, siteB, siteC"),
    process = c("Synaptic transmission", "Cell cycle", NA, "Circadian regulation"),
    interaction = c("Dopamine receptor binding", "Insulin receptor", NA, "Serotonin transporter")
  )
  
  scores <- score_david_annotations(raw_df, id_col = "ID")
  
  expect_equal(scores$expre.brain, c(1, 0, 0, 1))
  expect_equal(scores$dise.brain, c(2, 1, 0, 2))
  expect_equal(scores$pathway.brain, c(2, 1, 0, 2))
  expect_equal(scores$tfbs, c(1, 0, 0, 1))
  expect_equal(scores$bio.brain, c(2, 1, 0, 2))
  expect_equal(scores$bio.inter, c(2, 1, 0, 2))
})

test_that("load_builtin_annotations works with both version options", {
  # Option 1: curated benchmark
  curated <- load_builtin_annotations("entrez", version = "curated")
  expect_true(is.data.frame(curated))
  expect_equal(nrow(curated), 19503)
  expect_true("dise.brain" %in% colnames(curated))
  
  # Option 2: raw_scored
  raw_df <- data.frame(
    ID = c("100", "200"),
    tissue = c("brain", "lung"),
    disease = c("autism", "cancer"),
    pathway = c("neuro", "metabolism"),
    tfbs_count = c(2, 0),
    process = c("brain dev", "mitosis"),
    interaction = c("synapse", "tubulin")
  )
  fresh <- load_builtin_annotations("entrez", version = "raw_scored", raw_df = raw_df)
  expect_equal(nrow(fresh), 2)
  expect_equal(fresh$expre.brain, c(1, 0))
  expect_equal(fresh$dise.brain, c(2, 1))
})
