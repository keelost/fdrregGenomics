# fdrregGenomics

Genomics-Integrated False Discovery Rate Regression Pipeline

## Disclaimer / Attribution

> **The underlying FDRreg algorithm is developed by Scott et al. (2016,
> JASA); this package provides the comprehensive genomics integration,
> preprocessing, and validation workflow.**

## Overview

`fdrregGenomics` provides a unified pipeline for integrating GWAS summary
statistics with biological annotations using False Discovery Rate
Regression (FDRreg). It supports three genomic analysis levels:

| Tier | Function | Description |
|------|----------|-------------|
| SNP | `run_fdrreg_snp()` | SNP-level GWAS with LDSC sample-overlap decorrelation |
| MAGMA Gene | `run_fdrreg_magma_gene()` | MAGMA gene-level results |
| TWAS | `run_fdrreg_twas()` | S-PrediXcan and S-MultiXcan results |

### Features

- **Sample-overlap decorrelation** via LDSC intercepts and matrix power
  transformation (SNP level)
- **Variable selection**: LASSO, Marginal Screening, Elastic Net, or none
- **Biological annotation integration** with flexible ID column matching
- **Both null types**: theoretical and empirical FDRreg models
- **Standardized S3 output** with `print()`, `summary()`, and
  `significant()` methods
- **Reproducibility logs** including seed, R version, package versions,
  and timestamps
- **Enhanced simulation**: Multiple simulation modes (full, summary, raw),
  complex signal models, and mixture effect distributions
- **Evaluation utilities**: Functions for FDR control and variable selection
  performance assessment

## Installation

```r
# Install from GitHub
# install.packages("devtools")
devtools::install_github("keelost/fdrregGenomics")
```

## Quick Start

### Basic Example (Original Mode)

```r
library(fdrregGenomics)

# Simulate example data (full mode - default)
sim <- simulate_example_data(n_snps = 2000, n_genes = 200, seed = 42)

# SNP-level analysis with LDSC decorrelation
result <- run_fdrreg_snp(
target = sim$snp_target,
aux = sim$snp_aux,
overlap_traits = sim$overlap_traits,
ldsc_intercepts = sim$ldsc_intercepts,
annotations = sim$annotations,
feature_transform = "signed",
var_select = "lasso",
seed = 42
)

# View results
print(result)

# Extract significant findings
sig <- significant(result, threshold = 0.05, type = "theoretical")
nrow(sig)
```

### Advanced Example (New Simulation Modes)

```r
library(fdrregGenomics)

# Simulate summary statistics with complex signal model
sim_summary <- simulate_example_data(
n_snps = 2000,
simulation_mode = "summary_only",
signal_model = "complex",
signal_function = function(x1, x2, x3) {
 -3.0 + 0.8*x1 + 1.0*x2 + 1.2*x3
},
n_annot = 3,
seed = 42
)

# SNP analysis with new simulation mode
result_summary <- run_fdrreg_snp(
target = sim_summary$snp,
id_col = "id",
z_col = "z",
p_col = "pval",
feature_transform = "abs",
seed = 42
)

# Evaluate FDR control performance
# Use the true signals from simulation
fdr_perf <- evaluate_fdr_performance(
fdr_values = result_summary$full_results$fdr_theo,
true_signals = sim_summary$true_info$snp$is_signal
)
print(fdr_perf)
```

### Example with Variable Selection Evaluation

```r
library(fdrregGenomics)

# Simulate data with known signals
sim <- simulate_example_data(
n_snps = 1000,
simulation_mode = "summary_only",
signal_model = "complex",
signal_function = function(x1, x2, x3) -3 + 0.8*x1 + 1.0*x2 + 1.2*x3,
n_annot = 3,
seed = 42
)

# Run analysis with variable selection
result <- run_fdrreg_snp(
target = sim$snp,
annotations = sim$annotations,
id_col = "id",
z_col = "z",
p_col = "pval",
feature_transform = "abs",
var_select = "lasso",
seed = 42
)

# Evaluate variable selection if performed
if (!is.null(result$varselect_model)) {
# Get the true signals (convert to logical)
true_signals <- which(sim$true_info$snp$is_signal == 1)

# Get the selected features from the model
selected_features <- result$covariates_used

# Find column indices of selected features in annotations
# Note: annotations has 'id' as first column, so features start from column 2
annot_colnames <- colnames(sim$annotations)
selected_indices <- which(annot_colnames %in% selected_features)

# Get total number of annotation features (excluding ID column)
n_annot_features <- ncol(sim$annotations) - 1

# Evaluate variable selection
var_eval <- evaluate_variable_selection(
 selected_vars = selected_indices,
 true_vars = true_signals,
 total_vars = n_annot_features
)

cat("Variable Selection Performance:\n")
cat("Precision:", round(var_eval$precision, 3), "\n")
cat("Recall:", round(var_eval$recall, 3), "\n")
cat("F1 Score:", round(var_eval$F1, 3), "\n")
cat("Selected features:", length(selected_features), "\n")
cat("True signals:", sum(sim$true_info$snp$is_signal), "\n")
}
```

## Input Data Formats

### Target GWAS (SNP-level)
```
snpid,chr,bp,a1,a2,z,pval
rs12345,1,100000,A,G,2.31,0.021
```

### S-PrediXcan
```
gene,gene_name,zscore,pvalue
ENSG00000167550.6,RHEBL1,5.32,1.1e-07
```

### S-MultiXcan
```
gene,gene_name,pvalue
ENSG00000167550.6,RHEBL1,8.9e-08
```

### MAGMA Gene Output (`.genes.out`)
```
GENE CHR START STOP NSNPS NPARAM N ZSTAT P
1001 1 100000 200000 50 10 100000 4.5 6.7e-06
```

## Simulation Modes

The `simulate_example_data()` function supports three simulation modes:

| Mode | Description | Use Case |
|------|-------------|----------|
| `full` | Complete simulated datasets (default) | Testing full pipeline |
| `summary_only` | Z-scores and p-values only | Method comparison studies |
| `raw_only` | Raw data with features and response | Traditional statistical analysis |

### Complex Signal Models

```r
# Define a complex signal model
complex_model <- function(x1, x2, x3, x4, x5) {
-2.5 + 0.5*x1 + 0.8*x2 + 1.2*x3 - 0.3*x4 + 0.1*x5
}

# Generate simulated data with complex signals
sim_complex <- simulate_example_data(
n_snps = 2000,
simulation_mode = "summary_only",
signal_model = "complex",
signal_function = complex_model,
n_annot = 5,
effect_distribution = list(
 type = "mixture",
 weights = c(0.4, 0.2, 0.4),
 means = c(-1.25, 0, 1.25),
 sds = c(1, 0.8, 1)
),
seed = 42
)

# Run analysis
result_complex <- run_fdrreg_snp(
target = sim_complex$snp,
id_col = "id",
z_col = "z",
p_col = "pval",
feature_transform = "abs",
seed = 42
)

# Evaluate
fdr_perf <- evaluate_fdr_performance(
fdr_values = result_complex$full_results$fdr_theo,
true_signals = sim_complex$true_info$snp$is_signal
)
print(fdr_perf)
```

## Evaluation Functions

```r
# Evaluate FDR control performance
fdr_eval <- evaluate_fdr_performance(
fdr_values = result$full_results$fdr_theo,
true_signals = sim$true_info$snp$is_signal, # Use actual true signals
thresholds = c(0.05, 0.1, 0.2)
)

# Evaluate variable selection performance
# Note: Adjust selected_vars and true_vars based on your data
var_eval <- evaluate_variable_selection(
selected_vars = selected_indices, # Actual selected indices
true_vars = true_signals, # Actual true signals
total_vars = n_annot_features # Total number of features
)
```

## Citation

If you use this package, please cite:

- Scott, J. G., Kelly, R. C., Smith, M. A., Zhou, P., & Kass, R. E.
  (2016). False discovery rate regression: an application to neural
  synchrony detection in primary visual cortex. *Journal of the American
  Statistical Association*, 110(510), 459-471.

## License

MIT License - see [LICENSE](LICENSE) file