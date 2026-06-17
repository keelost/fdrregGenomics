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
- **Evaluation utilities**: Functions to assess FDR control and variable selection
- **Intelligent parameter handling**: Automatic adjustment of signal function arguments

## Installation

```r
# Install from GitHub
# install.packages("devtools")
devtools::install_github("keelost/fdrregGenomics")
```

## Quick Start

```r
library(fdrregGenomics)

# Simulate summary statistics with complex signal model
sim_complex <- simulate_example_data(
  n_snps = 2000,
  simulation_mode = "summary_only",
  signal_model = "complex",
  signal_function = function(x1, x2, x3) {
    -3.0 + 0.8*x1 + 1.0*x2 + 1.2*x3
  },
  effect_distribution = list(
    type = "mixture",
    weights = c(0.4, 0.2, 0.4),
    means = c(-1.25, 0, 1.25),
    sds = c(1, 0.8, 1)
  ),
  seed = 42
)


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

# Evaluate FDR performance
fdr_perf <- evaluate_fdr_performance(
  fdr_values = result$full_results$fdr_theo,
  true_signals = sim_complex$true_info$snp$is_signal
)
print(fdr_perf)
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

## Citation

If you use this package, please cite:

- Scott, J. G., Kelly, R. C., Smith, M. A., Zhou, P., & Kass, R. E.
  (2016). False discovery rate regression: an application to neural
  synchrony detection in primary visual cortex. *Journal of the American
  Statistical Association*, 110(510), 459-471.

## License

MIT License - see [LICENSE](LICENSE) file
