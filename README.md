# fdrregGenomics

A comprehensive R package for integrating GWAS summary statistics with biological annotations using False Discovery Rate Regression (FDRreg).

## Installation

```r
# Install the development version
devtools::install_github("keelost/fdrregGenomics")

# Or install from source code
devtools::install("path/to/fdrregGenomics")
```

## Quick Start

```r
library(fdrregGenomics)

# Load sample data
sim <- simulate_example_data(n_snps = 1000, n_genes = 100)

# SNP level analysis
result <- run_fdrreg_snp(
target = sim$snp_target,
aux = sim$snp_aux,
feature_transform = "signed",
var_select = "lasso"
)

# View results
print(result)
significant(result, threshold = 0.05)
```

## Features

- Supports three analysis levels: SNP, MAGMA gene, and TWAS (S-PrediXcan/S-MultiXcan)
- Optional sample overlap removal
- Biological annotation integration
- Variable selection methods: LASSO, elasticity network, marginal screening
- Complete results output and visualization

## Documentation

- [Quick Start](vignettes/01-quickstart.Rmd)
- [SNP Level Analysis Tutorial](vignettes/02-snp-analysis.Rmd)
- [Real Data Analysis Workflow](vignettes/03-real-data-workflow.Rmd)

## License

MIT License - see [LICENSE](LICENSE) file
