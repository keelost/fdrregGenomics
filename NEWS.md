# NEWS / Changelog for fdrregGenomics

## fdrregGenomics 0.3.0 (Major Release)

### Major Features & Data Additions
* **Embedded Psychiatric Biological Annotations**: Added full, frozen gene-level annotation dataset comprising 23 biological features across 19,503 human genes (`psychiatric_annotations.rda`, Supplementary Table S1.21).
* **Automated Deterministic Annotation Scoring (`score_david_annotations()`)**: Implemented deterministic substring keyword matching engine translating raw functional database queries (DAVID 6.8, OMIM, Reactome) into the exact 6 tiered and binary numerical features specified in Supplementary Methods Section 1.
* **Built-in Annotation Loader (`load_builtin_annotations()`)**: Seamless extraction of annotation matrices keyed by Entrez ID (`"entrez"`), Ensembl ID (`"ensembl"`), or HGNC Gene Symbol (`"symbol"`), with support for both frozen benchmark (`version = "curated"`) and freshly calculated (`version = "raw_scored"`) matrices.
* **Extensible Custom Annotation Pipeline (`prepare_custom_annotations()`)**: Added general-purpose utility to format, validate, zero-fill, and merge arbitrary user-supplied annotation matrices (e.g. single-cell markers, custom pathway scores, or ChIP-seq peaks) with automatic zero-variance column pruning.
* **Command-Line Interface (CLI)**: Added executable CLI utility `inst/scripts/score_david_annotations.r` for reproducible batch scoring from the terminal.

### Documentation & Vignettes
* **Comprehensive Data Dictionary**: Fully documented all 23 annotation columns in `?psychiatric_annotations` detailing the 6 DAVID 6.8 functional tiers, Open Targets Platform scores, and denovo-db mutations.
* **Expanded Getting Started Vignette**: Added dedicated section demonstrating end-to-end usage with built-in and user-customized annotations on independent GWAS datasets.
* **Updated README**: Documented new built-in annotation access and extensibility guidelines.

---

## fdrregGenomics 0.2.0

* **Realistic Simulation Framework**: Added `simulate_example_data()` generating synthetic multi-trait GWAS summary statistics (SNP level with LDSC overlap, MAGMA gene statistics, and S-PrediXcan/S-MultiXcan TWAS data) for benchmarking and tutorials.
* **Simulation Modes**: Supported `full`, `summary_only`, and `raw_only` simulation schemes with complex nonlinear signal functions and mixture effect distributions.
* **Evaluation Utilities**: Added `evaluate_fdr_performance()` and `evaluate_variable_selection()` for benchmarking error control and variable selection metrics (Precision, Recall, F1).
* **Enhanced Diagnostic Outputs**: Added S3 methods (`print`, `summary`, `significant`) and standardized discovery reporting.

---

## fdrregGenomics 0.1.0

* Initial release of the unified genomics FDR regression pipeline supporting SNP, MAGMA gene, and TWAS analysis tiers.
