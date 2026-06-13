#' Unified FDRreg Pipeline Dispatcher
#'
#' Routes to the appropriate tier-specific function.
#'
#' @param tier Character, "snp", "magma_gene", or "twas".
#' @param ... Arguments passed to the tier-specific function.
#' @return An object of class \code{fdrreg_result}.
#' @export
run_fdrreg <- function(tier = c("snp", "magma_gene", "twas"), ...) {
  tier <- match.arg(tier)
  switch(tier,
    snp        = run_fdrreg_snp(...),
    magma_gene = run_fdrreg_magma_gene(...),
    twas       = run_fdrreg_twas(...)
  )
}
