#' @keywords internal
#' @aliases fdrregGenomics-package
#'
#' @title fdrregGenomics: Genomics-Integrated FDR Regression Pipeline
#'
#' @description A comprehensive pipeline for integrating GWAS summary statistics
#'   with biological annotations using False Discovery Rate Regression (FDRreg).
#'
#' @section Disclaimer / Attribution:
#' The underlying FDRreg algorithm is developed by Scott et al. (2016);
#' this package provides the comprehensive genomics integration,
#' preprocessing, and validation workflow.
#'
#' @references
#' Scott, J. G., Kelly, R. C., Smith, M. A., Zhou, P., & Kass, R. E. (2016).
#' False discovery rate regression: an application to neural synchrony
#' detection in primary visual cortex.
#' \emph{Journal of the American Statistical Association}, 110(510), 459-471.
#'
#' @importFrom FDRreg FDRreg
#' @importFrom powerplus Matpow
#' @importFrom glmnet cv.glmnet
#' @importFrom HelpersMG SEfromHessian
#' @importFrom stats qnorm pnorm rnorm cor sd coef runif rbinom
#' @importFrom utils packageVersion head
"_PACKAGE"