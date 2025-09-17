#' FS4Clustering: Feature screening for clustering analysis
#'
#' Feature screening for clustering analysis of count data with an application to single-cell RNA-sequencing. See Wang, Chen and Xi (2023) <https://doi.org/10.48550/arXiv.2306.12671> for more details.

#' @name FS4Clustering
#' @docType package
#' @useDynLib FS4Clustering
## usethis namespace: start
#' @import RcppEigen
#' @import dplyr
#' @importFrom nloptr nl.grad
#' @importFrom nloptr nloptr
#' @importFrom Rcpp evalCpp
#' @importFrom stats dnbinom
#' @importFrom stats lm
#' @importFrom stats p.adjust
#' @importFrom stats pchisq
#' @importFrom stats dpois
#' @importFrom stats quantile
#' @importFrom stats rpois
#' @importFrom stats var
#' @importFrom grDevices colorRampPalette
#' @importFrom rlang .data
## usethis namespace: end
NULL
