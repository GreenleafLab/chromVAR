#' chromVAR: A package for computing variability across sets of peaks.
#'
#' Determine variation in chromatin accessibility across sets of
#' annotations or peaks. Designed primarily for single-cell or sparse chromatin
#' accessibility, e.g. from scATAC-seq or sparse ATAC or DNAse-seq experiments.
#'
#' @import Matrix
#' @import GenomicRanges
#' @import Biostrings
#' @import methods
#' @importFrom Rcpp sourceCpp
#' @useDynLib chromVAR
#' @docType package
#' @name chromVAR
NULL
#> NULL


.onUnload <- function (libpath) {
  library.dynam.unload("mypackage", libpath)
}
