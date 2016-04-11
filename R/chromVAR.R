#' chromVAR: A package for computing variability across sets of peaks.
#'
#' The package provides three categories of important functions:
#' functions!
#'
#' @section Functions:
#' The functions ...
#'
#' @section Classes:
#' the classes ...
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
