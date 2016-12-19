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
#' @importFrom graphics text
#' @importFrom stats anova approx as.dist cor cov dist dnorm hclust kruskal.test lm median
#' @importFrom stats na.omit oneway.test p.adjust pchisq pnorm prcomp quantile runif sd t.test var wilcox.test
#' @importFrom utils installed.packages read.delim setTxtProgressBar txtProgressBar
#' @importFrom Rcpp sourceCpp
#' @useDynLib chromVAR
#' @docType package
#' @name chromVAR
NULL
#> NULL


.onUnload <- function (libpath) {
  library.dynam.unload("chromVAR", libpath)
}
