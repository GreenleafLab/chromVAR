#' chromVAR: A package for computing variability across sets of peaks.
#'
#' Determine variation in chromatin accessibility across sets of
#' annotations or peaks. Designed primarily for single-cell or sparse chromatin
#' accessibility, e.g. from scATAC-seq or sparse ATAC or DNAse-seq experiments.
#'
#' @import Matrix
#' @import methods
#' @import GenomicRanges
#' @import SummarizedExperiment
#' @import Biostrings
#' @import ggplot2
#' @import shiny
#' @import miniUI
#' @import BSgenome.Hsapiens.UCSC.hg19 
#' @importFrom plotly ggplotly plotlyOutput renderPlotly
#' @importFrom graphics text
#' @importFrom stats anova approx as.dist cor cov dist dnorm hclust kruskal.test lm median
#' @importFrom stats na.omit oneway.test p.adjust pchisq pnorm prcomp quantile runif sd t.test var wilcox.test
#' @importFrom utils installed.packages read.delim setTxtProgressBar txtProgressBar
#' @importFrom Rcpp sourceCpp
#' @importFrom S4Vectors queryHits subjectHits DataFrame elementLengths isSorted
#' @importFrom GenomeInfoDb seqlevels seqlevels<- sortSeqlevels seqnames
#' @importFrom BiocParallel bplapply
#' @importFrom Rsamtools ScanBamParam scanBamFlag scanBam countBam
#' @importFrom IRanges IRanges overlapsAny
#' @importFrom BiocParallel bplapply
#' @importFrom TFBSTools PWMatrixList bg name
#' @importMethodsFrom GenomicRanges sort start end
#' @importClassesFrom SummarizedExperiment SummarizedExperiment
#' @useDynLib chromVAR
#' @docType package
#' @name chromVAR
NULL
# > NULL

#' human_pwms_v1
#'
#' Collection of human pwms
#' @docType data
#' @keywords datasets
#' @name human_pwms_v1
#' @usage data(human_pwms_v1)
#' @return \code{\link[TFBSTools]{XMatrixList}} of length 1764
#' @examples 
#' data(human_pwms_v1)
NULL

#' mouse_pwms_v1
#'
#' Collection of mouse pwms
#' @docType data
#' @keywords datasets
#' @name mouse_pwms_v1
#' @usage data(mouse_pwms_v1)
#' @return \code{\link[TFBSTools]{XMatrixList}} of length 1346
#' @examples 
#' data(mouse_pwms_v1)
NULL

#' example_counts
#' 
#' Very small sample data set for trying out chromVAR
#' @docType data
#' @keywords datasets
#' @name example_counts
#' @usage data(example_counts)
#' @return \code{\link[SummarizedExperiment]{RangedSummarizedExperiment}}
#' @examples 
#' data(example_counts)
NULL

#' mini_counts
#' 
#' Tiny sample data set for chromVAR funtion examples
#' @docType data
#' @keywords datasets
#' @name mini_counts
#' @usage data(mini_counts)
#' @return \code{\link[SummarizedExperiment]{RangedSummarizedExperiment}}
#' @examples 
#' data(mini_counts)
NULL

.onUnload <- function(libpath) {
  library.dynam.unload("chromVAR", libpath)
}
