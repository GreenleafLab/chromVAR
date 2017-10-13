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
#' @importFrom plotly ggplotly plotlyOutput renderPlotly
#' @importFrom graphics text
#' @importFrom stats anova approx as.dist cor cov dist dnorm hclust
#' kruskal.test lm median runif sd t.test var wilcox.test
#' na.omit oneway.test p.adjust pchisq pnorm prcomp quantile 
#' @importFrom utils installed.packages read.delim setTxtProgressBar 
#' txtProgressBar
#' @importFrom Rcpp sourceCpp
#' @importFrom S4Vectors queryHits subjectHits DataFrame elementNROWS isSorted
#' @importFrom GenomeInfoDb seqlevels seqlevels<- sortSeqlevels seqnames
#' @importFrom Rsamtools ScanBamParam scanBamFlag scanBam countBam
#' @importFrom IRanges IRanges overlapsAny
#' @importFrom BiocParallel bplapply
#' @importFrom TFBSTools PWMatrixList bg name
#' @importFrom BSgenome getBSgenome
#' @importMethodsFrom GenomicRanges sort start end
#' @importClassesFrom SummarizedExperiment SummarizedExperiment
#' @useDynLib chromVAR, .registration = TRUE
#' @docType package
#' @name chromVAR
NULL
# > NULL

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
#' @seealso \code{\link{mini_dev}}, \code{\link{mini_ix}}
#' @examples 
#' data(mini_counts)
NULL

#' mini_ix
#' 
#' Tiny sample annotation object for use in chromVAR examples
#' Result from running matchMotifs(example_motifs,mini_counts,"hg19) on 
#' example_motifs from motifmatchr package and mini_counts from this package
#' @docType data
#' @keywords datasets
#' @name mini_ix
#' @usage data(mini_ix)
#' @return \code{\link[SummarizedExperiment]{RangedSummarizedExperiment}}
#' @seealso \code{\link{mini_counts}}, \code{\link{mini_dev}}
#' @examples 
#' data(mini_ix)
NULL


#' mini_dev
#' 
#' Tiny sample chromVARDeviations object resulting from computeDeviations
#' Result from running computeDeviations(mini_counts, mini_ix) on 
#' mini_ix and mini_counts data from this package
#' @docType data
#' @keywords datasets
#' @name mini_dev
#' @usage data(mini_dev)
#' @return \code{\link{chromVARDeviations-class}}
#' @seealso \code{\link{computeDeviations}}, \code{\link{mini_counts}},
#' \code{\link{mini_ix}}
#' @examples 
#' data(mini_dev)
NULL




.onUnload <- function(libpath) {
  library.dynam.unload("chromVAR", libpath)
}




