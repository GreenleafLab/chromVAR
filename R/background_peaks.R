#' addGCBias
#'
#' Computes GC content for peaks
#' @param object (Ranged)SummarizedExperiment
#' @param ... additional arguments
#' @return (Ranged)SummarizedExperiment object with new column in row metadata
#' with the gc content of the peak in question
#' @export
#' @examples 
#' 
#' data(example_counts, package = "chromVAR")
#' # show example on small part of data 
#' subset_counts <- example_counts[1:500,]
#' library(BSgenome.Hsapiens.UCSC.hg19)
#' example_counts <- addGCBias(subset_counts, 
#'                               genome = BSgenome.Hsapiens.UCSC.hg19)
#' 
setGeneric("addGCBias", function(object, ...) standardGeneric("addGCBias"))

#' @describeIn addGCBias method for RangedSummarizedExperiment
#' @param genome BSgenome object, by defualt hg19
#' @export
setMethod(addGCBias, c(object = "RangedSummarizedExperiment"), 
          function(object, 
                   genome = GenomeInfoDb::genome(object)) {
            genome <- validate_genome_input(genome)
            peaks <- rowRanges(object)
            seqs <- getSeq(genome, peaks)
            nucfreqs <- letterFrequency(seqs, c("A", "C", "G", "T"))
            gc <- rowSums(nucfreqs[, 2:3]) / rowSums(nucfreqs)
            rowRanges(object)$bias <- gc
            return(object)
          })

#' @param peaks GenomicRanges with peaks, needed if object is 
#' SummarizedExperiment and not RangedSummarizedExperiment
#' @describeIn addGCBias method for SummarizedExperiment
#' @export
setMethod(addGCBias, c(object = "SummarizedExperiment"), 
          function(object, peaks, 
                   genome = GenomeInfoDb::genome(peaks)) {
            genome <- validate_genome_input(genome)
            seqs <- getSeq(genome, peaks)
            nucfreqs <- letterFrequency(seqs, c("A", "C", "G", "T"))
            gc <- rowSums(nucfreqs[, 2:3]) / rowSums(nucfreqs)
            rowData(object)$bias <- gc
            return(object)
          })

#' getBackgroundPeaks
#'
#' Function to get a set of background peaks for each peak based on GC content 
#' and # of fragments across all samples
#' @param object fragment counts as SummarizedExperiment, RangedSummarized, 
#' Matrix, or matrix
#' @param bias vector of values for some bias signal for each row of object
#' @param niterations number of background peaks to sample
#' @param w parameter controlling similarity of background peaks
#' @param bs bin size parameter
#' @param ... addtional arguments
#' @return matrix with one row per peak and one column per iteration.  values in
#'  a row represent indices of background peaks for the peak with that index
#' @details Background peaks are chosen by sampling peaks based on similarity in
#' GC content and # of fragments across samples using the Mahalanobis distance.
#' The w paramter controls how similar background peaks should be. The bs 
#' parameter controls the precision with which the similarity is computed;
#' increasing bs will make the function run slower. Sensible default parameters 
#' are chosen for both.
#' @export
#' @examples
#' 
#' # Load very small example counts (already filtered)
#' data(mini_counts, package = "chromVAR")
#' 
#' # get background peaks
#' bgpeaks <- getBackgroundPeaks(mini_counts)
#' 
setGeneric("getBackgroundPeaks", 
           function(object, ...) standardGeneric("getBackgroundPeaks"))

#' @describeIn getBackgroundPeaks method for SummarizedExperiment
#' @export
setMethod(getBackgroundPeaks, c(object = "SummarizedExperiment"), 
          function(object, 
                   bias = rowData(object)$bias, 
                   niterations = 50, 
                   w = 0.1, 
                   bs = 50) {
            object <- counts_check(object)
            get_background_peaks_core(object, bias, niterations, w, bs)
          })

#' @describeIn getBackgroundPeaks method for RangedSummarizedExperiment
#' @export
setMethod(getBackgroundPeaks, c(object = "RangedSummarizedExperiment"), 
          function(object, 
                   bias = rowRanges(object)$bias, 
                   niterations = 50, 
                   w = 0.1, 
                   bs = 50) {
            object <- counts_check(object)
            get_background_peaks_core(object, bias, niterations, w, bs)
          })

#' @describeIn getBackgroundPeaks method for Matrix or matrix
#' @export
setMethod(getBackgroundPeaks, c(object = "MatrixOrmatrix"), 
          function(object, 
                   bias, 
                   niterations = 50, 
                   w = 0.1, 
                   bs = 50) {
            get_background_peaks_core(object, bias, niterations, w, bs)
          })

get_background_peaks_core <- function(object, 
                                      bias, 
                                      niterations = 50, 
                                      w = 0.1, 
                                      bs = 50) {
  
  fragments_per_peak <- getFragmentsPerPeak(object)
  stopifnot(length(bias) == length(fragments_per_peak))
  if (min(fragments_per_peak) <= 0) 
    stop("All peaks must have at least one fragment in one sample")
  
  intensity <- log10(fragments_per_peak)
  norm_mat <- matrix(c(intensity, bias), ncol = 2, byrow = FALSE)
  
  chol_cov_mat <- chol(cov(norm_mat))
  trans_norm_mat <- t(forwardsolve(t(chol_cov_mat), t(norm_mat)))
  
  # make bins
  bins1 <- seq(min(trans_norm_mat[, 1]), max(trans_norm_mat[, 1]), 
               length.out = bs)
  bins2 <- seq(min(trans_norm_mat[, 2]), max(trans_norm_mat[, 2]), 
               length.out = bs)
  
  bin_data <- do.call(rbind, lapply(seq_len(bs), 
                                    function(x) matrix(c(rep(bins1[x], bs), 
                                                         bins2), ncol = 2, 
                                                       byrow = FALSE)))
  
  bin_dist <- euc_dist(bin_data)
  bin_p <- dnorm(bin_dist, 0, w)
  
  bin_membership <- nabor::knn(bin_data, query = trans_norm_mat, k = 1)$nn.idx
  
  bin_density <- tabulate2(bin_membership, min_val = 1, max_val = bs^2)
  
  background_peaks <- bg_sample_helper(bin_membership - 1, bin_p, bin_density, 
                                       niterations)
  
  return(background_peaks)
  
}

#' getPermutedData
#'
#' Function to get permuted data while maintaining biases
#' @param object SummarizedExperiment
#' @param niterations number of background peaks to sample
#' @param w parameter controlling similarity of background peaks
#' @param bs bin size parameter
#' @return new SummarizedExperiment with addition assays representing permuted
#' version of counts
#' @details Replaces the counts at a given peak with the count from another peak
#' with similar GC content and average accessibility
#' @export
#' @examples 
#' 
#' # Load very small example counts (already filtered)
#' data(mini_counts, package = "chromVAR")
#' 
#' # get background peaks
#' perm_counts <- getPermutedData(mini_counts, niterations = 2)
#' 
getPermutedData <- function(object, niterations = 10, w = 0.1, bs = 50) {
  
  out <- bplapply(seq_len(niterations), 
                                function(x) get_permuted_data_helper(object, 
                                                                     w = w, 
                                                                     bs = bs))
  names(out) <- paste("perm", seq_len(niterations), sep = "_")
  
  assays(object) <- c(assays(object), out)
  
  return(object)
}


get_permuted_data_helper <- function(object, w, bs) {
  bgpeaks <- getBackgroundPeaks(object, niterations = ncol(object), 
                                  w = w, bs = bs)
  reorder_columns(assays(object)$counts, bgpeaks)
}

reorder_columns <- function(mat, colixmat) {
  reordered <- lapply(seq_len(ncol(mat)), 
                      function(x) 
                        mat[colixmat[, x], x, drop = FALSE])
  return(do.call(cbind, reordered))
}




