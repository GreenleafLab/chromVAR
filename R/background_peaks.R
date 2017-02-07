#' add_gc_bias
#'
#' Computes GC content for peaks
#' @param object SummarizedExperiment
#' @param ... additional arguments
#' @export
setGeneric("add_gc_bias", function(object, ...) standardGeneric("add_gc_bias"))

#' @describeIn add_gc_bias method for RangedSummarizedExperiment
#' @param genome BSgenome object, by defualt hg19
#' @export
setMethod(add_gc_bias, c(object = "RangedSummarizedExperiment"), 
          function(object, 
                   genome = BSgenome.Hsapiens.UCSC.hg19) {
            peaks <- rowRanges(object)
            seqs <- getSeq(genome, peaks)
            nucfreqs <- letterFrequency(seqs, c("A", "C", "G", "T"))
            gc <- apply(nucfreqs, 1, function(x) sum(x[2:3])/sum(x))
            rowRanges(object)$bias <- gc
            return(object)
          })

#' @param peaks GenomicRanges with peaks, needed if object is SummarizedExperiment
#' and not RangedSummarizedExperiment
#' @describeIn add_gc_bias method for SummarizedExperiment
#' @export
setMethod(add_gc_bias, c(object = "SummarizedExperiment"), 
          function(object, peaks, 
                   genome = BSgenome.Hsapiens.UCSC.hg19) {
            seqs <- getSeq(genome, peaks)
            nucfreqs <- letterFrequency(seqs, c("A", "C", "G", "T"))
            gc <- apply(nucfreqs, 1, function(x) sum(x[2:3])/sum(x))
            rowData(object)$bias <- gc
            return(object)
          })

#' get_background_peaks
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
#' @return matrix with one row per peak and one column per iteration.  values in a row
#' represent indices of background peaks for the peak with that index
#' @details Background peaks are chosen by sampling peaks based on similarity in
#' GC content and # of fragments across samples using the Mahalanobis distance.
#' The w paramter controls how similar background peaks should be. The bs parameter
#' controls the precision with which the similarity is computed; increasing bs will
#' make the function run slower. Sensible default parameters are chosen for both.
#' @export
setGeneric("get_background_peaks", 
           function(object, ...) standardGeneric("get_background_peaks"))

#' @describeIn get_background_peaks method for SummarizedExperiment
#' @export
setMethod(get_background_peaks, c(object = "SummarizedExperiment"), 
          function(object, 
                   bias = rowData(object)$bias, 
                   niterations = 50, 
                   w = 0.1, 
                   bs = 50) {
            object <- counts_check(object)
            get_background_peaks_core(object, bias, niterations, w, bs)
          })

#' @describeIn get_background_peaks method for RangedSummarizedExperiment
#' @export
setMethod(get_background_peaks, c(object = "RangedSummarizedExperiment"), 
          function(object, 
                   bias = rowRanges(object)$bias, 
                   niterations = 50, 
                   w = 0.1, 
                   bs = 50) {
            object <- counts_check(object)
            get_background_peaks_core(object, bias, niterations, w, bs)
          })

#' @describeIn get_background_peaks method for Matrix or matrix
#' @export
setMethod(get_background_peaks, c(object = "MatrixOrmatrix"), 
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
  
  fragments_per_peak <- get_fragments_per_peak(object)
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
  
  bin_data <- do.call(rbind, lapply(1:bs, 
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

#' get_permuted_data
#'
#' Function to get permuted data while maintaining biases
#' @param object SummarizedExperiment
#' @param niterations number of background peaks to sample
#' @param w parameter controlling similarity of background peaks
#' @param bs bin size parameter
#' @return matrix with one row per peak and one column per iteration.  values in a row
#' represent indices of background peaks for the peak with that index
#' @details Background peaks are chosen by sampling peaks based on similarity in
#' GC content and # of fragments across samples using the Mahalanobis distance.
#' The w paramter controls how similar background peaks should be.  The bs parameter
#' controls the precision with which the similarity is computed; increasing bs will
#' make the function run slower. Sensible default parameters are chosen for both.
#' @export
get_permuted_data <- function(object, niterations = 10, w = 0.1, bs = 50) {
  
  out <- BiocParallel::bplapply(1:niterations, 
                                function(x) get_permuted_data_helper(object, 
                                                                     w, bs))
  names(out) <- paste("perm", 1:niterations, sep = "_")
  
  assays(object) <- c(assays(object), out)
  
  return(object)
}


get_permuted_data_helper <- function(object, w, bs) {
  bgpeaks <- get_background_peaks(object, ncol(object), w, bs)
  reorder_columns(assays(object)$counts, bgpeaks)
}

reorder_columns <- function(mat, colixmat) {
  reordered <- lapply(1:ncol(mat), function(x) mat[colixmat[, x], x, 
                                                   drop = FALSE])
  return(do.call(cBind, reordered))
}




