#' backgroundPeaks
#'
#' backgroundPeaks is a class to store background peaks. The associated peaks are also stored
#' @slot peaks
#' @slot background_peaks
backgroundPeaks <- setClass("backgroundPeaks",
                            slots = c(peaks = 'GenomicRanges',
                                      background_peaks = 'matrix'
                            ))

#' sampleBackgroundPeaks
#' 
#' gets sample counts for background sets for a given peak set
#' @param object backgroundPeaks object
#' @param peak_set list of peak indices
#' @param counts_mat fragmentCounts object
#' @param niterations number of background sets to sample
#' @return matrix with sampled counts for background sets
#' @export
sampleBackgroundPeaks <- function(object, peak_set, counts_mat, niterations = 50){
  
  stopifnot(inherits(object, "backgroundPeaks"))
  stopifnot(niterations <= ncol(object@background_peaks))

  sample_mat = sparseMatrix(j = as.vector(object@background_peaks[peak_set,1:niterations]), i = rep(1:niterations, each = length(peak_set)), x=1, dims = c(niterations, length(counts_mat@peaks)))

  sampled_counts =  as.matrix(t(sample_mat %*% counts_mat@counts))

  return(sampled_counts)
}


#' getBackgroundPeakSets
#' 
#' Function to get a set of background peaks for each peak based on GC content and # of counts
#' across all samples
#' @param counts_mat a fragmentCounts object
#' @param niterations number of background peaks to sample
#' @param window window size around peak in which to sample a background peak
#' @param BPPARAM multiprocessing parameter for use by BiocParallel
#' @param for_set computing background peaks for a set of peaks? affect whether sampling is done with replacement
#' @return backgroundPeaks object
#' @details Background peaks are chosen by finding the window nearest neighbors to a peak in terms of 
#' GC content and # of fragments across samples using the Mahalanobis distance.  From those nearest
#' neighbors, niterations peaks are sampled from that background.
#' @export
getBackgroundPeakSets <- function(counts_mat, niterations = 50, window = 1000, BPPARAM = BiocParallel::bpparam(), for_set = TRUE){
  stopifnot(inherits(counts_mat,"fragmentCounts"))
  #if bias not available...
  if ("bias" %ni% colnames(S4Vectors::mcols(counts_mat@peaks))){
    stop("Peaks must have metadata column named 'bias'. Compute using compute_bias.")
  }
  if (!for_set && niterations > window){
    stop("If for_set is FALSE, then niterations must be less than window")
  }
  norm_mat <- cbind(log10(counts_mat@fragments_per_peak + 1 ), counts_mat@peaks$bias)
  reflected <- add_reflections(norm_mat, window = window)
  cov_mat = cov(norm_mat)
  sample_peaks <- function(index){
    mdist <- mahalanobis(reflected$data, norm_mat[index,], cov = cov_mat)
    mdist[reflected$identify(index)] = Inf
    closest <- which(mdist <= quantile(mdist, window/nrow(norm_mat), type=1))
    s <- sample(closest, niterations, replace = for_set)
    s <- reflected$replace(s)
    return(s)
  }
  if (niterations > 1){
    sampled_peaks = do.call(rbind,BiocParallel::bplapply(1:nrow(norm_mat), sample_peaks, BPPARAM = BPPARAM))
  }
  else if (niterations == 1){
    sampled_peaks = matrix(do.call(c,BiocParallel::bplapply(1:nrow(norm_mat), sample_peaks, BPPARAM = BPPARAM)), nrow = nrow(norm_mat), ncol =1)
  }
  background_peaks = backgroundPeaks(background_peaks = sampled_peaks, peaks = counts_mat@peaks)
  return(background_peaks)
}

# internal function, not exported
add_reflections <- function(x, window = 1000, as_ranks = FALSE){
  if (as_ranks){
    r = x
  } else{
    r <- apply(x,2, rank, ties.method="random")
  }
  half_window = window %/% 2
  n = nrow(x)
  mapping = c()
  out = x
  for (i in 1:ncol(x)){
    low = which(r[,i] <= half_window)
    high = which(r[,i] > n - half_window)
    low_add = x[low,]
    low_add[,i] = min(x[,i]) - (x[low,i] - min(x[,i]))
    high_add = x[high,]
    high_add[,i] = max(x[,i]) + (max(x[,i]) - x[high,i])
    out = rbind(out, low_add, high_add)
    mapping = c(mapping, low, high)
  }
  replace_reflected <- function(y){
    tmp <- which(y %in% (n+1):(n + half_window * ncol(x) *2))
    y <- replace(y, tmp, mapping[y[tmp]-n])
    return(y)
  }
  identify_same <- function(y){
    c(y, which(mapping == y) + n)
  }
  return(list("data" = out, "replace" = replace_reflected, "identify" = identify_same))
}






