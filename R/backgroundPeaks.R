#' sampleBackgroundPeaks
#' 
#' gets sample counts for background sets for a given peak set
#' @param object fragmentCounts object
#' @param peak_set list of peak indices
#' @param niterations number of background sets to sample
#' @details input fragmentCounts object must include set of background peaks, computed
#' using \code{\link{getBackgroundPeakSets}}
#' @return matrix with sampled counts for background sets
#' @export
sampleBackgroundPeaks <- function(object, peak_set, niterations = 50){
  
  stopifnot(inherits(object, "fragmentCounts"))
  stopifnot(niterations <= ncol(object@background_peaks))

  sample_mat = sparseMatrix(j = unlist(object@background_peaks[peak_set,1:niterations],use.names=FALSE), 
                            i = rep(1:niterations, each = length(peak_set)), 
                            x=1, 
                            dims = c(niterations, object@npeak))

  sampled_counts =  as.matrix(t(sample_mat %*% object@counts))

  return(sampled_counts)
}


#' getBackgroundPeakSets
#' 
#' Function to get a set of background peaks for each peak based on GC content and # of counts
#' across all samples
#' @param counts_mat a fragmentCounts object
#' @param niterations number of background peaks to sample
#' @param window window size around peak in which to sample a background peak
#' @param with_replacement sample peaks with replacement? Default is TRUE
#' @return fragmentCounts object
#' @details Background peaks are chosen by finding the window nearest neighbors to a peak in terms of 
#' GC content and # of fragments across samples using the Mahalanobis distance.  From those nearest
#' neighbors, niterations peaks are sampled from that background.
#' @export
getBackgroundPeakSets <- function(counts_mat, niterations = 50, window = 500, with_replacement = TRUE, count = TRUE){
  stopifnot(inherits(counts_mat,"fragmentCounts"))
  #if bias not available...
  if ("bias" %ni% colnames(S4Vectors::mcols(counts_mat@peaks))){
    stop("Peaks in fragmentCounts must have metadata column named 'bias'. Compute using compute_bias.")
  }
  if (!with_replacement && niterations > window){
    stop("If with_replacement is FALSE, then niterations must be less than window")
  }
  if (count){
    norm_mat <- cbind(log10(counts_mat@fragments_per_peak + 1 ), counts_mat@peaks$bias)
  } else{
    norm_mat <- cbind(counts_mat@fragments_per_peak, counts_mat@peaks$bias)
  }
  reflected <- add_reflections(norm_mat, window = window)
  chol_cov_mat <- chol(cov(norm_mat))
  tmp_vals <- t(forwardsolve(t(chol_cov_mat),t(reflected$data)))
  
  grpsize <- 10000
  grps <- lapply(1:(counts_mat@npeak %/% grpsize + ((counts_mat@npeak %% grpsize)!=0)), function(x) ((x-1)*grpsize +1):(min(x*grpsize,counts_mat@npeak)))
  
  bghelper <- function(grp, reflected, tmp_vals, with_replacement, niterations){
    in1 = tmp_vals
    in2 = tmp_vals[grp,]
    tmp_nns <- FNN::get.knnx(in1, query = in2, k = window)$nn.index
    if (niterations == 1){
        return(matrix(apply(tmp_nns, 1, function(x) reflected$replace(sample(x, niterations, replace = with_replacement))),
                                                                 ncol = 1))
    } else{
        return(t(apply(tmp_nns,1, function(x) reflected$replace(sample(x, niterations, replace = with_replacement)))))
    }
  }
  
  counts_mat@background_peaks <- do.call(rbind, BiocParallel::bplapply(grps, bghelper, reflected, tmp_vals, with_replacement, niterations))
  return(counts_mat)
}
  
#   
#   nns <- FNN::get.knn(tmp_vals, k = window)$nn.index
#   if (niterations == 1){
#     counts_mat@background_peaks <- matrix(apply(nns[1:length(counts_mat@peaks),], 1, function(x) reflected$replace(sample(x, niterations, replace = with_replacement))),
#                                           ncol = 1)
#   } else{
#     counts_mat@background_peaks <- t(apply(nns[1:length(counts_mat@peaks),], 1, function(x) reflected$replace(sample(x, niterations, replace = with_replacement))))
#   }
#   return(counts_mat)
#}


# helper function, not exported
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






