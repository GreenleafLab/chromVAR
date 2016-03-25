#' get_gc
#' 
#' Computes GC content for peaks
#' @param peaks GenomicRanges
#' @param genome BSgenome object
#' @export
get_gc <- function(peaks, 
                   genome = BSgenome.Hsapiens.UCSC.hg19::BSgenome.Hsapiens.UCSC.hg19){
  seqs = Biostrings::getSeq(genome, peaks)
  nucfreqs <- Biostrings::letterFrequency(seqs, c("A","C","G","T"))
  gc <- apply(nucfreqs, 1, function(x) sum(x[2:3])/sum(x))
  return(gc)
}

#' get_background_peaks
#' 
#' Function to get a set of background peaks for each peak based on GC content and # of counts
#' across all samples
#' @param counts_mat counts matrix
#' @param bias vector with bias values for peaks, typically GC content
#' @param niterations number of background peaks to sample
#' @param window window size around peak in which to sample a background peak
#' @param with_replacement sample peaks with replacement? Default is TRUE
#' @return matrix with one row per peak and one column per iteration.  values in a row
#' represent indices of background peaks for the peak with that index
#' @details Background peaks are chosen by finding the window nearest neighbors to a peak in terms of 
#' GC content and # of fragments across samples using the Mahalanobis distance.  From those nearest
#' neighbors, niterations peaks are sampled from that background.
#' @export
get_background_peaks <- function(counts_mat, bias, niterations = 50, window = 500, with_replacement = TRUE, count = TRUE){

  countsum <- counts_summary(counts_mat)
  if (min(countsum$fragments_per_peaks <= 0)) stop("All peaks must have at least one fragment in one sample")
  if (!with_replacement && niterations > window){
    stop("If with_replacement is FALSE, then niterations must be less than window")
  }
  if (count){
    norm_mat <- cbind(log10(countsum$fragments_per_peak + 1 ), bias)
  } else{
    norm_mat <- cbind(countsum$fragments_per_peak, bias)
  }
  reflected <- add_reflections(norm_mat, window = window)
  chol_cov_mat <- chol(cov(norm_mat))
  tmp_vals <- t(forwardsolve(t(chol_cov_mat),t(reflected$data)))
  
  grpsize <- 10000
  grps <- lapply(1:(countsum$npeak %/% grpsize + ((countsum$npeak %% grpsize)!=0)), function(x) ((x-1)*grpsize +1):(min(x*grpsize,countsum$npeak)))
  
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
  
  background_peaks <- do.call(rbind, BiocParallel::bplapply(grps, bghelper, reflected, tmp_vals, with_replacement, niterations))

  return(background_peaks)
}


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

