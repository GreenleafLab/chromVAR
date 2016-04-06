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
  if (min(countsum$fragments_per_peak)<=0) stop("All peaks must have at least one fragment in one sample")
  if (!with_replacement && niterations > window){
    stop("If with_replacement is FALSE, then niterations must be less than window")
  }
  if (count){
    norm_mat <- cbind(log10(countsum$fragments_per_peak), bias)
  } else{
    norm_mat <- cbind(countsum$fragments_per_peak, bias)
  }
  chol_cov_mat <- chol(cov(norm_mat))
  tmp_vals <- t(forwardsolve(t(chol_cov_mat),t(norm_mat)))
  
  grpsize <- 2000
  grps <- lapply(1:(countsum$npeak %/% grpsize + ((countsum$npeak %% grpsize)!=0)), function(x) ((x-1)*grpsize +1):(min(x*grpsize,countsum$npeak)))
  
  bghelper <- function(grp, w, tmp_vals, with_replacement, niterations){
    in2 = tmp_vals[grp,]
    tmp_nns <- w$query(in2, k = window + 1, eps = 0)$nn.idx
    if (niterations == 1){
      return(matrix(sapply(1:nrow(tmp_nns), function(x) sample(tmp_nns[x,][tmp_nns[x,] != grp[x]], niterations, replace = with_replacement)),
                    ncol = 1))
    } else{
      return(t(sapply(1:nrow(tmp_nns), function(x) sample(tmp_nns[x,][tmp_nns[x,] != grp[x]], niterations, replace = with_replacement))))
    }
  }
  
  knn_kd_tree = nabor::WKNND(tmp_vals)
  background_peaks <- do.call(rbind, BiocParallel::bplapply(grps, bghelper, knn_kd_tree, tmp_vals, with_replacement, niterations))
  
  return(background_peaks)
}
