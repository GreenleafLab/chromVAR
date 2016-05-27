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
  mcols(peaks)$gc <- gc 
  return(peaks)
}


#' get_background_peaks
#' 
#' Function to get a set of background peaks for each peak based on GC content and # of counts
#' across all samples
#' @param counts_mat counts matrix
#' @param bias vector with bias values for peaks, typically GC content
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
get_background_peaks <- function(counts_mat, bias, niterations = 50, w = 0.1, bs = 50){
  
  if (inherits(bias, "GenomicRanges")){
    bias = mcols(bias)$gc
  }
  
  countsum <- counts_summary(counts_mat)
  if (min(countsum$fragments_per_peak)<=0) stop("All peaks must have at least one fragment in one sample")
  
  intensity = log10(countsum$fragments_per_peak)
  norm_mat = matrix(c(intensity, bias), ncol = 2, byrow = FALSE)
  
  chol_cov_mat <- chol(cov(norm_mat))
  trans_norm_mat <- t(forwardsolve(t(chol_cov_mat),t(norm_mat)))  
  
  #make bins
  bins1 = seq(min(trans_norm_mat[,1]), max(trans_norm_mat[,1]), length.out = bs)
  bins2 = seq(min(trans_norm_mat[,2]), max(trans_norm_mat[,2]), length.out = bs)
  
  bin_data = do.call(rbind,lapply(1:50, function(x) matrix(c(rep(bins1[x], bs),bins2), ncol = 2, byrow = FALSE)))
  
  bin_dist = euc_dist(bin_data)
  bin_p = dnorm(bin_dist, 0, w)
  
  bin_membership <- nabor::knn(bin_data, query = trans_norm_mat, k = 1)$nn.idx
  
  bin_density <- tabulate2(bin_membership, min_val = 1, max_val = bs**2)
  
  background_peaks <- bg_sample_helper(bin_membership-1, bin_p, bin_density, niterations)
  
  return(background_peaks)
}


get_background_peaks_old <- function(counts_mat, bias, niterations = 50, window = 500, with_replacement = TRUE){
  
  if (inherits(bias, "GenomicRanges")){
    bias = mcols(bias)$gc
  }
  
  countsum <- counts_summary(counts_mat)
  if (min(countsum$fragments_per_peak)<=0) stop("All peaks must have at least one fragment in one sample")
  if (!with_replacement && niterations > window){
    stop("If with_replacement is FALSE, then niterations must be less than window")
  }
  
  norm_mat <- cbind(countsum$fragments_per_peak, bias)
  
  chol_cov_mat <- chol(cov(norm_mat))
  tmp_vals <- t(forwardsolve(t(chol_cov_mat),t(norm_mat)))
  
  grpsize <- 2000
  grps <- lapply(1:(countsum$npeak %/% grpsize + ((countsum$npeak %% grpsize)!=0)), function(x) ((x-1)*grpsize +1):(min(x*grpsize,countsum$npeak)))
  
  bghelper <- function(grp, tmp_vals, with_replacement, niterations){
    tmp_nns <- nabor::knn(tmp_vals, tmp_vals[grp,], window + 1, eps = 0)$nn.idx
    if (niterations == 1){
      return(matrix(sapply(1:nrow(tmp_nns), function(x) sample(tmp_nns[x,][tmp_nns[x,] != grp[x]], niterations, replace = with_replacement)),
                    ncol = 1))
    } else{
      return(t(sapply(1:nrow(tmp_nns), function(x) sample(tmp_nns[x,][tmp_nns[x,] != grp[x]], niterations, replace = with_replacement))))
    }
  }
  
  background_peaks <- do.call(rbind, BiocParallel::bplapply(grps, bghelper, tmp_vals, with_replacement, niterations))
  
  return(background_peaks)
}
  
